#include <mpi.h>
#include <vector>
#include <algorithm>
#include "ParallelFW.hpp"
// On inclut Distribution.hpp mais on ne l'utilise plus dans cette version
#include "Distribution.hpp"

using namespace std;

// On choisit un INF suffisamment grand pour représenter "pas d'arête"
static const int INF = 1000000000;

// Petit helper pour distribuer les lignes entre les processus
// rowsPerRank[r] = nombre de lignes pour le rang r
// startRow[r]    = indice de première ligne globale pour le rang r
static void computeRowDistribution(int n, int size,
                                   vector<int>& rowsPerRank,
                                   vector<int>& startRow) {
    rowsPerRank.resize(size);
    startRow.resize(size);

    int base = n / size;
    int rest = n % size;   // les "rest" premiers rangs auront une ligne en plus

    int offset = 0;
    for (int r = 0; r < size; ++r) {
        rowsPerRank[r] = base + (r < rest ? 1 : 0);
        startRow[r] = offset;
        offset += rowsPerRank[r];
    }
}

// =============================
//   Version MPI par bandes
// =============================
// - n      : nb de sommets
// - mat    : matrice d'adjacence (0 = pas d'arête, sinon poids)
// Retour : sur rank 0, pointeur vers matrice des distances n×n
//          sur les autres rangs : nullptr
int* ParallelFloydWarshall(int n, int* mat) {

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // -----------------------------
    // 1) Distribuer les lignes
    // -----------------------------
    vector<int> rowsPerRank, startRow;
    computeRowDistribution(n, size, rowsPerRank, startRow);

    int localRows = rowsPerRank[rank];       // nombre de lignes gérées par ce processus
    int firstRow  = startRow[rank];          // indice global de la première ligne locale

    // localD : petit sous-ensemble des lignes de la matrice D (distance)
    // dimension = localRows × n
    int* localD = new int[localRows * n];

    // -----------------------------
    // 2) Initialiser D^(-1) localement
    // -----------------------------
    // On part de la matrice d'adjacence "mat" (déjà présente chez tout le monde,
    // car main_mpi l'a diffusée avec un MPI_Bcast).
    //
    // Rappel de la règle (comme InitDk dans le code de la prof) :
    // - D[i][j] = 0        si i == j
    // - D[i][j] = INF      si pas d'arête (mat[i][j] == 0 et i != j)
    // - D[i][j] = poids    sinon
    for (int li = 0; li < localRows; ++li) {
        int i = firstRow + li;   // indice global de la ligne
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                localD[li * n + j] = 0;
            } else if (mat[i * n + j] == 0) {
                localD[li * n + j] = INF;
            } else {
                localD[li * n + j] = mat[i * n + j];
            }
        }
    }

    // Buffer pour contenir la ligne k (dimension n)
    // Cette ligne est diffusée à tous les processus à chaque itération k
    int* row_k = new int[n];

    // -----------------------------
    // 3) Boucle principale de Floyd–Warshall
    // -----------------------------
    // Version classique :
    // for k in 0..n-1:
    //   for i in 0..n-1:
    //     for j in 0..n-1:
    //       D[i][j] = min(D[i][j], D[i][k] + D[k][j])
    //
    // Ici : chaque processus a un paquet de lignes i, donc il met à jour
    // uniquement ces lignes, en utilisant :
    //  - localD[li][k] (distance i->k)
    //  - row_k[j]      (distance k->j) diffusée depuis le processus qui possède la ligne k
    // -----------------------------
    for (int k = 0; k < n; ++k) {

        // 3.1) Trouver qui possède la ligne k
        int owner = 0;
        for (int r = 0; r < size; ++r) {
            int start = startRow[r];
            int end   = startRow[r] + rowsPerRank[r] - 1;
            if (k >= start && k <= end) {
                owner = r;
                break;
            }
        }

        // 3.2) Le propriétaire remplit row_k avec la ligne k
        if (rank == owner) {
            int localIndex = k - firstRow;  // indice local (dans localD) de la ligne k
            for (int j = 0; j < n; ++j) {
                row_k[j] = localD[localIndex * n + j];
            }
        }

        // 3.3) Diffusion de la ligne k à tous les processus
        MPI_Bcast(row_k, n, MPI_INT, owner, MPI_COMM_WORLD);

        // 3.4) Mise à jour des lignes locales avec la ligne k
        for (int li = 0; li < localRows; ++li) {
            int idx_ik = localD[li * n + k];   // D[i][k] pour la ligne locale i
            if (idx_ik == INF) continue;       // si i->k est INF, inutile de tester i->k->j

            for (int j = 0; j < n; ++j) {
                int dkj = row_k[j];            // D[k][j]
                int& dij = localD[li * n + j]; // D[i][j]

                if (dkj == INF) continue;      // k->j impossible, pas d'amélioration

                int via_k = idx_ik + dkj;
                if (via_k < dij) {
                    dij = via_k;
                }
            }
        }
    }

    // -----------------------------
    // 4) Rassembler la matrice finale sur le rang 0
    // -----------------------------
    int* D_final = nullptr;
    vector<int> recvCounts(size), displs(size);

    for (int r = 0; r < size; ++r) {
        recvCounts[r] = rowsPerRank[r] * n;
    }
    displs[0] = 0;
    for (int r = 1; r < size; ++r) {
        displs[r] = displs[r-1] + recvCounts[r-1];
    }

    if (rank == 0) {
        D_final = new int[n * n];
    }

    MPI_Gatherv(localD,
                localRows * n, MPI_INT,
                D_final,
                recvCounts.data(), displs.data(), MPI_INT,
                0, MPI_COMM_WORLD);

    delete[] localD;
    delete[] row_k;

    // Sur rank 0, on renvoie la matrice complète D
    // Sur les autres rangs, on renvoie nullptr (inutile d'avoir une copie partout)
    if (rank == 0) {
        return D_final;
    } else {
        return nullptr;
    }
}
