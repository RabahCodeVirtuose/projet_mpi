#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <omp.h>
#include "ParallelFWBlocks.hpp"
#include "Distribution.hpp"

/**
 * @file ParallelFWBlocks.cpp
 * @brief Implémentation du Floyd–Warshall parallèle par blocs (MPI + OpenMP).
 */

using namespace std;

static const int INF = 1000000000;

// ============================================================================
// HELPERS
// ============================================================================

/**
 * @brief Applique Floyd–Warshall sur le bloc pivot Dkk (taille b×b).
 *
 * @param Dkk Pointeur vers le bloc pivot local.
 * @param bs  Taille réelle du bloc (peut être < b en bord).
 * @param b   Taille logique du bloc.
 */
static void fw_block(int* Dkk, int bs, int b) {
    #pragma omp parallel
    {
        for (int kk = 0; kk < bs; ++kk) {
            #pragma omp for schedule(static)
            for (int i = 0; i < bs; ++i) {
                int dik = Dkk[i * b + kk];
                if (dik == INF) continue;
                for (int j = 0; j < bs; ++j) {
                    int kkj = Dkk[kk * b + j];
                    if (kkj == INF) continue;
                    int via = dik + kkj;
                    int& dij = Dkk[i * b + j];
                    if (via < dij) dij = via;
                }
            }
        }
    }
}

/**
 * @brief Met à jour un bloc D(k,J) de la ligne du pivot à partir de D(k,k).
 *
 * @param Dkk Bloc pivot D(k,k).
 * @param DkJ Bloc cible D(k,J).
 * @param bs  Taille réelle du pivot.
 * @param wJ  Largeur réelle du bloc D(k,J).
 * @param b   Taille logique du bloc.
 */
static void fw_row(const int* Dkk, int* DkJ, int bs, int wJ, int b) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < bs; ++i) {
        for (int kk = 0; kk < bs; ++kk) {
            int dik = Dkk[i * b + kk];
            if (dik == INF) continue;
            for (int j = 0; j < wJ; ++j) {
                int kkj = DkJ[kk * b + j];
                if (kkj == INF) continue;
                int via = dik + kkj;
                int& dij = DkJ[i * b + j];
                if (via < dij) dij = via;
            }
        }
    }
}

/**
 * @brief Met à jour un bloc D(I,k) de la colonne du pivot à partir de D(k,k).
 *
 * @param Dik Bloc cible D(I,k).
 * @param Dkk Bloc pivot D(k,k).
 * @param hI  Hauteur réelle du bloc D(I,k).
 * @param bs  Taille réelle du pivot.
 * @param b   Taille logique du bloc.
 */
static void fw_col(int* Dik, const int* Dkk, int hI, int bs, int b) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < hI; ++i) {
        for (int kk = 0; kk < bs; ++kk) {
            int ik = Dik[i * b + kk];
            if (ik == INF) continue;
            for (int j = 0; j < bs; ++j) {
                int kj = Dkk[kk * b + j];
                if (kj == INF) continue;
                int via = ik + kj;
                int& ij = Dik[i * b + j];
                if (via < ij) ij = via;
            }
        }
    }
}

/**
 * @brief Met à jour un bloc interne D(I,J) à partir de D(I,k) et D(k,J).
 *
 * @param Dik Bloc D(I,k).
 * @param DkJ Bloc D(k,J).
 * @param Dij Bloc D(I,J) à mettre à jour.
 * @param hI  Hauteur réelle du bloc D(I,J).
 * @param wJ  Largeur réelle du bloc D(I,J).
 * @param bs  Taille réelle du pivot.
 * @param b   Taille logique du bloc.
 */
static void fw_inner(const int* Dik, const int* DkJ, int* Dij, 
                     int hI, int wJ, int bs, int b) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < hI; ++i) {
        for (int kk = 0; kk < bs; ++kk) {
            int ik = Dik[i * b + kk];
            if (ik == INF) continue;
            for (int j = 0; j < wJ; ++j) {
                int kj = DkJ[kk * b + j];
                if (kj == INF) continue;
                int via = ik + kj;
                int& ij = Dij[i * b + j];
                if (via < ij) ij = via;
            }
        }
    }
}

// ============================================================================
// FONCTION PRINCIPALE
// ============================================================================

/**
 * @brief Applique Floyd–Warshall par blocs sur une matrice d'adjacence.
 *
 * @param n   Taille de la matrice (n × n).
 * @param mat Matrice d'adjacence d'entrée (row-major).
 * @return Matrice finale des distances sur le rang 0, nullptr sinon.
 */
int* ParallelFloydWarshallBlocks(int n, int* mat) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Calcul de la taille des blocs
    int sqrtp = (int)lround(sqrt((double)size));
    bool grilleCarree = (sqrtp * sqrtp == size) && (sqrtp > 0) && (n % sqrtp == 0);

    int b;
    if (grilleCarree) {
        // Cas idéal : grille carrée
        b = n / sqrtp;
        // On limite b pour éviter des blocs trop gros qui ne rentrent pas en cache
        if (b > 128) b = 128;
    } else {
        // Cas général : on calcule un b raisonnable
        int denom = max(1, sqrtp);
        b = (n + denom - 1) / denom;
        // On garde b entre 32 et 128 pour de bonnes perfs
        if (b < 32) b = 32;
        if (b > 128) b = 128;
    }

    int nb = (n + b - 1) / b;

    if (rank == 0) {
        cout << "[INFO] Taille matrice : " << n << "x" << n << endl;
        cout << "[INFO] Taille bloc    : " << b << "x" << b << endl;
        cout << "[INFO] Nombre blocs   : " << nb << "x" << nb << endl;
        cout << "[INFO] Processus      : " << size << endl;
    }

    // Construction de la grille de processus
    int dims[2] = {0, 0};
    if (grilleCarree) {
        dims[0] = sqrtp;
        dims[1] = sqrtp;
    }
    MPI_Dims_create(size, 2, dims);
    int Pr = dims[0];
    int Pc = dims[1];

    // Distribution des blocs
    vector<BlockInfo> localBlocks = computeLocalBlocks(n, b, Pr, Pc, rank);
    int numLocal = (int)localBlocks.size();
    int blockArea = b * b;

    vector<int> localData(numLocal * blockArea);
    vector<int> localIndex(nb * nb, -1);
    
    for (int idx = 0; idx < numLocal; ++idx) {
        const BlockInfo& info = localBlocks[idx];
        localIndex[info.bi * nb + info.bj] = idx;
    }

    // Initialisation des blocs locaux
    for (int idx = 0; idx < numLocal; ++idx) {
        const BlockInfo& info = localBlocks[idx];
        int i0 = info.offset_i;
        int j0 = info.offset_j;
        int* blk = &localData[idx * blockArea];

        for (int ii = 0; ii < b; ++ii) {
            int gi = i0 + ii;
            for (int jj = 0; jj < b; ++jj) {
                int gj = j0 + jj;
                int val;
                if (gi >= n || gj >= n) {
                    val = INF;
                } else if (gi == gj) {
                    val = 0;
                } else if (mat[gi * n + gj] == 0) {
                    val = INF;
                } else {
                    val = mat[gi * n + gj];
                }
                blk[ii * b + jj] = val;
            }
        }
    }

    // Buffers pour les communications
    vector<vector<int>> rowBlocks(nb, vector<int>(blockArea));
    vector<vector<int>> colBlocks(nb, vector<int>(blockArea));

    // Boucle principale sur les blocs pivots
    for (int kk = 0; kk < nb; ++kk) {
        int pivotOwner = ownerOf(kk, kk, Pr, Pc);
        int pivotLocalIdx = localIndex[kk * nb + kk];

        vector<int> pivotBlock(blockArea, INF);
        int bs = min(b, n - kk * b);

        // Phase A : Traitement du bloc pivot
        if (rank == pivotOwner && pivotLocalIdx != -1) {
            int* DkkLocal = &localData[pivotLocalIdx * blockArea];
            fw_block(DkkLocal, bs, b);
            copy(DkkLocal, DkkLocal + blockArea, pivotBlock.begin());
        }

        MPI_Bcast(pivotBlock.data(), blockArea, MPI_INT, pivotOwner, MPI_COMM_WORLD);

        rowBlocks[kk] = pivotBlock;
        colBlocks[kk] = pivotBlock;

        // Phase B : Mise à jour ligne et colonne
        // J'ai remarqué qu'on peut lancer les Ibcast de la ligne et de la colonne
        // en même temps au lieu de faire deux Waitall séparés
        vector<MPI_Request> requests;

        // Mise à jour de la ligne
        for (int jb = 0; jb < nb; ++jb) {
            if (jb == kk) continue;

            int ownerRow = ownerOf(kk, jb, Pr, Pc);
            int wJ = min(b, n - jb * b);

            if (rank == ownerRow) {
                int localIdx = localIndex[kk * nb + jb];
                if (localIdx != -1) {
                    int* DkJ = &localData[localIdx * blockArea];
                    fw_row(pivotBlock.data(), DkJ, bs, wJ, b);
                    copy(DkJ, DkJ + blockArea, rowBlocks[jb].begin());
                }
            }

            MPI_Request req;
            MPI_Ibcast(rowBlocks[jb].data(), blockArea, MPI_INT, 
                      ownerRow, MPI_COMM_WORLD, &req);
            requests.push_back(req);
        }

        // Mise à jour de la colonne (en parallèle de la ligne)
        for (int ib = 0; ib < nb; ++ib) {
            if (ib == kk) continue;

            int ownerCol = ownerOf(ib, kk, Pr, Pc);
            int hI = min(b, n - ib * b);

            if (rank == ownerCol) {
                int localIdx = localIndex[ib * nb + kk];
                if (localIdx != -1) {
                    int* Dik = &localData[localIdx * blockArea];
                    fw_col(Dik, pivotBlock.data(), hI, bs, b);
                    copy(Dik, Dik + blockArea, colBlocks[ib].begin());
                }
            }

            MPI_Request req;
            MPI_Ibcast(colBlocks[ib].data(), blockArea, MPI_INT,
                      ownerCol, MPI_COMM_WORLD, &req);
            requests.push_back(req);
        }

        // On attend que tous les Ibcast soient finis (ligne + colonne)
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

        // Phase C : Mise à jour des blocs internes
        for (int idx = 0; idx < numLocal; ++idx) {
            BlockInfo& info = localBlocks[idx];
            int ib = info.bi;
            int jb = info.bj;

            if (ib == kk || jb == kk) continue;

            int hI = min(b, n - ib * b);
            int wJ = min(b, n - jb * b);

            int* Dij = &localData[idx * blockArea];
            const int* Dik = colBlocks[ib].data();
            const int* DkJ = rowBlocks[jb].data();

            fw_inner(Dik, DkJ, Dij, hI, wJ, bs, b);
        }
    }

    // Rassemblement final
    int* D_final = nullptr;
    if (rank == 0) {
        D_final = new int[n * n];
        for (int i = 0; i < n * n; ++i) D_final[i] = INF;
    }

    // Chaque processus envoie ses blocs au rang 0
    for (int r = 0; r < size; ++r) {
        vector<BlockInfo> blocks_r = computeLocalBlocks(n, b, Pr, Pc, r);

        for (size_t idx_r = 0; idx_r < blocks_r.size(); ++idx_r) {
            const BlockInfo& info = blocks_r[idx_r];
            vector<int> buf(blockArea);
            
            if (rank == r) {
                int localIdx = localIndex[info.bi * nb + info.bj];
                if (localIdx != -1) {
                    int* src = &localData[localIdx * blockArea];
                    copy(src, src + blockArea, buf.begin());
                }
            }

            if (r == 0) {
                // Rien à faire, déjà local
            } else {
                if (rank == r) {
                    MPI_Send(buf.data(), blockArea, MPI_INT, 0, 0, MPI_COMM_WORLD);
                }
                if (rank == 0) {
                    MPI_Recv(buf.data(), blockArea, MPI_INT, r, 0, 
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }

            if (rank == 0) {
                int i0 = info.offset_i;
                int j0 = info.offset_j;
                for (int ii = 0; ii < b; ++ii) {
                    int gi = i0 + ii;
                    if (gi >= n) break;
                    for (int jj = 0; jj < b; ++jj) {
                        int gj = j0 + jj;
                        if (gj >= n) break;
                        D_final[gi * n + gj] = buf[ii * b + jj];
                    }
                }
            }
        }
    }

    return (rank == 0) ? D_final : nullptr;
}
