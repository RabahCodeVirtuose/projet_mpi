#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <vector>
#include <algorithm>
#include <cmath>

#include "ParallelFWBlocks.hpp"
#include "Distribution.hpp"   // BlockInfo, ownerOf, computeLocalBlocks

using namespace std;

// "Infini" suffisamment grand
static const int INF = 1000000000;

// -------------------------------------------------------------------
//  Helpers bloc → mini-Floyd sur des blocs (b comme leading dimension)
// -------------------------------------------------------------------

// Bloc diagonal : Dkk (bs x bs), stocké dans un tableau de taille b*b
// bs = taille réelle du bloc (utile sur la dernière diagonale si n n'est
//     pas multiple de b)
// b  = leading dimension (nombre de colonnes allouées)
static void fw_block(int* Dkk, int bs, int b) {
    for (int kk = 0; kk < bs; ++kk) {
        for (int i = 0; i < bs; ++i) {
            int dik = Dkk[i * b + kk];
            if (dik == INF) continue;
            for (int j = 0; j < bs; ++j) {
                int kkj = Dkk[kk * b + j];
                if (kkj == INF) continue;

                int via = dik + kkj;
                int& dij = Dkk[i * b + j];
                if (via < dij) {
                    dij = via;
                }
            }
        }
    }
}

// Bloc de ligne : DkJ (bs x wJ) mis à jour avec le pivot Dkk (bs x bs)
// Dkk : bs x bs
// DkJ : bs x wJ
// b   : leading dimension
static void fw_row(const int* Dkk, int* DkJ, int bs, int wJ, int b) {
    for (int i = 0; i < bs; ++i) {
        for (int kk = 0; kk < bs; ++kk) {
            int dik = Dkk[i * b + kk];
            if (dik == INF) continue;

            for (int j = 0; j < wJ; ++j) {
                int kkj = DkJ[kk * b + j];
                if (kkj == INF) continue;

                int via = dik + kkj;
                int& dij = DkJ[i * b + j];
                if (via < dij) {
                    dij = via;
                }
            }
        }
    }
}

// Bloc de colonne : Dik (hI x bs) mis à jour avec le pivot Dkk (bs x bs)
// Dik : hI x bs
// Dkk : bs x bs
// b   : leading dimension
static void fw_col(int* Dik, const int* Dkk, int hI, int bs, int b) {
    for (int i = 0; i < hI; ++i) {
        for (int kk = 0; kk < bs; ++kk) {
            int ik = Dik[i * b + kk];
            if (ik == INF) continue;

            for (int j = 0; j < bs; ++j) {
                int kj = Dkk[kk * b + j];
                if (kj == INF) continue;

                int via = ik + kj;
                int& ij = Dik[i * b + j];
                if (via < ij) {
                    ij = via;
                }
            }
        }
    }
}

// Bloc "interne" : Dij (hI x wJ) mis à jour avec Dik (hI x bs) et DkJ (bs x wJ)
// Dik : hI x bs
// DkJ : bs x wJ
// Dij : hI x wJ
// b   : leading dimension
static void fw_inner(const int* Dik, const int* DkJ,
                     int* Dij, int hI, int wJ, int bs, int b) {
    for (int i = 0; i < hI; ++i) {
        for (int kk = 0; kk < bs; ++kk) {
            int ik = Dik[i * b + kk];
            if (ik == INF) continue;

            for (int j = 0; j < wJ; ++j) {
                int kj = DkJ[kk * b + j];
                if (kj == INF) continue;

                int via = ik + kj;
                int& ij = Dij[i * b + j];
                if (via < ij) {
                    ij = via;
                }
            }
        }
    }
}

// -------------------------------------------------------------------
//  Version Floyd–Warshall MPI par blocs B×B (2D)
// -------------------------------------------------------------------
int* ParallelFloydWarshallBlocks(int n, int* mat) {

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // -------------------------------------------------------
    // 0) Paramètres de la grille de processus + taille de bloc
    // -------------------------------------------------------
    int dims[2] = {0, 0};
    MPI_Dims_create(size, 2, dims);  // calcule Pr, Pc "équilibrés"
    int Pr = dims[0];
    int Pc = dims[1];

    // Choix "simple" d'une taille de bloc :
    // approx ~ n / sqrt(p)
    int approxBlocks = (int)std::sqrt((double)size);
    if (approxBlocks < 1) approxBlocks = 1;
    int b = (n + approxBlocks - 1) / approxBlocks;  // taille théorique du bloc

    if (b < 1) b = 1;

    // Nombre de blocs par dimension
    int nb = (n + b - 1) / b;  // ceil(n/b)

    // -------------------------------------------------------
    // 1) Quels blocs appartiennent à ce rang ?
    // -------------------------------------------------------
    vector<BlockInfo> localBlocks = computeLocalBlocks(n, b, Pr, Pc, rank);
    int numLocal = (int)localBlocks.size();
    int blockArea = b * b;

    // Tableau local : concaténation de tous les blocs locaux
    // bloc idx → pointeur = &localData[idx * blockArea]
    vector<int> localData(numLocal * blockArea);

    // Table de correspondance (bi,bj) → index local ou -1
    vector<int> localIndex(nb * nb, -1);
    for (int idx = 0; idx < numLocal; ++idx) {
        const BlockInfo& info = localBlocks[idx];
        localIndex[info.bi * nb + info.bj] = idx;
    }

    // -------------------------------------------------------
    // 2) Initialisation D^(-1) dans les blocs locaux
    // -------------------------------------------------------
    for (int idx = 0; idx < numLocal; ++idx) {
        const BlockInfo& info = localBlocks[idx];
        int i0 = info.offset_i;
        int j0 = info.offset_j;

        int* blk = &localData[idx * blockArea];

        for (int ii = 0; ii < b; ++ii) {
            int gi = i0 + ii;        // indice global i
            for (int jj = 0; jj < b; ++jj) {
                int gj = j0 + jj;    // indice global j

                int val;
                if (gi >= n || gj >= n) {
                    // En dehors de la matrice globale
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

    // Buffers pour la ligne et la colonne de blocs du pivot k
    // rowBlocks[j] contiendra le bloc (k, j) pour tous j
    // colBlocks[i] contiendra le bloc (i, k) pour tous i
    vector<vector<int>> rowBlocks(nb, vector<int>(blockArea));
    vector<vector<int>> colBlocks(nb, vector<int>(blockArea));

    // -------------------------------------------------------
    // 3) Boucle principale sur les blocs de pivot kk
    // -------------------------------------------------------
    for (int kk = 0; kk < nb; ++kk) {

        int pivotOwner = ownerOf(kk, kk, Pr, Pc);
        int pivotLocalIdx = localIndex[kk * nb + kk];
        vector<int> pivotBlock(blockArea, INF);

        // Taille réelle du bloc diagonal (peut être < b au bord)
        int bs = std::min(b, n - kk * b);

        // ---- Phase A : bloc pivot D[kk,kk] ----
        if (rank == pivotOwner && pivotLocalIdx != -1) {
            int* DkkLocal = &localData[pivotLocalIdx * blockArea];
            fw_block(DkkLocal, bs, b);              // met à jour Dkk local
            std::copy(DkkLocal, DkkLocal + blockArea, pivotBlock.begin());
        }

        // Diffusion du bloc pivot à tous
        MPI_Bcast(pivotBlock.data(), blockArea, MPI_INT, pivotOwner, MPI_COMM_WORLD);

        // rowBlocks[kk] et colBlocks[kk] pourront référencer le pivot
        rowBlocks[kk] = pivotBlock;
        colBlocks[kk] = pivotBlock;

        // ---- Phase B.1 : mise à jour + diffusion des blocs de ligne (kk, j) ----
        for (int jb = 0; jb < nb; ++jb) {
            if (jb == kk) continue;  // déjà le pivot

            int ownerRow = ownerOf(kk, jb, Pr, Pc);
            vector<int> tempBlock(blockArea, INF);

            int wJ = std::min(b, n - jb * b); // largeur réelle de ce bloc

            if (rank == ownerRow) {
                int localIdx = localIndex[kk * nb + jb];
                if (localIdx != -1) {
                    int* DkJ = &localData[localIdx * blockArea];
                    fw_row(pivotBlock.data(), DkJ, bs, wJ, b);
                    std::copy(DkJ, DkJ + blockArea, tempBlock.begin());
                }
            }

            MPI_Bcast(tempBlock.data(), blockArea, MPI_INT, ownerRow, MPI_COMM_WORLD);
            rowBlocks[jb] = tempBlock;
        }

        // ---- Phase B.2 : mise à jour + diffusion des blocs de colonne (i, kk) ----
        for (int ib = 0; ib < nb; ++ib) {
            if (ib == kk) continue;  // déjà le pivot

            int ownerCol = ownerOf(ib, kk, Pr, Pc);
            vector<int> tempBlock(blockArea, INF);

            int hI = std::min(b, n - ib * b); // hauteur réelle de ce bloc

            if (rank == ownerCol) {
                int localIdx = localIndex[ib * nb + kk];
                if (localIdx != -1) {
                    int* Dik = &localData[localIdx * blockArea];
                    fw_col(Dik, pivotBlock.data(), hI, bs, b);
                    std::copy(Dik, Dik + blockArea, tempBlock.begin());
                }
            }

            MPI_Bcast(tempBlock.data(), blockArea, MPI_INT, ownerCol, MPI_COMM_WORLD);
            colBlocks[ib] = tempBlock;
        }

        // ---- Phase C : mise à jour de tous les autres blocs locaux (i,j) ----
        for (int idx = 0; idx < numLocal; ++idx) {
            BlockInfo& info = localBlocks[idx];
            int ib = info.bi;
            int jb = info.bj;

            if (ib == kk || jb == kk) {
                // Ligne ou colonne du pivot : déjà traitée
                continue;
            }

            int hI = std::min(b, n - ib * b);
            int wJ = std::min(b, n - jb * b);

            int* Dij = &localData[idx * blockArea];
            const int* Dik = colBlocks[ib].data();
            const int* DkJ = rowBlocks[jb].data();

            fw_inner(Dik, DkJ, Dij, hI, wJ, bs, b);
        }
    }

    // -------------------------------------------------------
    // 4) Rapatrier la matrice finale sur le rang 0
    // -------------------------------------------------------
    int* D_final = nullptr;
    if (rank == 0) {
        D_final = new int[n * n];
        // initialiser à INF par sécurité
        for (int i = 0; i < n * n; ++i) D_final[i] = INF;
    }

    // Pour chaque rang, le rang 0 sait quels blocs il possède,
    // car computeLocalBlocks(...) est déterministe.
    for (int r = 0; r < size; ++r) {

        vector<BlockInfo> blocks_r = computeLocalBlocks(n, b, Pr, Pc, r);

        for (size_t idx_r = 0; idx_r < blocks_r.size(); ++idx_r) {
            const BlockInfo& info = blocks_r[idx_r];

            // Buffer pour recevoir (ou lire) un bloc complet b*b
            vector<int> buf(blockArea);

            if (rank == r) {
                // On a déjà le bloc localement
                int localIdx = localIndex[info.bi * nb + info.bj];
                if (localIdx != -1) {
                    int* src = &localData[localIdx * blockArea];
                    std::copy(src, src + blockArea, buf.begin());
                }
            }

            // Transfert vers le rang 0
            if (r == 0) {
                // Rang 0 lit directement son propre bloc : pas de MPI
                if (rank == 0) {
                    // rien de spécial, buf contient déjà les données
                }
            } else {
                if (rank == r) {
                    MPI_Send(buf.data(), blockArea, MPI_INT, 0, 0, MPI_COMM_WORLD);
                }
                if (rank == 0) {
                    MPI_Recv(buf.data(), blockArea, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }

            // Rang 0 intègre le bloc dans D_final
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

    // Sur le rang 0, on renvoie la matrice complète
    // Sur les autres rangs, on renvoie nullptr
    if (rank == 0) {
        return D_final;
    } else {
        return nullptr;
    }
}
