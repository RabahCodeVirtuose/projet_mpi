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



/* dkk est un bloc diagonal de la matrice des distances un sous tableaux 2d mais stocké en 1D on fait 
une petite version de floyd warshal uniquemet sur ce bloc la on applique la formule à l'intérieur du bloc 

première boucle : 
bs : taille réelle du bloc par ex si le bloc fait 10 * 1O , bs = 10 
kk : joue le role de k dans floyd warshal


intuition : pour chaquue sommet k du bloc je regarde ......... est ce que passer par k ameliore les distances 
entre les couples i j du bloc 

deuxième boucle : 
i indice de la lige du bloc 
kk[i * b + kk ] 

*/
static void fw_block(int* Dkk, int bs, int b) {
    // ---------------------------------------------------------------------------
// Boucle principale de Floyd–Warshall à l'intérieur d'un bloc diagonal Dkk.
//
// Objectif : appliquer la formule classique de Floyd–Warshall
//     D[i][j] = min(D[i][j], D[i][k] + D[k][j])
// mais uniquement sur le sous-bloc diagonal courant.
//
// - Dkk est stocké en 1D avec un "leading dimension" b :
//       élément (i, j) du bloc → Dkk[i * b + j]
// - bs est la taille réelle du bloc (peut être < b en bord de matrice).
// - INF représente "pas de chemin" (distance infinie).
// ---------------------------------------------------------------------------
for (int kk = 0; kk < bs; ++kk) {
    // kk joue le rôle du pivot k dans l’algorithme de Floyd–Warshall.
    // À chaque itération, on teste si passer par kk améliore les distances
    // entre tous les couples (i, j) du bloc.

    for (int i = 0; i < bs; ++i) {
        // On récupère D[i][k] local dans le bloc diagonal.
        // Comme Dkk est en 1D, l’élément (i, kk) est à l’index i * b + kk.
        int dik = Dkk[i * b + kk];

        // Si la distance i → k est infinie, il n’existe aucun chemin i→k,
        // donc il est inutile d’essayer de passer par kk pour aller vers j.
        if (dik == INF) continue;

        for (int j = 0; j < bs; ++j) {
            // On récupère D[k][j] local dans le bloc diagonal.
            // L’élément (kk, j) est à l’index kk * b + j.
            int kkj = Dkk[kk * b + j];

            // Si la distance k → j est infinie, il n’y a pas de chemin k→j,
            // donc i→k→j est impossible et ne peut pas améliorer D[i][j].
            if (kkj == INF) continue;

            // via = distance i → j en passant par k :
            // via = D[i][k] + D[k][j]
            int via = dik + kkj;

            // Référence vers D[i][j] dans le bloc diagonal.
            // Toute modification de dij met directement à jour Dkk.
            int& dij = Dkk[i * b + j];

            // Si le chemin i → k → j est plus court que la valeur actuelle
            // stockée dans D[i][j], on met à jour D[i][j].
            if (via < dij) {
                dij = via;
            }
        }
    }
}

}

// fw_row : met à jour un bloc de ligne DkJ (bs x wJ) en utilisant le bloc pivot Dkk (bs x bs)
// ----------------------------------------------------------------------------------------
// Contexte : Floyd–Warshall par blocs
// On se place sur la ligne de blocs k (les blocs (k, j) pour tout j).
// - Dkk : bloc diagonal (k, k) déjà mis à jour par fw_block, taille logique bs x bs
// - DkJ : bloc de ligne (k, jBloc), taille logique bs x wJ
// - b   : leading dimension (nombre de colonnes allouées, sert pour l’indexation 1D)
//
// L’objectif est d’appliquer la formule de Floyd–Warshall à l’intérieur du bloc de ligne :
//    D[i][j] = min( D[i][j], D[i][k] + D[k][j] )
// mais restreinte aux indices appartenant à ces deux blocs :
//    - i et k sont dans le bloc pivot Dkk
//    - j est dans le bloc de ligne DkJ
//
// Dkk est indexé en 1D par : Dkk[i * b + kk]  (i, kk)
// DkJ est indexé en 1D par : DkJ[kk * b + j]  (kk, j) et DkJ[i * b + j] (i, j)
//
// Paramètres :
//   - Dkk : pointeur constant vers le bloc pivot (bs x bs), déjà optimisé par fw_block
//   - DkJ : pointeur vers le bloc de ligne à mettre à jour (bs x wJ)
//   - bs  : taille logique du bloc le long de i et kk (peut être < b en bordure)
//   - wJ  : largeur logique du bloc DkJ (peut être < b en bordure droite)
//   - b   : leading dimension utilisée pour l’indexation (taille d’une ligne allouée)
//
// Déroulement :
//   1) On parcourt toutes les lignes i du bloc (0 .. bs-1).
//   2) Pour chaque colonne pivot locale kk (0 .. bs-1), on lit D[i][k] dans Dkk.
//      - Si D[i][k] == INF : il n’existe pas de chemin i→k, inutile d’essayer i→k→j.
//   3) Pour chaque colonne j du bloc DkJ (0 .. wJ-1), on lit D[k][j] dans DkJ.
//      - Si D[k][j] == INF : il n’existe pas de chemin k→j, inutile de tester i→k→j.
//   4) Sinon, on calcule via = D[i][k] + D[k][j].
//      - Si via < D[i][j], on met à jour D[i][j] = via.
//   5) Au final, tout le bloc de ligne DkJ est optimisé en tenant compte des chemins passant
//      par tout k local au bloc pivot Dkk.
//
static void fw_row(const int* Dkk, int* DkJ, int bs, int wJ, int b) {
    for (int i = 0; i < bs; ++i) {                  // pour chaque ligne locale i du bloc
        for (int kk = 0; kk < bs; ++kk) {          // pour chaque pivot local kk dans Dkk
            int dik = Dkk[i * b + kk];             // D[i][k] dans le bloc pivot
            if (dik == INF) continue;              // pas de chemin i->k, on saute

            for (int j = 0; j < wJ; ++j) {         // pour chaque colonne j du bloc de ligne
                int kkj = DkJ[kk * b + j];         // D[k][j] dans le bloc de ligne
                if (kkj == INF) continue;          // pas de chemin k->j, on saute

                int via = dik + kkj;               // coût du chemin i->k->j
                int& dij = DkJ[i * b + j];         // référence vers D[i][j] dans DkJ
                if (via < dij) {                   // si i->k->j est meilleur
                    dij = via;                     // on met à jour D[i][j]
                }
            }
        }
    }
}


// fw_col met à jour un **bloc de colonne** Dik à l’aide du bloc pivot diagonal Dkk.
// Contexte : on applique la formule de Floyd–Warshall D[i][j] = min(D[i][j], D[i][k] + D[k][j])
// mais à l’échelle d’un bloc (sous-matrice).
//
// Paramètres :
//   - Dik : bloc de colonne (i,k) de taille logique hI × bs, stocké dans un tableau de taille b*b.
//           Il contient des distances D[i][j] pour un certain ensemble de lignes i
//           et pour les colonnes du bloc k. Ce bloc est **modifié** par la fonction.
//   - Dkk : bloc pivot diagonal (k,k) de taille logique bs × bs, stocké dans un tableau de taille b*b.
//           Il contient des distances D[k][j] pour k et j dans le même bloc pivot.
//   - hI  : nombre de lignes utiles dans le bloc Dik (peut être < b si on est sur le bord bas de la matrice).
//   - bs  : taille logique du bloc pivot (nombre de lignes/colonnes utiles dans Dkk et colonnes utiles dans Dik).
//           Peut être < b si on est sur le bord droit ou bas de la matrice.
//   - b   : leading dimension, c’est-à-dire le nombre de colonnes **allouées** pour chaque bloc en mémoire.
//           Pour accéder à l’élément (ligne L, colonne C) dans un bloc, on utilise l’index : L * b + C.
//
// Idée générale :
//   On parcourt toutes les lignes i du bloc Dik, puis tous les "pivots internes" kk du bloc k,
//   et enfin toutes les colonnes j du bloc k. À chaque triple (i, kk, j),
//   on tente d’améliorer la distance D[i][j] via le chemin i → kk → j.
//   Concrètement :
//      - ik = Dik[i][kk] = distance i → kk
//      - kj = Dkk[kk][j] = distance kk → j
//      - via = ik + kj   = distance i → kk → j
//      - si via < Dik[i][j], on met à jour Dik[i][j] = via.
//
// Détail de la triple boucle :
//   for (int i = 0; i < hI; ++i)
//     -> on itère sur toutes les lignes locales du bloc Dik (tous les sommets i de ce bloc).
//
//     for (int kk = 0; kk < bs; ++kk)
//       -> kk joue le rôle du sommet intermédiaire "k" dans la formule de Floyd–Warshall,
//          mais restreint à l’intérieur du bloc pivot (0..bs-1).
//       -> on lit ik = Dik[i * b + kk], la distance i → kk dans ce bloc de colonne.
//       -> si ik == INF, aucun chemin i → kk n’existe, donc i → kk → j ne pourra jamais améliorer i → j.
//
//       for (int j = 0; j < bs; ++j)
//         -> on balaye toutes les colonnes j du bloc pivot (0..bs-1),
//            qui correspondent aux destinations j dans le même bloc k.
//         -> on lit kj = Dkk[kk * b + j], la distance kk → j dans le bloc pivot.
//         -> si kj == INF, aucun chemin kk → j n’existe, donc i → kk → j n’existe pas.
//
//         -> on calcule via = ik + kj, c’est le coût du chemin i → kk → j.
//         -> on récupère une référence sur ij = Dik[i * b + j], la distance actuelle i → j dans le bloc Dik.
//         -> si via < ij, on met à jour ij = via.
//
// Au final, après l’exécution de fw_col, toutes les distances du bloc de colonne Dik
// ont été éventuellement raccourcies en utilisant le bloc pivot Dkk comme ensemble
// de sommets intermédiaires possibles.

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

// fw_inner met à jour un bloc "interne" Dij de la matrice des distances.
// Ce bloc Dij ne se trouve ni sur la ligne du pivot k, ni sur la colonne du pivot k.
// On utilise pour cela :
//   - Dik : bloc de colonne (i, k) de taille logique hI × bs
//   - DkJ : bloc de ligne   (k, j) de taille logique bs × wJ
//   - Dij : bloc interne    (i, j) de taille logique hI × wJ
// Tous ces blocs sont stockés dans des tableaux 1D avec un leading dimension = b.
//
// Pour chaque ligne locale i du bloc Dij (0 ≤ i < hI) :
//   on parcourt chaque "pivot interne" kk dans le bloc k (0 ≤ kk < bs) et on lit :
//     ik = Dik[i, kk] : distance i → k (dans le bloc de colonne)
//   si ik == INF, inutile de continuer sur ce kk (aucun chemin i → k).
//
//   Sinon, pour chaque colonne j du bloc Dij (0 ≤ j < wJ) :
//     kj = DkJ[kk, j] : distance k → j (dans le bloc de ligne)
//     si kj == INF, aucun chemin k → j via ce kk, on passe au suivant.
//
//     Sinon on calcule la distance via k : via = ik + kj,
//     et on compare avec la valeur actuelle ij = Dij[i, j].
//     Si via < ij, on met à jour Dij[i, j] avec cette nouvelle distance plus courte.
//
// Au final, fw_inner applique la relation de Floyd–Warshall D[i,j] = min(D[i,j], D[i,k] + D[k,j])
// mais limitée à un bloc interne Dij, en utilisant les blocs Dik (colonne) et DkJ (ligne).
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

       /*
     * Récupération des informations MPI de base :
     * - rank : identifiant du processus courant dans MPI_COMM_WORLD (0..size-1)
     * - size : nombre total de processus lancés.
     */
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // -------------------------------------------------------
    // 0) Construction d’une "grille" 2D de processus + choix de la taille de bloc
    // -------------------------------------------------------

    /*
     * Je demande à MPI de décomposer automatiquement mes 'size' processus
     * en une grille 2D de dimensions (Pr x Pc), aussi équilibrée que possible.
     * Exemple : size = 6  →  Pr=2, Pc=3 (2x3 = 6).
     */
    int dims[2] = {0, 0};
    MPI_Dims_create(size, 2, dims);  // calcule Pr, Pc "équilibrés"
    int Pr = dims[0];                // nombre de processus en "ligne"
    int Pc = dims[1];                // nombre de processus en "colonne"

    /*
     * Choix d’une taille de bloc b :
     * - approxBlocks ≈ sqrt(size) = nombre "souhaité" de blocs par dimension.
     * - b ≈ n / approxBlocks, avec un ceil fait "à la main".
     * L’idée : avoir un nombre de blocs cohérent avec le nombre de processus,
     * pour une répartition raisonnablement équilibrée.
     */
    int approxBlocks = (int)std::sqrt((double)size);
    if (approxBlocks < 1) approxBlocks = 1;
    int b = (n + approxBlocks - 1) / approxBlocks;  // taille théorique du bloc (ceil(n/approxBlocks))

    if (b < 1) b = 1;  // sécurité : ne jamais avoir de bloc de taille 0

    /*
     * nb = nombre de blocs par dimension (en ligne et en colonne).
     * C’est un ceil(n / b) : même si n n’est pas multiple de b, on couvre toute la matrice.
     */
    int nb = (n + b - 1) / b;  // ceil(n/b)

    // -------------------------------------------------------
    // 1) Déterminer quels blocs appartiennent au processus courant (rank)
    // -------------------------------------------------------

    /*
     * computeLocalBlocks(...) parcourt tous les blocs (bi, bj) de la matrice globale
     * et applique la stratégie de distribution (ownerOf) sur la grille Pr x Pc.
     * Il renvoie, pour ce rank, la liste des BlockInfo qu’il possède :
     *  - bi, bj      : indices de bloc (ligne/colonne de bloc)
     *  - owner       : rang MPI propriétaire
     *  - offset_i/j  : coordonnées (ligne/colonne) du coin haut-gauche dans la matrice globale.
     */
    vector<BlockInfo> localBlocks = computeLocalBlocks(n, b, Pr, Pc, rank);
    int numLocal = (int)localBlocks.size();  // nombre de blocs réellement gérés par ce rang
    int blockArea = b * b;                   // nombre de cases dans un bloc (taille max)

    /*
     * Je stocke tous les blocs locaux dans un seul grand tableau 1D :
     *  - bloc k → &localData[k * blockArea]
     * C’est plus simple à gérer et à envoyer en MPI que plein de allocations séparées.
     */
    vector<int> localData(numLocal * blockArea);

    /*
     * Table de correspondance (bi, bj) → index local dans localData.
     * - taille nb x nb (tous les blocs possibles),
     * - initialisée à -1 : valeur "bloc non présent sur ce rang".
     *
     * Ensuite, pour chaque bloc que je possède (localBlocks[idx]),
     * je marque localIndex[bi * nb + bj] = idx.
     *
     * Plus tard, si je veux savoir si je possède le bloc (bi, bj) :
     *  - int idx = localIndex[bi * nb + bj];
     *  - idx == -1  → je n’ai pas ce bloc,
     *  - idx >= 0   → je peux accéder à ce bloc via &localData[idx * blockArea].
     */
    vector<int> localIndex(nb * nb, -1);
    for (int idx = 0; idx < numLocal; ++idx) {
        const BlockInfo& info = localBlocks[idx];
        localIndex[info.bi * nb + info.bj] = idx;
    }


       // -------------------------------------------------------
    // 2) Initialisation de D^(-1) dans les blocs locaux
    // -------------------------------------------------------
    // Pour chaque bloc que ce rang MPI possède :
    //  - on récupère sa position dans la matrice globale (offset_i, offset_j)
    //  - on construit la matrice de distances initiale D^(-1) à l'intérieur de ce bloc
    //    à partir de la matrice d'adjacence globale `mat`.
    //
    // Rappel de la règle d'initialisation de Floyd–Warshall :
    //   - D[i][i] = 0                                 (distance d’un sommet à lui-même)
    //   - D[i][j] = poids de l’arête (i,j)            s’il existe une arête directe
    //   - D[i][j] = INF                              s’il n’y a pas d’arête directe
    //
    // Ici on applique cette règle mais bloc par bloc :
    //   * `localBlocks[idx]` décrit le bloc numéro idx :
    //       - info.offset_i = ligne globale du coin haut-gauche du bloc
    //       - info.offset_j = colonne globale du coin haut-gauche du bloc
    //   * `blk` pointe vers la zone mémoire du bloc dans `localData`
    //       (on a numLocal blocs à la suite, chacun de taille b*b).
    //
    // Pour chaque case locale (ii, jj) du bloc :
    //   - on calcule l’indice global correspondant (gi, gj)
    //   - si (gi, gj) est en dehors de la matrice n×n (bord quand n n’est pas multiple de b),
    //     on met INF
    //   - sinon, si gi == gj          → D[gi][gj] = 0
    //   - sinon, si mat[gi*n + gj] = 0 (et gi != gj) → pas d’arête, donc INF
    //   - sinon                        → on met le poids mat[gi*n + gj]
    //
    // Au final, après cette boucle, `localData` contient la partie locale de D^(-1),
    // organisée en blocs, prête à être mise à jour par les phases de l’algorithme
    // de Floyd–Warshall parallélisé.
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
        /*
     * Étape 3 : boucle principale de Floyd–Warshall par blocs
     * -------------------------------------------------------
     * On parcourt kk = 0..nb-1, où kk est l’indice du bloc pivot sur la diagonale.
     *
     * Pour chaque kk, on réalise trois phases :
     *
     *   A) Bloc pivot D[kk, kk]
     *      - On détermine d’abord quel processus possède le bloc diagonal (kk, kk)
     *        via ownerOf(kk, kk, Pr, Pc) et localIndex.
     *      - Ce processus applique fw_block(...) sur ce bloc : c’est un Floyd–Warshall
     *        complet mais limité à l’intérieur de ce bloc b×b.
     *      - Le bloc pivot mis à jour est copié dans pivotBlock, puis diffusé à
     *        tous les processus avec MPI_Bcast. Ainsi, tout le monde connaît D[kk, kk].
     *      - On stocke aussi ce bloc pivot dans rowBlocks[kk] et colBlocks[kk], car
     *        il appartient à la fois à la ligne de blocs kk et à la colonne de blocs kk.
     *
     *   B1) Blocs de ligne (kk, j) : mise à jour + diffusion
     *      - Pour chaque bloc de la même ligne de blocs que le pivot (kk, jb),
     *        sauf le pivot lui-même (jb != kk) :
     *          * ownerRow = ownerOf(kk, jb, Pr, Pc) indique quel rang possède ce bloc.
     *          * Sur ce rang, on récupère le bloc local DkJ = (kk, jb) dans localData,
     *            on l’actualise avec fw_row(pivotBlock, DkJ, ...) en utilisant le pivot.
     *          * Le bloc de ligne mis à jour est copié dans tempBlock.
     *          * tempBlock est diffusé à tous les processus via MPI_Bcast.
     *      - Après diffusion, rowBlocks[jb] contient le bloc (kk, jb) mis à jour
     *        sur tous les rangs. On a donc la ligne de blocs kk complète en mémoire locale.
     *
     *   B2) Blocs de colonne (i, kk) : mise à jour + diffusion
     *      - De façon symétrique, pour chaque bloc de la même colonne de blocs
     *        que le pivot (ib, kk), avec ib != kk :
     *          * ownerCol = ownerOf(ib, kk, Pr, Pc) donne le rang propriétaire.
     *          * Sur ce rang, on récupère Dik = bloc (ib, kk), on l’actualise avec
     *            fw_col(Dik, pivotBlock, ...) en utilisant le pivot.
     *          * On copie le bloc mis à jour dans tempBlock.
     *          * tempBlock est diffusé à tous via MPI_Bcast.
     *      - Après diffusion, colBlocks[ib] contient le bloc (ib, kk) mis à jour
     *        sur tous les rangs. On a donc la colonne de blocs kk complète en mémoire locale.
     *
     *   C) Mise à jour des blocs internes (i, j), avec i != kk et j != kk
     *      - Pour chaque bloc local (ib, jb) détenu par ce processus :
     *          * Si ib == kk ou jb == kk, on saute : ces blocs de ligne/colonne
     *            ont déjà été traités en B1/B2.
     *          * Sinon, c’est un bloc "interne" (ib, jb) :
     *              - hI = hauteur réelle du bloc (cas bord éventuel),
     *                wJ = largeur réelle du bloc.
     *              - Dij = pointeur vers le bloc (ib, jb) dans localData.
     *              - Dik = bloc (ib, kk) correspondant dans colBlocks[ib].
     *              - DkJ = bloc (kk, jb) correspondant dans rowBlocks[jb].
     *              - On appelle fw_inner(Dik, DkJ, Dij, ...) qui applique la formule
     *                de Floyd–Warshall au niveau bloc :
     *                   Dij = min(Dij, Dik + DkJ)
     *                pour tous les couples de sommets dans ces blocs.
     *      - À la fin de cette phase C, tous les blocs locaux de ce rang sont
     *        à jour pour ce pivot kk.
     *
     * En résumé :
     *   - On traite le bloc diagonal (pivot),
     *   - on propage les mises à jour sur la ligne et la colonne du pivot,
     *   - puis on met à jour tous les autres blocs (i, j) en utilisant
     *     les blocs de colonne (i, kk) et de ligne (kk, j).
     * C’est exactement la logique de Floyd–Warshall, mais appliquée
     * à une matrice découpée en blocs et distribuée entre plusieurs processus MPI.
     */

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
    // Objectif : reconstruire la matrice globale des distances D (n×n)
    //            à partir des blocs distribués sur tous les processus MPI,
    //            et ne conserver cette matrice complète que sur le rang 0.
    //
    // Étape 1 : allocation de D_final uniquement sur le rang 0
    //  - Tous les rangs déclarent le pointeur D_final, mais seul rank 0
    //    alloue un tableau de taille n×n.
    //  - On initialise toutes les cases à INF pour éviter les valeurs
    //    indéfinies dans les zones qui ne seront jamais écrasées
    //    (par exemple sur les bords quand n n’est pas multiple de b).
    //
    // Étape 2 : parcours de tous les rangs r = 0..size-1
    //  - computeLocalBlocks(n, b, Pr, Pc, r) est déterministe : le rang 0
    //    peut donc recalculer quels blocs appartiennent à chaque rang r
    //    sans communication supplémentaire.
    //  - Pour chaque bloc (bi, bj) appartenant au rang r :
    //       * Sur le rang r :
    //           - on retrouve l’index local du bloc via localIndex[bi * nb + bj]
    //           - on copie les données du bloc depuis localData[...] dans
    //             un buffer temporaire buf (taille b×b).
    //       * Transfert vers le rang 0 :
    //           - si r == 0 : pas de MPI, le rang 0 lit directement ses blocs
    //             depuis buf (déjà rempli localement).
    //           - si r != 0 : le rang r envoie son buf avec MPI_Send,
    //             et le rang 0 le reçoit avec MPI_Recv.
    //       * Intégration dans la matrice globale (sur rank 0 uniquement) :
    //           - on utilise info.offset_i / info.offset_j pour retrouver
    //             les coordonnées globales (gi, gj) dans D_final.
    //           - on recopie chaque case de buf[ii * b + jj] dans
    //             D_final[gi * n + gj], en ignorant les indices (gi, gj)
    //             qui sortent de la matrice (tests gi >= n et gj >= n),
    //             ce qui gère les blocs partiellement en dehors de la matrice.
    //
    // Étape 3 : valeur de retour
    //  - Sur le rang 0 : on renvoie le pointeur D_final (matrice complète).
    //  - Sur les autres rangs : on renvoie nullptr, car ils n’ont pas besoin
    //    de stocker la matrice globale des distances.
    //
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
