#ifndef PARALLEL_FW_BLOCKS_HPP
#define PARALLEL_FW_BLOCKS_HPP

/**
 * @file ParallelFWBlocks.hpp
 * @brief Déclaration de l'algorithme parallèle de Floyd–Warshall par blocs (MPI).
 *
 * Cette fonction constitue le cœur du calcul parallèle : la matrice des distances
 * est découpée en blocs de taille B×B, distribués sur une grille 2D de processus,
 * puis mise à jour itérativement selon le schéma bloc-cyclique :
 *
 *   - Phase A : calcul du bloc pivot Dkk (diagonale)
 *   - Phase B : mise à jour et diffusion des blocs de la ligne k (DkJ)
 *               et de la colonne k (DIk)
 *   - Phase C : mise à jour des blocs internes (DIJ)
 *
 * Les communications font appel à MPI_Bcast et MPI_Ibcast afin de recouvrir
 * calculs et communications lorsque cela est possible. La matrice finale est
 * rassemblée sur le processus 0.
 */

/**
 * @brief Algorithme parallèle de Floyd–Warshall utilisant une distribution en blocs.
 *
 * Cette fonction applique l'algorithme complet de Floyd–Warshall à une matrice
 * d'adjacence initiale en la distribuant sur plusieurs processus MPI selon une
 * grille 2D bloc-cyclique. Chaque processus stocke uniquement certains blocs
 * de taille B×B de la matrice globale, les met à jour localement et échange
 * les blocs nécessaires lors des phases pivot/ligne/colonne.
 *
 * @param n   Taille de la matrice initiale (n × n).
 * @param mat Matrice d'adjacence initiale stockée à plat (row-major) :
 *            l'élément (i, j) est à l'indice i * n + j.
 *            Cette matrice doit être identique sur tous les processus au moment
 *            de l'appel (diffusée préalablement).
 *
 * @return Sur le rang 0 : un pointeur vers la matrice finale des distances
 *         (n × n), allouée avec new[] et devant être libérée par l'appelant.
 *         Sur les autres rangs : nullptr.
 *
 * @note La fonction doit être appelée après MPI_Init et avant MPI_Finalize.
 */
int* ParallelFloydWarshallBlocks(int n, int* mat);

#endif // PARALLEL_FW_BLOCKS_HPP
