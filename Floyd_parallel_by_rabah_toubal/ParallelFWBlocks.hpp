#ifndef PARALLEL_FW_BLOCKS_HPP
#define PARALLEL_FW_BLOCKS_HPP

// Version Floyd–Warshall parallèle avec décomposition en blocs B×B
// Retour : sur rank 0, pointeur vers matrice des distances n×n
//          sur les autres rangs : nullptr
int* ParallelFloydWarshallBlocks(int n, int* mat);

#endif
