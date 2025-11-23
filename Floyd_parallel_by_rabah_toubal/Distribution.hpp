#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include <vector>

struct BlockInfo {
    int bi, bj;      // Coordonnées bloc
    int owner;       // Rang du propriétaire MPI
    int offset_i;    // Position dans la matrice globale
    int offset_j;
};

int ownerOf(int bi, int bj, int Pr, int Pc);
std::vector<BlockInfo> computeLocalBlocks(int nb_nodes, int b, int Pr, int Pc, int rank);

#endif
