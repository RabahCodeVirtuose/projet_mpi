#include "Distribution.hpp"

int ownerOf(int bi, int bj, int Pr, int Pc) {
    int pr = bi % Pr;
    int pc = bj % Pc;
    return pr * Pc + pc;
}

std::vector<BlockInfo> computeLocalBlocks(int nb_nodes, int b, int Pr, int Pc, int rank) {

    int nb = (nb_nodes + b - 1) / b; // ceil(n/b)
    std::vector<BlockInfo> list;

    for (int bi = 0; bi < nb; bi++) {
        for (int bj = 0; bj < nb; bj++) {

            int own = ownerOf(bi, bj, Pr, Pc);
            if (own == rank) {
                BlockInfo info;
                info.bi = bi;
                info.bj = bj;
                info.owner = own;
                info.offset_i = bi * b;
                info.offset_j = bj * b;
                list.push_back(info);
            }
        }
    }
    return list;
}
