#define OMPI_SKIP_MPICXX 1
#include "Distribution.hpp"

using namespace std;

/**
 * @file Distribution.cpp
 * @brief Implémentation des fonctions de distribution des blocs.
 */

/**
 * @brief Renvoie le rang MPI propriétaire du bloc (bi, bj).
 *
 * @param bi Indice de ligne du bloc.
 * @param bj Indice de colonne du bloc.
 * @param Pr Nombre de processus sur l'axe des lignes.
 * @param Pc Nombre de processus sur l'axe des colonnes.
 * @return Rang MPI propriétaire.
 */
int ownerOf(int bi, int bj, int Pr, int Pc) {
    int pr = bi % Pr;
    int pc = bj % Pc;
    return pr * Pc + pc;
}

/**
 * @brief Construit la liste des blocs locaux appartenant au processus courant.
 *
 * @param nb_nodes Taille de la matrice globale (n × n).
 * @param b        Taille d'un bloc.
 * @param Pr       Nombre de processus sur l'axe des lignes.
 * @param Pc       Nombre de processus sur l'axe des colonnes.
 * @param rank     Rang MPI du processus courant.
 * @return Vecteur de BlockInfo pour les blocs locaux.
 */
vector<BlockInfo> computeLocalBlocks(int nb_nodes, int b, int Pr, int Pc, int rank) {
    int nb = (nb_nodes + b - 1) / b; // ceil(n/b)
    vector<BlockInfo> list;

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
     
