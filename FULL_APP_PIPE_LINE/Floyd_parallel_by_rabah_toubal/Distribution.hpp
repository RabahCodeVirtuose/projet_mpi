#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include <vector>

/**
 * @file Distribution.hpp
 * @brief Structures et fonctions permettant de distribuer les blocs de la matrice
 *        dans une grille de processus MPI selon une distribution bloc-cyclique 2D.
 *
 * Ce module gère :
 *  - la description des blocs d'une matrice globale (structure BlockInfo),
 *  - l'affectation d'un bloc (bi, bj) à un processus MPI,
 *  - le calcul de tous les blocs locaux d'un processus donné.
 *
 * Il est utilisé dans l'implémentation parallèle de l'algorithme
 * de Floyd–Warshall par blocs.
 */

/**
 * @struct BlockInfo
 * @brief Décrit un bloc de taille b × b dans la matrice globale.
 *
 * Un bloc est identifié par :
 *  - ses coordonnées dans la grille de blocs (bi, bj),
 *  - le rang MPI du processus propriétaire (`owner`),
 *  - les offsets (offset_i, offset_j) représentant la position du coin
 *    supérieur gauche du bloc dans la matrice globale.
 *
 * Ces informations permettent à chaque processus de reconstruire ses blocs locaux
 * et de connaître leur emplacement exact dans la matrice distribuée.
 */
struct BlockInfo {
    int bi;       /**< Indice de ligne du bloc dans la grille des blocs. */
    int bj;       /**< Indice de colonne du bloc dans la grille des blocs. */
    int owner;    /**< Rang MPI du processus propriétaire du bloc. */
    int offset_i; /**< Indice global de la ligne correspondant au début du bloc. */
    int offset_j; /**< Indice global de la colonne correspondant au début du bloc. */
};

/**
 * @brief Détermine le rang MPI propriétaire d'un bloc (bi, bj).
 *
 * La distribution utilisée est une distribution bloc-cyclique 2D :
 * - la matrice globale est découpée en blocs de taille b × b,
 * - les processus sont organisés en une grille de dimensions Pr × Pc,
 * - le bloc (bi, bj) est affecté au processus :
 *       (bi % Pr, bj % Pc)
 *   converti en rang linéaire : pr * Pc + pc.
 *
 * @param bi  Indice de ligne du bloc dans la grille des blocs.
 * @param bj  Indice de colonne du bloc dans la grille des blocs.
 * @param Pr  Nombre de processus dans la dimension verticale (lignes).
 * @param Pc  Nombre de processus dans la dimension horizontale (colonnes).
 * @return Le rang MPI du propriétaire du bloc (bi, bj).
 */
int ownerOf(int bi, int bj, int Pr, int Pc);

/**
 * @brief Calcule la liste de tous les blocs locaux possédés par un processus MPI.
 *
 * Cette fonction :
 *  - parcourt l'ensemble des blocs de la matrice globale,
 *  - détermine le propriétaire de chacun via `ownerOf`,
 *  - construit un vecteur contenant uniquement les blocs appartenant au processus `rank`.
 *
 * Chaque bloc retourné inclut :
 *  - ses coordonnées (bi, bj),
 *  - le rang propriétaire,
 *  - les offsets dans la matrice globale.
 *
 * @param nb_nodes Taille de la matrice globale (n × n).
 * @param b        Taille d'un bloc (nombre de lignes/colonnes).
 * @param Pr       Nombre de processus dans la dimension des lignes.
 * @param Pc       Nombre de processus dans la dimension des colonnes.
 * @param rank     Rang MPI du processus courant.
 * @return Un vecteur contenant les BlockInfo correspondant aux blocs locaux du processus.
 */
std::vector<BlockInfo> computeLocalBlocks(int nb_nodes, int b, int Pr, int Pc, int rank);

#endif // DISTRIBUTION_HPP
