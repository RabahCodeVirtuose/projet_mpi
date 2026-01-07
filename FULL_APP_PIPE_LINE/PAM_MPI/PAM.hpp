#ifndef PAM_HPP
#define PAM_HPP

#include <vector>

using namespace std;

/**
 * @file PAM.hpp
 * @brief Définition des structures et de l'interface de l'algorithme PAM en version MPI+OpenMP.
 *
 * Ce module définit :
 *  - la structure PAMResult qui stocke le résultat de l'algorithme,
 *  - la fonction runPAM_MPI qui exécute l'algorithme PAM en parallèle (MPI + OpenMP).
 */

/**
 * @brief Résultat de l'algorithme PAM.
 *
 * - medoids : indices des k médioïdes choisis (entre 0 et n-1)
 * - clusterOf[i] : pour chaque sommet i, indice du cluster (0..k-1)
 *                  auquel i est affecté (c'est l'indice du méd(oïde) dans le tableau medoids)
 * - distToMedoid[i] : distance de i à son médioïde le plus proche
 * - totalCost : somme des distances de tous les sommets à leur médioïde
 */
struct PAMResult {
    vector<int> medoids;       /**< Indices des médioïdes choisis. */
    vector<int> clusterOf;     /**< Pour chaque sommet i, indice du cluster (0..k-1). */
    vector<int> distToMedoid;  /**< Distance entre chaque sommet et son médioïde. */
    long long totalCost = 0;        /**< Coût total (somme des distances au médioïde). */
};

/**
 * @brief Version MPI+OpenMP de l'algorithme PAM (Partitioning Around Medoids).
 *
 * Le principe :
 *  - La matrice de distances est répliquée sur tous les processus.
 *  - Le coût pour un ensemble de médioïdes est calculé en parallèle :
 *    chaque processus traite un sous-ensemble de sommets (MPI)
 *    et parallélise localement la boucle sur les sommets (OpenMP).
 *  - Les décisions d’amélioration (échanges médoïde / non-médoïde)
 *    sont prises par le rang 0 et diffusées à tous.
 *
 * @param dist Matrice des distances de taille n*n, stockée à plat :
 *             dist[i * n + j] contient la distance entre i et j.
 * @param n    Nombre total de sommets.
 * @param k    Nombre de groupes / médioïdes.
 *
 * @return Sur le rang 0 : résultat complet (médioïdes, clusters, coût).
 *         Sur les autres rangs : seul totalCost est rempli, le reste n’est pas utilisé.
 */
PAMResult runPAM_MPI(const vector<int>& dist, int n, int k);

#endif 
