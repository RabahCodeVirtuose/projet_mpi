#define OMPI_SKIP_MPICXX 1

/**
 * @file PAM.cpp
 * @brief Implémentation de l'algorithme PAM en version parallèle MPI.
 *
 * Ce fichier contient :
 *  - des fonctions utilitaires internes pour calculer le coût d'un ensemble
 *    de médioïdes (en local ou en distribué),
 *  - l'implémentation de la fonction publique runPAM_MPI().
 *
 * Le coût global est calculé en parallèle en répartissant les sommets entre
 * les processus. Les échanges médoïde / non-médoïde sont pilotés par le
 * rang 0, qui choisit les meilleurs médioïdes et les diffuse aux autres.
 */

#include "PAM.hpp"
#include <mpi.h>
#include <limits>
#include <random>
#include <algorithm>
#include <ctime>   

static const int INF = 100000;

/**
 * @brief Calcule la somme des distances au médioïde le plus proche,
 *        et remplit clusterOf & distToMedoid.
 *
 * Cette fonction est utilisée uniquement sur le rang 0 à la fin
 * pour reconstruire la partition détaillée :
 * pour chaque sommet i, on mémorise :
 *   - l'indice du cluster (indice du médioïde dans medoids),
 *   - la distance au médioïde le plus proche.
 *
 * @param dist         Matrice de distances (n × n), stockée à plat.
 * @param n            Nombre de sommets.
 * @param medoids      Liste des indices de médioïdes.
 * @param clusterOf    Vecteur de sortie : clusterOf[i] = indice du médioïde associé.
 * @param distToMedoid Vecteur de sortie : distToMedoid[i] = distance au médioïde.
 *
 * @return La somme totale des distances de chaque sommet à son médioïde.
 */
static long long computeCostAndAssign(const std::vector<int>& dist,
                                      int n,
                                      const std::vector<int>& medoids,
                                      std::vector<int>& clusterOf,
                                      std::vector<int>& distToMedoid)
{
    int k = (int)medoids.size();
    long long totalCost = 0;
    // Pour chaque sommet i, je cherche le medioide le plus proche.

    for (int i = 0; i < n; ++i) {
        int bestMedoidIdx = 0;   // indice du meilleur medioide dans [0..k-1]
        int bestDist = INF;      // distance minimale pour ce sommet

                // Je teste tous les medioides possibles pour ce sommet.
        for (int m = 0; m < k; ++m) {
            int med = medoids[m];        // indice global du medioide
            int d   = dist[i * n + med]; // distance de i à ce medioide
            if (d < bestDist) {
                bestDist = d;
                bestMedoidIdx = m;
            }
        }
  // A la fin de la boucle sur m, j'ai trouvé le medioide le plus proche.
        // Je mémorise :
        //  - dans quel cluster i tombe (bestMedoidIdx),
        //  - la distance de i à son medioide.
        clusterOf[i]    = bestMedoidIdx;
        distToMedoid[i] = bestDist;
        totalCost      += bestDist;
    }
    // Je renvoie la somme des distances de tous les sommets
    // à leur medioide respectif.

    return totalCost;
}

/**
 * @brief Calcule le coût global pour un ensemble de médioïdes, en parallèle MPI.
 *
 * Chaque processus :
 *   - traite l’intervalle de sommets [start, end),
 *   - pour chaque sommet i il cherche le médioïde le plus proche,
 *   - accumule sa contribution dans localCost.
 *
 * Un MPI_Allreduce(MPI_SUM) donne le coût total global.
 */
static long long computeCostDistributed(const std::vector<int>& dist,
                                        int n,
                                        const std::vector<int>& medoids)
{
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int k = (int)medoids.size();

    // Découpage simple des sommets entre les processus :
    // chaque rang traite un intervalle [start, end[.
    // chunk = nombre de sommets par processus (en arrondissant vers le haut)

    // Découpage simple des sommets entre les processus
    int chunk = (n + size - 1) / size;  // ceil(n / size)
    int start = rank * chunk;
    int end   = std::min(n, start + chunk);

    long long localCost = 0;
    // Ici, chaque rang parcourt seulement ses sommets à lui

    for (int i = start; i < end; ++i) {
        int bestDist = INF;

        // Comme dans PAM classique : pour ce sommet i,
        // je cherche la distance au medioide le plus proche
        for (int m = 0; m < k; ++m) {
            int med = medoids[m];
            int d   = dist[i * n + med];
            if (d < bestDist) {
                bestDist = d;
            }
        }
        // J'ajoute la meilleure distance pour ce sommet à mon coût local.

        localCost += bestDist;
    }

     // Réduction pour obtenir le coût global
    // Ici chaque rang envoie son localCost, et MPI_Allreduce fait la somme.
    // Au final, globalCost contient la somme de tous les localCost,
    // donc le coût total comme si j'avais tout fait sur un seul processus.
    long long globalCost = 0;
    MPI_Allreduce(&localCost, &globalCost, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    return globalCost;
}

/* ===================================================================== */
/*  Version MPI de PAM                                                   */
/* ===================================================================== */
// Cette fonction c'est ma version parallèle de l'algo PAM.
// L'idée c'est :
//   - j'ai une matrice de distances dist (taille n*n, identique sur tous les rangs),
//   - je choisis k medioides,
//   - je calcule le coût de la partition en parallèle avec MPI,
//   - et j'améliore les medioides en testant des échanges.
//
// Le calcul de coût est distribué (computeCostDistributed),
// mais la logique "qui est le meilleur échange ?" est centralisée sur le rang 0.
//
PAMResult runPAM_MPI(const std::vector<int>& dist, int n, int k) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    PAMResult res;
    res.medoids.resize(k);

    // -----------------------------
    // 1) Initialisation des médioïdes (rang 0)
    // -----------------------------
        // Pour commencer, je laisse le rang 0 choisir des medioides initiaux au hasard.
    if (rank == 0) {
        std::vector<int> allIndices(n);
            for (int i = 0; i < n; ++i) allIndices[i] = i;

        // Je seed mon générateur aléatoire avec le temps.
            unsigned seed = (unsigned) time(NULL);
// mt19937 en gros c’est un générateur aléatoire moderne, plus stable que le vieux rand()   
// Je l’utilise parce que quand je mélange les indices pour choisir mes medioides,
// j’ai besoin d’un truc un minimum fiable sinon le tirage sera éclaté.
            std::mt19937 gen(seed);

        // Je mélange la liste des indices
            std::shuffle(allIndices.begin(), allIndices.end(), gen);

        // Et je prends les k premiers comme medioides initiaux
            for (int m = 0; m < k; ++m) {
                res.medoids[m] = allIndices[m];
            }

    }

    // Diffuser les médioïdes initiaux à tous
    MPI_Bcast(res.medoids.data(), k, MPI_INT, 0, MPI_COMM_WORLD);

    // Je calcule le coût initial de ces medioides avec la version distribuée.
    long long bestCost = computeCostDistributed(dist, n, res.medoids);

    // -----------------------------
    // 2) Boucles d'amélioration : on teste tous les échanges (m, h)
   //    On teste tous les échanges (m, h) :
    //    remplacer le medioide m par un sommet h qui n'est pas encore medioide.
    //    Tous les rangs calculent le coût pour chaque candidat,
    //    mais seul le rang 0 garde le meilleur et décide.
    // -----------------------------
    while (true) {
        bool improvedLocal= false;  // seulement sur rang 0
        long long bestCostThisPass = bestCost;
        std::vector<int> bestMedoidsThisPass = res.medoids;

        // Je parcours tous les medioides m possibles
        for (int m = 0; m < k; ++m) {
             // Et pour chaque sommet h du graphe,
            // je regarde ce que ça donnerait si je le mettais medioide à la place.
            for (int h = 0; h < n; ++h) {
                // D'abord je vérifie si h est déjà un medioide, sinon ça ne sert à rien.
                bool isMedoid = false;
                for (int mm = 0; mm < k; ++mm) {
                    if (res.medoids[mm] == h) {
                        isMedoid = true;
                        break;
                    }
                }
                if (isMedoid) continue;

                 // Je construis un ensemble de medioides candidat :
                // c'est la même liste, sauf que je remplace le medioide m par h
                std::vector<int> newMedoids = res.medoids;
                newMedoids[m] = h;

               // Je demande le coût global de cette configuration newMedoids.
                // computeCostDistributed utilise tous les rangs pour calculer ce coût.
                long long newCost = computeCostDistributed(dist, n, newMedoids);
               
                // Seul le rang 0 décide si c'est mieux ou pas.
                if (rank == 0 && newCost < bestCostThisPass) {
                    bestCostThisPass = newCost;
                    bestMedoidsThisPass = newMedoids;
                    improvedLocal = true;
                }
            }
        }

        // Rang 0 dit aux autres si une amélioration a été trouvée
        int flag = improvedLocal ? 1 : 0;
                // Maintenant le rang 0 doit dire aux autres si on a trouvé une amélioration.
        MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!flag) {
            // Si flag == 0, ça veut dire qu'aucun échange (m,h) n'améliore le coût.
            // Donc on arrête la boucle : PAM est arrivé à un optimum local.
            break;
        }

        // Sinon, il y a eu une amélioration.
        // Le rang 0 met à jour les medioides et le coût courant.
        if (rank == 0) {
            res.medoids = bestMedoidsThisPass;
            bestCost    = bestCostThisPass;
        }

        // Puis on rediffuse les nouveaux medioides à tous les rangs,
        // et on repart pour un tour de boucle while.
        MPI_Bcast(res.medoids.data(), k, MPI_INT, 0, MPI_COMM_WORLD);
    }

  // -----------------------------
    // 3) Affectation finale (rang 0 uniquement)
    // -----------------------------
    // Ici, je veux reconstruire les clusters finaux proprement.
    // Je le fais juste sur le rang 0, parce que c'est lui qui va écrire le résultat.
    if (rank == 0) {
        res.clusterOf.resize(n);
        res.distToMedoid.resize(n);

        long long finalCost = computeCostAndAssign(dist, n,
                                                   res.medoids,
                                                   res.clusterOf,
                                                   res.distToMedoid);
        res.totalCost = finalCost;
    } 

    return res;
}
