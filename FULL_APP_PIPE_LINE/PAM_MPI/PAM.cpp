#define OMPI_SKIP_MPICXX 1

/**
 * @file PAM.cpp
 * @brief Implémentation de l'algorithme PAM en version hybride MPI+OpenMP.
 *
 * Stratégie de parallélisation :
 *  - MPI : distribution des sommets entre les processus.
 *  - OpenMP : parallélisation du coût local (boucle sur i).
 *
 * Important : les appels MPI (Allreduce, Bcast) restent hors des régions OpenMP.
 */

#include "PAM.hpp"
#include <mpi.h>
#include <omp.h>
#include <limits>
#include <random>
#include <algorithm>
#include <ctime>
#include <iostream>

using namespace std;

static const int INF = 100000;

/**
 * @brief Calcule le coût et assigne les sommets aux médioïdes (version séquentielle).
 *
 * Utilisée uniquement par le rang 0 pour construire le résultat final.
 *
 * @param dist         Matrice des distances (n × n) stockée à plat.
 * @param n            Nombre de sommets.
 * @param medoids      Liste des indices de médioïdes.
 * @param clusterOf    Cluster associé à chaque sommet (sortie).
 * @param distToMedoid Distance au médioïde le plus proche (sortie).
 * @return Coût total (somme des distances aux médioïdes).
 */
static long long computeCostAndAssign(const vector<int>& dist,
                                      int n,
                                      const vector<int>& medoids,
                                      vector<int>& clusterOf,
                                      vector<int>& distToMedoid)
{
    int k = (int)medoids.size();
    long long totalCost = 0;

    for (int i = 0; i < n; ++i) {
        int bestMedoidIdx = 0;
        int bestDist = INF;

        for (int m = 0; m < k; ++m) {
            int med = medoids[m];
            int d   = dist[i * n + med];
            if (d < bestDist) {
                bestDist = d;
                bestMedoidIdx = m;
            }
        }

        clusterOf[i]    = bestMedoidIdx;
        distToMedoid[i] = bestDist;
        totalCost      += bestDist;
    }

    return totalCost;
}

/**
 * @brief Calcule le coût global (MPI) avec coût local parallélisé (OpenMP).
 *
 * @param dist    Matrice des distances (n × n) stockée à plat.
 * @param n       Nombre de sommets.
 * @param medoids Liste des indices de médioïdes.
 * @return Coût total global (somme sur tous les sommets).
 */
static long long computeCostDistributed(const vector<int>& dist,
                                        int n,
                                        const vector<int>& medoids)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int k = (int)medoids.size();

    // Découpage des sommets entre les processus MPI
    int chunk = (n + size - 1) / size;
    int start = rank * chunk;
    int end   = min(n, start + chunk);

    long long localCost = 0;

    // Calcul local parallélisé (OpenMP), puis réduction MPI.
    #pragma omp parallel for reduction(+:localCost) schedule(static)
    for (int i = start; i < end; ++i) {
        int bestDist = INF;

        for (int m = 0; m < k; ++m) {
            int med = medoids[m];
            int d   = dist[i * n + med];
            if (d < bestDist) {
                bestDist = d;
            }
        }

        localCost += bestDist;
    }

    // Réduction MPI pour obtenir le coût global
    long long globalCost = 0;
    MPI_Allreduce(&localCost, &globalCost, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    return globalCost;
}

/**
 * @brief Version hybride MPI+OpenMP de l'algorithme PAM.
 *
 * L'algorithme fonctionne en trois phases :
 *  1. Initialisation aléatoire des médioïdes (rang 0).
 *  2. Boucle d'amélioration : test de tous les échanges possibles.
 *  3. Affectation finale des sommets aux clusters.
 *
 * Parallélisation à deux niveaux :
 *  - MPI : distribution du calcul de coût entre processus.
 *  - OpenMP : parallélisation du coût local dans computeCostDistributed.
 *
 * @param dist Matrice des distances (n × n) stockée à plat.
 * @param n    Nombre de sommets.
 * @param k    Nombre de clusters (médioïdes).
 * @return Résultat PAM (médioïdes, clusters, coût).
 */
PAMResult runPAM_MPI(const vector<int>& dist, int n, int k) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int nthreads = omp_get_max_threads();
    
    if (rank == 0) {
        cout << "[INFO] PAM hybride : " << size << " processus MPI, "
             << nthreads << " threads OpenMP/processus" << endl;
        cout << "[INFO] Configuration : n=" << n << ", k=" << k << endl;
    }

    PAMResult res;
    res.medoids.resize(k);

    // =========================================================================
    // Phase 1 : Initialisation des médioïdes
    // =========================================================================
    
    if (rank == 0) {
        vector<int> allIndices(n);
        for (int i = 0; i < n; ++i) allIndices[i] = i;

        unsigned seed = (unsigned) time(NULL);
        mt19937 gen(seed);
        shuffle(allIndices.begin(), allIndices.end(), gen);

        for (int m = 0; m < k; ++m) {
            res.medoids[m] = allIndices[m];
        }
    }

    // Diffusion des médioïdes initiaux
    MPI_Bcast(res.medoids.data(), k, MPI_INT, 0, MPI_COMM_WORLD);

    // Calcul du coût initial
    long long bestCost = computeCostDistributed(dist, n, res.medoids);

    if (rank == 0) {
        cout << "[INFO] Coût initial : " << bestCost << endl;
    }

    // =========================================================================
    // Phase 2 : Boucle d'amélioration
    // =========================================================================
    
    int iteration = 0;
    while (true) {
        bool improved = false;
        long long bestCostThisPass = bestCost;
        vector<int> bestMedoidsThisPass = res.medoids;

        // Pour chaque médioïde m, tester tous les remplacements possibles
        for (int m = 0; m < k; ++m) {
            long long localBestCost = bestCostThisPass;
            int localBestH = -1;

            // Boucle séquentielle sur h (MPI_Allreduce reste hors OpenMP)
            for (int h = 0; h < n; ++h) {
                bool isMedoid = false;
                for (int mm = 0; mm < k; ++mm) {
                    if (res.medoids[mm] == h) {
                        isMedoid = true;
                        break;
                    }
                }
                if (isMedoid) continue;

                vector<int> newMedoids = res.medoids;
                newMedoids[m] = h;

                long long newCost = computeCostDistributed(dist, n, newMedoids);
                if (newCost < localBestCost) {
                    localBestCost = newCost;
                    localBestH = h;
                }
            }
            
            // Le rang 0 décide de la meilleure amélioration
            if (rank == 0 && localBestH != -1 && localBestCost < bestCostThisPass) {
                bestCostThisPass = localBestCost;
                bestMedoidsThisPass = res.medoids;
                bestMedoidsThisPass[m] = localBestH;
                improved = true;
            }
        }

        // Diffuser si une amélioration a été trouvée
        int flag = improved ? 1 : 0;
        MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!flag) {
            // Convergence : aucune amélioration trouvée
            if (rank == 0) {
                cout << "[INFO] Convergence après " << iteration << " itérations" << endl;
            }
            break;
        }

        // Mise à jour des médioïdes
        if (rank == 0) {
            res.medoids = bestMedoidsThisPass;
            bestCost = bestCostThisPass;
            cout << "[INFO] Itération " << iteration
                 << " : nouveau coût = " << bestCost << endl;
        }

        MPI_Bcast(res.medoids.data(), k, MPI_INT, 0, MPI_COMM_WORLD);
        iteration++;
    }

    // =========================================================================
    // Phase 3 : Affectation finale
    // =========================================================================
    
    if (rank == 0) {
        res.clusterOf.resize(n);
        res.distToMedoid.resize(n);

        long long finalCost = computeCostAndAssign(dist, n,
                                                   res.medoids,
                                                   res.clusterOf,
                                                   res.distToMedoid);
        res.totalCost = finalCost;
        
        cout << "[INFO] Coût final : " << finalCost << endl;
    }

    return res;
}
