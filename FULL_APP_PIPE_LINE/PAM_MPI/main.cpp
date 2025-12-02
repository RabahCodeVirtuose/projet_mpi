#define OMPI_SKIP_MPICXX 1
/**
 * @file main.cpp
 * @brief Point d'entrée du programme MPI pour l'exécution de l'algorithme PAM.
 *
 * Le rang 0 lit une matrice de distances depuis un fichier texte, diffuse cette
 * matrice à tous les processus, puis lance la version MPI de PAM (runPAM_MPI).
 * À la fin, le rang 0 affiche le coût final, les médioïdes, le temps d'exécution,
 * et écrit un fichier de sortie contenant la partition.
 */

#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>

#include "Utils.hpp"
#include "PAM.hpp"

/**
 * @brief Point d'entrée du programme.
 *
 * Usage :
 * @code
 *   mpirun -np <nb_processus> ./pam_mpi fichier_distances.txt
 * @endcode
 *
 * @param argc Nombre d'arguments de la ligne de commande.
 * @param argv Tableau d'arguments (argv[1] doit contenir le fichier de distances).
 * @return 0 en cas de succès, une valeur non nulle en cas d'erreur.
 */
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // ---------- NOMS EN DUR ICI ----------
    const std::string distFile = argv[1];
    const std::string outFile  = "../../DATA/resultat_pam_parallel.txt";                  // résultat PAM
    const int k = 4;                                                                   // nombre de clusters
    // -------------------------------------

    int n = 0;
    std::vector<int> dist;

    if (rank == 0) {
        try {
            dist = readDistanceMatrix(distFile, n);
            std::cout << "\n\nLecture de la matrice de distances: n = " << n << "\n";
            std::cout << "Execution de PAM MPI avec k = " << k << " ...\n";
        } catch (const std::exception& e) {
            std::cerr << "Erreur (rang 0) : " << e.what() << "\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Diffuser n à tous
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allouer la matrice de distances sur les autres rangs
    if (rank != 0) {
        dist.resize(n * n);
    }

    // Diffuser la matrice complète
    MPI_Bcast(dist.data(), n * n, MPI_INT, 0, MPI_COMM_WORLD);

    // ------------------------------------------------------------------
    // Mesure du temps d'exécution de PAM MPI
    // ------------------------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);            // on synchronise tous les processus
    double t0 = MPI_Wtime();                // temps de départ (en secondes)

    // Exécuter l'algorithme PAM en parallèle
    PAMResult res = runPAM_MPI(dist, n, k);

    MPI_Barrier(MPI_COMM_WORLD);            // on attend que tout le monde ait fini
    double t1 = MPI_Wtime();                // temps de fin
    double elapsed_ms = (t1 - t0) * 1000.0; // conversion en millisecondes

    // Rang 0 : afficher + écrire le résultat
    if (rank == 0) {
        std::cout << "Cout final = " << res.totalCost << "\n";
        std::cout << "Medioides : ";
        for (int m = 0; m < (int)res.medoids.size(); ++m) {
            std::cout << res.medoids[m] << " ";
        }
        std::cout << "\n";

        std::cout << "\n\n>>> Temps d'execution (PAM MPI) : "
                  << elapsed_ms << " ms\n\n";

        try {
            writePAMResult(outFile, res);
            std::cout << "Resultats ecrits dans " << outFile << "\n";
        } catch (const std::exception& e) {
            std::cerr << "Erreur d'ecriture du resultat: " << e.what() << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}
