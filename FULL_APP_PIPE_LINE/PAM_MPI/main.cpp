#define OMPI_SKIP_MPICXX 1

/**
 * @file main.cpp
 * @brief Point d'entrée du programme PAM hybride (MPI + OpenMP).
 *
 * Ce binaire :
 *  - lit une matrice de distances,
 *  - diffuse les données à tous les processus MPI,
 *  - exécute PAM en mode hybride,
 *  - écrit le résultat sur le rang 0.
 */

#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>

#include "Utils.hpp"
#include "PAM.hpp"

using namespace std;

/**
 * @brief Point d'entrée principal.
 *
 * @param argc Nombre d'arguments.
 * @param argv Arguments : argv[1] = fichier des distances, argv[2] = k (optionnel).
 * @return Code retour du programme.
 */
int main(int argc, char** argv) {
    // =========================================================================
    // INITIALISATION MPI EN MODE THREAD-SAFE
    // =========================================================================
    
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    
    if (provided < MPI_THREAD_FUNNELED) {
        cerr << "[WARN] MPI ne supporte pas MPI_THREAD_FUNNELED\n";
        cerr << "       Le programme peut crasher avec OpenMP+MPI\n";
    }

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Vérification des arguments
    if (argc < 2) {
        if (rank == 0) {
            cerr << "\n[ERREUR] Argument manquant !\n";
            cerr << "Usage: " << (argv[0] ? argv[0] : "./pam_mpi")
                 << " <fichier_distances.txt> [k]\n\n";
            cerr << "Exemples:\n";
            cerr << "  mpirun -np 4 ./pam_mpi ../../DATA/matrice_finale.txt\n";
            cerr << "  mpirun -np 4 ./pam_mpi ../../DATA/matrice_finale.txt 5\n\n";
        }
        MPI_Finalize();
        return 1;
    }

    // Paramètres
    string distFile = argv[1];
    
    int k = 4;
    if (argc >= 3) {
        k = atoi(argv[2]);
        if (k <= 0) {
            if (rank == 0) {
                cerr << "[ERREUR] k doit être > 0, valeur reçue : " << argv[2] << "\n";
            }
            MPI_Finalize();
            return 1;
        }
    }
    
    const string outFile = "../../DATA/resultat_pam_parallel.txt";

    // Lecture et diffusion de la matrice
    int n = 0;
    vector<int> dist;

    if (rank == 0) {
        cout << "\n========================================\n";
        cout << "  PAM Parallèle (MPI+OpenMP)\n";
        cout << "========================================\n";
        cout << "Processus MPI      : " << size << "\n";
        cout << "Threads OpenMP/proc: " << omp_get_max_threads() << "\n";
        cout << "Support MPI threads: " << (provided >= MPI_THREAD_FUNNELED ? "OUI" : "NON") << "\n";
        cout << "Fichier d'entrée   : " << distFile << "\n";
        cout << "Nombre de clusters : " << k << "\n";
        cout << "========================================\n\n";

        try {
            dist = readDistanceMatrix(distFile, n);
            cout << "[INFO] Matrice de distances chargée : " << n << "x" << n << "\n";
            
            if (k > n) {
                cerr << "[ERREUR] k (" << k << ") ne peut pas être > n (" << n << ")\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
        } catch (const exception& e) {
            cerr << "[ERREUR] Lecture de la matrice : " << e.what() << "\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Diffusion
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        dist.resize(n * n);
    }

    MPI_Bcast(dist.data(), n * n, MPI_INT, 0, MPI_COMM_WORLD);

    // Exécution de PAM
    if (rank == 0) {
        cout << "[INFO] Démarrage de PAM...\n\n";
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    PAMResult res = runPAM_MPI(dist, n, k);

    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    double elapsed_s = t1 - t0;
    double elapsed_ms = elapsed_s * 1000.0;

    // Affichage des résultats
    if (rank == 0) {
        cout << "\n========================================\n";
        cout << "  Résultats PAM\n";
        cout << "========================================\n";
        cout << "Coût total    : " << res.totalCost << "\n";
        cout << "Médioïdes     : ";
        for (size_t m = 0; m < res.medoids.size(); ++m) {
            cout << res.medoids[m];
            if (m < res.medoids.size() - 1) cout << ", ";
        }
        cout << "\n";
        cout << "========================================\n";
        cout << "Temps d'exécution : " << elapsed_s << " s (" << elapsed_ms << " ms)\n";
        cout << "========================================\n\n";

        try {
            writePAMResult(outFile, res);
            cout << "[INFO] Résultats sauvegardés dans : " << outFile << "\n\n";
        } catch (const exception& e) {
            cerr << "[ERREUR] Écriture du résultat : " << e.what() << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}
