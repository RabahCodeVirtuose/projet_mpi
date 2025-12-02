#define OMPI_SKIP_MPICXX 1

/**
 * @file main_mpi.cpp
 * @brief Point d'entrée du programme MPI pour le calcul des plus courts chemins (Floyd–Warshall par blocs).
 *
 * Le rang 0 lit un graphe au format Graphviz (.dot), construit la matrice d'adjacence
 * correspondante, puis diffuse cette matrice à tous les processus MPI. L'algorithme
 * parallèle de Floyd–Warshall par blocs est ensuite lancé via ParallelFloydWarshallBlocks().
 *
 * À la fin du calcul, le rang 0 :
 *  - récupère la matrice finale des distances,
 *  - la sauvegarde dans un fichier texte utilisable par l'algorithme PAM,
 *  - affiche le temps d'exécution de la partie parallèle.
 */

#include <mpi.h>
#include <iostream>
#include <string>
#include <map>

#include "Utils.hpp"
#include "ForGraphMPI.hpp"
#include "ParallelFWBlocks.hpp"

using namespace std;

/**
 * @brief Point d'entrée du programme MPI.
 *
 * Usage :
 * @code
 *   mpirun -np <nb_processus> ./main_mpi fichier.dot
 * @endcode
 *
 * @param argc 
 * @param argv 
 *
 * @return 
 */
int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    // Vérification des arguments
    if (argc != 2) {
        if (rank == 0)
            cout << "Usage : mpirun -np X<=6 ./main_mpi fichier.dot\n";
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    char* file_name = argv[1];

    int nb_nodes;
    int* mat_adjacence = nullptr;
    map<string,int> my_nodes;  

    // -----------------------------
    // Rank 0 lit le graphe
    // -----------------------------
    if (rank == 0) {
        mat_adjacence = lectureGrapheMPI(file_name, &nb_nodes, &my_nodes);

        // Debug éventuel : affichage de la matrice d'adjacence
        // cout << "=== Matrice d'adjacence ===\n";
        // affichage(mat_adjacence, nb_nodes, nb_nodes, 3);
        cout << endl;
    }

    // Diffusion du nombre de sommets à tous les processus
    MPI_Bcast(&nb_nodes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocation de la matrice d'adjacence pour les autres rangs
    if (rank != 0)
        mat_adjacence = new int[nb_nodes * nb_nodes];

    // Diffusion de la matrice brute à tous les processus
    MPI_Bcast(mat_adjacence, nb_nodes * nb_nodes, MPI_INT, 0, MPI_COMM_WORLD);

    // -------------------------------------------------
    //   Mesure du temps de l'algorithme parallèle
    //   (lecture + Bcast initiaux exclus)
    // -------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);              // Synchronisation de tous les rangs
    double t_start = MPI_Wtime();

    int* Dk_final = ParallelFloydWarshallBlocks(nb_nodes, mat_adjacence);

    MPI_Barrier(MPI_COMM_WORLD);              // On attend que tout le monde ait fini
    double t_end = MPI_Wtime();

    double elapsed_ms = (t_end - t_start) * 1000.0;

    // Rank 0 : sauvegarde et affichage des résultats
    if (rank == 0) {
        // cout << "\n=== Matrice de distances (MPI) ===\n";
        // affichage(Dk_final, nb_nodes, nb_nodes, 5);

        // Sauvegarde dans un fichier texte pour PAM
        writeMatrixToFile("../../DATA/matrice_finale_sortie_de_floyd_warshal.txt",
                          Dk_final, nb_nodes, nb_nodes, 5);

        cout << "\n>>> Temps d'exécution (parallèle MPI) : "
             << elapsed_ms << " ms" << endl;
        
        cout << "\nRésultat sauvegardé dans "
             << "'../../DATA/matrice_finale_sortie_de_floyd_warshal.txt'\n"
             << endl;     
    }

    delete[] mat_adjacence;
    delete[] Dk_final;   // sur les autres rangs, c'est nullptr => OK

    MPI_Finalize();
    return 0;
}
