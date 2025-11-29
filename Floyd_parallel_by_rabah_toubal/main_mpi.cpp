#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <iostream>
#include <string>
#include <map>

#include "Utils.hpp"
#include "ForGraphMPI.hpp"
#include "ParallelFWBlocks.hpp"

using namespace std;

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

  //  if (size != 6 && rank == 0) {
    //    cout << "[ERREUR] Ce programme doit être lancé avec 6 processus MPI.\n";
      //  MPI_Finalize();
       // return EXIT_FAILURE;
   // }

    if (argc != 2) {
        if (rank == 0)
            cout << "Usage : mpirun -np 6 ./main_mpi fichier.dot\n";
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

        cout << "=== Matrice d'adjacence ===\n";
        affichage(mat_adjacence, nb_nodes, nb_nodes, 3);
        cout << endl;
    }

    // Diffusion du nb_nodes
    MPI_Bcast(&nb_nodes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocation de mat_adjacence pour les autres ranks
    if (rank != 0)
        mat_adjacence = new int[nb_nodes * nb_nodes];

    // Diffusion de la matrice brute
    MPI_Bcast(mat_adjacence, nb_nodes * nb_nodes, MPI_INT, 0, MPI_COMM_WORLD);

    // -------------------------------------------------
    //   Mesure du temps de l'algorithme parallèle
    //   (on exclut la lecture + les MPI_Bcast initiaux)
    // -------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);              // Synchronisation de tous les rangs
    double t_start = MPI_Wtime();

    int* Dk_final = ParallelFloydWarshallBlocks(nb_nodes, mat_adjacence);

    MPI_Barrier(MPI_COMM_WORLD);              // On attend que tout le monde ait fini
    double t_end = MPI_Wtime();

    double elapsed_ms = (t_end - t_start) * 1000.0;

    // Rank 0 affichage
    if (rank == 0) {
        cout << "\n=== Matrice de distances (MPI) ===\n";
        affichage(Dk_final, nb_nodes, nb_nodes, 5);

        cout << "\n>>> Temps d'exécution (parallèle MPI) : "
             << elapsed_ms << " ms" << endl;
    }

    delete[] mat_adjacence;
    delete[] Dk_final;   // sur les autres rangs, c'est nullptr => OK

    MPI_Finalize();
    return 0;
}
