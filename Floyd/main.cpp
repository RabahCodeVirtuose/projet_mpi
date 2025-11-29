#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <chrono>   // ajout chrono

#include "Utils.hpp"
#include "ForGraph.hpp"

using namespace std;

int main(int argc, char* argv[]) {

    if (argc != 2) {
        cout << "Usage : ./main fichier.dot (graphe au format dot)";
        return EXIT_FAILURE;
    }

    char* file_name = argv[1];

    map<string,int> my_nodes; // Dans le .dot les sommets sont désignés par un nom 
    
    int nb_nodes;

    int* mat_adjacence = lectureGraphe(file_name,&nb_nodes,&my_nodes);
   
    cout << "matrice d'adjacence" << endl;
    affichage(mat_adjacence,nb_nodes,nb_nodes,2);
    cout << endl;

    // ========================
    // CHRONOMÉTRAGE SANS ARRONDI
    // ========================
    auto start = chrono::high_resolution_clock::now();

    int* Dk = MatDistance(nb_nodes, mat_adjacence);

    auto end = chrono::high_resolution_clock::now();

    // durée précise, en millisecondes, sans arrondi
    chrono::duration<double, milli> duration_ms = end - start;

   
    // ========================

    cout << "La matrice de distances" << endl;
    affichage(Dk,nb_nodes,nb_nodes,3);
 cout << "\n>>> Temps d'exécution (séquentiel) : "
         << duration_ms.count() << " ms" << endl;

    return 0;
}
