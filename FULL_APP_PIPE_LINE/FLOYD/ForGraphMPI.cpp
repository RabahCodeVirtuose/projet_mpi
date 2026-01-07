#define OMPI_SKIP_MPICXX 1
#include "ForGraphMPI.hpp"

/**
 * @file ForGraphMPI.cpp
 * @brief Implémentation de la lecture d'un graphe .dot et construction de la matrice.
 */

using namespace std;

/**
 * @brief Lit un graphe .dot et construit la matrice d'adjacence.
 *
 * @param f        Chemin du fichier .dot.
 * @param nb_nodes Pointeur pour récupérer le nombre de sommets.
 * @param my_nodes Table de correspondance nom -> indice.
 * @return Matrice d'adjacence allouée dynamiquement (new[]).
 */
int* lectureGrapheMPI(char* f, int* nb_nodes, map<string,int>* my_nodes) {

    FILE *fp = fopen(f, "r");
    if (!fp) {
        cout << "Impossible d'ouvrir le fichier " << f << endl;
        exit(1);
    }




    
    Agraph_t *g = agread(fp, NULL);
    fclose(fp);

    int nn = agnnodes(g);
    (*nb_nodes) = nn;
    
    int t = 0;
    for (Agnode_t *n = agfstnode(g); n; n = agnxtnode(g, n)) {
        (*my_nodes)[agnameof(n)] = t;
        t++;
    }

    int* mat = new int[nn * nn]();
    
    for (Agnode_t *n = agfstnode(g); n; n = agnxtnode(g, n)) {
        int i = (*my_nodes)[agnameof(n)];
        for (Agedge_t *e = agfstout(g, n); e; e = agnxtout(g, e)) {
            int j = (*my_nodes)[agnameof(aghead(e))];
            int w = stoi(agget(e, (char*)"weight"));
            mat[i*nn + j] = w;
            mat[j*nn + i] = w;
        }
    }

    agclose(g);
    return mat;
}
