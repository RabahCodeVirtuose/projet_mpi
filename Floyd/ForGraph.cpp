#include "ForGraph.hpp"

int* lectureGraphe(char* f, int* nb_nodes, map<string, int>* my_nodes) {
    FILE *fp = fopen(f, "r");
    Agraph_t *g = agread(fp, NULL);
    fclose(fp);

    int nn = agnnodes(g);

    (*nb_nodes) = nn;

    int t = 0;
    for (Agnode_t *n = agfstnode(g); n; n = agnxtnode(g, n)) {
        (*my_nodes)[agnameof(n)]=t;
        t++;
    }


    /* La matrice d'adjacence est "normale" sans l'infini
     * C'est dans la construction de la matrice de distance initiale
     * qu'on initialise à l'infini une case non diagonale qui ne 
     * correspond pas à une arête du graphe. 
     */
    int* mat_adjacence = new int[nn*nn]();
    
    int i, j;
    int weight;
    for (Agnode_t *n = agfstnode(g); n; n = agnxtnode(g, n)) {
        i = (*my_nodes)[agnameof(n)];
        for (Agedge_t *e = agfstout(g, n); e; e = agnxtout(g, e)) {
            j = (*my_nodes)[agnameof(aghead(e))];
            weight = std::stoi(agget(e, (char*)"weight"));
            mat_adjacence[i*nn+j]=weight;
            mat_adjacence[j*nn+i]=weight;
     }
    }
    agclose(g);
    return mat_adjacence;
}

int* InitDk(int nb_nodes, int* mat_adjacence) {
    int* Dk = new int[nb_nodes*nb_nodes];
    for (int i=0; i<nb_nodes; i++)
        for (int j=0; j<nb_nodes; j++)
            if (i==j)
                Dk[i*nb_nodes+j] = 0;
            else if (mat_adjacence[i*nb_nodes+j]==0)
                Dk[i*nb_nodes+j] = INF;
            else
                Dk[i*nb_nodes+j] = mat_adjacence[i*nb_nodes+j];
    return Dk;
}

int* MatDistance(int nb_nodes, int* mat_adjacence) {
    int* Dk = InitDk(nb_nodes, mat_adjacence);
    int* Dk1 = new int[nb_nodes*nb_nodes]; 

    for (int k=0; k<nb_nodes; k++) {
        for (int i=0; i<nb_nodes; i++)
            for (int j=0; j<nb_nodes; j++) {
                if (Dk[i*nb_nodes+j] <= Dk[i*nb_nodes+k]+Dk[k*nb_nodes+j]) 
                    Dk1[i*nb_nodes+j] = Dk[i*nb_nodes+j];
                else 
                    Dk1[i*nb_nodes+j] = Dk[i*nb_nodes+k]+Dk[k*nb_nodes+j];
                
            }
        for (int i=0; i<nb_nodes*nb_nodes; i++)
            Dk[i] = Dk1[i];
    }
    delete[] Dk1;
    return Dk;
}

int* MatDistance(int nb_nodes, int* mat_adjacence, int* next) {
    int* Dk = InitDk(nb_nodes, mat_adjacence);
    int* Dk1 = new int[nb_nodes*nb_nodes]; 

    for (int k=0; k<nb_nodes; k++) {
        for (int i=0; i<nb_nodes; i++)
            for (int j=0; j<nb_nodes; j++) {
                if (Dk[i*nb_nodes+j] <= Dk[i*nb_nodes+k]+Dk[k*nb_nodes+j]) 
                    Dk1[i*nb_nodes+j] = Dk[i*nb_nodes+j];
                else {
                    Dk1[i*nb_nodes+j] = Dk[i*nb_nodes+k]+Dk[k*nb_nodes+j];
                    next[i*nb_nodes+j] = next[i*nb_nodes+k];
                }
            }
        for (int i=0; i<nb_nodes*nb_nodes; i++)
            Dk[i] = Dk1[i];
    }
    delete[] Dk1;
    return Dk;
}

int* InitNext(int nb_nodes, int* mat_adjacence) {
    int* next = new int[nb_nodes*nb_nodes];
    for (int i=0; i<nb_nodes; i++)
        for (int j=0; j<nb_nodes; j++)
            if (mat_adjacence[i*nb_nodes+j]!=0)
                next[i*nb_nodes+j]=j;
            else
                next[i*nb_nodes+j]=-1;

    return next;
}


