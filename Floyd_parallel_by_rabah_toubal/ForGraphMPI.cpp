#include "ForGraphMPI.hpp"

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
