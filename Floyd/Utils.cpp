#include "Utils.hpp"

void affichage (int* tab, int n, int m, int format) {
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++)
            cout << setw(format) << tab[i*m+j] << " ";
        cout << endl;
    }
}