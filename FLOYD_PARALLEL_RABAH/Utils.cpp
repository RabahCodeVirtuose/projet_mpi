#include "Utils.hpp"
#include <fstream>   // pour std::ofstream

void affichage(int* tab, int n, int m, int format) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            std::cout << std::setw(format) << tab[i * m + j] << " ";
        std::cout << std::endl;
    }
}

// Écrit une matrice n × m dans un fichier texte.
//
// Format du fichier :
//   ligne 1 : n m
//   lignes suivantes : la matrice, une ligne par i, éléments séparés par des espaces.
//
// Exemple de contenu :
//   4 4
//   0 1 3 9
//   1 0 2 4
//   3 2 0 5
//   9 4 5 0
//
void writeMatrixToFile(const std::string& filename,
                       const int* tab,
                       int n, int m,
                       int format)
{
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "[ERREUR] Impossible d'ouvrir le fichier " << filename << " en écriture.\n";
        return;
    }

    // On écrit d'abord les dimensions
    out << n << " " << m << "\n";

    // Puis la matrice, même logique que dans affichage()
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (format > 0)
                out << std::setw(format) << tab[i * m + j] << " ";
            else
                out << tab[i * m + j] << " ";
        }
        out << "\n";
    }
}
