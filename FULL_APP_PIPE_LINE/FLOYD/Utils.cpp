#include "Utils.hpp"
#include <fstream>   // pour ofstream

using namespace std;

/**
 * @file Utils.cpp
 * @brief Implémentation des fonctions utilitaires d'affichage et d'écriture.
 */

/**
 * @brief Affiche une matrice d'entiers sur la sortie standard.
 *
 * @param tab    Pointeur vers la matrice stockée à plat.
 * @param n      Nombre de lignes.
 * @param m      Nombre de colonnes.
 * @param format Largeur minimale d'affichage.
 */
void affichage(int* tab, int n, int m, int format) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            cout << setw(format) << tab[i * m + j] << " ";
        cout << endl;
    }
}

/**
 * @brief Écrit une matrice n × m dans un fichier texte.
 *
 * Format du fichier :
 *  - ligne 1 : "n m"
 *  - lignes suivantes : la matrice, une ligne par i.
 *
 * @param filename Nom du fichier de sortie.
 * @param tab      Pointeur vers la matrice.
 * @param n        Nombre de lignes.
 * @param m        Nombre de colonnes.
 * @param format   Largeur minimale d'affichage (setw).
 */
void writeMatrixToFile(const string& filename,
                       const int* tab,
                       int n, int m,
                       int format)
{
    ofstream out(filename);
    if (!out) {
        cerr << "[ERREUR] Impossible d'ouvrir le fichier " << filename << " en écriture.\n";
        return;
    }

    // On écrit d'abord les dimensions
    out << n << " " << m << "\n";

    // Puis la matrice, même logique que dans affichage()
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (format > 0)
                out << setw(format) << tab[i * m + j] << " ";
            else
                out << tab[i * m + j] << " ";
        }
        out << "\n";
    }
}
