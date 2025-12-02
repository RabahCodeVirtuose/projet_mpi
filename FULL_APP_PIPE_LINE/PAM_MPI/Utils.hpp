#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <string>
#include "PAM.hpp"

/**
 * @file Utils.hpp
 * @brief Fonctions utilitaires pour l'affichage, la lecture/écriture
 *        de matrices de distances et de résultats PAM.
 */

/**
 * @brief Affiche une matrice n x m stockée à plat dans tab.
 *
 * @param tab    Tableau de taille n*m (accès tab[i*m + j]).
 * @param n      Nombre de lignes.
 * @param m      Nombre de colonnes.
 * @param format Largeur d'affichage pour chaque entier (par défaut 4).
 */
void affichage(const int* tab, int n, int m, int format = 4);

/**
 * @brief Lit une matrice de distances carrée dans un fichier texte.
 *
 * Format attendu (comme ce que produit le Floyd) :
 * @code
 *   n m
 *   d00 d01 ... d0(n-1)
 *   d10 d11 ... d1(n-1)
 *   ...
 * @endcode
 *
 * On suppose que n == m. On utilise seulement n.
 *
 * @param filename Nom du fichier.
 * @param n_out    Paramètre de sortie : nombre de sommets n.
 *
 * @return Un vecteur de taille n*n contenant la matrice ligne par ligne.
 */
std::vector<int> readDistanceMatrix(const std::string& filename, int& n_out);

/**
 * @brief Écrit le résultat de PAM dans un fichier texte.
 *
 * Format (exemple) :
 * @code
 *   # PAM results
 *   # n = 20
 *   # k = 4
 *   # total_cost = 67
 *
 *   # medoids:
 *   6 1 4 7
 *
 *   # columns: vertex cluster medoid dist
 *   0 1 1 1
 *   1 1 1 0
 *   ...
 * @endcode
 *
 * @param filename Nom du fichier de sortie.
 * @param res      Résultat PAM à écrire.
 */
void writePAMResult(const std::string& filename, const PAMResult& res);

#endif // UTILS_HPP
