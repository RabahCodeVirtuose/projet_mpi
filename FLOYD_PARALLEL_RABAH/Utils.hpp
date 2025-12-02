#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <iomanip>
#include <string>

/**
 * @file Utils.hpp
 * @brief Fonctions utilitaires pour l'affichage et l'écriture de matrices.
 *
 * Ce module regroupe deux fonctions :
 *  - affichage d'une matrice stockée à plat sur la sortie standard,
 *  - écriture d'une matrice dans un fichier texte, au format utilisé
 *    par le projet (première ligne : "n m", puis n lignes de m valeurs).
 *
 * Ces fonctions sont utilisées à la fois pour le débogage et pour
 * sauvegarder les résultats du Floyd–Warshall afin d'alimenter PAM.
 */

/**
 * @brief Affiche une matrice d'entiers sur la sortie standard.
 *
 * La matrice est supposée stockée en mémoire de manière linéaire
 * (row-major). L'élément (i, j) correspond donc à tab[i * m + j].
 *
 * @param tab    Pointeur vers le tableau contenant la matrice.
 * @param n      Nombre de lignes de la matrice.
 * @param m      Nombre de colonnes de la matrice.
 * @param format Largeur minimale utilisée lors de l'affichage de chaque valeur
 *               (paramètre utilisé par std::setw). Cela permet d'aligner
 *               les colonnes pour une meilleure lisibilité.
 */
void affichage(int* tab, int n, int m, int format);

/**
 * @brief Écrit une matrice d'entiers dans un fichier texte.
 *
 * Format du fichier généré :
 *   - Première ligne : "n m"
 *   - Puis n lignes, chacune contenant m entiers séparés par un espace.
 *
 * La matrice est supposée être stockée à plat (row-major).
 *
 * @param filename Nom du fichier de sortie (ex : "distances.txt").
 * @param tab      Pointeur vers la matrice à écrire.
 * @param n        Nombre de lignes de la matrice.
 * @param m        Nombre de colonnes de la matrice.
 * @param format   Largeur minimale utilisée pour l'affichage de chaque entier.
 *                 Si format > 0, std::setw(format) est appliqué.
 *                 Si format == 0, les valeurs sont écrites sans largeur fixe.
 */
void writeMatrixToFile(const std::string& filename,
                       const int* tab,
                       int n, int m,
                       int format);

#endif // UTILS_HPP
