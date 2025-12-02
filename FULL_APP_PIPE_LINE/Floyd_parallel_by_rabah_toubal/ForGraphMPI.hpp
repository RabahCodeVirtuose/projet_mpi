#ifndef FORGRAPHMPI_HPP
#define FORGRAPHMPI_HPP

#include <iostream>
#include <map>
#include <string>
#include <graphviz/cgraph.h>

/**
 * @file ForGraphMPI.hpp
 * @brief Fonctions de lecture d'un graphe au format Graphviz (.dot) et
 *        construction de la matrice d'adjacence associée.
 *
 * Ce module fournit une fonction utilitaire permettant :
 *  - de lire un fichier de graphe au format .dot grâce à la bibliothèque Graphviz,
 *  - de numéroter les sommets de 0 à n-1,
 *  - de construire une matrice d'adjacence symétrique (graphe non orienté),
 *  - de remplir une table de correspondance entre les noms textuels des sommets
 *    et leurs indices entiers.
 *
 * Cette matrice d'adjacence est ensuite utilisée comme entrée de l'algorithme
 * parallèle de Floyd–Warshall.
 */

/**
 * @brief Lit un graphe au format .dot et construit sa matrice d'adjacence.
 *
 * La fonction effectue les opérations suivantes :
 *  - ouverture du fichier .dot,
 *  - lecture du graphe avec la bibliothèque Graphviz (cgraph),
 *  - numérotation des sommets de 0 à n-1,
 *  - construction d'une matrice d'adjacence dense de taille n × n,
 *    initialisée à 0 et remplie avec les poids des arêtes.
 *
 * Le graphe est supposé non orienté : pour chaque arête (u, v) de poids w,
 * la matrice est mise à jour en position (i, j) et (j, i).
 *
 * @param f        Chemin du fichier .dot à lire.
 * @param nb_nodes Pointeur vers un entier dans lequel sera stocké
 *                 le nombre total de sommets du graphe.
 * @param my_nodes Pointeur vers une map qui sera remplie avec la correspondance
 *                 entre le nom textuel du sommet (std::string) et son indice entier.
 *
 * @return Un pointeur vers un tableau d'entiers de taille nb_nodes * nb_nodes
 *         représentant la matrice d'adjacence (stockée à plat, ligne par ligne).
 *         La mémoire est allouée avec new[] et doit être libérée par l'appelant.
 */
int* lectureGrapheMPI(char* f, int* nb_nodes, std::map<std::string,int>* my_nodes);

#endif // FORGRAPHMPI_HPP
