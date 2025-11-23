#ifndef FORGRAPH_HPP
#define FORGRAPH_HPP

#include <iostream>
#include <graphviz/cgraph.h>
#include <map>
#include <vector>

/**
 * @brief Pour représenter l'infini 
 * 
 */
#define INF 1000

using namespace std;

/**
 * @brief Lecture du graphe et génération de la matrice d'adjacence adaptée
 *        La librairie libcgraph.so de GraphViz est utilisée pour lire le graphe .dot
 * @param f le fichier .dot contenant le graphe
 * @param nb_nodes le nombre de sommets du graphe (sortie)
 * @param my_nodes ll
 * @return int* en retour la matrice d'adjacence non adaptée
 */
int* lectureGraphe(char* f, int* nb_nodes, map<string, int>* my_nodes);

/**
 * @brief Pour initialiser D^(-1) avant les nb_nodes itérations du calcul
 *        de la matrice de distance
 * 
 * @param nb_nodes le nombre de noeuds
 * @param mat_adjacence la matrice d'adjacence normale
 * @return int* la matrice D^(-1) 
 */
int* InitDk(int nb_nodes, int* mat_adjacence);

/**
 * @brief Calcul de la matrice de distance seule
 * 
 * @param nb_nodes le nombre de noeuds
 * @param mat_adjacence la matrice d'adjacence normale
 * @return int* la matrice des distances autrement dit D à l'itération nb_nodes.
 */
int* MatDistance(int nb_nodes, int* mat_adjacence);

/**
 * @brief 
 * 
 * @param nb_nodes 
 * @param mat_adjacence 
 * @return int* 
 */

/**
 * @brief Si on souhaite connaître le chemin le plus court
 *        il est nécessaire de mémoriser les noeuds visités
 *        qui permettent de construire la distance la plus courte
 *        
 * 
 * @param nb_nodes le nombre de noeuds
 * @param mat_adjacence la matrice d'adjacence
 * @return int* une matrice next telle que s'il existe une arête entre le sommet i et le sommet j
 *        alors next[i] = j
 */
int* InitNext(int nb_nodes, int* mat_adjacence);


/**
 * @brief Calcul de la matrice de distance + la matrice next pour reconstruire 
 *        les chemins les plus courts entre 2 paires de sommets
 * 
 * @param nb_nodes le nombre de noeuds
 * @param mat_adjacence la matrice d'adjacence
 * @param next la matrice next initialisée par la fonction précédente
 * @return int* la matrice des distances autrement dit D à l'itération nb_nodes.
 */
int* MatDistance(int nb_nodes, int* mat_adjacence, int* next);



#endif