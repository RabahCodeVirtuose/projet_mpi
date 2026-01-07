Master 1 ARIAS - Programmation parallele
Parallele MPI + OpenMP : Needleman-Wunsch, Floyd-Warshall, PAM
Rabah Toubal


1. Introduction

L'objectif de ce projet est de construire un pipeline complet pour traiter des sequences ARN. Dans la premiere version, la distance de Hamming etait utilisee pour construire le graphe. Dans cette version, je remplace Hamming par le score Needleman-Wunsch (alignement global avec penalites de trous affines) et j'introduis de l'hybridation MPI + OpenMP dans les etapes qui s'y pretent.

Pipeline final :
FASTA -> Needleman-Wunsch -> Graphe DOT -> Floyd-Warshall -> PAM

Je travaille avec des sequences de meme longueur. Le graphe est filtre avec un seuil epsilon (70 dans mes tests), comme dans la partie 1. Le rapport decrit (1) la nouvelle mesure de similarite et sa parallelisation, et (2) les optimisations hybrides dans le reste du pipeline.


2. Remplacement de Hamming par Needleman-Wunsch

2.1 Score Needleman-Wunsch (gaps affines)

Le score est calcule avec les parametres imposes :
- match = +1
- mismatch = -1
- gap opening = -3
- gap extension = -1

Pour gerer les gaps affines, on utilise 3 matrices de programmation dynamique :
- M[i][j] : meilleur score si A[i] et B[j] sont alignees (match/mismatch)
- X[i][j] : meilleur score si A[i] est aligne avec un trou dans B (gap vertical)
- Y[i][j] : meilleur score si B[j] est aligne avec un trou dans A (gap horizontal)

Recurrences (idee) :
- M[i][j] = max(M, X, Y sur la diagonale) + score(match/mismatch)
- X[i][j] = max( M[i-1][j] + gap_open, X[i-1][j] + gap_extend )
- Y[i][j] = max( M[i][j-1] + gap_open, Y[i][j-1] + gap_extend )

Le score final est : max(M[n][m], X[n][m], Y[n][m]).

2.2 Initialisation

Les bords sont initialises avec la penalite d'ouverture puis l'extension :
- ligne 0 : gap_open + j * gap_extend
- colonne 0 : gap_open + i * gap_extend

Cela force un alignement global et respecte bien la difference entre ouverture et extension de gap.

2.3 Conversion score -> distance

Pour garder la meme logique que dans la premiere partie (comparaison avec un epsilon), je transforme le score en distance :

d = (L - score) / 2

Explication rapide : si on avait seulement match/mismatch, le score vaut L - 2 * (nb_mismatch). Donc d correspond a un equivalent de distance. Avec les gaps, c'est une approximation mais elle reste coherente pour classer les sequences et filtrer les aretes du graphe.

2.4 Complexite

Pour deux sequences de longueur L, Needleman-Wunsch coute O(L^2) en temps et O(L^2) en memoire. Pour N sequences, on calcule toutes les paires, donc O(N^2 * L^2). C'est la partie la plus lourde du pipeline.

2.5 Parallelisation OpenMP

Dans cette version, le coeur de Needleman-Wunsch reste sequentiel (les 3 matrices M/X/Y sont remplies par lignes). La parallelisation se fait au niveau des paires de sequences : on distribue les lignes i entre les threads OpenMP (boucle externe), et chaque thread calcule tous les alignements (i, j) associes. C'est exactement ce qui est implemente dans le code via un `#pragma omp parallel for` sur i.

Cette approche est simple, evite le sur-parallellisme et reste efficace car le nombre de paires est tres grand. Elle respecte aussi les dependances internes de l'algorithme, qui rendent difficile une parallelisation fine du remplissage de la matrice.

Pour gagner du temps, je calcule seulement la partie triangulaire (i < j) puis je recopie la matrice symetrique. Cela divise presque par 2 le temps de calcul, sans changer le resultat.

2.6 Generation du graphe DOT

Une fois la matrice de distances obtenue, on ecrit un graphe pond ere. Une arete (i,j) est ajoutee si la distance est inferieure a epsilon. Cette etape permet de filtrer les aretes faibles et de limiter la taille du graphe.


3. Floyd-Warshall par blocs (MPI + OpenMP)

3.1 Rappel du principe

La matrice des distances est decoupee en blocs. Une grille 2D de processus MPI se partage ces blocs. A chaque etape k :
- bloc pivot (k,k) mis a jour puis diffuse
- mise a jour de la ligne k et colonne k, puis diffusion
- mise a jour des blocs internes avec les blocs de ligne/colonne

L'idee est d'eviter les dependances globales et de travailler sur des blocs locaux. Cela reduit les communications et permet un meilleur cache.

3.2 Distribution et proprietaires

Une grille 2D (Pr x Pc) est creee avec MPI_Dims_create. Chaque bloc (bi,bj) est assigne a un processus. Les fonctions ownerOf et computeLocalBlocks determinent les blocs locaux. Chaque processus ne stocke que ses blocs, et les autres sont reconstruits via communications.

3.3 Optimisation de la taille des blocs

J'ai ajuste la taille de bloc b pour eviter des blocs trop gros qui cassent le cache et reduisent le parallele. La logique est :
- si la grille est carree (p carre parfait), b = n / sqrt(p) mais b <= 128
- sinon b = ceil(n / sqrt(p)) puis b borne entre 32 et 128

Cette regle evite des blocs enormes (ex: 1000x1000 quand p=4) et donne une meilleure repartition du travail. Cela stabilise les temps et reduit les variations selon le nombre de processus.

3.4 OpenMP dans Floyd

A l'interieur de chaque processus MPI, les mises a jour locales de blocs (fw_block, fw_row, fw_col, fw_inner) sont parallelisees avec OpenMP. L'idee est de garder MPI pour le gros grain (distribution des blocs) et OpenMP pour les boucles internes. Les communications utilisent des broadcast (souvent non bloquants) pour diffuser les blocs pivot, puis on calcule en local.


4. PAM hybride (MPI + OpenMP)

4.1 Rappel de PAM

PAM choisit k medioids parmi n points. A chaque iteration, on teste des echanges (m, h) : on remplace un medioid par un point non-medoid, puis on recalcule le cout global. L'algorithme s'arrete quand aucun echange ne reduit le cout.

4.2 Parallelisation

Le cout total est la somme, pour chaque point, de la distance au medioid le plus proche.
- MPI : chaque processus calcule le cout sur une partie des lignes, puis MPI_Allreduce somme tout.
- OpenMP : le cout local est parallelise avec une reduction sur les lignes.

4.3 Hybridation MPI + OpenMP

Pour autoriser l'usage de threads avec MPI, j'initialise avec MPI_Init_thread et je verifie la valeur returned (provided >= MPI_THREAD_FUNNELED). Dans ce mode, les appels MPI restent en dehors des regions OpenMP, ce qui garde un comportement correct.


5. Resultats (machine personnelle)

5.1 Configuration

Machine : Intel i5-13420H, 6 coeurs / 12 threads.
Jeu de test : 2000 sequences, longueur L=100.

Configuration finale testee :
- Needleman : OpenMP 12 threads
- Floyd     : 4 processus MPI x 3 threads OpenMP
- PAM       : 4 processus MPI x 3 threads OpenMP

5.2 Temps observes (meilleur cas)

- Needleman-Wunsch (OpenMP) : 4.62 s
- Floyd-Warshall (MPI+OMP)  : 0.813 s
- PAM (MPI+OMP)             : 0.422 s

Le temps total du pipeline est d'environ 5.85 s. On voit que la partie la plus couteuse reste le calcul des distances par Needleman-Wunsch, ce qui est logique car on calcule n^2 paires.

5.3 Observations

- Augmenter le nombre de threads au-dela des coeurs physiques ne donne pas toujours un gain a cause de l'overhead et de la memoire.
- La taille des blocs dans Floyd-Warshall a un impact tres fort. Sans bloc limite, les performances chutent.
- PAM passe bien a l'echelle jusqu'a quelques processus, puis le gain ralentit a cause des reductions MPI et du cout fixe de communication.


6. Conclusion

Cette deuxieme version du projet remplace la distance de Hamming par un score d'alignement global plus realiste. J'ai implemente Needleman-Wunsch avec gaps affines (3 matrices) et j'ai parallelise le calcul en distribuant les paires de sequences entre les threads OpenMP. Le coeur de l'algorithme reste sequentiel, ce qui est plus simple et stable pour cette etape.

Dans le pipeline complet, l'approche hybride MPI + OpenMP permet de mieux exploiter la machine : MPI distribue les blocs ou les lignes, OpenMP accelere les calculs locaux. L'optimisation de la taille des blocs dans Floyd-Warshall a un impact important sur les performances et donne une execution plus stable.

Ce projet m'a permis de comprendre concretement comment combiner MPI et OpenMP sur un pipeline complet, et de voir les limites liees au cout des calculs et des communications.
