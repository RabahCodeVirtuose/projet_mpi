
## 1. Description rapide

Ce dossier contient ma version parallèle de l’algorithme de Floyd-Warshall.
J’ai utilisé **MPI** pour découper la grosse matrice en blocs et distribuer le calcul entre plusieurs processus.

L’idée générale est la même que Floyd classique mais, au lieu de tout faire sur un seul processus, chaque bloc est envoyé au processus qui doit le traiter.
On utilise aussi des `MPI_Bcast` (et `MPI_Ibcast`) pour partager les blocs pivots.

Je n’ai pas cherché à faire plus compliqué que ce qu’on a vu en cours, juste la version avec grille de processus et mises à jour par blocs.

---

## 2. Fichiers importants

* `main_mpi.cpp` → le programme principal
* `ParallelFWBlocks.cpp` + `.hpp` → l’algorithme en lui-même (version blocs)
* `Distribution.cpp` → qui gère la répartition des blocs entre les processus
* `Utils.cpp` → quelques fonctions pratiques (affichage, écriture fichier…)
* `Makefile` → pour compiler

---

## 3. Compilation

Normalement il suffit d’être dans ce dossier (`FLOYD_PARALLEL_RABAH/`) et de faire :

```
make
```

Le Makefile compile tout et génère un exécutable :

```
./main_mpi
```

---

## 4. Exécution

Le programme prend en entrée un fichier contenant **la matrice d’adjacence** (format simple : n, puis n lignes de n entiers).

Pour lancer le programme avec MPI :

```
mpirun -np <nb_processus> ./main_mpi <fichier_matrice>
```

Exemple (4 processus) :

```
mpirun -np 4 ./main_mpi ../../DATA/PetitExemple.dot
```

Ou pour la grande matrice :

```
mpirun -np 9 ./main_mpi ../../DATA/exemplemassi.dot
```

*(Les fichiers sont dans le dossier `DATA` à la racine du projet.)*

---

## 5. Sortie du programme

À la fin, le programme :

* calcule la matrice des plus courts chemins,
* rassemble tout sur le rang 0,
* écrit le résultat dans `DATA/matrice_finale_sortie_de_floyd_warshal.txt`

Le rang 0 affiche aussi le temps total d’exécution.

---

## 6. Remarque

* Si le nombre de processus n’est **pas carré**, le programme ajuste la taille des blocs pour quand même faire tourner l’algo (grâce à MPI_Dims_create).
* Si la taille n’est pas divisible par la taille des blocs, les cases en trop sont juste remplies avec **INF** (padding classique).

---

## 7. Exemple complet

```
make
mpirun -np 4 ./main_mpi ../DATA/PetitExemple.dot
```

---

## 8. Nettoyage

Pour supprimer les `.o` et tout recompiler propre :

```
make clean
```

---
