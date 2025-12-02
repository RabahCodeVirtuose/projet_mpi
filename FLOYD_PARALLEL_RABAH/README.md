-

#  Version Parallèle de Floyd-Warshall (MPI)

## 1. Description rapide

Ce dossier contient **ma version parallèle** de l’algorithme de Floyd-Warshall, écrite en C++ avec **MPI**.
L’idée générale est la même que l’algo séquentiel : on cherche les plus courts chemins entre tous les sommets.
Sauf qu’ici, **la matrice est découpée en blocs**, et chaque processus s’occupe d’une partie différente.

À chaque itération `k`, le bloc pivot est envoyé aux autres processus (via `MPI_Bcast` ou `MPI_Ibcast`), et chacun met à jour ses blocs locaux. C’est exactement ce qu’on a vu en cours sur la **parallélisation par blocs**.

Je ne suis pas allé chercher quelque chose d’extra : j’ai juste appliqué la méthode classique avec une grille de processus et une mise à jour bloc-par-bloc.

 **Inspiration utilisée**
En travaillant dessus, j’ai aussi regardé un document universitaire qui décrit la même stratégie (découpage 2D, broadcasts, mises à jour locales). Ça m’a aidé à comprendre le schéma global.
Référence : Asmita Gautam, *Parallel Floyd-Warshall Algorithm*, University at Buffalo, 2019.

---

## 2. Fichiers importants

* **`main_mpi.cpp`** → programme principal (lecture, distribution, lancer l’algo)
* **`ParallelFWBlocks.cpp / .hpp`** → cœur de l’algorithme (version blocs)
* **`Distribution.cpp`** → répartit chaque bloc à son processus
* **`Utils.cpp`** → outils : affichage, écriture fichier…
* **`Makefile`** → compilation automatique

---

## 3. Compilation

Place-toi simplement dans le dossier :

```
FLOYD_PARALLEL_RABAH/
```

Puis lance :

```
make
```

Ça génère l’exécutable :

```
./main_mpi
```

---

## 4. Exécution

Le programme prend en entrée un fichier contenant **la matrice d’adjacence**.

Pour le lancer :

```
mpirun -np <nb_processus> ./main_mpi <chemin_fichier_matrice>
```

Exemple simple :

```
mpirun -np 4 ./main_mpi ../../DATA/PetitExemple.dot
```

Ou pour un plus gros fichier :

```
mpirun -np 9 ./main_mpi ../../DATA/exemplemassi.dot
```

> ⚠️ Les fichiers d’entrée sont dans le dossier `DATA` à la racine.

---

## 5. Sortie du programme

À la fin, le programme :

* calcule la matrice des plus courts chemins,
* rassemble le résultat sur le **rang 0**,
* écrit la matrice finale dans :

```
DATA/matrice_finale_sortie_de_floyd_warshal.txt
```

Le rang 0 affiche aussi le **temps total d’exécution**.

---

## 6. Remarques utiles

* Si le nombre de processus n’est **pas carré**, le programme adapte automatiquement la grille avec `MPI_Dims_create`.
* Si la taille n’est **pas divisible par la taille des blocs**, les cases “en trop” sont remplies avec **INF** (c’est juste du padding, ça ne gêne pas les calculs).
* La version utilise des **communications non bloquantes** pour accélérer la propagation des blocs (ça évite que tout le monde attende).

---

## 7. Exemple complet

```
make
mpirun -np 4 ./main_mpi ../../DATA/PetitExemple.dot
```

---

## 8. Nettoyage

Pour supprimer les `.o` et repartir propre :

```
make clean
```

---

