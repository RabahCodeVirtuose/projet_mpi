

# Version Parallèle de Floyd-Warshall (MPI)

## 1. Description rapide

Ce dossier contient **ma version parallèle** de l’algorithme de Floyd-Warshall, écrite en C++ avec **MPI**.
L’idée reste la même que pour la version séquentielle : on veut les plus courts chemins entre tous les sommets.
La différence est que **la grande matrice est divisée en blocs**, et chaque processus traite les blocs dont il est responsable.

À chaque itération `k`, le bloc pivot est transmis aux autres processus (`MPI_Bcast` ou `MPI_Ibcast`), ce qui leur permet de mettre à jour leurs blocs locaux.
C’est la méthode classique de **parallélisation 2D par blocs**, comme vue en TP.

📌 **Référence consultée**
Pendant la réalisation, j’ai aussi regardé un document externe qui explique une approche proche (découpage 2D, broadcasts, etc.).
Cela m’a aidé à organiser mon code.

> Asmita Gautam, *Parallel Floyd-Warshall Algorithm*, University at Buffalo, 2019.

---

## 2. Fichiers importants

* **`main_mpi.cpp`** – programme principal
* **`ParallelFWBlocks.cpp` / `.hpp`** – implémentation de Floyd-Warshall par blocs
* **`Distribution.cpp`** – répartition des blocs entre les processus
* **`Utils.cpp`** – affichage, écriture dans un fichier, etc.
* **`Makefile`** – compilation automatique

---

## 3. Compilation

Se placer dans le dossier :

```
FLOYD_PARALLEL_RABAH/
```

Puis compiler :

```
make
```

Un exécutable apparaît :

```
./main_mpi
```

---

## 4. Exécution

Le programme attend un fichier contenant **une matrice d’adjacence**.

Commande générale :

```
mpirun -np <nb_processus> ./main_mpi <chemin_fichier_matrice>
```

Exemples :

```
mpirun -np 4 ./main_mpi ../../DATA/PetitExemple.dot
```

```
mpirun -np 9 ./main_mpi ../../DATA/exemplemassi.dot
```

> Les fichiers d’entrée se trouvent dans le dossier `DATA`.

---

## 5. Sortie du programme

Le programme :

* calcule la matrice des plus courts chemins,
* rassemble tout sur le **rang 0**,
* écrit le résultat dans :

```
DATA/matrice_finale_sortie_de_floyd_warshal.txt
```

Le rang 0 affiche aussi le **temps d’exécution MPI**.

---

## 6. Remarques utiles

* Si le nombre de processus n’est **pas un carré**, la grille est adaptée automatiquement (`MPI_Dims_create`).
* Si la taille de la matrice n’est **pas un multiple de la taille des blocs**, les endroits “qui dépassent” sont remplis avec **INF** (padding).
* Des **communications non bloquantes** sont utilisées pour éviter que les processus attendent inutilement.

---

## 7. Exemple complet

```
make
mpirun -np 4 ./main_mpi ../../DATA/PetitExemple.dot
```

---

## 8. Nettoyage

Pour repartir propre :

```
make clean
```



