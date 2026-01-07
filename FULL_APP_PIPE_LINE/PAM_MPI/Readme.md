
# Version hybride MPI + OpenMP de l’algorithme PAM (k-médoides)

## 1. Description rapide

Ce dossier contient **ma version hybride** de l’algorithme **PAM (k-médoides)**.

* En entrée : une **matrice de distances** `n × n` (par exemple celle calculée par Floyd-Warshall).
* En sortie : les **médoines** choisis, l’affectation de chaque sommet à un cluster, et le coût total.

Le calcul du coût est réparti entre les processus MPI : chacun traite une partie des lignes de la matrice, puis on fait une réduction pour avoir le coût global.
À l’intérieur de chaque processus, la boucle locale est parallélisée avec **OpenMP**.

Par défaut, j’utilise `k = 4` groupes.

---

## 2. Fichiers importants

* `main.cpp` → programme principal (lecture du fichier de distances, appel de `runPAM_MPI`)
* `PAM.cpp / PAM.hpp` → implémentation de l’algorithme PAM en version hybride (MPI + OpenMP)
* `Utils.cpp / Utils.hpp` → fonctions utilitaires (lecture de la matrice, écriture du résultat)
* `Makefile` → compilation automatique

---

## 3. Compilation

Depuis le dossier `PAM_MPI` :

```bash
cd PAM_MPI
make
```

Ça génère l’exécutable :

```bash
./pam_mpi
```

---

## 4. Format du fichier d’entrée

Le programme attend un fichier texte contenant une **matrice de distances carrée**, au format :

```txt
n m
d00 d01 d02 ... d0(n-1)
d10 d11 d12 ... d1(n-1)
...
d(n-1)0 ...          d(n-1)(n-1)
```

Dans mon projet, ce fichier est en général :

```bash
../../DATA/matrice_finale_sortie_de_floyd_warshal.txt
```

…qui est justement la sortie de l’algorithme de Floyd parallèle.

---

## 5. Exécution

Commande typique :

```bash
OMP_NUM_THREADS=2 mpirun -np 6 ./pam_mpi ../../DATA/matrice_finale_sortie_de_floyd_warshal.txt
```

* `OMP_NUM_THREADS` → nombre de threads OpenMP par processus
* `-np 6` → nombre de processus MPI
* dernier argument → chemin du fichier de distances

---

## 6. Sortie du programme

Le rang 0 :

* affiche le **coût final** et la liste des **médoides** choisis,
* écrit le résultat détaillé dans :

```bash
../../DATA/resultat_pam_parallel.txt
```

Ce fichier contient :

* la taille `n`,
* la valeur de `k`,
* le coût total,
* la liste des médoines,
* puis, pour chaque sommet :
  `sommet  cluster  medoid  dist`.

---

## 7. Nettoyage

Pour supprimer les fichiers objets et l’exécutable :

```bash
make clean
```

---
