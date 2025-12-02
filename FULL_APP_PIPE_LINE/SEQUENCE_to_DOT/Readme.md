Voilà un README simple pour le dossier **SEQUENCE_to_DOT** (BuildMatrixMPI.cpp).

---

# Construction de la matrice de distances / fichier DOT (MPI)

## 1. Description rapide

Ce programme prend un fichier **FASTA** avec des séquences d’ARN et produit un **graphe pondéré au format DOT**.

En résumé, il fait :

1. lecture des séquences dans le fichier FASTA (sur le rang 0),
2. vérification que toutes les séquences ont la même longueur,
3. diffusion des séquences à tous les processus MPI,
4. calcul en parallèle de toutes les **distances de Hamming** entre les séquences,
5. rassemblement de la matrice des distances sur le rang 0,
6. génération d’un fichier DOT avec un graphe non orienté, où
   le poids de l’arête entre deux sommets = distance de Hamming entre les deux séquences.

Une arête n’est créée que si la distance est **strictement inférieure** à un seuil `ε` (dans mon code, `epsilon = 70`).

Le fichier DOT généré sert ensuite d’entrée à l’algorithme de Floyd–Warshall parallèle.

---

## 2. Compilation

Depuis le dossier `SEQUENCE_to_DOT` :

```bash
cd SEQUENCE_to_DOT
make
```

Le `Makefile` génère l’exécutable :

```bash
./build_dot
```

---

## 3. Fichier d’entrée (FASTA)

Le fichier FASTA est supposé être “simple” et toutes les séquences ont la même taille 

Exemples de fichiers :

* `../../DATA/dataset_500seq.fa`
* `../../DATA/dataset_2000seq.fa`

---

## 4. Exécution

Pour lancer le programme avec MPI :

```bash
mpirun -np 4 ./build_dot ../../DATA/dataset_2000seq.fa
```

* `-np 4` : nombre de processus MPI,
* dernier argument : chemin vers le fichier FASTA.

Le programme :

* lit le FASTA sur le rang 0,
* calcule la matrice de distances de Hamming en parallèle,
* mesure le temps total (calcul + rassemblement),
* écrit le graphe DOT dans :

```bash
../../DATA/Resulat_sequence_by_premier_algo.dot
```

C’est ce fichier DOT qui sera utilisé après par **Floyd–Warshall**.

---

## 5. Paramètre epsilon

Dans le code, le seuil est fixé à :

```cpp
const int epsilon = 70;
```

Pour chaque paire de séquences `(i, j)` :

* on calcule `d = distance_de_Hamming(i, j)`
* si `d < epsilon`, on crée une arête `Ai -- Aj` dans le fichier DOT,
* sinon, aucune arête n’est écrite entre ces deux sommets.

---

## 6. Nettoyage

Pour supprimer les fichiers objets / recompiler propre :

```bash
make clean
```

---
