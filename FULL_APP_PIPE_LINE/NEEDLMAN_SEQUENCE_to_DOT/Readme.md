# Construction de la matrice de distances / fichier DOT (OpenMP)

## 1. Description rapide

Ce programme prend un fichier **FASTA** avec des séquences d’ARN et produit un **graphe pondéré au format DOT**.

En résumé, il fait :

1. lecture des séquences dans le fichier FASTA,
2. vérification que toutes les séquences ont la même longueur,
3. calcul en parallèle (OpenMP) de toutes les distances **Needleman–Wunsch**,
4. génération d’un fichier DOT avec un graphe non orienté, où
   le poids de l’arête entre deux sommets = distance dérivée du score Needleman–Wunsch.

Une arête n’est créée que si la distance est **strictement inférieure** à un seuil `ε` (dans mon code, `epsilon = 70`).

Le fichier DOT généré sert ensuite d’entrée à l’algorithme de Floyd–Warshall parallèle.

---

## 2. Compilation

Depuis le dossier `NEEDLMAN_SEQUENCE_to_DOT` :

```bash
cd NEEDLMAN_SEQUENCE_to_DOT
make
```

Le `Makefile` génère l’exécutable :

```bash
./build_matrix_needleman
```

---

## 3. Fichier d’entrée (FASTA)

Le fichier FASTA est supposé être “simple” et toutes les séquences ont la même taille 

Exemples de fichiers :

* `../../DATA/dataset_500seq.fa`
* `../../DATA/dataset_2000seq.fa`

---

## 4. Exécution

Pour lancer le programme avec OpenMP :

```bash
OMP_NUM_THREADS=12 ./build_matrix_needleman ../DATA/dataset_2000seq.fa
```

* `OMP_NUM_THREADS` : nombre de threads OpenMP,
* dernier argument : chemin vers le fichier FASTA.

Le programme :

* lit le FASTA,
* calcule la matrice de distances Needleman–Wunsch en parallèle,
* mesure le temps total (calcul),
* écrit le graphe DOT dans :

```bash
../../DATA/Resulat_sequence_by_premier_algo.dot
```

C’est ce fichier DOT qui sera utilisé après par **Floyd–Warshall**.

---

## 5. Benchmark (1 thread vs plusieurs)

Le Makefile contient une cible `benchmark` qui lance plusieurs executions
avec des nombres de threads differents, pour comparer le temps sequentiel
(1 thread) et parallele (plusieurs threads).

Commande :

```bash
make benchmark
```

Par defaut, la liste des threads testes est : `1 2 4 8 12`.
Vous pouvez modifier cette liste directement dans le Makefile si besoin.

---

## 6. Paramètre epsilon

Dans le code, le seuil est fixé à :

```cpp
const int epsilon = 70;
```

Pour chaque paire de séquences `(i, j)` :

* on calcule un score Needleman–Wunsch,
* on convertit en distance avec `d = (L - score) / 2`,
* si `d < epsilon`, on crée une arête `Ai -- Aj` dans le fichier DOT,
* sinon, aucune arête n’est écrite entre ces deux sommets.

---

## 7. Nettoyage

Pour supprimer les fichiers objets / recompiler propre :

```bash
make clean
```

---
