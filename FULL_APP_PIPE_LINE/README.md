# Pipeline complet : séquences ARN → Floyd–Warshall → PAM (MPI)

Ce dossier correspond au **pipeline complet** de l’application :

1. Lecture d’un fichier FASTA contenant des séquences d’ARN.
2. Calcul des distances de Hamming entre toutes les séquences et génération d’un graphe pondéré au format **DOT**.
3. Calcul parallèle des plus courts chemins avec l’algorithme de **Floyd–Warshall** (version blocs, MPI).
4. Clustering des sommets avec l’algorithme **PAM** en parallèle.

Le but de ce dossier est de pouvoir **enchaîner automatiquement toutes les étapes**, du FASTA jusqu’au fichier de résultats PAM.

---

## 1. Organisation

Le dossier contient trois sous-parties principales :

- `SEQUENCE_to_DOT` : conversion FASTA → distances de Hamming → graphe DOT.
- `Floyd_parallel_by_rabah_toubal` : algorithme de Floyd–Warshall parallèle (MPI, par blocs).
- `PAM_MPI` : algorithme PAM parallèle appliqué à la matrice de distances.

Un `Makefile` à la racine permet de **compiler** ces trois modules et de **lancer le pipeline complet**.

---

## 2. Compilation globale

Depuis `FULL_APP_PIPE_LINE/`, vous pouvez compiler tous les modules avec :

```bash
make
````

Cela appelle les `Makefile` présents dans :

* `SEQUENCE_to_DOT/`
* `Floyd_parallel_by_rabah_toubal/`
* `PAM_MPI/`

et génère les exécutables suivants :

* `SEQUENCE_to_DOT/build_dot`
* `Floyd_parallel_by_rabah_toubal/main_mpi`
* `PAM_MPI/pam_mpi`

---

## 3. Variables importantes

Dans le `Makefile` racine, quelques variables peuvent être ajustées si besoin :

* `FASTA` : chemin du fichier FASTA (jeu de séquences ARN).
* `DOT_FILE` : fichier DOT généré à partir des distances de Hamming.
* `DIST_FILE` : fichier contenant la matrice de distances finale (sortie de Floyd–Warshall).
* `PAM_OUT` : fichier de sortie pour les résultats de PAM.
* `NP_SEQ`, `NP_FLOYD`, `NP_PAM` : nombre de processus MPI utilisés pour chaque étape.

Par défaut, ces variables sont définies au début du `Makefile`, mais vous pouvez les surcharger à l’appel (voir plus bas).

---

## 4. Exécution du pipeline complet

Une fois la compilation effectuée (`make`), vous pouvez lancer **toute la chaîne** (FASTA → DOT → Floyd → PAM) avec :

```bash
make run
```

Le `Makefile` exécute alors successivement :

1. `build_dot` sur le fichier FASTA
2. `main_mpi` sur le fichier DOT généré
3. `pam_mpi` sur la matrice de distances calculée par Floyd–Warshall

Si vous souhaitez modifier le nombre de processus MPI utilisés pour chaque étape, vous pouvez faire par exemple :

```bash
make run NP_SEQ=4 NP_FLOYD=6 NP_PAM=6
```

Les valeurs par défaut sont fixées dans le `Makefile` via :

```make
NP_SEQ   ?= 6
NP_FLOYD ?= 6
NP_PAM   ?= 6
```

---

## 5. Nettoyage

Pour nettoyer les fichiers objets et exécutables dans tous les sous-dossiers, vous pouvez utiliser :

```bash
make clean
```

Cela appelle la cible `clean` dans :

* `SEQUENCE_to_DOT/`
* `Floyd_parallel_by_rabah_toubal/`
* `PAM_MPI/`

---

## 6. Test des parties séparément

Si vous souhaitez tester chaque étape **séparément** (par exemple seulement Floyd ou seulement PAM), chaque sous-dossier contient :

* son propre `Makefile`,
* un `Readme.md` qui explique comment **compiler et exécuter cette partie de manière isolée**.

Il suffit d’ouvrir le dossier correspondant et de suivre le `Readme.md` local.

