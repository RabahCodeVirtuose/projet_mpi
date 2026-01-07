# Pipeline complet : séquences ARN → Floyd–Warshall → PAM (hybride)

Ce dossier correspond au **pipeline complet** de l’application :

1. Lecture d’un fichier FASTA contenant des séquences d’ARN.
2. Calcul des distances **Needleman–Wunsch** et génération d’un graphe pondéré au format **DOT** (OpenMP).
3. Calcul parallèle des plus courts chemins avec l’algorithme de **Floyd–Warshall** (version blocs, MPI + OpenMP).
4. Clustering des sommets avec l’algorithme **PAM** en mode hybride (MPI + OpenMP).

Le but de ce dossier est de pouvoir **enchaîner automatiquement toutes les étapes**, du FASTA jusqu’au fichier de résultats PAM.

---

## 1. Organisation

Le dossier contient trois sous-parties principales :

- `NEEDLMAN_SEQUENCE_to_DOT` : conversion FASTA → distances Needleman–Wunsch → graphe DOT (OpenMP).
- `FLOYD` : algorithme de Floyd–Warshall parallèle (MPI + OpenMP, par blocs).
- `PAM_MPI` : algorithme PAM hybride appliqué à la matrice de distances.

Un `Makefile` à la racine permet de **compiler** ces trois modules et de **lancer le pipeline complet**.

---

## 2. Compilation globale

Depuis `FULL_APP_PIPE_LINE/`, vous pouvez compiler tous les modules avec :

```bash
make
````

Cela appelle les `Makefile` présents dans :

* `NEEDLMAN_SEQUENCE_to_DOT/`
* `FLOYD/`
* `PAM_MPI/`

et génère les exécutables suivants :

* `NEEDLMAN_SEQUENCE_to_DOT/build_matrix_needleman`
* `FLOYD/main_mpi`
* `PAM_MPI/pam_mpi`

---

## 3. Variables importantes

Dans le `Makefile` racine, quelques variables peuvent être ajustées si besoin :

* `FASTA` : chemin du fichier FASTA (jeu de séquences ARN).
* `DOT_FILE` : fichier DOT généré à partir des distances Needleman–Wunsch.
* `DIST_FILE` : fichier contenant la matrice de distances finale (sortie de Floyd–Warshall).
* `PAM_OUT` : fichier de sortie pour les résultats de PAM.
* `OMP_SEQ`, `OMP_FLOYD` : nombre de threads OpenMP pour Needleman et Floyd.
* `NP_FLOYD`, `NP_PAM` : nombre de processus MPI pour Floyd et PAM.

Par défaut, ces variables sont définies au début du `Makefile`, mais vous pouvez les surcharger à l’appel (voir plus bas).

---

## 4. Exécution du pipeline complet

Une fois la compilation effectuée (`make`), vous pouvez lancer **toute la chaîne** (FASTA → DOT → Floyd → PAM) avec :

```bash
make run
```

Le `Makefile` exécute alors successivement :

1. `build_matrix_needleman` sur le fichier FASTA
2. `main_mpi` sur le fichier DOT généré
3. `pam_mpi` sur la matrice de distances calculée par Floyd–Warshall

Si vous souhaitez modifier le nombre de processus MPI utilisés pour chaque étape, vous pouvez faire par exemple :

```bash
make run OMP_SEQ=12 NP_FLOYD=4 OMP_FLOYD=3 NP_PAM=4 OMP_PAM=3
```

Les valeurs par défaut sont fixées dans le `Makefile` via :

```make
OMP_SEQ   ?= 12
NP_FLOYD  ?= 4
OMP_FLOYD ?= 3
NP_PAM    ?= 4
OMP_PAM   ?= 3 
```

---

## 5. Nettoyage

Pour nettoyer les fichiers objets et exécutables dans tous les sous-dossiers, vous pouvez utiliser :

```bash
make clean
```

Cela appelle la cible `clean` dans :

* `NEEDLMAN_SEQUENCE_to_DOT/`
* `FLOYD/`
* `PAM_MPI/`

---

## 6. Test des parties séparément

Si vous souhaitez tester chaque étape **séparément** (par exemple seulement Floyd ou seulement PAM), chaque sous-dossier contient :

* son propre `Makefile`,
* un `Readme.md` qui explique comment **compiler et exécuter cette partie de manière isolée**.

Il suffit d’ouvrir le dossier correspondant et de suivre le `Readme.md` local.
