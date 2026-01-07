# ============================
# Makefile RACINE - PIPELINE COMPLET
# ============================

# Compilateur MPI
CXX      = mpic++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pedantic

# Dossiers
SEQ_DIR    = NEEDLMAN_SEQUENCE_to_DOT
FLOYD_DIR  = FLOYD
PAM_DIR    = PAM_MPI

# Executables dans chaque dossier
SEQ_EXE    = build_matrix_needleman
FLOYD_EXE  = main_mpi
PAM_EXE    = pam_mpi

# Fichiers de données (chemins vus depuis la RACINE)
FASTA      = ../../DATA/dataset_2000seq.fa         
DOT_FILE   = ../../DATA/Resulat_sequence_by_premier_algo.dot
DIST_FILE  = ../../DATA/matrice_finale_sortie_de_floyd_warshal.txt
PAM_OUT    = ../../DATA/resultat_pam_parallel.txt

# ============================
# CONFIGURATION PARALLÉLISME
# ============================

# Nombre de processus MPI (modifiable à l'appel)
NP_FLOYD ?= 4
NP_PAM   ?= 4

# Nombre de threads OpenMP par processus
OMP_SEQ   ?= 12    # Needleman : OpenMP pur (pas de MPI)
OMP_FLOYD ?= 3    # Floyd : MPI+OpenMP hybride (6 proc × 2 threads = 12 threads)
OMP_PAM   ?= 3    # PAM   : MPI+OpenMP hybride (6 proc × 2 threads = 12 threads)

# ============================
# 1) Compilation
# ============================

all: seq floyd pam

seq:
	@echo "=== Compilation Needleman-Wunsch ==="
	$(MAKE) -C $(SEQ_DIR)

floyd:
	@echo "=== Compilation Floyd-Warshall ==="
	$(MAKE) -C $(FLOYD_DIR)

pam:
	@echo "=== Compilation PAM ==="
	$(MAKE) -C $(PAM_DIR)

# ============================
# 2) Exécution du pipeline
# ============================

run: all
	@echo
	@echo "=========================================="
	@echo "     PIPELINE COMPLET - DÉMARRAGE"
	@echo "=========================================="
	@echo "Configuration :"
	@echo "  - Needleman : $(OMP_SEQ) threads OpenMP"
	@echo "  - Floyd     : $(NP_FLOYD) processus MPI × $(OMP_FLOYD) threads OpenMP = $$(( $(NP_FLOYD) * $(OMP_FLOYD) )) threads"
	@echo "  - PAM       : $(NP_PAM) processus MPI × $(OMP_PAM) threads OpenMP = $$(( $(NP_PAM) * $(OMP_PAM) )) threads"
	@echo "=========================================="
	@echo
	@echo
	@echo "=== ÉTAPE 1/3 : Construction du graphe (Needleman-Wunsch) ==="
	@echo "Input  : $(FASTA)"
	@echo "Output : $(DOT_FILE)"
	@echo
	cd $(SEQ_DIR) && OMP_NUM_THREADS=$(OMP_SEQ) ./$(SEQ_EXE) $(FASTA)
	@echo
	@echo
	@echo "=== ÉTAPE 2/3 : Calcul des plus courts chemins (Floyd-Warshall) ==="
	@echo "Input  : $(DOT_FILE)"
	@echo "Output : $(DIST_FILE)"
	@echo
	cd $(FLOYD_DIR) && OMP_NUM_THREADS=$(OMP_FLOYD) mpirun -np $(NP_FLOYD) ./$(FLOYD_EXE) $(DOT_FILE)
	@echo
	@echo
	@echo "=== ÉTAPE 3/3 : Clustering (PAM) ==="
	@echo "Input  : $(DIST_FILE)"
	@echo "Output : $(PAM_OUT)"
	@echo
	cd $(PAM_DIR) && OMP_NUM_THREADS=$(OMP_PAM) mpirun -np $(NP_PAM) ./$(PAM_EXE) $(DIST_FILE)
	@echo
	@echo
	@echo "=========================================="
	@echo "     PIPELINE TERMINÉ AVEC SUCCÈS"
	@echo "=========================================="
	@echo "Fichiers générés :"
	@echo "  1. Graphe     : $(DOT_FILE)"
	@echo "  2. Distances  : $(DIST_FILE)"
	@echo "  3. Clusters   : $(PAM_OUT)"
	@echo "=========================================="

# ============================
# 3) Exécution étape par étape
# ============================

run-seq:
	@echo "=== Exécution : Needleman-Wunsch ==="
	cd $(SEQ_DIR) && OMP_NUM_THREADS=$(OMP_SEQ) ./$(SEQ_EXE) $(FASTA)

run-floyd:
	@echo "=== Exécution : Floyd-Warshall ==="
	cd $(FLOYD_DIR) && OMP_NUM_THREADS=$(OMP_FLOYD) mpirun -np $(NP_FLOYD) ./$(FLOYD_EXE) $(DOT_FILE)

run-pam:
	@echo "=== Exécution : PAM ==="
	cd $(PAM_DIR) && OMP_NUM_THREADS=$(OMP_PAM) mpirun -np $(NP_PAM) ./$(PAM_EXE) $(DIST_FILE)

# ============================
# 4) Tests de scalabilité
# ============================

# Test de scalabilité complète
benchmark:
	@echo "=== BENCHMARK DE SCALABILITÉ ==="
	@echo
	@echo "Test 1 : Configuration minimale (1×1)"
	$(MAKE) run NP_FLOYD=1 OMP_FLOYD=1 NP_PAM=1 OMP_PAM=1 OMP_SEQ=1
	@echo
	@echo "Test 2 : Configuration moyenne (4×2)"
	$(MAKE) run NP_FLOYD=4 OMP_FLOYD=2 NP_PAM=4 OMP_PAM=2 OMP_SEQ=8
	@echo
	@echo "Test 3 : Configuration maximale (6×2)"
	$(MAKE) run NP_FLOYD=6 OMP_FLOYD=2 NP_PAM=6 OMP_PAM=2 OMP_SEQ=12

# ============================
# 5) Nettoyage
# ============================

clean:
	@echo "=== Nettoyage des sous-projets ==="
	$(MAKE) -C $(SEQ_DIR) clean
	$(MAKE) -C $(FLOYD_DIR) clean
	$(MAKE) -C $(PAM_DIR) clean

clean-data:
	@echo "=== Nettoyage des fichiers de données ==="
	rm -f $(DOT_FILE) $(DIST_FILE) $(PAM_OUT)

clean-all: clean clean-data

# ============================
# 6) Aide
# ============================

help:
	@echo "=========================================="
	@echo "  Makefile Pipeline - Aide"
	@echo "=========================================="
	@echo
	@echo "Cibles principales :"
	@echo "  make              - Compile tous les projets"
	@echo "  make run          - Exécute le pipeline complet"
	@echo "  make clean        - Nettoie les executables"
	@echo "  make clean-all    - Nettoie tout (executables + données)"
	@echo
	@echo "Cibles individuelles :"
	@echo "  make run-seq      - Exécute seulement Needleman"
	@echo "  make run-floyd    - Exécute seulement Floyd"
	@echo "  make run-pam      - Exécute seulement PAM"
	@echo
	@echo "Tests :"
	@echo "  make benchmark    - Test de scalabilité"
	@echo
	@echo "Configuration (variables) :"
	@echo "  OMP_SEQ           - Threads OpenMP pour Needleman (défaut: 12)"
	@echo "  NP_FLOYD          - Processus MPI pour Floyd (défaut: 6)"
	@echo "  OMP_FLOYD         - Threads OpenMP/proc pour Floyd (défaut: 2)"
	@echo "  NP_PAM            - Processus MPI pour PAM (défaut: 6)"
	@echo "  OMP_PAM           - Threads OpenMP/proc pour PAM (défaut: 2)"
	@echo
	@echo "Exemples :"
	@echo "  make run"
	@echo "  make run NP_FLOYD=4 OMP_FLOYD=3 NP_PAM=4 OMP_PAM=3"
	@echo "  make run OMP_SEQ=16 NP_FLOYD=8 OMP_FLOYD=2"
	@echo
	@echo "=========================================="

.PHONY: all seq floyd pam run run-seq run-floyd run-pam benchmark clean clean-data clean-all help