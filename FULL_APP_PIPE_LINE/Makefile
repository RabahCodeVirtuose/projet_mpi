# ============================
# Makefile RACINE - PIPELINE
# ============================

# Compilateur MPI
CXX      = mpic++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pedantic

# Dossiers
SEQ_DIR    = SEQUENCE_to_DOT
FLOYD_DIR  = Floyd_parallel_by_rabah_toubal
PAM_DIR    = PAM_MPI

# Executables dans chaque dossier
SEQ_EXE    = build_dot
FLOYD_EXE  = main_mpi
PAM_EXE    = pam_mpi

# Fichiers de données (chemins vus depuis la RACINE)
FASTA      = ../../DATA/dataset_500seq.fa         
DOT_FILE   = ../../DATA/Resulat_sequence_by_premier_algo.dot
DIST_FILE  = ../../DATA/matrice_finale_sortie_de_floyd_warshal.txt
PAM_OUT    = ../../DATA/resultat_pam_parallel.txt

# Nombre de processus MPI pour chaque étape (modifiable à l'appel)
NP_SEQ   ?= 6
NP_FLOYD ?= 6
NP_PAM   ?= 6

# Nombre de threads OpenMP (Needleman) pour l'étape SEQ
OMP_SEQ  ?= 1
OMP_FLOYD ?= 12
OMP_PAM   ?= 12

# ============================
# 1) Compilation
# ============================

# "all" compile les 3 projets en appelant leurs Makefile
all: seq floyd pam

seq:
	$(MAKE) -C $(SEQ_DIR)

floyd:
	$(MAKE) -C $(FLOYD_DIR)

pam:
	$(MAKE) -C $(PAM_DIR)

# ============================
# 2) Exécution du pipeline
# ============================

# Pipeline complet : FASTA -> DOT -> Floyd -> PAM
run: all
	@echo
	@echo
	@echo
	@echo
	@echo "=== Étape 1 : FASTA -> DOT (build_dot) ==="
	@echo
	cd $(SEQ_DIR) && OMP_NUM_THREADS=$(OMP_SEQ) mpirun -np $(NP_SEQ) ./$(SEQ_EXE) $(FASTA) $(OMP_SEQ)
	@echo
	@echo
	@echo
	@echo
	@echo "=== Étape 2 : DOT -> matrice (Floyd) ==="
	@echo
	cd $(FLOYD_DIR) && OMP_NUM_THREADS=$(OMP_FLOYD) mpirun -np $(NP_FLOYD) ./$(FLOYD_EXE) $(DOT_FILE)
	@echo
	@echo
	@echo
	@echo
	@echo "=== Étape 3 : matrice -> PAM (clustering) ==="
	@echo
	cd $(PAM_DIR) && OMP_NUM_THREADS=$(OMP_PAM) mpirun -np $(NP_PAM) ./$(PAM_EXE) $(DIST_FILE)
	@echo
	@echo
	@echo
	@echo
	@echo "=== Pipeline terminé ==="
	@echo "  DOT      : $(DOT_FILE)"
	@echo "  Matrice  : $(DIST_FILE)"
	@echo "  PAM out  : $(PAM_OUT)"


# ============================
# 3) Nettoyage global
# ============================

clean:
	$(MAKE) -C $(SEQ_DIR) clean
	$(MAKE) -C $(FLOYD_DIR) clean
	$(MAKE) -C $(PAM_DIR) clean

.PHONY: all seq floyd pam run run-seq run-floyd run-pam clean
