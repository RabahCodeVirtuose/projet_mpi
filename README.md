# projet_mpi
Projet de Rabah TOUBAL

# FloydParallel — Version MPI du Floyd–Warshall

Structure :
- Lecture du graphe : ForGraphMPI.cpp
- Distribution des blocs : Distribution.cpp
- Algorithme parallèle : ParallelFW.cpp
- Main MPI : main_mpi.cpp

Lancer :
séquentiel : ./main ../Exemple2.dot
parallel : mpirun -np 6 ./main_mpi ../Exemple2.dot

Étapes :
1. Rank 0 lit le graphe (.dot).
2. Diffusion de la matrice brute.
3. Création des blocs b×b.
4. Répartition bloc-cyclique sur une grille 2×3.
5. Trois phases par itération :
   - Phase A : pivot
   - Phase B : ligne k + colonne k
   - Phase C : mise à jour restante
6. Recomposition finale sur rank 0.
