// BuildMatrixMPI.cpp (ou BuildDotMPI.cpp)
#define OMPI_SKIP_MPICXX 1

/**
 * @file BuildMatrixMPI.cpp
 * @brief Construction parallèle de la matrice de distances de Hamming entre séquences d'ARN
 *        et génération d'un graphe pondéré au format DOT.
 *
 * Ce programme :
 *   - lit un fichier FASTA contenant des séquences d'ARN (rang 0),
 *   - vérifie que toutes les séquences ont la même longueur,
 *   - diffuse les séquences à tous les processus MPI,
 *   - calcule en parallèle toutes les distances de Hamming entre les séquences,
 *   - rassemble la matrice de distances sur le rang 0,
 *   - génère un graphe pondéré non orienté au format DOT (poids = distance de Hamming),
 *   - n'écrit une arête que si la distance est strictement inférieure à un seuil ε.
 *
 * Le fichier DOT généré sert ensuite d'entrée à l'algorithme de Floyd–Warshall parallèle.
 */

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>

/**
 * @brief Lit un fichier FASTA "simple" et renvoie la liste des séquences.
 *
 * Format supposé :
 *   - les lignes qui commencent par '>' sont des en-têtes (identifiants),
 *   - les lignes suivantes (jusqu'au prochain '>') contiennent la séquence,
 *   - les éventuels retours à la ligne dans une séquence sont concaténés.
 *
 * @param filename Chemin du fichier FASTA à lire.
 * @return Un vecteur de chaînes, chaque élément étant une séquence complète.
 *
 * @throw std::runtime_error si le fichier ne peut pas être ouvert.
 */
static std::vector<std::string> readFasta(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("Impossible d'ouvrir " + filename);
    }

    std::vector<std::string> seqs;
    std::string line;
    std::string current;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            // Nouvelle séquence : on enregistre la précédente si elle existe
            if (!current.empty()) {
                seqs.push_back(current);
                current.clear();
            }
        } else {
            // Fragment de séquence : on concatène
            current += line;
        }
    }
    if (!current.empty()) {
        seqs.push_back(current);
    }

    return seqs;
}

/**
 * @brief Calcule la distance de Hamming entre deux séquences de même longueur.
 *
 * La distance de Hamming est le nombre de positions i telles que a[i] != b[i].
 *
 * @param a Pointeur vers la première séquence de longueur L.
 * @param b Pointeur vers la deuxième séquence de longueur L.
 * @param L Longueur des deux séquences.
 *
 * @return La distance de Hamming entre a et b.
 */
static int hamming(const char* a, const char* b, int L) {
    int d = 0;
    for (int i = 0; i < L; ++i) {
        if (a[i] != b[i]) ++d;
    }
    return d;
}

/**
 * @brief Écrit un graphe pondéré non orienté au format DOT à partir d'une matrice de distances.
 *
 * On génère un graphe de la forme :
 * @code
 *   graph graphe_pondere {
 *     A1 [label="0"];
 *     A2 [label="1"];
 *     ...
 *     A1 -- A2 [label="d", weight=d];
 *     ...
 *   }
 * @endcode
 *
 * Pour chaque paire (i, j) avec i < j, une arête est créée uniquement si
 * la distance d(i, j) est strictement inférieure à epsilon.
 *
 * @param filename Nom du fichier DOT à générer.
 * @param dist     Matrice des distances de taille n × n, stockée à plat (row-major).
 *                 L'élément (i, j) est à l'indice i * n + j.
 * @param n        Nombre de séquences / sommets du graphe.
 * @param epsilon  Seuil sur la distance de Hamming : on ne met une arête que si d < epsilon.
 *
 * @throw std::runtime_error si le fichier ne peut pas être ouvert en écriture.
 */
static void writeDotGraph(const std::string& filename,
                          const std::vector<int>& dist,
                          int n,
                          int epsilon)
{
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Impossible d'ouvrir " + filename + " en écriture");
    }

    out << "graph graphe_pondere {\n";
    out << "    node [shape=circle, style=filled, color=lightyellow, fontcolor=black];\n";
    out << "    edge [color=black, fontcolor=blue];\n\n";

    // Déclaration des sommets : A1, A2, ..., An
    // label = indice de la séquence (0,1,2,...)
    for (int i = 0; i < n; ++i) {
        out << "    A" << (i + 1) << " [label=\"" << i << "\"];\n";
    }
    out << "\n    // Les aretes avec poids (distance de Hamming < epsilon)\n";

    // Arêtes non orientées : i < j
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int d = dist[i * n + j];
            if (d < epsilon) {
                out << "    A" << (i + 1) << " -- A" << (j + 1)
                    << " [label=\"" << d << "\", weight=" << d << "];\n";
            }
        }
    }

    out << "}\n";
}

/**
 * @brief Programme principal MPI : construction de la matrice de distances et du graphe DOT.
 *
 * Étapes principales :
 *   - rang 0 lit un fichier FASTA et vérifie que toutes les séquences ont la même longueur,
 *   - n (nombre de séquences) et L (longueur des séquences) sont diffusés à tous,
 *   - les séquences sont diffusées à tous les rangs sous forme de tableau contigu,
 *   - chaque rang calcule un sous-ensemble de lignes de la matrice des distances
 *     de Hamming (d(i, j) pour ses lignes i),
 *   - le rang 0 rassemble les sous-matrices pour reconstruire la matrice n × n complète,
 *   - le rang 0 écrit un fichier DOT pondéré (utilisé ensuite par l'algorithme de Floyd–Warshall),
 *   - le temps total (calcul + rassemblement) est mesuré avec MPI_Wtime().
 *
 * Usage typique :
 * @code
 *   mpirun -np <nb_processus> ./build_matrix_mpi dataset_500seq.fa
 * @endcode
 *
 * @param argc Nombre d'arguments de la ligne de commande.
 * @param argv Tableau d'arguments (argv[1] doit être le fichier FASTA).
 *
 * @return 0 en cas de succès, une valeur non nulle si une erreur survient.
 */
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

     // Ici je récupère le fichier FASTA en argument.
    // C'est lui qui contient toutes les séquences 
    const std::string fastaFile = argv[1];
        // Et je fixe le nom du fichier DOT de sortie, que je vais donner ensuite à Floyd.
    const std::string dotFile   = "../../DATA/Resulat_sequence_by_premier_algo.dot";

   // Paramètre epsilon (voir énoncé, genre epsilon = 70).
    // Si la distance de Hamming entre deux séquences est < epsilon,
    // alors je crée une arête entre elles dans le graphe.    const int epsilon = 70;
    const int epsilon = 70;

    int n = 0;       // nombre de séquences
    int L = 0;       // longueur d'une séquence
     // Je vais stocker toutes les séquences à la suite dans un seul gros tableau de chars,
    // de taille n * L. La séquence i commence à l'offset i * L.
    std::vector<char> allSeqs;  // tableau contigu de taille n * L

    // ----------------------------------------------------------
    // Rang 0 : lecture du FASTA + contrôle des longueurs
    // ----------------------------------------------------------
                // Je fais toute la lecture disque uniquement sur le rang 0

    if (rank == 0) {
        try {
            // Je lis le FASTA et je récupère un vecteur de chaînes,
            // chaque string correspond à une séquence.
            std::vector<std::string> seqs = readFasta(fastaFile);
            n = (int)seqs.size();
            if (n == 0) {
                throw std::runtime_error("Aucune sequence lue dans " + fastaFile);
            }
                        // Je prends la longueur de la première séquence comme référence
            L = (int)seqs[0].size();
         

            std::cout << "\n\nLecture FASTA: n = " << n
                      << ", longueur L = " << L << "\n";

            // Copie dans un tableau contigu n * L
              // Maintenant je recopie tout dans un grand tableau contigu n * L.
            // Comme ça après, chaque processus peut calculer les distances
            // juste avec un &allSeqs[i * L].
            allSeqs.resize(n * L);
            for (int i = 0; i < n; ++i) {
                std::copy(seqs[i].begin(), seqs[i].end(),
                          allSeqs.begin() + i * L);
            }
        } catch (const std::exception& e) {
            // Si je tombe sur un problème (fichier introuvable, séquences pas de même taille, etc.),
            // j'affiche l'erreur sur le rang 0 et j'arrête tout le monde proprement
            std::cerr << "Erreur (rang 0) : " << e.what() << "\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // ----------------------------------------------------------
    // Diffuser n et L à tous les processus
    // ----------------------------------------------------------
      // Ici j’envoie à tout le monde le nombre de séquences (n)
    // et la longueur d’une séquence (L), que seul le rang 0 connaît au début.
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);



    // Tous les rangs allouent le tableau allSeqs,
    // sauf le rang 0 qui l’a déjà rempli plus haut.
    if (rank != 0) {
        allSeqs.resize(n * L);
    }

    // Ici je diffuse toutes les séquences en une seule fois :
    // le rang 0 envoie son gros tableau n*L vers tout le monde.
    MPI_Bcast(allSeqs.data(), n * L, MPI_CHAR, 0, MPI_COMM_WORLD);

    // ----------------------------------------------------------
    // Mesure du temps MPI (calcul + rassemblement)
    // ----------------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    // ----------------------------------------------------------
    // Chaque rang calcule un sous-ensemble de lignes de la matrice
    // ----------------------------------------------------------
        // Je découpe les n séquences en "tranches" entre les processus.
    // chunk = nombre de lignes (séquences) par processus, en arrondissant vers le haut.
    int chunk = (n + size - 1) / size;  // ceil(n / size)
    int start = rank * chunk; // première ligne (incluse) pour ce rang
    int end   = std::min(n, start + chunk);  // dernière ligne (exclue) pour ce rang on fait ça pour éviter de dépasser n

     // Nombre de lignes réellement calculées par ce rang
    int localRows = std::max(0, end - start);


    
    // Chaque rang va calculer localRows lignes de la matrice des distances,
    // donc au total localRows * n entiers. 
    std::vector<int> localDist(localRows * n);

    for (int i = start; i < end; ++i) {
          // seq_i = séquence de la ligne i
        const char* seq_i = &allSeqs[i * L];
                // rowOffset = où commence la ligne i dans le tableau localDist

        int rowOffset = (i - start) * n;

        for (int j = 0; j < n; ++j) {
            const char* seq_j = &allSeqs[j * L];
              // Distance de Hamming entre i et j.
            // Par convention je mets 0 sur la diagonale (i == j).
            int d = (i == j) ? 0 : hamming(seq_i, seq_j, L);
            localDist[rowOffset + j] = d;
        }
    }

    // ----------------------------------------------------------
    // Rassemblement de la matrice des distances sur le rang 0
    // ----------------------------------------------------------
    std::vector<int> fullDist;
    if (rank == 0) {
                // Rang 0 alloue la matrice complète n*n.

        fullDist.resize(n * n);
    }
    // Je repasse sur tous les rangs r pour reconstruire la matrice globale
    // dans le bon ordre des lignes.
    for (int r = 0; r < size; ++r) {
        int r_start = r * chunk;
        int r_end   = std::min(n, r_start + chunk);
        int r_rows  = std::max(0, r_end - r_start);
        int count   = r_rows * n;

        if (r_rows == 0) continue; // ce rang n'a pas de lignes

        if (rank == 0) {
            // Côté rang 0 : je récupère ce que le rang r a calculé.
            if (r == 0) {
                 // Pour r == 0, je copie directement mon propre localDist
                // au bon endroit dans fullDist.

                std::copy(localDist.begin(), localDist.end(),
                          fullDist.begin() + r_start * n);
            } else {
                // Recevoir le bloc du rang r
 // Pour les autres rangs, je reçois leurs lignes dans un buffer,
                // puis je copie ça dans la zone correspondant à r_start, r_end.

                std::vector<int> buf(count);
                MPI_Recv(buf.data(), count, MPI_INT, r, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::copy(buf.begin(), buf.end(),
                          fullDist.begin() + r_start * n);
            }
        } else if (rank == r) {
            // Rang r (/= 0) : envoie son bloc local au rang 0
            MPI_Send(localDist.data(), count, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();

    // ----------------------------------------------------------
    // Rang 0 : écriture du fichier .dot + affichage du temps
    // ----------------------------------------------------------
    if (rank == 0) {
        std::cout << "\n\n>>> Temps total calcul distances + rassemblement = "
                  << (t1 - t0) * 1000 << " millisecondes\n\n";

        try {
            writeDotGraph(dotFile, fullDist, n, epsilon);
            std::cout << "Graphe .dot ecrit dans " << dotFile << "\n";
        } catch (const std::exception& e) {
            std::cerr << "Erreur d'ecriture du fichier .dot : " << e.what() << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}
