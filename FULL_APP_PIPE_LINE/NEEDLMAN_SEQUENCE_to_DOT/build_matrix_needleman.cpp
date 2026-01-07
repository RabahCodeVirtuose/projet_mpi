/**
 * @file build_matrix_needleman.cpp
 * @brief Construction parallèle OpenMP de la matrice de distances Needleman-Wunsch
 * 
 * Stratégie de parallélisation :
 * - Parallélisation au niveau de la boucle sur les paires (i, j)
 * - Chaque calcul Needleman-Wunsch est séquentiel
 * - Calcul triangulaire (i < j) puis recopie symétrique
 * - Une seule région parallèle pour toute la matrice
 * 
 * Compilation : g++ -O3 -std=c++17 -fopenmp build_matrix_needleman.cpp -o build_matrix_needleman
 * Usage : OMP_NUM_THREADS=12 ./build_matrix_needleman sequences.fasta
 */

#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <stdexcept>

using namespace std;

// ============================================================================
// PARAMÈTRES DE L'ALGORITHME NEEDLEMAN-WUNSCH
// ============================================================================

/** @brief Score d'un match (caractères identiques). */
const int MATCH = 1;
/** @brief Pénalité d'un mismatch (caractères différents). */
const int MISMATCH = -1;
/** @brief Pénalité d'ouverture d'un trou (gap opening). */
const int GAP_OPENING = -3;
/** @brief Pénalité d'extension d'un trou (gap extension). */
const int GAP_EXTENSION = -1;

// ============================================================================
// LECTURE DU FICHIER FASTA
// ============================================================================

/**
 * @brief Lit un fichier FASTA et retourne les séquences.
 * 
 * Format supposé :
 * - Les lignes commençant par '>' sont des en-têtes (ignorées)
 * - Les lignes suivantes contiennent la séquence
 * - Les séquences multi-lignes sont concaténées
 * 
 * @param filename Chemin du fichier FASTA
 * @return Vector de séquences (strings).
 *
 * @throw runtime_error si le fichier ne peut pas être ouvert.
 */
vector<string> readFasta(const string& filename) {
    ifstream in(filename);
    if (!in) {
        throw runtime_error("Impossible d'ouvrir " + filename);
    }

    vector<string> seqs;
    string line, current;

    while (getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!current.empty()) {
                seqs.push_back(current);
                current.clear();
            }
        } else {
            current += line;
        }
    }
    if (!current.empty()) {
        seqs.push_back(current);
    }

    return seqs;
}

// ============================================================================
// ALGORITHME NEEDLEMAN-WUNSCH SÉQUENTIEL
// ============================================================================

/**
 * @brief Calcule le score d'alignement entre deux séquences.
 * 
 * Implémentation Needleman-Wunsch avec pénalités affines et trois matrices :
 * - M  : scores de match/mismatch
 * - IX : scores avec gap dans seq2 (mouvement vertical)
 * - IY : scores avec gap dans seq1 (mouvement horizontal)
 * 
 * @param seq1 Première séquence
 * @param seq2 Deuxième séquence
 * @return Score d'alignement optimal.
 */
int calculer_score_needleman(const string& seq1, const string& seq2) {
    int len1 = seq1.length();
    int len2 = seq2.length();
    int largeur = len2 + 1;
    int taille = (len1 + 1) * largeur;

    // Allocation des trois matrices
    int* m = new int[taille];
    int* ix = new int[taille];
    int* iy = new int[taille];

    // Initialisation à zéro
    memset(m, 0, taille * sizeof(int));
    memset(ix, 0, taille * sizeof(int));
    memset(iy, 0, taille * sizeof(int));

    // Initialisation des bords
    m[0] = 0;
    ix[0] = GAP_OPENING;
    iy[0] = GAP_OPENING;

    // Première ligne (gaps dans seq1)
    for (int j = 1; j <= len2; j++) {
        m[j] = GAP_OPENING + j * GAP_EXTENSION;
        ix[j] = GAP_OPENING + j * GAP_EXTENSION;
        iy[j] = GAP_OPENING + j * GAP_EXTENSION;
    }

    // Première colonne (gaps dans seq2)
    for (int i = 1; i <= len1; i++) {
        m[i * largeur] = GAP_OPENING + i * GAP_EXTENSION;
        ix[i * largeur] = GAP_OPENING + i * GAP_EXTENSION;
        iy[i * largeur] = GAP_OPENING + i * GAP_EXTENSION;
    }

    // Remplissage de la matrice
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            int idx = i * largeur + j;
            int idx_haut = (i - 1) * largeur + j;
            int idx_gauche = i * largeur + (j - 1);
            int idx_diag = (i - 1) * largeur + (j - 1);

            // Score de match ou mismatch
            int score_match = (seq1[i - 1] == seq2[j - 1]) ? MATCH : MISMATCH;

            // Calcul des trois matrices
            ix[idx] = max(m[idx_haut] + GAP_OPENING, 
                          ix[idx_haut] + GAP_EXTENSION);
            iy[idx] = max(m[idx_gauche] + GAP_OPENING, 
                          iy[idx_gauche] + GAP_EXTENSION);
            m[idx] = max({m[idx_diag] + score_match, ix[idx], iy[idx]});
        }
    }

    // Score final (coin inférieur droit)
    int score_final = m[len1 * largeur + len2];

    // Libération mémoire
    delete[] m;
    delete[] ix;
    delete[] iy;

    return score_final;
}

// ============================================================================
// ÉCRITURE DU GRAPHE DOT
// ============================================================================

/**
 * @brief Écrit un graphe pondéré au format DOT.
 * 
 * @param filename Nom du fichier DOT à générer
 * @param dist Matrice des distances (stockée à plat, row-major).
 * @param n Nombre de séquences.
 * @param epsilon Seuil : une arête est créée si distance < epsilon.
 */
void writeDotGraph(const string& filename, const vector<int>& dist, int n, int epsilon) {
    ofstream out(filename);
    if (!out) {
        throw runtime_error("Impossible d'ouvrir " + filename + " en écriture");
    }

    out << "graph graphe_pondere {\n";
    out << "    node [shape=circle, style=filled, color=lightyellow, fontcolor=black];\n";
    out << "    edge [color=black, fontcolor=blue];\n\n";

    // Déclaration des sommets
    for (int i = 0; i < n; ++i) {
        out << "    A" << (i + 1) << " [label=\"" << i << "\"];\n";
    }
    out << "\n";

    // Arêtes non orientées (i < j uniquement pour éviter les doublons)
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

// ============================================================================
// PROGRAMME PRINCIPAL
// ============================================================================

/**
 * @brief Programme principal : construction de la matrice de distances.
 * 
 * Étapes :
 * 1. Lecture du fichier FASTA
 * 2. Calcul parallèle de toutes les paires de distances
 * 3. Génération du fichier DOT
 * 
 * Stratégie OpenMP :
 * - Une seule région parallèle (#pragma omp parallel for)
 * - Parallélisation au niveau de la boucle sur les lignes i
 * - Schedule dynamique (5 lignes par chunk) pour équilibrer la charge
 * - Calcul triangulaire, puis recopie symétrique
 */
int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <fichier.fasta>\n";
        return 1;
    }

    const string fastaFile = argv[1];
    const string dotFile = "../../DATA/Resulat_sequence_by_premier_algo.dot";
    const int epsilon = 70;
    int nthreads = omp_get_max_threads();

    cout << "Threads: " << nthreads << "\n";

    // Lecture du fichier FASTA
    vector<string> sequences;
    int n = 0, L = 0;

    try {
        sequences = readFasta(fastaFile);
        n = sequences.size();
        if (n == 0) {
            throw runtime_error("Aucune sequence lue");
        }
        L = sequences[0].length();
        cout << "Sequences: " << n << " (L=" << L << ")\n";
    } catch (const exception& e) {
        cerr << "Erreur: " << e.what() << "\n";
        return 1;
    }

    // Allocation de la matrice de distances
    vector<int> dist(n * n, 0);
    
    // Début du chronomètre
    double t0 = omp_get_wtime();

    // Calcul parallèle de la matrice de distances
    #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 5)
    for (int i = 0; i < n; ++i) {
        // Diagonale : distance nulle
        dist[i * n + i] = 0;
        
        // Calcul de la ligne i (partie triangulaire supérieure uniquement)
        for (int j = i + 1; j < n; ++j) {
            // Calcul du score Needleman-Wunsch
            int score = calculer_score_needleman(sequences[i], sequences[j]);
            
            // Conversion score -> distance
            int d = (L - score) / 2;
            
            // Matrice symétrique
            dist[i * n + j] = d;
            dist[j * n + i] = d;
        }
    }

    // Fin du chronomètre
    double t1 = omp_get_wtime();
    cout << "Temps: " << (t1 - t0) << " s\n";

    // Écriture du fichier DOT
    try {
        writeDotGraph(dotFile, dist, n, epsilon);
        cout << "Fichier: " << dotFile << "\n";
    } catch (const exception& e) {
        cerr << "Erreur: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
