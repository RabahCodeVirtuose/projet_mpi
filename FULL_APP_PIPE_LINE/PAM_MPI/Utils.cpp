#define OMPI_SKIP_MPICXX 1
/**
 * @file Utils.cpp
 * @brief Implémentation des fonctions utilitaires pour PAM (affichage,
 *        lecture / écriture de matrices de distances et de résultats).
 */
#include "Utils.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <sstream>

/**
 * @brief Affiche une matrice n × m sur la sortie standard.
 */
void affichage(const int* tab, int n, int m, int format) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            std::cout << std::setw(format) << tab[i * m + j] << " ";
        }
        std::cout << "\n";
    }
}

/**
 * @brief Lit une matrice de distances carrée depuis un fichier texte.
 */
std::vector<int> readDistanceMatrix(const std::string& filename, int& n_out) {
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("Impossible d'ouvrir le fichier de distances: " + filename);
    }

    int n, m;
    if (!(in >> n >> m)) {
        throw std::runtime_error("Lecture de n, m echouee dans le fichier de distances");
    }
    if (n != m) {
        std::ostringstream oss;
        oss << "La matrice de distances n'est pas carree : n=" << n << ", m=" << m;
        throw std::runtime_error(oss.str());
    }
    n_out = n;

    std::vector<int> dist(n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(in >> dist[i * n + j])) {
                std::ostringstream oss;
                oss << "Erreur de lecture de la distance (" << i << "," << j << ")";
                throw std::runtime_error(oss.str());
            }
        }
    }

    return dist;
}

/**
 * @brief Écrit un résultat PAM détaillé dans un fichier texte.
 */
void writePAMResult(const std::string& filename, const PAMResult& res) {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Impossible d'ouvrir le fichier de resultat PAM: " + filename);
    }

    int n = (int)res.clusterOf.size();
    int k = (int)res.medoids.size();

    out << "# PAM results\n";
    out << "# n = " << n << "\n";
    out << "# k = " << k << "\n";
    out << "# total_cost = " << res.totalCost << "\n\n";

    out << "# medoids:\n";
    for (int m = 0; m < k; ++m) {
        out << res.medoids[m] << (m + 1 < k ? ' ' : '\n');
    }
    out << "\n";

    out << "# columns: vertex cluster medoid dist\n";

    for (int i = 0; i < n; ++i) {
        int cluster = res.clusterOf[i];      // indice du cluster (0..k-1)
        int medoid  = res.medoids[cluster];  // indice du sommet médioïde
        int d       = res.distToMedoid[i];

        out << i << " "
            << cluster << " "
            << medoid << " "
            << d << "\n";
    }
}
