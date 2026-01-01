#include "needleman_common.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

static void afficherSeparateur(int cols, int cell_width) {
    cout << '+';
    for (int i = 0; i < cols; ++i) {
        cout << string(cell_width + 2, '-') << '+';
    }
    cout << "\n";
}

static void afficherMatrice(const string& name,
                            const vector<vector<int>>& mat,
                            const char* A,
                            const char* B,
                            int n,
                            int m,
                            int neg_inf) {
    const int cell_width = 6;
    int data_cols = m + 1;
    int total_cols = data_cols + 1;

    cout << "\n" << name << "\n";
    afficherSeparateur(total_cols, cell_width);

    cout << "| " << setw(cell_width) << ' ' << " |";
    cout << " " << setw(cell_width) << '-' << " |";
    for (int j = 0; j < m; ++j) {
        cout << " " << setw(cell_width) << B[j] << " |";
    }
    cout << "\n";
    afficherSeparateur(total_cols, cell_width);

    for (int i = 0; i <= n; ++i) {
        char row = (i == 0) ? '-' : A[i - 1];
        cout << "| " << setw(cell_width) << row << " |";
        for (int j = 0; j <= m; ++j) {
            int v = mat[i][j];
            string cell = (v <= neg_inf / 2) ? "-INF" : to_string(v);
            cout << " " << setw(cell_width) << cell << " |";
        }
        cout << "\n";
        afficherSeparateur(total_cols, cell_width);
    }
}

void afficherMatricesDebug(const char* a,
                           const char* b,
                           int n,
                           int m,
                           const ParametresNW& params) {
    const int neg_inf = -1000000000;

    vector<vector<int>> M(n + 1, vector<int>(m + 1, neg_inf));
    vector<vector<int>> X(n + 1, vector<int>(m + 1, neg_inf));
    vector<vector<int>> Y(n + 1, vector<int>(m + 1, neg_inf));

    M[0][0] = 0;
    for (int i = 1; i <= n; ++i) {
        X[i][0] = params.gap_open + (i - 1) * params.gap_extend;
    }
    for (int j = 1; j <= m; ++j) {
        Y[0][j] = params.gap_open + (j - 1) * params.gap_extend;
    }

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int s = (a[i - 1] == b[j - 1]) ? params.match : params.mismatch;
            int diag = max({M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]});
            M[i][j] = diag + s;
            X[i][j] = max(M[i - 1][j] + params.gap_open,
                          X[i - 1][j] + params.gap_extend);
            Y[i][j] = max(M[i][j - 1] + params.gap_open,
                          Y[i][j - 1] + params.gap_extend);
        }
    }

    cout << "Matrices M/X/Y (debug):\n";
    cout << "A = " << a << "\n";
    cout << "B = " << b << "\n";
    cout << "match=" << params.match
         << " mismatch=" << params.mismatch
         << " gap_open=" << params.gap_open
         << " gap_extend=" << params.gap_extend << "\n";

    afficherMatrice("M (match/mismatch)", M, a, b, n, m, neg_inf);
    afficherMatrice("X (gap in B, vertical)", X, a, b, n, m, neg_inf);
    afficherMatrice("Y (gap in A, horizontal)", Y, a, b, n, m, neg_inf);
}
