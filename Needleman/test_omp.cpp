#include "needleman_common.hpp"

#include <cstring>
#include <iostream>

using namespace std;

int main() {
    const char* A = "AGCT";
    const char* B = "ATGC";
    int nthreads = 4;

    int n = static_cast<int>(strlen(A));
    int m = static_cast<int>(strlen(B));

    ParametresNW params;
    int score = scoreNeedlemanOMP(A, B, n, m, params, nthreads);

    cout << "Sequence A: " << A << "\n";
    cout << "Sequence B: " << B << "\n";
    cout << "Threads: " << nthreads << "\n";
    cout << "Score Needleman OMP: " << score << "\n";

    // Affichages supplementaires (a commenter si besoin)
    cout << "Longueur A: " << n << "\n";
    cout << "Longueur B: " << m << "\n";
    cout << "match=" << params.match
         << " mismatch=" << params.mismatch
         << " gap_open=" << params.gap_open
         << " gap_extend=" << params.gap_extend << "\n";

    // Matrices M/X/Y (a commenter si besoin)
    afficherMatricesDebug(A, B, n, m, params);

    return 0;
}
