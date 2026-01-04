#include "needleman_common.hpp"

#include <cstring>
#include <iostream>
#include <omp.h>

int main() {
    const char* A = "AGCT";
    const char* B = "ATGC";

    int n = static_cast<int>(std::strlen(A));
    int m = static_cast<int>(std::strlen(B));

    ParametresNW params;
    double t0 = omp_get_wtime();
    int score = scoreNeedleman(A, B, n, m, params);
    double t1 = omp_get_wtime();

    std::cout << "Sequence A: " << A << "\n";
    std::cout << "Sequence B: " << B << "\n";
    std::cout << "Score Needleman (seq): " << score << "\n";
    std::cout << "Temps: " << (t1 - t0) << " s\n";

    return 0;
}
