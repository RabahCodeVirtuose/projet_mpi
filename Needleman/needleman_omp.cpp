#include "needleman_common.hpp"

#include <vector>
#include <omp.h>

int scoreNeedlemanOMP(const char* a,
                      const char* b,
                      int n,
                      int m,
                      const ParametresNW& params,
                      int nthreads) {
    if (nthreads <= 0) {
        nthreads = 1;
    }

    const int neg_inf = -1000000000;
    const int cols = m + 1;
    const int size = (n + 1) * (m + 1);

    std::vector<int> M(size, neg_inf);
    std::vector<int> X(size, neg_inf);
    std::vector<int> Y(size, neg_inf);

    M[0] = 0;
    for (int i = 1; i <= n; ++i) {
        X[i * cols] = params.gap_open + (i - 1) * params.gap_extend;
    }
    for (int j = 1; j <= m; ++j) {
        Y[j] = params.gap_open + (j - 1) * params.gap_extend;
    }

    #pragma omp parallel num_threads(nthreads) shared(M, X, Y, a, b, n, m, params, cols)
    {
        for (int d = 2; d <= n + m; ++d) {
            int i_start = d - m;
            if (i_start < 1) i_start = 1;
            int i_end = d - 1;
            if (i_end > n) i_end = n;

            #pragma omp for schedule(static)
            for (int i = i_start; i <= i_end; ++i) {
                int j = d - i;
                int row = i * cols;
                int row_prev = (i - 1) * cols;
                int idx = row + j;
                int diag_idx = row_prev + j - 1;
                int up_idx = row_prev + j;
                int left_idx = row + j - 1;

                int s = (a[i - 1] == b[j - 1]) ? params.match : params.mismatch;
                M[idx] = max3(M[diag_idx], X[diag_idx], Y[diag_idx]) + s;
                X[idx] = max2(M[up_idx] + params.gap_open,
                              X[up_idx] + params.gap_extend);
                Y[idx] = max2(M[left_idx] + params.gap_open,
                              Y[left_idx] + params.gap_extend);
            }
        }
    }

    int last = n * cols + m;
    return max3(M[last], X[last], Y[last]);
}
