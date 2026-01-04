#include "needleman_common.hpp"

#include <vector>

int scoreNeedleman(const char* a,
                   const char* b,
                   int n,
                   int m,
                   const ParametresNW& params) {
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

    for (int i = 1; i <= n; ++i) {
        int row = i * cols;
        int row_prev = (i - 1) * cols;
        for (int j = 1; j <= m; ++j) {
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

    int last = n * cols + m;
    return max3(M[last], X[last], Y[last]);
}
