    #include "needleman_common.hpp"

#include <algorithm>
#include <vector>

using namespace std;

int scoreNeedleman(const char* a,
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

    return max({M[n][m], X[n][m], Y[n][m]});
}
