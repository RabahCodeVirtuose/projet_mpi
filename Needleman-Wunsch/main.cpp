#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

static std::string centerText(const std::string& text, int width) {
    if (static_cast<int>(text.size()) >= width) {
        return text;
    }
    int pad = width - static_cast<int>(text.size());
    int left = pad / 2;
    int right = pad - left;
    return std::string(left, ' ') + text + std::string(right, ' ');
}

static std::string rightText(const std::string& text, int width) {
    if (static_cast<int>(text.size()) >= width) {
        return text;
    }
    return std::string(width - text.size(), ' ') + text;
}

static void printSeparator(std::size_t cols, int cell_width) {
    std::cout << '+';
    for (std::size_t i = 0; i < cols; ++i) {
        std::cout << std::string(cell_width + 2, '-') << '+';
    }
    std::cout << "\n";
}

static void printMatrix(const std::string& name,
                        const std::vector<std::vector<int>>& mat,
                        const std::string& A,
                        const std::string& B,
                        int neg_inf) {
    const int cell_width = 6;
    const std::size_t data_cols = B.size() + 1;
    const std::size_t total_cols = data_cols + 1;

    std::cout << "\n" << name << "\n";
    printSeparator(total_cols, cell_width);

    std::cout << '|';
    std::cout << ' ' << centerText("", cell_width) << ' ' << '|';
    std::cout << ' ' << centerText("-", cell_width) << ' ' << '|';
    for (char c : B) {
        std::cout << ' ' << centerText(std::string(1, c), cell_width) << ' ' << '|';
    }
    std::cout << "\n";
    printSeparator(total_cols, cell_width);

    for (std::size_t i = 0; i < mat.size(); ++i) {
        char row = (i == 0) ? '-' : A[i - 1];
        std::cout << '|';
        std::cout << ' ' << centerText(std::string(1, row), cell_width) << ' ' << '|';
        for (std::size_t j = 0; j < mat[0].size(); ++j) {
            int v = mat[i][j];
            std::string cell = (v <= neg_inf / 2) ? "-INF" : std::to_string(v);
            std::cout << ' ' << rightText(cell, cell_width) << ' ' << '|';
        }
        std::cout << "\n";
        printSeparator(total_cols, cell_width);
    }
}

int main() {
    std::string A = "AGCTATGC";
    std::string B = "ATGCAGCT";

    const int MATCH = 1;
    const int MISMATCH = -1;
    const int GAP_OPEN = -3;
    const int GAP_EXT = -1;
    const int NEG_INF = -1000000000;

    auto score = [&](char a, char b) {
        return (a == b) ? MATCH : MISMATCH;
    };

    int n = static_cast<int>(A.size());
    int m = static_cast<int>(B.size());

    std::vector<std::vector<int>> M(n + 1, std::vector<int>(m + 1, NEG_INF));
    std::vector<std::vector<int>> X(n + 1, std::vector<int>(m + 1, NEG_INF));
    std::vector<std::vector<int>> Y(n + 1, std::vector<int>(m + 1, NEG_INF));

    M[0][0] = 0;

    for (int i = 1; i <= n; ++i) { 
        X[i][0] = GAP_OPEN + (i - 1) * GAP_EXT;
    }
    for (int j = 1; j <= m; ++j) {
        Y[0][j] = GAP_OPEN + (j - 1) * GAP_EXT;
    }

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int s = score(A[i - 1], B[j - 1]);
            M[i][j] = std::max({M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]}) + s;
            X[i][j] = std::max(M[i - 1][j] + GAP_OPEN, X[i - 1][j] + GAP_EXT);
            Y[i][j] = std::max(M[i][j - 1] + GAP_OPEN, Y[i][j - 1] + GAP_EXT);
        }
    }

    std::cout << "Sequences:\n";
    std::cout << "A = " << A << "\n";
    std::cout << "B = " << B << "\n";
    std::cout << "Scoring: match=+" << MATCH
              << ", mismatch=" << MISMATCH
              << ", gap_open=" << GAP_OPEN
              << ", gap_extend=" << GAP_EXT << "\n";

    printMatrix("M (match/mismatch)", M, A, B, NEG_INF);
    printMatrix("X (gap in B, vertical)", X, A, B, NEG_INF);
    printMatrix("Y (gap in A, horizontal)", Y, A, B, NEG_INF);

    int best = std::max({M[n][m], X[n][m], Y[n][m]});
    std::cout << "\nFinal score = " << best << "\n";

    return 0;
}
