#ifndef NEEDLEMAN_COMMON_HPP
#define NEEDLEMAN_COMMON_HPP

struct ParametresNW {
    int match = 1;
    int mismatch = -1;
    int gap_open = -3;
    int gap_extend = -1;
};

int scoreNeedleman(const char* a,
                   const char* b,
                   int n,
                   int m,
                   const ParametresNW& params);

int scoreNeedlemanOMP(const char* a,
                      const char* b,
                      int n,
                      int m,
                      const ParametresNW& params,
                      int nthreads);

void afficherMatricesDebug(const char* a,
                           const char* b,
                           int n,
                           int m,
                           const ParametresNW& params);

#endif // NEEDLEMAN_COMMON_HPP
