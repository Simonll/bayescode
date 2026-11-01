#include "doctest.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "AACodonMutSelOmegaCodonSubMatrix.hpp"
#include "CodonStateSpace.hpp"
#include "GTRSubMatrix.hpp"

// ---------------------------------------------------------------------------
// Reversibility / detailed-balance of the AACodonMutSel codon matrix under the
// SENCA parameterisation (mutselsenca branch).
//
// In SENCA-hybrid, the matrix receives codonweight_c = n_aa * theta_c (where
// theta sums to 1 within each amino acid; n_aa = degeneracy). We verify the
// three classical checks a valid CTMC codon matrix must pass:
//   1. rows/cols sum to 0 on the diagonal (generator, trace 0)
//   2. detailed balance: pi_i * Q_ij == pi_j * Q_ji   (reversibility)
//   3. pi is the stationary distribution: pi * Q == 0  (approx)
// Under normalized pi, (2) implies (3).
// ---------------------------------------------------------------------------

static void BuildSencaCodonWeights(CodonStateSpace &ss, std::vector<double> &codonweight,
    std::vector<double> &aa) {
    // 61-vector; codonweight = n_aa * theta, theta normalized to sum 1 per AA.
    int nstate = ss.GetNstate();
    codonweight.assign(nstate, 0.0);
    aa.assign(20, 0.0);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> u(0.5, 2.0);
    for (int a = 0; a < 20; a++) { aa[a] = u(rng); }

    // accumulate theta per amino acid then build n_aa * theta
    std::vector<double> theta(nstate, 0.0);
    std::vector<double> tot(20, 0.0);
    // assign random positive theta, then normalize per AA in a second pass
    for (int c = 0; c < nstate; c++) { theta[c] = u(rng); }
    // first pass: totals per AA
    std::vector<int> deg(20, 0);
    for (int c = 0; c < nstate; c++) {
        int a = ss.Translation(c);
        tot[a] += theta[c];
        deg[a]++;
    }
    // normalize per AA and multiply by degeneracy
    for (int c = 0; c < nstate; c++) {
        int a = ss.Translation(c);
        theta[c] /= tot[a];
        codonweight[c] = deg[a] * theta[c];
    }
}

TEST_CASE("AACodonMutSel SENCA: detailed balance and reversibility") {
    CodonStateSpace ss(Universal);
    int nstate = ss.GetNstate();
    REQUIRE(nstate == 61);

    // nucleotide GTR (reversible), parameters = reasonable
    std::vector<double> nucrelrate{0.5, 0.4, 0.3, 0.2, 0.1, 0.05};
    std::vector<double> nucstat{0.3, 0.2, 0.2, 0.3};
    GTRSubMatrixCodon nuc(4, nucrelrate, nucstat, true);

    std::vector<double> codonweight, aa;
    BuildSencaCodonWeights(ss, codonweight, aa);

    for (double omega : {1.0, 0.7, 1.3}) {
        AACodonMutSelOmegaCodonSubMatrix m(&ss, &nuc, codonweight, aa, omega, false);
        m.CorruptMatrix();
        m.GetStationary();

        double maxviol = 0.0;
        double sumrow = 0.0;
        for (int i = 0; i < nstate; i++) {
            double rowsum = 0.0;
            for (auto j : ss.GetNeighbors(i)) {
                double qij = m(i, j);
                double pi_i_qij = m.Stationary(i) * qij;
                double pi_j_qji = m.Stationary(j) * m(j, i);
                double rel = std::fabs(pi_i_qij - pi_j_qji) /
                             (std::fabs(pi_i_qij) + std::fabs(pi_j_qji) + 1e-12);
                if (rel > maxviol) { maxviol = rel; }
                rowsum += qij;
            }
            // row sum (off-diagonal) should equal -Q_ii
            REQUIRE(std::fabs(rowsum + m(i, i)) < 1e-9 * (1 + std::fabs(rowsum)));
            sumrow += std::fabs(rowsum);
        }

        std::printf("[omega %.2f] max relative detailed-balance violation = %.3g\n", omega, maxviol);
        CHECK(maxviol < 1e-9);

        // pi * Q ~ 0 for each column j
        for (int j = 0; j < nstate; j++) {
            double drift = 0.0;
            for (int i = 0; i < nstate; i++) { drift += m.Stationary(i) * m(i, j); }
            REQUIRE(std::fabs(drift) < 1e-9 * (1 + std::fabs(drift)));
        }
    }
}

TEST_CASE("AACodonMutSel SENCA: diagonal sums to -offdiag, rates nonneg") {
    CodonStateSpace ss(Universal);
    std::vector<double> nucrelrate{0.5, 0.4, 0.3, 0.2, 0.1, 0.05};
    std::vector<double> nucstat{0.3, 0.2, 0.2, 0.3};
    GTRSubMatrixCodon nuc(4, nucrelrate, nucstat, true);
    std::vector<double> codonweight, aa;
    BuildSencaCodonWeights(ss, codonweight, aa);
    AACodonMutSelOmegaCodonSubMatrix m(&ss, &nuc, codonweight, aa, 1.0, false);
    m.CorruptMatrix();
    int n = ss.GetNstate();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) { REQUIRE(m(i, j) >= 0); }
        }
        REQUIRE(m(i, i) <= 0);
    }
}
