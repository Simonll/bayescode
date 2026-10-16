#include "doctest.h"

#include "CodonMutSelCodonMatrixBidimArray.hpp"
#include "CodonStateSpace.hpp"
#include "GTRSubMatrix.hpp"

// Regression test for C1: the bidim array must build its codon matrices with
// normalise=false, so omega scales the overall substitution rate instead of
// cancelling out of it (and contaminating synonymous rates).
//
// With a GC-favoring fitness profile, doubling omega must raise the total rate
// by a factor strictly greater than 1 (only the non-synonymous part doubles),
// and the synonymous entry Q(0,1) must be identical at omega=1 and omega=2.
static double total_rate(const CodonMutSelOmegaCodonSubMatrix &m, CodonStateSpace &ss) {
    const auto &stat = m.GetStationary();
    double tot = 0;
    for (int i = 0; i < ss.GetNstate(); i++) {
        double row = 0;
        for (auto j : ss.GetNeighbors(i)) { row += m(i, j); }
        tot += stat[i] * row;
    }
    return tot;
}

TEST_CASE("CodonMutSelCodonMatrixBidimArray: omega scales the total rate") {
    CodonStateSpace statespace(Universal);
    int Nstate = statespace.GetNstate();

    std::vector<double> rr(6, 1.0);
    std::vector<double> stat(4, 0.25);
    GTRSubMatrix nucmatrix(4, rr, stat, true);

    // GC-favoring fitness profile (non-uniform, realistic)
    std::vector<double> fitness(Nstate, 0.0);
    for (int i = 0; i < Nstate; i++) {
        int a = statespace.GetCodonPosition(0, i);
        int b = statespace.GetCodonPosition(1, i);
        int c = statespace.GetCodonPosition(2, i);
        double gc = (a == 1 || a == 2) + (b == 1 || b == 2) + (c == 1 || c == 2);
        fitness[i] = std::exp(0.5 * gc);
    }
    double s = 0; for (auto &f : fitness) s += f;
    for (auto &f : fitness) f /= s;
    SimpleArray<std::vector<double>> fitnessarray(1, fitness);

    // two omega columns: delta=1 and delta=2, shift=0 => omega=1 and omega=2
    SimpleArray<double> deltaarray(2);
    deltaarray[0] = 1.0;
    deltaarray[1] = 2.0;
    double omega_shift = 0.0;

    CodonMutSelCodonMatrixBidimArray bidim(&statespace, &nucmatrix, &fitnessarray, &deltaarray,
        omega_shift);
    bidim.CorruptCodonMatrices();

    double t1 = total_rate(bidim.GetVal(0, 0), statespace);
    double t2 = total_rate(bidim.GetVal(0, 1), statespace);
    double ratio = t2 / t1;

    // omega must scale the rate: ratio strictly > 1 (here ~1.72 with this profile).
    // If normalise=true were used, omega would cancel and ratio would be ~1.0.
    CHECK(ratio > 1.1);
    CHECK(ratio < 2.5);

    // synonymous entry must be identical at omega=1 and omega=2
    int i0 = 0;
    auto nb = statespace.GetNeighbors(i0);
    double q1 = bidim.GetVal(0, 0)(i0, nb[0]);
    double q2 = bidim.GetVal(0, 1)(i0, nb[0]);
    CHECK(std::fabs(q1 - q2) < 1e-6);
}
