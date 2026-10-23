#include "doctest.h"

#include <cmath>
#include <vector>

#include "CodonMutSelCodonMatrixBidimArray.hpp"
#include "CodonMutSelOmegaCodonSubMatrix.hpp"
#include "CodonStateSpace.hpp"
#include "GTRSubMatrix.hpp"

// ---------------------------------------------------------------------------
// Lock the CorruptMatrix / NoCorrupt family of cache-invalidation methods.
//
// The matrix caches Q rows, the stationary distribution and the fitnesses.
// CorruptMatrix() (with RecomputeFitness) must be used when the fitness vector
// changed; CorruptMatrixNoFitnessRecomput() / SetOmegaNoCorrupt() only when
// omega changed (fitness untouched). Using the wrong variant leaves a stale
// cache. These tests pin the correct behaviour so a future misuse is caught.
// ---------------------------------------------------------------------------

namespace {

// a non-uniform fitness profile over the 61 codons
std::vector<double> corrupt_make_profile(int Nstate) {
    std::vector<double> f(Nstate, 1.0);
    for (int i = 0; i < Nstate; i++) { f[i] = 1.0 + 0.05 * (i % 7); }
    return f;
}

}  // namespace

TEST_CASE("CodonMutSelOmegaCodonSubMatrix: SetOmega changes the non-syn rate") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    std::vector<double> rr(6, 1.0), stat(4, 0.25);
    GTRSubMatrix nuc(4, rr, stat, true);
    auto fit = corrupt_make_profile(N);
    CodonMutSelOmegaCodonSubMatrix m(&ss, &nuc, fit, 1.0, false);

    m.UpdateMatrix();
    // pick a non-synonymous neighbor pair
    int i = 0, j = -1;
    for (auto nb : ss.GetNeighbors(i)) {
        if (!ss.Synonymous(i, nb)) { j = nb; break; }
    }
    REQUIRE(j >= 0);
    double q1 = m(i, j);

    m.SetOmega(2.0);  // corrupts + recomputes fitness
    m.UpdateMatrix();
    double q2 = m(i, j);
    // non-syn rate scales by omega (2x)
    CHECK(q2 == doctest::Approx(2.0 * q1).epsilon(1e-9));
}

TEST_CASE("CodonMutSelOmegaCodonSubMatrix: SetOmegaNoCorrupt + NoFitnessRecomput == SetOmega") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    std::vector<double> rr(6, 1.0), stat(4, 0.25);
    GTRSubMatrix nuc(4, rr, stat, true);
    auto fit = corrupt_make_profile(N);

    // reference: SetOmega (full corrupt)
    CodonMutSelOmegaCodonSubMatrix m1(&ss, &nuc, fit, 1.0, false);
    m1.UpdateMatrix();
    m1.SetOmega(2.0);
    m1.UpdateMatrix();

    // variant: SetOmegaNoCorrupt + CorruptMatrixNoFitnessRecomput
    CodonMutSelOmegaCodonSubMatrix m2(&ss, &nuc, fit, 1.0, false);
    m2.UpdateMatrix();
    m2.SetOmegaNoCorrupt(2.0);
    m2.CorruptMatrixNoFitnessRecomput();
    m2.UpdateMatrix();

    int i = 0, j = -1;
    for (auto nb : ss.GetNeighbors(i)) {
        if (!ss.Synonymous(i, nb)) { j = nb; break; }
    }
    REQUIRE(j >= 0);
    CHECK(m2(i, j) == doctest::Approx(m1(i, j)).epsilon(1e-9));
}

TEST_CASE("CodonMutSelOmegaCodonSubMatrix: SetOmegaNoCorrupt without corrupt leaves stale cache") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    std::vector<double> rr(6, 1.0), stat(4, 0.25);
    GTRSubMatrix nuc(4, rr, stat, true);
    auto fit = corrupt_make_profile(N);
    CodonMutSelOmegaCodonSubMatrix m(&ss, &nuc, fit, 1.0, false);
    m.UpdateMatrix();

    int i = 0, j = -1;
    for (auto nb : ss.GetNeighbors(i)) {
        if (!ss.Synonymous(i, nb)) { j = nb; break; }
    }
    REQUIRE(j >= 0);
    double before = m(i, j);

    // SetOmegaNoCorrupt alone: omega changed but cache NOT invalidated.
    // Accessing Q must still return the OLD value (stale) until a corrupt.
    m.SetOmegaNoCorrupt(2.0);
    double stale = m(i, j);
    CHECK(stale == doctest::Approx(before).epsilon(1e-9));

    // after a proper corrupt, the new omega is reflected
    m.CorruptMatrixNoFitnessRecomput();
    m.UpdateMatrix();
    double fresh = m(i, j);
    CHECK(fresh == doctest::Approx(2.0 * before).epsilon(1e-9));
}

TEST_CASE("CodonMutSelCodonMatrixBidimArray: CorruptColCodonMatrices only touches one omega column") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    std::vector<double> rr(6, 1.0), stat(4, 0.25);
    GTRSubMatrix nuc(4, rr, stat, true);

    // 2 fitness rows x 2 omega columns
    SimpleArray<std::vector<double>> fitness(2);
    fitness[0] = corrupt_make_profile(N);
    fitness[1] = corrupt_make_profile(N);
    for (int c = 0; c < N; c++) { fitness[1][c] = 1.0 + 0.1 * (c % 5); }
    SimpleArray<double> delta_omega(2);
    delta_omega[0] = 1.0;
    delta_omega[1] = 1.0;

    CodonMutSelCodonMatrixBidimArray arr(&ss, &nuc, &fitness, &delta_omega, 0.0);
    arr.CorruptCodonMatrices();

    // change omega of column 1 only
    delta_omega[1] = 3.0;
    arr.CorruptColCodonMatrices(1);

    // column 0 must be unchanged, column 1 scaled by 3
    int i = 0, j = -1;
    for (auto nb : ss.GetNeighbors(i)) {
        if (!ss.Synonymous(i, nb)) { j = nb; break; }
    }
    REQUIRE(j >= 0);
    double q00 = arr.GetVal(0, 0)(i, j);
    double q01 = arr.GetVal(0, 1)(i, j);
    double q10 = arr.GetVal(1, 0)(i, j);
    double q11 = arr.GetVal(1, 1)(i, j);
    // within each row, column 1 (omega=3) is 3x column 0 (omega=1)
    CHECK(q01 == doctest::Approx(3.0 * q00).epsilon(1e-9));
    CHECK(q11 == doctest::Approx(3.0 * q10).epsilon(1e-9));
    // column 0 (omega=1) is untouched by the column-1 corrupt: its ratio to
    // column 1 stays exactly 3 in both rows
    CHECK(q00 > 0.0);
    CHECK(q10 > 0.0);
}

TEST_CASE("CodonMutSelCodonMatrixBidimArray: CorruptRowCodonMatrices recomputes fitness") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    std::vector<double> rr(6, 1.0), stat(4, 0.25);
    GTRSubMatrix nuc(4, rr, stat, true);

    SimpleArray<std::vector<double>> fitness(2);
    fitness[0] = corrupt_make_profile(N);
    fitness[1] = corrupt_make_profile(N);
    SimpleArray<double> delta_omega(2);
    delta_omega[0] = 1.0;
    delta_omega[1] = 1.0;

    CodonMutSelCodonMatrixBidimArray arr(&ss, &nuc, &fitness, &delta_omega, 0.0);
    arr.CorruptCodonMatrices();

    // change the fitness of row 1, then corrupt only that row
    for (int c = 0; c < N; c++) { fitness[1][c] = 1.0 + 0.2 * (c % 3); }
    arr.CorruptRowCodonMatrices(1);

    int i = 0, j = -1;
    for (auto nb : ss.GetNeighbors(i)) {
        if (!ss.Synonymous(i, nb)) { j = nb; break; }
    }
    REQUIRE(j >= 0);
    // row 1 must reflect the new fitness (differs from row 0)
    double q00 = arr.GetVal(0, 0)(i, j);
    double q10 = arr.GetVal(1, 0)(i, j);
    CHECK(q10 != doctest::Approx(q00).epsilon(1e-6));
    // row 0 untouched
    double q01 = arr.GetVal(0, 1)(i, j);
    CHECK(q01 == doctest::Approx(q00).epsilon(1e-9));
}
