#include "doctest.h"

#include "CodonMutSelOmegaCodonSubMatrix.hpp"
#include "CodonStateSpace.hpp"
#include "GTRSubMatrix.hpp"

// ---------------------------------------------------------------------------
// FixationFactor: kimura-style fixation probability shared by ComputeArray and
// the dN/dS flow statistics. Fire-proof against the cases that were pinned down
// in the review (exploding exp() near 0, large positive / negative deltas).
// ---------------------------------------------------------------------------
TEST_CASE("CodonMutSelOmegaCodonSubMatrix: FixationFactor is monotone and bounded") {
    // neutral: factor ~= 1
    CHECK(FixationFactor(0.0) == doctest::Approx(1.0).epsilon(1e-6));
    // small deltas use the linear approximation 1 + delta/2
    CHECK(FixationFactor(1e-9) == doctest::Approx(1.0).epsilon(1e-6));
    // beneficial mutation: factor > 1, increases with deltaS
    CHECK(FixationFactor(0.5) > 1.0);
    CHECK(FixationFactor(1.0) > FixationFactor(0.5));
    // strongly beneficial: linear regime (deltaS)
    CHECK(FixationFactor(100.0) == doctest::Approx(100.0).epsilon(1e-9));
    // deleterious: factor in (0,1), decreasing
    CHECK(FixationFactor(-0.5) < 1.0);
    CHECK(FixationFactor(-0.5) > 0.0);
    CHECK(FixationFactor(-1.0) < FixationFactor(-0.5));
    // strongly deleterious: collapses to 0 (never negative)
    CHECK(FixationFactor(-100.0) == doctest::Approx(0.0).epsilon(1e-6));
    CHECK(FixationFactor(-100.0) >= 0.0);
    // symmetry identity: FixationFactor(a) = exp(a) * FixationFactor(-a)
    for (double a : {0.1, 0.5, 1.0, 2.0}) {
        CHECK(FixationFactor(a) == doctest::Approx(std::exp(a) * FixationFactor(-a)).epsilon(1e-6));
    }
}

TEST_CASE("CodonMutSelOmegaCodonSubMatrix: kFitnessFloor prevents log(0)") {
    // The fitness vector is floored at kFitnessFloor before log is taken, so a
    // zero codon in the input must not produce -inf. Exercise through a matrix.
    CodonStateSpace statespace(Universal);
    int Nstate = statespace.GetNstate();
    std::vector<double> rr(6, 1.0);
    std::vector<double> stat(4, 0.25);
    GTRSubMatrix nucmatrix(4, rr, stat, true);

    std::vector<double> fitness(Nstate, 0.0);
    fitness[0] = 1.0;  // only one codon occupied on purpose
    CodonMutSelOmegaCodonSubMatrix m(&statespace, &nucmatrix, fitness, 1.0, false);
    m.UpdateMatrix();

    // None of the log-fitnesses may be infinite / NaN, and fitness is floored.
    for (int i = 0; i < Nstate; i++) {
        CHECK(std::isfinite(m.GetLogFitness(i)));
        CHECK(m.GetFitness(i) >= kFitnessFloor);
    }
}

// ---------------------------------------------------------------------------
// GetRelativeFlowDNDS: the exact statistic reported in the trace. With a flat
// (neutral) profile and omega=1 the relative dN/dS must be ~1 (degenerate model:
// no selection), and never negative / NaN. This guards the computation used by
// GetPredictedRelativedNdS against regressions.
// ---------------------------------------------------------------------------
TEST_CASE("CodonMutSelOmegaCodonSubMatrix: relative dN/dS is sane on flat profile") {
    CodonStateSpace statespace(Universal);
    int Nstate = statespace.GetNstate();
    std::vector<double> rr(6, 1.0);
    std::vector<double> stat(4, 0.25);
    GTRSubMatrix nucmatrix(4, rr, stat, true);

    std::vector<double> fitness(Nstate, 1.0);  // flat => no selection
    CodonMutSelOmegaCodonSubMatrix m(&statespace, &nucmatrix, fitness, 1.0, false);

    double dn, ds;
    std::tie(dn, ds) = m.GetRelativeFlowDNDS();
    CHECK(std::isfinite(dn));
    CHECK(std::isfinite(ds));
    CHECK(ds > 0.0);
    double dnds = m.GetPredictedRelativeDNDS();
    CHECK(std::isfinite(dnds));
    // Neutral profile, omega=1: flows are exactly the neutral ones, so the
    // relative dN/dS is 1 up to rounding.
    CHECK(dnds == doctest::Approx(1.0).epsilon(1e-9));
}