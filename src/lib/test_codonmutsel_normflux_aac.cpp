#include "doctest.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "AACodonMutSelOmegaCodonSubMatrix.hpp"
#include "CodonStateSpace.hpp"
#include "GTRSubMatrix.hpp"

// ---------------------------------------------------------------------------
// Normalisation of the AACodonMutSel codon matrix (mutselaac branch).
//
// The nucleotide GTR as used by LFP (GetRateBayesCode) is normalised to the
// *mutation* flux 2*norm*3 = 1 per codon (norm = sum_{i<j} pi_i pi_j rr_ij).
// The mutselc model builds its codon matrix with innormalise=false on top of
// that GTRSubMatrixCodon, so its branch lengths are in *mutation* units and LFP
// recovers exactly the simulated substitution count (Nsub = blength * flux).
//
// The AAC model instead builds its codon matrix with innormalise=true
// (AACodonMutSelCodonMatrixBidimArray.hpp): Normalise() divides by
// SubMatrix::GetRate() = 2 * sum pi_c (-Q_cc) over the *substitution* flux
// (selection included), so its branch lengths are in *substitution* units.
// This test measures that scale difference.
// ---------------------------------------------------------------------------

namespace {

double aacfx_flux(const SubMatrix &m, int Nstate) {
    double f = 0;
    for (int c = 0; c < Nstate; c++) { f += m.Stationary(c) * -m(c, c); }
    return f;
}

std::vector<double> aacfx_aa20(double concentration) {
    std::mt19937 gen(12345);
    std::gamma_distribution<double> g(concentration, 1.0);
    std::vector<double> f(20);
    double tot = 0;
    for (auto &x : f) {
        x = g(gen);
        tot += x;
    }
    for (auto &x : f) { x /= tot; }
    return f;
}

std::vector<double> aacfx_codon61(double concentration) {
    std::mt19937 gen(54321);
    std::gamma_distribution<double> g(concentration, 1.0);
    std::vector<double> f(61);
    double tot = 0;
    for (auto &x : f) {
        x = g(gen);
        tot += x;
    }
    for (auto &x : f) { x /= tot; }
    return f;
}

struct aacfx_Nuc {
    std::vector<double> rr{0.149, 0.249, 0.080, 0.273, 0.163, 0.086};
    std::vector<double> stat{0.230, 0.231, 0.231, 0.308};
    GTRSubMatrixCodon codon;  // 2*norm*3 : used by mutselc and LFP
    GTRSubMatrix plain;       // 2*norm     : used by the AAC model
    aacfx_Nuc() : codon(4, rr, stat, true), plain(4, rr, stat, true) {}
};

}  // namespace

TEST_CASE("normflux AAC: codon flux of the AAC matrix vs LFP mutation scale") {
    CodonStateSpace ss(Universal);
    aacfx_Nuc nuc;
    auto aa = aacfx_aa20(0.3);
    auto codon = aacfx_codon61(0.3);

    // The AAC model now builds its codon matrix with innormalise=false on a
    // GTRSubMatrixCodon (2*norm*3 mutation scale), exactly as mutselc and LFP:
    // the flux is the real substitution flux (< 1 under purifying selection).
    AACodonMutSelOmegaCodonSubMatrix sel_mut(&ss, &nuc.codon, codon, aa, 1.0, false);
    // innormalise=true would force the substitution flux to 1 (substitution units).
    AACodonMutSelOmegaCodonSubMatrix sel_subs(&ss, &nuc.codon, codon, aa, 1.0, true);

    double f_mut = aacfx_flux(sel_mut, 61);
    double f_subs = aacfx_flux(sel_subs, 61);

    std::printf("[normflux AAC] flux false/mutation-scale (model) %.4f | true/substitution-scale %.4f\n",
        f_mut, f_subs);

    // model convention (innormalise=false on GTRSubMatrixCodon): real substitution flux < 1
    CHECK(f_mut != doctest::Approx(1.0).epsilon(1e-3));
    CHECK(f_mut < 1.0);
    CHECK(f_mut > 0.1);
    // innormalise=true would force flux = 1
    CHECK(f_subs == doctest::Approx(1.0).epsilon(1e-6));

    // In mutation units (LFP/mutselc), a branch of length blength carries blength*f_mut
    // substitutions per codon. This is the scale LFP recovers.
    std::printf("[normflux AAC] substitution/mutation scale factor = %.4f\n", 1.0 / f_mut);
}
