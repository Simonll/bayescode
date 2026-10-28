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

    // Innnormalise=false built on the *mutation*-normalised GTR (mutselc/LFP convention):
    // the flux is the real substitution flux (< 1 under selection).
    AACodonMutSelOmegaCodonSubMatrix sel_mut(&ss, &nuc.codon, codon, aa, 1.0, false);
    // Innnormalise=true (the value used by the AAC model): Normalise() forces the
    // substitution flux to 1 -> branch lengths in substitutions.
    AACodonMutSelOmegaCodonSubMatrix sel_subs(&ss, &nuc.codon, codon, aa, 1.0, true);
    // Also on the plain GTR (the AAC model's nucmatrix): innormalise=true.
    AACodonMutSelOmegaCodonSubMatrix sel_plain(&ss, &nuc.plain, codon, aa, 1.0, true);

    double f_mut = aacfx_flux(sel_mut, 61);
    double f_subs = aacfx_flux(sel_subs, 61);
    double f_plain = aacfx_flux(sel_plain, 61);

    std::printf("[normflux AAC] flux false/mutation-scale %.4f | true/substitution-scale %.4f | "
                "true on plain GTR %.4f\n",
        f_mut, f_subs, f_plain);

    // innormalise=true always forces flux = 1 (substitution units)
    CHECK(f_subs == doctest::Approx(1.0).epsilon(1e-6));
    // on the plain GTR it is also 1 (the AAC model normalises per-codon)
    CHECK(f_plain == doctest::Approx(1.0).epsilon(1e-6));
    // on the mutation-normalised GTR *without* extra normalisation, the flux is the real
    // substitution flux, < 1 under purifying selection
    CHECK(f_mut != doctest::Approx(1.0).epsilon(1e-3));
    CHECK(f_mut < 1.0);

    // The scale factor between the two conventions: LFP (mutation scale) reads branch
    // lengths as blength*flux_mut, so if the model uses substitution units the simulated
    // Nsub per codon on a branch of length blength would be blength*1 instead of blength*f_mut.
    std::printf("[normflux AAC] substitution/mutation scale factor = %.4f\n", 1.0 / f_mut);
}
