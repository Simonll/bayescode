#include "doctest.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "CodonMutSelOmegaCodonSubMatrix.hpp"
#include "CodonStateSpace.hpp"
#include "GTRSubMatrix.hpp"
#include "Random.hpp"

// ---------------------------------------------------------------------------
// Matrix normalisation (innormalise) and its impact on the branch length.
//
// blength is a TIME: the expected number of substitutions per codon along a branch
// is blength * flux, with flux = sum_c pi_c * (-Q_cc) the stationary flux of the codon
// matrix. The nucleotide matrix (GTRSubMatrixCodon) is normalised to a neutral flux of
// 1 per codon; the codon matrix of the model is built with innormalise=false, so its
// flux is the *real* flux (< 1 under purifying selection), whereas innormalise=true
// forces flux = 1.
//
// All helpers are prefixed nrmfx_ (all_tests.cpp is a single translation unit).
// ---------------------------------------------------------------------------

namespace {

//! nucleotide matrix as built by the model (CodonMutSelMultipleOmegaModel): GTRSubMatrixCodon, normalised
struct nrmfx_Nuc {
    std::vector<double> rr{0.149, 0.249, 0.080, 0.273, 0.163, 0.086};
    std::vector<double> stat{0.230, 0.231, 0.231, 0.308};
    GTRSubMatrixCodon codon;
    GTRSubMatrix plain;
    nrmfx_Nuc() : codon(4, rr, stat, true), plain(4, rr, stat, true) {}
};

//! flux sum_c pi_c * (-Q_cc) of an already-built matrix (works for any SubMatrix)
double nrmfx_flux(const SubMatrix &m, int Nstate) {
    double f = 0;
    for (int c = 0; c < Nstate; c++) { f += m.Stationary(c) * -m(c, c); }
    return f;
}

std::vector<double> nrmfx_uniform_fitness() { return std::vector<double>(61, 1.0 / 61); }

//! peaked, deterministic fitness profile (Dirichlet(0.3)-like), summing to 1
std::vector<double> nrmfx_selective_fitness() {
    std::mt19937 gen(12345);
    std::gamma_distribution<double> g(0.3, 1.0);
    std::vector<double> f(61);
    double tot = 0;
    for (auto &x : f) {
        x = g(gen);
        tot += x;
    }
    for (auto &x : f) { x /= tot; }
    return f;
}

//! jump-chain simulation along one branch of duration T (same algorithm as simumutselc),
//! each site starting from the stationary distribution; returns Nsub per codon
double nrmfx_simulate_nsub_per_codon(const CodonMutSelOmegaCodonSubMatrix &m, int Nsite, double T) {
    const int Nstate = 61;
    long long nsub = 0;
    for (int s = 0; s < Nsite; s++) {
        int cur = Random::DrawFromDiscreteDistribution(m.GetStationary(), Nstate);
        double time = 0;
        while (true) {
            double total = -m(cur, cur);
            if (total <= 0) { break; }
            time += Random::sExpo() / total;
            if (time >= T) { break; }
            double u = Random::Uniform() * total, acc = 0;
            int next = -1;
            for (int j = 0; j < Nstate; j++) {
                if (j == cur) { continue; }
                double q = m(cur, j);
                if (q > 0) {
                    acc += q;
                    next = j;
                    if (u < acc) { break; }
                }
            }
            cur = next;
            nsub++;
        }
    }
    return (double)nsub / Nsite;
}

}  // namespace

TEST_CASE("normflux A: stationary flux of the codon matrix vs innormalise and fitness") {
    CodonStateSpace ss(Universal);
    nrmfx_Nuc nuc;
    auto uni = nrmfx_uniform_fitness();
    auto sel = nrmfx_selective_fitness();

    // codon matrix: RecomputeFitness() is called by the constructor -> usable right away
    CodonMutSelOmegaCodonSubMatrix uni_f(&ss, &nuc.codon, uni, 1.0, false);
    CodonMutSelOmegaCodonSubMatrix uni_t(&ss, &nuc.codon, uni, 1.0, true);
    CodonMutSelOmegaCodonSubMatrix sel_f(&ss, &nuc.codon, sel, 1.0, false);
    CodonMutSelOmegaCodonSubMatrix sel_t(&ss, &nuc.codon, sel, 1.0, true);

    double f_uni_f = nrmfx_flux(uni_f, 61), f_uni_t = nrmfx_flux(uni_t, 61);
    double f_sel_f = nrmfx_flux(sel_f, 61), f_sel_t = nrmfx_flux(sel_t, 61);
    std::printf("[normflux A] codon flux  uniform: innormalise=false %.4f  true %.4f\n", f_uni_f,
        f_uni_t);
    std::printf("[normflux A] codon flux selective: innormalise=false %.4f  true %.4f\n", f_sel_f,
        f_sel_t);

    // innormalise=true -> flux 1 whatever the fitness
    CHECK(f_uni_t == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(f_sel_t == doctest::Approx(1.0).epsilon(1e-9));

    // innormalise=false -> real flux. Uniform: neutral flux minus the edges to stop codons (< 1);
    // selective: purifying selection lowers it further.
    CHECK(f_uni_f != doctest::Approx(1.0).epsilon(1e-3));
    CHECK(f_uni_f < 1.0);
    CHECK(f_uni_f > 0.9);
    CHECK(f_sel_f != doctest::Approx(1.0).epsilon(1e-3));
    CHECK(f_sel_f < f_uni_f);

    // with the *plain* GTRSubMatrix (no x3, flux 3 per codon) the uniform flux is ~2.9: this is
    // the "2.87" obtained when the nuc matrix is not the GTRSubMatrixCodon used by the model
    CodonMutSelOmegaCodonSubMatrix uni_plain(&ss, &nuc.plain, uni, 1.0, false);
    double f_plain = nrmfx_flux(uni_plain, 61);
    std::printf("[normflux A] uniform, plain GTRSubMatrix (x3 missing), innormalise=false %.4f\n",
        f_plain);
    CHECK(f_plain == doctest::Approx(3.0 * f_uni_f).epsilon(1e-6));
}

TEST_CASE("normflux A2: AA matrix needs CorruptMatrix() to fill its fitness cache") {
    CodonStateSpace ss(Universal);
    nrmfx_Nuc nuc;
    std::vector<double> aa(20, 1.0 / 20);
    // the AA constructor does NOT populate the fitness cache (unlike the codon one)
    AAMutSelOmegaCodonSubMatrix aa_t(&ss, &nuc.codon, aa, 1.0, 1.0, true);
    AAMutSelOmegaCodonSubMatrix aa_f(&ss, &nuc.codon, aa, 1.0, 1.0, false);
    aa_t.CorruptMatrix();
    aa_f.CorruptMatrix();
    double ft = nrmfx_flux(aa_t, 61), ff = nrmfx_flux(aa_f, 61);
    std::printf("[normflux A2] AA uniform flux: innormalise=true %.4f  false %.4f\n", ft, ff);
    CHECK(ft == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(ff != doctest::Approx(1.0).epsilon(1e-3));
}

TEST_CASE("normflux B: simulated Nsub = blength x flux; recovered length depends on innormalise") {
    CodonStateSpace ss(Universal);
    nrmfx_Nuc nuc;
    auto sel = nrmfx_selective_fitness();
    CodonMutSelOmegaCodonSubMatrix sel_f(&ss, &nuc.codon, sel, 1.0, false);
    CodonMutSelOmegaCodonSubMatrix sel_t(&ss, &nuc.codon, sel, 1.0, true);
    double flux_f = nrmfx_flux(sel_f, 61), flux_t = nrmfx_flux(sel_t, 61);

    const double T = 1.0;
    const int Nsite = 10000;
    Random::InitRandom(20240601);
    double nsub_f = nrmfx_simulate_nsub_per_codon(sel_f, Nsite, T);
    double nsub_t = nrmfx_simulate_nsub_per_codon(sel_t, Nsite, T);

    // moment estimator of the branch length given the matrix used by the analysis
    // (the Gamma(alpha+N, beta+B) conditioning of GammaWhiteNoise concentrates on N / (Nsite*flux))
    // nsub_* is already Nsub/codon (divided by Nsite inside the simulator)
    double rec_f_f = nsub_f / flux_f;  // simulated false, analysed false -> T
    double rec_f_t = nsub_f / flux_t;  // simulated false, analysed true  -> T*flux
    double rec_t_t = nsub_t / flux_t;  // simulated true,  analysed true  -> T
    std::printf("[normflux B] T=%.2f flux(false)=%.4f flux(true)=%.4f\n", T, flux_f, flux_t);
    std::printf("[normflux B] Nsub/codon: sim innormalise=false %.4f (expected %.4f) | true %.4f "
                "(expected %.4f)\n", nsub_f, T * flux_f, nsub_t, T * flux_t);
    std::printf("[normflux B] recovered length: false/false %.4f | false->true %.4f (T*flux=%.4f) | "
                "true/true %.4f\n", rec_f_f, rec_f_t, T * flux_f, rec_t_t);

    CHECK(nsub_f == doctest::Approx(T * flux_f).epsilon(0.04));
    CHECK(nsub_t == doctest::Approx(T * flux_t).epsilon(0.04));
    CHECK(nsub_t > nsub_f);  // same blength, more substitutions when normalised
    CHECK(rec_f_f == doctest::Approx(T).epsilon(0.04));
    CHECK(rec_t_t == doctest::Approx(T).epsilon(0.04));
    CHECK(rec_f_t == doctest::Approx(T * flux_f).epsilon(0.04));
    CHECK(rec_f_t < 0.95 * T);  // mismatched normalisation biases the length by the factor flux
}
