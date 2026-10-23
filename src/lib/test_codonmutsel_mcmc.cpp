#include "doctest.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "CodonMutSelMultipleOmegaModel.hpp"
#include "CodonStateSpace.hpp"
#include "GTRSubMatrix.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"

// ---------------------------------------------------------------------------
// QC lot 4: gaps of the review (T1).
//   1. MCMC moves of CodonMutSelMultipleOmegaModel run, keep logprior/lnL
//      finite, keep the parameters on their simplex and leave the model state
//      coherent (accepted OR rejected).
//   2. The state read back by ReadCodonMutSelMultipleOmega (.param + .chain)
//      matches the state that was written, and what the reader derives from it
//      is coherent. The reader is a main() and cannot be linked in, so the
//      model-level reading path it uses (istream constructor + ChainReader) is
//      exercised directly and its source is scanned for the parts that cannot.
//   3. Non-regression of the matrices of the other models (MG omega, AA
//      mut-sel omega).
// ---------------------------------------------------------------------------
namespace {

const char *const mcmc_ali = "mutselc_mcmc_ali.phy";
const char *const mcmc_tree = "mutselc_mcmc_tree.tre";

std::string mcmc_src_dir() {
    std::string here = __FILE__;
    return here.substr(0, here.find_last_of("/\\") + 1) + "../";
}

std::string mcmc_slurp(const std::string &path) {
    std::ifstream in(path);
    REQUIRE(in.good());
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

void mcmc_write_dataset() {
    std::vector<std::string> seq{
        "ATGGCCGATATTGTCAGCGTGAAAGCCTGCCCGCCG",
        "ATGGCCGATTTTGTCAGCGTGAAAGCCTGCCCGCCG",
        "ATGGCCGATATTGTCAGCGTGAAAGCCTGCCCGCCA",
        "ATGGCCGATATTGTCAGCATGAAAGCCTGCCCGCCG",
        "ATGGCCGATATTGTCAGCGTGAAAGCCTGCCCGCTG",
    };
    std::ofstream out(mcmc_ali);
    out << seq.size() << " " << seq[0].size() << "\n";
    for (size_t i = 0; i < seq.size(); i++) { out << "taxon" << i << "\t" << seq[i] << "\n"; }
    out.close();
    std::ofstream t(mcmc_tree);
    t << "(((taxon0:1,taxon1:1):1,(taxon2:1,taxon3:1):1):1,taxon4:1);\n";
}

void mcmc_cleanup_dataset() {
    std::remove(mcmc_ali);
    std::remove(mcmc_tree);
}

//! fully free model: non-flat nuc stat / rel rates / fitness, 3 profile
//! components, 2 omega components, omega estimated
std::unique_ptr<CodonMutSelMultipleOmegaModel> mcmc_make_free_model() {
    return std::unique_ptr<CodonMutSelMultipleOmegaModel>(new CodonMutSelMultipleOmegaModel(
        mcmc_ali, mcmc_tree, "Null", "Null", "Null",
        0,      // omegamode: free
        3,      // Ncat
        1,      // baseNcat
        2,      // omegaNcat
        0.0,    // omegashift
        false,  // flatfitness
        "Null",
        false,    // flatnucstat
        false));  // flatnucrelrate
}

void mcmc_check_simplex(const std::vector<double> &v) {
    double s = 0;
    for (double x : v) {
        CHECK(std::isfinite(x));
        CHECK(x > 0.0);
        s += x;
    }
    CHECK(s == doctest::Approx(1.0).epsilon(1e-9));
}

std::vector<double> mcmc_nucstat(CodonMutSelMultipleOmegaModel &m) {
    std::vector<double> v;
    for (int i = 0; i < Nnuc; i++) { v.push_back(m.GetNucStat(i)); }
    return v;
}

std::vector<double> mcmc_nucrr(CodonMutSelMultipleOmegaModel &m) {
    std::vector<double> v;
    for (int i = 0; i < Nrr; i++) { v.push_back(m.GetNucRR(i)); }
    return v;
}

//! the full state must be finite and on its simplices
void mcmc_check_model_state(CodonMutSelMultipleOmegaModel &m) {
    CHECK(std::isfinite(m.GetLogPrior()));
    CHECK(std::isfinite(m.GetLogLikelihood()));
    CHECK(std::isfinite(m.GetLogProb()));
    mcmc_check_simplex(mcmc_nucstat(m));
    mcmc_check_simplex(mcmc_nucrr(m));
    for (int k = 0; k < m.GetNcat(); k++) {
        std::vector<double> p;
        for (int c = 0; c < 61; c++) { p.push_back(m.GetProfileCodon(k, c)); }
        mcmc_check_simplex(p);
    }
    for (int s = 0; s < m.GetNsite(); s++) {
        CHECK(m.GetProfileAlloc(s) >= 0);
        CHECK(m.GetProfileAlloc(s) < m.GetNcat());
        CHECK(m.GetOmegaAlloc(s) >= 0);
        CHECK(m.GetOmegaAlloc(s) < m.GetOmegaNcat());
        CHECK(std::isfinite(m.GetSiteOmega(s)));
        CHECK(m.GetSiteOmega(s) >= 0.0);
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// 1. MCMC moves
// ---------------------------------------------------------------------------
TEST_CASE("mcmc: MoveNucStat keeps nucstat on the simplex and the state coherent") {
    mcmc_write_dataset();
    {
        auto model = mcmc_make_free_model();
        model->ResampleSub(1.0);
        model->CollectSitePathSuffStat();
        mcmc_check_model_state(*model);

        std::vector<double> initial = mcmc_nucstat(*model);
        bool moved = false;
        for (int rep = 0; rep < 20; rep++) {
            model->MoveNucStat();
            mcmc_check_model_state(*model);
            // the rate matrices must be consistent with the (accepted or
            // restored) parameters: the target of the move, recomputed from
            // scratch after an explicit matrix refresh, must not change.
            double before = model->NucRatesLogProb();
            model->UpdateMatrices();
            CHECK(model->NucRatesLogProb() == doctest::Approx(before).epsilon(1e-9));
            if (mcmc_nucstat(*model) != initial) { moved = true; }
        }
        CHECK(moved);  // 200 proposals: at least one must be accepted
    }
    mcmc_cleanup_dataset();
}

TEST_CASE("mcmc: MoveNucRates keeps nucrelrate on the simplex and the state coherent") {
    mcmc_write_dataset();
    {
        auto model = mcmc_make_free_model();
        model->ResampleSub(1.0);
        model->CollectSitePathSuffStat();

        std::vector<double> initial = mcmc_nucrr(*model);
        bool moved = false;
        for (int rep = 0; rep < 20; rep++) {
            model->MoveNucRates();
            mcmc_check_model_state(*model);
            double before = model->NucRatesLogProb();
            model->UpdateMatrices();
            CHECK(model->NucRatesLogProb() == doctest::Approx(before).epsilon(1e-9));
            if (mcmc_nucrr(*model) != initial) { moved = true; }
        }
        CHECK(moved);
    }
    mcmc_cleanup_dataset();
}

TEST_CASE("mcmc: MoveCodonProfiles keeps profiles on the simplex and the state coherent") {
    mcmc_write_dataset();
    {
        auto model = mcmc_make_free_model();
        model->ResampleSub(1.0);
        model->CollectSitePathSuffStat();
        model->UpdateProfileOccupancies();

        std::vector<std::vector<double>> initial;
        for (int k = 0; k < model->GetNcat(); k++) {
            std::vector<double> p;
            for (int c = 0; c < 61; c++) { p.push_back(model->GetProfileCodon(k, c)); }
            initial.push_back(p);
        }

        bool moved = false;
        for (int rep = 0; rep < 10; rep++) {
            model->MoveCodonProfiles();
            mcmc_check_model_state(*model);
            double before = model->CodonLogPrior() + model->PathSuffStatLogProb();
            CHECK(std::isfinite(before));
            // recomputing the matrices from the (accepted or restored) profiles
            // must give the same path log probability
            model->UpdateMatrices();
            double after = model->CodonLogPrior() + model->PathSuffStatLogProb();
            CHECK(after == doctest::Approx(before).epsilon(1e-9));
            for (int k = 0; k < model->GetNcat(); k++) {
                for (int c = 0; c < 61; c++) {
                    if (model->GetProfileCodon(k, c) != initial[k][c]) { moved = true; }
                }
            }
        }
        CHECK(moved);
    }
    mcmc_cleanup_dataset();
}

TEST_CASE("mcmc: MoveOmegaMixture keeps omega finite and non-negative") {
    mcmc_write_dataset();
    {
        auto model = mcmc_make_free_model();
        model->ResampleSub(1.0);
        model->CollectSitePathSuffStat();
        for (int rep = 0; rep < 10; rep++) {
            model->CollectSiteOmegaPathSuffStat();
            model->MoveOmegaMixture(1);
            mcmc_check_model_state(*model);
            for (int k = 0; k < model->GetOmegaNcat(); k++) {
                CHECK(std::isfinite(model->GetComponentOmega(k)));
                CHECK(model->GetComponentOmega(k) >= 0.0);
            }
        }
    }
    mcmc_cleanup_dataset();
}

TEST_CASE("mcmc: full Move() cycles keep logprior and lnL finite") {
    mcmc_write_dataset();
    {
        auto model = mcmc_make_free_model();
        mcmc_check_model_state(*model);
        for (int it = 0; it < 5; it++) {
            model->Move();
            mcmc_check_model_state(*model);
            // joint log prob is the sum of its parts
            CHECK(model->GetLogProb() ==
                  doctest::Approx(model->GetLogPrior() + model->GetLogLikelihood()));
        }
    }
    mcmc_cleanup_dataset();
}

TEST_CASE("mcmc: clamped parameters are left untouched by Move()") {
    mcmc_write_dataset();
    {
        // flat nuc stat / rel rate => clamped; flat fitness => profiles clamped
        CodonMutSelMultipleOmegaModel model(mcmc_ali, mcmc_tree, "Null", "Null", "Null", 0, 1, 1,
            1, 0.0, true, "Null", true, true);
        std::vector<double> stat = mcmc_nucstat(model);
        std::vector<double> rr = mcmc_nucrr(model);
        std::vector<double> prof;
        for (int c = 0; c < 61; c++) { prof.push_back(model.GetProfileCodon(0, c)); }
        for (int it = 0; it < 3; it++) { model.Move(); }
        CHECK(mcmc_nucstat(model) == stat);
        CHECK(mcmc_nucrr(model) == rr);
        for (int c = 0; c < 61; c++) { CHECK(model.GetProfileCodon(0, c) == prof[c]); }
        CHECK(std::isfinite(model.GetLogPrior()));
        CHECK(std::isfinite(model.GetLogLikelihood()));
    }
    mcmc_cleanup_dataset();
}

// ---------------------------------------------------------------------------
// 2. Reader
// ---------------------------------------------------------------------------
TEST_CASE("reader: .param round trip (ChainDriver::fake_read + istream constructor)") {
    mcmc_write_dataset();
    {
        auto model = mcmc_make_free_model();
        for (int it = 0; it < 3; it++) { model->Move(); }

        // layout of a .param file: ChainDriver stamp (4 tokens), then the model
        std::stringstream param;
        param << "chain_name 1 -1 0\n";
        model->ToStream(param);

        ChainDriver::fake_read(param);
        CodonMutSelMultipleOmegaModel back(param);

        // stamp (everything the constructor needs to rebuild the model)
        CHECK(back.GetNsite() == model->GetNsite());
        CHECK(back.GetNcat() == model->GetNcat());
        CHECK(back.GetOmegaNcat() == model->GetOmegaNcat());
        CHECK(back.FlatFitness() == model->FlatFitness());
        CHECK(back.FlatNucStat() == model->FlatNucStat());
        CHECK(back.FlatNucRelRate() == model->FlatNucRelRate());

        // state (written with the default stream precision: 6 digits)
        for (int i = 0; i < Nnuc; i++) {
            CHECK(back.GetNucStat(i) == doctest::Approx(model->GetNucStat(i)).epsilon(1e-4));
        }
        for (int i = 0; i < Nrr; i++) {
            CHECK(back.GetNucRR(i) == doctest::Approx(model->GetNucRR(i)).epsilon(1e-4));
        }
        for (int k = 0; k < model->GetNcat(); k++) {
            for (int c = 0; c < 61; c++) {
                CHECK(back.GetProfileCodon(k, c) ==
                      doctest::Approx(model->GetProfileCodon(k, c)).epsilon(1e-4));
            }
        }
        for (int s = 0; s < model->GetNsite(); s++) {
            CHECK(back.GetProfileAlloc(s) == model->GetProfileAlloc(s));
            CHECK(back.GetOmegaAlloc(s) == model->GetOmegaAlloc(s));
            CHECK(back.GetSiteOmega(s) ==
                  doctest::Approx(model->GetSiteOmega(s)).epsilon(1e-4));
        }
        for (int k = 0; k < model->GetOmegaNcat(); k++) {
            CHECK(back.GetComponentOmega(k) ==
                  doctest::Approx(model->GetComponentOmega(k)).epsilon(1e-4));
        }
        CHECK(std::isfinite(back.GetLogPrior()));
        CHECK(std::isfinite(back.GetLogLikelihood()));
    }
    mcmc_cleanup_dataset();
}

TEST_CASE("reader: .chain points are read by ChainReader and are coherent") {
    mcmc_write_dataset();
    const std::string chain = "mutselc_mcmc_chain.chain";
    {
        auto model = mcmc_make_free_model();
        // write a 3-point chain (header + 3 lines) like ChainDriver does
        std::vector<std::vector<double>> nucstats;
        std::vector<std::vector<double>> omegas;
        {
            Tracer tracer(*model);
            std::ofstream os(chain);
            tracer.write_header(os);
            for (int it = 0; it < 3; it++) {
                model->Move();
                tracer.write_line(os);
                nucstats.push_back(mcmc_nucstat(*model));
                std::vector<double> o;
                for (int k = 0; k < model->GetOmegaNcat(); k++) {
                    o.push_back(model->GetComponentOmega(k));
                }
                omegas.push_back(o);
            }
        }

        // a fresh model, as the reader builds it from the .param stamp
        auto fresh = mcmc_make_free_model();
        ChainReader cr(*fresh, chain);
        for (int pt = 0; pt < 3; pt++) {
            cr.skip(1);
            fresh->UpdateMatrices();

            // what ReadCodonMutSelMultipleOmega computes at each point
            for (int i = 0; i < Nnuc; i++) {
                CHECK(fresh->GetNucStat(i) ==
                      doctest::Approx(nucstats[pt][i]).epsilon(1e-4));
            }
            for (int k = 0; k < fresh->GetOmegaNcat(); k++) {
                CHECK(fresh->GetComponentOmega(k) ==
                      doctest::Approx(omegas[pt][k]).epsilon(1e-4));
            }
            for (int s = 0; s < fresh->GetNsite(); s++) {
                // per-site omega is the omega of the site component
                CHECK(fresh->GetSiteOmega(s) ==
                      fresh->GetComponentOmega(fresh->GetOmegaAlloc(s)));
                // site-specific profile is the profile of the site component
                const std::vector<double> &p = fresh->GetProfile(s);
                REQUIRE(p.size() == 61);
                for (int c = 0; c < 61; c++) {
                    CHECK(p[c] == fresh->GetProfileCodon(fresh->GetProfileAlloc(s), c));
                }
                CHECK(std::isfinite(fresh->GetPredictedSiteOmegaKnot(s)));
            }
            for (int i = 0; i < Nnuc; i++) {
                for (int j = 0; j < Nnuc; j++) {
                    if (i != j) {
                        CHECK(std::isfinite(fresh->GetNucRate(i, j)));
                        CHECK(fresh->GetNucRate(i, j) > 0.0);
                    }
                }
            }
            CHECK(std::isfinite(fresh->GetPredictedRelativedS()));
            CHECK(std::isfinite(fresh->GetPredictedRelativedN()));
            CHECK(fresh->GetPredictedRelativedS() > 0.0);
            CHECK(fresh->GetPredictedRelativedN() > 0.0);
            // dN/dS is a mean of per-site ratios (scaled by omega), not the ratio of
            // the means, so it is not asserted equal to dN/dS here.
        }
    }
    std::remove(chain.c_str());
    mcmc_cleanup_dataset();
}

// The nucleotide-pair -> relative-rate index used by the reader (selection
// coefficient distribution) must follow the GTR convention AC AG AT CG CT GT.
TEST_CASE("reader: nucRRIndex formula matches the GTR relative-rate ordering") {
    std::vector<double> rr{1.0, 2.0, 3.0, 5.0, 7.0, 11.0};  // distinct, AC AG AT CG CT GT
    std::vector<double> stat{0.1, 0.2, 0.3, 0.4};
    GTRSubMatrixCodon nuc(4, rr, stat, false);
    nuc.UpdateMatrix();

    // expected pairs for the six indices
    const int pa[6] = {0, 0, 0, 1, 1, 2};
    const int pb[6] = {1, 2, 3, 2, 3, 3};
    for (int a = 0; a < Nnuc; a++) {
        for (int b = 0; b < Nnuc; b++) {
            if (a == b) continue;
            int lo = std::min(a, b), hi = std::max(a, b);
            // same expression as ReadCodonMutSelMultipleOmega.cpp
            int idx = (2 * Nnuc - lo - 1) * lo / 2 + hi - lo - 1;
            REQUIRE(idx >= 0);
            REQUIRE(idx < 6);
            CHECK(pa[idx] == lo);
            CHECK(pb[idx] == hi);
            // and the matrix agrees: Q(a,b) = rr[idx] * stat[b] * common scale
            double ref = nuc(0, 1) / (rr[0] * stat[1]);
            CHECK(nuc(a, b) == doctest::Approx(ref * rr[idx] * stat[b]).epsilon(1e-9));
        }
    }
    // the formula is still the one in the reader
    std::string src = mcmc_slurp(mcmc_src_dir() + "ReadCodonMutSelMultipleOmega.cpp");
    CHECK(src.find("(2 * Nnuc - nucFrom - 1) * nucFrom / 2 + nucTo - nucFrom - 1") !=
          std::string::npos);
    CHECK(src.find("(2 * Nnuc - nucTo - 1) * nucTo / 2 + nucFrom - nucTo - 1") !=
          std::string::npos);
}

// Parts of the reader that cannot be linked: lock them by scanning the source.
TEST_CASE("reader: source reads .param/.chain in the writer's order and validates input") {
    std::string src = mcmc_slurp(mcmc_src_dir() + "ReadCodonMutSelMultipleOmega.cpp");
    auto pos = [&](const char *s) { return src.find(s); };

    // .param: skip ChainDriver stamp, then build the model from the stream, then chain
    size_t p_param = pos("chain_name + \".param\"");
    size_t p_fake = pos("ChainDriver::fake_read");
    size_t p_model = pos("CodonMutSelMultipleOmegaModel model(is)");
    size_t p_chain = pos("ChainReader cr{model, chain_name + \".chain\"}");
    REQUIRE(p_param != std::string::npos);
    REQUIRE(p_fake != std::string::npos);
    REQUIRE(p_model != std::string::npos);
    REQUIRE(p_chain != std::string::npos);
    CHECK(p_param < p_fake);
    CHECK(p_fake < p_model);
    CHECK(p_model < p_chain);

    // the confidence interval is validated before use
    CHECK(pos("ci < 0.0 || ci > 1.0") != std::string::npos);
    // upper index is clamped inside the sample
    CHECK(pos("pt_up >= gene_omega.size()") != std::string::npos);
    CHECK(pos("pt_up >= omega[i].size()") != std::string::npos);
    // model-side matrices are refreshed at every point before derived stats
    CHECK(pos("model.UpdateMatrices()") != std::string::npos);
    // selection coefficients use the same fitness floor as the matrix
    CHECK(pos("kFitnessFloor") != std::string::npos);
}

// ---------------------------------------------------------------------------
// 3. Non-regression of the other models' matrices
// ---------------------------------------------------------------------------
namespace {

struct mcmc_Nuc {
    std::vector<double> rr{1.0, 2.0, 0.5, 1.5, 3.0, 0.7};
    std::vector<double> stat{0.1, 0.2, 0.3, 0.4};
    GTRSubMatrixCodon m{4, rr, stat, true};
    mcmc_Nuc() { m.UpdateMatrix(); }
};

}  // namespace

TEST_CASE("matrices: MGOmegaCodonSubMatrix rows, omega scaling, stationary") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    mcmc_Nuc nuc;
    for (double w : {0.0, 0.3, 1.0, 4.0}) {
        MGOmegaCodonSubMatrix m(&ss, &nuc.m, w);
        m.UpdateMatrix();
        double sumpi = 0;
        for (int i = 0; i < N; i++) {
            double out = 0;
            for (auto j : ss.GetNeighbors(i)) {
                int pos = ss.GetDifferingPosition(i, j);
                int a = ss.GetCodonPosition(pos, i);
                int b = ss.GetCodonPosition(pos, j);
                double expected = nuc.m(a, b) * (ss.Synonymous(i, j) ? 1.0 : w);
                CHECK(m(i, j) == doctest::Approx(expected).epsilon(1e-12).scale(1.0));
                CHECK(m(i, j) >= 0.0);
                out += m(i, j);
                // reversibility: pi_i Q_ij = pi_j Q_ji
                CHECK(m.Stationary(i) * m(i, j) ==
                      doctest::Approx(m.Stationary(j) * m(j, i)).epsilon(1e-9).scale(1e-3));
            }
            CHECK(m(i, i) == doctest::Approx(-out).epsilon(1e-12).scale(1.0));
            sumpi += m.Stationary(i);
            // stationary = product of nucleotide stationary frequencies, renormalised
        }
        CHECK(sumpi == doctest::Approx(1.0).epsilon(1e-9));
    }
}

TEST_CASE("matrices: AAMutSelOmegaCodonSubMatrix rows, omega, detailed balance") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    mcmc_Nuc nuc;

    // 20 amino-acid fitnesses summing to 1, kept alive (the matrix stores a reference)
    std::vector<double> aa(Naa, 0.0);
    double tot = 0;
    for (int a = 0; a < Naa; a++) {
        aa[a] = 1.0 + 0.3 * a;
        tot += aa[a];
    }
    for (auto &x : aa) { x /= tot; }

    for (double w : {0.2, 1.0, 3.0}) {
        AAMutSelOmegaCodonSubMatrix m(&ss, &nuc.m, aa, w, 1.0, false);
        m.CorruptMatrix();  // computes the fitnesses from aa and Ne
        m.UpdateMatrix();
        double sumpi = 0;
        for (int i = 0; i < N; i++) {
            double out = 0;
            for (auto j : ss.GetNeighbors(i)) {
                int pos = ss.GetDifferingPosition(i, j);
                int a = ss.GetCodonPosition(pos, i);
                int b = ss.GetCodonPosition(pos, j);
                double expected = nuc.m(a, b);
                if (!ss.Synonymous(i, j)) {
                    double ds = m.GetLogFitness(ss.Translation(j)) - m.GetLogFitness(ss.Translation(i));
                    expected *= w * (std::fabs(ds) < 1e-30 ? 1.0 : ds / (1.0 - std::exp(-ds)));
                }
                CHECK(m(i, j) == doctest::Approx(expected).epsilon(1e-9).scale(1.0));
                CHECK(m(i, j) >= 0.0);
                out += m(i, j);
                CHECK(m.Stationary(i) * m(i, j) ==
                      doctest::Approx(m.Stationary(j) * m(j, i)).epsilon(1e-6).scale(1e-3));
            }
            CHECK(m(i, i) == doctest::Approx(-out).epsilon(1e-12).scale(1.0));
            sumpi += m.Stationary(i);
        }
        CHECK(sumpi == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(m.GetOmega() == doctest::Approx(w));
        double dnds = m.GetPredictedDNDS();
        CHECK(std::isfinite(dnds));
        CHECK(dnds > 0.0);
    }
}

TEST_CASE("matrices: AA mut-sel with flat fitness reduces to MG omega") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    mcmc_Nuc nuc;
    std::vector<double> aa(Naa, 1.0 / Naa);
    const double w = 0.4;
    AAMutSelOmegaCodonSubMatrix aam(&ss, &nuc.m, aa, w, 1.0, false);
    aam.CorruptMatrix();
    aam.UpdateMatrix();
    MGOmegaCodonSubMatrix mg(&ss, &nuc.m, w);
    mg.UpdateMatrix();
    for (int i = 0; i < N; i++) {
        for (auto j : ss.GetNeighbors(i)) {
            CHECK(aam(i, j) == doctest::Approx(mg(i, j)).epsilon(1e-9).scale(1.0));
        }
    }
    // with flat fitness the predicted dN/dS is the pure mutational ratio, i.e.
    // independent of omega (omega is a multiplier applied afterwards)
    CHECK(aam.GetPredictedDNDS() == doctest::Approx(1.0).epsilon(1e-6));
}
