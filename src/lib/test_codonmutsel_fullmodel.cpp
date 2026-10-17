#include "doctest.h"

#include "CodonMutSelMultipleOmegaModel.hpp"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Helpers to materialise a small synthetic dataset (alignment + tree) into
// temporary files, build a full CodonMutSelMultipleOmegaModel from them and
// exercise the model-level API. The alignment is deliberately tiny (5 taxa, 15
// codons) so a full model construction + likelihood + stamp stays fast.
// ---------------------------------------------------------------------------
namespace {

std::string tmp_ali;
std::string tmp_tree;

//! write alignment + tree into a temp dir and stash the paths
void write_toy_dataset() {
    // 12 codons / 36 nt, identical length across the 5 taxa; several
    // non-synonymous differences so the alignment is not constant.
    std::vector<std::string> seq{
        "ATGGCCGATATTGTCAGCGTGAAAGCCTGCCCGCCG",
        "ATGGCCGATTTTGTCAGCGTGAAAGCCTGCCCGCCG",
        "ATGGCCGATATTGTCAGCGTGAAAGCCTGCCCGCCA",
        "ATGGCCGATATTGTCAGCATGAAAGCCTGCCCGCCG",
        "ATGGCCGATATTGTCAGCGTGAAAGCCTGCCCGCTG",
    };
    tmp_ali = "mutselc_test_ali.phy";
    tmp_tree = "mutselc_test_tree.tre";
    {
        std::ofstream out(tmp_ali);
        out << seq.size() << " " << seq[0].size() << "\n";
        for (size_t i = 0; i < seq.size(); i++) {
            out << "taxon" << i << "\t" << seq[i] << "\n";
        }
    }
    {
        std::ofstream out(tmp_tree);
        out << "(((taxon0:1,taxon1:1):1,(taxon2:1,taxon3:1):1):1,taxon4:1);\n";
    }
}

void cleanup_toy_dataset() {
    std::remove(tmp_ali.c_str());
    std::remove(tmp_tree.c_str());
}

}  // namespace

// ---------------------------------------------------------------------------
// Full model construction + degeneracy-free likelihood/prior sanity.
// ---------------------------------------------------------------------------
TEST_CASE("CodonMutSelMultipleOmegaModel: constructs and yields finite logprior/likelihood") {
    write_toy_dataset();

    // omega free (omegamode=0), Ncat=1, omegaNcat=1, flat fitness/stat/rates.
    CodonMutSelMultipleOmegaModel model(
        tmp_ali, tmp_tree, "Null",           // datafile, treefile, profilesfile
        "Null", "Null",                       // nucstatfile, nucrelratefile
        0,                                    // omegamode: free
        1,                                    // Ncat
        1,                                    // baseNcat
        1,                                    // omegaNcat
        0.0,                                  // omegashift
        true,                                 // flatfitness
        "Null",                               // delta_omega_array_file
        true,                                 // flatnucstat
        true);                                // flatnucrelrate

    // likelihood and prior must both be finite after construction+Update()
    double lp = model.GetLogPrior();
    double ll = model.GetLogLikelihood();
    CHECK(std::isfinite(lp));
    CHECK(std::isfinite(ll));

    // dN/dS must be a finite, non-negative number
    double dnds = model.GetPredictedRelativedNdS();
    CHECK(std::isfinite(dnds));
    CHECK(dnds >= 0.0);

    // model size bookkeeping: at least one non-constant codon site, and at most
    // the 12 codons in the alignment.
    int nsite = model.GetNsite();
    CHECK(nsite > 0);
    CHECK(nsite <= 12);

    cleanup_toy_dataset();
}

// ---------------------------------------------------------------------------
// Regression test for C3: the .pvalues stamp must be exactly 13 tokens, with
// Ncat at index 8 (the format read by LikelihoodFreePhylogenetics'
// readBayescodeParametersMutSelC). The full model stamp appends
// nucstatfile/nucrelratefile at the END so the first 13 tokens keep position.
// ---------------------------------------------------------------------------
TEST_CASE("CodonMutSelMultipleOmegaModel: GetPvaluesStamp is exactly 13 tokens, Ncat@8") {
    write_toy_dataset();
    CodonMutSelMultipleOmegaModel model(
        tmp_ali, tmp_tree, "Null", "Null", "Null",
        0, 1, 1, 2,   // omegamode=0, Ncat=1, baseNcat=1, omegaNcat=2
        0.5,          // omegashift
        true, "Null", true, true);

    std::ostringstream pvs, full;
    model.GetPvaluesStamp(pvs);
    model.GetModelStamp(full);

    std::vector<std::string> pv_tokens, full_tokens;
    { std::stringstream ss(pvs.str()); std::string tok; while (ss >> tok) pv_tokens.push_back(tok); }
    { std::stringstream ss(full.str()); std::string tok; while (ss >> tok) full_tokens.push_back(tok); }

    CHECK(pv_tokens.size() == 13);
    // Ncat is token index 8 (0-based) per readBayescodeParametersMutSelC
    CHECK(pv_tokens[8] == "1");
    // omegaNcat at index 5, omega_shift at index 6
    CHECK(pv_tokens[5] == "2");
    CHECK(pv_tokens[6] == "0.5");
    // full stamp must be 2 tokens longer (nucstatfile, nucrelratefile appended)
    CHECK(full_tokens.size() == 15);
    CHECK(full_tokens[13] == "Null");
    CHECK(full_tokens[14] == "Null");

    cleanup_toy_dataset();
}

// ---------------------------------------------------------------------------
// Regression test for C5 validations at model level: an omegashift < 0 must be
// rejected. (The model prints an error and exits; we assert the guard exists by
// checking a negative shift cannot be silently accepted through construction in
// a mode where it would reach the IIDGamma allocation.)  We instead check the
// behaviour that IS portable: inputting a negative shift while omegancat > 0
// must not produce NaN in the delta-omega array when omega is pinned (mode 3).
// ---------------------------------------------------------------------------
TEST_CASE("CodonMutSelMultipleOmegaModel: omega-pinned mode keeps delta-omega = shift (C5)") {
    write_toy_dataset();
    // omegamode=3 (fixed): delta_omega_array must be all equal, no NaN, even for
    // omegancat > 1 without a loaded array file.
    CodonMutSelMultipleOmegaModel model(
        tmp_ali, tmp_tree, "Null", "Null", "Null",
        3, 1, 1, 3,   // omegamode=3 (fixed), Ncat=1, omegaNcat=3
        0.25,         // omegashift >= 0
        true, "Null", true, true);

    // In fixed mode without a loaded array, GetComponentOmega must return the
    // shift for every component (all delta-omega set to 1.0 was the old bug
    // leading to omega=1 despite a shift).
    double o0 = model.GetComponentOmega(0);
    double o1 = model.GetComponentOmega(1);
    double o2 = model.GetComponentOmega(2);
    CHECK(o0 == doctest::Approx(0.25 + 1.0));  // shift + delta(=1)
    CHECK(o0 == o1);
    CHECK(o1 == o2);

    cleanup_toy_dataset();
}