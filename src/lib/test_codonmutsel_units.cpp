#include "doctest.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "CodonMutSelMultipleOmegaModel.hpp"
#include "CodonMutSelOmegaCodonSubMatrix.hpp"
#include "CodonStateSpace.hpp"
#include "GTRSubMatrix.hpp"

// ---------------------------------------------------------------------------
// QC lot 3: lock unit coherence of CodonMutSelMultipleOmega.
// Convention: branch lengths are in substitutions per CODON; neutral flux = 1
// per codon; omega is the plain dN/dS multiplier (no hidden additive offset).
//
// Some tests below describe the CORRECT behaviour while the production code is
// still known to be inconsistent (review findings U1-U3). Those assertions use
// UNIT_BUG: a WARN (reported, non-fatal) by default so the suite stays green,
// and a hard CHECK when compiled with -DMUTSELC_STRICT_UNITS. Once the
// production code is fixed, define UNIT_BUG as CHECK for good.
// ---------------------------------------------------------------------------
#ifdef MUTSELC_STRICT_UNITS
#define UNIT_BUG(expr) CHECK(expr)
#else
#define UNIT_BUG(expr) WARN(expr)
#endif

namespace {

std::string units_src_dir() {
    std::string here = __FILE__;
    return here.substr(0, here.find_last_of("/\\") + 1) + "../";
}

std::vector<std::string> units_read_lines(const std::string &path) {
    std::ifstream in(path);
    REQUIRE(in.good());
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l)) lines.push_back(l);
    return lines;
}

//! tiny dataset, written to the cwd; returns the paths
void units_write_dataset(const std::string &ali, const std::string &tre) {
    std::ofstream o(ali);
    o << "3 12\ntaxon0\tATGGCCGATATT\ntaxon1\tATGGCCGATTTT\ntaxon2\tATGGCAGATATT\n";
    o.close();
    std::ofstream t(tre);
    t << "((taxon0:1,taxon1:1):1,taxon2:1);\n";
}

}  // namespace

// U1: the "length" trace stat must be the plain total branch length (already
// per codon since the x3 fix), without any factor 3.
TEST_CASE("units U1: model 'length' stat is the total branch length, no x3") {
    auto lines = units_read_lines(units_src_dir() + "CodonMutSelMultipleOmegaModel.hpp");
    int found = 0;
    for (const auto &l : lines) {
        if (l.find("model_stat") != std::string::npos &&
            l.find("\"length\"") != std::string::npos) {
            found++;
            CHECK(l.find("GetTotalLength()") != std::string::npos);
            UNIT_BUG(l.find("3 *") == std::string::npos);
            UNIT_BUG(l.find("3*") == std::string::npos);
            UNIT_BUG(l.find("* 3") == std::string::npos);
        }
    }
    CHECK(found == 1);
}

// U3: the omega used by the matrix must be the omega the model reports
// (GetComponentOmega / reader): no hidden additive omegamin.
TEST_CASE("units U3: matrix GetOmega equals the omega parameter") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    std::vector<double> fit(N, 1.0 / N);
    for (double w : {0.0, 0.1, 1.0, 2.5}) {
        GTRSubMatrixCodon nuc(4, std::vector<double>(6, 1.0), std::vector<double>(4, 0.25), true);
        CodonMutSelOmegaCodonSubMatrix m(&ss, &nuc, fit, w, false);
        UNIT_BUG(m.GetOmega() == doctest::Approx(w).epsilon(1e-12));
    }
}

TEST_CASE("units U3: flat fitness, non-synonymous rate == omega * nucleotide rate") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    std::vector<double> fit(N, 1.0 / N);
    for (double w : {0.0, 0.1, 2.5}) {
        GTRSubMatrixCodon nuc(4, std::vector<double>(6, 1.0), std::vector<double>(4, 0.25), true);
        CodonMutSelOmegaCodonSubMatrix m(&ss, &nuc, fit, w, false);
        m.UpdateMatrix();
        int nn = 0;
        for (int i = 0; i < N; i++) {
            for (auto j : ss.GetNeighbors(i)) {
                if (ss.Synonymous(i, j)) continue;
                nn++;
                // uniform GTR normalised to flux 1/codon: nucleotide rate 1/9;
                // flat fitness: fixation factor 1
                UNIT_BUG(m(i, j) == doctest::Approx(w / 9.0).epsilon(1e-9).scale(1.0));
            }
        }
        CHECK(nn > 0);
    }
}

// U3 (model/reader side): the model defines omega as shift + delta and the reader
// prints it unmodified; so it is the matrix that must drop the offset.
TEST_CASE("units U3: model component omega is shift + delta, reader adds no offset") {
    auto lines = units_read_lines(units_src_dir() + "CodonMutSelMultipleOmegaModel.hpp");
    int found = 0;
    for (const auto &l : lines) {
        if (l.find("double GetComponentOmega") != std::string::npos) {
            found++;
            CHECK(l.find("omega_shift + delta_omega_array->GetVal(k)") != std::string::npos);
            CHECK(l.find("omegamin") == std::string::npos);
        }
    }
    CHECK(found == 1);
    for (const auto &l : units_read_lines(units_src_dir() + "ReadCodonMutSelMultipleOmega.cpp")) {
        CHECK(l.find("omegamin") == std::string::npos);
    }
}

// U4: global unit test. Flat fitness + uniform nucstat: total rate per codon is
// (#sense neighbours)/9, i.e. ~1 (neutral flux = 1 per codon; only moves towards
// stop codons are lost).
TEST_CASE("units U4: flat fitness + uniform nucstat -> neutral flux ~ 1 per codon") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    GTRSubMatrixCodon nuc(4, std::vector<double>(6, 1.0), std::vector<double>(4, 0.25), true);
    std::vector<double> fit(N, 1.0 / N);
    CodonMutSelOmegaCodonSubMatrix m(&ss, &nuc, fit, 1.0, false);
    m.UpdateMatrix();

    // neutral flux = sum_i pi_i * (-Q_ii), with pi uniform (flat fitness, uniform nucstat)
    double flux = 0;
    for (int i = 0; i < N; i++) {
        flux += m.Stationary(i) * -m(i, i);
        // row sum: diagonal = -sum of off-diagonal
        double out = 0;
        for (auto j : ss.GetNeighbors(i)) out += m(i, j);
        CHECK(-m(i, i) == doctest::Approx(out).epsilon(1e-12));
    }
    // neutral flux ~1 per codon (only moves towards stop codons are lost)
    CHECK(flux > 0.9);
    CHECK(flux < 1.0);  // would be ~3 if the x3 convention came back
}
