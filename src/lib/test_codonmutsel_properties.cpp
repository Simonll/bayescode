#include "doctest.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "CodonMutSelOmegaCodonSubMatrix.hpp"
#include "CodonStateSpace.hpp"
#include "GTRSubMatrix.hpp"

// ---------------------------------------------------------------------------
// Structural properties of the mut-sel codon matrix, plus the unit convention
// (neutral mutation flux = 1 substitution per CODON per unit branch length).
// The x3 bug in the reader came from violating that convention, so it is
// locked here at the matrix level (flux) and at the reader level (source scan).
// ---------------------------------------------------------------------------
namespace {

//! GTRSubMatrix exposing the protected GetRate()
class ExposedGTR : public GTRSubMatrix {
  public:
    ExposedGTR(int inNstate, const std::vector<double> &rr, const std::vector<double> &stat,
            bool innormalise = false)
            : SubMatrix(inNstate, innormalise), GTRSubMatrix(inNstate, rr, stat, innormalise) {}
    double Rate() const { return GetRate(); }
};

//! non-trivial, deterministic, strictly positive profile summing to 1
std::vector<double> make_profile(int n) {
    std::vector<double> f(n);
    double tot = 0;
    for (int i = 0; i < n; i++) {
        f[i] = 0.2 + std::fabs(std::sin(1.7 * i + 0.3)) * (1 + (i % 5));
        tot += f[i];
    }
    for (auto &x : f) x /= tot;
    return f;
}

const std::vector<double> nonuniform_stat{0.1, 0.2, 0.3, 0.4};
// relative rates in order AC AG AT CG CT GT
const std::vector<double> nonuniform_rr{0.5, 2.0, 1.3, 0.7, 3.1, 1.1};

}  // namespace

TEST_CASE("CodonMutSelOmegaCodonSubMatrix: detailed balance pi_i Q_ij == pi_j Q_ji") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    GTRSubMatrix nuc(4, nonuniform_rr, nonuniform_stat, true);
    auto fit = make_profile(N);
    CodonMutSelOmegaCodonSubMatrix m(&ss, &nuc, fit, 0.37, false);
    m.UpdateMatrix();

    int nonzero = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            double fwd = m.Stationary(i) * m(i, j);
            double bwd = m.Stationary(j) * m(j, i);
            CHECK(fwd == doctest::Approx(bwd).epsilon(1e-9));
            if (fwd > 0) nonzero++;
        }
    }
    CHECK(nonzero > 0);
}

TEST_CASE("CodonMutSelOmegaCodonSubMatrix: rows sum to zero, off-diagonals >= 0") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    GTRSubMatrix nuc(4, nonuniform_rr, nonuniform_stat, true);
    auto fit = make_profile(N);
    CodonMutSelOmegaCodonSubMatrix m(&ss, &nuc, fit, 2.5, false);
    m.UpdateMatrix();

    for (int i = 0; i < N; i++) {
        double out = 0;
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            CHECK(m(i, j) >= 0.0);
            out += m(i, j);
        }
        CHECK(m(i, i) <= 0.0);
        CHECK(out == doctest::Approx(-m(i, i)).epsilon(1e-12));
    }
}

TEST_CASE("CodonMutSelOmegaCodonSubMatrix: non-synonymous rates are linear in omega") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    auto fit = make_profile(N);
    GTRSubMatrix nuc1(4, nonuniform_rr, nonuniform_stat, true);
    GTRSubMatrix nuc2(4, nonuniform_rr, nonuniform_stat, true);
    GTRSubMatrix nuc7(4, nonuniform_rr, nonuniform_stat, true);
    CodonMutSelOmegaCodonSubMatrix m1(&ss, &nuc1, fit, 1.0, false);
    CodonMutSelOmegaCodonSubMatrix m2(&ss, &nuc2, fit, 2.0, false);
    CodonMutSelOmegaCodonSubMatrix m7(&ss, &nuc7, fit, 0.7, false);
    m1.UpdateMatrix();
    m2.UpdateMatrix();
    m7.UpdateMatrix();

    int nsyn = 0, nnonsyn = 0;
    for (int i = 0; i < N; i++) {
        for (auto j : ss.GetNeighbors(i)) {
            if (ss.Synonymous(i, j)) {
                nsyn++;
                CHECK(m2(i, j) == doctest::Approx(m1(i, j)).epsilon(1e-12));
                CHECK(m7(i, j) == doctest::Approx(m1(i, j)).epsilon(1e-12));
            } else {
                nnonsyn++;
                CHECK(m2(i, j) == doctest::Approx(2.0 * m1(i, j)).epsilon(1e-5));
                CHECK(m7(i, j) == doctest::Approx(0.7 * m1(i, j)).epsilon(1e-5));
            }
        }
    }
    CHECK(nsyn > 0);
    CHECK(nnonsyn > 0);
}

// Unit convention: the nucleotide GTR is normalised so that the neutral flux is
// 1 substitution per codon (3 positions x 1/3 per nt). With a flat profile and
// uniform nucstat the stationary law is uniform over the 61 sense codons, and
// the mean total rate equals (mean number of sense neighbours) / 9 exactly,
// i.e. just below 1 (only the moves towards stop codons are lost). A x3
// violation (flux ~3 per codon) is caught here.
TEST_CASE("CodonMutSelOmegaCodonSubMatrix: neutral flux is 1 substitution per codon") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    std::vector<double> rr(6, 1.0);
    std::vector<double> stat(4, 0.25);
    GTRSubMatrix nuc(4, rr, stat, true);
    std::vector<double> fit(N, 1.0 / N);  // flat profile
    CodonMutSelOmegaCodonSubMatrix m(&ss, &nuc, fit, 1.0, false);
    m.UpdateMatrix();

    double flux = 0, mean_neighbors = 0;
    for (int i = 0; i < N; i++) {
        flux += m.Stationary(i) * -m(i, i);
        mean_neighbors += static_cast<double>(ss.GetNeighbors(i).size()) / N;
    }
    CHECK(m.Stationary(0) == doctest::Approx(1.0 / N).epsilon(1e-6));
    CHECK(flux == doctest::Approx(mean_neighbors / 9.0).epsilon(1e-4));
    CHECK(flux < 1.0);
    CHECK(flux > 0.9);

    // each single-nucleotide change has rate exactly 1/9 under uniform GTR
    for (int i = 0; i < N; i++) {
        for (auto j : ss.GetNeighbors(i)) {
            CHECK(m(i, j) == doctest::Approx(1.0 / 9.0).epsilon(1e-6));
        }
    }
}

// Same convention with a non-uniform nucleotide process: neutral flux
// (mean of -Q_ii under the neutral stationary law) stays <= 1 and close to it.
TEST_CASE("CodonMutSelOmegaCodonSubMatrix: neutral flux ~1 per codon, non-uniform nucleotides") {
    CodonStateSpace ss(Universal);
    int N = ss.GetNstate();
    GTRSubMatrix nuc(4, nonuniform_rr, nonuniform_stat, true);
    std::vector<double> fit(N, 1.0 / N);
    CodonMutSelOmegaCodonSubMatrix m(&ss, &nuc, fit, 1.0, false);
    m.UpdateMatrix();
    double flux = 0;
    for (int i = 0; i < N; i++) flux += m.Stationary(i) * -m(i, i);
    CHECK(flux < 1.0);
    CHECK(flux > 0.8);
}

TEST_CASE("GTRSubMatrix::GetRate == 2 * 3 * sum_{i<j} pi_i pi_j rr_ij") {
    for (int variant = 0; variant < 2; variant++) {
        std::vector<double> rr = variant ? nonuniform_rr : std::vector<double>(6, 1.0);
        std::vector<double> stat = variant ? nonuniform_stat : std::vector<double>(4, 0.25);
        ExposedGTR g(4, rr, stat, false);

        double norm = 0;
        int k = 0;  // (0,1)(0,2)(0,3)(1,2)(1,3)(2,3)
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 4; j++) norm += stat[i] * stat[j] * rr[k++];
        }
        CHECK(g.Rate() == doctest::Approx(2 * norm * 3).epsilon(1e-12));
    }
    // uniform case in closed form: 2 * (6/16) * 3 = 2.25
    std::vector<double> rr(6, 1.0), stat(4, 0.25);
    ExposedGTR u(4, rr, stat, false);
    CHECK(u.Rate() == doctest::Approx(2.25).epsilon(1e-12));
}

// The reader is a main(), so it cannot be linked in. Lock the x3 fix by scanning
// its source: the --for_simulation branch must write GetBranchLength(node) as is
// (subst/codon), never multiplied by 3 (subst/nt).
TEST_CASE("ReadCodonMutSelMultipleOmega: --for_simulation writes branch lengths unscaled") {
    std::string here = __FILE__;
    std::string dir = here.substr(0, here.find_last_of("/\\") + 1);
    std::ifstream in(dir + "../ReadCodonMutSelMultipleOmega.cpp");
    REQUIRE(in.good());

    std::string line;
    int found = 0;
    while (std::getline(in, line)) {
        if (line.find("set_tag") != std::string::npos &&
            line.find("\"length\"") != std::string::npos) {
            found++;
            CHECK(line.find("to_string(model.GetBranchLength(node))") != std::string::npos);
            CHECK(line.find("3 *") == std::string::npos);
            CHECK(line.find("* 3") == std::string::npos);
            CHECK(line.find("3*") == std::string::npos);
        }
    }
    CHECK(found == 1);
}
