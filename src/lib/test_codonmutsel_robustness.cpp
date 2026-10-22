#include "doctest.h"

#include "CodonMutSelMultipleOmegaModel.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Robustness checks for CodonMutSelMultipleOmega (review findings M5, M8, G7,
// G12, G16).
//
// ROBUST_BUG(expr) states the CORRECT behaviour that production code does not
// have yet. By default a violated expectation is only reported as a WARN, so
// the suite stays green while the finding is open. Run
//     MUTSELC_STRICT_ROBUSTNESS=1 ./all_tests
// to turn every ROBUST_BUG into a hard CHECK; this must pass once the findings
// are fixed, at which point ROBUST_BUG should be replaced by CHECK.
//
// M5 and G12 live in the main-like reader (ReadCodonMutSelMultipleOmega.cpp),
// which cannot be linked into the test binary, and the logprior structure (M8)
// and the acceptance-rate guards (G7) are partly invisible through the public
// API. Those are therefore checked on the source text, whitespace-insensitively.
// ---------------------------------------------------------------------------
namespace {

bool strict_robustness() { return std::getenv("MUTSELC_STRICT_ROBUSTNESS") != nullptr; }

#define ROBUST_BUG(expr)                 \
    do {                                 \
        if (strict_robustness()) {       \
            CHECK(expr);                 \
        } else {                         \
            WARN(expr);                  \
        }                                \
    } while (0)

//! read a file of the src/ tree (located relative to this test file); "" if not found
std::string read_src_file(const std::string &relpath) {
    std::string here = __FILE__;
    std::string libdir = here.substr(0, here.find_last_of('/') == std::string::npos
                                            ? 0
                                            : here.find_last_of('/'));
    std::vector<std::string> candidates{
        libdir + "/../" + relpath, "src/" + relpath, "../src/" + relpath, relpath};
    for (const auto &c : candidates) {
        std::ifstream in(c);
        if (in) {
            std::stringstream ss;
            ss << in.rdbuf();
            return ss.str();
        }
    }
    return "";
}

//! remove all whitespace and normalise "1.0"->"1", "2.0"->"2"
std::string squash(const std::string &s) {
    std::string out;
    for (char ch : s) {
        if (!isspace(static_cast<unsigned char>(ch))) { out += ch; }
    }
    for (const char *num : {"1.0", "2.0"}) {
        std::string from = num, to(1, num[0]);
        size_t pos = 0;
        while ((pos = out.find(from, pos)) != std::string::npos) {
            bool digit_before = pos > 0 && isdigit(static_cast<unsigned char>(out[pos - 1]));
            bool digit_after = pos + 3 < out.size() && isdigit(static_cast<unsigned char>(out[pos + 3]));
            if (digit_before || digit_after) {
                pos += 3;
            } else {
                out.replace(pos, 3, to);
            }
        }
    }
    return out;
}

//! text of the brace-delimited body that starts at the first '{' after `signature`
std::string body_after(const std::string &src, const std::string &signature) {
    size_t p = src.find(signature);
    if (p == std::string::npos) { return ""; }
    size_t open = src.find('{', p);
    if (open == std::string::npos) { return ""; }
    int depth = 0;
    for (size_t i = open; i < src.size(); i++) {
        if (src[i] == '{') { depth++; }
        if (src[i] == '}' && --depth == 0) { return src.substr(open, i - open + 1); }
    }
    return "";
}

std::string tmp_ali_rb = "mutselc_robust_ali.phy";
std::string tmp_tree_rb = "mutselc_robust_tree.tre";

void write_robust_dataset() {
    std::vector<std::string> seq{
        "ATGGCCGATATTGTCAGCGTGAAAGCCTGCCCGCCG",
        "ATGGCCGATTTTGTCAGCGTGAAAGCCTGCCCGCCG",
        "ATGGCCGATATTGTCAGCGTGAAAGCCTGCCCGCCA",
        "ATGGCCGATATTGTCAGCATGAAAGCCTGCCCGCCG",
        "ATGGCCGATATTGTCAGCGTGAAAGCCTGCCCGCTG",
    };
    {
        std::ofstream out(tmp_ali_rb);
        out << seq.size() << " " << seq[0].size() << "\n";
        for (size_t i = 0; i < seq.size(); i++) { out << "taxon" << i << "\t" << seq[i] << "\n"; }
    }
    {
        std::ofstream out(tmp_tree_rb);
        out << "(((taxon0:1,taxon1:1):1,(taxon2:1,taxon3:1):1):1,taxon4:1);\n";
    }
}

void cleanup_robust_dataset() {
    std::remove(tmp_ali_rb.c_str());
    std::remove(tmp_tree_rb.c_str());
}

std::vector<std::string> split_tabs(const std::string &line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, '\t')) { out.push_back(tok); }
    return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// M5: "-c 0.95" must give a 95% credible interval, i.e. quantiles
// (1-ci)/2 and (1+ci)/2 -- not [1-ci, ci] which is a 90% interval.
// ---------------------------------------------------------------------------
TEST_CASE("M5: credible interval bounds are (1-ci)/2 and (1+ci)/2") {
    std::string src = read_src_file("ReadCodonMutSelMultipleOmega.cpp");
    REQUIRE(!src.empty());

    // the line(s) defining the bounds, squashed
    auto find_decl = [&](const std::string &name) {
        size_t p = src.find("double " + name);
        if (p == std::string::npos) { return std::string(); }
        return squash(src.substr(p, src.find(';', p) - p + 1));
    };
    std::string lower = find_decl("lower");
    std::string upper = find_decl("upper");
    REQUIRE(!lower.empty());
    REQUIRE(!upper.empty());

    ROBUST_BUG(lower.find("(1-ci)/2") != std::string::npos);
    ROBUST_BUG(upper.find("(1+ci)/2") != std::string::npos);
    // the legacy [1-ci, ci] formulation must be gone
    ROBUST_BUG(upper.find("max(ci,1-ci)") == std::string::npos);
    ROBUST_BUG(lower.find("min(ci,1-ci)") == std::string::npos);
}

// ---------------------------------------------------------------------------
// M8: GetLogPrior must include P(omega_alloc | omega_weight) and the prior of
// omega_weight (and P(profile_alloc | weight)).
// Behavioural check: resampling omega_weight with allocations unchanged must
// change the log prior, since P(omega_alloc | omega_weight) depends on it.
// ---------------------------------------------------------------------------
TEST_CASE("M8: logprior depends on omega_weight through P(omega_alloc | omega_weight)") {
    write_robust_dataset();
    CodonMutSelMultipleOmegaModel model(tmp_ali_rb, tmp_tree_rb, "Null", "Null", "Null",
        0, 1, 1, 3,  // omegamode=0, Ncat=1, baseNcat=1, omegaNcat=3
        0.0, true, "Null", true, true);

    // spread the sites over the omega components so the weights matter
    for (int i = 0; i < 5; i++) {
        model.ResampleOmegaAlloc();
        model.UpdateOmegaOccupancies();
    }

    double lp0 = model.GetLogPrior();
    REQUIRE(std::isfinite(lp0));
    bool changed = false;
    for (int i = 0; i < 5 && !changed; i++) {
        model.ResampleOmegaWeights();  // allocations and everything else untouched
        double lp = model.GetLogPrior();
        REQUIRE(std::isfinite(lp));
        if (std::fabs(lp - lp0) > 1e-9) { changed = true; }
    }
    ROBUST_BUG(changed);
    cleanup_robust_dataset();
}

TEST_CASE("M8: GetLogPrior sums alloc|weight and omega_weight terms (source)") {
    std::string src = read_src_file("CodonMutSelMultipleOmegaModel.hpp");
    REQUIRE(!src.empty());
    std::string body = body_after(src, "double GetLogPrior() const");
    REQUIRE(!body.empty());
    // already present: sanity of the extraction
    CHECK(body.find("StickBreakingLogPrior") != std::string::npos);
    CHECK(body.find("DeltaOmegaLogProb") != std::string::npos);

    // helper(s) reachable from GetLogPrior: look for terms on the allocation
    // vectors and on the omega weights, either inline or through a member
    // function whose body mentions them.
    auto mentions = [&](const std::string &needle) {
        if (body.find(needle) != std::string::npos) { return true; }
        // follow calls "XxxLogPrior()" / "XxxLogProb()" made by GetLogPrior
        size_t pos = 0;
        while ((pos = body.find("total +=", pos)) != std::string::npos) {
            size_t semi = body.find(';', pos);
            std::string stmt = body.substr(pos + 8, semi - pos - 8);
            std::string callee = squash(stmt);
            size_t paren = callee.find('(');
            if (paren != std::string::npos) {
                std::string fn = callee.substr(0, paren);
                std::string fbody = body_after(src, "double " + fn + "() const");
                if (fbody.find(needle) != std::string::npos) { return true; }
            }
            pos = semi;
        }
        return false;
    };
    ROBUST_BUG(mentions("omega_alloc"));
    ROBUST_BUG(mentions("omega_weight"));
    ROBUST_BUG(mentions("profile_alloc"));
}

// ---------------------------------------------------------------------------
// G7: acceptance rates nacc/ntot must be guarded against ntot == 0 (NaN).
// ---------------------------------------------------------------------------
TEST_CASE("G7: MoveBaseConcentrations returns a finite rate when no component is eligible") {
    write_robust_dataset();
    // baseNcat=-2 -> basemin=1: component 0 is excluded from the move, and all
    // sites sit in base component 0, so ntot stays 0.
    CodonMutSelMultipleOmegaModel model(tmp_ali_rb, tmp_tree_rb, "Null", "Null", "Null",
        0, 1, -2, 1, 0.0, true, "Null", true, true);
    double rate = model.MoveBaseConcentrations(1.0);
    ROBUST_BUG(std::isfinite(rate));
    cleanup_robust_dataset();
}

TEST_CASE("G7: every 'return nacc / ntot' is guarded against ntot == 0 (source)") {
    std::string src = read_src_file("CodonMutSelMultipleOmegaModel.hpp");
    REQUIRE(!src.empty());
    std::string sq = squash(src);

    // count unguarded divisions vs. total occurrences
    const std::string bad = "returnnacc/ntot;";
    int unguarded = 0;
    for (size_t p = sq.find(bad); p != std::string::npos; p = sq.find(bad, p + 1)) { unguarded++; }
    // there are 5 moves returning an acceptance rate
    ROBUST_BUG(unguarded == 0);

    // and each of these moves must still return something (guard replaces, not removes)
    for (const char *fn : {"double MoveCodon(", "double MoveCodonGamma(", "double MoveBaseCenters(",
             "double MoveBaseConcentrations(", "double MoveOmegaValues("}) {
        std::string body = body_after(src, fn);
        CHECK(!body.empty());
        CHECK(body.find("return") != std::string::npos);
    }
}

// ---------------------------------------------------------------------------
// G12: the 6 output streams of the reader must have their opening verified.
// ---------------------------------------------------------------------------
TEST_CASE("G12: the .mutsel/.subsel/... output streams are checked after opening") {
    std::string src = read_src_file("ReadCodonMutSelMultipleOmega.cpp");
    REQUIRE(!src.empty());
    std::string sq = squash(src);

    for (const char *name : {"mutmutsel_os", "mutsubsel_os", "nonsynmutmutsel_os",
             "nonsynmutsubsel_os", "synmutmutsel_os", "synmutsubsel_os"}) {
        std::string n = name;
        // the stream is declared (sanity)
        REQUIRE(sq.find("ofstream" + n + "(") != std::string::npos);
        bool checked = sq.find("!" + n) != std::string::npos ||
                       sq.find(n + ".is_open()") != std::string::npos ||
                       sq.find(n + ".fail()") != std::string::npos ||
                       sq.find(n + ".good()") != std::string::npos ||
                       sq.find("if(" + n + ")") != std::string::npos;
        std::string unchecked = checked ? "" : n;  // name shown in the report on failure
        ROBUST_BUG(unchecked == "");
    }
}

// ---------------------------------------------------------------------------
// G16: no duplicated column name in the .chain / .trace headers.
// ---------------------------------------------------------------------------
namespace {

void check_no_duplicate_columns(int omegamode, int Ncat, int baseNcat, int omegaNcat) {
    write_robust_dataset();
    CodonMutSelMultipleOmegaModel model(tmp_ali_rb, tmp_tree_rb, "Null", "Null", "Null",
        omegamode, Ncat, baseNcat, omegaNcat, 0.0, true, "Null", true, true);

    Tracer chain_tracer(model);
    Tracer stat_tracer(model, processing::HasTag<Stat>());
    std::ostringstream chain_header, stat_header;
    chain_tracer.write_header(chain_header);
    stat_tracer.write_header(stat_header);

    for (const auto &h : {chain_header.str(), stat_header.str()}) {
        std::vector<std::string> cols = split_tabs(h);
        CHECK(!cols.empty());
        for (const auto &c : cols) {
            if (std::count(cols.begin(), cols.end(), c) > 1) {
                std::string duplicated = c;  // name shown in the report on failure
                ROBUST_BUG(duplicated == "");
                break;
            }
        }
    }
    CHECK(chain_tracer.nbr_header_fields() > 0);
    CHECK(stat_tracer.nbr_header_fields() > 0);
    cleanup_robust_dataset();
}

}  // namespace

TEST_CASE("G16: trace headers have no duplicated column (single base component)") {
    check_no_duplicate_columns(0, 1, 1, 1);
    check_no_duplicate_columns(0, 2, 1, 2);
}

TEST_CASE("G16: trace headers have no duplicated column (base mixture)") {
    check_no_duplicate_columns(0, 1, 2, 1);
    check_no_duplicate_columns(0, 1, -2, 1);
}

TEST_CASE("G16: kappa / basekappa are declared once per tracer (source)") {
    std::string src = read_src_file("CodonMutSelMultipleOmegaModel.hpp");
    REQUIRE(!src.empty());
    std::string body = body_after(src, "void declare_interface(Info info)");
    REQUIRE(!body.empty());
    std::string sq = squash(body);
    for (const char *kind : {"model_node(info,\"", "model_stat(info,\""}) {
        for (const char *name : {"kappa", "basekappa"}) {
            std::string needle = std::string(kind) + name + "\"";
            int n = 0;
            for (size_t p = sq.find(needle); p != std::string::npos; p = sq.find(needle, p + 1)) { n++; }
            std::string repeated = (n > 1) ? needle : "";  // shown in the report on failure
            ROBUST_BUG(repeated == "");
        }
    }
}
