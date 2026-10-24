#include "doctest.h"

#include "CodonMutSelMultipleOmegaModel.hpp"
#include "test_tmp_names.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Lot 6 QC: regression tests for the fixes of lots A (input validation),
// B (logprior of free allocations with clamped profiles) and C (precision of
// omega_shift in the checkpoint stamp).
//
// Same convention as test_codonmutsel_parsing.cpp, whose helpers (KNOWN_BUG,
// run_in_child, is_clean_rejection, TmpFile, read_file, source_dir, ...) are
// reused: this file must be included after it in all_tests.cpp. KNOWN_BUG is a
// WARN by default and a hard CHECK with MUTSELC_STRICT_PARSING=1; once the
// fixes are merged, replace KNOWN_BUG by CHECK.
// All helpers defined here are prefixed lot6_ (single translation unit).
// ---------------------------------------------------------------------------
namespace {

std::string lot6_ali = "mutselc_lot6_ali.phy";
std::string lot6_tree = "mutselc_lot6_tree.tre";

void lot6_write_dataset() {
    std::vector<std::string> seq{
        "ATGGCCGATATTGTCAGCGTGAAAGCCTGCCCGCCG",
        "ATGGCCGATTTTGTCAGCGTGAAAGCCTGCCCGCCG",
        "ATGGCCGATATTGTCAGCGTGAAAGCCTGCCCGCCA",
        "ATGGCCGATATTGTCAGCATGAAAGCCTGCCCGCCG",
        "ATGGCCGATATTGTCAGCGTGAAAGCCTGCCCGCTG",
    };
    {
        std::ofstream out(lot6_ali);
        out << seq.size() << " " << seq[0].size() << "\n";
        for (size_t i = 0; i < seq.size(); i++) { out << "taxon" << i << "\t" << seq[i] << "\n"; }
    }
    std::ofstream(lot6_tree) << "(((taxon0:1,taxon1:1):1,(taxon2:1,taxon3:1):1):1,taxon4:1);\n";
}

void lot6_cleanup_dataset() {
    std::remove(lot6_ali.c_str());
    std::remove(lot6_tree.c_str());
}

std::string lot6_squash(std::string const &s) {
    std::string out;
    for (char c : s) {
        if (!std::isspace(static_cast<unsigned char>(c))) { out += c; }
    }
    return out;
}

bool lot6_contains_any(std::string const &s, std::vector<std::string> const &needles) {
    for (auto const &n : needles) {
        if (s.find(n) != std::string::npos) { return true; }
    }
    return false;
}

//! status of open_delta_omega_array(file) run in a child
int lot6_delta_status(std::string const &file) {
    return run_in_child([&] {
        open_delta_omega_array(file);
        return true;
    });
}

//! a 61-column preference file (header + one line per row, rows of 61 values)
std::string lot6_pref_text(std::vector<std::vector<double>> const &rows) {
    std::ostringstream os;
    os.precision(17);
    for (int i = 0; i < 61; i++) { os << (i ? "," : "") << "c" << i; }
    os << "\n";
    for (auto const &r : rows) {
        for (size_t i = 0; i < r.size(); i++) { os << (i ? "," : "") << r[i]; }
        os << "\n";
    }
    return os.str();
}

std::vector<double> lot6_uniform_profile() { return std::vector<double>(61, 1.0 / 61); }

std::vector<double> lot6_skewed_profile() {
    std::vector<double> f(61);
    double tot = 0;
    for (int i = 0; i < 61; i++) {
        f[i] = 1.0 + (i % 7);
        tot += f[i];
    }
    for (auto &x : f) { x /= tot; }
    return f;
}

//! 61-value profile with entry idx = val; the others share the rest so that the sum is 1
std::vector<double> lot6_profile_with(int idx, double val) {
    std::vector<double> f(61, (1.0 - val) / 60);
    f[idx] = val;
    return f;
}

int lot6_pref61_status(std::string const &file) {
    return run_in_child([&] {
        open_preferences(file, 61);
        return true;
    });
}

//! text of the reader between two markers ("" if one is missing)
std::string lot6_reader_block(std::string const &from, std::string const &to) {
    std::string src = read_file(source_dir() + "../ReadCodonMutSelMultipleOmega.cpp");
    auto b = src.find(from);
    if (b == std::string::npos) { return ""; }
    auto e = src.find(to, b);
    if (e == std::string::npos) { return ""; }
    return src.substr(b, e - b);
}

std::vector<std::string> lot6_split_tabs(std::string const &s) {
    std::vector<std::string> out;
    std::string tok;
    std::istringstream is(s);
    while (std::getline(is, tok, '\t')) { out.push_back(tok); }
    return out;
}

}  // namespace

// ===========================================================================
// LOT A -- open_delta_omega_array
// ===========================================================================
TEST_CASE("lot6 A: open_delta_omega_array reads a well-formed file") {
    TmpFile f("lot6_dom_ok.tmp", "0.5\n1.5\n\n2.25\n");
    int s = run_in_child([&] {
        auto r = open_delta_omega_array(f.name);
        return r.size() == 3 && std::abs(r[0] - 0.5) < 1e-12 && std::abs(r[2] - 2.25) < 1e-12;
    });
    CHECK(s == kReturnedOk);
}

TEST_CASE("lot6 A: open_delta_omega_array rejects a non-numeric line cleanly") {
    TmpFile f("lot6_dom_abc.tmp", "0.5\nabc\n1.5\n");
    int s = lot6_delta_status(f.name);
    KNOWN_BUG(s != kRawStodException);
    KNOWN_BUG(is_clean_rejection(s));
}

TEST_CASE("lot6 A: open_delta_omega_array rejects a header line") {
    TmpFile f("lot6_dom_hdr.tmp", "delta_omega\n0.5\n1.5\n");
    KNOWN_BUG(is_clean_rejection(lot6_delta_status(f.name)));
}

TEST_CASE("lot6 A: open_delta_omega_array rejects trailing garbage after a number") {
    TmpFile f("lot6_dom_tail.tmp", "0.5\n1.5abc\n");
    KNOWN_BUG(is_clean_rejection(lot6_delta_status(f.name)));
}

TEST_CASE("lot6 A: open_delta_omega_array rejects nan and inf") {
    TmpFile n("lot6_dom_nan.tmp", "0.5\nnan\n");
    TmpFile i("lot6_dom_inf.tmp", "0.5\ninf\n");
    TmpFile m("lot6_dom_minf.tmp", "0.5\n-inf\n");
    KNOWN_BUG(is_clean_rejection(lot6_delta_status(n.name)));
    KNOWN_BUG(is_clean_rejection(lot6_delta_status(i.name)));
    KNOWN_BUG(is_clean_rejection(lot6_delta_status(m.name)));
}

TEST_CASE("lot6 A: open_delta_omega_array rejects zero and negative values (log(0) in the prior)") {
    TmpFile z("lot6_dom_zero.tmp", "0.5\n0\n");
    TmpFile n("lot6_dom_neg.tmp", "0.5\n-1.5\n");
    KNOWN_BUG(is_clean_rejection(lot6_delta_status(z.name)));
    KNOWN_BUG(is_clean_rejection(lot6_delta_status(n.name)));
}

TEST_CASE("lot6 A: open_delta_omega_array rejects an empty file (omegaNcat would be 0)") {
    TmpFile f("lot6_dom_empty.tmp", "");
    KNOWN_BUG(is_clean_rejection(lot6_delta_status(f.name)));
}

// ===========================================================================
// LOT A -- open_preferences: non-positive fitness
// ===========================================================================
TEST_CASE("lot6 A: open_preferences reads a well-formed 61-column file") {
    TmpFile f("lot6_pref_ok.tmp",
        lot6_pref_text({lot6_uniform_profile(), lot6_skewed_profile(), lot6_uniform_profile()}));
    int s = run_in_child([&] {
        auto r = open_preferences(f.name, 61);
        return std::get<0>(r).size() == 2 && std::get<1>(r) == std::vector<size_t>({0, 1, 0});
    });
    CHECK(s == kReturnedOk);
}

TEST_CASE("lot6 A: open_preferences rejects a zero fitness (profile still sums to 1)") {
    TmpFile f("lot6_pref_zero.tmp", lot6_pref_text({lot6_profile_with(5, 0.0)}));
    KNOWN_BUG(is_clean_rejection(lot6_pref61_status(f.name)));
}

TEST_CASE("lot6 A: open_preferences rejects a negative fitness (profile still sums to 1)") {
    TmpFile f("lot6_pref_neg.tmp", lot6_pref_text({lot6_profile_with(5, -0.5)}));
    KNOWN_BUG(is_clean_rejection(lot6_pref61_status(f.name)));
}

TEST_CASE("lot6 A: open_preferences rejects nan and inf fitness") {
    std::string good = lot6_pref_text({lot6_uniform_profile()});
    // replace the first value of the data line
    auto bad = [&](std::string const &tok) {
        auto nl = good.find('\n');
        auto comma = good.find(',', nl);
        return good.substr(0, nl + 1) + tok + good.substr(comma);
    };
    TmpFile n("lot6_pref_nan.tmp", bad("nan"));
    TmpFile i("lot6_pref_inf.tmp", bad("inf"));
    KNOWN_BUG(is_clean_rejection(lot6_pref61_status(n.name)));
    KNOWN_BUG(is_clean_rejection(lot6_pref61_status(i.name)));
}

TEST_CASE("lot6 A: open_preferences rejects a header-only file") {
    TmpFile f("lot6_pref_hdr.tmp", lot6_pref_text({}));
    KNOWN_BUG(is_clean_rejection(lot6_pref61_status(f.name)));
}

// ===========================================================================
// LOT A -- --ci / -c option of the reader (a main(): source scan, as for M4)
// ===========================================================================
TEST_CASE("lot6 A: --ci block is found in the reader source") {
    std::string blk =
        lot6_reader_block("confidence_interval.getValue().empty()", "string filename");
    CHECK(blk.find("confidence_interval") != std::string::npos);
    CHECK(blk.find("[0,1]") != std::string::npos);
}

TEST_CASE("lot6 A: --ci value is parsed without an unhandled stod exception ('-c abc')") {
    std::string blk = lot6_squash(
        lot6_reader_block("confidence_interval.getValue().empty()", "string filename"));
    REQUIRE(!blk.empty());
    // either a try/catch around stod, or a strict parser / strtod
    KNOWN_BUG(lot6_contains_any(blk, {"catch", "strtod", "parse_number", "from_chars"}));
}

TEST_CASE("lot6 A: --ci range check rejects NaN ('-c nan')") {
    std::string blk = lot6_squash(
        lot6_reader_block("confidence_interval.getValue().empty()", "string filename"));
    REQUIRE(!blk.empty());
    // 'ci < 0 || ci > 1' is false for NaN: need isfinite/isnan or a negated '>=' test
    KNOWN_BUG(lot6_contains_any(blk, {"isfinite", "isnan", "!(ci>=0", "!(ci>0"}));
}

TEST_CASE("lot6 A: the number of points to read is checked (size == 0 must not abort)") {
    std::string blk = lot6_squash(
        lot6_reader_block("read_args.GetSize()", "confidence_interval.getValue().empty()"));
    REQUIRE(!blk.empty());
    KNOWN_BUG(lot6_contains_any(blk, {"size<=0", "size<1", "size==0", "size>0"}));
}

// ===========================================================================
// LOT B -- logprior includes the allocation terms when profiles are clamped
// but the allocation is free (profiles file with fewer profiles than sites)
// ===========================================================================
TEST_CASE("lot6 B: logprior reacts to the profile weights when profiles are clamped and the "
          "allocation is free") {
    lot6_write_dataset();
    TmpFile prefs("lot6_prefs2.tmp",
        lot6_pref_text({lot6_uniform_profile(), lot6_skewed_profile()}));  // 2 profiles < 12 sites
    {
        // omegaNcat = 1: every omega term is constant, only the profile terms can move
        CodonMutSelMultipleOmegaModel model(lot6_ali, lot6_tree, prefs.name, "Null", "Null",
            0, 1, 1, 1, 0.0, false, "Null", true, true);

        // spread the sites over the two profiles, then redraw the weights given the occupancies
        for (int i = 0; i < 5; i++) {
            model.ResampleProfileAlloc();
            model.UpdateProfileOccupancies();
        }
        double lp0 = model.GetLogPrior();
        REQUIRE(std::isfinite(lp0));

        bool changed = false;
        for (int i = 0; i < 10 && !changed; i++) {
            model.ResampleProfileWeights();  // allocations, profiles and omega untouched
            double lp = model.GetLogPrior();
            REQUIRE(std::isfinite(lp));
            if (std::fabs(lp - lp0) > 1e-9) { changed = true; }
        }
        KNOWN_BUG(changed);
    }
    lot6_cleanup_dataset();
}

TEST_CASE("lot6 B: logprior is finite for a model with clamped profiles and free allocation") {
    lot6_write_dataset();
    TmpFile prefs("lot6_prefs2b.tmp",
        lot6_pref_text({lot6_uniform_profile(), lot6_skewed_profile()}));
    {
        CodonMutSelMultipleOmegaModel model(lot6_ali, lot6_tree, prefs.name, "Null", "Null",
            0, 1, 1, 1, 0.0, false, "Null", true, true);
        CHECK(std::isfinite(model.GetLogPrior()));
        CHECK(std::isfinite(model.GetLogProb()));
    }
    lot6_cleanup_dataset();
}

// ===========================================================================
// LOT C -- omega_shift survives the stamp round trip
// ===========================================================================
TEST_CASE("lot6 C: GetModelStamp writes omega_shift with enough digits") {
    lot6_write_dataset();
    {
        CodonMutSelMultipleOmegaModel model(lot6_ali, lot6_tree, "Null", "Null", "Null",
            0, 1, 1, 1, 0.123456789, true, "Null", true, true);
        std::ostringstream stamp;
        model.GetModelStamp(stamp);
        auto fields = lot6_split_tabs(stamp.str());
        REQUIRE(fields.size() > 6);
        // fields: name data tree profiles deltafile omegaNcat omega_shift ...
        KNOWN_BUG(std::abs(std::stod(fields[6]) - 0.123456789) < 1e-12);
    }
    lot6_cleanup_dataset();
}

TEST_CASE("lot6 C: GetPvaluesStamp writes omega_shift with enough digits") {
    lot6_write_dataset();
    {
        CodonMutSelMultipleOmegaModel model(lot6_ali, lot6_tree, "Null", "Null", "Null",
            0, 1, 1, 1, 0.123456789, true, "Null", true, true);
        std::ostringstream stamp;
        model.GetPvaluesStamp(stamp);
        auto fields = lot6_split_tabs(stamp.str());
        REQUIRE(fields.size() > 6);
        KNOWN_BUG(std::abs(std::stod(fields[6]) - 0.123456789) < 1e-12);
    }
    lot6_cleanup_dataset();
}

TEST_CASE("lot6 C: omega_shift round-trips through ToStream / istream constructor") {
    lot6_write_dataset();
    {
        CodonMutSelMultipleOmegaModel model(lot6_ali, lot6_tree, "Null", "Null", "Null",
            0, 1, 1, 1, 0.123456789, true, "Null", true, true);
        std::ostringstream os;
        model.ToStream(os);
        std::istringstream is(os.str());
        CodonMutSelMultipleOmegaModel resumed(is);
        double shift = resumed.GetComponentOmega(0) - resumed.GetDeltaOmega(0);
        KNOWN_BUG(std::abs(shift - 0.123456789) < 1e-12);
        // the stamp must be a fixed point: writing it again gives the same omega_shift field
        std::ostringstream a, b;
        model.GetModelStamp(a);
        resumed.GetModelStamp(b);
        KNOWN_BUG(lot6_split_tabs(a.str())[6] == lot6_split_tabs(b.str())[6]);
    }
    lot6_cleanup_dataset();
}

TEST_CASE("lot6 C: round shifts (0, 0.5, 2) still round-trip exactly") {
    lot6_write_dataset();
    {
        for (double shift : {0.0, 0.5, 2.0}) {
            CodonMutSelMultipleOmegaModel model(lot6_ali, lot6_tree, "Null", "Null", "Null",
                0, 1, 1, 1, shift, true, "Null", true, true);
            std::ostringstream os;
            model.ToStream(os);
            std::istringstream is(os.str());
            CodonMutSelMultipleOmegaModel resumed(is);
            CHECK(std::abs(resumed.GetComponentOmega(0) - resumed.GetDeltaOmega(0) - shift) <
                  1e-12);
        }
    }
    lot6_cleanup_dataset();
}
