#include "doctest.h"

#include "CodonMutSelMultipleOmegaModel.hpp"
#include "test_tmp_names.hpp"

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Parsing of the input files of CodonMutSelMultipleOmega (profiles, nucstat,
// nucrelrate) and numerical precision of the --for_simulation output.
//
// open_* report fatal errors with exit(1), so every malformed-input scenario is
// executed in a forked child and judged on its exit status.
//
// KNOWN_BUG(expr) states the CORRECT behaviour that production code does not
// have yet (review findings M1, M3, M4). By default a violated expectation is
// only reported as a WARN, so the suite stays green while the bug is open. Run
//     MUTSELC_STRICT_PARSING=1 ./all_tests
// to turn every KNOWN_BUG into a hard CHECK; this must pass once M1/M3/M4 are
// fixed, at which point KNOWN_BUG should be replaced by CHECK.
// (the vendored doctest 1.1.4 has no should_fail decorator)
// ---------------------------------------------------------------------------
namespace {

bool strict_parsing() { return std::getenv("MUTSELC_STRICT_PARSING") != nullptr; }

#define KNOWN_BUG(expr)                  \
    do {                                 \
        if (strict_parsing()) {          \
            CHECK(expr);                 \
        } else {                         \
            WARN(expr);                  \
        }                                \
    } while (0)

enum ChildStatus {
    kReturnedOk = 0,        // function returned and the callback accepted the result
    kExit1 = 1,             // clean diagnostic + exit(1)
    kRawStodException = 2,  // std::invalid_argument / out_of_range leaked from stod
    kOtherException = 3,    // any other exception (e.g. a deliberate runtime_error)
    kWrongResult = 4,       // function returned but the result is wrong
    kOtherFailure = 99
};

//! run f in a child process (stderr silenced); f returns true when the result is correct
int run_in_child(const std::function<bool()> &f) {
    fflush(nullptr);
    pid_t pid = fork();
    if (pid == 0) {
        int devnull = open("/dev/null", O_WRONLY);
        if (devnull >= 0) { dup2(devnull, 2); }
        try {
            _exit(f() ? kReturnedOk : kWrongResult);
        } catch (std::invalid_argument const &) {
            _exit(kRawStodException);
        } catch (std::out_of_range const &) {
            _exit(kRawStodException);
        } catch (...) { _exit(kOtherException); }
    }
    int status = 0;
    waitpid(pid, &status, 0);
    return WIFEXITED(status) ? WEXITSTATUS(status) : kOtherFailure;
}

//! a deliberate, readable rejection of the input
bool is_clean_rejection(int s) { return s == kExit1 || s == kOtherException; }

//! rejected cleanly OR parsed to exactly the right value; never silently wrong
bool is_rejected_or_correct(int s) { return s == kReturnedOk || is_clean_rejection(s); }

struct TmpFile {
    std::string name;
    TmpFile(std::string const &n, std::string const &content) : name(unique_tmp_name(n)) {
        std::ofstream out(name);
        out << content;
    }
    ~TmpFile() { std::remove(name.c_str()); }
};

bool near(double a, double b) { return std::abs(a - b) < 1e-12; }

bool same(std::vector<double> const &a, std::vector<double> const &b) {
    if (a.size() != b.size()) { return false; }
    for (size_t i = 0; i < a.size(); i++) {
        if (!near(a[i], b[i])) { return false; }
    }
    return true;
}

bool profiles_are(std::string const &file, int n, std::vector<std::vector<double>> const &expected,
        std::vector<size_t> const &expected_alloc) {
    auto res = open_preferences(file, n);
    auto const &p = std::get<0>(res);
    if (p.size() != expected.size() || std::get<1>(res) != expected_alloc) { return false; }
    for (size_t i = 0; i < p.size(); i++) {
        if (!same(p[i], expected[i])) { return false; }
    }
    return true;
}

//! status of open_preferences(file, 3) run in a child
int pref_status(std::string const &file) {
    return run_in_child([&] {
        open_preferences(file, 3);
        return true;
    });
}

int nucstat_status(std::string const &file) {
    return run_in_child([&] {
        open_nucstat(file, 4);
        return true;
    });
}

int nucrelrate_status(std::string const &file) {
    return run_in_child([&] {
        open_nucrelrate(file, 6);
        return true;
    });
}

std::string read_file(std::string const &path) {
    std::ifstream in(path);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

std::string source_dir() {
    std::string here = __FILE__;
    return here.substr(0, here.find_last_of("/\\") + 1);
}

}  // namespace

// ===========================================================================
// M1 -- open_preferences
// ===========================================================================
TEST_CASE("open_preferences: well-formed file, each separator, dedup of identical profiles") {
    std::vector<std::vector<double>> prof{{0.2, 0.3, 0.5}, {0.6, 0.3, 0.1}};
    std::vector<size_t> alloc{0, 1, 0};
    TmpFile sp("mutselc_pref_sp.tmp", "a b c\n0.2 0.3 0.5\n0.6 0.3 0.1\n0.2 0.3 0.5\n");
    TmpFile co("mutselc_pref_co.tmp", "a,b,c\n0.2,0.3,0.5\n0.6,0.3,0.1\n0.2,0.3,0.5\n");
    TmpFile ta("mutselc_pref_ta.tmp", "a\tb\tc\n0.2\t0.3\t0.5\n0.6\t0.3\t0.1\n0.2\t0.3\t0.5\n");
    CHECK(profiles_are(sp.name, 3, prof, alloc));
    CHECK(profiles_are(co.name, 3, prof, alloc));
    CHECK(profiles_are(ta.name, 3, prof, alloc));
}

TEST_CASE("open_preferences: leading extra columns (site index) are skipped, blank lines ignored") {
    TmpFile f("mutselc_pref_lead.tmp", "site,a,b,c\n1,0.2,0.3,0.5\n\n2,0.6,0.3,0.1\n");
    CHECK(profiles_are(f.name, 3, {{0.2, 0.3, 0.5}, {0.6, 0.3, 0.1}}, {0, 1}));
}

TEST_CASE("open_preferences: header with fewer columns than Nstate is rejected") {
    TmpFile f("mutselc_pref_hdr.tmp", "a,b\n0.5,0.5\n");
    CHECK(is_clean_rejection(pref_status(f.name)));
}

TEST_CASE("open_preferences: row with more fitness columns than Nstate is rejected") {
    TmpFile f("mutselc_pref_long.tmp", "a,b,c\n0.2,0.3,0.3,0.2\n");
    CHECK(is_clean_rejection(pref_status(f.name)));
}

TEST_CASE("open_preferences: trailing separator on rows only never shifts the profile") {
    TmpFile f("mutselc_pref_trow.tmp", "a,b,c\n0.2,0.3,0.5,\n");
    int s = run_in_child([&] { return profiles_are(f.name, 3, {{0.2, 0.3, 0.5}}, {0}); });
    CHECK(is_rejected_or_correct(s));
}

TEST_CASE("open_preferences: trailing separator in header only does not shift the profile") {
    // M1
    TmpFile f("mutselc_pref_thdr.tmp", "a,b,c,\n0.2,0.3,0.5\n");
    int s = run_in_child([&] { return profiles_are(f.name, 3, {{0.2, 0.3, 0.5}}, {0}); });
    KNOWN_BUG(is_rejected_or_correct(s));
}

TEST_CASE("open_preferences: trailing separator in header and rows does not shift the profile") {
    // M1
    TmpFile f("mutselc_pref_tboth.tmp", "a,b,c,\n0.2,0.3,0.5,\n");
    int s = run_in_child([&] { return profiles_are(f.name, 3, {{0.2, 0.3, 0.5}}, {0}); });
    KNOWN_BUG(is_rejected_or_correct(s));
}

TEST_CASE("open_preferences: a row with fewer than Nstate columns is rejected") {
    // M1
    TmpFile f("mutselc_pref_short.tmp", "a,b,c\n0.2,0.3,0.5\n0.4,0.6\n");
    KNOWN_BUG(is_clean_rejection(pref_status(f.name)));
}

TEST_CASE("open_preferences: empty field 'a,,b' gives a clean error, not a raw stod exception") {
    // M1
    TmpFile f("mutselc_pref_empty.tmp", "a,b,c\n0.2,,0.8\n");
    int s = pref_status(f.name);
    KNOWN_BUG(s != kRawStodException);
    KNOWN_BUG(is_clean_rejection(s));
}

TEST_CASE("open_preferences: non-numeric field gives a clean error, not a raw stod exception") {
    // M1
    TmpFile f("mutselc_pref_nan.tmp", "a,b,c\n0.2,x,0.8\n");
    KNOWN_BUG(is_clean_rejection(pref_status(f.name)));
}

TEST_CASE("open_preferences: a profile not summing to 1 is rejected") {
    // M1 (only a warning today)
    TmpFile f("mutselc_pref_sum.tmp", "a,b,c\n0.5,0.5,0.5\n");
    KNOWN_BUG(is_clean_rejection(pref_status(f.name)));
}

// ===========================================================================
// M3 -- open_nucstat / open_nucrelrate
// ===========================================================================
TEST_CASE("open_nucstat / open_nucrelrate: well-formed files, each separator") {
    TmpFile ns("mutselc_ns_ok.tmp", "A C G T\n0.1 0.2 0.3 0.4\n");
    TmpFile nsc("mutselc_ns_okc.tmp", "A,C,G,T\n0.1,0.2,0.3,0.4\n");
    TmpFile rr("mutselc_rr_ok.tmp", "AC\tAG\tAT\tCG\tCT\tGT\n0.1\t0.2\t0.1\t0.2\t0.2\t0.2\n");
    CHECK(same(open_nucstat(ns.name, 4), {0.1, 0.2, 0.3, 0.4}));
    CHECK(same(open_nucstat(nsc.name, 4), {0.1, 0.2, 0.3, 0.4}));
    CHECK(same(open_nucrelrate(rr.name, 6), {0.1, 0.2, 0.1, 0.2, 0.2, 0.2}));
}

TEST_CASE("open_nucstat: leading label column is skipped") {
    TmpFile ns("mutselc_ns_lead.tmp", "id,A,C,G,T\n1,0.1,0.2,0.3,0.4\n");
    CHECK(same(open_nucstat(ns.name, 4), {0.1, 0.2, 0.3, 0.4}));
}

TEST_CASE("open_nucstat / open_nucrelrate: missing file or too-short header is rejected") {
    CHECK(is_clean_rejection(nucstat_status("mutselc_does_not_exist.tmp")));
    CHECK(is_clean_rejection(nucrelrate_status("mutselc_does_not_exist.tmp")));
    TmpFile ns("mutselc_ns_hdr.tmp", "A C G\n0.2 0.3 0.5\n");
    CHECK(is_clean_rejection(nucstat_status(ns.name)));
}

TEST_CASE("open_nucstat: row with too many columns is rejected") {
    TmpFile ns("mutselc_ns_long.tmp", "A C G T\n0.1 0.2 0.3 0.2 0.2\n");
    CHECK(is_clean_rejection(nucstat_status(ns.name)));
}

TEST_CASE("open_nucstat / open_nucrelrate: empty file is rejected") {
    TmpFile ns("mutselc_ns_empty.tmp", "");
    TmpFile rr("mutselc_rr_empty.tmp", "");
    CHECK(is_clean_rejection(nucstat_status(ns.name)));
    CHECK(is_clean_rejection(nucrelrate_status(rr.name)));
}

TEST_CASE("open_nucstat: header only (no data row) is rejected, not a zero vector") {
    // M3
    TmpFile f("mutselc_ns_nodata.tmp", "A C G T\n");
    KNOWN_BUG(is_clean_rejection(nucstat_status(f.name)));
}

TEST_CASE("open_nucstat: file without header (data line only) is rejected") {
    // M3
    TmpFile f("mutselc_ns_nohdr.tmp", "0.1 0.2 0.3 0.4\n");
    KNOWN_BUG(is_clean_rejection(nucstat_status(f.name)));
}

TEST_CASE("open_nucrelrate: header only (no data row) is rejected, not a zero vector") {
    // M3
    TmpFile f("mutselc_rr_nodata.tmp", "AC AG AT CG CT GT\n");
    KNOWN_BUG(is_clean_rejection(nucrelrate_status(f.name)));
}

TEST_CASE("open_nucrelrate: file without header (data line only) is rejected") {
    // M3
    TmpFile f("mutselc_rr_nohdr.tmp", "0.1 0.2 0.1 0.2 0.2 0.2\n");
    KNOWN_BUG(is_clean_rejection(nucrelrate_status(f.name)));
}

TEST_CASE("open_nucstat: data row with fewer than 4 values is rejected") {
    // M3
    TmpFile f("mutselc_ns_short.tmp", "A C G T\n0.5 0.5\n");
    KNOWN_BUG(is_clean_rejection(nucstat_status(f.name)));
}

TEST_CASE("open_nucrelrate: data row with fewer than 6 values is rejected") {
    // M3
    TmpFile f("mutselc_rr_short.tmp", "AC AG AT CG CT GT\n0.5 0.5\n");
    KNOWN_BUG(is_clean_rejection(nucrelrate_status(f.name)));
}

TEST_CASE("open_nucstat: sum different from 1 is rejected") {
    // M3 (only a warning today)
    TmpFile f("mutselc_ns_sum.tmp", "A C G T\n0.5 0.5 0.5 0.5\n");
    KNOWN_BUG(is_clean_rejection(nucstat_status(f.name)));
}

TEST_CASE("open_nucrelrate: sum different from 1 is rejected") {
    // M3 (only a warning today)
    TmpFile f("mutselc_rr_sum.tmp", "AC AG AT CG CT GT\n0.5 0.5 0.5 0.5 0.5 0.5\n");
    KNOWN_BUG(is_clean_rejection(nucrelrate_status(f.name)));
}

TEST_CASE("open_nucstat: a zero frequency is rejected (would give NaN downstream)") {
    // M3
    TmpFile f("mutselc_ns_zero.tmp", "A C G T\n0.0 0.3 0.3 0.4\n");
    KNOWN_BUG(is_clean_rejection(nucstat_status(f.name)));
}

TEST_CASE("open_nucstat / open_nucrelrate: negative values are rejected") {
    // M3
    TmpFile ns("mutselc_ns_neg.tmp", "A C G T\n-0.1 0.4 0.4 0.3\n");
    TmpFile rr("mutselc_rr_neg.tmp", "AC AG AT CG CT GT\n-0.1 0.3 0.2 0.2 0.2 0.2\n");
    KNOWN_BUG(is_clean_rejection(nucstat_status(ns.name)));
    KNOWN_BUG(is_clean_rejection(nucrelrate_status(rr.name)));
}

TEST_CASE("open_nucstat: empty field gives a clean error, not a raw stod exception") {
    // M3 (same stod issue as M1); double space -> empty field
    TmpFile f("mutselc_ns_stod.tmp", "A C G T\n0.1  0.5 0.4\n");
    KNOWN_BUG(is_clean_rejection(nucstat_status(f.name)));
}

// ===========================================================================
// M4 -- precision of the --for_simulation output. ReadCodonMutSelMultipleOmega.cpp
// is a main(), so (as for the x3 fix) the source of the branch is scanned.
// ===========================================================================
namespace {

//! text of the `prepare_files_for_simulation` branch of the reader
std::string for_simulation_block() {
    std::string src = read_file(source_dir() + "../ReadCodonMutSelMultipleOmega.cpp");
    auto b = src.find("prepare_files_for_simulation.getValue()");
    auto e = src.find("selection_coef_distribution.getValue()", b);
    if (b == std::string::npos || e == std::string::npos) { return ""; }
    return src.substr(b, e - b);
}

//! setprecision(N) with N >= 10, max_digits10, or a printf-style %.Ng with N >= 10
bool has_high_precision(std::string const &s) {
    for (size_t pos = s.find("setprecision("); pos != std::string::npos;
            pos = s.find("setprecision(", pos + 1)) {
        if (std::atoi(s.c_str() + pos + 13) >= 10) { return true; }
    }
    if (s.find("max_digits10") != std::string::npos) { return true; }
    for (size_t pos = s.find("%."); pos != std::string::npos; pos = s.find("%.", pos + 1)) {
        if (std::atoi(s.c_str() + pos + 2) >= 10) { return true; }
    }
    return false;
}

}  // namespace

TEST_CASE("std::to_string(double) loses small branch lengths (why M4 matters)") {
    CHECK(std::to_string(4e-7) == "0.000000");
    CHECK(std::to_string(0.123456789) == "0.123457");
}

TEST_CASE("--for_simulation: block is found in the reader source") {
    CHECK(for_simulation_block().find("GetBranchLength") != std::string::npos);
    CHECK(for_simulation_block().find("GetProfileCodon") != std::string::npos);
}

TEST_CASE("--for_simulation: branch lengths are not written with std::to_string") {
    // M4
    std::string blk = for_simulation_block();
    REQUIRE(!blk.empty());
    for (size_t pos = blk.find("GetBranchLength"); pos != std::string::npos;
            pos = blk.find("GetBranchLength", pos + 1)) {
        size_t line_begin = blk.rfind('\n', pos);
        std::string line = blk.substr(line_begin + 1, blk.find('\n', pos) - line_begin - 1);
        KNOWN_BUG(line.find("to_string(") == std::string::npos);
    }
}

TEST_CASE("--for_simulation: numbers are written with at least 10 significant digits") {
    // M4: setprecision(17) or %.10g
    KNOWN_BUG(has_high_precision(for_simulation_block()));
}
