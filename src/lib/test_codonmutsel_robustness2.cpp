#include "doctest.h"

#include "CodonMutSelMultipleOmegaModel.hpp"

#include <sys/stat.h>
#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Second batch of robustness checks for CodonMutSelMultipleOmega (review
// findings G7, G9, G11, G17; G12 and G16 are covered in
// test_codonmutsel_robustness.cpp).
//
// Same convention as test_codonmutsel_robustness.cpp: ROBUST_BUG(expr) states
// the CORRECT behaviour that production code does not have yet. It is a WARN
// by default and a hard CHECK with MUTSELC_STRICT_ROBUSTNESS=1. The helpers of
// that file (anonymous namespace, same translation unit since all_tests.cpp
// includes both) are reused here, so this file must be included after it.
//
// G9 and G11 are checked on the source text: the faulty paths abort() (assert)
// or crash inside the parser, which would kill the whole test binary.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// G7 (behavioural): no NaN acceptance rate when a component is empty
// ---------------------------------------------------------------------------
TEST_CASE("G7: MoveBaseCenters returns a finite rate when no component is eligible") {
    write_robust_dataset();
    // baseNcat=-2 -> basemin=1: component 0 (the only occupied one) is excluded
    CodonMutSelMultipleOmegaModel model(tmp_ali_rb, tmp_tree_rb, "Null", "Null", "Null",
        0, 1, -2, 1, 0.0, true, "Null", true, true);
    double rate = model.MoveBaseCenters(1.0, 1);
    ROBUST_BUG(std::isfinite(rate));
    cleanup_robust_dataset();
}

TEST_CASE("G7: no unguarded 'nacc / ntot' remains, whatever its formatting (source)") {
    std::string src = read_src_file("CodonMutSelMultipleOmegaModel.hpp");
    REQUIRE(!src.empty());
    std::string sq = squash(src);
    int unguarded = 0;
    for (size_t p = sq.find("nacc/ntot"); p != std::string::npos; p = sq.find("nacc/ntot", p + 1)) {
        unguarded++;
    }
    ROBUST_BUG(unguarded == 0);
}

// ---------------------------------------------------------------------------
// G9: Ncat <= Nsite must be an explicit runtime check, not only an assert()
// (asserts vanish with -DNDEBUG, i.e. in Release).
// ---------------------------------------------------------------------------
TEST_CASE("G9: Ncat > Nsite is rejected by an explicit check, not only by assert (source)") {
    std::string src = read_src_file("CodonMutSelMultipleOmegaModel.hpp");
    REQUIRE(!src.empty());
    std::string sq = squash(src);

    // the constraint is still stated (sanity)
    REQUIRE((sq.find("Ncat<=Nsite") != std::string::npos ||
             sq.find("Ncat>Nsite") != std::string::npos));

    bool explicit_check = sq.find("if(Ncat>Nsite") != std::string::npos ||
                          sq.find("if(!(Ncat<=Nsite") != std::string::npos ||
                          sq.find("if(Nsite<Ncat") != std::string::npos ||
                          sq.find("if(!(Nsite>=Ncat") != std::string::npos;
    ROBUST_BUG(explicit_check);
}

TEST_CASE("G9: Ncat == Nsite is still accepted") {
    write_robust_dataset();
    // 12 codon sites: the boundary must remain valid
    CHECK_NOTHROW(CodonMutSelMultipleOmegaModel(tmp_ali_rb, tmp_tree_rb, "Null", "Null", "Null",
        0, 12, 1, 1, 0.0, true, "Null", true, true));
    cleanup_robust_dataset();
}

TEST_CASE("G9: profile allocation bound is not enforced by assert alone (source)") {
    std::string src = read_src_file("CodonMutSelMultipleOmegaModel.hpp");
    REQUIRE(!src.empty());
    std::string sq = squash(src);
    const std::string needle = "assert(alloc[alloc.size()-1]<fitness_profiles.size())";
    size_t p = sq.find(needle);
    if (p != std::string::npos) {
        // the assert, if kept, must be accompanied by an explicit test of the same bound
        std::string rest = sq.substr(0, p) + sq.substr(p + needle.size());
        bool also_explicit = rest.find("alloc.size()-1]>=fitness_profiles.size()") !=
                             std::string::npos;
        ROBUST_BUG(also_explicit);
    }
}

// ---------------------------------------------------------------------------
// G11: the tree file must be checked before parsing.
// ---------------------------------------------------------------------------
TEST_CASE("G11: a missing tree file is detected before parsing (source)") {
    std::string src = read_src_file("CodonMutSelMultipleOmegaModel.hpp");
    REQUIRE(!src.empty());
    std::string sq = squash(src);

    size_t open = sq.find("std::ifstreamtree_stream{treefile};");
    if (open == std::string::npos) { open = sq.find("std::ifstreamtree_stream(treefile);"); }
    REQUIRE(open != std::string::npos);
    size_t parse = sq.find("NHXParserparser{tree_stream};", open);
    if (parse == std::string::npos) { parse = sq.find("NHXParserparser(tree_stream);", open); }
    REQUIRE(parse != std::string::npos);

    std::string between = sq.substr(open, parse - open);
    bool checked = between.find("!tree_stream") != std::string::npos ||
                   between.find("tree_stream.is_open()") != std::string::npos ||
                   between.find("tree_stream.fail()") != std::string::npos ||
                   between.find("tree_stream.good()") != std::string::npos ||
                   between.find("if(tree_stream)") != std::string::npos;
    ROBUST_BUG(checked);
}

// ---------------------------------------------------------------------------
// G17: the model stamp is whitespace-delimited on read; a path with a space
// breaks the resume (the istream constructor reads paths with operator>>).
// ---------------------------------------------------------------------------
TEST_CASE("G17: the stamp survives data/tree paths containing a space") {
    const std::string dir = "mutselc robust dir";
    mkdir(dir.c_str(), 0755);
    const std::string ali = dir + "/ali.phy", tre = dir + "/tree.tre";
    write_robust_dataset();
    {
        std::ifstream a(tmp_ali_rb), t(tmp_tree_rb);
        std::ofstream(ali) << a.rdbuf();
        std::ofstream(tre) << t.rdbuf();
    }
    cleanup_robust_dataset();

    {
        CodonMutSelMultipleOmegaModel model(ali, tre, "Null", "Null", "Null",
            0, 1, 1, 1, 0.0, true, "Null", true, true);
        std::ostringstream stamp;
        model.GetModelStamp(stamp);

        // the stamp must still carry the paths verbatim (tab-separated fields)
        std::vector<std::string> fields = split_tabs(stamp.str());
        REQUIRE(fields.size() >= 3);
        CHECK(fields[1] == ali);
        CHECK(fields[2] == tre);

        // Emulate the current reader (operator>>). A correct round trip needs
        // a space-safe reader, in which case the plain `is >> datafile >>
        // treefile` disappears from the source.
        std::istringstream is(stamp.str());
        std::string name, datafile, treefile;
        is >> name >> datafile >> treefile;
        bool roundtrip = (datafile == ali && treefile == tre);

        std::string src = read_src_file("CodonMutSelMultipleOmegaModel.hpp");
        REQUIRE(!src.empty());
        bool plain_reader =
            squash(src).find("is>>datafile>>treefile>>profilesfile;") != std::string::npos;
        std::string unsafe = (plain_reader && !roundtrip) ? "plain operator>> reader" : "";
        ROBUST_BUG(unsafe == "");
    }

    std::remove(ali.c_str());
    std::remove(tre.c_str());
    rmdir(dir.c_str());
}
