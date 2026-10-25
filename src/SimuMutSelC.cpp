// Standalone sequence simulator for CodonMutSelMultipleOmega.
//
// Reads samples of a `.pvalues` file (readmutselcomega --for_simulation) and simulates
// codon sequences along the tree by a hand-written jump-chain (CTMC), drawing the rates
// from the BayesCode matrix CodonMutSelOmegaCodonSubMatrix itself (same matrix as the
// inference). PhyloProcess::PostPredSample is deliberately NOT used.
//
// .pvalues layout (see ReadCodonMutSelMultipleOmega.cpp, --for_simulation):
//   line 0 : tab-separated 13-token stamp (omegaNcat at index 5, Ncat at index 8)
//   per sample, Ncat+7 lines:
//     "ds \t dn"; newick tree with [&&NHX:length=...]; nucstat (4); nucrelrate (6);
//     Ncat lines of 61 codon fitnesses; Nsite profile allocations; omegaNcat omegas
//     (already including omega_shift); Nsite omega allocations.
//   Nsite is deduced from the number of tokens of the profile-allocation line.
//
// Note on gaps: the .pvalues does not store the original alignment mask, so simulated
// sequences are gap-free (every taxon has all Nsite codons).
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include "CodonMutSelOmegaCodonSubMatrix.hpp"
#include "CodonStateSpace.hpp"
#include "GTRSubMatrix.hpp"
#include "Random.hpp"
#include "tree/implem.hpp"

using namespace std;

static vector<double> parse_doubles(const string &line) {
    vector<double> v;
    istringstream is(line);
    string tok;
    while (is >> tok) { v.push_back(stod(tok)); }
    return v;
}

static string next_line(istream &is, const char *what) {
    string line;
    if (!getline(is, line)) {
        cerr << "error: .pvalues truncated while reading " << what << "\n";
        exit(1);
    }
    return line;
}

static void usage(const char *prog, ostream &os) {
    os << "usage: " << prog << " <file.pvalues> <out_prefix> [options]\n"
       << "\n"
       << "Simulates codon sequences under the mut-sel model (CodonMutSelOmegaCodonSubMatrix,\n"
       << "the BayesCode matrix) along the tree of each selected .pvalues sample (jump-chain CTMC).\n"
       << "\n"
       << "options:\n"
       << "  -b, --start N     first sample to simulate (0-based, default 0)\n"
       << "  -e, --every N     step between simulated samples (default 1)\n"
       << "  -u, --until N     last sample, inclusive (default = start, i.e. a single sample;\n"
       << "                    -1 = up to the end of the file)\n"
       << "  -n, --nrep N      independent replicates per sample (default 1)\n"
       << "  -s, --seed N      random seed (default -1); set once, samples and replicates\n"
       << "                    are consecutive draws of the same stream\n"
       << "      --normalise   build matrices with innormalise=true (flux 1 per unit of time)\n"
       << "      --fasta       also write <prefix>[_s<k>]_<rep>.fasta\n"
       << "  -h, --help        this help\n"
       << "\n"
       << "outputs (<tag> is empty if one sample is simulated, '_s<k>' (k = sample index) otherwise):\n"
       << "  <prefix><tag>_<rep>.ali    phylip alignment of replicate <rep> (0-based)\n"
       << "  <prefix><tag>_<rep>.fasta  same in fasta (with --fasta)\n"
       << "  <prefix><tag>.tree         newick tree of the sample\n"
       << "  <prefix><tag>.stats        per replicate Nsub, Nsynsub, Nnon (= Nsub - Nsynsub), plus\n"
       << "                             means, ratio over tree length and stationary flux\n"
       << "\n"
       << "note: the .pvalues has no original alignment mask, so simulated sequences have no gaps.\n";
}

struct Options {
    string pvalues, prefix;
    int start = 0, every = 1, until = -2;  // -2: not given (= start)
    int nrep = 1, seed = -1;
    bool normalise = false, fasta = false;
};

static Options parse_args(int argc, char **argv) {
    Options o;
    vector<string> pos;
    auto need = [&](int &i) -> string {
        if (i + 1 >= argc) {
            cerr << "error: option " << argv[i] << " needs a value\n";
            exit(1);
        }
        return argv[++i];
    };
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "-h" || a == "--help") {
            usage(argv[0], cout);
            exit(0);
        } else if (a == "-b" || a == "--start") {
            o.start = stoi(need(i));
        } else if (a == "-e" || a == "--every") {
            o.every = stoi(need(i));
        } else if (a == "-u" || a == "--until") {
            o.until = stoi(need(i));
        } else if (a == "-n" || a == "--nrep") {
            o.nrep = stoi(need(i));
        } else if (a == "-s" || a == "--seed") {
            o.seed = stoi(need(i));
        } else if (a == "--normalise") {
            o.normalise = true;
        } else if (a == "--fasta") {
            o.fasta = true;
        } else if (a.size() > 1 && a[0] == '-') {
            cerr << "error: unknown option " << a << "\n";
            usage(argv[0], cerr);
            exit(1);
        } else {
            pos.push_back(a);
        }
    }
    if (pos.size() != 2) {
        usage(argv[0], cerr);
        exit(1);
    }
    o.pvalues = pos[0];
    o.prefix = pos[1];
    if (o.until == -2) { o.until = o.start; }
    if (o.start < 0 || o.every < 1 || o.nrep < 1) {
        cerr << "error: need start >= 0, every >= 1, nrep >= 1\n";
        exit(1);
    }
    return o;
}

int main(int argc, char **argv) {
    Options opt = parse_args(argc, argv);
    Random::InitRandom(opt.seed);

    ifstream is(opt.pvalues);
    if (!is) {
        cerr << "error: cannot open " << opt.pvalues << "\n";
        return 1;
    }
    // stamp
    vector<string> stamp;
    {
        istringstream ss(next_line(is, "stamp"));
        string f;
        while (getline(ss, f, '\t')) { stamp.push_back(f); }
    }
    if (stamp.size() < 13) {
        cerr << "error: stamp has " << stamp.size() << " tokens, expected 13\n";
        return 1;
    }
    int omegaNcat = stoi(stamp[5]);
    int Ncat = stoi(stamp[8]);

    // samples start, start+every, ... up to `until` (inclusive; -1 = end of file)
    bool multi = opt.until < 0 || opt.until > opt.start;
    int nsimulated = 0;
    for (int sample = 0; opt.until < 0 || sample <= opt.until; sample++) {
        // read the Ncat+7 lines of this sample (EOF on its first line = end of the file)
        vector<string> lines(1);
        if (!getline(is, lines[0])) { break; }
        for (int l = 1; l < Ncat + 7; l++) { lines.push_back(next_line(is, "sample")); }
        if (sample < opt.start || (sample - opt.start) % opt.every != 0) { continue; }
        string tag = multi ? "_s" + to_string(sample) : "";

        int li = 0;
        vector<double> dsdn = parse_doubles(lines[li++]);
        string treeline = lines[li++];
        vector<double> nucstat = parse_doubles(lines[li++]);
        vector<double> nucrr = parse_doubles(lines[li++]);
        vector<vector<double>> profiles;
        for (int i = 0; i < Ncat; i++) { profiles.push_back(parse_doubles(lines[li++])); }
        vector<double> palloc_d = parse_doubles(lines[li++]);
        vector<double> omegas = parse_doubles(lines[li++]);
        vector<double> oalloc_d = parse_doubles(lines[li++]);
        int Nsite = palloc_d.size();
        if (nucstat.size() != 4 || nucrr.size() != 6 || (int)omegas.size() != omegaNcat ||
            (int)oalloc_d.size() != Nsite) {
            cerr << "error: inconsistent .pvalues sample " << sample << " (nucstat "
                 << nucstat.size() << ", nucrr " << nucrr.size() << ", omegas " << omegas.size()
                 << "/" << omegaNcat << ", omega alloc " << oalloc_d.size() << "/" << Nsite
                 << ")\n";
            return 1;
        }

        // tree (with branch lengths from the "length" tag)
        istringstream tree_stream(treeline);
        NHXParser parser{tree_stream};
        const AnnotatedTree &atree = parser.get_tree();
        unique_ptr<const Tree> tree = make_from_parser(parser);
        int Nnode = tree->nb_nodes();
        vector<double> blength(Nnode, 0.0);
        double totlength = 0;
        for (int n = 0; n < Nnode; n++) {
            if (!tree->is_root(n)) {
                blength[n] = stod(atree.tag(n, "length"));
                totlength += blength[n];
            }
        }
        vector<int> leaves;
        for (int n : tree->root_to_leaves_iter()) {
            if (tree->is_leaf(n)) { leaves.push_back(n); }
        }

        // BayesCode objects: genetic code, nucleotide matrix, one codon matrix per (profile, omega)
        CodonStateSpace statespace(Universal);
        int Nstate = statespace.GetNstate();
        GTRSubMatrixCodon nucmatrix(4, nucrr, nucstat, true);
        for (auto &p : profiles) {
            if ((int)p.size() != Nstate) {
                cerr << "error: profile of size " << p.size() << " != " << Nstate << "\n";
                return 1;
            }
        }
        map<pair<int, int>, unique_ptr<CodonMutSelOmegaCodonSubMatrix>> matrices;
        vector<const CodonMutSelOmegaCodonSubMatrix *> sitemat(Nsite);
        for (int i = 0; i < Nsite; i++) {
            int p = (int)palloc_d[i], o = (int)oalloc_d[i];
            if (p < 0 || p >= Ncat || o < 0 || o >= omegaNcat) {
                cerr << "error: allocation out of range at site " << i << "\n";
                return 1;
            }
            auto &m = matrices[make_pair(p, o)];
            if (!m) {
                m.reset(new CodonMutSelOmegaCodonSubMatrix(
                    &statespace, &nucmatrix, profiles[p], omegas[o], opt.normalise));
            }
            sitemat[i] = m.get();
        }

        // expected substitutions per site and unit of time under stationarity, for diagnostics
        double flux = 0;
        for (int i = 0; i < Nsite; i++) {
            const auto &m = *sitemat[i];
            for (int c = 0; c < Nstate; c++) { flux += m.Stationary(c) * -m(c, c); }
        }
        flux /= Nsite;

        ofstream ot(opt.prefix + tag + ".tree");
        ot << treeline << "\n";

        vector<long long> nsubs, nsyns;
        vector<vector<int>> state(Nnode, vector<int>(Nsite));
        for (int rep = 0; rep < opt.nrep; rep++) {
            long long nsub = 0, nsyn = 0;
            for (int site = 0; site < Nsite; site++) {
                const auto &m = *sitemat[site];
                state[tree->root()][site] =
                    Random::DrawFromDiscreteDistribution(m.GetStationary(), Nstate);
            }
            for (int node : tree->root_to_leaves_iter()) {
                if (tree->is_root(node)) { continue; }
                int par = tree->parent(node);
                for (int site = 0; site < Nsite; site++) {
                    const auto &m = *sitemat[site];
                    int cur = state[par][site];
                    double time = 0;
                    while (true) {
                        double total = -m(cur, cur);
                        if (total <= 0) { break; }
                        time += Random::sExpo() / total;
                        if (time >= blength[node]) { break; }
                        // draw the target codon by CDF inversion over the outgoing rates
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
                        if (statespace.Synonymous(cur, next)) { nsyn++; }
                        cur = next;
                        nsub++;
                    }
                    state[node][site] = cur;
                }
            }
            nsubs.push_back(nsub);
            nsyns.push_back(nsyn);
            cout << "sample " << sample << " rep " << rep << " Nsub " << nsub << " Nsynsub " << nsyn
                 << " Nnon " << nsub - nsyn << " Nsub_per_codon " << (double)nsub / Nsite << "\n";

            string base = opt.prefix + tag + "_" + to_string(rep);
            ofstream os(base + ".ali");
            os << leaves.size() << "\t" << 3 * Nsite << "\n";
            for (int leaf : leaves) {
                os << tree->node_name(leaf) << "\t";
                for (int site = 0; site < Nsite; site++) {
                    os << statespace.GetState(state[leaf][site]);
                }
                os << "\n";
            }
            if (opt.fasta) {
                ofstream of(base + ".fasta");
                for (int leaf : leaves) {
                    of << ">" << tree->node_name(leaf) << "\n";
                    for (int site = 0; site < Nsite; site++) {
                        of << statespace.GetState(state[leaf][site]);
                    }
                    of << "\n";
                }
            }
        }

        double sum_nsub = 0, sum_syn = 0;
        for (int r = 0; r < opt.nrep; r++) {
            sum_nsub += nsubs[r];
            sum_syn += nsyns[r];
        }
        double mean = sum_nsub / opt.nrep / Nsite;
        ofstream st(opt.prefix + tag + ".stats");
        st << setprecision(10) << "sample\t" << sample << "\nNsite\t" << Nsite << "\nNtaxa\t"
           << leaves.size() << "\nnrep\t" << opt.nrep << "\ntree_length\t" << totlength << "\n";
        for (int r = 0; r < opt.nrep; r++) {
            st << "rep\t" << r << "\tNsub\t" << nsubs[r] << "\tNsynsub\t" << nsyns[r] << "\tNnon\t"
               << nsubs[r] - nsyns[r] << "\n";
        }
        st << "mean_Nsub\t" << sum_nsub / opt.nrep << "\nmean_Nsynsub\t" << sum_syn / opt.nrep
           << "\nmean_Nnon\t" << (sum_nsub - sum_syn) / opt.nrep << "\nmean_Nsub_per_codon\t"
           << mean << "\nratio_over_tree_length\t" << mean / totlength << "\nstationary_flux\t"
           << flux << "\nexpected_Nsub_per_codon\t" << flux * totlength << "\n";
        cout << setprecision(8) << "sample " << sample << ": tree_length " << totlength
             << " mean_Nsub_per_codon " << mean << " ratio " << mean / totlength << " flux "
             << flux << " expected " << flux * totlength << "\n";
        nsimulated++;
    }
    if (nsimulated == 0) {
        cerr << "error: no sample simulated (range outside the file?)\n";
        return 1;
    }
    return 0;
}
