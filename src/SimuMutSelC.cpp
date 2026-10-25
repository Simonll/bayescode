// Standalone sequence simulator for CodonMutSelMultipleOmega.
//
// Reads one sample of a `.pvalues` file (readmutselcomega --for_simulation) and simulates
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

int main(int argc, char **argv) {
    if (argc < 3) {
        cerr << "usage: " << argv[0]
             << " <file.pvalues> <out_prefix> [sample_index=0] [seed=-1] [nrep=1] [normalise=0]\n";
        return 1;
    }
    string pvalues = argv[1], prefix = argv[2];
    int sample = argc > 3 ? stoi(argv[3]) : 0;
    int seed = argc > 4 ? stoi(argv[4]) : -1;
    int nrep = argc > 5 ? stoi(argv[5]) : 1;
    // normalise=1 builds the codon matrices with innormalise=true (flux 1 per unit of time)
    bool normalise = argc > 6 && stoi(argv[6]) != 0;
    Random::InitRandom(seed);

    ifstream is(pvalues);
    if (!is) {
        cerr << "error: cannot open " << pvalues << "\n";
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
    for (int s = 0; s < sample; s++) {
        for (int l = 0; l < Ncat + 7; l++) { next_line(is, "skipped sample"); }
    }

    vector<double> dsdn = parse_doubles(next_line(is, "ds dn"));
    string treeline = next_line(is, "tree");
    vector<double> nucstat = parse_doubles(next_line(is, "nucstat"));
    vector<double> nucrr = parse_doubles(next_line(is, "nucrelrate"));
    vector<vector<double>> profiles;
    for (int i = 0; i < Ncat; i++) { profiles.push_back(parse_doubles(next_line(is, "profile"))); }
    vector<double> palloc_d = parse_doubles(next_line(is, "profile allocation"));
    vector<double> omegas = parse_doubles(next_line(is, "omegas"));
    vector<double> oalloc_d = parse_doubles(next_line(is, "omega allocation"));
    int Nsite = palloc_d.size();
    if (nucstat.size() != 4 || nucrr.size() != 6 || (int)omegas.size() != omegaNcat ||
        (int)oalloc_d.size() != Nsite) {
        cerr << "error: inconsistent .pvalues sample (nucstat " << nucstat.size() << ", nucrr "
             << nucrr.size() << ", omegas " << omegas.size() << "/" << omegaNcat
             << ", omega alloc " << oalloc_d.size() << "/" << Nsite << ")\n";
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
                &statespace, &nucmatrix, profiles[p], omegas[o], normalise));
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

    vector<vector<int>> state(Nnode, vector<int>(Nsite));
    double sum_nsub = 0;
    for (int rep = 0; rep < nrep; rep++) {
        long long nsub = 0;
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
                    cur = next;
                    nsub++;
                }
                state[node][site] = cur;
            }
        }
        double per_codon = (double)nsub / Nsite;
        sum_nsub += per_codon;
        cout << "rep " << rep << " Nsub " << nsub << " Nsub_per_codon " << per_codon << "\n";
        if (rep == nrep - 1) {
            ofstream os(prefix + ".ali");
            os << leaves.size() << "\t" << 3 * Nsite << "\n";
            for (int leaf : leaves) {
                os << tree->node_name(leaf) << "\t";
                for (int site = 0; site < Nsite; site++) {
                    os << statespace.GetState(state[leaf][site]);
                }
                os << "\n";
            }
            ofstream ot(prefix + ".tree");
            ot << treeline << "\n";
        }
    }
    double mean = sum_nsub / nrep;
    cout << setprecision(8) << "tree_length " << totlength << "\nmean_Nsub_per_codon " << mean
         << "\nratio_mean_Nsub_over_tree_length " << mean / totlength
         << "\nstationary_flux_per_unit_time " << flux
         << "\nexpected_Nsub_per_codon_(flux*length) " << flux * totlength << "\n";
    ofstream st(prefix + ".stats");
    st << setprecision(10) << "Nsite\t" << Nsite << "\nNtaxa\t" << leaves.size() << "\nnrep\t" << nrep
       << "\ntree_length\t" << totlength << "\nmean_Nsub_per_codon\t" << mean
       << "\nratio_over_tree_length\t" << mean / totlength << "\nstationary_flux\t" << flux
       << "\nexpected_Nsub_per_codon\t" << flux * totlength << "\n";
    return 0;
}
