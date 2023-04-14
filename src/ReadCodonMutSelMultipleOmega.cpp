#include <cmath>
#include <fstream>
#include "CodonMutSelMultipleOmegaModel.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"
#include "components/ReadArgParse.hpp"
#include "components/stats_posterior.hpp"
#include "tclap/CmdLine.h"

using namespace std;
using namespace TCLAP;

class ReadCodonMutSelMultipleOmegaArgParse : public ReadArgParse {
  public:
    explicit ReadCodonMutSelMultipleOmegaArgParse(CmdLine &cmd) : ReadArgParse(cmd) {}

    SwitchArg nuc{"n", "nuc", "Mean posterior nucleotide matrix.", cmd};
    ValueArg<string> confidence_interval{"c", "confidence_interval",
        "Posterior credible interval for ω (per site and at the gene level).", false, "", "string",
        cmd};
    SwitchArg omega_knot{"", "omega_0",
        "Posterior credible interval for ω0 predicted at the mutation-selection "
        "equilibrium from the fitness profiles (instead of ω). "
        "To use combined with the option `confidence_interval`.",
        cmd};
    SwitchArg ss{"s", "ss",
        "Computes the mean posterior site-specific amino-acid equilibrium frequencies"
        "(amino-acid fitness profiles).",
        cmd};
    ValueArg<double> omega_pp{"", "omega_threshold",
        "Threshold to compute the mean posterior probability that ω⁎ "
        "(or ω if option `flatfitness` is used in `mutselomega`) is greater than a given value.",
        false, 1.0, "double", cmd};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"CodonMutSelMultipleOmega", ' ', "0.1"};
    ReadCodonMutSelMultipleOmegaArgParse read_args(cmd);
    cmd.parse(argc, argv);

    string chain_name = read_args.GetChainName();
    int burnin = read_args.GetBurnIn();
    int every = read_args.GetEvery();
    int size = read_args.GetSize();

    ifstream is{chain_name + ".param"};
    ChainDriver::fake_read(is);  // We're not interested in the ChainDriver of the param file
    CodonMutSelMultipleOmegaModel model(is);
    ChainReader cr{model, chain_name + ".chain"};

    cr.skip(burnin);
    cerr << size << " points to read\n";

    if (read_args.GetPpred()) {
        for (int i = 0; i < size; i++) {
            cerr << '.';
            cr.skip(every);
            model.PostPred("ppred_" + chain_name + "_" + to_string(i) + ".ali");
        }
        cerr << '\n';
    } else if (read_args.ss.getValue()) {
        vector<vector<double>> sitestat(model.GetNsite(), {0});

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int i = 0; i < model.GetNsite(); i++) {
                vector<double> const &profile = model.GetProfile(i);
                if (sitestat[i].size() != profile.size()) {
                    sitestat[i].resize(profile.size(), 0);
                };
                for (unsigned k{0}; k < profile.size(); k++) { sitestat[i][k] += profile[k]; }
            }
        }
        cerr << '\n';

        ofstream os((chain_name + ".siteprofiles").c_str());

        os << "site\tTTT\tTTC\tTTA\tTTG\tTCT\tTCC\tTCA\tTCG\tTAT\tTAC\tTGT\tTGC\tTGG\tCTT\tCTC\tCTA"
              "\tCTG\tCCT\tCCC\tCCA\tCCG\tCAT\tCAC\tCAA\tCAG\tCGT\tCGC\tCGA\tCGG\tATT\tATC\tATA\tAT"
              "G\tACT\tACC\tACA\tACG\tAAT\tAAC\tAAA\tAAG\tAGT\tAGC\tAGA\tAGG\tGTT\tGTC\tGTA\tGTG\tG"
              "CT\tGCC\tGCA\tGCG\tGAT\tGAC\tGAA\tGAG\tGGT\tGGC\tGGA\tGGG\n";
        for (int i = 0; i < model.GetNsite(); i++) {
            os << i + 1;
            for (auto &codon : sitestat[i]) {
                codon /= size;
                os << '\t' << codon;
            }
            os << '\n';
        }
        cerr << "mean site-specific profiles in " << chain_name << ".siteprofiles\n";
        cerr << '\n';
    } else if (!read_args.confidence_interval.getValue().empty()) {
        double ci = stod(read_args.confidence_interval.getValue());
        vector<vector<double>> omega(model.GetNsite());
        vector<double> gene_omega{};
        double upper = max(ci, 1.0 - ci);
        double lower = min(ci, 1.0 - ci);

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            double mean{0.0};
            for (int site = 0; site < model.GetNsite(); site++) {
                double val = read_args.omega_knot.getValue() ? model.GetPredictedSiteOmegaKnot(site)
                                                             : model.GetSiteOmega(site);
                omega[site].push_back(val);
                mean += val;
            }
            gene_omega.push_back(mean / model.GetNsite());
        }
        cerr << '\n';

        string filename{chain_name + ".ci" + read_args.confidence_interval.getValue() + ".tsv"};
        ofstream os(filename.c_str());
        os << "#site\tomega_lower\tomega\tomega_upper\n";

        double mean = accumulate(gene_omega.begin(), gene_omega.end(), 0.0) / size;
        sort(gene_omega.begin(), gene_omega.end());
        auto pt_up = static_cast<size_t>(upper * gene_omega.size());
        if (pt_up >= gene_omega.size()) { pt_up = gene_omega.size() - 1; }
        double up = gene_omega.at(pt_up);
        double down = gene_omega.at(static_cast<size_t>(lower * gene_omega.size()));
        os << "#Mean\t" << down << '\t' << mean << '\t' << up << '\n';

        for (int i = 0; i < model.GetNsite(); i++) {
            mean = accumulate(omega[i].begin(), omega[i].end(), 0.0) / size;
            sort(omega[i].begin(), omega[i].end());
            pt_up = static_cast<size_t>(upper * omega[i].size());
            if (pt_up >= omega[i].size()) { pt_up = omega[i].size() - 1; }
            up = omega[i].at(pt_up);
            down = omega[i].at(static_cast<size_t>(lower * omega[i].size()));
            os << i + 1 << '\t' << down << '\t' << mean << '\t' << up << '\n';
        }
        cerr << '\n';
    } else if (read_args.nuc.getValue()) {
        vector<vector<double>> rates(Nnuc * Nnuc);
        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int i = 0; i < Nnuc; i++) {
                for (int j = 0; j < Nnuc; j++) {
                    if (i != j) {
                        int r = i * Nnuc + j;
                        rates[r].push_back(model.GetNucRate(i, j));
                    }
                }
            }
        }
        cerr << '\n';
        string filename{chain_name + ".nucmatrix.tsv"};
        ofstream os(filename.c_str());
        os << "Name\tRate\n";
        for (int i = 0; i < Nnuc; i++) {
            for (int j = 0; j < Nnuc; j++) {
                if (i != j) {
                    int r = i * Nnuc + j;
                    double q_mean = accumulate(rates.at(r).begin(), rates.at(r).end(), 0.0) / size;
                    os << "q_" << DNAletters[i] << "_" << DNAletters[j] << "\t" << q_mean << '\n';
                }
            }
        }
        cerr << '\n';
    } else {
        vector<double> omegappgto(model.GetNsite(), 0);
        vector<double> omega(model.GetNsite(), 0);

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int site = 0; site < model.GetNsite(); site++) {
                omega[site] += model.GetSiteOmega(site);
                if (model.GetSiteOmega(site) > read_args.omega_pp.getValue()) {
                    omegappgto[site]++;
                }
            }
        }
        cerr << '\n';

        string filename{chain_name + ".omegappgt" + to_string(read_args.omega_pp.getValue())};
        ofstream os(filename.c_str());
        if (model.FlatFitness()) {
            os << "#site\tp(ω>" << read_args.omega_pp.getValue() << ")\tω\n";
        } else {
            os << "#site\tp(ω*>" << read_args.omega_pp.getValue() << ")\tω*\n";
        }

        for (int i = 0; i < model.GetNsite(); i++) {
            os << i + 1 << '\t' << omegappgto[i] / size << '\t' << omega[i] / size << '\n';
        }
        cerr << "Posterior prob of omega greater than " << read_args.omega_pp.getValue() << " in "
             << filename << "\n";
        cerr << '\n';
    }
}