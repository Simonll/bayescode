#include <cmath>
#include <fstream>
#include "AACodonMutSelMultipleOmegaModel.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"
#include "components/ReadArgParse.hpp"
#include "components/stats_posterior.hpp"
#include "tclap/CmdLine.h"
#include "tree/export.hpp"

using namespace std;
using namespace TCLAP;

class ReadAACodonMutSelDSBDPOmegaArgParse : public ReadArgParse {
  public:
    explicit ReadAACodonMutSelDSBDPOmegaArgParse(CmdLine& cmd) : ReadArgParse(cmd) {}

    SwitchArg nuc{"n", "nuc", "Mean posterior nucleotide matrix.", cmd};
    SwitchArg codonfitness{"f", "codon",
        "Computes the mean posterior codon equilibrium frequencies"
        "(codon fitness profile).",
        cmd};
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
    SwitchArg sel{"d", "distribution", "Computes scaled selection coefficients", cmd};
    SwitchArg simu{"", "for_simulation", "Prepare files for jump chain simulations", cmd};
    ValueArg<double> omega_pp{"", "omega_threshold",
        "Threshold to compute the mean posterior probability that ω⁎ "
        "(or ω if option `flatfitness` is used in `mutselomega`) is greater than a given value.",
        false, 1.0, "double", cmd};
};

int main(int argc, char* argv[]) {
    CmdLine cmd{AACodonMutSelMultipleOmega::GetModelName(), ' ', "0.1"};
    ReadAACodonMutSelDSBDPOmegaArgParse read_args(cmd);
    cmd.parse(argc, argv);

    string chain_name = read_args.GetChainName();
    int burnin = read_args.GetBurnIn();
    int every = read_args.GetEvery();
    int size = read_args.GetSize();


    ifstream is{chain_name + ".param"};
    ChainDriver::fake_read(is);  // We're not interested in the ChainDriver of the param file
    AACodonMutSelMultipleOmega model(is);
    ChainReader cr{model, chain_name + ".chain"};
    int Nstate = model.GetCodonStateSpace()->GetNstate();
    cr.skip(burnin);
    cerr << size << " points to read\n";

    if (read_args.GetPpred()) {
        for (int i = 0; i < size; i++) {
            cerr << '.';
            cr.skip(every);
            model.PostPred(chain_name + "_" + "ppred_" + to_string(i) + ".ali");
        }
        cerr << '\n';
    } else if (read_args.ss.getValue()) {
        vector<vector<double>> sitestat(model.GetNsite(), {0});

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int i = 0; i < model.GetNsite(); i++) {
                vector<double> const& profile = model.GetProfile(i);
                if (sitestat[i].size() != profile.size()) {
                    sitestat[i].resize(profile.size(), 0);
                };
                for (unsigned k{0}; k < profile.size(); k++) { sitestat[i][k] += profile[k]; }
            }
        }
        cerr << '\n';

        ofstream os((chain_name + ".siteprofiles").c_str());
        os << "site\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\n";
        for (int i = 0; i < model.GetNsite(); i++) {
            os << i + 1;
            for (auto& aa : sitestat[i]) {
                aa /= size;
                os << '\t' << aa;
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
    } else if (read_args.codonfitness.getValue()) {
        vector<double> codonfitness(Nstate, {0});
        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int i = 0; i < Nstate; i++) { codonfitness[i] += model.GetCodonFitness(i); }
        }
        cerr << '\n';
        string filename{chain_name + ".codonfitness.tsv"};
        ofstream os(filename.c_str());
        os << "TTT\tTTC\tTTA\tTTG\tTCT\tTCC\tTCA\tTCG\tTAT\tTAC\tTGT\tTGC\tTGG\tCTT\tCTC\tCTA"
              "\tCTG\tCCT\tCCC\tCCA\tCCG\tCAT\tCAC\tCAA\tCAG\tCGT\tCGC\tCGA\tCGG\tATT\tATC\tATA\tAT"
              "G\tACT\tACC\tACA\tACG\tAAT\tAAC\tAAA\tAAG\tAGT\tAGC\tAGA\tAGG\tGTT\tGTC\tGTA\tGTG\tG"
              "CT\tGCC\tGCA\tGCG\tGAT\tGAC\tGAA\tGAG\tGGT\tGGC\tGGA\tGGG\n";
        for (int i = 0; i < Nstate; i++) { os << codonfitness[i] / size; }
        cerr << '\n';
    } else if (read_args.simu.getValue()) {
        string filename{chain_name + ".pvalues"};
        std::ofstream os(filename.c_str());
        model.GetModelStamp(os);
        os << '\n';
        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            ExportTree export_tree(model.GetTree());
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(model.GetTree().nb_nodes());
                 node++) {
                if (!model.GetTree().is_root(node)) {
                    export_tree.set_tag(node, "length", to_string(model.GetBranchLength(node)));
                }
            }
            os << export_tree.as_string() << '\n';
            for (int i = 0; i < Nnuc; i++) {
                if (i != Nnuc - 1) {
                    os << model.GetNucStat(i) << '\t';
                } else {
                    os << model.GetNucStat(i) << '\n';
                }
            }
            for (int i = 0; i < Nrr; i++) {
                if (i != Nrr - 1) {
                    os << model.GetNucRR(i) << '\t';
                } else {
                    os << model.GetNucRR(i) << '\n';
                }
            }
            for (int i = 0; i < Nstate; i++) {
                if (i != Nstate - 1) {
                    os << model.GetCodonFitness(i) << '\t';
                } else {
                    os << model.GetCodonFitness(i) << '\n';
                }
            }
            for (int i = 0; i < model.GetNsite(); i++) {
                for (int j = 0; j < Naa; j++) {
                    if (j != Naa - 1) {
                        os << model.GetAASiteFitness(i, j) << '\t';
                    } else {
                        os << model.GetAASiteFitness(i, j) << '\n';
                    }
                }
            }
            for (int i = 0; i < model.GetNsite(); i++) {
                if (i != model.GetNsite() - 1) {
                    os << model.GetSiteOmega(i) << '\t';
                } else {
                    os << model.GetSiteOmega(i) << '\n';
                }
            }
        }
        os.close();
    }

    else if (read_args.sel.getValue()) {
        int Ncat = 241;
        double min = -30;
        double max = 30;
        double bin = 0.25;

        vector<double> ghistoMut(Ncat, {0});
        vector<double> ghistoSub(Ncat, {0});
        vector<double> ghistoNonsynMut(Ncat, {0});
        vector<double> ghistoNonsynSub(Ncat, {0});
        vector<double> ghistoSynMut(Ncat, {0});
        vector<double> ghistoSynSub(Ncat, {0});


        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            int c = 0;
            double statMutRate, deltaS, statSubRate;
            double totalMut = 0;
            double totalSub = 0;
            double totalNonsynMut = 0;
            double totalNonsynSub = 0;
            double totalSynMut = 0;
            double totalSynSub = 0;
            vector<double> stat(Nstate, {0});
            vector<double> shistoMut(Ncat, {0});
            vector<double> shistoSub(Ncat, {0});
            vector<double> shistoNonsynMut(Ncat, {0});
            vector<double> shistoNonsynSub(Ncat, {0});
            vector<double> shistoSynMut(Ncat, {0});
            vector<double> shistoSynSub(Ncat, {0});

            for (int site = 0; site < model.GetNsite(); site++) {
                double Z = 0;
                for (int s = 0; s < Nstate; s++) {
                    stat[s] =
                        model.GetNucStat(model.GetCodonStateSpace()->GetCodonPosition(0, s)) *
                        model.GetNucStat(model.GetCodonStateSpace()->GetCodonPosition(1, s)) *
                        model.GetNucStat(model.GetCodonStateSpace()->GetCodonPosition(2, s)) *
                        model.GetCodonFitness(s) *
                        model.GetAASiteFitness(site, model.GetCodonStateSpace()->Translation(s));
                    Z += stat[s];
                }
                for (int s = 0; s < Nstate; s++) { stat[s] /= Z; }
                int nonsyncount = 0;
                int syncount = 0;
                for (int codonFrom = 0; codonFrom < Nstate; codonFrom++) {
                    for (auto codonTo : model.GetCodonStateSpace()->GetNeighbors(codonFrom)) {
                        double pos =
                            model.GetCodonStateSpace()->GetDifferingPosition(codonFrom, codonTo);
                        double nucFrom =
                            model.GetCodonStateSpace()->GetCodonPosition(pos, codonFrom);
                        double nucTo = model.GetCodonStateSpace()->GetCodonPosition(pos, codonTo);
                        double nucRRIndex;
                        if (nucFrom < nucTo) {
                            nucRRIndex =
                                (2 * Nnuc - nucFrom - 1) * nucFrom / 2 + nucTo - nucFrom - 1;
                        } else {
                            nucRRIndex = (2 * Nnuc - nucTo - 1) * nucTo / 2 + nucFrom - nucTo - 1;
                        }
                        statMutRate =
                            model.GetNucRR(nucRRIndex) * model.GetNucStat(nucTo) * stat[codonFrom];

                        if (!model.GetCodonStateSpace()->Synonymous(codonFrom, codonTo)) {
                            int aaFrom = model.GetCodonStateSpace()->Translation(codonFrom);
                            int aaTo = model.GetCodonStateSpace()->Translation(codonTo);
                            deltaS = log(model.GetAASiteFitness(site, aaTo)) -
                                     log(model.GetAASiteFitness(site, aaFrom)) +
                                     log(model.GetCodonFitness(codonTo)) -
                                     log(model.GetCodonFitness(codonFrom));
                        } else {
                            deltaS = log(model.GetCodonFitness(codonTo)) -
                                     log(model.GetCodonFitness(codonFrom));
                        }

                        if ((fabs(deltaS)) < 1e-30) {
                            statSubRate = statMutRate * 1.0 / (1.0 - (deltaS / 2));
                        } else {
                            statSubRate = statMutRate * (deltaS / (1.0 - exp(-deltaS)));
                        }

                        if (!model.GetCodonStateSpace()->Synonymous(codonFrom, codonTo)) {
                            statSubRate *= model.GetSiteOmega(site);
                            nonsyncount++;
                        } else {
                            syncount++;
                        }

                        if (deltaS < min) {
                            c = 0;
                        } else if (deltaS > max) {
                            c = Ncat - 1;
                        } else {
                            c = 0;
                            double tmp = min + ((double)c * bin) - bin / 2 + bin;
                            do {
                                c++;
                                tmp = min + ((double)c * bin) - bin / 2 + bin;
                            } while (tmp < deltaS);
                        }
                        if (c == Ncat) {
                            cout << "error, c==Ncat.\n";
                            cout.flush();
                        }

                        if (!model.GetCodonStateSpace()->Synonymous(codonFrom, codonTo)) {
                            shistoNonsynMut[c] += statMutRate;
                            shistoNonsynSub[c] += statSubRate;
                            totalNonsynMut += statMutRate;
                            totalNonsynSub += statSubRate;
                        } else {
                            shistoSynMut[c] += statMutRate;
                            shistoSynSub[c] += statSubRate;
                            totalSynMut += statMutRate;
                            totalSynSub += statSubRate;
                        }
                        shistoMut[c] += statMutRate;
                        shistoSub[c] += statSubRate;
                        totalMut += statMutRate;
                        totalSub += statSubRate;
                    }
                }
            }

            for (c = 0; c < Ncat; c++) {
                ghistoMut[c] += shistoMut[c] / totalMut;
                ghistoSub[c] += shistoSub[c] / totalSub;
                ghistoNonsynMut[c] += shistoNonsynMut[c] / totalNonsynMut;
                ghistoNonsynSub[c] += shistoNonsynSub[c] / totalNonsynSub;
                ghistoSynMut[c] += shistoSynMut[c] / totalSynMut;
                ghistoSynSub[c] += shistoSynSub[c] / totalSynSub;
            }
        }
        ofstream mutmutsel_os((chain_name + ".mutsel").c_str(), std::ios::out);
        ofstream mutsubsel_os((chain_name + ".subsel").c_str(), std::ios::out);
        ofstream nonsynmutmutsel_os((chain_name + ".nonsynmutsel").c_str(), std::ios::out);
        ofstream nonsynmutsubsel_os((chain_name + ".nonsynsubsel").c_str(), std::ios::out);
        ofstream synmutmutsel_os((chain_name + ".synmutsel").c_str(), std::ios::out);
        ofstream synmutsubsel_os((chain_name + ".synsubsel").c_str(), std::ios::out);

        for (int c = 0; c < Ncat; c++) {
            mutmutsel_os << (min + (c * bin)) << "\t" << (ghistoMut[c] / size) << '\n';
            mutsubsel_os << (min + (c * bin)) << "\t" << (ghistoSub[c] / size) << '\n';
            nonsynmutmutsel_os << (min + (c * bin)) << "\t" << (ghistoNonsynMut[c] / size) << '\n';
            nonsynmutsubsel_os << (min + (c * bin)) << "\t" << (ghistoNonsynSub[c] / size) << '\n';
            synmutmutsel_os << (min + (c * bin)) << "\t" << (ghistoSynMut[c] / size) << '\n';
            synmutsubsel_os << (min + (c * bin)) << "\t" << (ghistoSynSub[c] / size) << '\n';
        }
    }

    else {
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