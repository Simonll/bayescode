#include <algorithm>
#include <cmath>
#include <fstream>
#include "CodonMutSelMultipleOmegaModel.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"
#include "components/ReadArgParse.hpp"
#include "components/stats_posterior.hpp"
#include "tclap/CmdLine.h"
#include "tree/export.hpp"

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
    SwitchArg selection_coef_distribution{"d", "distribution", "Computes selection coefficients", cmd};
    SwitchArg prepare_files_for_simulation{"", "for_simulation", "Prepare files for jump chain simulations", cmd};
    SwitchArg site_specific{"s", "ss",
        "Computes the mean posterior site-specific codon equilibrium frequencies"
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
    } else if (read_args.site_specific.getValue()) {
        vector<vector<double>> sitestat(model.GetNsite(), vector<double>(Nstate, 0.0));      
        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int i = 0; i < model.GetNsite(); i++) {
                vector<double> const &profile = model.GetProfile(i);
                assert(profile.size() == Nstate);               
                for (unsigned k{0}; k < profile.size(); k++) { sitestat[i][k] += profile[k]; }
            }
        }
        cerr << '\n';

        ofstream os(chain_name + ".siteprofiles");
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
        ofstream os(filename);
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
        ofstream os(filename);
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
    } else if (read_args.prepare_files_for_simulation.getValue()) {
        string filename{chain_name + ".pvalues"};
        std::ofstream os(filename);
        model.GetModelStamp(os);
        os << '\n';
        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            os << model.GetPredictedRelativedS() << "\t" << model.GetPredictedRelativedN() << "\n";
            ExportTree export_tree(model.GetTree());
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(model.GetTree().nb_nodes());
                 node++) {
                if (!model.GetTree().is_root(node)) {
                    export_tree.set_tag(node, "length", to_string(3 * model.GetBranchLength(node)));
                }
            }
            os << export_tree.as_string() << '\n';
            for (int i = 0; i < Nnuc; i++) {
                os << model.GetNucStat(i) << (i == Nnuc - 1 ? '\n' : '\t');
            }
            for (int i = 0; i < Nrr; i++) {
                os << model.GetNucRR(i) << (i == Nrr - 1 ? '\n' : '\t');
            }
            for (int i = 0; i < model.GetNcat(); i++) {
                for (int j = 0; j < Nstate; j++) {
                    os << model.GetProfileCodon(i, j) << (j == Nstate - 1 ? '\n' : '\t');
                }
            }
            for (int i = 0; i < model.GetNsite(); i++) {
                os << model.GetProfileAlloc(i) << (i == model.GetNsite() - 1 ? '\n' : '\t');
                
            }
            for (int i = 0; i < model.GetOmegaNcat(); i++) {
                os << model.GetOmega(i) << (i == model.GetOmegaNcat() - 1 ? '\n' : '\t');
            }
            for (int i = 0; i < model.GetNsite(); i++) {
                os << model.GetOmegaAlloc(i) << (i == model.GetNsite() - 1 ? '\n' : '\t');
            }
        }
        os.close();
    } else if (read_args.selection_coef_distribution.getValue()) {
        int Ncat = 241;
        double min = -30;
        double max = 30;
        double bin = 0.25;
        vector<double> ghistoMut(Ncat, 0.0);
        vector<double> ghistoSub(Ncat, 0.0);
        vector<double> ghistoNonsynMut(Ncat, 0.0);
        vector<double> ghistoNonsynSub(Ncat, 0.0);
        vector<double> ghistoSynMut(Ncat, 0.0);
        vector<double> ghistoSynSub(Ncat, 0.0);
        vector<double> stat(Nstate, 0.0);
        vector<double> shistoMut(Ncat, 0.0);
        vector<double> shistoSub(Ncat, 0.0);
        vector<double> shistoNonsynMut(Ncat, 0.0);
        vector<double> shistoNonsynSub(Ncat, 0.0);
        vector<double> shistoSynMut(Ncat, 0.0);
        vector<double> shistoSynSub(Ncat, 0.0);

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            
            double statMutRate, deltaS, statSubRate;
            double totalMut = 0;
            double totalSub = 0;
            double totalNonsynMut = 0;
            double totalNonsynSub = 0;
            double totalSynMut = 0;
            double totalSynSub = 0;
            std::fill(stat.begin(), stat.end(), 0.0);
            std::fill(shistoMut.begin(), shistoMut.end(), 0.0);
            std::fill(shistoSub.begin(), shistoSub.end(), 0.0);
            std::fill(shistoNonsynMut.begin(), shistoNonsynMut.end(), 0.0);
            std::fill(shistoNonsynSub.begin(), shistoNonsynSub.end(), 0.0);
            std::fill(shistoSynMut.begin(), shistoSynMut.end(), 0.0);
            std::fill(shistoSynSub.begin(), shistoSynSub.end(), 0.0);

            for (int site = 0; site < model.GetNsite(); site++) {
                double Z = 0;
                for (int state = 0; state < Nstate; state++) {
                    stat[state] = model.GetNucStat(model.GetCodonStateSpace()->GetCodonPosition(0, state)) *
                                  model.GetNucStat(model.GetCodonStateSpace()->GetCodonPosition(1, state)) *
                                  model.GetNucStat(model.GetCodonStateSpace()->GetCodonPosition(2, state)) *
                                  model.GetSiteCodonFitness(site, state);
                    Z += stat[state];
                }
                for (int state = 0; state < Nstate; state++) { stat[state] /= Z; }

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
                        
                        deltaS = log(model.GetSiteCodonFitness(site, codonTo)) -
                                 log(model.GetSiteCodonFitness(site, codonFrom));

                        if ((fabs(deltaS)) < 1e-30) {
                            statSubRate = statMutRate / (1.0 - deltaS / 2.0);
                        } else {
                            statSubRate = statMutRate * (deltaS / (1.0 - exp(-deltaS)));
                        }

                        if (!model.GetCodonStateSpace()->Synonymous(codonFrom, codonTo)) {
                            statSubRate *= model.GetSiteOmega(site);
                        }
                        int c = std::max(0, std::min(Ncat - 1, static_cast<int>(std::floor((deltaS - min) / bin))));
                        // int c = std::clamp(static_cast<int>(std::floor((deltaS - min) / bin)), 0, Ncat - 1);

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
            for (int c = 0; c < Ncat; c++) {
                ghistoMut[c] += shistoMut[c] / totalMut;
                ghistoSub[c] += shistoSub[c] / totalSub;
                ghistoNonsynMut[c] += shistoNonsynMut[c] / totalNonsynMut;
                ghistoNonsynSub[c] += shistoNonsynSub[c] / totalNonsynSub;
                ghistoSynMut[c] += shistoSynMut[c] / totalSynMut;
                ghistoSynSub[c] += shistoSynSub[c] / totalSynSub;
            }
        }
        ofstream mutmutsel_os(chain_name + ".mutsel", std::ios::out);
        ofstream mutsubsel_os(chain_name + ".subsel", std::ios::out);
        ofstream nonsynmutmutsel_os(chain_name + ".nonsynmutsel", std::ios::out);
        ofstream nonsynmutsubsel_os(chain_name + ".nonsynsubsel", std::ios::out);
        ofstream synmutmutsel_os(chain_name + ".synmutsel", std::ios::out);
        ofstream synmutsubsel_os(chain_name + ".synsubsel", std::ios::out);

        for (int c = 0; c < Ncat; c++) {
            mutmutsel_os << (min + (c * bin)) << "\t" << (ghistoMut[c] / size) << '\n';
            mutsubsel_os << (min + (c * bin)) << "\t" << (ghistoSub[c] / size) << '\n';
            nonsynmutmutsel_os << (min + (c * bin)) << "\t" << (ghistoNonsynMut[c] / size) << '\n';
            nonsynmutsubsel_os << (min + (c * bin)) << "\t" << (ghistoNonsynSub[c] / size) << '\n';
            synmutmutsel_os << (min + (c * bin)) << "\t" << (ghistoSynMut[c] / size) << '\n';
            synmutsubsel_os << (min + (c * bin)) << "\t" << (ghistoSynSub[c] / size) << '\n';
        }
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
        ofstream os(filename);
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