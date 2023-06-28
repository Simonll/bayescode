#include <cmath>
#include <fstream>
#include "AACodonMutSelMultipleOmegaModel.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"
#include "components/ReadArgParse.hpp"
#include "components/stats_posterior.hpp"
#include "tclap/CmdLine.h"

using namespace std;
using namespace TCLAP;

class ReadAACodonMutSelDSBDPOmegaArgParse : public ReadArgParse {
  public:
    explicit ReadAACodonMutSelDSBDPOmegaArgParse(CmdLine& cmd) : ReadArgParse(cmd) {}

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
    SwitchArg sel{"d", "distribution", "Computes selection coefficients", cmd};
    ValueArg<double> omega_pp{"", "omega_threshold",
        "Threshold to compute the mean posterior probability that ω⁎ "
        "(or ω if option `flatfitness` is used in `mutselomega`) is greater than a given value.",
        false, 1.0, "double", cmd};
};

int main(int argc, char* argv[]) {
    CmdLine cmd{"AACodonMutSelMultipleOmega", ' ', "0.1"};
    ReadAACodonMutSelDSBDPOmegaArgParse read_args(cmd);
    cmd.parse(argc, argv);

    string chain_name = read_args.GetChainName();
    int burnin = read_args.GetBurnIn();
    int every = read_args.GetEvery();
    int size = read_args.GetSize();

    ifstream is{chain_name + ".param"};
    ChainDriver::fake_read(is);  // We're not interested in the ChainDriver of the param file
    AACodonMutSelMultipleOmegaModel model(is);
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
    } else if (read_args.sel.getValue()) {
        int Ncat = 241;
        double min = -30;
        double max = 30;
        double bin = 0.25;
        // vector<vector<double>> sshistoMut(model.GetNsite(), {0});
        // vector<vector<double>> sshistoSub(model.GetNsite(), {0});
        // vector<vector<double>> sshistoNonsynMut(model.GetNsite(), {0});
        // vector<vector<double>> sshistoNonsynSub(model.GetNsite(), {0});
        // vector<vector<double>> sshistoSynMut(model.GetNsite(), {0});
        // vector<vector<double>> sshistoSynSub(model.GetNsite(), {0});


        // vector<double> ssStatNonsynSubRate(model.GetNsite(), {0});
        // vector<double> ssStatSynSubRate(model.GetNsite(), {0});
        // vector<double> ssStatNonsynMutRate(model.GetNsite(), {0});
        // vector<double> ssStatSynMutRate(model.GetNsite(), {0});

        vector<double> ghistoMut(Ncat, {0});
        vector<double> ghistoSub(Ncat, {0});
        vector<double> ghistoNonsynMut(Ncat, {0});
        vector<double> ghistoNonsynSub(Ncat, {0});
        vector<double> ghistoSynMut(Ncat, {0});
        vector<double> ghistoSynSub(Ncat, {0});
        vector<double> shistoMut(Ncat, {0});
        vector<double> shistoSub(Ncat, {0});
        vector<double> shistoNonsynMut(Ncat, {0});
        vector<double> shistoNonsynSub(Ncat, {0});
        vector<double> shistoSynMut(Ncat, {0});
        vector<double> shistoSynSub(Ncat, {0});
        // vector<double> tsshistoMut(Ncat, {0});
        // vector<double> tsshistoSub(Ncat, {0});
        // vector<double> tsshistoNonsynMut(Ncat, {0});
        // vector<double> tsshistoNonsynSub(Ncat, {0});
        // vector<double> tsshistoSynMut(Ncat, {0});
        // vector<double> tsshistoSynSub(Ncat, {0});

        vector<double> stat(model.GetCodonStateSpace()->GetNstate(), {0});
        int pos, nucFrom, nucTo, nucRRIndex, c;
        double statMutRate, deltaS, statSubRate, totalMut, totalSub, totalNonsynMut, totalNonsynSub,
            totalSynMut, totalSynSub, siteTotalMut, siteTotalSub, siteTotalNonsynMut,
            siteTotalNonsynSub, siteTotalSynMut, siteTotalSynSub, Z;
        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            totalMut = 0;
            totalSub = 0;
            totalNonsynMut = 0;
            totalNonsynSub = 0;
            totalSynMut = 0;
            totalSynSub = 0;

            shistoMut.resize(Ncat, {0});
            shistoSub.resize(Ncat, {0});
            shistoNonsynMut.resize(Ncat, {0});
            shistoNonsynSub.resize(Ncat, {0});
            shistoSynMut.resize(Ncat, {0});
            shistoSynSub.resize(Ncat, {0});

            for (int site = 0; site < model.GetNsite(); site++) {
                // siteTotalMut = 0;
                // siteTotalSub = 0;
                // siteTotalNonsynMut = 0;
                // siteTotalNonsynSub = 0;
                // siteTotalSynMut = 0;
                // siteTotalSynSub = 0;
                // tsshistoMut.resize(Ncat, {0});
                // tsshistoSub.resize(Ncat, {0});
                // tsshistoNonsynMut.resize(Ncat, {0});
                // tsshistoNonsynSub.resize(Ncat, {0});
                // tsshistoSynMut.resize(Ncat, {0});
                // tsshistoSynSub.resize(Ncat, {0});

                Z = 0;
                for (int s = 0; s < model.GetCodonStateSpace()->GetNstate(); s++) {
                    stat[s] =
                        model.GetNucStat(model.GetCodonStateSpace()->GetCodonPosition(0, s)) *
                        model.GetNucStat(model.GetCodonStateSpace()->GetCodonPosition(1, s)) *
                        model.GetNucStat(model.GetCodonStateSpace()->GetCodonPosition(2, s)) *
                        model.GetCodonFitness(s) *
                        model.GetAASiteFitness(site, model.GetCodonStateSpace()->Translation(s));
                    Z += stat[s];
                }
                for (int s = 0; s < model.GetCodonStateSpace()->GetNstate(); s++) { stat[s] /= Z; }
                int nonsyncount = 0;
                int syncount = 0;
                for (int codonFrom = 0; codonFrom < model.GetCodonStateSpace()->GetNstate();
                     codonFrom++) {
                    for (int codonTo = 0; codonTo < model.GetCodonStateSpace()->GetNstate();
                         codonTo++) {
                        pos = model.GetCodonStateSpace()->GetDifferingPosition(codonFrom, codonTo);
                        if ((pos != -1) && (pos != 3)) {
                            nucFrom = model.GetCodonStateSpace()->GetCodonPosition(pos, codonFrom);
                            nucTo = model.GetCodonStateSpace()->GetCodonPosition(pos, codonTo);
                            if (nucFrom < nucTo) {
                                nucRRIndex =
                                    (2 * Nnuc - nucFrom - 1) * nucFrom / 2 + nucTo - nucFrom - 1;
                            } else {
                                nucRRIndex =
                                    (2 * Nnuc - nucTo - 1) * nucTo / 2 + nucFrom - nucTo - 1;
                            }
                            statMutRate = model.GetNucRR(nucRRIndex) * model.GetNucStat(nucTo) *
                                          stat[codonFrom];
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
                                    tmp = min + ((double)(c)*bin) - bin / 2 + bin;
                                } while (tmp < deltaS);
                            }
                            if (c == Ncat) {
                                cout << "error, c==Ncat.\n";
                                cout.flush();
                            }

                            shistoMut[c] += statMutRate;
                            shistoSub[c] += statSubRate;
                            // tsshistoMut[c] += statMutRate;
                            // tsshistoSub[c] += statSubRate;
                            if (!model.GetCodonStateSpace()->Synonymous(codonFrom, codonTo)) {
                                shistoNonsynMut[c] += statMutRate;
                                shistoNonsynSub[c] += statSubRate;
                                totalNonsynMut += statMutRate;
                                totalNonsynSub += statSubRate;

                                // tsshistoNonsynMut[c] += statMutRate;
                                // tsshistoNonsynSub[c] += statSubRate;
                                // siteTotalNonsynMut += statMutRate;
                                // siteTotalNonsynSub += statSubRate;
                            } else {
                                shistoSynMut[c] += statMutRate;
                                shistoSynSub[c] += statSubRate;
                                totalSynMut += statMutRate;
                                totalSynSub += statSubRate;

                                // tsshistoSynMut[c] += statMutRate;
                                // tsshistoSynSub[c] += statSubRate;
                                // siteTotalSynMut += statMutRate;
                                // siteTotalSynSub += statSubRate;
                            }
                            totalMut += statMutRate;
                            totalSub += statSubRate;
                            // siteTotalMut += statMutRate;
                            // siteTotalSub += statSubRate;
                        }
                    }
                }

                // for (c = 0; c < Ncat; c++) {
                //     sshistoMut[site][c] += tsshistoMut[c] / siteTotalMut;
                //     sshistoSub[site][c] += tsshistoSub[c] / siteTotalSub;
                //     sshistoNonsynMut[site][c] += tsshistoNonsynMut[c] / siteTotalNonsynMut;
                //     sshistoNonsynSub[site][c] += tsshistoNonsynSub[c] / siteTotalNonsynSub;
                //     sshistoSynMut[site][c] += tsshistoSynMut[c] / siteTotalSynMut;
                //     sshistoSynSub[site][c] += tsshistoSynSub[c] / siteTotalSynSub;
                // }

                // ssStatNonsynSubRate[site] += siteTotalNonsynSub;
                // ssStatNonsynMutRate[site] += siteTotalNonsynMut;
                // ssStatSynSubRate[site] += siteTotalSynSub;
                // ssStatSynMutRate[site] += siteTotalSynMut;
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

        for (c = 0; c < Ncat; c++) {
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