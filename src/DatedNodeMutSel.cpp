#include <cmath>
#include <fstream>
#include "DatedNodeMutSelModel.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/InferenceAppArgParse.hpp"
#include "components/StandardTracer.hpp"
#include "components/restart_check.hpp"

using namespace std;

class DatedNodeMutselArgParse : public BaseArgParse {
  public:
    DatedNodeMutselArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}

    ValueArg<int> ncat{
        "", "ncat", "truncation of the first-level stick-breaking process", false, 100, "int", cmd};
    ValueArg<int> basencat{"", "basencat", "truncation of the second-level stick-breaking process",
        false, 1, "int", cmd};
    SwitchArg condition_aware{"b", "condition_aware",
        "One Ne per condition, if the tree doesn't have condition, then one Ne per branch", cmd,
        false};
    ValueArg<std::string> traitsfile{
        "", "traitsfile", "Traits file for taxon at the leaves", false, "Null", "string", cmd};
    ValueArg<std::string> profiles{
        "c", "profiles", "Preferences profiles (to clamp)", false, "Null", "string", cmd};
    SwitchArg move_root_pop_size{"", "move_root_pop_size", "Clamp the branch mutation rate", cmd, false};
    SwitchArg clamp_pop_sizes{
        "", "clamp_pop_sizes", "Clamp the branch population size", cmd, false};
    SwitchArg clamp_nuc_matrix{"", "clamp_nuc_matrix", "Clamp the nucleotide matrix", cmd, false};
    SwitchArg clamp_corr_matrix{
        "", "clamp_corr_matrix", "Clamp the correlation matrix", cmd, false};
    SwitchArg polymorphism_aware{"p", "polymorphism_aware", "Use polymorphic data", cmd, false};
    ValueArg<unsigned> precision{
        "", "precision", "The precision of PRF computation", false, 6, "unsigned", cmd};
    SwitchArg arithmetic{"d", "arithmetic", "Use arithmetic mean instead of arithmetic", cmd, false};

    void check() {
        if (condition_aware.getValue()) {
            cerr << "The switch parameter ([-b] or [--condition_aware]) is not yet implemented."
                 << endl;
        }
        if (!profiles.getValue().empty()) {
            cout << "Preferences are clamped (option [--profiles <string>]), thus options [--ncat] "
                    "and [-basencat] are not used"
                 << endl;
        }
    }
};

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "DatedMutSel", ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    unique_ptr<DatedNodeMutSelModel> model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        is >> model;
        check_restart(*model, cmd.chain_name() + ".trace");
    } else {
        InferenceAppArgParse args(cmd);
        DatedNodeMutselArgParse datedmutsel_args(cmd);
        cmd.parse();
        datedmutsel_args.check();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.every.getValue(), args.until.getValue());
        model = unique_ptr<DatedNodeMutSelModel>(new DatedNodeMutSelModel(args.alignment.getValue(),
            args.treefile.getValue(), datedmutsel_args.traitsfile.getValue(), datedmutsel_args.profiles.getValue(),
            datedmutsel_args.ncat.getValue(), datedmutsel_args.basencat.getValue(),
            datedmutsel_args.condition_aware.getValue(),
            datedmutsel_args.polymorphism_aware.getValue(), datedmutsel_args.precision.getValue(),
            datedmutsel_args.arithmetic.getValue(), datedmutsel_args.move_root_pop_size.getValue(),
            datedmutsel_args.clamp_pop_sizes.getValue(),
            datedmutsel_args.clamp_nuc_matrix.getValue(),
            datedmutsel_args.clamp_corr_matrix.getValue()));
        model->Update();
    }
    model->ResampleSub(1.0);
    model->MoveParameters(10);
    ConsoleLogger console_logger;
    ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    StandardTracer trace(*model, cmd.chain_name());
    chain_driver->add(*model);
    chain_driver->add(console_logger);
    chain_driver->add(chain_checkpoint);
    chain_driver->add(trace);
    chain_driver->go();
}