#include <cmath>
#include <fstream>
#include "InferenceAppArgParse.hpp"
#include "SingleOmegaModel.hpp"
#include "components/ChainCheckpoint.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/StandardTracer.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "SingleOmega", ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    SingleOmegaModel *model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        model = new SingleOmegaModel(is);
    } else {
        InferenceAppArgParse args(cmd);
        cmd.parse();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.every.getValue(), args.until.getValue());
        model = new SingleOmegaModel(args.alignment.getValue(), args.treefile.getValue());
    }

    ConsoleLogger console_logger;
    ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    StandardTracer trace(*model, cmd.chain_name());
    chain_driver->add(*model);
    chain_driver->add(console_logger);
    chain_driver->add(chain_checkpoint);
    chain_driver->add(trace);
    chain_driver->go();

    delete chain_driver;
    delete model;
}
