#include "mpi.h"
#include <stk_balance/balance.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/Inputs.hpp>
#include <stk_balance/internal/balanceCommandLine.hpp>
#include <stk_balance/internal/balanceDefaults.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/command_line/CommandLineParser.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/environment/FileUtils.hpp>

#include <string>
#include <iostream>
#include <fstream>

#include "stk_util/environment/EnvData.hpp"

int main(int argc, const char**argv)
{
    MPI_Init(&argc, const_cast<char***>(&argv));
    MPI_Comm comm = MPI_COMM_WORLD;

    if(stk::parallel_machine_rank(comm) != 0)
    {
        stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
    }

    std::string execName = stk::util::tailname(argv[0]);
    stk::balance::CommandLineOptions cmdLineOptions;
    std::string quickExample = stk::balance::get_quick_example(execName, cmdLineOptions.infile.name, cmdLineOptions.outputDirectory.name, comm);
    stk::balance::ParsedOptions balanceOptions;
    try
    {
        balanceOptions = stk::balance::parse_balance_command_line(argc,
                                                                  argv,
                                                                  execName,
                                                                  comm);
    }
    catch(const std::exception &e)
    {
        std::string errorMessage = e.what() + stk::get_quick_error(execName, quickExample);
        stk::parallel::print_and_exit(errorMessage, comm);
    }

    stk::parallel::require_file_exists(balanceOptions.inFile, execName, quickExample, comm);

    stk::balance::print_running_msg(execName, balanceOptions, comm);
    try {
        stk::balance::run_stk_rebalance(balanceOptions.outputDirectory, balanceOptions.inFile, balanceOptions.appTypeDefaults, comm);
    }
    catch(std::exception& e)
    {
        std::cerr<<e.what()<<std::endl;
    }

    MPI_Finalize();
    return 0;
}


