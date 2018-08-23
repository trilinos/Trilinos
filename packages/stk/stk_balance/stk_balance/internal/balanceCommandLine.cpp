#include <internal/balanceCommandLine.hpp>
#include <internal/Inputs.hpp>

#include "mpi.h"

#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/command_line/CommandLineParser.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/environment/FileUtils.hpp>

#include <string>
#include <iostream>

namespace
{

stk::balance::CommandLineOptions cmdLineOptions;

}


namespace stk {
namespace balance {

enum { OK=0, NOT_OK = 1 };

bool am_proc0(MPI_Comm comm)
{
    return stk::parallel_machine_rank(comm)==0;
}

int handle_output_directory(const std::string& outputDirectory, MPI_Comm comm)
{
    bool path_ok = true;
    if(am_proc0(comm))
        path_ok = stk::balance::create_path(outputDirectory);

    if(!stk::is_true_on_all_procs(comm, path_ok))
    {
        return NOT_OK;
    }
    return OK;
}

std::string get_quick_example(const std::string &execName,
                              const std::string &infileName,
                              const std::string &outputDirectoryName,
                              MPI_Comm comm)
{
    std::string mpiCmd = "  > mpirun -n <numProcsDecomp> " + execName + " ";
    std::string usage = "Usage:\n" + mpiCmd + stk::angle_it(infileName) + " " + stk::bracket_it(outputDirectoryName)
                            + "\n" + mpiCmd + stk::dash_it(infileName) + " " + stk::angle_it(infileName)
                            + " " + stk::dash_it(outputDirectoryName) + " " + stk::bracket_it(outputDirectoryName)
                  + "\n\n";
    return usage;
}

std::string get_examples(const std::string &executableName)
{
    std::string examples = "Examples:\n\n";
    std::string tab = "  ";
    examples += tab + "To decompose for 16 processors:\n";
    examples += tab + tab + "> mpirun -n 16 " + executableName + " file.exo\n";
    examples += "\n";
    examples += tab + "To decompose for 512 processors and put the decomposition into a directory named 'temp1':\n";
    examples += tab + tab + "> mpirun -n 512 " + executableName + " file.exo temp1\n";
    examples += "\n";
    examples += tab + "To decompose for 64 processors and use settings suitable for solving Solid Mechanics problems':\n";
    examples += tab + tab + "> mpirun -n 64 " + executableName + " file.exo --sm\n";
    return examples;
}

void parse_single_app_defaults(const stk::CommandLineParserParallel& commandLine,
                               const std::string &optionName,
                               enum AppTypeDefaults appType,
                               ParsedOptions &balanceOptions)
{
    if(commandLine.is_option_provided(optionName))
    {
        if(balanceOptions.appTypeDefaults != NO_DEFAULTS)
            throw std::invalid_argument("Can't set default settings for multiple apps at the same time.");
        balanceOptions.appTypeDefaults = appType;
    }
}

void parse_app_defaults(const stk::CommandLineParserParallel& commandLine, ParsedOptions& balanceOptions)
{
    parse_single_app_defaults(commandLine, cmdLineOptions.smDefaults.name, SM_DEFAULTS, balanceOptions);
    parse_single_app_defaults(commandLine, cmdLineOptions.sdDefaults.name, SD_DEFAULTS, balanceOptions);
}

ParsedOptions parse_balance_command_line(int argc,
                                 const char**argv,
                                 const std::string &execName,
                                 MPI_Comm comm)
{
    stk::CommandLineParserParallel commandLine(comm);
    commandLine.add_required_positional<std::string>(cmdLineOptions.infile);
    commandLine.add_optional_positional<std::string>(cmdLineOptions.outputDirectory, ".");
    commandLine.add_flag(cmdLineOptions.smDefaults.name, cmdLineOptions.smDefaults.description);
    commandLine.add_flag(cmdLineOptions.sdDefaults.name, cmdLineOptions.sdDefaults.description);

    std::string quickExample = stk::balance::get_quick_example(execName, cmdLineOptions.infile.name, cmdLineOptions.outputDirectory.name, comm);

    stk::parse_command_line(argc, argv, quickExample, get_examples(execName), commandLine, comm);
    ParsedOptions balanceOptions{commandLine.get_option_value<std::string>(cmdLineOptions.infile.name),
                                 commandLine.get_option_value<std::string>(cmdLineOptions.outputDirectory.name),
                                 NO_DEFAULTS};
    parse_app_defaults(commandLine, balanceOptions);
    stk::parallel::require(handle_output_directory(balanceOptions.outputDirectory, comm) == OK, "Unable to create output directory.", comm);

    return balanceOptions;
}

void print_running_msg(const std::string &execName, const ParsedOptions &balanceOptions, MPI_Comm comm)
{
    if(am_proc0(comm))
    {
        std::cerr << "\n";
        std::cerr << "\tRunning: " << execName << " " << balanceOptions.inFile << " " << balanceOptions.outputDirectory << std::endl;
        std::cerr << "\n";
    }
}

}
}
