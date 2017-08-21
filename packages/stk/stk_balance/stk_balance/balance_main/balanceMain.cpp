#include "mpi.h"
#include <stk_balance/balance.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/Inputs.hpp>

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

namespace
{

struct CommandLineOptions
{
    stk::CommandLineOption infile{"infile", "i", "undecomposed serial input mesh file"};
    stk::CommandLineOption outputDirectory{"outputDirectory", "o", "output directory for decomposition"};
};

CommandLineOptions cmdLineOptions;

struct ParsedOptions
{
    std::string inFile;
    std::string outputDirectory;
};

}

bool create_path(const std::string &filename);
bool am_proc0();
void print_running_msg(const std::string &execName, const ParsedOptions &balanceOptions);

enum { OK=0, NOT_OK = 1 };

int handle_output_directory(const std::string& outputDirectory)
{
    bool path_ok = true;
    if(am_proc0())
        path_ok = stk::balance::create_path(outputDirectory);

    if(!stk::is_true_on_all_procs(MPI_COMM_WORLD, path_ok))
    {
        return NOT_OK;
    }
    return OK;
}

std::string get_quick_example(const std::string &execName, MPI_Comm comm)
{
    std::string infileName = cmdLineOptions.infile.name;
    std::string outputDirectoryName = cmdLineOptions.outputDirectory.name;
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
    return examples;
}

ParsedOptions parse_balance_command_line(int argc,
                                 const char**argv,
                                 const std::string &quickExample,
                                 stk::CommandLineParserParallel &commandLine,
                                 MPI_Comm comm)
{
    std::string execName = stk::util::tailname(argv[0]);
    stk::parse_command_line(argc, argv, quickExample, get_examples(execName), commandLine, comm);
    ParsedOptions balanceOptions{commandLine.get_option_value<std::string>(cmdLineOptions.infile.name),
                                 commandLine.get_option_value<std::string>(cmdLineOptions.outputDirectory.name)};
    stk::parallel::require(handle_output_directory(balanceOptions.outputDirectory) == OK, "Unable to create output directory.", comm);
    stk::parallel::require_file_exists(balanceOptions.inFile, execName, get_quick_example(execName, comm), comm);

    return balanceOptions;
}


int main(int argc, const char**argv)
{
    MPI_Init(&argc, const_cast<char***>(&argv));
    MPI_Comm comm = MPI_COMM_WORLD;

    stk::CommandLineParserParallel commandLine(comm);
    commandLine.add_required_positional<std::string>(cmdLineOptions.infile);
    commandLine.add_optional_positional<std::string>(cmdLineOptions.outputDirectory, ".");

    std::string execName = stk::util::tailname(argv[0]);
    ParsedOptions balanceOptions = parse_balance_command_line(argc, argv, get_quick_example(execName, comm), commandLine, comm);

    print_running_msg(execName, balanceOptions);
    stk::balance::run_stk_rebalance(balanceOptions.outputDirectory, balanceOptions.inFile, comm);

    MPI_Finalize();
    return OK;
}

bool am_proc0()
{
    return stk::parallel_machine_rank(MPI_COMM_WORLD)==0;
}

void print_running_msg(const std::string &execName, const ParsedOptions &balanceOptions)
{
    if(am_proc0())
    {
        std::cerr << "\n";
        std::cerr << "\tRunning: " << execName << " " << balanceOptions.inFile << " " << balanceOptions.outputDirectory << std::endl;
        std::cerr << "\n";
    }
}



