#include "mpi.h"
#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/balanceMtoN.hpp>
#include <stk_balance/internal/Inputs.hpp>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_util/command_line/CommandLineParser.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/environment/FileUtils.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/util/string_utils.hpp>

#include <string>
#include <iostream>
#include <fstream>

namespace
{

struct CommandLineOptions
{
    stk::CommandLineOption infile{"infile", "i", "input file decomposed for old number of processors"};
    stk::CommandLineOption nprocs{"nprocs", "n", "new number of processors"};
};

CommandLineOptions m2nOptions;

struct ParsedOptions
{
    std::string inFile;
    int targetNumProcs;
};

std::string get_quick_example(const std::string &execName, MPI_Comm comm)
{
    std::string nProcs = std::to_string(stk::parallel_machine_size(comm));
    std::string infileName = m2nOptions.infile.name;
    std::string nProcName = m2nOptions.nprocs.name;
    std::string mpiCmd = "  > mpirun -n " + nProcs + " " + execName + " ";
    std::string usage = "Usage:\n" + mpiCmd + stk::angle_it(infileName) + " " + stk::angle_it(nProcName)
                            + "\n" + mpiCmd + stk::dash_it(infileName) + " " + stk::angle_it(infileName)
                            + " " + stk::dash_it(nProcName) + " " + stk::angle_it(nProcName)
                  + "\n\n";
    return usage;
}

std::string get_examples(const std::string &executableName)
{
    std::string examples = "Examples:\n\n";
    std::string tab = "  ";
    examples += tab + "To change from 16 processor decomposition to 128 processor decomposition:\n";
    examples += tab + tab + "> mpirun -n 16 " + executableName + " file.exo 128\n";
    examples += "\n";
    examples += tab + "To change from 16 processor decomposition to 8 processor decomposition:\n";
    examples += tab + tab + "> mpirun -n 16 " + executableName + " file.exo 8\n";
    return examples;
}

ParsedOptions parse_m2n_command_line(int argc, const char**argv, stk::CommandLineParserParallel &commandLine, MPI_Comm comm)
{
    std::string execName = stk::tailname(argv[0]);
    stk::parse_command_line(argc, argv, get_quick_example(execName, comm), get_examples(execName), commandLine, comm);

    std::string inFile = commandLine.get_option_value<std::string>(m2nOptions.infile.name);
    int targetNumProcs = commandLine.get_option_value<int>(m2nOptions.nprocs.name);
    ThrowRequireMsg(targetNumProcs > 0, "Please specify a valid target processor count.");

    return ParsedOptions{inFile, targetNumProcs};
}

void rebalance_m_to_n(ParsedOptions &parsedOptions, MPI_Comm comm)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, comm);

    stk::mesh::Field<double> &field = meta.declare_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK, "TargetDecomp", 1);
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), (double*)nullptr);

    stk::io::StkMeshIoBroker ioBroker;
    stk::io::fill_mesh_preexisting(ioBroker, parsedOptions.inFile, bulk);

    stk::balance::internal::rebalanceMtoN(ioBroker, field, parsedOptions.targetNumProcs, parsedOptions.inFile);
}

}

int main(int argc, const char**argv)
{
    MPI_Init(&argc, const_cast<char***>(&argv));
    MPI_Comm comm = MPI_COMM_WORLD;

    stk::CommandLineParserParallel commandLine(comm);
    commandLine.add_required_positional<std::string>(m2nOptions.infile);
    commandLine.add_required_positional<int>(m2nOptions.nprocs);

    ParsedOptions balanceOptions = parse_m2n_command_line(argc, argv, commandLine, comm);
    rebalance_m_to_n(balanceOptions, comm);

    MPI_Finalize();
    return 0;
}
