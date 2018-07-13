#include "mpi.h"
#include <iosfwd>
#include <vector>
#include <stk_tools/block_extractor/ParseCsv.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include "ExtractBlocks.hpp"

int main(int argc, const char**argv)
{
    MPI_Init(&argc, const_cast<char***>(&argv));
    MPI_Comm comm = MPI_COMM_WORLD;

    stk::CommandLineParserParallel commandLine(comm);
    commandLine.add_required<std::string>({"infile", "i", "input file"});
    commandLine.add_required<std::string>({"outfile", "o", "output file"});
    commandLine.add_required<std::string>({"extract-blocks", "b", "Comma-separated list of blocks to extract"});

    stk::CommandLineParser::ParseState state = commandLine.parse(argc, argv);
    stk::parallel::require(state == stk::CommandLineParser::ParseComplete, commandLine.get_usage(), comm);

    std::string inFile = commandLine.get_option_value<std::string>("infile");
    std::string outFile = commandLine.get_option_value<std::string>("outfile");
    std::string csvBlocks = commandLine.get_option_value<std::string>("extract-blocks");
    std::vector<std::string> blockNames = stk::tools::get_block_names_given_ids(stk::tools::get_csv(csvBlocks));
    stk::tools::extract_blocks_from_file(inFile, outFile, blockNames, comm);

    MPI_Finalize();
    return 0;
}


