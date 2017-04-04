#include <mpi.h>
#include <iosfwd>
#include <vector>
#include <stk_util/environment/CommandLineParser.hpp>
#include <stk_tools/block_extractor/ParseCsv.hpp>
#include "ExtractBlocks.hpp"

int main(int argc, const char**argv)
{
    MPI_Init(&argc, const_cast<char***>(&argv));
    MPI_Comm comm = MPI_COMM_WORLD;

    stk::CommandLineParser commandLine;
    commandLine.add_option<std::string>("infile,i", "input file");
    commandLine.add_option<std::string>("outfile,o", "output file");
    commandLine.add_option<std::string>("extract-blocks,b", "Comma-separated list of blocks to extract");
    commandLine.parse(argc, argv);

    std::string inFile = commandLine.get_option_value<std::string>("infile");
    std::string outFile = commandLine.get_option_value<std::string>("outfile");
    std::string csvBlocks = commandLine.get_option_value<std::string>("extract-blocks");
    std::vector<std::string> blockNames = stk::tools::get_block_names_given_ids(stk::tools::get_csv(csvBlocks));
    stk::tools::extract_blocks_from_file(inFile, outFile, blockNames, comm);

    MPI_Finalize();
    return 0;
}


