#include "mpi.h"
#include <iosfwd>
#include <vector>
#include <stk_util/stk_config.h>

#if defined(STK_HAVE_BOOSTLIB)

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
    commandLine.add_optional<std::string>({"extract-blocks", "b", "Comma-separated list of blocks to extract"},{});
    commandLine.add_optional<std::string>({"extract-nodesets","n","Comma-separated list of nodesets to extract"},{});

    stk::CommandLineParser::ParseState state = commandLine.parse(argc, argv);
    stk::parallel::require(state == stk::CommandLineParser::ParseComplete, commandLine.get_usage(), comm);

    std::string inFile = commandLine.get_option_value<std::string>("infile");
    std::string outFile = commandLine.get_option_value<std::string>("outfile");
    std::string csvBlocks;
    if (commandLine.is_option_provided("extract-blocks")) {
      csvBlocks = commandLine.get_option_value<std::string>("extract-blocks");
    }
    std::vector<std::string> blockNames = stk::tools::get_block_names_given_ids(stk::tools::get_csv(csvBlocks));

    std::string csvNodesets;
    if (commandLine.is_option_provided("extract-nodesets")) {
      csvNodesets = commandLine.get_option_value<std::string>("extract-nodesets");
    }
    std::vector<std::string> nodesetNames = stk::tools::get_nodeset_names_given_ids(stk::tools::get_csv(csvNodesets));

    if (nodesetNames.empty())
    {
    	std::cout << "No nodesets specified, extracting only desired blocks" << "\n";
        stk::tools::extract_blocks_from_file(inFile, outFile, blockNames, comm);

    }
    else
    {
    	std::cout << "Nodesets specified, extracting combination of blocks and nodesets" << "\n";
    	stk::tools::extract_blocks_and_ns_from_file(inFile,outFile,blockNames, nodesetNames, comm);

    }

    MPI_Finalize();
    return 0;
}

#else

int main(int argc, const char**argv)
{
  std::cerr<<"ERROR, stk_block_extractor requires that Trilinos was configured with -DTPL_ENABLE_BoostLib:BOOL=ON"
           << std::endl;
  return -1;
}

#endif

