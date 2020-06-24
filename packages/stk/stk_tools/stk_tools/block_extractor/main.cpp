#include "mpi.h"
#include <iosfwd>
#include <vector>
#include <stk_util/stk_config.h>

#include <stk_tools/block_extractor/ParseCsv.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include "ExtractBlocks.hpp"

void print_ids_to_extract(const std::string & idType, const std::vector<int> & ids)
{
  std::cout << "Extracting " << idType << ":";
  for (int id : ids) {
    std::cout << " " << id;
  }
  std::cout << std::endl;
}

int main(int argc, const char**argv)
{
    MPI_Init(&argc, const_cast<char***>(&argv));
    MPI_Comm comm = MPI_COMM_WORLD;

    stk::CommandLineParserParallel commandLine(comm);
    commandLine.add_required<std::string>({"infile", "i", "input file"});
    commandLine.add_required<std::string>({"outfile", "o", "output file"});
    commandLine.add_optional<std::string>({"extract-blocks", "b", "Comma-separated list of block IDs or ID ranges of the "
                                           "form first:last:stride. (Example: \"1,10:15,20:30:2,100\")"},{});
    commandLine.add_optional<std::string>({"extract-nodesets","n","Comma-separated list of nodeset IDs or ID ranges"},{});

    stk::CommandLineParser::ParseState state = commandLine.parse(argc, argv);
    ThrowRequireMsg(state == stk::CommandLineParser::ParseComplete, commandLine.get_usage());

    std::string inFile = commandLine.get_option_value<std::string>("infile");
    std::string outFile = commandLine.get_option_value<std::string>("outfile");

    std::string csvBlocks;
    if (commandLine.is_option_provided("extract-blocks")) {
      csvBlocks = commandLine.get_option_value<std::string>("extract-blocks");
    }
    std::vector<int> block_ids = stk::tools::get_ids_from_strings(stk::tools::get_csv(csvBlocks));

    std::string csvNodesets;
    if (commandLine.is_option_provided("extract-nodesets")) {
      csvNodesets = commandLine.get_option_value<std::string>("extract-nodesets");
    }
    std::vector<int> nodeset_ids = stk::tools::get_ids_from_strings(stk::tools::get_csv(csvNodesets));

    print_ids_to_extract("blocks", block_ids);

    if (!nodeset_ids.empty()) {
      print_ids_to_extract("nodesets", nodeset_ids);
    }

    stk::tools::extract_blocks_and_ns_from_file(inFile, outFile, block_ids, nodeset_ids, comm);

    MPI_Finalize();
    return 0;
}

