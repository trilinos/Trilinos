#include "mpi.h"
#include <iosfwd>
#include <vector>
#include <stk_util/stk_config.h>

#include <stk_tools/block_extractor/ParseCsv.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/string_utils.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include "ExtractBlocks.hpp"

void print_ids_to_extract(const std::string & idType, const std::vector<int> & ids)
{
  std::cout << "Extracting " << idType << ":";
  for (int id : ids) {
    std::cout << " " << id;
  }
  std::cout << std::endl;
}

std::string get_short_example(std::string execPath)
{
  std::string execName = stk::tailname(execPath);
  std::string mpiCmd = "  > mpirun -n <numProcsDecomp> " + execName + " ";
  std::string quickExample = "Usage\n" + mpiCmd
                             + " --infile " + stk::angle_it("infile")
                             + " --outfile " + stk::angle_it("outfile")
                             + " --extract-blocks " + stk::bracket_it("block IDs")
                             + " --extract-nodesets " + stk::bracket_it("nodeset IDs") + "\n\n";
  return quickExample;
}

std::string get_long_example(std::string execPath)
{
  std::string execName = stk::tailname(execPath);
  std::string mpiCmd = "  > mpirun -n <numProcsDecomp> " + execName + " ";
  std::string example = "Example\n" + mpiCmd
                        + " -i infile.exo"
                        + " -o outfile.exo"
                        + " -b 1,10:15,20:30:2,100"
                        + " -n 1,2,100\n\n";
  return example;
}

int main(int argc, char** argv)
{
  stk::initialize(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  stk::CommandLineParserParallel commandLine(comm);
  commandLine.add_required<std::string>({"infile", "i", "input file"});
  commandLine.add_required<std::string>({"outfile", "o", "output file"});
  commandLine.add_optional<std::string>({"extract-blocks", "b", "Comma-separated list of block IDs or ID ranges of the "
                                         "form first:last:stride. (Example: \"1,10:15,20:30:2,100\")"},{});
  commandLine.add_optional<std::string>({"extract-nodesets","n","Comma-separated list of nodeset IDs or ID ranges"},{});

  std::string quickExample = get_short_example(argv[0]);
  std::string longExample = get_long_example(argv[0]);

  stk::parse_command_line(argc, const_cast<const char**>(argv), quickExample, longExample, commandLine, comm);

  std::string inFile = commandLine.get_option_value<std::string>("infile");
  std::string outFile = commandLine.get_option_value<std::string>("outfile");

  std::string csvBlocks;
  if (commandLine.is_option_provided("extract-blocks")) {
    csvBlocks = commandLine.get_option_value<std::string>("extract-blocks");
  }
  std::vector<int> block_ids = stk::tools::get_ids_from_strings(stk::split_csv_string(csvBlocks));

  std::string csvNodesets;
  if (commandLine.is_option_provided("extract-nodesets")) {
    csvNodesets = commandLine.get_option_value<std::string>("extract-nodesets");
  }
  std::vector<int> nodeset_ids = stk::tools::get_ids_from_strings(stk::split_csv_string(csvNodesets));

  print_ids_to_extract("blocks", block_ids);

  if (!nodeset_ids.empty()) {
    print_ids_to_extract("nodesets", nodeset_ids);
  }

  stk::tools::extract_blocks_and_ns_from_file(inFile, outFile, block_ids, nodeset_ids, comm);

  stk::finalize();
  return 0;
}

