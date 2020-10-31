
#include <iostream>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>

#include "src/test_stk_simd.hpp"
#include "src/test_stk_io.hpp"

int main(int argc, char** argv)
{
  if (MPI_SUCCESS != MPI_Init(&argc, &argv)) {
    std::cout << "MPI_Init failed." << std::endl;
    return -1;
  }

  const bool proc0 = (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0);

  stk::CommandLineParserParallel clp(MPI_COMM_WORLD);

  clp.add_optional<std::string>({"mesh", "m", "mesh file or generated-mesh description"},
                                "generated:5x20x100");

  clp.add_flag({"autodecomp", "a", "read serial mesh file into parallel run"});

  stk::CommandLineParser::ParseState parseResult = clp.parse(argc, (const char**)argv);

  if (stk::CommandLineParser::ParseError == parseResult) {
    std::cout << "Exiting early due to command-line-parser error." << std::endl;
    return -1;
  }

  std::string meshSource = clp.get_option_value<std::string>("mesh");
  const bool useAutoDecomp = clp.is_option_provided("autodecomp");

  if (proc0) {
    std::cout << "Test-STK-App" << std::endl;
  }

  test_stk_lib::test_stk_simd(MPI_COMM_WORLD);

  test_stk_lib::test_stk_io(MPI_COMM_WORLD, meshSource, useAutoDecomp);

  MPI_Finalize();

  if (proc0) {
    std::cout << "... exiting." << std::endl;
  }
  return 0;
}

