
#include <stk_util/stk_config.h>
#ifdef STK_HAVE_KOKKOS
#include <Kokkos_Core.hpp>
#endif
#include <iostream>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>

#include "src/test_stk_coupling.hpp"
#include "src/test_stk_search.hpp"
#include "src/test_stk_simd.hpp"
#include "src/test_stk_io.hpp"
#include "src/test_stk_tools.hpp"

int main(int argc, char** argv)
{
  stk::parallel_machine_init(&argc, &argv);

#ifdef STK_HAVE_KOKKOS
  Kokkos::initialize(argc, argv);
  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    std::cout << "Kokkos::DefaultExecutionSpace: " << Kokkos::DefaultExecutionSpace::device_type::execution_space::name() << std::endl;
  }
#endif

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

  test_stk_lib::test_stk_coupling();

  test_stk_lib::test_stk_search();

  test_stk_lib::test_stk_simd(MPI_COMM_WORLD);

  test_stk_lib::test_stk_io(MPI_COMM_WORLD, meshSource, useAutoDecomp);

  test_stk_lib::test_stk_tools();

#ifdef STK_HAVE_KOKKOS
  Kokkos::finalize();
#endif

  stk::parallel_machine_finalize();

  if (proc0) {
    std::cout << "... exiting." << std::endl;
  }
  return 0;
}

