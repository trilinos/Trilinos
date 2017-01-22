#include <Kokkos_OpenMP.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

#include "ShyLUTacho_config.h"

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::OpenMP exec_space;

#ifdef HAVE_SHYLUTACHO_MKL
#include "Tacho_ExampleCholPardiso.hpp"
using namespace Tacho;
#endif

int main (int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program measure the performance of Pardiso Chol algorithms on Kokkos::OpenMP execution space.\n");

  int nthreads = 0;
  clp.setOption("nthreads", &nthreads, "Number of threads");

  bool verbose = false;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  std::string file_input = "test.mtx";
  clp.setOption("file-input", &file_input, "Input file (MatrixMarket SPD matrix)");

  int nrhs = 1;
  clp.setOption("nrhs", &nrhs, "Number of RHS vectors");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  int r_val = 0;
  {
    exec_space::initialize(nthreads);
    exec_space::print_configuration(cout, true);

#ifdef HAVE_SHYLUTACHO_MKL
    const bool skip_factorize = false;
    const bool skip_solve     = false;
    r_val = examplePardisoChol<exec_space>
      (file_input,
       nrhs, 
       nthreads,
       skip_factorize,
       skip_solve,
       verbose);
#else
    r_val = -1;
    cout << "MKL is NOT configured in Trilinos" << endl;
#endif

    exec_space::finalize();
  }

  return r_val;
}
