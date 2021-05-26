#include "Tacho_CommandLineParser.hpp" 
#include "Tacho_ExampleDeviceDenseCholesky.hpp"
#include "Tacho_ExampleDeviceDenseLDL.hpp"

int main (int argc, char *argv[]) {
  Tacho::CommandLineParser opts("This example program measure the Tacho on Kokkos::OpenMP");

  int m = 10;
  bool verbose = false;
  bool test_chol = false;
  bool test_ldl = false;
  //bool test_lu = false;

  opts.set_option<int>("m", "Dense problem size", &m);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<bool>("test-chol", "Flag for testing Cholesky", &test_chol);
  opts.set_option<bool>("test-ldl", "Flag for testing LDL", &test_ldl);
  //  opts.set_option<bool>("test-lu", "Flag for testing LU", &test_lu);


  const bool r_parse = opts.parse(argc, argv);
  if (r_parse) return 0; // print help return

  int r_val(0);
  Kokkos::initialize(argc, argv);
  {
    if (test_chol) {
      const int r_val_chol = driver_chol<double>(m, verbose); r_val += r_val_chol;
    }
    if (test_ldl) {
      const int r_val_ldl  = driver_ldl <double>(m, verbose); r_val += r_val_ldl;
    }
  }
  Kokkos::finalize();
  
  return r_val;
}
