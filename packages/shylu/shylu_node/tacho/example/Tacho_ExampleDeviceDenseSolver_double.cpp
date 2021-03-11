#include "Tacho_CommandLineParser.hpp" 
#include "Tacho_ExampleDeviceDenseCholesky.hpp"
#include "Tacho_ExampleDeviceDenseLDL.hpp"

int main (int argc, char *argv[]) {
  Tacho::CommandLineParser opts("This example program measure the Tacho on Kokkos::OpenMP");

  int m = 10;
  bool verbose = false;

  opts.set_option<int>("m", "Dense problem size", &m);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);


  const bool r_parse = opts.parse(argc, argv);
  if (r_parse) return 0; // print help return

  int r_val(0);
  Kokkos::initialize(argc, argv);
  {
    const int r_val_chol = 0; //driver_chol<double>(m, verbose);
    const int r_val_ldl  = driver_ldl <double>(m, verbose);
    r_val = r_val_chol + r_val_ldl;
  }
  Kokkos::finalize();
  
  return r_val;
}
