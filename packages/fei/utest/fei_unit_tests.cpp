/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_utils.hpp"
#include "fei_CommUtils.hpp"
#include "fei_test_utils.hpp"

#include "fei_unit_test_runner.hpp"

int main(int argc, char** argv) {

  int numProcs, localProc;
  if (fei_test_utils::initialize_mpi(argc, argv, localProc, numProcs) != 0) {
    return(-1);
  }

  if (fei::localProc(MPI_COMM_WORLD) == 0) {
    FEI_COUT << "\nFEI version: " << fei::utils::version() << "\n\n"<<FEI_ENDL;
  }

  fei::unit::test_runner runner;

  int return_value = runner.run_tests(numProcs, localProc,  MPI_COMM_WORLD);

#ifndef FEI_SER
  MPI_Finalize();
#endif

  return(return_value);
}

