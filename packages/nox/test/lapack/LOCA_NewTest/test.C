
#include "LOCA.H"  // LOCA headers
#include "LOCA_LAPACK.H" // LOCA LAPACK Interface headers
#include "nox_test_err.h" // common file for testing 

#ifdef HAVE_MPI
#include <mpi.h>
#else 
#endif

int main(int argc, char *argv[]) {

  // Set up the printing utilities
  NOX::Parameter::List noxParams;
  NOX::Parameter::List& printParams = noxParams.sublist("Printing");
  printParams.setParameter("Output Precision", 5);
  if (argc > 1) { 
    if (argv[1][0]=='-' && argv[1][1]=='v')
       printParams.setParameter("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Parameters + 
			NOX::Utils::Details + 
			NOX::Utils::Warning +
			NOX::Utils::TestDetails);
    else
       printParams.setParameter("Output Information", NOX::Utils::Error);
  }
  NOX::Utils printing(printParams);

  // This test (and lapack routines in general) is only for SERIAL
  // Exit if trying to run in parallel
#ifdef HAVE_MPI
  return 0;
#endif
  
  // Identify the test
  if (printing.isPrintProcessAndType(NOX::Utils::TestDetails)) {
    cout << "Starting lapack/LOCA_NewTest/LOCA_NewTest.exe" << endl;
  }

  // Final return value (0 = successfull, non-zero = failure)
  int status = 0;

  // *** Insert your testing here! ***



  // Summarize test results  
  if (printing.isPrintProcessAndType(NOX::Utils::TestDetails)) {
    if (status == 0) 
      cout << "Test was successful!" << endl;
    else 
      cout << "Test Failed!" << endl;
  }

  // Final return value (0 = succefull, non-zero = failure)
  return status;
}

/*
  end of file main.cc
*/
