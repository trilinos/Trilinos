
#include "NOX.H"  // NOX headers
#include "NOX_LAPACK.H" // NOX LAPACK Interface headers
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

  // This test is only for SERIAL
  // Exit if trying to run in parallel
#ifdef HAVE_MPI
  return 0;
#endif
  
  if (printing.isPrintProcessAndType(NOX::Utils::TestDetails)) {
    cout << "Starting lapack/NOX_Group/NOX_Group.exe" << endl;
  }

  int status = 0;

  // Begin real testing here!
  
  if (printing.isPrintProcessAndType(NOX::Utils::TestDetails)) {
    if (status == 0) 
      cout << "Test was successful!" << endl;
    else 
      cout << "Test Failed!" << endl;
  }

  // 0 is success
  return status;
}

/*
  end of file main.cc
*/
