// Trilinos headers 
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

// ParaCont headers
#include "ContinuationManager.H"
#include "LinearSystem.H"

// Main driver
int main( int argc, char **argv )
{

  try {
  
#ifdef HAVE_MPI
    // Initialise MPI
    MPI_Init(&argc, &argv);

    // Create a communicator
    Teuchos::RefCountPtr <Epetra_MpiComm> comm = 
      Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    // Create a communicator
    Teuchos::RefCountPtr <Epetra_SerialComm> comm = 
      Teuchos::rcp(new Epetra_SerialComm);
#endif

    std::string fileName = "task.xml";
    if (argc>1) 
       fileName = argv[1];

    // Instantiate the continuation manager
    Teuchos::RefCountPtr <ContinuationManager> contManager = 
      Teuchos::rcp(new ContinuationManager(comm,fileName));

    // Instantiate the problem
    Teuchos::RefCountPtr <LinearSystem> problem = 
      Teuchos::rcp(new LinearSystem(comm)); 

    // Set the problem in the continuation manager
    contManager->SetLOCAProblem(problem);

    // Prepare to run LOCA
    contManager->BuildLOCAStepper();

    // Run LOCA
    bool status = contManager->RunLOCAStepper();

  if (status)
    std::cout << "\nAll tests passed" << std::endl;

  }

  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  catch (const char *s) {
    std::cout << s << std::endl;
  }

  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
  }

#ifdef HAVE_MPI
  // Finalise MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);

}
