// Trilinos headers 
#ifdef HAVE_MPI
#include "mpi.h"
#include "EpetraExt_MultiMpiComm.h"
#else
#include "EpetraExt_MultiSerialComm.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

// ParaCont headers
#include "ContinuationManager.H"
#include "PeriodicLinearSystem.H"

// Main driver
int main( int argc, char **argv )
{

  try {
  
#ifdef HAVE_MPI
    // Initialise MPI
    MPI_Init(&argc, &argv);

    // Create a two-level communicator for space&time parallelism
    int numSpatialProcs = 1;
    int numTimeSteps = 2;
    Teuchos::RefCountPtr <EpetraExt::MultiMpiComm> globalComm = 
      Teuchos::rcp(new EpetraExt::MultiMpiComm(MPI_COMM_WORLD,
                                              numSpatialProcs,
                                              numTimeSteps));

    // Get the communicator for the spatial sub-problem
    //Teuchos::RefCountPtr <Epetra_MpiComm> comm = 
    //   Teuchos::rcp(&(globalComm->SubDomainComm()), false);
#else
    // Create a two-level communicator for space&time parallelism
    int numTimeSteps = 2;
    Teuchos::RefCountPtr <EpetraExt::MultiSerialComm> globalComm = 
      Teuchos::rcp(new EpetraExt::MultiSerialComm(numTimeSteps));

    // Get the communicator for the spatial sub-problem
   // Teuchos::RefCountPtr <Epetra_SerialComm> comm = 
    //   Teuchos::rcp(&(globalComm->SubDomainComm()), false);
#endif
    Teuchos::RefCountPtr <Epetra_Comm> comm = 
       Teuchos::rcp(&(globalComm->SubDomainComm()), false);

    std::string fileName = "task.xml";
    if (argc>1) 
       fileName = argv[1];

    // Instantiate the continuation manager
    Teuchos::RefCountPtr <ContinuationManager> contManager = 
      Teuchos::rcp(new ContinuationManager(comm,fileName));

    // Instantiate the problem
    Teuchos::RefCountPtr <PeriodicLinearSystem> problem = 
      Teuchos::rcp(new PeriodicLinearSystem(comm)); 

    // Set the problem in the continuation manager
    contManager->SetLOCAProblem(problem);

    // Prepare to run LOCA
    contManager->BuildLOCAPeriodicStepper(globalComm);

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
