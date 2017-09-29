#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include "Teuchos_Assert.hpp"
#include "Teuchos_LocalTestingHelpers.hpp"
#include <Teuchos_StandardCatchMacros.hpp>
#include "Teuchos_ParameterList.hpp"

// Xpetra
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_IO.hpp"
#include "Xpetra_RegionManager_impl.hpp"
#include "Xpetra_RegionMatrix_impl.hpp"
#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraCrsMatrix.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraCrsMatrix.hpp"
#endif

// Epetra
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

// MueLu
#include "MueLu_Utilities.hpp"


// =========== //
// main driver //
// =========== //

int main(int argc, char* argv[])
{
  typedef double                                      scalar_type;
  typedef int                                         local_ordinal_type;
  typedef int                                         global_ordinal_type;
  typedef scalar_type         												Scalar;
  typedef local_ordinal_type  												LocalOrdinal;
  typedef global_ordinal_type 												GlobalOrdinal;
  typedef KokkosClassic::DefaultNode::DefaultNodeType Node;
  typedef Scalar SC;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Node NO;

  typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Xpetra::RegionManager<SC,LO,GO,NO> RegionManager;
  typedef Xpetra::RegionMatrix<SC,LO,GO,NO,Xpetra::UseTpetra,Xpetra::region_split> RegionMatrix;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm CommEpetra(MPI_COMM_WORLD);
#else
  Epetra_SerialComm CommEpetra;
#endif

  // wrap communicator into Teuchos::Comm
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

//  //  MueLu::Utilities<SC,LO,GO,NO>::PauseForDebugger();
//  if (comm->getRank() == 0) {
//    std::cout << "** Enter a character to continue > " << std::endl;
//    std::cin.get();
//  }
//  comm->barrier();

  if (comm->getSize() > 1)
    MueLu::Utilities<SC,LO,GO,NO>::PauseForDebugger();


  // wrap the output stream
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  TEUCHOS_TEST_FOR_EXCEPT_MSG(argc<4, "\nNot enough command line arguments!\n");

  // process input arguments
  std::string xmlFileName = argv[1]; // xml-file to configure MueLu
  std::string matrixFileName = argv[2]; // file with matrix entries
  std::string mappingFileName = argv[3]; // file with node-to-region mapping

  // create the RegionManager to deal with node-to-region mappings
  Teuchos::RCP<RegionManager> regionManager = Teuchos::rcp(new RegionManager(mappingFileName, comm));
  regionManager->printNodeRegionPairs(*out);
  regionManager->printNodesToRegionMapping(*out);
  regionManager->printInterfaceNodesToRegionMapping(*out);
  regionManager->printInactiveProcs(*out);
//  regionManager->printNumRegionsPerProc(*out);
//  regionManager->printProcsPerRegion(*out);


  // create the RegionMatrix to access the assembled, the composite, and the regional matrix
  Teuchos::RCP<RegionMatrix> regionMatrix = Teuchos::rcp(new RegionMatrix(matrixFileName, regionManager, comm));

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

