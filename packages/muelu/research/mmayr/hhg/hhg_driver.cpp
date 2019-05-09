#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include "Teuchos_Assert.hpp"
#include "Teuchos_LocalTestingHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

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

/* Read xmlFile with all possible problem definitions and extract the one to be
 * solved.
 *
 * All parameters are stored in \c xmlParams, while the problem-specific
 * parameters are stored in \c probParams.
 *
 * \return probParams
 */
Teuchos::RCP<const Teuchos::ParameterList>
readProblemParams(const std::string xmlFileName,
    Teuchos::RCP<const Teuchos::Comm<int> > comm)
{
  // check for reasonable filename of xml-file to be read
  Teuchos::RCP<Teuchos::ParameterList> xmlParams = Teuchos::null;
  if (xmlFileName.length() && xmlFileName.rfind(".xml"))
  {
    // read the parameter list
    xmlParams = Teuchos::getParametersFromXmlFile(xmlFileName.c_str());
//    xmlParams->print(std::cout);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "The file name " << xmlFileName << " is not a valid XML file name.");

  // Access sublist with data of current problem to be solved
  const std::string probName = xmlParams->get<std::string>("problem type");
  if (comm->getRank() == 0)
    std::cout << "Reading parameter list for problem '"
        << probName << "'." << std::endl;

  const Teuchos::ParameterList& probParams = xmlParams->sublist(probName);

  // Perform some sanity checks
  TEUCHOS_TEST_FOR_EXCEPT_MSG(probParams.get<int>("number of processors") != comm->getSize(),
      "Number of processes defined in input file (" << probParams.get<int>("number of processors")
      << ") does not match number of MPI ranks (" << comm->getSize() << ").");

  // Wrap into RCP to return
  Teuchos::RCP<const Teuchos::ParameterList> probParamsRCP = Teuchos::rcp(new Teuchos::ParameterList(probParams));

  return probParamsRCP;
}

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

//  TEUCHOS_TEST_FOR_EXCEPT_MSG(argc != 1, "Do only provide an input xml-file, nothing else.");

  // process input arguments
  std::string xmlFileName = argv[1]; // xml-file to configure MueLu

  // check for reasonable filename of xml-file to be read
  Teuchos::RCP<const Teuchos::ParameterList> probParams = readProblemParams(xmlFileName, comm);
  if (comm->getRank() == 0)
    probParams->print(*out);

  const std::string& mappingFileName = probParams->get<std::string>("mapping file name");
  // const std::string& matrixFileName = probParams->get<std::string>("matrix file name");

  // create the RegionManager to deal with node-to-region mappings
  Teuchos::RCP<RegionManager> regionManager = Teuchos::rcp(new RegionManager(mappingFileName, comm));
//  regionManager->printNodeRegionPairs(*out);
//  regionManager->printNodesToRegionMapping(*out);
//  regionManager->printInterfaceNodesToRegionMapping(*out);
//  regionManager->printInactiveProcs(*out);
////  regionManager->printNumRegionsPerProc(*out);
////  regionManager->printProcsPerRegion(*out);
//

//  // create the RegionMatrix to access the assembled, the composite, and the regional matrix
//  Teuchos::RCP<RegionMatrix> regionMatrix = Teuchos::rcp(new RegionMatrix(matrixFileName, regionManager, comm));
//  regionMatrix->printCompositeMatrix(*out, Teuchos::VERB_EXTREME);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

