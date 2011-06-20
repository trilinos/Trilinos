#include <iostream>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp> // TODO: rename MueLu Gallery

/*
  This driver simply generates a Tpetra matrix, prints it to screen, and exits.

  THIS EXAMPLE DOES NOT USE CTHULHU.

  Use the "--help" option to get verbose help.
*/

int main(int argc, char** argv) 
{
  using Teuchos::RCP;

  typedef int    LO; // LocalOrdinal
  typedef int    GO; // GlobalOrdinal
  typedef double SC; // Scalar

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);
  
  MueLu::Gallery::Parameters<GO> matrixParameters(clp);   // manage parameters of the test case
  
  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }
  
  matrixParameters.check();
  matrixParameters.print();

  /**********************************************************************************/
  /* CREATE INITAL MATRIX                                                           */
  /**********************************************************************************/
  RCP<const Tpetra::Map<LO,GO> > map = rcp( new Tpetra::Map<LO,GO>(matrixParameters.GetNumGlobalElements(), 0, comm) );
  RCP<Tpetra::CrsMatrix<SC,LO,GO> > A = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Tpetra::Map<LO,GO>, Tpetra::CrsMatrix<SC,LO,GO> >(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());

  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  if (comm->getRank() == 0)
    std::cout << "\n================ MAP =====================================================\n" << std::endl;
  map->describe(*out, Teuchos::VERB_EXTREME);
  comm->barrier();
  sleep(1);

  if (comm->getRank() == 0)
    std::cout << "\n================ MATRIX ==================================================\n" << std::endl;
  A->describe(*out, Teuchos::VERB_EXTREME);

  return EXIT_SUCCESS;
} 
