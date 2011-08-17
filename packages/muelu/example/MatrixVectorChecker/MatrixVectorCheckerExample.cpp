// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Using Galeri:
// #include <Epetra_Maps.h>
// #include <Epetra_CrsMatrix.h>
// #include <Galeri_Maps.h>
// #include <Galeri_CrsMatrices.h>

// Using MueLu:
#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_CrsMatrix.hpp>

#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>

// MueLu Gallery :
#define XPETRA_ENABLED
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

//
#include "MatrixVectorChecker.hpp"

int main(int argc, char *argv[])
{
  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv, &blackhole);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);
  
  MueLu::Gallery::Parameters<GO> matrixParameters(clp);   // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);         // manage parameters of xpetra
  
  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }
  
  matrixParameters.check();
  xpetraParameters.check();

  matrixParameters.print();
  xpetraParameters.print();

  /**********************************************************************************/
  /* CREATE INITAL MATRIX                                                           */
  /**********************************************************************************/
  const RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<const Operator>   Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Operator vs. CrsOperator

  // Using Galeri:
  //
  // RCP<Tpetra_Map> map = rcp( Galeri::CreateMap("Cartesian2D", comm, paramList) );
  // RCP<Tpetra_CrsMatrix> Op = rcp( Galeri::CreateCrsMatrix("Laplace2D", map.get(), paramList) );
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  if ( MueLu::MatrixVectorChecker<SC,LO,GO,NO>(Op) ) {
    std::cout << "OK !" << std::endl;
    return EXIT_SUCCESS;
  } else {
    std::cout << "FAILURE !" << std::endl;
    return EXIT_FAILURE;
  }

}

// JG TODO:
// - add a method CreateMap for xpetra/gallery (as Galeri)
// - wrap galeri matrix for the new MatrixVectorChecker
