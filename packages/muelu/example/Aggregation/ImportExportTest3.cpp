#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include "Teuchos_VerboseObject.hpp"
#include <Teuchos_FancyOStream.hpp>

// Cthulhu
#include <Cthulhu_Parameters.hpp>
#include <Cthulhu_Map.hpp>
#include <Cthulhu_MapFactory.hpp>
#include <Cthulhu_CrsOperator.hpp>
#include <Cthulhu.hpp>

// Gallery
#define CTHULHU_ENABLED // == Gallery have to be build with the support of Cthulhu matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

#include <Cthulhu_VectorFactory.hpp>
#include <Cthulhu_ImportFactory.hpp>
#include <Cthulhu_MapFactory.hpp>
#include <MueLu_Graph.hpp>

#include "MueLu_UseShortNames.hpp"

using Teuchos::RCP;
using Teuchos::ArrayRCP;

int main(int argc, char *argv[]) {
  
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor cmdp(false);
  
  MueLu::Gallery::Parameters matrixParameters(cmdp); // manage parameters of the test case
  Cthulhu::Parameters cthulhuParameters(cmdp);       // manage parameters of cthulhu
  
  switch (cmdp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }
  
  matrixParameters.check();
  cthulhuParameters.check();

  if (comm->getRank() == 0) {
    matrixParameters.print();
    cthulhuParameters.print();
  }

  /**********************************************************************************/
  /* CREATE INITAL MATRIX                                                           */
  /**********************************************************************************/
  const RCP<const Map> map = MapFactory::Build(cthulhuParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<CrsOperator> Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Operator vs. CrsOperator
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  RCP<Graph> graph = rcp(new Graph(Op->getCrsGraph(), "Uncoupled"));
  
  // Mimic comm of ERR3
  // must be run with --nx=2
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/
  int MyPid = comm->getRank();
  RCP<Teuchos::FancyOStream> out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));
      
  const RCP<const Map> nonUniqueMap = graph->GetImportMap();

  RCP<LOVector> companion = LOVectorFactory::Build(nonUniqueMap);
  companion->putScalar(-1);

  int numMyWinners = 0;
  ArrayRCP<LO> myWinners(numMyWinners);

  Cthulhu::global_size_t g = -1; // TODO for Tpetra -1 == ??
  RCP<Map> winnerMap = MapFactory::Build(map->lib(), g, myWinners(), 0, comm);
  std::cout << MyPid << ": winnerMap(source)" << *winnerMap << std::endl;

  std::cout << MyPid << ": nonUniqueMap(dest)=" << *nonUniqueMap << std::endl;

  RCP<LOVector> justWinners = LOVectorFactory::Build(winnerMap);
  //justWinners->putScalar(-2);
  
  std::cout << MyPid << ": justWinners (input)=" << std::endl;
  justWinners->describe(*out, Teuchos::VERB_EXTREME);
  std::cout << std::endl << std::endl;

  const RCP<Import> pushWinners = ImportFactory::Build(winnerMap, nonUniqueMap); // error here





  companion->doImport(*justWinners, *pushWinners, Cthulhu::INSERT);

  std::cout << MyPid << ": nonUniqueMap=" << *nonUniqueMap << std::endl;

  std::cout << MyPid << ": companion (output)=" << std::endl;
  companion->describe(*out, Teuchos::VERB_EXTREME);
  std::cout << std::endl << std::endl;

  return EXIT_SUCCESS;
}
