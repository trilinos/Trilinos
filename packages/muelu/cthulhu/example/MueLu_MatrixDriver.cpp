#include <iostream>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Cthulhu_Map.hpp>
#include <Cthulhu_MapFactory.hpp>
#include <Cthulhu_CrsMatrix.hpp>
#include <Cthulhu_Operator.hpp>
#include <Cthulhu_Example.hpp>
#include <Cthulhu_Parameters.hpp>

// Gallery
#define CTHULHU_ENABLED // == Gallery have to be build with the support of Cthulhu matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

//#define CTHULHU_MAPFACTORY_SHORT
#include <Cthulhu_UseShortNamesOrdinal.hpp>



/*
  This driver simply generates a Cthulhu matrix, prints it to screen, and exits.

  Use the "--help" option to get verbose help.
*/

int main(int argc, char* argv[]) 
{
  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);
  
  MueLu::Gallery::Parameters<GO> matrixParameters(clp); // manage parameters of the test case
  Cthulhu::Parameters cthulhuParameters(clp);       // manage parameters of cthulhu
  
  switch (clp.parse(argc,argv)) {
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

//   std::cout << "#threads = " << numThreads << std::endl;
//   std::cout << "problem size = " << nx*ny << std::endl;
//   std::cout << "matrix type = " << matrixType << std::endl;

//   Teuchos::ParameterList pl;
//   pl.set("Num Threads",numThreads);

  //typedef Cthulhu::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactory;

  const RCP<const Map> map = MapFactory::Build(cthulhuParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);

  {
    RCP<Operator> A = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, Operator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
    
    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    if (comm->getRank() == 0)
      std::cout << "\n================ MAP =====================================================\n" << std::endl;
    map->describe(*out, Teuchos::VERB_EXTREME);
    comm->barrier();
    //    sleep(1);
    
    if (comm->getRank() == 0)
      std::cout << "\n================ MATRIX ==================================================\n" << std::endl;
    A->describe(*out, Teuchos::VERB_EXTREME);
  }
  
  return EXIT_SUCCESS;
} 
