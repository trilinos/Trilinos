#include <iostream>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>

#if !defined(CTHULHU_USE_TPETRA) && !defined(CTHULHU_USE_EPETRA)
#define CTHULHU_USE_TPETRA
#endif
#include <Cthulhu_Map.hpp>
#include <Cthulhu_CrsMatrix.hpp>
#include <Cthulhu_Operator.hpp>
#include <Cthulhu.hpp>

#define CTHULHU_ENABLED //TODO
#include <MueLu_MatrixFactory.hpp>

using Teuchos::RCP;

/*
  This driver simply generates a Cthulhu matrix, prints it to screen, and exits.

  Use the "--help" option to get verbose help.
*/

int main(int argc, char* argv[]) 
{

#ifdef CTHULHU_USE_TPETRA
  std::cout << "CTHULHU_USE_TPETRA" << std::endl;
#else
  std::cout << "CTHULHU_USE_EPETRA" << std::endl;
#endif

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  LO numThreads=1;
  GO nx=4;
  GO ny=4;
  GO nz=4;
  Teuchos::CommandLineProcessor cmdp(false,true);
  std::string matrixType("Laplace1D");
  cmdp.setOption("nt",&numThreads,"number of threads.");
  cmdp.setOption("nx",&nx,"mesh points in x-direction.");
  cmdp.setOption("ny",&ny,"mesh points in y-direction.");
  cmdp.setOption("nz",&nz,"mesh points in z-direction.");
  cmdp.setOption("matrixType",&matrixType,"matrix type: Laplace1D, Laplace2D, Star2D, Laplace3D");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return EXIT_FAILURE;
  }

  std::cout << "#threads = " << numThreads << std::endl;
  std::cout << "problem size = " << nx*ny << std::endl;
  std::cout << "matrix type = " << matrixType << std::endl;

  Teuchos::ParameterList pl;
  pl.set("Num Threads",numThreads);

  GO numGlobalElements = nx*ny;
  if (matrixType == "Laplace3D")
    numGlobalElements *= nz;
  LO indexBase = 0;

  RCP<const Map > map;
  map = rcp( new MyMap(numGlobalElements, indexBase, comm) );

  Teuchos::ParameterList matrixList;
  matrixList.set("nx",nx);
  matrixList.set("ny",ny);
  matrixList.set("nz",nz);

  {
    // RCP<CrsMatrix> A = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsMatrix>(matrixType,map,matrixList);
    RCP<Operator> A = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, Operator>(matrixType,map,matrixList);
    
    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    if (comm->getRank() == 0)
      std::cout << "\n================ MAP =====================================================\n" << std::endl;
    map->describe(*out, Teuchos::VERB_EXTREME);
    comm->barrier();
    sleep(1);
    
    if (comm->getRank() == 0)
      std::cout << "\n================ MATRIX ==================================================\n" << std::endl;
    A->describe(*out, Teuchos::VERB_EXTREME);
  }
  
  return EXIT_SUCCESS;
} 
