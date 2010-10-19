#include <iostream>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>

#include "Cthulhu_DefaultPlatform.hpp"
#include "Cthulhu_Map.hpp"
#include "Cthulhu_CrsMatrix.hpp"

#include "Cthulhu_TpetraMap.hpp"

#include "MueLu_MatrixFactory.hpp"
#include "MueLu_MatrixTypes.hpp"

using Teuchos::RCP;

/*
  This driver simply generates a Cthulhu matrix, prints it to screen, and exits.

  Use the "--help" option to get verbose help.
*/

int main(int argc, char* argv[]) 
{
  typedef int                                                       Ordinal;
  typedef double                                                    Scalar;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  Ordinal numThreads=1;
  Ordinal nx=4;
  Ordinal ny=4;
  Ordinal nz=4;
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

  Ordinal numGlobalElements = nx*ny;
  if (matrixType == "Laplace3D")
    numGlobalElements *= nz;
  Ordinal indexBase = 0;

  RCP<const Cthulhu::Map<Ordinal,Ordinal> > map;
  map = rcp( new Cthulhu::TpetraMap<Ordinal,Ordinal>(numGlobalElements, indexBase, comm) );

  RCP<Cthulhu::CrsMatrix<Scalar,Ordinal,Ordinal> > A;
  Teuchos::ParameterList matrixList;
  matrixList.set("nx",nx);
  matrixList.set("ny",ny);
  matrixList.set("nz",nz);
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  if (comm->getRank() == 0)
    std::cout << "\n================ MAP =====================================================\n" << std::endl;
  map->describe(*out, Teuchos::VERB_EXTREME);
  comm->barrier();
  sleep(1);

  if (comm->getRank() == 0)
    std::cout << "\n================ MATRIX ==================================================\n" << std::endl;
  A = CreateCrsMatrix<Scalar,Ordinal,Ordinal>(matrixType,map,matrixList);
  A->describe(*out, Teuchos::VERB_EXTREME);

  return EXIT_SUCCESS;
} // main
