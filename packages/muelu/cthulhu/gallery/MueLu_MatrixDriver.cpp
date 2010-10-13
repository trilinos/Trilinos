/*
  This driver simply generates a Cthulhu matrix, prints it to screen, and exits.
*/
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_FancyOStream.hpp"

#include "Cthulhu_DefaultPlatform.hpp"
#include "Cthulhu_Map.hpp"
#include "Cthulhu_CrsMatrix.hpp"

#include "Cthulhu_TpetraMap.hpp"

#include "MueLu_MatrixFactory.hpp"
#include "MueLu_MatrixTypes.hpp"

#include <iostream>

#include <Kokkos_SerialNode.hpp>
#ifdef KOKKOS_HAVE_TBB
#include <Kokkos_TBBNode.hpp>
#endif
#ifdef KOKKOS_HAVE_THREADPOOL
#include <Kokkos_TPINode.hpp>
#endif

//** ***************************************************************************
//**                               main program
//** ***************************************************************************

int main(int argc, char** argv) 
{
  typedef int                                                       Ordinal;
  typedef double                                                    Scalar;
  //typedef Cthulhu::DefaultPlatform::DefaultPlatformType              Platform;
#ifdef HAVE_MPI
  typedef Cthulhu::MpiPlatform<Kokkos::SerialNode>                   Platform;
#else
  typedef Cthulhu::SerialPlatform<Kokkos::SerialNode>                   Platform;
#endif
  //Change this next typedef to get different sorts of Kokkos nodes.
  typedef Platform::NodeType                                        NodeType;


  using namespace Teuchos;

  oblackholestream blackhole;
  GlobalMPISession mpiSession(&argc,&argv,&blackhole);


  //use the "--help" option to get verbose help
  Ordinal numThreads=1;
  Ordinal nx=4;
  Ordinal ny=4;
  Ordinal nz=4;
  CommandLineProcessor cmdp(false,true);
  std::string matrixType("Laplace1D");
  cmdp.setOption("nt",&numThreads,"number of threads.");
  cmdp.setOption("nx",&nx,"mesh points in x-direction.");
  cmdp.setOption("ny",&ny,"mesh points in y-direction.");
  cmdp.setOption("nz",&nz,"mesh points in z-direction.");
  cmdp.setOption("matrixType",&matrixType,"matrix type: Laplace1D, Laplace2D, Star2D, Laplace3D");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  std::cout << "#threads = " << numThreads << std::endl;
  std::cout << "problem size = " << nx*ny << std::endl;
  std::cout << "matrix type = " << matrixType << std::endl;

  ParameterList pl;
  pl.set("Num Threads",numThreads);

  RCP<NodeType> node = rcp(new NodeType(pl));

  Platform myplat(node);
  RCP<const Comm<int> > comm = myplat.getComm();

  Ordinal numGlobalElements = nx*ny;
  if (matrixType == "Laplace3D")
    numGlobalElements *= nz;
  Ordinal indexBase = 0;

  RCP<const Cthulhu::Map<Ordinal,Ordinal,NodeType> > map;
  map = rcp( new Cthulhu::TpetraMap<Ordinal,Ordinal,NodeType>(numGlobalElements, indexBase, comm, Tpetra::GloballyDistributed, node) );

  RCP<Cthulhu::CrsMatrix<Scalar,Ordinal,Ordinal,NodeType> > A;
  ParameterList matrixList;
  matrixList.set("nx",nx);
  matrixList.set("ny",ny);
  matrixList.set("nz",nz);
  RCP<FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  if (comm->getRank() == 0)
    std::cout << "\n================ MAP =====================================================\n" << std::endl;
  map->describe(*out, VERB_EXTREME);
  comm->barrier();
  sleep(1);

  if (comm->getRank() == 0)
    std::cout << "\n================ MATRIX ==================================================\n" << std::endl;
  A = CreateCrsMatrix<Scalar,Ordinal,Ordinal,NodeType>(matrixType,map,matrixList);
  A->describe(*out, VERB_EXTREME);

  return(0);
} //main
