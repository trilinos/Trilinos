/*
  This driver simply generates a Tpetra matrix, prints it to screen, and exits.
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

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "MueLu_MatrixFactory.hpp"
#include "MueLu_MatrixTypes.hpp"

#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaLevel.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"

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
  //typedef Tpetra::DefaultPlatform::DefaultPlatformType              Platform;
#ifdef HAVE_MPI
  typedef Tpetra::MpiPlatform<Kokkos::SerialNode>                   Platform;
#else
  typedef Tpetra::SerialPlatform<Kokkos::SerialNode>                   Platform;
#endif
  //Change this next typedef to get different sorts of Kokkos nodes.
  //typedef Platform::NodeType                                        NodeType;
  typedef Platform::NodeType                                        NO;

  typedef Ordinal  LO;
  typedef Ordinal  GO;
  typedef Scalar   SC;

  typedef Tpetra::Map<LO,GO,NO>              Map;
  typedef Tpetra::CrsMatrix<SC,LO,GO,NO>     CrsMatrix;
  typedef Tpetra::Vector<SC,LO,GO,NO>        Vector;

  typedef MueLu::SaPFactory<SC,LO,GO,NO>        SaPFactory;
  typedef MueLu::TransPFactory<SC,LO,GO,NO>     TransPFactory;
  typedef MueLu::RAPFactory<SC,LO,GO,NO>        RAPFactory;
  typedef MueLu::SmootherFactory<SC,LO,GO,NO>   SmootherFactory;

  using namespace Teuchos;

  oblackholestream blackhole;
  GlobalMPISession mpiSession(&argc,&argv,&blackhole);


  //use the "--help" option to get verbose help
  Ordinal numThreads=1;
  Ordinal nx=4;
  Ordinal ny=4;
  Ordinal nz=4;
  CommandLineProcessor cmdp(false,true);
  Ordinal maxLevels = 10;
  std::string matrixType("Laplace1D");
  cmdp.setOption("nt",&numThreads,"number of threads.");
  cmdp.setOption("nx",&nx,"mesh points in x-direction.");
  cmdp.setOption("ny",&ny,"mesh points in y-direction.");
  cmdp.setOption("nz",&nz,"mesh points in z-direction.");
  cmdp.setOption("matrixType",&matrixType,"matrix type: Laplace1D, Laplace2D, Star2D, Laplace3D");
  cmdp.setOption("maxLevels",&maxLevels,"maximum number of levels allowed");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  std::cout << "#threads = " << numThreads << std::endl;
  std::cout << "problem size = " << nx*ny << std::endl;
  std::cout << "matrix type = " << matrixType << std::endl;

  ParameterList pl;
  pl.set("Num Threads",numThreads);

  RCP<NO> node = rcp(new NO(pl));

  Platform myplat(node);
  RCP<const Comm<int> > comm = myplat.getComm();

  GO numGlobalElements = nx*ny;
  if (matrixType == "Laplace3D")
    numGlobalElements *= nz;
  GO indexBase = 0;

  RCP<const Map > map;
  map = rcp( new Map(numGlobalElements, indexBase, comm, Tpetra::GloballyDistributed, node) );

  RCP<CrsMatrix> A;
  ParameterList matrixList;
  matrixList.set("nx",nx);
  matrixList.set("ny",ny);
  matrixList.set("nz",nz);

  A = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Tpetra::Map<LO,GO>, Tpetra::CrsMatrix<Scalar,LO,GO> >(matrixType,map,matrixList);
  //RCP<FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  //A->describe(*out, VERB_EXTREME);

  RCP<Vector> nullSpace = rcp(new Vector(map) );
  nullSpace->putScalar( (SC) 1.0);

//Teuchos::VerboseObject<AlgorithmA>::setDefaultVerbLevel(Teuchos::VERB_NONE);

  MueLu::Hierarchy<SC,LO,GO,NO> H;
  H.setDefaultVerbLevel(Teuchos::VERB_NONE); //FIXME does nothing right now, must check inside Hierarchy
  RCP<MueLu::SaLevel<SC,LO,GO,NO> > Finest = rcp( new MueLu::SaLevel<SC,LO,GO,NO>() );
  Finest->setDefaultVerbLevel(Teuchos::VERB_NONE); //FIXME does nothing right now, must check inside SaLevel
  Finest->SetA(A);
  Finest->SetNullSpace(nullSpace);
  H.SetLevel(Finest);

  RCP<SaPFactory>         Pfact = rcp( new SaPFactory() );
  RCP<TransPFactory>      Rfact = rcp( new TransPFactory() );
  RCP<RAPFactory>         Acfact = rcp( new RAPFactory() );
  RCP<SmootherFactory>    SmooFact = Teuchos::null;

  //H.FillHierarchy(Pfact,Rfact,Acfact);
  H.FullPopulate(Pfact,Rfact,Acfact,SmooFact,0,maxLevels);

  return(0);
} //main
