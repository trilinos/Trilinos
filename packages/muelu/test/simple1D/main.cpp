#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaLevel.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
//#include "MueLu_GaussSeidel.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_GenericPRFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"

/**********************************************************************************/
/* CREATE INITAL MATRIX                                                           */
/**********************************************************************************/
#include <Cthulhu_Map.hpp>
#include <Cthulhu_CrsMatrix.hpp>
#include <Cthulhu_EpetraCrsMatrix.hpp>
#include <Cthulhu_CrsOperator.hpp>
#include <Cthulhu.hpp>
#include <Cthulhu_Vector.hpp>
#include <Cthulhu_VectorFactory.hpp>
#include <Cthulhu_MultiVectorFactory.hpp>

#include <MueLu_MatrixFactory.hpp>

#include "MueLu_Utilities.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"
/**********************************************************************************/

int main(int argc, char *argv[]) {
  
#ifdef HAVE_MUELU_IFPACK //TODO

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  LO numThreads=1;
  LO its=10;
  GO nx=9;
  GO ny=9;
  GO nz=9;
  LO maxLevels = 3;
  Teuchos::CommandLineProcessor cmdp(false,true);
  std::string matrixType("Laplace1D");
  cmdp.setOption("nt",&numThreads,"number of threads.");
  cmdp.setOption("nx",&nx,"mesh points in x-direction.");
  cmdp.setOption("ny",&ny,"mesh points in y-direction.");
  cmdp.setOption("nz",&nz,"mesh points in z-direction.");
  cmdp.setOption("matrixType",&matrixType,"matrix type: Laplace1D, Laplace2D, Star2D, Laplace3D, Identity");
  cmdp.setOption("maxLevels",&maxLevels,"maximum number of levels allowed");
  cmdp.setOption("its",&its,"number of multigrid cycles");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return EXIT_FAILURE;
  }

  Teuchos::ParameterList pl;
  pl.set("Num Threads",numThreads);

  GO numGlobalElements = nx;
  if (matrixType == "Laplace2D" || matrixType == "Star2D")
    numGlobalElements *= ny;
  if (matrixType == "Laplace3D")
    numGlobalElements *= nz;
  if (numGlobalElements - (numGlobalElements/3)*3 != 0)
    throw(MueLu::Exceptions::RuntimeError("problem size must be divisible by 3"));
  LO indexBase = 0;

  std::cout << "#threads = " << numThreads << std::endl;
  std::cout << "problem size = " << numGlobalelements << std::endl;
  std::cout << "matrix type = " << matrixType << std::endl;

  RCP<const Map > map;
  map = rcp( new MyMap(numGlobalElements, indexBase, comm) );

  Teuchos::ParameterList matrixList;
  matrixList.set("nx",nx);
  matrixList.set("ny",ny);
  matrixList.set("nz",nz);

  RCP<CrsOperator> Op = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>(matrixType,map,matrixList); //TODO: Operator vs. CrsOperator

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  RCP<Epetra_MultiVector> foo = Utils::MV2NonConstEpetraMV(nullSpace);
  double n;
  foo->Norm1(&n);
  std::cout << "||NS|| = " << n << std::endl;

  MueLu::Hierarchy<SC,LO,GO,NO,LMO> H;
  H.setDefaultVerbLevel(Teuchos::VERB_HIGH);
  RCP<MueLu::Level<SC,LO,GO,NO,LMO> > Finest = rcp( new MueLu::Level<SC,LO,GO,NO,LMO>() );
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  Finest->SetA(Op);
  Finest->Save("Nullspace",nullSpace);
  Finest->Request("Nullspace"); //FIXME putting this in to avoid error until Merge needs business
                                          //FIXME is implemented

  Finest->Save("NullSpace",nullSpace);
  H.SetLevel(Finest);

  RCP<SaPFactory>         Pfact = rcp( new SaPFactory() );
  RCP<GenericPRFactory>   PRfact = rcp( new GenericPRFactory(Pfact));
  RCP<RAPFactory>         Acfact = rcp( new RAPFactory() );
  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smooProto = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );

  RCP<SmootherFactory>    SmooFact = rcp( new SmootherFactory(smooProto) );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  Teuchos::ParameterList status;
  status = H.FullPopulate(PRfact,Acfact,SmooFact,0,maxLevels);
  std::cout  << "======================\n Multigrid statistics \n======================" << std::endl;
  status.print(std::cout,Teuchos::ParameterList::PrintOptions().indent(2));

  //FIXME we should be able to just call smoother->SetNIts(50) ... but right now an exception gets thrown
  Teuchos::ParameterList amesosList;
  amesosList.set("PrintTiming",true);
  RCP<SmootherPrototype> coarseProto = rcp( new AmesosSmoother("Amesos_Klu",amesosList) );
  SmootherFactory coarseSolveFact(coarseProto);
  H.SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  RCP<Epetra_MultiVector> epX = Utils::MV2NonConstEpetraMV(X);
  epX->SetSeed(846930886);
  X->randomize();
  Op->multiply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  epX->Norm2(&n);
  std::cout << "||X_true|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << n << std::endl;

  X->putScalar( (SC) 0.0);

  H.PrintResidualHistory(true);
  H.Iterate(RHS,its,X);

  epX->Norm2(&n);
  std::cout << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << n << std::endl;


#endif // HAVE_MUELU_IFPACK
  
  return EXIT_SUCCESS;

}
