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

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"
/**********************************************************************************/

int main(int argc, char *argv[]) {
  
  std::cout << "Hello World !" << std::endl;

  /**********************************************************************************/
  /* CREATE INITAL MATRIX                                                           */
  /**********************************************************************************/
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  LO numThreads=1;
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
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return EXIT_FAILURE;
  }

  std::cout << "#threads = " << numThreads << std::endl;
  std::cout << "problem size = " << nx*ny << std::endl;
  std::cout << "matrix type = " << matrixType << std::endl;

  Teuchos::ParameterList pl;
  pl.set("Num Threads",numThreads);

  GO numGlobalElements = nx*ny;
  if (nx*nx - (nx*ny/3)*3 != 0)
    throw(std::logic_error("problem size must be divisible by 3"));
  if (matrixType == "Laplace3D")
    numGlobalElements *= nz;
  LO indexBase = 0;

  RCP<const Map > map;
  map = rcp( new MyMap(numGlobalElements, indexBase, comm) );

  Teuchos::ParameterList matrixList;
  matrixList.set("nx",nx);
  matrixList.set("ny",ny);
  matrixList.set("nz",nz);

  RCP<CrsOperator> Op = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>(matrixType,map,matrixList); //TODO: Operator vs. CrsOperator

  RCP<const Epetra_CrsMatrix> A;

  { // Get the underlying Epetra Mtx (Wow ! It's paintful ! => I should create a function to do that)
    RCP<const CrsMatrix> tmp_CrsMtx = Op->getCrsMatrix();
    const RCP<const Cthulhu::EpetraCrsMatrix> &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Cthulhu::EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null) { std::cout << "Error !" << std::endl; return 1; }
    A = tmp_ECrsMtx->getEpetra_CrsMatrix();
  }

  //RCP<Vector> nullSpace = VectorFactory::Build(map);
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);

  MueLu::Hierarchy<SC,LO,GO,NO,LMO> H;
  H.setDefaultVerbLevel(Teuchos::VERB_HIGH);
  RCP<MueLu::Level<SC,LO,GO,NO,LMO> > Finest = rcp( new MueLu::Level<SC,LO,GO,NO,LMO>() );
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  Finest->SetA(Op);
  Finest->Save("Nullspace",nullSpace);

  Teuchos::RCP< Operator > cOp = Finest->GetA();
  //GO nFineDofs = cOp->getGlobalNumRows();
  //prolongator is nFineDofs by nCoarseDofs
  //Teuchos::RCP<Cthulhu::CrsOperator> Ptent = Teuchos::rcp( new Cthulhu::CrsOperator(cOp->Rowmap(), 2) );
  Teuchos::RCP< Operator > Ptent = Teuchos::rcp( new CrsOperator(cOp->getRowMap(), 2) );




  Finest->Save("NullSpace",nullSpace);
  H.SetLevel(Finest);

  RCP<SaPFactory>         Pfact = rcp( new SaPFactory() );
  RCP<TransPFactory>      Rfact = rcp( new TransPFactory() );
  RCP<RAPFactory>         Acfact = rcp( new RAPFactory() );
  RCP<SmootherFactory>    SmooFact = Teuchos::null;
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  //H.FillHierarchy(Pfact,Rfact,Acfact);
  H.FullPopulate(Pfact,Rfact,Acfact,SmooFact,0,maxLevels);
  
  return EXIT_SUCCESS;

}
