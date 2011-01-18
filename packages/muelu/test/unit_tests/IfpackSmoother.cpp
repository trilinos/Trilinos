#include "Teuchos_UnitTestHarness.hpp"

#include "Tpetra_DefaultPlatform.hpp"

#include "MueLu_Version.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Utilities.hpp"

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

namespace {

bool testMpi = true;

Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
{
  Teuchos::RCP<const Teuchos::Comm<int> > ret;
  if (testMpi) {
    ret = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  }
  else {
    ret = Teuchos::rcp(new Teuchos::SerialComm<int>());
  }
  return ret;
}

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(IfpackSmoother, NotSetup)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  GO nx,ny,nz;
  nx = ny = nz = 5;
  GO numGlobalElements = nx*ny*nz;
  LO indexBase = 0;
  RCP<const Map > map;
  map = rcp( new Cthulhu::EpetraMap(numGlobalElements, indexBase, comm) );

  Teuchos::ParameterList  ifpackList;
  RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );

  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  //try applying without setting up
  TEST_THROW( smoother->Apply(X,RHS) , MueLu::Exceptions::RuntimeError );
  TEST_THROW( smoother->SetNIts(5), MueLu::Exceptions::RuntimeError );

}

TEUCHOS_UNIT_TEST(IfpackSmoother, GaussSeidelApply)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  GO nx,ny,nz;
  nx = ny = nz = 5;
  GO numGlobalElements = nx*ny*nz;
  LO indexBase = 0;
  RCP<const Map > map;
  map = rcp( new Cthulhu::EpetraMap(numGlobalElements, indexBase, comm) );

  Teuchos::ParameterList matrixList;
  matrixList.set("nx",nx);
  matrixList.set("ny",ny);
  matrixList.set("nz",nz);

  std::string matrixType("Laplace1D");
  RCP<CrsOperator> Op = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>(matrixType,map,matrixList);

  //Epetra_Util blech;

  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (int) 1);
  ifpackList.set("relaxation: damping factor", (double) 1.0);
  ifpackList.set("relaxation: zero starting solution", false);
  RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  //RCP<Level> aLevel = rcp( new Level() );
  //aLevel->SetA(Op);
  Level aLevel;
  aLevel.SetA(Op);

  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  smoother->Setup(aLevel);

  RCP<Epetra_MultiVector> epX = Utils::MV2NonConstEpetraMV(X);
  epX->SetSeed(846930886);
  X->randomize();
  Op->multiply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  //double n;
  //epX->Norm2(&n);
  //std::cout << "||X_true|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << n << std::endl;

  X->putScalar( (SC) 0.0);

  out << "Applying one GS sweep" << std::endl;
  smoother->Apply(X,RHS);
  Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms;
  norms = Utils::ResidualNorm(*Op,*X,*RHS);
  TEUCHOS_TEST_FLOATING_EQUALITY(norms[0],6.04555396884098,1e-12,out,success)

  int numIts = 50;
  smoother->SetNIts(numIts);
  TEUCHOS_TEST_EQUALITY(smoother->GetNIts(),50,out,success);
  out << "Applying " << numIts << " GS sweeps" << std::endl;
  X->putScalar( (SC) 0.0);
  smoother->Apply(X,RHS);
  norms = Utils::ResidualNorm(*Op,*X,*RHS);
  TEUCHOS_TEST_FLOATING_EQUALITY(norms[0],0.00912675857196253,1e-12,out,success)

}

}//namespace <anonymous>

