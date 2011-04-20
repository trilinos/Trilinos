#include "Teuchos_UnitTestHarness.hpp"
#include "Cthulhu_DefaultPlatform.hpp"
#include "MueLu_Version.hpp"
#include "test_helpers.hpp"

#include <Cthulhu_Map.hpp>
#include <Cthulhu_Operator.hpp>
#include <Cthulhu_MultiVectorFactory.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

bool testMpi = true;

Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
{
  Teuchos::RCP<const Teuchos::Comm<int> > ret;
  if (testMpi) {
    ret = Cthulhu::DefaultPlatform::getDefaultPlatform().getComm();
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
  TEST_THROW( smoother->Apply(*X,*RHS) , MueLu::Exceptions::RuntimeError );
  TEST_THROW( smoother->SetNIts(5), MueLu::Exceptions::RuntimeError );

}

TEUCHOS_UNIT_TEST(IfpackSmoother, GaussSeidel)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using Teuchos::RCP;
  using Teuchos::rcp;

  //FIXME this will probably fail in parallel b/c it becomes block Jacobi

  out << "version: " << MueLu::Version() << std::endl;
  out << "Tests interface to Ifpack's Gauss-Seidel preconditioner." << std::endl;

  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  RCP<Operator> Op = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO,NO,LMO>(125);

  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (int) 1);
  ifpackList.set("relaxation: damping factor", (double) 1.0);
  ifpackList.set("relaxation: zero starting solution", false);
  RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  Level aLevel;
  aLevel.SetA(Op);

  RCP<MultiVector> X = MultiVectorFactory::Build(Op->getDomainMap(),1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getDomainMap(),1);

  smoother->Setup(aLevel);

  RCP<Epetra_MultiVector> epX = Utils::MV2NonConstEpetraMV(X);
  epX->SetSeed(846930886);
  X->randomize();
  Op->multiply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  X->putScalar( (SC) 0.0);

  out << "Applying one GS sweep" << std::endl;
  smoother->Apply(*X,*RHS);
  Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms;
  norms = Utils::ResidualNorm(*Op,*X,*RHS);
  TEUCHOS_TEST_FLOATING_EQUALITY(norms[0],6.04555396884098,1e-12,out,success)

  int numIts = 50;
  smoother->SetNIts(numIts);
  TEUCHOS_TEST_EQUALITY(smoother->GetNIts(),50,out,success);
  out << "Applying " << numIts << " GS sweeps" << std::endl;
  X->putScalar( (SC) 0.0);
  smoother->Apply(*X,*RHS);
  norms = Utils::ResidualNorm(*Op,*X,*RHS);
  TEUCHOS_TEST_FLOATING_EQUALITY(norms[0],0.00912675857196253,1e-12,out,success)

} //GaussSeidel

TEUCHOS_UNIT_TEST(IfpackSmoother, Chebyshev)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  out << "Tests interface to Ifpack's Chebyshev preconditioner." << std::endl;

  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  RCP<Operator> Op = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO,NO,LMO>(125);

  Teuchos::ParameterList  ifpackList;
  ifpackList.set("chebyshev: degree", (int) 1);
  ifpackList.set("chebyshev: max eigenvalue", (double) 2.0);
  ifpackList.set("chebyshev: min eigenvalue", (double) 1.0);
  ifpackList.set("chebyshev: zero starting solution", false);
  RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("Chebyshev",ifpackList) );
  Level aLevel;
  aLevel.SetA(Op);

  RCP<MultiVector> X = MultiVectorFactory::Build(Op->getDomainMap(),1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getDomainMap(),1);
  RHS->putScalar( (SC) 0.0);

  smoother->Setup(aLevel);

  RCP<Epetra_MultiVector> epX = Utils::MV2NonConstEpetraMV(X);
  epX->SetSeed(846930886);
  X->randomize();
  double n;
  epX->Norm2(&n);
  X->scale(1/n);
  epX->Norm2(&n);
  out << "||X_initial|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << n << std::endl;
  //Op->multiply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);


  int numIts = 1;
  out << "Applying degree " << numIts << " Chebyshev smoother" << std::endl;
  smoother->SetNIts(numIts);
  smoother->Apply(*X,*RHS);
  epX->Norm2(&n);
  out << "||X_1|| = " << std::setiosflags(ios::fixed) << std::setprecision(25) << n << std::endl;
  TEUCHOS_TEST_EQUALITY(n<0.7,true,out,success);  //FIXME should calculate reduction analytically

  numIts = 10;
  smoother->SetNIts(numIts);
  TEUCHOS_TEST_EQUALITY(smoother->GetNIts(),10,out,success);
  out << "Applying degree " << numIts << " Chebyshev smoother" << std::endl;
  epX->SetSeed(846930886);
  X->randomize();
  epX->Norm2(&n);
  X->scale(1/n);
  epX->Norm2(&n);
  out << "||X_initial|| = " << std::setiosflags(ios::fixed) << std::setprecision(25) << n << std::endl;
  smoother->Apply(*X,*RHS);
  epX->Norm2(&n);
  out << "||X_" << std::setprecision(2) << numIts << "|| = " << std::setiosflags(ios::fixed) <<
std::setprecision(20) << n << std::endl;
  TEUCHOS_TEST_EQUALITY(n<0.25,true,out,success);  //FIXME should calculate reduction analytically

} //Chebyshev

TEUCHOS_UNIT_TEST(IfpackSmoother, ILU)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  //FIXME this will probably fail in parallel b/c it becomes block Jacobi

  out << "version: " << MueLu::Version() << std::endl;

  out << "Tests interface to Ifpack's ILU(0) preconditioner." << std::endl;

  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  RCP<Operator> Op = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO,NO,LMO>(125);

  Teuchos::ParameterList  ifpackList;
  RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("ILU",ifpackList) );
  Level aLevel;
  aLevel.SetA(Op);

  RCP<MultiVector> Xtrue = MultiVectorFactory::Build(Op->getDomainMap(),1);
  RCP<MultiVector> X = MultiVectorFactory::Build(Op->getDomainMap(),1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getDomainMap(),1);

  smoother->Setup(aLevel);

  RCP<Epetra_MultiVector> epXtrue = Utils::MV2NonConstEpetraMV(Xtrue);
  epXtrue->SetSeed(846930886);
  Xtrue->randomize();
  double n;
  epXtrue->Norm2(&n);
  Xtrue->scale(1/n);
  Op->multiply(*Xtrue,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
  RCP<Epetra_MultiVector> epRHS = Utils::MV2NonConstEpetraMV(RHS);
  epRHS->Norm2(&n);
  out << "||RHS|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << n << std::endl;
  X->putScalar( (SC) 0.0);


  smoother->Apply(*X,*RHS);
  RCP<Epetra_MultiVector> epX = Utils::MV2NonConstEpetraMV(X);
  epX->Norm2(&n);
  out << "||X_final|| = " << std::setiosflags(ios::fixed) << std::setprecision(25) << n << std::endl;
  Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms;
  norms = Utils::ResidualNorm(*Op,*X,*RHS);
  out << "||residual|| = " << norms[0] << std::endl;
  TEUCHOS_TEST_EQUALITY(norms[0]<1e-10,true,out,success)

} //ILU

}//namespace <anonymous>

