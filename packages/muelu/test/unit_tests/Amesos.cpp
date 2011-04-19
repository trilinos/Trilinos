#include "Teuchos_UnitTestHarness.hpp"
#include "Cthulhu_DefaultPlatform.hpp"
#include "MueLu_Version.hpp"
#include "test_helpers.hpp"

#include <Cthulhu_Map.hpp>
#include <Cthulhu_CrsOperator.hpp>
#include <Cthulhu_MultiVectorFactory.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"
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
TEUCHOS_UNIT_TEST(Amesos, NotSetup)
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

  Teuchos::ParameterList  amesosList;
  RCP<AmesosSmoother>  smoother = rcp( new AmesosSmoother("Amesos_Klu",amesosList) );

  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  //try applying without setting up
  TEST_THROW( smoother->Apply(*X,*RHS) , MueLu::Exceptions::RuntimeError );
  TEST_THROW( smoother->SetNIts(5), MueLu::Exceptions::RuntimeError );

}

TEUCHOS_UNIT_TEST(Amesos, KLUSolve)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  RCP<CrsOperator> Op = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(125);

  Teuchos::ParameterList  amesosList;
  amesosList.set("PrintTiming",false);
  amesosList.set("OutputLevel",0);
  RCP<AmesosSmoother>  smoother = rcp( new AmesosSmoother("Amesos_Klu",amesosList) );
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
  Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> res;
  res = Utils::ResidualNorm(*Op,*X,*RHS);
  out << "||initial residual|| = " << res[0] << std::endl;

  smoother->Apply(*X,*RHS);
  res = Utils::ResidualNorm(*Op,*X,*RHS);
  out << "||final residual|| = " << res[0] << std::endl;
  TEUCHOS_TEST_EQUALITY(res[0] < 1e-12,true,out,success)
}

}//namespace <anonymous>

