#include "Teuchos_UnitTestHarness.hpp"

#include "MueLu_Version.hpp"

#include "test_helpers.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_IfpackSmoother.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(SmootherFactory, DefaultCtor_Exception)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using namespace Teuchos;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<SmootherPrototype>  smoother = Teuchos::null;
  RCP<SmootherFactory> SmooFactory;


  TEST_THROW( SmooFactory = rcp(new SmootherFactory(smoother) ) , MueLu::Exceptions::RuntimeError );

}

TEUCHOS_UNIT_TEST(SmootherFactory, DefaultCtor_OneArg)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using namespace Teuchos;

  out << "version: " << MueLu::Version() << std::endl;


  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  RCP<SmootherFactory> SmooFactory = rcp(new SmootherFactory(smoother) );
  TEUCHOS_TEST_EQUALITY(SmooFactory != Teuchos::null, true, out, success);

}

TEUCHOS_UNIT_TEST(SmootherFactory, DefaultCtor_TwoArgs)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using namespace Teuchos;

  out << "version: " << MueLu::Version() << std::endl;


  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  RCP<SmootherFactory> SmooFactory = rcp(new SmootherFactory(smoother,smoother) );
  TEUCHOS_TEST_EQUALITY(SmooFactory != Teuchos::null, true, out, success);

}

TEUCHOS_UNIT_TEST(SmootherFactory, SetSmootherPrototypes)
{

  using namespace Teuchos;

  out << "version: " << MueLu::Version() << std::endl;

  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  RCP<SmootherFactory> SmooFactory = rcp(new SmootherFactory(smoother) );

  RCP<SmootherPrototype>  newSmoo1 = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  RCP<SmootherPrototype>  newSmoo2 = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  SmooFactory->SetSmootherPrototypes(newSmoo1,newSmoo2);
  RCP<SmootherPrototype>  checkSmoo1;
  RCP<SmootherPrototype>  checkSmoo2;
  SmooFactory->GetSmootherPrototypes(checkSmoo1,checkSmoo2);

  TEUCHOS_TEST_EQUALITY(checkSmoo1 == newSmoo1, true, out, success);
  TEUCHOS_TEST_EQUALITY(checkSmoo2 == newSmoo2, true, out, success);

}

TEUCHOS_UNIT_TEST(SmootherFactory, Build)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Testing SmootherFactory::Build method" << std::endl;

  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  RCP<SmootherFactory> smooFactory = rcp(new SmootherFactory(smoother) );

  RCP<Level> aLevel = rcp(new Level() );
  aLevel->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);

  RCP<SmootherPrototype>  preSmoo, postSmoo;
  //Check for exception if matrix is not set in Level.
  //FIXME once Level-GetA() doesn't throw an exception, this must be changed
  //TEST_THROW(smooFactory->Build(aLevel,preSmoo,postSmoo),MueLu::Exceptions::RuntimeError);
  TEST_THROW(smooFactory->Build(aLevel,preSmoo,postSmoo),std::logic_error);

  aLevel->SetA(A);
  smooFactory->Build(aLevel,preSmoo,postSmoo);
}

}//namespace <anonymous>

