#include "Teuchos_UnitTestHarness.hpp"

#include "MueLu_Version.hpp"

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

}//namespace <anonymous>

