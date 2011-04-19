#include "Teuchos_UnitTestHarness.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_TransPFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"  
#include "MueLu_UseShortNames.hpp"  

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(TransPFactory, Test0)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  out << "version: " << MueLu::Version() << std::endl;

  RCP<TransPFactory> transPFact = rcp(new TransPFactory);
  TEUCHOS_TEST_EQUALITY(transPFact != Teuchos::null, true, out, success);

  out << *transPFact << std::endl;

}

}//namespace <anonymous>

