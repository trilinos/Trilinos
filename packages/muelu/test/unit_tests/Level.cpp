#include "Teuchos_UnitTestHarness.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Level.hpp"

//TODO: some headers are missing
#include <Cthulhu_CrsOperator.hpp>

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(Level, SetCoreData)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

#ifdef HAVE_CTHULHU_TPETRA

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Operator> A = MueLu::UnitTest::create_1d_poisson_matrix<ScalarType, LO, GO>(10);

  Level firstLevel;
  firstLevel.SetA(A);
  TEUCHOS_TEST_EQUALITY(firstLevel.GetA(), A, out, success);
  firstLevel.SetR(A);
  TEUCHOS_TEST_EQUALITY(firstLevel.GetR(), A, out, success);
  firstLevel.SetP(A);
  TEUCHOS_TEST_EQUALITY(firstLevel.GetP(), A, out, success);
  firstLevel.SetLevelID(42);
  TEUCHOS_TEST_EQUALITY(firstLevel.GetLevelID(), 42, out, success);

/*
  RCP<Smoother> preSmoo = Smoother<ScalarType, LO, GO, Node, LMO>();
  TEUCHOS_TEST_EQUALITY(firstLevel.GetPreSmoother(), preSmoo, out, success);
  //RCP<Smoother> postSmoo = Smoother<ScalarType, LO, GO, Map, CrsOperator>();
*/


  //out << firstLevel << std::endl;
  /*
  out << "Testing copy ctor" << std::endl;
  Level secondLevel(firstLevel);
  //out << secondLevel << std::endl;
  */

#endif
}

}//namespace <anonymous>

