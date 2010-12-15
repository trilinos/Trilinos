#include "Teuchos_UnitTestHarness.hpp"
//#include "Teuchos_ParameterList.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(Hierarchy,Test0)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using namespace Teuchos;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<Level> levelTwo = rcp(new Level() );
  levelTwo->SetLevelID(2);

  Hierarchy H;

  H.SetLevel(levelOne);
  H.SetLevel(levelTwo);

  std::cout << H << std::endl;

  /* TODO
  Test set/get of R & P matrices.
  */
  /*
  RCP<Vector> x = rcp(new Cthulhu::TpetraVector<Scalar,LO,GO,Node>(map,nx) );
  RCP<Vector> y = rcp(new Cthulhu::TpetraVector<Scalar,LO,GO,Node>(map,nx) );
  x->putScalar(1);
  */

} //TEST0

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy1)
{

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP< SaPFactory >    PFact = Teuchos::null;
  TEST_THROW(H.FillHierarchy(PFact), std::logic_error);
}

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy2)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory>    PFact = rcp(new SaPFactory());

  out << "Providing just prolongator factory to FillHierarchy." << std::endl;
  H.FillHierarchy(PFact);
}

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy3)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory>    PFact = rcp(new SaPFactory());
  RCP<TransPFactory> RFact = rcp(new TransPFactory());
  RCP<RAPFactory>    AcFact= rcp(new RAPFactory());

  out << "Providing all three factories to FillHierarchy." << std::endl;
  H.FillHierarchy(PFact,RFact,AcFact);
}

TEUCHOS_UNIT_TEST(Hierarchy,SetSmoothers)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);

  Hierarchy H;
  H.SetLevel(levelOne);
  H.SetSmoothers();
}

}//namespace <anonymous>

