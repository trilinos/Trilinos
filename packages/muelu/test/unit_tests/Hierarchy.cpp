#include "Teuchos_UnitTestHarness.hpp"
//#include "Teuchos_ParameterList.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_PRFactory.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_RFactory.hpp"
#include "MueLu_GenericPRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

#ifdef SKIP_THIS_FOR_NOW
//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(Hierarchy,Test0)
{

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
#endif
TEUCHOS_UNIT_TEST(Hierarchy,Constructor)
{

  using namespace Teuchos;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Hierarchy> H = rcp(new Hierarchy);

  TEUCHOS_TEST_INEQUALITY(H, Teuchos::null, out, success);

} //Constructor

TEUCHOS_UNIT_TEST(Hierarchy,SetAndGetLevel)
{

  using namespace Teuchos;

  out << "version: " << MueLu::Version() << std::endl;

  Hierarchy H;
  RCP<Level> level = rcp(new Level() );
  H.SetLevel(level);
  RCP<Level> dupLevel = H.GetLevel(0);
  TEUCHOS_TEST_EQUALITY(level, dupLevel, out, success);

}//SetAndGetLevel

TEUCHOS_UNIT_TEST(Hierarchy,NumberOfLevels)
{

  using namespace Teuchos;

  out << "version: " << MueLu::Version() << std::endl;

  Hierarchy H;
  RCP<Level> levelOne = rcp(new Level() );
  RCP<Level> levelTwo = rcp(new Level() );
  RCP<Level> levelThree = rcp(new Level() );
  H.SetLevel(levelOne);
  H.SetLevel(levelTwo);
  H.SetLevel(levelThree);
  TEUCHOS_TEST_EQUALITY(H.GetNumberOfLevels(), 3, out, success);

}//NumberOfLevels

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_PFact_Only)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  /*
  TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  SaPFactory must inherit from PFactory ... dude
  TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  */

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  GenericPRFactory PRFact(PFact);

  out << "Providing just prolongator factory to FillHierarchy." << std::endl;
  //FIXME until Hierarchy.FillHierarchy takes a generic (?) RP factory, this
  //FIXME won't work with with more than 2 levels
  H.FillHierarchy(PRFact);
}

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_NoFactoriesGiven)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  out << "Providing no factories to FillHierarchy." << std::endl;
  //FIXME until Hierarchy.FillHierarchy takes a generic (?) RP factory, this
  //FIXME test will ALWAYS fail because it is trying to produce more than two levels.
  //FIXME Without an R factory, Acoarse can't be made, and the next time through,
  //FIXME GetA fails.
  H.FillHierarchy();
}

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy1)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP< PRFactory >    PFact = Teuchos::null;
  //FIXME this won't throw anymore ...
  //TEST_THROW(H.FillHierarchy(PFact), std::logic_error); FIXME -- do we need this test anymore?
}

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_PRFactoryOnly)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  GenericPRFactory PRFact(PFact);

  out << "Providing just PR factory to FillHierarchy." << std::endl;
  H.FillHierarchy(PRFact);
}

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_BothFactories)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  GenericPRFactory PRFact(PFact);
  RAPFactory    AcFact;

  out << "Providing both factories to FillHierarchy." << std::endl;
  H.FillHierarchy(PRFact,AcFact);
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

