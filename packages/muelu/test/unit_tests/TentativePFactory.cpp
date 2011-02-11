#include "Teuchos_UnitTestHarness.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(TentativePFactory, Constructor)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<TentativePFactory> tentPFact = rcp(new TentativePFactory);
  TEUCHOS_TEST_EQUALITY(tentPFact != Teuchos::null, true, out, success);

} //Constructor

TEUCHOS_UNIT_TEST(TentativePFactory, SetGetMethods)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  TentativePFactory tentPFact;

  bool flag = tentPFact.TentativeWithQR();
  TEUCHOS_TEST_EQUALITY(flag, false, out, success);
  tentPFact.TentativeWithQR(true);
  flag = tentPFact.TentativeWithQR();
  TEUCHOS_TEST_EQUALITY(flag, true, out, success);
} //SetGetMethods

//TODO test BuildP
//TODO test MakeTentative

TEUCHOS_UNIT_TEST(TentativePFactory, MakeTentativeWithQR)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel;
  fineLevel.SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(12);
  fineLevel.SetA(A);

  Level coarseLevel;

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(),1);
  nullSpace->putScalar( (SC) 1.0);
  fineLevel.Save("Nullspace",nullSpace);
  fineLevel.Request("Nullspace"); //FIXME putting this in to avoid error until Merge needs business
                                  //FIXME is implemented

  TentativePFactory::MakeTentativeWithQR(fineLevel,coarseLevel);

  RCP<MultiVector> coarseNullSpace; 
  coarseLevel.Examine("Nullspace",coarseNullSpace);
  out << *coarseNullSpace << std::endl;


}


}//namespace <anonymous>

