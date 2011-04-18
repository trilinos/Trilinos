#include "Teuchos_UnitTestHarness.hpp"

#include "Cthulhu_DefaultPlatform.hpp"

#define HAVE_CTHULHU_EPETRA
#define CTHULHU_USE_EPETRA
#include "MueLu_config.hpp"

//#include "Cthulhu_config.hpp"

#include <Cthulhu_Map.hpp>
#include <Cthulhu_Example.hpp>
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"

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
TEUCHOS_UNIT_TEST(RAPFactory, Constructor)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<RAPFactory> rapFactory = rcp(new RAPFactory);
  TEUCHOS_TEST_EQUALITY(rapFactory != Teuchos::null, true, out, success);

  out << *rapFactory << std::endl;
} //Constructor test

TEUCHOS_UNIT_TEST(RAPFactory, Correctness)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  RCP<CrsOperator> Op = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(27);
  MueLu::Level<SC,LO,GO,NO,LMO> fineLevel;
  fineLevel.SetA(Op);
  MueLu::Level<SC,LO,GO,NO,LMO> coarseLevel;

  MueLu::SaPFactory<SC,LO,GO,NO,LMO> sapFactory;
  sapFactory.BuildP(fineLevel,coarseLevel);

  RCP<Operator> P = coarseLevel.GetP();
  RCP<Operator> A = fineLevel.GetA();

  MueLu::TransPFactory<SC,LO,GO,NO,LMO> transPFactory;
  transPFactory.BuildR(fineLevel,coarseLevel);

  RCP<Operator> R = coarseLevel.GetR();

  RCP<MultiVector> workVec1 = MultiVectorFactory::Build(P->getRangeMap(),1);
  RCP<MultiVector> workVec2 = MultiVectorFactory::Build(Op->getRangeMap(),1);
  RCP<MultiVector> result1 = MultiVectorFactory::Build(R->getRangeMap(),1);
  RCP<MultiVector> X = MultiVectorFactory::Build(P->getDomainMap(),1);
  X->randomize();

  //Calculate result1 = R*(A*(P*X))
  P->multiply(*X,*workVec1,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
  Op->multiply(*workVec1,*workVec2,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
  R->multiply(*workVec2,*result1,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  MueLu::RAPFactory<SC,LO,GO,NO,LMO> rap;
  rap.Build(fineLevel,coarseLevel);

  RCP<Operator> coarseOp = coarseLevel.GetA();

  //Calculate result2 = (R*A*P)*X
  RCP<MultiVector> result2 = MultiVectorFactory::Build(Op->getRangeMap(),1);
  coarseOp->multiply(*X,*result2,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
  
  Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> normX(1), normResult1(1),normResult2(1);
  X->norm2(normX);
  out << "This test checks the correctness of the Galerkin triple "
      << "matrix product by comparing (RAP)*X to R(A(P*X))." << std::endl;
  out << "||X||_2 = " << normX << std::endl; 
  result1->norm2(normResult1);
  result2->norm2(normResult2);
  TEUCHOS_TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 1e-12, out, success);
} //Correctness test

}//namespace <anonymous>

