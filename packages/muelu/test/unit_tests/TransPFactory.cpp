#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_TransPFactory.hpp"
#include "MueLu_SaPFactory.hpp"  

#include "MueLu_UseDefaultTypes.hpp"  
#include "MueLu_UseShortNames.hpp"  

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(TransPFactory, Constructor)
  {

    out << "version: " << MueLu::Version() << std::endl;

    RCP<TransPFactory> transPFact = rcp(new TransPFactory);
    TEST_EQUALITY(transPFact != Teuchos::null, true);

    out << *transPFact << std::endl;

  }

  TEUCHOS_UNIT_TEST(TransPFactory, Correctness)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();

    Level fineLevel, coarseLevel;
    fineLevel.SetupPhase(true); coarseLevel.SetupPhase(true);
    MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    // Test of createTwoLevelHierarchy: to be moved...
    TEST_EQUALITY(fineLevel.GetLevelID(), 1);
    TEST_EQUALITY(coarseLevel.GetLevelID(), 2);
    //compilation warning TEST_EQUALITY(fineLevel.GetPreviousLevel().get(), NULL);
    //TEST_EQUALITY(coarseLevel.GetPreviousLevel().get(), &fineLevel);
    // --

    RCP<Operator> Op = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(27*comm->getSize());
    fineLevel.Set("A",Op,NULL);

    SaPFactory sapFactory;
    sapFactory.BuildP(fineLevel,coarseLevel);

    TransPFactory transPFact(rcpFromRef(sapFactory)); //todo:rcpFromRef
    transPFact.BuildR(fineLevel,coarseLevel);

    RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", &sapFactory);
    RCP<Operator> R = coarseLevel.Get< RCP<Operator> >("R", &transPFact);

    RCP<MultiVector> result1 = MultiVectorFactory::Build(P->getDomainMap(),1);
    RCP<MultiVector> result2  = MultiVectorFactory::Build(R->getRangeMap(),1);
    RCP<MultiVector> X = MultiVectorFactory::Build(P->getRangeMap(),1);
    X->randomize();

    //Calculate P^T * X
    P->apply(*X,*result1,Teuchos::TRANS,(SC)1.0,(SC)0.0);
    //Calculate R * X
    R->apply(*X,*result2,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

    Teuchos::Array<ST::magnitudeType> normX(1), normResult1(1),normResult2(1);
    X->norm2(normX);
    out << "This test checks the correctness of the R created by TransPFactory." << std::endl;
    out << "||X||_2 = " << normX << std::endl; 
    result1->norm2(normResult1);
    result2->norm2(normResult2);
    TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 1e-12);

  } //Correctness test

}//namespace MueLuTests
