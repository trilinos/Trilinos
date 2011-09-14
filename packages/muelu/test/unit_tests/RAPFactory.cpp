#include "Teuchos_UnitTestHarness.hpp"

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Utilities.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {
  
  TEUCHOS_UNIT_TEST(RAPFactory, Constructor)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<RAPFactory> rapFactory = rcp(new RAPFactory);
    TEST_EQUALITY(rapFactory != Teuchos::null, true);

    out << *rapFactory << std::endl;
  } //Constructor test

  TEUCHOS_UNIT_TEST(RAPFactory, Correctness)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel, coarseLevel; TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    RCP<Operator> Op = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(27*comm->getSize());
    fineLevel.Request("A",NULL);
    fineLevel.Set("A",Op,NULL);

    TentativePFactory tentpFactory;
    SaPFactory sapFactory(rcpFromRef(tentpFactory));
    TransPFactory transPFactory(rcpFromRef(sapFactory)); //todo:rcpFromRef

    coarseLevel.Request("P",&sapFactory);
    coarseLevel.Request("R",&transPFactory);
    coarseLevel.Request("A",NULL,false); // don't call DeclareInput for default factory of "A"

    RCP<GenericPRFactory> PRFac = rcp(new GenericPRFactory(rcpFromRef(sapFactory),rcpFromRef(transPFactory)));
    PRFac->SetMaxCoarseSize(1);
    PRFac->DeclareInput(fineLevel,coarseLevel);
    PRFac->Build(fineLevel,coarseLevel);
    RAPFactory rap(PRFac);
    rap.DeclareInput(fineLevel,coarseLevel);
    rap.Build(fineLevel,coarseLevel);

    RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", &sapFactory);
    RCP<Operator> A = fineLevel.Get< RCP<Operator> >("A",NULL);
    RCP<Operator> R = coarseLevel.Get< RCP<Operator> >("R", &transPFactory);

    RCP<MultiVector> workVec1 = MultiVectorFactory::Build(P->getRangeMap(),1);
    RCP<MultiVector> workVec2 = MultiVectorFactory::Build(Op->getRangeMap(),1);
    RCP<MultiVector> result1 = MultiVectorFactory::Build(R->getRangeMap(),1);
    RCP<MultiVector> X = MultiVectorFactory::Build(P->getDomainMap(),1);
    X->randomize();

    //Calculate result1 = R*(A*(P*X))
    P->apply(*X,*workVec1,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    Op->apply(*workVec1,*workVec2,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    R->apply(*workVec2,*result1,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

    RCP<Operator> coarseOp = coarseLevel.Get< RCP<Operator> >("A",NULL);

    //Calculate result2 = (R*A*P)*X
    RCP<MultiVector> result2 = MultiVectorFactory::Build(R->getRangeMap(),1);
    coarseOp->apply(*X,*result2,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
  
    Teuchos::Array<ST::magnitudeType> normX(1), normResult1(1),normResult2(1);
    X->norm2(normX);
    out << "This test checks the correctness of the Galerkin triple "
        << "matrix product by comparing (RAP)*X to R(A(P*X))." << std::endl;
    out << "||X||_2 = " << normX << std::endl; 
    result1->norm2(normResult1);
    result2->norm2(normResult2);
    TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 1e-12);

  } //Correctness test

  TEUCHOS_UNIT_TEST(RAPFactory, ImplicitTranspose)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel, coarseLevel;
    TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    RCP<Operator> Op = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(19*comm->getSize());
    fineLevel.Request("A",NULL);
    fineLevel.Set("A",Op,NULL);

    TentativePFactory tentpFactory;
    SaPFactory sapFactory(rcpFromRef(tentpFactory));
    TransPFactory transPFactory(rcpFromRef(sapFactory));

    coarseLevel.Request("P", &sapFactory);
    coarseLevel.Request("R", &transPFactory);
    coarseLevel.Request("A",NULL,false);
    fineLevel.Request("A",NULL);

    RCP<GenericPRFactory> PRFac = rcp(new GenericPRFactory(rcpFromRef(sapFactory),rcpFromRef(transPFactory)));
    PRFac->SetMaxCoarseSize(1);
    PRFac->DeclareInput(fineLevel,coarseLevel);
    PRFac->Build(fineLevel,coarseLevel);
    RAPFactory rap(PRFac);
    rap.SetImplicitTranspose(true);
    rap.DeclareInput(fineLevel,coarseLevel);
    rap.Build(fineLevel,coarseLevel);

//    sapFactory.DeclareInput(fineLevel,coarseLevel);
//    transPFactory.DeclareInput(fineLevel,coarseLevel);
//
//    sapFactory.BuildP(fineLevel,coarseLevel);
//    transPFactory.BuildR(fineLevel,coarseLevel);

    RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", &sapFactory);
    RCP<Operator> A = fineLevel.Get< RCP<Operator> >("A",NULL);
    RCP<Operator> R = coarseLevel.Get< RCP<Operator> >("R", &transPFactory);

    //std::string filename = "A.dat";
    //Utils::Write(filename,Op);
    //filename = "P.dat";
    //Utils::Write(filename,P);

    RCP<MultiVector> workVec1 = MultiVectorFactory::Build(P->getRangeMap(),1);
    RCP<MultiVector> workVec2 = MultiVectorFactory::Build(Op->getRangeMap(),1);
    RCP<MultiVector> result1 = MultiVectorFactory::Build(P->getDomainMap(),1);
    RCP<MultiVector> X = MultiVectorFactory::Build(P->getDomainMap(),1);
    X->randomize();
    //out.precision(12);
    //out.setOutputToRootOnly(-1);
    //X->describe(out,Teuchos::VERB_EXTREME);

    //Calculate result1 = P^T*(A*(P*X))
    P->apply(*X,*workVec1,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    Op->apply(*workVec1,*workVec2,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    P->apply(*workVec2,*result1,Teuchos::TRANS,(SC)1.0,(SC)0.0);



    RCP<Operator> coarseOp = coarseLevel.Get< RCP<Operator> >("A",NULL);

    //Calculate result2 = (R*A*P)*X
    RCP<MultiVector> result2 = MultiVectorFactory::Build(P->getDomainMap(),1);
    coarseOp->apply(*X,*result2,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
  
    Teuchos::Array<ST::magnitudeType> normX(1), normResult1(1),normResult2(1);
    X->norm2(normX);
    out << "This test checks the correctness of the Galerkin triple "
        << "matrix product by comparing (RAP)*X to R(A(P*X)), where R is the implicit tranpose of P." << std::endl;
    out << "||X||_2 = " << normX << std::endl; 
    result1->norm2(normResult1);
    result2->norm2(normResult2);
    TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 1e-12);
  } //Correctness test

}//namespace MueLuTests

