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
  } // Constructor test

  TEUCHOS_UNIT_TEST(RAPFactory, Correctness)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel, coarseLevel; TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    RCP<Operator> Op = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(27*comm->getSize());
    fineLevel.Set("A",Op);

    TentativePFactory tentpFactory;
    SaPFactory sapFactory(rcpFromRef(tentpFactory));
    TransPFactory transPFactory(rcpFromRef(sapFactory)); //todo:rcpFromRef

    coarseLevel.Request("P",&sapFactory);
    coarseLevel.Request("R",&transPFactory);

    coarseLevel.Request(sapFactory);
    coarseLevel.Request(transPFactory);
    sapFactory.Build(fineLevel,coarseLevel);
    transPFactory.Build(fineLevel,coarseLevel);

    RAPFactory rap(rcpFromRef(sapFactory), rcpFromRef(transPFactory));
    coarseLevel.Request(rap);

    coarseLevel.Request("A",&rap);
    rap.Build(fineLevel,coarseLevel);

    RCP<Operator> A = fineLevel.Get< RCP<Operator> >("A");
    RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", &sapFactory);
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

    RCP<Operator> coarseOp = coarseLevel.Get< RCP<Operator> >("A", &rap);

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

  } // Correctness test

  TEUCHOS_UNIT_TEST(RAPFactory, ImplicitTranspose)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    if (comm->getSize() > 1 && TestHelpers::Parameters::getLib() == Xpetra::UseEpetra ) {
      out << "Skipping ImplicitTranspose test for Epetra and #proc>1" << std::endl;
      return;
    }

    // build test-specific default factory manager
    RCP<MueLu::FactoryManagerBase> defManager = rcp(new MueLu::FactoryManagerBase());
    defManager->SetDefaultFactory("A", rcp(MueLu::NoFactory::get(),false));         // dummy factory for A
    defManager->SetDefaultFactory("Nullspace", rcp(new NullspaceFactory()));        // real null space factory for Ptent
    defManager->SetDefaultFactory("Graph", rcp(new CoalesceDropFactory()));         // real graph factory for Ptent
    defManager->SetDefaultFactory("Aggregates", rcp(new UCAggregationFactory()));   // real aggregation factory for Ptent

    Level fineLevel, coarseLevel;
    TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    // overwrite default factory manager
    fineLevel.SetDefaultFactoryHandler(defManager);
    coarseLevel.SetDefaultFactoryHandler(defManager);

    RCP<Operator> Op = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(19*comm->getSize());
    fineLevel.Set("A",Op);

    TentativePFactory tentpFactory;
    SaPFactory sapFactory(rcpFromRef(tentpFactory));
    TransPFactory transPFactory(rcpFromRef(sapFactory));

    coarseLevel.Request("P", &sapFactory);
    coarseLevel.Request("R", &transPFactory);

    coarseLevel.Request(sapFactory);
    coarseLevel.Request(transPFactory);
    sapFactory.Build(fineLevel, coarseLevel);
    transPFactory.Build(fineLevel,coarseLevel);
    RAPFactory rap(rcpFromRef(sapFactory), rcpFromRef(transPFactory));

    coarseLevel.Request("A", &rap);

    rap.SetImplicitTranspose(true);
    coarseLevel.Request(rap);
    rap.Build(fineLevel,coarseLevel);

    RCP<Operator> A = fineLevel.Get< RCP<Operator> >("A");
    RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", &sapFactory);
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

    RCP<Operator> coarseOp = coarseLevel.Get< RCP<Operator> >("A", &rap);

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

  } // Correctness test

} // namespace MueLuTests

