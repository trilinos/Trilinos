#include "Teuchos_UnitTestHarness.hpp"

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Utilities.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_MultiVectorTransferFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_SmootherFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"

#define XPETRA_ENABLED // == Gallery have to be build with the support of Xpetra matrices.
#include "MueLu_GalleryUtils.hpp"

#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {
  
  TEUCHOS_UNIT_TEST(MultiVectorTransferFactory, Constructor)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<TentativePFactory>    PtentFact = rcp(new TentativePFactory());
    RCP<MueLu::MultiVectorTransferFactory<SC, LO, GO, NO, LMO> > mvtf = rcp(new MueLu::MultiVectorTransferFactory<SC, LO, GO, NO, LMO>("Coordinates","P",PtentFact));
    TEST_EQUALITY(mvtf != Teuchos::null, true);
  } // Constructor test

  //------------------------------------------------------------------------------------------
  
  TEUCHOS_UNIT_TEST(MultiVectorTransferFactory, Build)
  {
    out << "version: " << MueLu::Version() << std::endl;

    out << "Tests the action of the transfer factory on a vector.  In this test, the transfer is the tentative" << std::endl;
    out << "prolongator, and the vector is all ones.  So the norm of the resulting coarse grid vector should be" << std::endl;
    out << "equal to the number of fine degrees of freedom." << std::endl;

    Level fineLevel, coarseLevel;
    TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
    GO nx = 199;
    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(nx);
    fineLevel.Set("A",A);

    RCP<MultiVector> fineOnes = MultiVectorFactory::Build(A->getRowMap(),1);
    fineOnes->putScalar(1.0);
    fineLevel.Set("onesVector",fineOnes);


    RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
    RCP<TentativePFactory>    PtentFact = rcp(new TentativePFactory(UCAggFact));
    RCP<TransPFactory>        RFact = rcp(new TransPFactory(PtentFact));

    RCP<MueLu::MultiVectorTransferFactory<SC, LO, GO, NO, LMO> > mvtf = rcp(new MueLu::MultiVectorTransferFactory<SC, LO, GO, NO, LMO>("onesVector","R",RFact));

    coarseLevel.Request("onesVector",mvtf.get()); 
    coarseLevel.Request("R",RFact.get());
    coarseLevel.Request("P",PtentFact.get());

    mvtf->Build(fineLevel,coarseLevel);

    RCP<MultiVector> coarseOnes = coarseLevel.Get<RCP<MultiVector> >("onesVector",mvtf.get());
    Teuchos::Array<ST::magnitudeType> vn(1);
    coarseOnes->norm2(vn);

    TEST_FLOATING_EQUALITY(vn[0]*vn[0],((SC)fineOnes->getGlobalLength()),1e-12);
  } // Build test

  //------------------------------------------------------------------------------------------
  
  TEUCHOS_UNIT_TEST(MultiVectorTransferFactory, ThreeLevels)
  {
    out << "version: " << MueLu::Version() << std::endl;

    out << "Tests usage on a three-level hierarchy." << std::endl;

    GO nx = 199;
    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(nx);


    // Set up three level hierarchy.
    RCP<Hierarchy> H = rcp( new Hierarchy() );
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

    RCP<Level> fineLevel = H->GetLevel();
    fineLevel->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    fineLevel->Set("A",A);                      // set fine level matrix
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(),1);
    nullSpace->putScalar( (SC) 1.0);
    fineLevel->Set("Nullspace",nullSpace);       // set null space information for finest level

    RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
    UCAggFact->SetMinNodesPerAggregate(3);
    UCAggFact->SetMaxNeighAlreadySelected(0);
    UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
    UCAggFact->SetPhase3AggCreation(0.5);

    RCP<TentativePFactory> PFact = rcp(new TentativePFactory(UCAggFact)); //just using plain aggregation
    RCP<RFactory>          RFact = rcp( new TransPFactory(PFact) );
    RCP<RAPFactory>        AcFact = rcp( new RAPFactory(PFact,RFact) );
    H->SetMaxCoarseSize(1);

    Teuchos::ParameterList smootherParamList;
    smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smootherParamList.set("relaxation: sweeps", (LO) 1);
    smootherParamList.set("relaxation: damping factor", (SC) 1.0);
    RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother("RELAXATION", smootherParamList) );
    RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
    AcFact->setVerbLevel(Teuchos::VERB_HIGH);

    //set up the transfer factory
    RCP<MultiVector> fineOnes = MultiVectorFactory::Build(A->getRowMap(),1);
    fineOnes->putScalar(1.0);
    fineLevel->Set("onesVector",fineOnes);
    RCP<MueLu::MultiVectorTransferFactory<SC, LO, GO, NO, LMO> > mvtf = rcp(new MueLu::MultiVectorTransferFactory<SC, LO, GO, NO, LMO>("onesVector","R",RFact));
    AcFact->AddTransferFactory(mvtf);

    int maxLevels = 3;
    Teuchos::ParameterList status = H->FullPopulate(*PFact,*RFact, *AcFact,*SmooFact,0,maxLevels);

/*
    //FIXME we probably need to do some requests....
    coarseLevel.Request("onesVector",mvtf.get()); 
    coarseLevel.Request("R",RFact.get());
    coarseLevel.Request("P",PtentFact.get());
*/

/*
    RCP<MultiVector> coarseOnes = coarseLevel.Get<RCP<MultiVector> >("onesVector",mvtf.get());
    Teuchos::Array<ST::magnitudeType> vn(1);
    coarseOnes->norm2(vn);

    TEST_FLOATING_EQUALITY(vn[0]*vn[0],((SC)fineOnes->getGlobalLength()),1e-12);
*/
  } // ThreeLevels

} // namespace MueLuTests

