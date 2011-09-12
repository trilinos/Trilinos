#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(TentativePFactory, Constructor)
  {

    out << "version: " << MueLu::Version() << std::endl;

    RCP<TentativePFactory> tentPFact = rcp(new TentativePFactory);
    TEST_EQUALITY(tentPFact != Teuchos::null, true);

  } //Constructor

  TEUCHOS_UNIT_TEST(TentativePFactory, SetGetMethods)
  {
    out << "version: " << MueLu::Version() << std::endl;

    TentativePFactory tentPFact;

    bool flag = tentPFact.TentativeWithQR();
    TEST_EQUALITY(flag, false);
    tentPFact.TentativeWithQR(true);
    flag = tentPFact.TentativeWithQR();
    TEST_EQUALITY(flag, true);
  } //SetGetMethods

  //TODO test BuildP

#if 0 // todo check me
  TEUCHOS_UNIT_TEST(TentativePFactory, MakeTentative)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Test QR with user-supplied nullspace" << std::endl;

    Level fineLevel, coarseLevel;
    MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(199);
    fineLevel.Request("A");
    fineLevel.Set("A",A);

    // first iteration calls LAPACK QR
    // second iteration (with only one NS vector) exercises manual orthogonalization
    for (LO NSdim = 2; NSdim >= 1; --NSdim) {
      RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(),NSdim);
      nullSpace->randomize();
      fineLevel.Request("Nullspace",NULL); //FIXME putting this in to avoid error until Merge needs business
                                      //FIXME is implemented
      fineLevel.Set("Nullspace",nullSpace,NULL);
      RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
      UCAggFact->SetMinNodesPerAggregate(3);
      UCAggFact->SetMaxNeighAlreadySelected(0);
      UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
      UCAggFact->SetPhase3AggCreation(0.5);

      //fineLevel.Request("Aggregates",&UCAggFact); //FIXME putting this in to avoid error until Merge needs business
      RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory(UCAggFact));

      coarseLevel.Request("P",TentativePFact);  // request Ptent
      coarseLevel.Request("Nullspace",NULL);    // request coarse nullspace
      TentativePFact->DeclareInput(fineLevel,coarseLevel);
      TentativePFact->Build(fineLevel,coarseLevel);

      RCP<Operator> Ptent; 
      coarseLevel.Get("P",Ptent,TentativePFact);

      RCP<MultiVector> coarseNullSpace; 
      coarseLevel.Get("Nullspace",coarseNullSpace,NULL);

      //check interpolation
      RCP<MultiVector> PtN = MultiVectorFactory::Build(Ptent->getRangeMap(),NSdim);
      Ptent->apply(*coarseNullSpace,*PtN,Teuchos::NO_TRANS,1.0,0.0);

      RCP<MultiVector> diff = MultiVectorFactory::Build(A->getRowMap(),NSdim);
      diff->putScalar(0.0);

      coarseLevel.Release("P",TentativePFact); // release Ptent
      coarseLevel.Release("Nullspace",NULL);   // release coarse nullspace

      //diff = fineNS + (-1.0)*(P*coarseNS) + 0*diff
      diff->update(1.0,*nullSpace,-1.0,*PtN,0.0);

      Teuchos::Array<ST::magnitudeType> norms(NSdim);
      diff->norm2(norms);
      for (LO i=0; i<NSdim; ++i) {
        out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
        TEST_EQUALITY(norms[i]<1e-12, true);
      }

      //TODO check normalization of columns
      //RCP<const Map > coarseMap = coarseNullSpace->
      //TODO check orthogonality of columns

    } //for (LO NSdim = 1; NSdim <= 2; ++NSdim)

  } //MakeTentative  */
#endif

  TEUCHOS_UNIT_TEST(TentativePFactory, MakeTentativeUsingDefaultNullSpace)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Test QR when nullspace isn't supplied by user" << std::endl;

    Level fineLevel, coarseLevel; MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
    RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(199);
    fineLevel.Request("A",NULL);
    fineLevel.Set("A",A,NULL);

    RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
    UCAggFact->SetMinNodesPerAggregate(3);
    UCAggFact->SetMaxNeighAlreadySelected(0);
    UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
    UCAggFact->SetPhase3AggCreation(0.5);

    //UCAggFact->DeclareInput(fineLevel);
    //UCAggFact->Build(fineLevel);
    //fineLevel.Request("Aggregates"); //FIXME putting this in to avoid error until Merge needs business
    RCP<TentativePFactory> tentativePFact = rcp(new TentativePFactory(UCAggFact));

    coarseLevel.Request("P",tentativePFact);  // request Ptent
    coarseLevel.Request("Nullspace",NULL);    // request coarse nullspace
    tentativePFact->DeclareInput(fineLevel,coarseLevel);
    tentativePFact->Build(fineLevel,coarseLevel);

    RCP<Operator> Ptent; 
    coarseLevel.Get("P",Ptent,tentativePFact);

    RCP<MultiVector> coarseNullSpace; 
    coarseLevel.Get("Nullspace",coarseNullSpace,NULL);

    coarseLevel.Release("P",tentativePFact); // release Ptent
    coarseLevel.Release("Nullspace",NULL);   // release coarse nullspace

    //grab default fine level nullspace (vector of all ones)
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), 1);
    nullSpace->putScalar(1.0);

    //check interpolation
    LO NSdim = 1;
    RCP<MultiVector> PtN = MultiVectorFactory::Build(A->getRowMap(),NSdim);
    Ptent->apply(*coarseNullSpace,*PtN,Teuchos::NO_TRANS,1.0,0.0);

    RCP<MultiVector> diff = MultiVectorFactory::Build(A->getRowMap(),NSdim);
    diff->putScalar(0.0);

    //diff = fineNS + (-1.0)*(P*coarseNS) + 0*diff
    diff->update(1.0,*nullSpace,-1.0,*PtN,0.0);

    Teuchos::Array<ST::magnitudeType> norms(NSdim);
    diff->norm2(norms);
    for (LO i=0; i<NSdim; ++i) {
      out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
      TEST_EQUALITY(norms[i]<1e-12, true);
    }

  } //MakeTentative


}//namespace MueLuTests
