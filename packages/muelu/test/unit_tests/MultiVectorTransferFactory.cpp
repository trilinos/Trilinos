#include "Teuchos_UnitTestHarness.hpp"

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Utilities.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_MultiVectorTransferFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"

#define XPETRA_ENABLED // == Gallery have to be build with the support of Xpetra matrices.
#include "MueLu_GalleryUtils.hpp"

#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {
  
  TEUCHOS_UNIT_TEST(MultiVectorTransferFactory, Constructor)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<TentativePFactory>    Ptentfact = rcp(new TentativePFactory());
    RCP<MultiVectorTransferFactory> mvtf = rcp(new MultiVectorTransferFactory("Coordinates","P",Ptentfact));
    TEST_EQUALITY(mvtf != Teuchos::null, true);
  } // Constructor test
  
  TEUCHOS_UNIT_TEST(MultiVectorTransferFactory, Build)
  {
    out << "version: " << MueLu::Version() << std::endl;

    Level fineLevel, coarseLevel;
    TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
    GO nx = 199;
    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(nx);
    fineLevel.Set("A",A);

    //build coordinates
    Teuchos::ParameterList list;
    list.set("nx",nx);
    RCP<MultiVector> coordVector = MueLu::GalleryUtils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D",A->getRowMap(),list);
    fineLevel.Set("Coordinates",coordVector);


    RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
    RCP<TentativePFactory>    Ptentfact = rcp(new TentativePFactory(UCAggFact));

    RCP<MultiVectorTransferFactory> mvtf = rcp(new MultiVectorTransferFactory("Coordinates","P",Ptentfact));

    coarseLevel.Request("Coordinates",mvtf.get()); 
    coarseLevel.Request("P",Ptentfact.get());

    mvtf->Build(fineLevel,coarseLevel);

    RCP<MultiVector> coarseCoords = coarseLevel.Get<RCP<MultiVector> >("Coordinates",mvtf.get());
    TEST_EQUALITY(coarseCoords != Teuchos::null, true);
  } // Constructor test

} // namespace MueLuTests

