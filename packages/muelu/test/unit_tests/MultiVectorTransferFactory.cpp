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

#include "MueLu_UseDefaultTypes.hpp"

#define XPETRA_ENABLED // == Gallery have to be build with the support of Xpetra matrices.
#include "MueLu_GalleryUtils.hpp"

#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {
  
  TEUCHOS_UNIT_TEST(MultiVectorTransferFactory, Constructor)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<TentativePFactory>    PtentFact = rcp(new TentativePFactory());
    RCP<MultiVectorTransferFactory> mvtf = rcp(new MultiVectorTransferFactory("Coordinates","P",PtentFact));
    TEST_EQUALITY(mvtf != Teuchos::null, true);
  } // Constructor test
  
  TEUCHOS_UNIT_TEST(MultiVectorTransferFactory, Build)
  {
    out << "version: " << MueLu::Version() << std::endl;

    out << "Tests the action of the transfer factory on a vector.  In this case, the transfer is the tentative" << std::endl;
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

    RCP<MultiVectorTransferFactory> mvtf = rcp(new MultiVectorTransferFactory("onesVector","R",RFact));

    coarseLevel.Request("onesVector",mvtf.get()); 
    coarseLevel.Request("R",RFact.get());
    coarseLevel.Request("P",PtentFact.get());

    mvtf->Build(fineLevel,coarseLevel);

    RCP<MultiVector> coarseOnes = coarseLevel.Get<RCP<MultiVector> >("onesVector",mvtf.get());
    Teuchos::Array<ST::magnitudeType> vn(1);
    coarseOnes->norm2(vn);

    TEST_FLOATING_EQUALITY(vn[0]*vn[0],((SC)fineOnes->getGlobalLength()),1e-12);
  } // Constructor test

} // namespace MueLuTests

