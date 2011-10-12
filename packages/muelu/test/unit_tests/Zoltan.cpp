#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_UseDefaultTypes.hpp"

#include "MueLu_Zoltan.hpp"

#include "MueLu_GalleryUtils.hpp"

#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(Zoltan, Constructor)
  {

    out << "version: " << MueLu::Version() << std::endl;

    RCP<Zoltan> zoltan = rcp(new Zoltan);
    TEST_EQUALITY(zoltan != Teuchos::null, true);

  } //Constructor

  TEUCHOS_UNIT_TEST(Zoltan, Build)
  {

    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    RCP<ZoltanInterface> zoltan = rcp(new ZoltanInterface(comm));
    Level level;
    RCP<DefaultFactoryHandler> factoryHandler = rcp(new DefaultFactoryHandler());
    level.SetDefaultFactoryHandler(factoryHandler);
    //int nx=199;
    int nx=7;
    int ny=nx;
    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build2DPoisson(nx,ny);
    level.Set("A",A);
    //build coordinates
    RCP<const Map> rowMap = A->getRowMap();
    Teuchos::ParameterList list;
    list.set("nx",nx);
    list.set("ny",ny);
    RCP<MultiVector> XYZ = MueLu::GalleryUtils::CreateCartesianCoordinates<SC,LO,GO,Map>("2D",rowMap,list);
    level.Set("coordinates",XYZ);
    RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));
    //XYZ->describe(*fos,Teuchos::VERB_EXTREME);

    zoltan->SetNumberOfPartitions(2);
    zoltan->Build(level);

    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = level.Get<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("partition");
    decomposition->describe(*fos,Teuchos::VERB_EXTREME);

  } //Constructor

}//namespace MueLuTests
