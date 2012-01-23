#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_UseDefaultTypes.hpp"

#include "MueLu_Repartition.hpp"

#define XPETRA_ENABLED // == Gallery have to be build with the support of Xpetra matrices.
#include "MueLu_GalleryUtils.hpp"

#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_ExportFactory.hpp"

#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(Repartition, Constructor)
  {

    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    RCP<Repartition> repart = rcp(new Repartition());
    TEST_EQUALITY(repart != Teuchos::null, true);

  } //Constructor

  TEUCHOS_UNIT_TEST(Repartition, Build)
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;

    out << "version: " << MueLu::Version() << std::endl;
    out << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    int mypid = comm->getRank();

    if (comm->getSize() != 4) {
      std::cout << std::endl;
      std::cout << "This test must be run on 4 processors!" << std::endl << std::endl;
      return;
    }

    Level level;
    RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
    level.SetFactoryManager(factoryHandler);
    int nx=3;
    int ny=5;
    GO numGlobalElements = nx*ny; //15

    // Describes the initial layout of unknowns across processors. Note that PID 2 initially has 0 unknowns.
    size_t numMyElements;
    switch(mypid) {
       case 0:
         numMyElements = 6;
         break;
       case 1:
         numMyElements = 5;
         break;
       case 2:
         numMyElements = 0;
         break;
       case 3:
         numMyElements = 4;
         break;
    } //switch
    GO indexBase = 0;
    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, numMyElements, indexBase, comm);
    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(map,false);
    
    Teuchos::ParameterList matrixList;
    matrixList.set("nx",nx);
    matrixList.set("ny",ny);
    RCP<Operator> Op = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>("Laplace2D",map,matrixList);
    level.Set<RCP<Operator> >("A",Op);

    Teuchos::ArrayRCP<GO> partitionThisDofBelongsTo;
    if (decomposition->getLocalLength() > 0)
      partitionThisDofBelongsTo = decomposition->getDataNonConst(0);

    // Indicate the number of partitions that there should be.
    level.Set<GO>("number of partitions",3);

    /* Assign the partition that each unknown belongs to.  In this case,

       partition 0 has 6 unknowns 
       partition 1 has 4 unknowns 
       partition 2 has 5 unknowns 
       there is no partition 3
    */

    switch (mypid)  {
      case 0:
        partitionThisDofBelongsTo[0] = 0;
        partitionThisDofBelongsTo[1] = 1;
        partitionThisDofBelongsTo[2] = 2;
        partitionThisDofBelongsTo[3] = 0;
        partitionThisDofBelongsTo[4] = 2;
        partitionThisDofBelongsTo[5] = 2;
        break;
      case 1:
        partitionThisDofBelongsTo[0] = 0;
        partitionThisDofBelongsTo[1] = 1;
        partitionThisDofBelongsTo[2] = 0;
        partitionThisDofBelongsTo[3] = 2;
        partitionThisDofBelongsTo[4] = 1;
        break;
      case 3:
        partitionThisDofBelongsTo[0] = 2;
        partitionThisDofBelongsTo[1] = 0;
        partitionThisDofBelongsTo[2] = 0;
        partitionThisDofBelongsTo[3] = 1;
        break;
      default:
        break;
    } //switch

    RCP<Vector> decompositionAsScalar = VectorFactory::Build(map,false);
    Teuchos::ArrayRCP<SC> das;
    if (decompositionAsScalar->getLocalLength() > 0)
      das = decompositionAsScalar->getDataNonConst(0);
    for (int i=0; i<das.size(); ++i)
      das[i] = partitionThisDofBelongsTo[i];
    das = Teuchos::null;

    partitionThisDofBelongsTo = Teuchos::null;

    level.Set<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("partition",decomposition);

    RCP<Repartition> repart = rcp(new Repartition());
    level.Request("permMat",repart.get());  // request permutation matrix

    repart->Build(level);

    RCP<Operator> permMat;
    level.Get("permMat",permMat,repart.get());
    RCP<Vector> result = VectorFactory::Build(permMat->getRangeMap(),false);
    permMat->apply(*decompositionAsScalar,*result,Teuchos::NO_TRANS,1,0);
    int thisPidFailed=-1;
    // local entries in the resulting vector should be equal to the PID
    Teuchos::ArrayRCP<SC> resultData;
    if (result->getLocalLength() > 0)
      resultData = result->getDataNonConst(0);
    for (int i=0; i<resultData.size(); ++i) {
      if (resultData[i] != mypid) {
        thisPidFailed = mypid; // error condition
        break;
      }
    }
    resultData = Teuchos::null;

    //RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    //fos->setOutputToRootOnly(-1);
    //result->describe(*fos,Teuchos::VERB_EXTREME);

    // Test passes if all pid's return -1. If a pid hits an error, it sets
    // thisPidFailed equal to comm->getRank().  In this way, you can tell with --details=ALL
    // which pid failed.
    int whichPidFailed;
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,thisPidFailed,Teuchos::outArg(whichPidFailed));
    TEST_EQUALITY(whichPidFailed,-1);

    //  Ptent->apply(*coarseNullSpace,*PtN,Teuchos::NO_TRANS,1.0,0.0);

  } //Build

}//namespace MueLuTests
