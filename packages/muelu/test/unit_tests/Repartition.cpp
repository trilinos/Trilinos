#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_UseDefaultTypes.hpp"

#include "MueLu_Repartition.hpp"

#include "MueLu_GalleryUtils.hpp"

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

    //FIXME switch to 4 processors to be consistent with other tests
    if (comm->getSize() != 5) {
      std::cout << std::endl;
      std::cout << "This test must be run on 5 processors!" << std::endl << std::endl;
      return;
    }

    Level level;
    RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
    level.SetFactoryManager(factoryHandler);
    int nx=3;
    int ny=5;
    GO numGlobalElements = nx*ny; //15

    // Describes the initial layout of unknowns across processors. Note that PID 3 initially has 0 unknowns.
    size_t numMyElements;
    switch(mypid) {
       case 0:
       case 1:
       case 2:
         numMyElements = 4;
         break;
       case 3:
         numMyElements = 0;
         break;
       case 4:
         numMyElements = 3;
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
    comm->barrier();

    /* Assign the partition that each unknown belongs to.  In this case,

       partition 0 has 2 unknowns 
       partition 1 has 3 unknowns 
       partition 2 has 4 unknowns 
       partition 3 has 6 unknowns .
    */
    switch (mypid)  {
      case 0:
        partitionThisDofBelongsTo[0] = 0;
        partitionThisDofBelongsTo[1] = 1;
        partitionThisDofBelongsTo[2] = 2;
        partitionThisDofBelongsTo[3] = 3;
        break;
      case 1:
        partitionThisDofBelongsTo[0] = 0;
        partitionThisDofBelongsTo[1] = 1;
        partitionThisDofBelongsTo[2] = 2;
        partitionThisDofBelongsTo[3] = 3;
        break;
      case 2:
        partitionThisDofBelongsTo[0] = 1;
        partitionThisDofBelongsTo[1] = 2;
        partitionThisDofBelongsTo[2] = 3;
        partitionThisDofBelongsTo[3] = 3;
        break;
      case 4:
        partitionThisDofBelongsTo[0] = 2;
        partitionThisDofBelongsTo[1] = 3;
        partitionThisDofBelongsTo[2] = 3;
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

    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    fos->setOutputToRootOnly(-1);
/*
    if (mypid == 0) std::cout << "\n==========\ninitial decomposition\n=========" << std::endl;
    sleep(1);
    decomposition->describe(*fos,Teuchos::VERB_EXTREME);

    sleep(1);
    comm->barrier();
*/

    level.Set<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("partition",decomposition);
    level.Set<GO>("number of partitions",4);

    RCP<Repartition> repart = rcp(new Repartition());
    level.Request("permMat",repart.get());  // request permutation matrix

    repart->Build(level);

    RCP<Operator> permMat;
    level.Get("permMat",permMat,repart.get());
    RCP<Vector> result = VectorFactory::Build(permMat->getRangeMap(),false);
    permMat->apply(*decompositionAsScalar,*result,Teuchos::NO_TRANS,1,0);
    bool testPassed=true;
    // the resulting vector should entries equal to the PID
    Teuchos::ArrayRCP<SC> resultData;
    if (result->getLocalLength() > 0)
      resultData = result->getDataNonConst(0);
    for (int i=0; i<resultData.size(); ++i) {
      if (resultData[i] != mypid) {
        testPassed = false;
        break;
      }
    }
    resultData = Teuchos::null;
    result->describe(*fos,Teuchos::VERB_EXTREME);

    // FIXME  do an all reduce on testPassed!!!
    TEST_EQUALITY(testPassed,true);

    //  Ptent->apply(*coarseNullSpace,*PtN,Teuchos::NO_TRANS,1.0,0.0);

  } //Build

}//namespace MueLuTests
