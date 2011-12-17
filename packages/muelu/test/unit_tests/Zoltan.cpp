#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_UseDefaultTypes.hpp"

#include "MueLu_Zoltan.hpp"

#include "MueLu_GalleryUtils.hpp"

#include "Xpetra_MultiVectorFactory.hpp"

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
    typedef Teuchos::ScalarTraits<Scalar> ST;

    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    Level level;
    RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
    level.SetFactoryManager(factoryHandler);
    //int nx=199;
    int nx=7;
    int ny=nx;
    GO numGlobalElements = nx*ny;
    size_t maxEntriesPerRow=10;

    // Populate CrsMatrix with random number of entries (up to maxEntriesPerRow) per row.
    RCP<const Map> map = MapFactory::createUniformContigMap(TestHelpers::Parameters::getLib(), numGlobalElements, comm);
    const size_t numMyElements = map->getNodeNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getNodeElementList();
    RCP<Operator> A = rcp(new CrsOperator(map, 1)); // Force underlying linear algebra library to allocate more
                                                    // memory on the fly.  While not super efficient, this
                                                    // ensures that no zeros are being stored.  Thus, from
                                                    // Zoltan's perspective the matrix is imbalanced.
    // Create a vector with random integer entries in [1,maxEntriesPerRow].
    ST::seedrandom(8675309);
    RCP<Xpetra::Vector<LO,LO,GO,NO> > entriesPerRow = Xpetra::VectorFactory<LO,LO,GO,NO>::Build(map,false);
    Teuchos::ArrayRCP<LO> eprData = entriesPerRow->getDataNonConst(0);
    for (Teuchos::ArrayRCP<LO>::iterator i=eprData.begin(); i!=eprData.end(); ++i) {
      *i = Teuchos::as<LO>(std::floor(((ST::random()+1)*0.5*maxEntriesPerRow)+1));
    }

    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    fos->setOutputToRootOnly(-1);
    //entriesPerRow->describe(*fos,Teuchos::VERB_EXTREME);

    Teuchos::Array<Scalar> vals(maxEntriesPerRow);
    Teuchos::Array<GO> cols(maxEntriesPerRow);
    for (size_t i = 0; i < numMyElements; ++i) {
      Teuchos::ArrayView<SC> av(&vals[0],eprData[i]);
      Teuchos::ArrayView<GO> iv(&cols[0],eprData[i]);
      //stick in ones for values
      for (LO j=0; j< eprData[i]; ++j) vals[j] = ST::one();
      //figure out valid column indices
      GO start = std::max(myGlobalElements[i]-eprData[i]+1,0);
      for (LO j=0; j< eprData[i]; ++j) cols[j] = start+j;
      A->insertGlobalValues(myGlobalElements[i], iv, av);
    }

    A->fillComplete();
    //A->describe(*fos, Teuchos::VERB_EXTREME);
    level.Set("A",A);

    //build coordinates
    RCP<const Map> rowMap = A->getRowMap();
    //RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    //fos->setOutputToRootOnly(-1);
    //rowMap->describe(*fos,Teuchos::VERB_EXTREME);
    Teuchos::ParameterList list;
    list.set("nx",nx);
    list.set("ny",ny);
    RCP<MultiVector> XYZ = MueLu::GalleryUtils::CreateCartesianCoordinates<SC,LO,GO,Map>("2D",rowMap,list);
    //std::cout << "pid " << comm->getRank() << ": XYZ local leng = " << XYZ->getLocalLength() << std::endl;
    level.Set("coordinates",XYZ);
    //XYZ->describe(*fos,Teuchos::VERB_EXTREME);

    RCP<ZoltanInterface> zoltan = rcp(new ZoltanInterface(comm));
    LO numPartitions = 5;
    zoltan->SetNumberOfPartitions(numPartitions);
    //zoltan->SetOutputLevel(0); //options are 0=none, 1=summary, 2=every pid prints
    zoltan->Build(level);

    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition = level.Get<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("partition");
    //decomposition->describe(*fos,Teuchos::VERB_EXTREME);

    //Create vector whose local length is the global number of partitions.
    //This vector will record the local number of nonzeros associated with each partition.
    Teuchos::Array<GO> parts(numPartitions);
    for (int i=0; i<numPartitions; ++i) parts[i] = i;
    Teuchos::ArrayView<GO> partsView(&parts[0],numPartitions);
    RCP<const Map> partitionMap = MapFactory::Build(TestHelpers::Parameters::getLib(),
                                                    Teuchos::OrdinalTraits<global_size_t>::invalid(), partsView,
                                                    map->getIndexBase(),comm);
    RCP<Xpetra::Vector<LO,LO,GO,NO> > localPartsVec = Xpetra::VectorFactory<LO,LO,GO,NO>::Build(partitionMap);

    //For the local rows in each partition, tally up the number of nonzeros.  This is what
    //Zoltan should be load-balancing.
    Teuchos::ArrayRCP<GO> lpvData = localPartsVec->getDataNonConst(0);
    Teuchos::ArrayRCP<const GO> decompData = decomposition->getData(0);
    for (size_t i=0; i<decomposition->getLocalLength();++i) {
      Teuchos::ArrayView<const LO> c;
      Teuchos::ArrayView<const SC> v;
      A->getLocalRowView(i,c,v);
      lpvData[decompData[i]] += v.size();
    }

    lpvData = Teuchos::null;
    decompData = Teuchos::null;

    //localPartsVec->describe(*fos,Teuchos::VERB_EXTREME);

    //Send the local nnz tallies to pid 0, which can report the global sums.
    size_t mysize=1;
    if (comm->getRank() == 0) mysize = numPartitions;
    RCP<const Map> globalTallyMap = MapFactory::Build(TestHelpers::Parameters::getLib(),
                                                Teuchos::OrdinalTraits<global_size_t>::invalid(),
                                                mysize,
                                                map->getIndexBase(),
                                                comm);
    RCP<Xpetra::Vector<LO,LO,GO,NO> > globalTallyVec = Xpetra::VectorFactory<LO,LO,GO,NO>::Build(globalTallyMap);
    RCP<const Export> exporter = ExportFactory::Build( partitionMap, globalTallyMap);
    globalTallyVec->doExport(*localPartsVec,*exporter,Xpetra::ADD);

    //globalTallyVec->describe(*fos,Teuchos::VERB_EXTREME);
    //FIXME cool ... this next line causes a hang if locally the globalyTallyVec has no data.
    //FIXME I get around this by making mysize (above) 1 instead of 0. Is this a bug or feature
    //FIXME in getData?
    ArrayRCP<const LO> gtvData = globalTallyVec->getData(0);
    ArrayRCP<GO> expectedResults(numPartitions);
    switch (comm->getSize()) {
       case 1:
         expectedResults[0] = 54;
         expectedResults[1] = 54;
         expectedResults[2] = 55;
         expectedResults[3] = 50;
         expectedResults[4] = 54;
         break;
       case 2:
         expectedResults[0] = 49;
         expectedResults[1] = 47;
         expectedResults[2] = 46;
         expectedResults[3] = 51;
         expectedResults[4] = 48;
         break;

       case 3:
         expectedResults[0] = 46;
         expectedResults[1] = 40;
         expectedResults[2] = 47;
         expectedResults[3] = 44;
         expectedResults[4] = 40;
         break;

       case 4:
         expectedResults[0] = 39;
         expectedResults[1] = 36;
         expectedResults[2] = 41;
         expectedResults[3] = 42;
         expectedResults[4] = 41;
         break;

       case 5:
         expectedResults[0] = 42;
         expectedResults[1] = 44;
         expectedResults[2] = 43;
         expectedResults[3] = 42;
         expectedResults[4] = 43;
         break;

       default:
         break;
    };

    for (int i=0; i<numPartitions; ++i) {
      // if pid 0, compare expectedResults[i] to gtvData[i].
      if (comm->getRank() == 0) TEST_EQUALITY( expectedResults[i], gtvData[i]);
    }

    //
    //Now write everything to a comma-separate list that ParaView can grok
    //
    Teuchos::ArrayRCP<const Scalar> X = XYZ->getData(0);
    Teuchos::ArrayRCP<const Scalar> Y = XYZ->getData(1);
    Teuchos::ArrayRCP<const GO> D = decomposition->getData(0);
    RCP<std::ofstream> outFile;
    std::string fileName = "zoltanResults.csv";

    //write header information
    if (comm->getRank() == 0) {
      outFile = rcp(new std::ofstream(fileName.c_str()));
      *outFile << "x coord, y coord, z coord, scalar" << std::endl;
    }
    comm->barrier();

    //append coordinates
    for (int j=0; j<comm->getSize(); ++j) {
      int mypid = comm->getRank();
      if (mypid == j) {
        outFile = rcp(new std::ofstream(fileName.c_str(),std::ios::app));
        for (int i=0; i < D.size(); ++i) {
          *outFile << X[i] << ", " << Y[i] << ", " << ST::zero() << ", " << D[i] << std::endl;
        }
      }
    } //for (int i=0; i<comm->getSize(); ++i)

    out << std::endl;
    out << "You can view the Zoltan decomposition in ParaView 3.10.1 or later:" << std::endl;
    out << "   1) Load the data file " << fileName << "." << std::endl;
    out << "   2) Run the filter Filters/ Alphabetical/ Table To Points." << std::endl;
    out << "   3) Tell ParaView what columns are the X, Y and Z coordinates." << std::endl;
    out << "   4) Split screen horizontally (Icon, top right)." << std::endl;
    out << "   5) Click on the eyeball in the Pipeline Browser to see the points." << std::endl;
    out << "   6) Under the Display tab, you can color points by scalar value and resize them." << std::endl;

  } //Build

}//namespace MueLuTests
