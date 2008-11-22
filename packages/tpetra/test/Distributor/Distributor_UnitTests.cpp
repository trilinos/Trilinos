#include "Teuchos_UnitTestHarness.hpp"

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"
#include "Tpetra_Distributor.hpp"
#include <Teuchos_Array.hpp>

// FINISH: test handling of null export in createFromSends

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Tpetra::Distributor;
  using Tpetra::DefaultPlatform;
  using Tpetra::Platform;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Array;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  template <typename T>
  T generateValue(T const x, T const y) {
    const T two = as<T>(2);
    // formula for z(x,y) = 0.5(x^2 + y^2 + 3x + y) + xy
    return(((x*x + y*y + x+x+x + y) / two) + (x*y));
  }

#define PRINT_VECTOR(v) \
   { \
     out << #v << ": "; \
     copy(v.begin(), v.end(), ostream_iterator<Ordinal>(out," ")); \
     out << endl; \
   }

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignord and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  template<class Ordinal>
  RCP<const Platform<Ordinal> > getDefaultPlatform()
  {
    if (testMpi) {
      return DefaultPlatform<Ordinal>::getPlatform();
    }
    return rcp(new Tpetra::SerialPlatform<Ordinal>());
  }

  //
  // UNIT TESTS
  // 

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, basic, Ordinal )
  {
    out << "sizeof(Ordinal): " << sizeof(Ordinal) << std::endl;
    out << "sizeof(std::size_t): " << sizeof(std::string::size_type) << std::endl;
    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    TEST_INEQUALITY_CONST( platform->createComm(), Teuchos::null );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromSendsContig, Ordinal )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;
    typedef typename std::vector<Ordinal>::size_type size_type;

    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    RCP<Teuchos::Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    const Ordinal ZERO = OT::zero();

    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    Teuchos_Ordinal numImports = 0;
    // fill exportImageIDs with {0,0, 1,1, 2,2, ... numImages-1,numImages-1}
    // two sends to each image, contiguous, in order
    vector<Ordinal> exportImageIDs(0);
    exportImageIDs.reserve(numImages*2);
    for(Ordinal i = ZERO; i < as<Ordinal>(numImages); ++i) {
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
    }

    // create from contiguous sends
    Distributor<Ordinal> distributor(comm);
    distributor.createFromSends(exportImageIDs, numImports);

    // tests
    TEST_EQUALITY(numImports, as<Teuchos_Ordinal>(2*numImages));
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<Ordinal>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<Ordinal>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), as<Ordinal>(numImages > 1 ? 2 : 0))
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<Ordinal>(2*numImages));
    {
      ArrayView<const Ordinal> imgFrom(distributor.getImagesFrom());
      ArrayView<const Ordinal> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      ArrayView<const Ordinal> lenFrom(distributor.getLengthsFrom());
      ArrayView<const Ordinal> lenTo(distributor.getLengthsTo());
      TEST_EQUALITY(lenFrom.size(),as<typename ArrayView<Ordinal>::Ordinal>(numImages));
      TEST_EQUALITY(lenTo.size()  ,as<typename ArrayView<Ordinal>::Ordinal>(numImages));
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], as<Ordinal>(2) );
        TEST_EQUALITY_CONST( lenTo[i],   as<Ordinal>(2) );
      }
    }
    TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), 0);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromSendsMixedContig, Ordinal )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;
    typedef typename std::vector<Ordinal>::size_type size_type;

    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    RCP<Teuchos::Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const Ordinal ZERO = OT::zero();

    if (numImages < 2) return;

    // even is black, odd is red
    bool even = ((myImageID % 2) == 0);

    // two exports to each image, including myself
    // on even imageIDs, send data contig
    // on odd imageIDs, send data non-contig
    Teuchos_Ordinal numImports = 0;
    vector<Ordinal> exportImageIDs(0);
    exportImageIDs.reserve(numImages*2);
    if (even) {
      // fill exportImageIDs with {0,0, 1,1, 2,2, ... numImages-1,numImages-1}
      for(Ordinal i = ZERO; i < as<Ordinal>(numImages); ++i) {
        exportImageIDs.push_back(i);
        exportImageIDs.push_back(i);
      }
    }
    else {
      // fill exportImageIDs with {0,0, 1,1, 2,2, ... numImages-1,numImages-1}
      for(Ordinal i = ZERO; i < as<Ordinal>(numImages); ++i) {
        exportImageIDs.push_back(i);
      }
      for(Ordinal i = ZERO; i < as<Ordinal>(numImages); ++i) {
        exportImageIDs.push_back(i);
      }
    }

    // create from sends, contiguous and non-contiguous
    Distributor<Ordinal> distributor(comm);
    distributor.createFromSends(exportImageIDs, numImports);

    // tests
    TEST_EQUALITY(numImports, as<Teuchos_Ordinal>(2*numImages));
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<Ordinal>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<Ordinal>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), as<Ordinal>(2));
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<Ordinal>(2*numImages));
    {
      ArrayView<const Ordinal> imgFrom(distributor.getImagesFrom());
      ArrayView<const Ordinal> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      ArrayView<const Ordinal> lenFrom(distributor.getLengthsFrom());
      ArrayView<const Ordinal> lenTo(distributor.getLengthsTo());
      TEST_EQUALITY(lenFrom.size(),as<typename ArrayView<Ordinal>::Ordinal>(numImages));
      TEST_EQUALITY(lenTo.size()  ,as<typename ArrayView<Ordinal>::Ordinal>(numImages));
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], as<Ordinal>(2) );
        TEST_EQUALITY_CONST( lenTo[i],   as<Ordinal>(2) );
      }
    }
    if (even) {
      TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), 0);
    }
    else {
      TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), as<typename ArrayView<Ordinal>::Ordinal>(2*numImages));
    }

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromSendsRedBlack, Ordinal )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;
    typedef typename std::vector<Ordinal>::size_type size_type;

    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    RCP<Teuchos::Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const Ordinal ZERO = OT::zero();
    const Ordinal  ONE = OT::one();

    if (numImages < 3) return;

    // partition world into red/black (according to imageID even/odd)
    // even is black, odd is red
    bool black = ((myImageID % 2) == 0);
    Ordinal numInMyPartition = ZERO;

    Teuchos_Ordinal numImports = 0;
    // fill exportImageIDs with all images from partition
    vector<Ordinal> exportImageIDs(0);
    if (black) {
      // evens
      for(Ordinal i = ZERO; i < as<Ordinal>(numImages); i+=2) {
        exportImageIDs.push_back(i);
        numInMyPartition++;
      }
    }
    else {
      // odds
      for(Ordinal i = ONE; i < as<Ordinal>(numImages); i+=2) {
        exportImageIDs.push_back(i);
        numInMyPartition++;
      }
    }

    // create from contiguous sends
    Distributor<Ordinal> distributor(comm);
    distributor.createFromSends(exportImageIDs, numImports);

    // tests
    TEST_EQUALITY(numImports, as<Teuchos_Ordinal>(numInMyPartition));
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), numInMyPartition-ONE);
    TEST_EQUALITY(distributor.getNumReceives(), numInMyPartition-ONE);
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), (numInMyPartition > 1 ? ONE : ZERO));
    TEST_EQUALITY(distributor.getTotalReceiveLength(), numInMyPartition);
    {
      ArrayView<const Ordinal> imgFrom(distributor.getImagesFrom());
      ArrayView<const Ordinal> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      ArrayView<const Ordinal> lenFrom(distributor.getLengthsFrom());
      ArrayView<const Ordinal> lenTo(distributor.getLengthsTo());
      TEST_EQUALITY(lenFrom.size(),as<typename ArrayView<Ordinal>::Ordinal>(numInMyPartition));
      TEST_EQUALITY(lenTo.size()  ,as<typename ArrayView<Ordinal>::Ordinal>(numInMyPartition));
      for (int i=0; i<numInMyPartition; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], as<Ordinal>(1) );
        TEST_EQUALITY_CONST( lenTo[i],   as<Ordinal>(1) );
      }
    }
    TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), 0);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromSendsContigNoself, Ordinal )
  {
    typedef typename std::vector<Ordinal>::size_type size_type;
    typedef Teuchos::OrdinalTraits<Ordinal> OT;

    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    RCP<Teuchos::Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const Ordinal ZERO = OT::zero();

    if (numImages < 2) return;

    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    Teuchos_Ordinal numImports = 0;
    // fill exportImageIDs with {0,1,...,myImageID-1,myImageID+1,...,numImages-1}
    // one send to each image, contiguous, in order, but not to myself
    vector<Ordinal> exportImageIDs(0);
    exportImageIDs.reserve(numImages-1);
    for(Ordinal i = ZERO; i < as<Ordinal>(myImageID); ++i) {
      exportImageIDs.push_back(i);
    }
    for(Ordinal i = as<Ordinal>(myImageID+1); i < as<Ordinal>(numImages); ++i) {
      exportImageIDs.push_back(i);
    }

    // create from contiguous sends
    Distributor<Ordinal> distributor(comm);
    distributor.createFromSends(exportImageIDs, numImports);

    // tests
    TEST_EQUALITY(numImports, as<Teuchos_Ordinal>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), false);
    TEST_EQUALITY(distributor.getNumSends(), as<Ordinal>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<Ordinal>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), as<Ordinal>(1));
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<Ordinal>(numImages-1));
    {
      ArrayView<const Ordinal> imgFrom(distributor.getImagesFrom());
      ArrayView<const Ordinal> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      ArrayView<const Ordinal> lenFrom(distributor.getLengthsFrom());
      ArrayView<const Ordinal> lenTo(distributor.getLengthsTo());
      TEST_EQUALITY(lenFrom.size(),as<typename ArrayView<Ordinal>::Ordinal>(numImages-1));
      TEST_EQUALITY(lenTo.size()  ,as<typename ArrayView<Ordinal>::Ordinal>(numImages-1));
      for (int i=0; i<numImages-1; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], as<Ordinal>(1) );
        TEST_EQUALITY_CONST( lenTo[i],   as<Ordinal>(1) );
      }
    }
    TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), 0);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromSendsContigUnordered, Ordinal )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;
    typedef typename std::vector<Ordinal>::size_type size_type;

    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    RCP<Teuchos::Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    const Ordinal ZERO = OT::zero();
    const Ordinal  ONE = OT::one();

    if (numImages < 3) return;

    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    Teuchos_Ordinal numImports = 0;
    // fill exportImageIDs with {0,0,0, 1,1,1, 2,2,2, ... numImages-1,numImages-1,numImages-1}
    // three sends to each image, out of order (even first, then odd)
    // only test if numImages > 2
    vector<Ordinal> exportImageIDs(0);
    exportImageIDs.reserve(numImages*3);
    // even first: {0,0,0, 2,2,2, 4,4,4, ...}
    for(Ordinal i = ZERO; i < as<Ordinal>(numImages); i+=2) {
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
    }
    // then odd: {1,1,1, 3,3,3, 5,5,5, ...}
    for(Ordinal i = ONE; i < as<Ordinal>(numImages); i+=2) {
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
    }

    // create from contiguous sends
    Distributor<Ordinal> distributor(comm);
    distributor.createFromSends(exportImageIDs, numImports);

    // tests
    TEST_EQUALITY(numImports, as<Teuchos_Ordinal>(3*numImages));
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<Ordinal>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<Ordinal>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), as<Ordinal>(3));
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<Ordinal>(3*numImages));
    {
      ArrayView<const Ordinal> imgFrom(distributor.getImagesFrom());
      ArrayView<const Ordinal> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      ArrayView<const Ordinal> lenFrom(distributor.getLengthsFrom());
      ArrayView<const Ordinal> lenTo(distributor.getLengthsTo());
      TEST_EQUALITY(lenFrom.size(),as<typename ArrayView<Ordinal>::Ordinal>(numImages));
      TEST_EQUALITY(lenTo.size()  ,as<typename ArrayView<Ordinal>::Ordinal>(numImages));
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], as<Ordinal>(3) );
        TEST_EQUALITY_CONST( lenTo[i],   as<Ordinal>(3) );
      }
    }
    TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), 0);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromSendsNonContig, Ordinal )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;
    typedef typename std::vector<Ordinal>::size_type size_type;

    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    RCP<Teuchos::Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    const Ordinal ZERO = OT::zero();

    if (numImages < 2) return;

    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    Teuchos_Ordinal numImports = 0;
    // fill exportImageIDs with {0, 1, 2, ... numImages-1,
    //                           0, 1, 2, ... numImages-1}
    vector<Ordinal> exportImageIDs(0);
    exportImageIDs.reserve(2*numImages);
    for(Ordinal i = ZERO; i < as<Ordinal>(numImages); ++i) {
      exportImageIDs.push_back(i);
    }
    for(Ordinal i = ZERO; i < as<Ordinal>(numImages); ++i) {
      exportImageIDs.push_back(i);
    }

    // create from non-contiguous sends
    Distributor<Ordinal> distributor(comm);
    distributor.createFromSends(exportImageIDs, numImports);

    // tests
    TEST_EQUALITY(numImports, as<Teuchos_Ordinal>(2*numImages));
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<Ordinal>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<Ordinal>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), as<Ordinal>(2));
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<Ordinal>(2*numImages));
    {
      ArrayView<const Ordinal> imgFrom(distributor.getImagesFrom());
      ArrayView<const Ordinal> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      ArrayView<const Ordinal> lenFrom(distributor.getLengthsFrom());
      ArrayView<const Ordinal> lenTo(distributor.getLengthsTo());
      TEST_EQUALITY(lenFrom.size(),as<typename ArrayView<Ordinal>::Ordinal>(numImages));
      TEST_EQUALITY(lenTo.size()  ,as<typename ArrayView<Ordinal>::Ordinal>(numImages));
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], as<Ordinal>(2) );
        TEST_EQUALITY_CONST( lenTo[i],   as<Ordinal>(2) );
      }
    }
    TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), as<typename ArrayView<Ordinal>::Ordinal>(2*numImages));

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Distributor, doPosts1, Ordinal, Packet )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;
    typedef Teuchos::ScalarTraits<Packet>   PT;
    typedef typename std::vector<Ordinal>::size_type size_type;

    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    RCP<Teuchos::Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const Ordinal ZERO = OT::zero();

    // send data to each image, including myself
    const Ordinal numExportIDs = as<Ordinal>(numImages); 
    Teuchos_Ordinal numRemoteIDs = 0;
    // fill exportImageIDs with {0, 1, 2, ... numImages-1}
    vector<Ordinal> exportImageIDs; 
    exportImageIDs.reserve(numExportIDs);
    for(Ordinal i = ZERO; i < numExportIDs; ++i) {
      exportImageIDs.push_back(i);
    }
    Distributor<Ordinal> distributor(comm);
    distributor.createFromSends(exportImageIDs, numRemoteIDs);

    // generate global random data set: each image sends 1 packet to each image
    // we need numImages*numImages "unique" values (we don't want redundant data allowing false positives)
    vector<Packet> exports(numImages*numImages);
    for (int i=0; i<numImages*numImages; i++) {
        exports[i] = PT::random();
    }
    // broadcast
    broadcast(*comm,0,arrayViewFromVector(exports));

    // pick a subset of entries to post
    Array<Packet> myExports(exports.begin()+myImageID*numImages,exports.begin()+(myImageID+1)*numImages);
    // do posts, one Packet to each image
    Array<Packet> imports(1*distributor.getTotalReceiveLength());
    distributor.doPostsAndWaits(myExports().getConst(), 1, imports());
    // imports[i] came from image i. it was element "myImageID" in his "myExports" vector. 
    // it corresponds to element i*numImages+myImageID in the global export vector
    // make a copy of the corresponding entries in the global vector, then compare these against the 
    // entries that I received
    vector<Packet> expectedImports(numImages);
    {
      typename vector<Packet>::iterator eI = expectedImports.begin(), 
                                         E = exports.begin()+myImageID;
      int left = numImages;
      while (true) {
        *eI = *E;
        if (--left > 0) {
          eI++;
          E += numImages;
        }
        else {
          break;
        }
      }
    }
    // check the values 
    TEST_COMPARE_ARRAYS(expectedImports,imports);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromReceives, Ordinal )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;
    typedef typename std::vector<Ordinal>::size_type size_type;

    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    RCP<Teuchos::Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const Ordinal ZERO = OT::zero();
    const Ordinal length = as<Ordinal>(numImages);

    // fill remoteImageIDs with {0, 1, 2, ... length-1}
    // we'll receive one GID from every image
    // 
    // fill remoteGIDs with row from generator
    // we'll receive generateValue(i,myImageID) from proc "i"
    // "i" sends us generateValue(i,myImageID)
    // similarly, we send generateValue(myImageID,i) to proc "i"
    vector<Ordinal> importImageIDs, importGIDs;
    importImageIDs.reserve(length);
    importGIDs.reserve(length);
    for(Ordinal i = ZERO; i < length; ++i) {
      importImageIDs.push_back(i);
      importGIDs.push_back( generateValue(i, as<Ordinal>(myImageID)) );
    }

    Distributor<Ordinal> distributor(comm);
    ArrayRCP<Ordinal> exportGIDs, exportImageIDs;
    distributor.createFromRecvs(importGIDs, importImageIDs, exportGIDs, exportImageIDs);
    
    TEST_EQUALITY(exportGIDs.size(), exportImageIDs.size());  // should *always* be the case

    vector<Ordinal> expectedGIDs;
    for (Ordinal i = ZERO; i < length; ++i) {
      expectedGIDs.push_back( generateValue(as<Ordinal>(myImageID),i) );
    }
    TEST_COMPARE_ARRAYS(importImageIDs, exportImageIDs);             
    TEST_COMPARE_ARRAYS(expectedGIDs, exportGIDs);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  //
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD


# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, basic, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromReceives, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsMixedContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsContigUnordered, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsContigNoself, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsRedBlack, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsNonContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Distributor, doPosts1, ORDINAL, double )

    UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

      // FINISH: add complex tests

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, basic, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromReceives, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsMixedContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsContigUnordered, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsContigNoself, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsRedBlack, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsNonContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Distributor, doPosts1, ORDINAL, char ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Distributor, doPosts1, ORDINAL, int ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Distributor, doPosts1, ORDINAL, double ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Distributor, doPosts1, ORDINAL, float )

    UNIT_TEST_GROUP_ORDINAL(int)
    typedef short int ShortInt;
    UNIT_TEST_GROUP_ORDINAL(ShortInt)
    typedef long int LongInt;
    UNIT_TEST_GROUP_ORDINAL(LongInt)
#   ifdef HAVE_TEUCHOS_LONG_LONG_INT
      typedef long long int LongLongInt;
      UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#   endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}

