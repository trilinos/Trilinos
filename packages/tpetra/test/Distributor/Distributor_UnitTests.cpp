#include "Teuchos_UnitTestHarness.hpp"

// #include "../tpetra_test_util.hpp"
#include "Tpetra_ConfigDefs.hpp" // for map, vector, string, and iostream 

#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"

#include "Tpetra_Distributor.hpp"

// FINISH: test that createFromSends picks up contiguous sends, whether they are
// sorted or not. test, of course, that it also detects non-contiguous sends.

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
     out << "#v: "; \
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
    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    out << "platform = " << *platform << std::endl;
    TEST_INEQUALITY_CONST( platform->createComm(), Teuchos::null );
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromSends, Ordinal )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;

    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    RCP<Teuchos::Comm<Ordinal> > comm = platform->createComm();
    const Ordinal numImages = comm->getSize();
    const Ordinal ZERO = OT::zero();

    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    const Ordinal numExportIDs = as<Ordinal>(numImages); 
    Ordinal numRemoteIDs = ZERO;
    // fill exportImageIDs with {0, 1, 2, ... numImages-1}
    vector<Ordinal> exportImageIDs; 
    exportImageIDs.reserve(numExportIDs);
    for(Ordinal i = ZERO; i < numExportIDs; ++i) {
      exportImageIDs.push_back(i);
    }

    Distributor<Ordinal> distributor(comm);
    distributor.createFromSends(exportImageIDs, numRemoteIDs);
    TEST_EQUALITY(numRemoteIDs, numImages);
    TEST_EQUALITY(distributor.getTotalReceiveLength(), numImages);
    TEST_EQUALITY(distributor.getNumReceives(), numImages-1);
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), numImages-1);
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), as<Ordinal>(numImages > 1 ? 1 : 0));
    {
      vector<Ordinal> imgFrom(distributor.getImagesFrom());
      vector<Ordinal> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      const vector<Ordinal> & lenFrom = distributor.getLengthsFrom();
      const vector<Ordinal> & lenTo   = distributor.getLengthsTo();
      TEST_EQUALITY(as<Ordinal>(lenFrom.size()),numImages);
      TEST_EQUALITY(as<Ordinal>(lenTo.size())  ,numImages);
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], as<Ordinal>(1) );
        TEST_EQUALITY_CONST( lenTo[i],   as<Ordinal>(1) );
      }
    }
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(),as<Ordinal>(numImages > 1 ? 1 : 0));
    TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), 0);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromSendsNonContig, Ordinal )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;

    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    RCP<Teuchos::Comm<Ordinal> > comm = platform->createComm();
    const Ordinal numImages = comm->getSize();
    const Ordinal ZERO = OT::zero();

    if (numImages == 1) return;

    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    Ordinal numRemoteIDs = ZERO;
    // fill exportImageIDs with {0, 1, 2, ... numImages-1}
    vector<Ordinal> exportImageIDs; 
    exportImageIDs.reserve(2*numImages);
    for(Ordinal i = ZERO; i < numImages; ++i) {
      exportImageIDs.push_back(i);
    }
    for(Ordinal i = ZERO; i < numImages; ++i) {
      exportImageIDs.push_back(i);
    }

    Distributor<Ordinal> distributor(comm);
    distributor.createFromSends(exportImageIDs, numRemoteIDs);
    TEST_EQUALITY(numRemoteIDs, 2*numImages);
    TEST_EQUALITY(distributor.getTotalReceiveLength(), 2*numImages);
    TEST_EQUALITY(distributor.getNumReceives(), 2*numImages-2);
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), 2*numImages-2);
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), as<Ordinal>(2));
    {
      vector<Ordinal> imgFrom(distributor.getImagesFrom());
      vector<Ordinal> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      const vector<Ordinal> & lenFrom = distributor.getLengthsFrom();
      const vector<Ordinal> & lenTo   = distributor.getLengthsTo();
      TEST_EQUALITY(as<Ordinal>(lenFrom.size()),2*numImages);
      TEST_EQUALITY(as<Ordinal>(lenTo.size())  ,2*numImages);
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], as<Ordinal>(1) );
        TEST_EQUALITY_CONST( lenTo[i],   as<Ordinal>(1) );
      }
    }
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(),as<Ordinal>(2));
    TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), as<unsigned int>(2*numImages));

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Distributor, doPosts1, Ordinal, Packet )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;
    typedef Teuchos::ScalarTraits<Packet>   PT;

    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    RCP<Teuchos::Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const Ordinal ZERO = OT::zero();

    // send data to each image, including myself
    const Ordinal numExportIDs = as<Ordinal>(numImages); 
    Ordinal numRemoteIDs = ZERO;
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
    vector<Packet> myExports(&exports[myImageID*numImages],&exports[(myImageID+1)*numImages]);
    // do posts, one Packet to each image
    vector<Packet> imports;
    distributor.doPostsAndWaits(myExports, 1, imports);
    // imports[i] came from image i. it was element "myImageID" in his "myExports" vector. 
    // it corresponds to element i*numImages+myImageID in the global export vector
    // make a copy of the corresponding entries in the global vector, then compare these against the 
    // entries that I received
    vector<Packet> expectedImports(numImages);
    {
      typename vector<Packet>::iterator eI = expectedImports.begin(), 
                                         E = exports.begin()+myImageID;
      for (; eI != expectedImports.end();) {
        *eI = *E;
        eI++;
        E += numImages;
      }
    }
    // check the values 
    TEST_COMPARE_ARRAYS(expectedImports,imports);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromReceives, Ordinal )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;

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
    vector<Ordinal> exportGIDs, exportImageIDs;
    distributor.createFromRecvs(importGIDs, importImageIDs, exportGIDs, exportImageIDs);
    
    TEST_EQUALITY(exportGIDs.size(), exportImageIDs.size());  // should *always* be the case

    std::vector<Ordinal> expectedGIDs;
    for (Ordinal i = ZERO; i < length; ++i) {
      expectedGIDs.push_back( generateValue(myImageID,i) );
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
# define FAST_DEVELOPMENT_UNIT_TEST_BUILD


# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, basic, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromReceives, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSends, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsNonContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Distributor, doPosts1, ORDINAL, double )

    UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, basic, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromReceives, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSends, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromSendsNonContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Distributor, doPosts1, ORDINAL, char ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Distributor, doPosts1, ORDINAL, int ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Distributor, doPosts1, ORDINAL, double ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Distributor, doPosts1, ORDINAL, float )

    UNIT_TEST_GROUP_ORDINAL(char)
    typedef short int ShortInt;
    UNIT_TEST_GROUP_ORDINAL(ShortInt)
    UNIT_TEST_GROUP_ORDINAL_WITH_PAIRS(int)
    typedef long int LongInt;
    UNIT_TEST_GROUP_ORDINAL_WITH_PAIRS(LongInt)
#   ifdef HAVE_TEUCHOS_LONG_LONG_INT
      typedef long long int LongLongInt;
      UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#   endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}

