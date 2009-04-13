#include "Teuchos_UnitTestHarness.hpp"

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"
#include "Tpetra_Distributor.hpp"
#include <Teuchos_Array.hpp>

// FINISH: test handling of null export in createFromSends
// FINISH: test efficiency warnings if HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS for non-contig 

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Tpetra::Distributor;
  using Tpetra::DefaultPlatform;
  using Tpetra::Platform;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::Comm;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  int generateValue(int x, int y) {
    // formula for z(x,y) = 0.5(x^2 + y^2 + 3x + y) + xy
    return(((x*x + y*y + x+x+x + y) / 2) + (x*y));
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

  RCP<const Comm<int> > getDefaultComm()
  {
    RCP<Platform<double> > plat;
    if (testMpi) {
      plat = DefaultPlatform<double>::getPlatform();
    }
    else {
      plat = rcp(new Tpetra::SerialPlatform<double>());
    }
    return plat->getComm();
  }

  //
  // UNIT TESTS
  // 


  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsContig)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    Teuchos_Ordinal numImports = 0;
    // fill exportImageIDs with {0,0, 1,1, 2,2, ... numImages-1,numImages-1}
    // two sends to each image, contiguous, in order
    Array<int> exportImageIDs(0);
    exportImageIDs.reserve(numImages*2);
    for(int i=0; i < numImages; ++i) {
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
    }

    // create from contiguous sends
    Distributor distributor(comm);
    distributor.createFromSends(exportImageIDs, numImports);

    // tests
    TEST_EQUALITY(numImports, as<Teuchos_Ordinal>(2*numImages));
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), numImages-1);
    TEST_EQUALITY(distributor.getNumReceives(), numImages-1);
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), (numImages > 1 ? 2 : 0))
    TEST_EQUALITY(distributor.getTotalReceiveLength(), 2*numImages);
    {
      ArrayView<const int> imgFrom(distributor.getImagesFrom());
      ArrayView<const int> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      ArrayView<const Teuchos_Ordinal> lenFrom(distributor.getLengthsFrom());
      ArrayView<const Teuchos_Ordinal> lenTo(distributor.getLengthsTo());
      TEST_EQUALITY(lenFrom.size(),as<Teuchos_Ordinal>(numImages));
      TEST_EQUALITY(lenTo.size()  ,as<Teuchos_Ordinal>(numImages));
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], 2);
        TEST_EQUALITY_CONST( lenTo[i],   2);
      }
    }
    TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), 0);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
#ifndef HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsMixedContig)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();

    if (numImages < 2) return;

    // even is black, odd is red
    bool even = ((myImageID % 2) == 0);

    // two exports to each image, including myself
    // on even imageIDs, send data contig
    // on odd imageIDs, send data non-contig
    Teuchos_Ordinal numImports = 0;
    Array<int> exportImageIDs(0);
    exportImageIDs.reserve(numImages*2);
    if (even) {
      // fill exportImageIDs with {0,0, 1,1, 2,2, ... numImages-1,numImages-1}
      for(int i = 0; i < numImages; ++i) {
        exportImageIDs.push_back(i);
        exportImageIDs.push_back(i);
      }
    }
    else {
      // fill exportImageIDs with {0,0, 1,1, 2,2, ... numImages-1,numImages-1}
      for(int i = 0; i < numImages; ++i) {
        exportImageIDs.push_back(i);
      }
      for(int i = 0; i < numImages; ++i) {
        exportImageIDs.push_back(i);
      }
    }

    // create from sends, contiguous and non-contiguous
    Distributor distributor(comm);
    distributor.createFromSends(exportImageIDs, numImports);

    // tests
    TEST_EQUALITY(numImports, as<Teuchos_Ordinal>(2*numImages));
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<Teuchos_Ordinal>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<Teuchos_Ordinal>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), 2);
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<Teuchos_Ordinal>(2*numImages));
    {
      ArrayView<const int> imgFrom(distributor.getImagesFrom());
      ArrayView<const int> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      ArrayView<const Teuchos_Ordinal> lenFrom(distributor.getLengthsFrom());
      ArrayView<const Teuchos_Ordinal> lenTo(distributor.getLengthsTo());
      TEST_EQUALITY(lenFrom.size(),as<Teuchos_Ordinal>(numImages));
      TEST_EQUALITY(lenTo.size()  ,as<Teuchos_Ordinal>(numImages));
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], 2);
        TEST_EQUALITY_CONST( lenTo[i],   2);
      }
    }
    if (even) {
      TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), 0);
    }
    else {
      TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), as<Teuchos_Ordinal>(2*numImages));
    }

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }
#endif


  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsRedBlack)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();

    if (numImages < 3) return;

    // partition world into red/black (according to imageID even/odd)
    // even is black, odd is red
    bool black = ((myImageID % 2) == 0);
    Teuchos_Ordinal numInMyPartition = 0;

    Teuchos_Ordinal numImports = 0;
    // fill exportImageIDs with all images from partition
    Array<int> exportImageIDs(0);
    if (black) {
      // evens
      for(int i=0; i < numImages; i+=2) {
        exportImageIDs.push_back(i);
        numInMyPartition++;
      }
    }
    else {
      // odds
      for(int i=1; i < numImages; i+=2) {
        exportImageIDs.push_back(i);
        numInMyPartition++;
      }
    }

    // create from contiguous sends
    Distributor distributor(comm);
    distributor.createFromSends(exportImageIDs, numImports);

    // tests
    TEST_EQUALITY(numImports, numInMyPartition);
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), numInMyPartition-1);
    TEST_EQUALITY(distributor.getNumReceives(), numInMyPartition-1);
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), (numInMyPartition > 1 ? 1 : 0));
    TEST_EQUALITY(distributor.getTotalReceiveLength(), numInMyPartition);
    {
      ArrayView<const int> imgFrom(distributor.getImagesFrom());
      ArrayView<const int> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      ArrayView<const Teuchos_Ordinal> lenFrom(distributor.getLengthsFrom());
      ArrayView<const Teuchos_Ordinal> lenTo(distributor.getLengthsTo());
      TEST_EQUALITY(lenFrom.size(),numInMyPartition);
      TEST_EQUALITY(lenTo.size()  ,numInMyPartition);
      for (int i=0; i<numInMyPartition; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], as<Teuchos_Ordinal>(1) );
        TEST_EQUALITY_CONST( lenTo[i],   as<Teuchos_Ordinal>(1) );
      }
    }
    TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), 0);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsContigNoself)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();

    if (numImages < 2) return;

    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    Teuchos_Ordinal numImports = 0;
    // fill exportImageIDs with {0,1,...,myImageID-1,myImageID+1,...,numImages-1}
    // one send to each image, contiguous, in order, but not to myself
    Array<int> exportImageIDs(0);
    exportImageIDs.reserve(numImages-1);
    for (int i=0; i < myImageID; ++i) {
      exportImageIDs.push_back(i);
    }
    for (int i = myImageID+1; i < numImages; ++i) {
      exportImageIDs.push_back(i);
    }

    // create from contiguous sends
    Distributor distributor(comm);
    distributor.createFromSends(exportImageIDs, numImports);

    // tests
    TEST_EQUALITY(numImports, as<Teuchos_Ordinal>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), false);
    TEST_EQUALITY(distributor.getNumSends(), as<Teuchos_Ordinal>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<Teuchos_Ordinal>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), 1);
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<Teuchos_Ordinal>(numImages-1));
    {
      ArrayView<const int> imgFrom(distributor.getImagesFrom());
      ArrayView<const int> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      ArrayView<const Teuchos_Ordinal> lenFrom(distributor.getLengthsFrom());
      ArrayView<const Teuchos_Ordinal> lenTo(distributor.getLengthsTo());
      TEST_EQUALITY(lenFrom.size(),as<Teuchos_Ordinal>(numImages-1));
      TEST_EQUALITY(lenTo.size()  ,as<Teuchos_Ordinal>(numImages-1));
      for (int i=0; i<numImages-1; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], 1);
        TEST_EQUALITY_CONST( lenTo[i],   1);
      }
    }
    TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), 0);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsContigUnordered)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

    if (numImages < 3) return;

    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    Teuchos_Ordinal numImports = 0;
    // fill exportImageIDs with {0,0,0, 1,1,1, 2,2,2, ... numImages-1,numImages-1,numImages-1}
    // three sends to each image, out of order (even first, then odd)
    // only test if numImages > 2
    Array<int> exportImageIDs(0);
    exportImageIDs.reserve(numImages*3);
    // even first: {0,0,0, 2,2,2, 4,4,4, ...}
    for(int i=0; i < numImages; i+=2) {
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
    }
    // then odd: {1,1,1, 3,3,3, 5,5,5, ...}
    for(int i=1; i < numImages; i+=2) {
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
    }

    // create from contiguous sends
    Distributor distributor(comm);
    distributor.createFromSends(exportImageIDs, numImports);

    // tests
    TEST_EQUALITY(numImports, as<Teuchos_Ordinal>(3*numImages));
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<Teuchos_Ordinal>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<Teuchos_Ordinal>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), 3);
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<Teuchos_Ordinal>(3*numImages));
    {
      ArrayView<const int> imgFrom(distributor.getImagesFrom());
      ArrayView<const int> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      ArrayView<const Teuchos_Ordinal> lenFrom(distributor.getLengthsFrom());
      ArrayView<const Teuchos_Ordinal> lenTo(distributor.getLengthsTo());
      TEST_EQUALITY(lenFrom.size(),as<Teuchos_Ordinal>(numImages));
      TEST_EQUALITY(lenTo.size()  ,as<Teuchos_Ordinal>(numImages));
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST(lenFrom[i], 3);
        TEST_EQUALITY_CONST(lenTo[i],   3);
      }
    }
    TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), 0);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
#ifndef HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsNonContig)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

    if (numImages < 2) return;

    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    Teuchos_Ordinal numImports = 0;
    // fill exportImageIDs with {0, 1, 2, ... numImages-1,
    //                           0, 1, 2, ... numImages-1}
    Array<int> exportImageIDs(0);
    exportImageIDs.reserve(2*numImages);
    for(int i=0; i < numImages; ++i) {
      exportImageIDs.push_back(i);
    }
    for(int i=0; i < numImages; ++i) {
      exportImageIDs.push_back(i);
    }

    // create from non-contiguous sends
    Distributor distributor(comm);
    distributor.createFromSends(exportImageIDs, numImports);

    // tests
    TEST_EQUALITY(numImports, as<Teuchos_Ordinal>(2*numImages));
    TEST_EQUALITY_CONST(distributor.getSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<Teuchos_Ordinal>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<Teuchos_Ordinal>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), 2);
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<Teuchos_Ordinal>(2*numImages));
    {
      ArrayView<const int> imgFrom(distributor.getImagesFrom());
      ArrayView<const int> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
    {
      ArrayView<const Teuchos_Ordinal> lenFrom(distributor.getLengthsFrom());
      ArrayView<const Teuchos_Ordinal> lenTo(distributor.getLengthsTo());
      TEST_EQUALITY(lenFrom.size(),as<Teuchos_Ordinal>(numImages));
      TEST_EQUALITY(lenTo.size()  ,as<Teuchos_Ordinal>(numImages));
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], 2);
        TEST_EQUALITY_CONST( lenTo[i],   2);
      }
    }
    TEST_EQUALITY_CONST(distributor.getIndicesTo().size(), as<Teuchos_Ordinal>(2*numImages));

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }
#endif


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, doPosts1, Packet )
  {
    typedef Teuchos::ScalarTraits<Packet>   PT;

    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();

    // send data to each image, including myself
    const Teuchos_Ordinal numExportIDs = numImages; 
    Teuchos_Ordinal numRemoteIDs = 0;
    // fill exportImageIDs with {0, 1, 2, ... numImages-1}
    Array<int> exportImageIDs; 
    exportImageIDs.reserve(numExportIDs);
    for(int i=0; i < numExportIDs; ++i) {
      exportImageIDs.push_back(i);
    }
    Distributor distributor(comm);
    distributor.createFromSends(exportImageIDs, numRemoteIDs);

    // generate global random data set: each image sends 1 packet to each image
    // we need numImages*numImages "unique" values (we don't want redundant data allowing false positives)
    Array<Packet> exports(numImages*numImages);
    for (int i=0; i<numImages*numImages; i++) {
        exports[i] = PT::random();
    }
    // broadcast
    broadcast(*comm,0,exports());

    // pick a subset of entries to post
    Array<Packet> myExports(exports.begin()+myImageID*numImages,exports.begin()+(myImageID+1)*numImages);
    // do posts, one Packet to each image
    Array<Packet> imports(1*distributor.getTotalReceiveLength());
    distributor.doPostsAndWaits(myExports().getConst(), 1, imports());
    // imports[i] came from image i. it was element "myImageID" in his "myExports" vector. 
    // it corresponds to element i*numImages+myImageID in the global export vector
    // make a copy of the corresponding entries in the global vector, then compare these against the 
    // entries that I received
    Array<Packet> expectedImports(numImages);
    {
      typename Array<Packet>::iterator eI = expectedImports.begin(), 
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

    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const int length = numImages;

    // fill remoteImageIDs with {0, 1, 2, ... length-1}
    // we'll receive one GID from every image
    // 
    // fill remoteGIDs with row from generator
    // we'll receive generateValue(i,myImageID) from proc "i"
    // "i" sends us generateValue(i,myImageID)
    // similarly, we send generateValue(myImageID,i) to proc "i"
    Array<int> importImageIDs;
    Array<Ordinal> importGIDs;
    importImageIDs.reserve(length);
    importGIDs.reserve(length);
    for(int i=0; i < length; ++i) {
      importImageIDs.push_back(i);
      importGIDs.push_back( as<Ordinal>(generateValue(i, myImageID)) );
    }

    Distributor distributor(comm);
    ArrayRCP<int> exportImageIDs;
    ArrayRCP<Ordinal> exportGIDs;
    distributor.createFromRecvs<Ordinal>(importGIDs, importImageIDs, exportGIDs, exportImageIDs);
    
    TEST_EQUALITY(exportGIDs.size(), exportImageIDs.size());  // should *always* be the case

    Array<Ordinal> expectedGIDs;
    for (int i=0; i < length; ++i) {
      expectedGIDs.push_back( as<Ordinal>(generateValue(myImageID,i)) );
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
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromReceives, ORDINAL )

    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, doPosts1, double )
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_FLOAT( Distributor, doPosts1 )
    UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromReceives, ORDINAL )

    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( Distributor, doPosts1 )

    typedef short int ShortInt;
    UNIT_TEST_GROUP_ORDINAL(ShortInt)
    UNIT_TEST_GROUP_ORDINAL(int)
    typedef long int LongInt;
    UNIT_TEST_GROUP_ORDINAL(LongInt)
#   ifdef HAVE_TEUCHOS_LONG_LONG_INT
      typedef long long int LongLongInt;
      UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#   endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
