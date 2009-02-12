#include "Teuchos_UnitTestHarness.hpp"

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"
#include "Tpetra_Directory.hpp"

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Tpetra::Map;
  using Tpetra::Directory;
  using Tpetra::DefaultPlatform;
  using Tpetra::Platform;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

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

  // test with a uniform, contiguous map of constant size
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, SmallUniformContig, LO, GO )
  {
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform map
    const GO numEntries = 1;
    M map(numEntries,0,comm);
    // create a directory
    D dir(map);

    // all GIDs
    const vector<GO> allGIDs(1,0);    // all GIDs (i.e., 0) are located on node 0
    const vector<int> expectedImageIDs(1,0); 
    const vector<LO> expectedLIDs(1,0);

    {
      vector<int> imageIDs(numEntries);              
      dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
    }

    {
      vector<int> imageIDs(numEntries);
      vector<LO> localIDs(numEntries); 
      dir.getDirectoryEntries(allGIDs,imageIDs,localIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
  }

  // test with a uniform, contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, UniformContig, LO, GO )
  {
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a uniform map
    const GO remainder = numImages/2;
    const GO numEntries = 2*numImages + remainder;
    M map(numEntries,0,comm);
    // create a directory
    D dir(map);

    // all GIDs
    vector<GO> allGIDs;
    allGIDs.reserve(numEntries);
    for (GO gid = 0; gid < numEntries; ++gid) {
      allGIDs.push_back(gid);
    }
    vector<int> expectedImageIDs;
    vector<LO> expectedLIDs;
    // expected image IDs and LIDs
    expectedImageIDs.reserve(numEntries);
    expectedLIDs.reserve(numEntries);
    GO remLeft = remainder;
    for (int id = 0; id < numImages; ++id) {
      expectedImageIDs.push_back(id);
      expectedImageIDs.push_back(id);
      expectedLIDs.push_back(0);
      expectedLIDs.push_back(1);
      if (remLeft) {
        expectedImageIDs.push_back(id);
        expectedLIDs.push_back(2);
        --remLeft;
      }
    }

    {
      vector<int> imageIDs(numEntries);                         
      dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
    }

    {
      vector<int> imageIDs(numEntries);
      vector<LO> localIDs(numEntries);
      dir.getDirectoryEntries(allGIDs,imageIDs,localIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
  }


  // test with a small contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, SmallContig, LO, GO )
  {
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const GO numEntries = as<GO>(numImages+1);
    // the last image gets two entries, others get none
    const LO numMyEntries = (myImageID == numImages-1 ? 2 : 1);
    M map(numEntries,numMyEntries,0,comm);
    // create a directory
    D dir(map);

    // all GIDs
    vector<GO> allGIDs;
    allGIDs.reserve(numEntries);
    for (GO gid = 0; gid < numEntries; ++gid) {
      allGIDs.push_back(gid);
    }
    vector<int> expectedImageIDs;
    vector<LO> expectedLIDs;
    // expected image IDs and LIDs
    expectedImageIDs.reserve(numEntries);
    expectedLIDs.reserve(numEntries);
    for (int id = 0; id < numImages; ++id) {
      expectedImageIDs.push_back(id);
      expectedLIDs.push_back(0);
      if (id == numImages-1) {
        expectedImageIDs.push_back(id);
        expectedLIDs.push_back(1);
      }
    }

    {
      vector<int> imageIDs(numEntries);                           
      dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
    }

    {
      vector<int> imageIDs(numEntries);
      vector<LO> localIDs(numEntries); 
      dir.getDirectoryEntries(allGIDs,imageIDs,localIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
  }


  // test with a contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, Contig, LO, GO )
  {
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // image i gets i+1 entries
    const LO numMyEntries = as<LO>(myImageID+1);
    // number of entries is (numImages+1)*numImages/2
    const GO numEntries = as<GO>((numImages*numImages+numImages)/2);
    M map(numEntries,numMyEntries,0,comm);
    // create a directory
    D dir(map);

    // all GIDs
    vector<GO> allGIDs;
    allGIDs.reserve(numEntries);
    for (GO gid = 0; gid < numEntries; ++gid) {
      allGIDs.push_back(gid);
    }
    vector<int> expectedImageIDs;
    vector<LO> expectedLIDs;
    // expected image IDs and LIDs
    expectedImageIDs.reserve(numEntries);
    expectedLIDs.reserve(numEntries);
    for (int id = 0; id < numImages; ++id) {
      for (Teuchos_Ordinal num = 0; num < id+1; ++num) {
        expectedImageIDs.push_back(id);
        expectedLIDs.push_back(as<LO>(num));
      }
    }

    {
      vector<int> imageIDs(numEntries);                           
      dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
    }

    {
      vector<int> imageIDs(numEntries);
      vector<LO> localIDs(numEntries); 
      dir.getDirectoryEntries(allGIDs,imageIDs,localIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
  }


  // test with a non-contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, NonContig, LO, GO )
  {
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // number of entries is 3*numImages
    // we will stripe the GIDs across images
    const GO numEntries = as<GO>(3*numImages);
    vector<GO> myGIDs(3);
    myGIDs[0] = as<GO>(myImageID);
    myGIDs[1] = as<GO>(myImageID + numImages);
    myGIDs[2] = as<GO>(myImageID + numImages*2);
    M map(numEntries,myGIDs,0,comm);
    // create a directory
    D dir(map);

    // all GIDs
    vector<GO> allGIDs;
    allGIDs.reserve(numEntries);
    for (GO gid = 0; gid < numEntries; ++gid) {
      allGIDs.push_back(gid);
    }
    vector<int> expectedImageIDs;
    vector<LO> expectedLIDs;
    // expected image IDs and LIDs
    expectedImageIDs.reserve(numEntries);
    expectedLIDs.reserve(numEntries);
    for (int i = 0; i < 3; ++i) {
      for (int id = 0; id < numImages; ++id) {
        expectedImageIDs.push_back(id);
        expectedLIDs.push_back(i);
      }
    }

    {
      vector<int> imageIDs(numEntries);                           
      dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
    }

    {
      vector<int> imageIDs(numEntries);
      vector<LO> localIDs(numEntries); 
      dir.getDirectoryEntries(allGIDs,imageIDs,localIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
  }


  //
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD


# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, SmallUniformContig, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, UniformContig, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, SmallContig, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, Contig, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, NonContig, LO, GO )

    UNIT_TEST_GROUP_ORDINAL( char , int )
    UNIT_TEST_GROUP_ORDINAL( int , int )

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, SmallUniformContig, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, UniformContig, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, SmallContig, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, Contig, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, NonContig, LO, GO )

    UNIT_TEST_GROUP_ORDINAL(char , int)
    UNIT_TEST_GROUP_ORDINAL(int , int)
    typedef short int ShortInt;
    UNIT_TEST_GROUP_ORDINAL(ShortInt , int)
    typedef long int LongInt;
    UNIT_TEST_GROUP_ORDINAL(int , LongInt)
#   ifdef HAVE_TEUCHOS_LONG_LONG_INT
      typedef long long int LongLongInt;
      UNIT_TEST_GROUP_ORDINAL(char , LongLongInt)
      UNIT_TEST_GROUP_ORDINAL(int , LongLongInt)
#   endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}

