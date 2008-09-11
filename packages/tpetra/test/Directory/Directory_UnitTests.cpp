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

  // test with a uniform, contiguous map of constant size
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Directory, SmallUniformContig, Ordinal )
  {
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    // create a platform  
    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform->createComm();
    // create a uniform map
    const Ordinal numEntries = 1;
    Map<Ordinal> map(numEntries,ZERO,*platform);
    // create a directory
    Directory<Ordinal> dir(map);

    // all GIDs
    const vector<Ordinal> allGIDs(1,ZERO);    // all GIDs (i.e., 0) are located on node 0
    const vector<Ordinal> expectedImageIDs(1,ZERO); 
    const vector<Ordinal> expectedLIDs(1,ZERO);

    {
      vector<Ordinal> imageIDs(numEntries);              
      dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
    }

    {
      vector<Ordinal> imageIDs(numEntries), localIDs(numEntries); 
      dir.getDirectoryEntries(allGIDs,imageIDs,localIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
  }

  // test with a uniform, contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Directory, UniformContig, Ordinal )
  {
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal  ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal  TWO = ONE + ONE;
    // create a platform  
    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    // create a uniform map
    const Ordinal remainder = as<Ordinal>(numImages/2);
    const Ordinal numEntries = as<Ordinal>(2*numImages + remainder);
    Map<Ordinal> map(numEntries,ZERO,*platform);
    // create a directory
    Directory<Ordinal> dir(map);

    // all GIDs
    vector<Ordinal> allGIDs;
    allGIDs.reserve(numEntries);
    for (Ordinal gid = ZERO; gid < numEntries; ++gid) {
      allGIDs.push_back(gid);
    }
    vector<Ordinal> expectedImageIDs, expectedLIDs;
    // expected image IDs and LIDs
    expectedImageIDs.reserve(numEntries);
    expectedLIDs.reserve(numEntries);
    Ordinal remLeft = remainder;
    for (int id = 0; id < numImages; ++id) {
      expectedImageIDs.push_back(as<Ordinal>(id));
      expectedImageIDs.push_back(as<Ordinal>(id));
      expectedLIDs.push_back(ZERO);
      expectedLIDs.push_back(ONE);
      if (remLeft) {
        expectedImageIDs.push_back(as<Ordinal>(id));
        expectedLIDs.push_back(TWO);
        --remLeft;
      }
    }

    {
      vector<Ordinal> imageIDs(numEntries);                         
      dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
    }

    {
      vector<Ordinal> imageIDs(numEntries), localIDs(numEntries);
      dir.getDirectoryEntries(allGIDs,imageIDs,localIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
  }


  // test with a small contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Directory, SmallContig, Ordinal )
  {
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal  ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal  TWO = ONE + ONE;
    // create a platform  
    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const Ordinal numEntries = as<Ordinal>(numImages+1);
    // the last image gets two entries, others get none
    const Ordinal numMyEntries = (myImageID == numImages-1 ? TWO : ONE);
    Map<Ordinal> map(numEntries,numMyEntries,ZERO,*platform);
    // create a directory
    Directory<Ordinal> dir(map);

    // all GIDs
    vector<Ordinal> allGIDs;
    allGIDs.reserve(numEntries);
    for (Ordinal gid = ZERO; gid < numEntries; ++gid) {
      allGIDs.push_back(gid);
    }
    vector<Ordinal> expectedImageIDs, expectedLIDs;
    // expected image IDs and LIDs
    expectedImageIDs.reserve(numEntries);
    expectedLIDs.reserve(numEntries);
    for (int id = 0; id < numImages; ++id) {
      expectedImageIDs.push_back(as<Ordinal>(id));
      expectedLIDs.push_back(ZERO);
      if (id == numImages-1) {
        expectedImageIDs.push_back(as<Ordinal>(id));
        expectedLIDs.push_back(ONE);
      }
    }

    {
      vector<Ordinal> imageIDs(numEntries);                           
      dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
    }

    {
      vector<Ordinal> imageIDs(numEntries), localIDs(numEntries); 
      dir.getDirectoryEntries(allGIDs,imageIDs,localIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
  }


  // test with a contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Directory, Contig, Ordinal )
  {
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    // create a platform  
    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // image i gets i+1 entries
    const Ordinal numMyEntries = myImageID+1;
    // number of entries is (numImages+1)*numImages/2
    const Ordinal numEntries = as<Ordinal>((numImages*numImages+numImages)/2);
    Map<Ordinal> map(numEntries,numMyEntries,ZERO,*platform);
    // create a directory
    Directory<Ordinal> dir(map);

    // all GIDs
    vector<Ordinal> allGIDs;
    allGIDs.reserve(numEntries);
    for (Ordinal gid = ZERO; gid < numEntries; ++gid) {
      allGIDs.push_back(gid);
    }
    vector<Ordinal> expectedImageIDs, expectedLIDs;
    // expected image IDs and LIDs
    expectedImageIDs.reserve(numEntries);
    expectedLIDs.reserve(numEntries);
    for (int id = 0; id < numImages; ++id) {
      for (int num = 0; num < id+1; ++num) {
        expectedImageIDs.push_back(as<Ordinal>(id));
        expectedLIDs.push_back(as<Ordinal>(num));
      }
    }

    {
      vector<Ordinal> imageIDs(numEntries);                           
      dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
    }

    {
      vector<Ordinal> imageIDs(numEntries), localIDs(numEntries); 
      dir.getDirectoryEntries(allGIDs,imageIDs,localIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
  }


  // test with a non-contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Directory, NonContig, Ordinal )
  {
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    // create a platform  
    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform->createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // number of entries is 3*numImages
    // we will stripe the GIDs across images
    const Ordinal numEntries = as<Ordinal>(3*numImages);
    vector<Ordinal> myGIDs(3);
    myGIDs[0] = as<Ordinal>(myImageID);
    myGIDs[1] = as<Ordinal>(myImageID + numImages);
    myGIDs[2] = as<Ordinal>(myImageID + numImages*2);
    Map<Ordinal> map(numEntries,myGIDs,ZERO,*platform);
    // create a directory
    Directory<Ordinal> dir(map);

    // all GIDs
    vector<Ordinal> allGIDs;
    allGIDs.reserve(numEntries);
    for (Ordinal gid = ZERO; gid < numEntries; ++gid) {
      allGIDs.push_back(gid);
    }
    vector<Ordinal> expectedImageIDs, expectedLIDs;
    // expected image IDs and LIDs
    expectedImageIDs.reserve(numEntries);
    expectedLIDs.reserve(numEntries);
    for (Ordinal i = ZERO; i < as<Ordinal>(3); ++i) {
      for (Ordinal id = ZERO; id < numImages; ++id) {
        expectedImageIDs.push_back(id);
        expectedLIDs.push_back(i);
      }
    }

    {
      vector<Ordinal> imageIDs(numEntries);                           
      dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
    }

    {
      vector<Ordinal> imageIDs(numEntries), localIDs(numEntries); 
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

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Directory, SmallUniformContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Directory, UniformContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Directory, SmallContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Directory, Contig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Directory, NonContig, ORDINAL )

    UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Directory, SmallUniformContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Directory, UniformContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Directory, SmallContig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Directory, Contig, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Directory, NonContig, ORDINAL )

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

