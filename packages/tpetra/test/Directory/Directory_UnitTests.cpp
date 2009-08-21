#include "Teuchos_UnitTestHarness.hpp"

#include <Teuchos_as.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Directory.hpp"

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Tpetra::Map;
  using Tpetra::Directory;
  using Tpetra::DefaultPlatform;
  using Teuchos::Array;
  using Teuchos::tuple;
  using std::sort;
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
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  // 

  // test with a uniform, contiguous map of constant size
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, BadSize, LO, GO )
  {
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform map
    const GO numEntries = 2;
    RCP<M> map = rcp(new M(numEntries,0,comm));
    // create a directory
    D dir(map);
    
    Array<int> imageIDs(2);
    Array<LO> localIDs(2);
    TEST_THROW( dir.getDirectoryEntries(tuple<GO>(0,1), imageIDs(0,1)), std::runtime_error );
    TEST_THROW( dir.getDirectoryEntries(tuple<GO>(0,1), imageIDs(0,1), localIDs(0,1)), std::runtime_error );
    TEST_THROW( dir.getDirectoryEntries(tuple<GO>(0,1), imageIDs(0,2), localIDs(0,1)), std::runtime_error );
    TEST_THROW( dir.getDirectoryEntries(tuple<GO>(0,1), imageIDs(0,1), localIDs(0,2)), std::runtime_error );
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  // test with a uniform, contiguous map of constant size
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, SmallUniformContig, LO, GO )
  {
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    const LO LINV = OrdinalTraits<LO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform map
    const GO numEntries = 1;
    RCP<M> map = rcp( new M(numEntries,0,comm) );
    // create a directory
    D dir(map);
    {
      bool invalid;
      Array<int> imageIDs(numEntries);
      Array<LO>  localIDs(numEntries); 
      invalid = dir.getDirectoryEntries(tuple<GO>(0),imageIDs);
      TEST_EQUALITY_CONST( invalid, false );
      TEST_COMPARE_ARRAYS( tuple<int>(0), imageIDs );
      invalid = dir.getDirectoryEntries(tuple<GO>(0),imageIDs,localIDs);
      TEST_EQUALITY_CONST( invalid, false );
      TEST_COMPARE_ARRAYS( tuple<int>(0), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0), localIDs );
    }
    {
      bool invalid;
      Array<int> imageIDs(numEntries+1);
      Array<LO>  localIDs(numEntries+1);
      invalid = dir.getDirectoryEntries(tuple<GO>(0,1), imageIDs);
      TEST_EQUALITY_CONST( invalid, true );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      invalid = dir.getDirectoryEntries(tuple<GO>(0,1),imageIDs,localIDs);
      TEST_EQUALITY_CONST( invalid, true );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0,LINV), localIDs );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  // test with a uniform, contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, UniformContig, LO, GO )
  {
    const LO LINV = OrdinalTraits<LO>::invalid();
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a uniform map
    const GO remainder = numImages/2;
    const GO numEntries = 2*numImages + remainder;
    RCP<M> map = rcp(new M(numEntries,0,comm));
    // create a directory
    D dir(map);
    // all GIDs
    Array<GO> allGIDs(numEntries);
    for (GO gid = 0; gid < numEntries; ++gid) {
      allGIDs[gid] = gid;
    }
    Array<int> expectedImageIDs;
    Array<LO>  expectedLIDs;
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
      bool invalid;
      Array<int> imageIDs(numEntries);
      Array<LO> localIDs(numEntries);
      invalid = dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_EQUALITY_CONST( invalid, false );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      invalid = dir.getDirectoryEntries(allGIDs,imageIDs,localIDs);
      TEST_EQUALITY_CONST( invalid, false );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
    {
      bool invalid;
      Array<int> imageIDs(2);
      Array<LO> localIDs(2);
      invalid = dir.getDirectoryEntries( tuple<GO>(0,numEntries),imageIDs);
      TEST_EQUALITY_CONST( invalid, true );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      invalid = dir.getDirectoryEntries( tuple<GO>(0,numEntries),imageIDs,localIDs);
      TEST_EQUALITY_CONST( invalid, true );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0,LINV), localIDs );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  // test with a small contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, SmallContig, LO, GO )
  {
    const LO LINV = OrdinalTraits<LO>::invalid();
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const GO numEntries = as<GO>(numImages+1);
    // the last image gets two entries, others get one
    const LO numMyEntries = (myImageID == numImages-1 ? 2 : 1);
    RCP<M> map = rcp(new M(numEntries,numMyEntries,0,comm));
    // create a directory
    D dir(map);
    // all GIDs
    Array<GO> allGIDs;
    allGIDs.reserve(numEntries);
    for (GO gid = 0; gid < numEntries; ++gid) 
    {
      allGIDs.push_back(gid);
    }
    Array<int> expectedImageIDs;
    Array<LO> expectedLIDs;
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
      bool invalid;
      Array<int> imageIDs(numEntries);
      Array<LO> localIDs(numEntries); 
      invalid = dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_EQUALITY_CONST( invalid, false );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      invalid = dir.getDirectoryEntries(allGIDs,imageIDs,localIDs);
      TEST_EQUALITY_CONST( invalid, false );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
    {
      bool invalid;
      Array<int> imageIDs(2);
      Array<LO> localIDs(2); 
      invalid = dir.getDirectoryEntries(tuple<GO>(0,numEntries),imageIDs);
      TEST_EQUALITY_CONST( invalid, true );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      invalid = dir.getDirectoryEntries(tuple<GO>(0,numEntries),imageIDs,localIDs);
      TEST_EQUALITY_CONST( invalid, true );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0,LINV), localIDs );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  // test with a contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, Contig, LO, GO )
  {
    const LO LINV = OrdinalTraits<LO>::invalid();
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
    RCP<M> map = rcp(new M(numEntries,numMyEntries,0,comm));
    // create a directory
    D dir(map);
    // all GIDs
    Array<GO> allGIDs(numEntries);
    for (GO gid = 0; gid < numEntries; ++gid) {
      allGIDs[gid] = gid;
    }
    Array<int> expectedImageIDs;
    Array<LO> expectedLIDs;
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
      bool invalid;
      Array<int> imageIDs(numEntries);                           
      Array<LO> localIDs(numEntries); 
      invalid = dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_EQUALITY_CONST( invalid, false );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      invalid = dir.getDirectoryEntries(allGIDs,imageIDs,localIDs);
      TEST_EQUALITY_CONST( invalid, false );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
    {
      bool invalid;
      Array<int> imageIDs(2);                           
      Array<LO> localIDs(2); 
      invalid = dir.getDirectoryEntries( tuple<GO>(0,numEntries) ,imageIDs);
      TEST_EQUALITY_CONST( invalid, true );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      invalid = dir.getDirectoryEntries( tuple<GO>(0,numEntries) ,imageIDs,localIDs);
      TEST_EQUALITY_CONST( invalid, true );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0,LINV), localIDs );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  // test with a non-contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, NonContig, LO, GO )
  {
    const LO LINV = OrdinalTraits<LO>::invalid();
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // number of entries is 3*numImages
    // we will stripe the GIDs across images
    const GO numEntries = as<GO>(3*numImages);
    RCP<M> map = rcp(new M(numEntries, tuple<GO>(myImageID, myImageID+numImages, myImageID+2*numImages) ,0,comm));
    // create a directory
    D dir(map);
    // all GIDs
    Array<GO> allGIDs;
    allGIDs.reserve(numEntries);
    for (GO gid = 0; gid < numEntries; ++gid) {
      allGIDs.push_back(gid);
    }
    Array<int> expectedImageIDs;
    Array<LO> expectedLIDs;
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
      bool invalid;
      Array<int> imageIDs(numEntries);                           
      Array<LO> localIDs(numEntries); 
      invalid = dir.getDirectoryEntries(allGIDs,imageIDs);
      TEST_EQUALITY_CONST( invalid, false );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      invalid = dir.getDirectoryEntries(allGIDs,imageIDs,localIDs);
      TEST_EQUALITY_CONST( invalid, false );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
    {
      bool invalid;
      Array<int> imageIDs(2);                           
      Array<LO> localIDs(2); 
      invalid = dir.getDirectoryEntries(tuple<GO>(0,numEntries), imageIDs);
      TEST_EQUALITY_CONST( invalid, true );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      invalid = dir.getDirectoryEntries(tuple<GO>(0,numEntries),imageIDs,localIDs);
      TEST_EQUALITY_CONST( invalid, true );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0,LINV), localIDs );
    }
    // All procs fail if any node fails
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

#   define UNIT_TEST_GROUP_ORDINAL( LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, SmallUniformContig, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, UniformContig, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, SmallContig, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, Contig, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, NonContig, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, BadSize, LO, GO )

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

    UNIT_TEST_GROUP_ORDINAL( char , int )
    UNIT_TEST_GROUP_ORDINAL( int , int )

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

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

