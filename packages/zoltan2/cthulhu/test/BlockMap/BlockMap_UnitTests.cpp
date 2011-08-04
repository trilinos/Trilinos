#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>

#include "Cthulhu_BlockMap.hpp"
#include "Cthulhu_DefaultPlatform.hpp"

#ifdef HAVE_CTHULHU_TPETRA
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"
//#include "Tpetra_Map.hpp"
//#include "Tpetra_BlockMap.hpp"
#include "Cthulhu_TpetraMap.hpp"
#include "Cthulhu_TpetraBlockMap.hpp"
#endif

#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraMap.hpp"
#include "Cthulhu_EpetraBlockMap.hpp"
#endif

namespace {

  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::Tuple;
  using Teuchos::tuple;
  //using Tpetra::Map;
  //using Tpetra::BlockMap;
#ifdef HAVE_CTHULHU_TPETRA
  using Cthulhu::TpetraMap;
  using Cthulhu::TpetraBlockMap;
  using Tpetra::global_size_t;
#endif
#ifdef HAVE_CTHULHU_EPETRA
  using Cthulhu::EpetraMap;
  using Cthulhu::EpetraBlockMap;
#endif
  using Cthulhu::DefaultPlatform;
  using std::sort;
  using std::find;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

#define TEST_IS_COMPATIBLE(m1,m2,is_compat)               \
{                                                         \
    TEST_EQUALITY_CONST(m1.isCompatible(m1), true);       \
    TEST_EQUALITY_CONST(m2.isCompatible(m2), true);       \
    TEST_EQUALITY_CONST(m1.isCompatible(m2), is_compat);  \
    TEST_EQUALITY_CONST(m2.isCompatible(m1), is_compat);  \
}

#define TEST_IS_SAME_AS(m1,m2,is_sameas)               \
{                                                      \
    TEST_EQUALITY_CONST(m1.isSameAs(m1), true);        \
    TEST_EQUALITY_CONST(m2.isSameAs(m2), true);        \
    TEST_EQUALITY_CONST(m1.isSameAs(m2), is_sameas);   \
    TEST_EQUALITY_CONST(m2.isSameAs(m1), is_sameas);   \
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
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockMap, ContigConstBlkSize, BM, M, LO, GO )
  {
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with four entries per node
    // this map will have the following entries:
    Array<GO> myGlobal( tuple<GO>(myImageID*4, myImageID*4+1, myImageID*4+2, myImageID*4+3) );
    Array<LO>  myLocal( tuple<LO>(0,1,2,3) );

    const size_t numGlobalEntries = numImages*4;
    const GO indexBase = 0;
    RCP<const M> tmap(new M(numGlobalEntries,indexBase,comm));

    // create a BlockMap with 2 blocks per node, each block has size 2
    LO blkSize = 2;
    Array<GO> blkIDs( tuple<GO>(myImageID*2, myImageID*2+1) );
    Array<LO> blkLIDs( tuple<LO>(0, 1) );
    Array<LO> blkSzs( tuple<LO>(blkSize,blkSize) );
    Array<LO> firstPt( tuple<LO>(0,2) );

    BM blkmap(tmap, blkIDs(), blkSzs());
    TEST_EQUALITY_CONST(blkmap.getNodeNumBlocks(), 2);
    TEST_EQUALITY_CONST(blkmap.getLocalBlockID(blkIDs[0]), blkLIDs[0]);
    TEST_EQUALITY_CONST(blkmap.getLocalBlockID(blkIDs[1]), blkLIDs[1]);
    TEST_EQUALITY_CONST(blkmap.getLocalBlockSize(blkLIDs[0]), blkSzs[0]);
    TEST_EQUALITY_CONST(blkmap.getLocalBlockSize(blkLIDs[1]), blkSzs[1]);
    TEST_EQUALITY_CONST(blkmap.getFirstLocalPointInLocalBlock(blkLIDs[0]), firstPt[0]);
    TEST_EQUALITY_CONST(blkmap.getFirstLocalPointInLocalBlock(blkLIDs[1]), firstPt[1]);

    // create the same BlockMap with a different constructor:
    BM blkmap2(numImages*2, blkSize, indexBase, comm);
    // this BlockMap should pass the same tests:
    TEST_EQUALITY_CONST(blkmap2.getNodeNumBlocks(), 2);
    TEST_EQUALITY_CONST(blkmap2.getLocalBlockID(blkIDs[0]), blkLIDs[0]);
    TEST_EQUALITY_CONST(blkmap2.getLocalBlockID(blkIDs[1]), blkLIDs[1]);
    TEST_EQUALITY_CONST(blkmap2.getLocalBlockSize(blkLIDs[0]), blkSzs[0]);
    TEST_EQUALITY_CONST(blkmap2.getLocalBlockSize(blkLIDs[1]), blkSzs[1]);
    TEST_EQUALITY_CONST(blkmap2.getFirstLocalPointInLocalBlock(blkLIDs[0]), firstPt[0]);
    TEST_EQUALITY_CONST(blkmap2.getFirstLocalPointInLocalBlock(blkLIDs[1]), firstPt[1]);

    //and this BlockMap should have the same point-map:
#ifdef CTHULHU_NOT_IMPLEMENTED
    const M& tmapref = *tmap;
    const M& tmap2ref = *(blkmap2.getPointMap());
    TEST_IS_SAME_AS(tmapref, tmap2ref, true);
#endif

    // All procs fail if any proc fails 
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockMap, ContigNonConstBlkSize, BM, M, LO, GO )
  {
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with four entries per node
    // this map will have the following entries:
    Array<GO> myGlobal( tuple<GO>(myImageID*4, myImageID*4+1, myImageID*4+2, myImageID*4+3) );
    Array<LO>  myLocal( tuple<LO>(0,1,2,3) );

    const size_t numGlobalEntries = numImages*4;
    const GO indexBase = 0;
    RCP<const M> tmap(new M(numGlobalEntries,indexBase,comm));

    // create a BlockMap with 2 blocks per node, one block size 1 , the other size 3
    Array<GO> blkIDs( tuple<GO>(myImageID*2, myImageID*2+1) );
    Array<GO> points( tuple<GO>(myImageID*4, myImageID*4+1) );
    Array<LO> blkLIDs( tuple<LO>(0, 1) );
    Array<LO> blkSzs( tuple<LO>(1,3) );
    Array<LO> firstPt( tuple<LO>(0,1) );

    BM blkmap(tmap, blkIDs(), blkSzs());
    TEST_EQUALITY_CONST(blkmap.getNodeNumBlocks(), 2);
    TEST_EQUALITY_CONST(blkmap.getLocalBlockID(blkIDs[0]), blkLIDs[0]);
    TEST_EQUALITY_CONST(blkmap.getLocalBlockID(blkIDs[1]), blkLIDs[1]);
    TEST_EQUALITY_CONST(blkmap.getLocalBlockSize(blkLIDs[0]), blkSzs[0]);
    TEST_EQUALITY_CONST(blkmap.getLocalBlockSize(blkLIDs[1]), blkSzs[1]);
    TEST_EQUALITY_CONST(blkmap.getFirstLocalPointInLocalBlock(blkLIDs[0]), firstPt[0]);
    TEST_EQUALITY_CONST(blkmap.getFirstLocalPointInLocalBlock(blkLIDs[1]), firstPt[1]);

    // create the same BlockMap with a different constructor:
    BM blkmap2(numImages*2, blkIDs(), points(), blkSzs(), indexBase, comm);
    // this BlockMap should pass the same tests:
    TEST_EQUALITY_CONST(blkmap2.getNodeNumBlocks(), 2);
    TEST_EQUALITY_CONST(blkmap2.getLocalBlockID(blkIDs[0]), blkLIDs[0]);
    TEST_EQUALITY_CONST(blkmap2.getLocalBlockID(blkIDs[1]), blkLIDs[1]);
    TEST_EQUALITY_CONST(blkmap2.getLocalBlockSize(blkLIDs[0]), blkSzs[0]);
    TEST_EQUALITY_CONST(blkmap2.getLocalBlockSize(blkLIDs[1]), blkSzs[1]);
    TEST_EQUALITY_CONST(blkmap2.getFirstLocalPointInLocalBlock(blkLIDs[0]), firstPt[0]);
    TEST_EQUALITY_CONST(blkmap2.getFirstLocalPointInLocalBlock(blkLIDs[1]), firstPt[1]);

    //and this BlockMap should have the same point-map:
#ifdef CTHULHU_NOT_IMPLEMENTED
    const M& tmapref = *tmap;
    const M& tmap2ref = *(blkmap2.getPointMap());
    TEST_IS_SAME_AS(tmapref, tmap2ref, true);
#endif

    // All procs fail if any proc fails 
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockMap, ConstructorBadLengths1, BM, M, LO, GO )
  {
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with four entries per node
    // this map will have the following entries:
    Array<GO> myGlobal( tuple<GO>(myImageID*4, myImageID*4+1, myImageID*4+2, myImageID*4+3) );
    Array<LO>  myLocal( tuple<LO>(0,1,2,3) );

    const size_t numGlobalEntries = numImages*4;
    const GO indexBase = 0;
    RCP<const M> map(new M(numGlobalEntries,indexBase,comm));

    // create a BlockMap with 3 blocks per node, each block has size 2
    // (to be compatible with map, should be 2 blocks each with size 2)
    Array<GO> blkIDs( tuple<GO>(myImageID*3, myImageID*3+1, myImageID*3+2) );
    Array<LO> blkSzs( tuple<LO>(2,2, 2) );
    TEST_THROW(BM blkmap(map, blkIDs(), blkSzs()), std::runtime_error);

    // All procs fail if any proc fails 
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockMap, ConstructorBadLengths2, BM, M, LO, GO )
  {
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with four entries per node
    // this map will have the following entries:
    Array<GO> myGlobal( tuple<GO>(myImageID*4, myImageID*4+1, myImageID*4+2, myImageID*4+3) );
    Array<LO>  myLocal( tuple<LO>(0,1,2,3) );

    const size_t numGlobalEntries = numImages*4;
    const GO indexBase = 0;
    RCP<const M> map(new M(numGlobalEntries,indexBase,comm));

    // create a BlockMap with 2 blocks per node, each block has size 3
    // (to be compatible with map, should be 2 blocks each with size 2)
    Array<GO> blkIDs( tuple<GO>(myImageID*2, myImageID*2+1) );
    Array<LO> blkSzs( tuple<LO>(3,3) );
    TEST_THROW(BM blkmap(map, blkIDs(), blkSzs()), std::runtime_error);

    // All procs fail if any proc fails 
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  // 
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD


# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL_( BM, M, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMap, ContigConstBlkSize, BM, M, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMap, ContigNonConstBlkSize, BM, M, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMap, ConstructorBadLengths1, BM, M, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMap, ConstructorBadLengths2, BM, M, LO, GO )

#  define UNIT_TEST_GROUP_ORDINAL( LO, GO ) \
      typedef Cthulhu::TpetraBlockMap<LO,GO> TpetraBlockMap ## LO ## GO; \
      UNIT_TEST_GROUP_ORDINAL_(TpetraBlockMap ## LO ## GO, LO, GO)

    UNIT_TEST_GROUP_ORDINAL( char , int )
    UNIT_TEST_GROUP_ORDINAL( int , int )

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL_( BM, M, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMap, ContigConstBlkSize, BM, M, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMap, ContigNonConstBlkSize, BM, M, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMap, ConstructorBadLengths1, BM, M, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMap, ConstructorBadLengths2, BM, M, LO, GO )

#ifdef HAVE_CTHULHU_TPETRA
#  define UNIT_TEST_GROUP_ORDINAL( LO, GO ) \
      typedef Cthulhu::TpetraMap<LO,GO> TpetraMap ## LO ## GO; \
      typedef Cthulhu::TpetraBlockMap<LO,GO> TpetraBlockMap ## LO ## GO; \
      UNIT_TEST_GROUP_ORDINAL_(TpetraBlockMap ## LO ## GO, TpetraMap ## LO ## GO, LO, GO)

    // UNIT_TEST_GROUP_ORDINAL(char , int)

    //TODO: missing constructor UNIT_TEST_GROUP_ORDINAL_(EpetraBlockMap, EpetraMap, int , int)

      UNIT_TEST_GROUP_ORDINAL(int , int)
#endif
    // typedef short int ShortInt;
    // UNIT_TEST_GROUP_ORDINAL(ShortInt, int)

    // typedef long int LongInt;
    // UNIT_TEST_GROUP_ORDINAL(int , LongInt)

#   ifdef HAVE_TEUCHOS_LONG_LONG_INT
      // typedef long long int LongLongInt;
      // UNIT_TEST_GROUP_ORDINAL(char , LongLongInt)
      // UNIT_TEST_GROUP_ORDINAL(int , LongLongInt)
#   endif



# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
