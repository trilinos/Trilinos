#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_TypeTraits.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_BlockCrsGraph.hpp"

#include "Kokkos_SerialNode.hpp"
#ifdef HAVE_KOKKOS_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif
#ifdef HAVE_KOKKOS_THRUST
#include "Kokkos_ThrustGPUNode.hpp"
#endif

namespace {

  using Teuchos::as;
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::arcpClone;
  using Tpetra::BlockMap;
  using Tpetra::OptimizeOption;
  using Tpetra::DoOptimizeStorage;
  using Tpetra::DoNotOptimizeStorage;
  using Tpetra::DefaultPlatform;
  using Tpetra::global_size_t;
  using std::sort;
  using Teuchos::arrayView;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;
  using std::endl;
  using std::swap;
  using std::min;
  using std::max;
  using Teuchos::Array;
  using Teuchos::TypeTraits::is_same;
  using Teuchos::ArrayView;
  using Tpetra::BlockCrsGraph;
  using Tpetra::global_size_t;
  using Teuchos::arcp;
  using std::string;
  using std::unique;
  using Teuchos::tuple;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using Tpetra::ProfileType;
  using Tpetra::StaticProfile;
  using Tpetra::DynamicProfile;
  using Tpetra::Array_size_type;

  typedef DefaultPlatform::DefaultPlatformType::NodeType Node;

  using Kokkos::SerialNode;
  RCP<SerialNode> snode;
#ifdef HAVE_KOKKOS_TBB
  using Kokkos::TBBNode;
  RCP<TBBNode> tbbnode;
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
  using Kokkos::TPINode;
  RCP<TPINode> tpinode;
#endif
#ifdef HAVE_KOKKOS_THRUST
  using Kokkos::ThrustGPUNode;
  RCP<ThrustGPUNode> thrustnode;
#endif

  bool testMpi = true;
  double errorTolSlack = 1e+1;
  string filedir;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption(
        "filedir",&filedir,"Directory of expected input files.");
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
    RCP<const Comm<int> > ret;
    if (testMpi) {
      ret = DefaultPlatform::getDefaultPlatform().getComm();
    }
    else {
      ret = rcp(new Teuchos::SerialComm<int>());
    }
    return ret;
  }

  template <class Node>
  RCP<Node> getNode() {
    assert(false);
  }

  template <>
  RCP<SerialNode> getNode<SerialNode>() {
    if (snode == null) {
      Teuchos::ParameterList pl;
      snode = rcp(new SerialNode(pl));
    }
    return snode;
  }

#ifdef HAVE_KOKKOS_TBB
  template <>
  RCP<TBBNode> getNode<TBBNode>() {
    if (tbbnode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      tbbnode = rcp(new TBBNode(pl));
    }
    return tbbnode;
  }
#endif

#ifdef HAVE_KOKKOS_THREADPOOL
  template <>
  RCP<TPINode> getNode<TPINode>() {
    if (tpinode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      tpinode = rcp(new TPINode(pl));
    }
    return tpinode;
  }
#endif

#ifdef HAVE_KOKKOS_THRUST
  template <>
  RCP<ThrustGPUNode> getNode<ThrustGPUNode>() {
    if (thrustnode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      pl.set<int>("Verbose",1);
      thrustnode = rcp(new ThrustGPUNode(pl));
    }
    return thrustnode;
  }
#endif

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( BlockCrsGraph, ColMap1, LO, GO)
  {
    //This test fills a (block-tri-diagonal) block-crs-graph such that in parallel
    //the column-map should have an overlapping set of entries (i.e.,
    //different than the row-map), and verifies that the column-map is correct.

    RCP<Node> node = getNode<Node>();
    typedef BlockCrsGraph<LO,GO,Node> BGRAPH;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;
    GO maxGlobalBlock = numLocalBlocks*comm->getSize();
    const LO blockSize = 2;
    const size_t maxEntriesPerRow = 3;
    RCP<BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );
    
    //now set up the list of block-column-ids that we expect the
    //column-map to contain after fillComplete:
    size_t numLocalColBlocks = numLocalBlocks;
    if (comm->getRank() != 0) ++numLocalColBlocks;
    if (comm->getRank() != comm->getSize()-1) ++numLocalColBlocks;
    Array<GO> blockColIDs(numLocalColBlocks);
    typedef typename Array<GO>::size_type Tsize_t;
    Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
    GO first_row = blk_rows[0];
    Tsize_t offset = 0;
    if (comm->getRank() != 0) {
      blockColIDs[offset++] = first_row - 1;
    }
    GO last_row = 0;
    for(LO i=0; i<blk_rows.size(); ++i) {
      blockColIDs[offset++] = blk_rows[i];
      last_row = blk_rows[i];
    }
    if (offset < blockColIDs.size()) blockColIDs[offset++] = last_row + 1;

    // create the graph
    RCP<BGRAPH> bgrph = rcp( new BGRAPH(rowmap,maxEntriesPerRow,DynamicProfile) );
    for(int i=0; i<blk_rows.size(); ++i) {
      GO row = blk_rows[i];
      {
        GO col = row;
        bgrph->insertGlobalIndices(row, Teuchos::arrayView(&col, 1));
      }
      if (row > indexBase) {
        GO col = row - 1;
        bgrph->insertGlobalIndices(row, Teuchos::arrayView(&col, 1));
      }
      if (row < maxGlobalBlock-1) {
        GO col = row + 1;
        bgrph->insertGlobalIndices(row, Teuchos::arrayView(&col, 1));
      }
    }

    bgrph->fillComplete();
    RCP<const BlockMap<LO,GO,Node> > colmap = bgrph->getBlockColMap();
    ArrayView<const GO> blk_cols = colmap->getNodeBlockIDs();
    TEST_EQUALITY(blk_cols.size(), blockColIDs.size());
    TEST_COMPARE_ARRAYS(blk_cols, blockColIDs() );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( BlockCrsGraph, Queries, LO, GO)
  {
    //This test fills a (block-tri-diagonal) block-crs-graph such that
    //in parallel the column-map should have an overlapping set of
    //entries (i.e., different than the row-map), and verifies that
    //some queries work for attributes such as row-lengths, etc.

    RCP<Node> node = getNode<Node>();
    typedef BlockCrsGraph<LO,GO,Node> BGRAPH;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;
    GO maxGlobalBlock = numLocalBlocks*comm->getSize();
    const LO blockSize = 2;
    const size_t maxEntriesPerRow = 3;
    RCP<BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );
    
    //now set up the list of block-column-ids that we expect the
    //column-map to contain after fillComplete:
    size_t numLocalColBlocks = numLocalBlocks;
    if (comm->getRank() != 0) ++numLocalColBlocks;
    if (comm->getRank() != comm->getSize()-1) ++numLocalColBlocks;
    Array<GO> blockColIDs(numLocalColBlocks);
    typedef typename Array<GO>::size_type Tsize_t;
    Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
    GO first_row = blk_rows[0];
    Tsize_t offset = 0;
    if (comm->getRank() != 0) {
      blockColIDs[offset++] = first_row - 1;
    }
    GO last_row = 0;
    for(LO i=0; i<blk_rows.size(); ++i) {
      blockColIDs[offset++] = blk_rows[i];
      last_row = blk_rows[i];
    }
    if (offset < blockColIDs.size()) blockColIDs[offset++] = last_row + 1;

    // create the graph
    RCP<BGRAPH> bgrph = rcp( new BGRAPH(rowmap,maxEntriesPerRow,DynamicProfile) );
    for(int i=0; i<blk_rows.size(); ++i) {
      GO row = blk_rows[i];
      {
        GO col = row;
        bgrph->insertGlobalIndices(row, Teuchos::arrayView(&col, 1));
      }
      if (row > indexBase) {
        GO col = row - 1;
        bgrph->insertGlobalIndices(row, Teuchos::arrayView(&col, 1));
      }
      if (row < maxGlobalBlock-1) {
        GO col = row + 1;
        bgrph->insertGlobalIndices(row, Teuchos::arrayView(&col, 1));
      }
    }

    bgrph->fillComplete();
    RCP<const BlockMap<LO,GO,Node> > colmap = bgrph->getBlockColMap();
    ArrayView<const GO> blk_cols = colmap->getNodeBlockIDs();
    TEST_EQUALITY(blk_cols.size(), blockColIDs.size());
    TEST_COMPARE_ARRAYS(blk_cols, blockColIDs() );

    size_t map_blk_elems = blk_rows.size();
    TEST_EQUALITY(bgrph->getNodeNumBlockRows(), map_blk_elems);
    TEST_EQUALITY(bgrph->getGlobalNumBlockRows(), rowmap->getGlobalNumBlocks());

    for(int i=0; i<blk_rows.size(); ++i) {
      GO row = blk_rows[i];
      if (row > indexBase && row < maxGlobalBlock-1) {
        size_t row_len = bgrph->getGlobalBlockRowLength(row);
        size_t expected_row_len = 3;
        TEST_EQUALITY(row_len, expected_row_len);
      }
      else {
        size_t row_len = bgrph->getGlobalBlockRowLength(row);
        size_t expected_row_len = 2;
        TEST_EQUALITY(row_len, expected_row_len);
      }
    }
  }

  // 
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_LO_GO( LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( BlockCrsGraph, ColMap1  , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( BlockCrsGraph, Queries  , LO, GO )

     UNIT_TEST_GROUP_LO_GO(int,int)
// #ifndef FAST_DEVELOPMENT_UNIT_TEST_BUILD
// #ifdef HAVE_TEUCHOS_LONG_LONG_INT
//         typedef long long int LongLongInt;
//         UNIT_TEST_GROUP_LO_GO(short,LongLongInt)
// #else
//         typedef long int LongInt;
//         UNIT_TEST_GROUP_LO_GO(short,LongInt)
// #endif
// # endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
