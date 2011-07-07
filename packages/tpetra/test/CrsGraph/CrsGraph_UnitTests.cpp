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
#include "Tpetra_CrsGraph.hpp"

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::arcpClone;
  using Tpetra::Map;
  using Tpetra::OptimizeOption;
  using Tpetra::DoOptimizeStorage;
  using Tpetra::DoNotOptimizeStorage;
  using Tpetra::DefaultPlatform;
  using Tpetra::createUniformContigMap;
  using Tpetra::createContigMap;
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
  using Tpetra::CrsGraph;
  using Tpetra::RowGraph;
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

  bool testMpi = true;
  double errorTolSlack = 1e+1;
  string filedir;

#define STD_TESTS(graph) \
  { \
    using Teuchos::outArg; \
    RCP<const Comm<int> > STCOMM = graph.getComm(); \
    ArrayView<const GO> STMYGIDS = graph.getRowMap()->getNodeElementList(); \
    size_t STMAX = 0; \
    for (size_t STR=0; STR<graph.getNodeNumRows(); ++STR) { \
      TEST_EQUALITY( graph.getNumEntriesInLocalRow(STR), graph.getNumEntriesInGlobalRow( STMYGIDS[STR] ) ); \
      STMAX = std::max( STMAX, graph.getNumEntriesInLocalRow(STR) ); \
    } \
    TEST_EQUALITY( graph.getNodeMaxNumRowEntries(), STMAX ); \
    global_size_t STGMAX; \
    Teuchos::reduceAll<int,global_size_t>( *STCOMM, Teuchos::REDUCE_MAX, STMAX, outArg(STGMAX) ); \
    TEST_EQUALITY( graph.getGlobalMaxNumRowEntries(), STGMAX ); \
  }


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

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, BadConst, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 10;
    RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,numLocal,comm);
    {
      // bad constructor
      TEST_THROW( GRAPH badgraph(map,INVALID), std::runtime_error ); // allocation hint must be >= 0
    }
    {
      // bad constructor
      ArrayRCP<size_t> hints = arcp<size_t>(numLocal+1);
      std::fill(hints.begin(),hints.end(),1);
      hints[0] = INVALID;
      TEST_THROW( GRAPH badgraph(map,hints.persistingView(0,numLocal+1)), std::runtime_error ); // too many
      TEST_THROW( GRAPH badgraph(map,hints.persistingView(0,numLocal-1)), std::runtime_error ); // too few
      TEST_THROW( GRAPH badgraph(map,hints.persistingView(0,numLocal)),   std::runtime_error ); // too few
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, BadGIDs, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    // create a Map
    const size_t numLocal = 10;
    RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,numLocal,comm);
    {
      Array<GO> gids(1);
      gids[0] = myImageID*numLocal+numLocal;    // off this node, except on the last proc, where it is off the map
      // bad gid on the last node (not in domain map), discovered at fillComplete()
      GRAPH goodgraph(map,1);
      goodgraph.insertGlobalIndices(map->getMinGlobalIndex(),gids());
      TEST_THROW( goodgraph.fillComplete(), std::runtime_error );
    }
    {
      Array<GO> gids(1);
      if (myImageID == numImages-1) {
        gids[0] = myImageID*numLocal+numLocal;   // not in the column map
      }
      else {
        gids[0] = myImageID*numLocal;            // in the column map
      }
      // bad gid on the last node (not in column map)
      // this gid doesn't throw an exception; it is ignored, because the column map acts as a filter
      GRAPH goodgraph(map,map,1);
      goodgraph.insertGlobalIndices(map->getMinGlobalIndex(),gids());
      goodgraph.fillComplete();
      TEST_EQUALITY( goodgraph.getNumEntriesInLocalRow(0),
                     (size_t)(myImageID == numImages-1 ? 0 : 1) );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, BadLIDs, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    // create a Map
    const size_t numLocal = 10;
    RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,numLocal,comm);
    {
      // bad lids on the last node, not in the column map, ignored
      Array<LO> lids(1);
      if (myImageID == numImages-1) {
        lids[0] = numLocal;
      }
      else {
        lids[0] = 0;
      }
      {
        GRAPH diaggraph(map,map,1);
        TEST_EQUALITY(diaggraph.hasColMap(), true);
        // insert on bad row
        TEST_THROW(diaggraph.insertLocalIndices(numLocal,lids()), std::runtime_error);
      }
      {
        GRAPH diaggraph(map,map,1);
        TEST_EQUALITY(diaggraph.hasColMap(), true);
        diaggraph.insertLocalIndices(0,lids());
        TEST_EQUALITY( diaggraph.getNumEntriesInLocalRow(0), (size_t)(myImageID == numImages-1 ? 0 : 1) );
      }
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, EmptyFillComplete, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map with numLocal entries per node
    const size_t numLocal = 10;
    RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,numLocal,comm);
    {
      // create static-profile graph, fill-complete without inserting (and therefore, without allocating)
      GRAPH graph(map,1,StaticProfile);
      graph.fillComplete(DoOptimizeStorage);
    }
    {
      // create dynamic-profile graph, fill-complete without inserting (and therefore, without allocating)
      GRAPH graph(map,1,DynamicProfile);
      graph.fillComplete(DoOptimizeStorage);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, insert_remove_LIDs, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,numLocal,comm);
    {
      Array<LO> lids(1);
      lids[0] = myImageID;
      {
        GRAPH diaggraph(map,map,1);
        TEST_EQUALITY(diaggraph.hasColMap(), true);
        // insert
        LO row = myImageID*numLocal;
        //insert a column-index:
        diaggraph.insertLocalIndices(row, lids());
        TEST_EQUALITY(as<Array_size_type>(diaggraph.getNumEntriesInLocalRow(row)), lids.size())
        //remove the column-index:
        diaggraph.removeLocalIndices(row);
        TEST_EQUALITY(diaggraph.getNumEntriesInLocalRow(row), 0)
        //now inserting the index again, should make the row-length be 1 again...
        diaggraph.insertLocalIndices(row, lids());
        TEST_EQUALITY(as<Array_size_type>(diaggraph.getNumEntriesInLocalRow(row)), lids.size())
      }
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, Bug20100622K, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // test filtering
    if (numImages > 1) {
      // we only need two procs to demonstrate this bug. ignore the others.
      Array<GO> mygids;
      if (myImageID == 0) {
        mygids.push_back(0);
        mygids.push_back(1);
        mygids.push_back(2);
      }
      else if (myImageID == 1) {
        mygids.push_back(2);
        mygids.push_back(3);
      }
      const RCP<const Map<LO,GO> > rmap = rcp(new Map<LO,GO>(INVALID,mygids(),0,comm));
      RCP< CrsGraph<LO,GO> > G = rcp(new CrsGraph<LO,GO>(rmap,2) );
      if (myImageID == 0) {
        G->insertGlobalIndices(0, tuple<GO>(0));
        G->insertGlobalIndices(1, tuple<GO>(0,1));
        G->insertGlobalIndices(2, tuple<GO>(2,1));
      }
      else if (myImageID == 1) {
        G->insertGlobalIndices(2, tuple<GO>(2, 1));
        G->insertGlobalIndices(3, tuple<GO>(3));
      }
      G->fillComplete();
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, WithColMap, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // test filtering
    if (numImages > 1) {
      const size_t numLocal = 1;
      const RCP<const Map<LO,GO> > rmap = createContigMap<LO,GO>(INVALID,numLocal,comm);
      const RCP<const Map<LO,GO> > cmap = createContigMap<LO,GO>(INVALID,numLocal,comm);
      // must allocate enough for all submitted indices, not accounting for filtering.
      RCP< CrsGraph<LO,GO> > G = rcp(new CrsGraph<LO,GO>(rmap,cmap,2,StaticProfile) );
      TEST_EQUALITY_CONST( G->hasColMap(), true );
      const GO myrowind = rmap->getGlobalElement(0);
      TEST_NOTHROW( G->insertGlobalIndices( myrowind, tuple<GO>(myrowind,myrowind+1) ) );
      TEST_NOTHROW( G->fillComplete(DoOptimizeStorage) );
      TEST_EQUALITY( G->getRowMap(), rmap );
      TEST_EQUALITY( G->getColMap(), cmap );
      TEST_EQUALITY( G->getNumEntriesInGlobalRow(myrowind), 1 );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, CopiesAndViews, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    if (numImages < 2) return;
    // create a Map, one row per processor
    const size_t numLocal = 1;
    RCP<const Map<LO,GO> > rmap = createContigMap<LO,GO>(INVALID,numLocal,comm);
    GO myrowind = rmap->getGlobalElement(0);
    // specify the column map to control ordering
    // construct tridiagonal graph
    Array<GO> ginds(0);
    if (myImageID==0) {
      ginds.resize(2);
      ginds[0] = myrowind;
      ginds[1] = myrowind+1;
    }
    else if (myImageID==numImages-1) {
      ginds.resize(2);
      ginds[0] = myrowind-1;
      ginds[1] = myrowind;
    }
    else {
      ginds.resize(3);
      ginds[0] = myrowind-1;
      ginds[1] = myrowind;
      ginds[2] = myrowind+1;
    }
    Array<LO> linds(ginds.size());
    for (typename Array<LO>::size_type i=0; i<linds.size(); ++i) linds[i] = i;
    RCP<const Map<LO,GO,Node> > cmap = rcp( new Map<LO,GO,Node>(INVALID,ginds,0,comm) );
    for (int T=0; T<4; ++T) {
      ProfileType pftype = ( (T & 1) == 1 ) ? StaticProfile : DynamicProfile;
      OptimizeOption os  = ( (T & 2) == 2 ) ? DoOptimizeStorage : DoNotOptimizeStorage;
      GRAPH trigraph(rmap,cmap, ginds.size(),pftype);   // only allocate as much room as necessary
      Array<GO> GCopy(4); Array<LO> LCopy(4);
      ArrayView<const GO> GView; 
      ArrayView<const LO> LView;
      size_t numindices;
      // at this point, there are no global or local indices, but views and copies should succeed
      trigraph.getLocalRowCopy(0,LCopy,numindices);
      trigraph.getLocalRowView(0,LView);
      trigraph.getGlobalRowCopy(myrowind,GCopy,numindices);
      trigraph.getGlobalRowView(myrowind,GView);
      // use multiple inserts: this illustrated an overwrite bug for column-map-specified graphs
      for (Array_size_type j=0; j < ginds.size(); ++j) {
        trigraph.insertGlobalIndices(myrowind,ginds(j,1));
      }
      TEST_EQUALITY( trigraph.getNumEntriesInLocalRow(0), trigraph.getNumAllocatedEntriesInLocalRow(0) ); // test that we only allocated as much room as necessary
      // if static graph, attempt to insert one additional entry on my row and verify that an exception is thrown
      if (pftype == StaticProfile) {
        TEST_THROW( trigraph.insertGlobalIndices(myrowind,tuple<GO>(myrowind)), std::runtime_error );
      }
      trigraph.fillComplete(os);
      // check that inserting global entries throws (inserting local entries is still allowed)
      {
        Array<GO> zero(0);
        TEST_THROW( trigraph.insertGlobalIndices(0,zero()), std::runtime_error );
      }
      // check for throws and no-throws/values
      TEST_THROW( trigraph.getGlobalRowView(myrowind,GView), std::runtime_error );
      TEST_THROW( trigraph.getLocalRowCopy(    0       ,LCopy(0,1),numindices), std::runtime_error );
      TEST_THROW( trigraph.getGlobalRowCopy(myrowind,GCopy(0,1),numindices), std::runtime_error );
      TEST_NOTHROW( trigraph.getLocalRowView(0,LView) );
      TEST_COMPARE_ARRAYS( LView, linds );
      TEST_NOTHROW( trigraph.getLocalRowCopy(0,LCopy,numindices) );
      TEST_COMPARE_ARRAYS( LCopy(0,numindices), linds );
      TEST_NOTHROW( trigraph.getGlobalRowCopy(myrowind,GCopy,numindices) );
      TEST_COMPARE_ARRAYS( GCopy(0,numindices), ginds );
      STD_TESTS(trigraph);
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, WithStaticProfile, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    // create a Map, one row per processor
    const size_t numLocal = 1;
    RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,numLocal,comm);
    {
      // add too many entries to a static graph
      // let node i contribute to row i+1, where node the last node contributes to row 0
      GRAPH diaggraph(map,1,StaticProfile);
      GO grow = myImageID;
      Array<GO> colinds(1);
      colinds[0] = grow;
      TEST_NOTHROW( diaggraph.insertGlobalIndices(grow,colinds()) );
      TEST_THROW( diaggraph.insertGlobalIndices(grow,colinds()), std::runtime_error );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, EmptyGraphAlloc0, LO, GO )
  {
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const size_t STINV = OrdinalTraits<size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = size(*comm);
    const int myImageID = comm->getRank();
    // create the empty graph
    // this one is empty (because of rows distribution) on only one node (demonstrating a previous bug)
    if (numImages > 1) {
      // create a Map
      const size_t numLocal = (myImageID == 1 ? 0 : 1);
      TEST_FOR_EXCEPTION( myImageID != 1 && numLocal != 1, std::logic_error, "numLocal assumed below to be 1.");
      RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,(myImageID == 1 ? 0 : numLocal),comm);
      RCP<RowGraph<LO,GO> > test_row;
      {
        // allocate with no space
        RCP<CrsGraph<LO,GO,Node> > test_crs = rcp(new CrsGraph<LO,GO,Node>(map,0));
        TEST_EQUALITY_CONST( test_crs->getNodeAllocationSize(), STINV ); // invalid, because none are allocated yet
        if (myImageID != 1) test_crs->insertGlobalIndices( map->getMinGlobalIndex(), tuple<GO>(map->getMinGlobalIndex()) );
        test_crs->fillComplete(DoOptimizeStorage);
        TEST_EQUALITY( test_crs->getNodeAllocationSize(), numLocal );
        test_row = test_crs;
      }
      RCP<const Map<LO,GO,Node> > cmap = test_row->getColMap();
      TEST_EQUALITY( cmap->getGlobalNumElements(), (size_t)numImages-1 );
      TEST_EQUALITY( test_row->getGlobalNumRows(), (size_t)numImages-1 );
      TEST_EQUALITY( test_row->getNodeNumRows(), numLocal );
      TEST_EQUALITY( test_row->getGlobalNumCols(), (size_t)numImages-1 );
      TEST_EQUALITY( test_row->getNodeNumCols(), numLocal );
      TEST_EQUALITY( test_row->getIndexBase(), 0 );
      TEST_EQUALITY( test_row->isUpperTriangular(), true );
      TEST_EQUALITY( test_row->isLowerTriangular(), true );
      TEST_EQUALITY( test_row->getGlobalNumDiags(), (size_t)numImages-1 );
      TEST_EQUALITY( test_row->getNodeNumDiags(), numLocal );
      TEST_EQUALITY( test_row->getGlobalNumEntries(), (size_t)numImages-1 );
      TEST_EQUALITY( test_row->getNodeNumEntries(), numLocal );
      TEST_EQUALITY( test_row->getGlobalMaxNumRowEntries(), 1 );
      TEST_EQUALITY( test_row->getNodeMaxNumRowEntries(), numLocal );
      STD_TESTS((*test_row));
    }
    // this one is empty on all nodes because of zero allocation size
    {
      // create a Map
      const size_t numLocal = 10;
      RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,numLocal,comm);
      RCP<RowGraph<LO,GO> > zero;
      {
        // allocate with no space
        RCP<CrsGraph<LO,GO,Node> > zero_crs = rcp(new CrsGraph<LO,GO,Node>(map,0));
        TEST_EQUALITY_CONST( zero_crs->getNodeAllocationSize(), STINV ); // invalid, because none are allocated yet
        zero_crs->fillComplete(DoOptimizeStorage);
        TEST_EQUALITY_CONST( zero_crs->getNodeAllocationSize(), 0     ); // zero, because none were allocated.
        zero = zero_crs;
      }
      RCP<const Map<LO,GO,Node> > cmap = zero->getColMap();
      TEST_EQUALITY( cmap->getGlobalNumElements(), 0 );
      TEST_EQUALITY( zero->getGlobalNumRows(), numImages*numLocal );
      TEST_EQUALITY( zero->getNodeNumRows(), numLocal );
      TEST_EQUALITY( zero->getGlobalNumCols(), 0 );
      TEST_EQUALITY( zero->getNodeNumCols(), 0 );
      TEST_EQUALITY( zero->getIndexBase(), 0 );
      TEST_EQUALITY( zero->isUpperTriangular(), true );
      TEST_EQUALITY( zero->isLowerTriangular(), true );
      TEST_EQUALITY( zero->getGlobalNumDiags(), 0 );
      TEST_EQUALITY( zero->getNodeNumDiags(), 0 );
      TEST_EQUALITY( zero->getGlobalNumEntries(), 0 );
      TEST_EQUALITY( zero->getNodeNumEntries(), 0 );
      TEST_EQUALITY( zero->getGlobalMaxNumRowEntries(), 0 );
      TEST_EQUALITY( zero->getNodeMaxNumRowEntries(), 0 );
      STD_TESTS((*zero));
    }
    // this one is empty on all nodes because of zero rows
    {
      // create a Map
      const size_t numLocal = 0;
      RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,numLocal,comm);
      RCP<RowGraph<LO,GO> > zero;
      {
        // allocate with no space
        RCP<CrsGraph<LO,GO,Node> > zero_crs = rcp(new CrsGraph<LO,GO,Node>(map,0));
        TEST_EQUALITY_CONST( zero_crs->getNodeAllocationSize(), STINV ); // invalid, because none are allocated yet
        zero_crs->fillComplete(DoOptimizeStorage);
        TEST_EQUALITY_CONST( zero_crs->getNodeAllocationSize(), 0     ); // zero, because none were allocated.
        zero = zero_crs;
      }
      RCP<const Map<LO,GO,Node> > cmap = zero->getColMap();
      TEST_EQUALITY( cmap->getGlobalNumElements(), 0 );
      TEST_EQUALITY( zero->getGlobalNumRows(), numImages*numLocal );
      TEST_EQUALITY( zero->getNodeNumRows(), numLocal );
      TEST_EQUALITY( zero->getGlobalNumCols(), 0 );
      TEST_EQUALITY( zero->getNodeNumCols(), 0 );
      TEST_EQUALITY( zero->getIndexBase(), 0 );
      TEST_EQUALITY( zero->isUpperTriangular(), true );
      TEST_EQUALITY( zero->isLowerTriangular(), true );
      TEST_EQUALITY( zero->getGlobalNumDiags(), 0 );
      TEST_EQUALITY( zero->getNodeNumDiags(), 0 );
      TEST_EQUALITY( zero->getGlobalNumEntries(), 0 );
      TEST_EQUALITY( zero->getNodeNumEntries(), 0 );
      TEST_EQUALITY( zero->getGlobalMaxNumRowEntries(), 0 );
      TEST_EQUALITY( zero->getNodeMaxNumRowEntries(), 0 );
      STD_TESTS((*zero));
    }

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, EmptyGraphAlloc1, LO, GO )
  {
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const size_t STINV = OrdinalTraits<size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    int numImages = size(*comm);
    // create a Map
    const size_t numLocal = 10;
    RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,numLocal,comm);
    // create the empty graph
    RCP<RowGraph<LO,GO> > zero;
    {
      // allocated with space for one entry per row
      RCP<CrsGraph<LO,GO,Node> > zero_crs = rcp(new CrsGraph<LO,GO,Node>(map,1));
      TEST_EQUALITY( zero_crs->getNodeAllocationSize(), STINV ); // zero, because none are allocated yet
      zero_crs->fillComplete(DoOptimizeStorage);
      zero = zero_crs;
    }
    RCP<const Map<LO,GO,Node> > cmap = zero->getColMap();
    TEST_EQUALITY( cmap->getGlobalNumElements(), 0 );
    TEST_EQUALITY( zero->getGlobalNumRows(), numImages*numLocal );
    TEST_EQUALITY( zero->getNodeNumRows(), numLocal );
    TEST_EQUALITY( zero->getGlobalNumCols(), 0 );
    TEST_EQUALITY( zero->getNodeNumCols(), 0 );
    TEST_EQUALITY( zero->getIndexBase(), 0 );
    TEST_EQUALITY( zero->isUpperTriangular(), true );
    TEST_EQUALITY( zero->isLowerTriangular(), true );
    TEST_EQUALITY( zero->getGlobalNumDiags(), 0 );
    TEST_EQUALITY( zero->getNodeNumDiags(), 0 );
    TEST_EQUALITY( zero->getGlobalNumEntries(), 0 );
    TEST_EQUALITY( zero->getNodeNumEntries(), 0 );
    TEST_EQUALITY( zero->getGlobalMaxNumRowEntries(), 0 );
    TEST_EQUALITY( zero->getNodeMaxNumRowEntries(), 0 );
    STD_TESTS((*zero));

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, DottedDiag, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a Map, three rows per processor
    const size_t numLocal = 3;
    RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,numLocal,comm);
    GO mymiddle = map->getGlobalElement(1);  // get my middle row
    for (int T=0; T<4; ++T) {
      ProfileType pftype = ( (T & 1) == 1 ) ? StaticProfile : DynamicProfile;
      OptimizeOption os  = ( (T & 2) == 2 ) ? DoOptimizeStorage : DoNotOptimizeStorage;
      {
        // create a diagonal graph, but where only my middle row has an entry
        ArrayRCP<size_t> toalloc = arcpClone<size_t>( tuple<size_t>(0,1,0) );
        GRAPH ddgraph(map,toalloc,pftype);
        ddgraph.insertGlobalIndices(mymiddle, tuple<GO>(mymiddle));
        // before globalAssemble(), there should be one local entry on middle, none on the others
        ArrayView<const GO> myrow_gbl;
        ddgraph.getGlobalRowView(mymiddle-1,myrow_gbl); TEST_EQUALITY( myrow_gbl.size(), 0 );
        ddgraph.getGlobalRowView(mymiddle  ,myrow_gbl); TEST_COMPARE_ARRAYS( myrow_gbl, tuple<GO>(mymiddle) );
        ddgraph.getGlobalRowView(mymiddle+1,myrow_gbl); TEST_EQUALITY( myrow_gbl.size(), 0 );
        if (pftype == StaticProfile) { // no room for more, on any row
          TEST_THROW( ddgraph.insertGlobalIndices(mymiddle-1,tuple<GO>(mymiddle)), std::runtime_error );
          TEST_THROW( ddgraph.insertGlobalIndices(mymiddle  ,tuple<GO>(mymiddle)), std::runtime_error );
          TEST_THROW( ddgraph.insertGlobalIndices(mymiddle+1,tuple<GO>(mymiddle)), std::runtime_error );
        }
        ddgraph.fillComplete(os);
        // after fillComplete(), there should be a single entry on my middle, corresponding to the diagonal, none on the others
        ArrayView<const LO> myrow_lcl;
        TEST_EQUALITY_CONST( ddgraph.getNumEntriesInLocalRow(0), 0 );
        TEST_EQUALITY_CONST( ddgraph.getNumEntriesInLocalRow(2), 0 );
        ddgraph.getLocalRowView(1,myrow_lcl);
        TEST_EQUALITY_CONST( myrow_lcl.size(), 1 );
        if (myrow_lcl.size() == 1) {
          TEST_EQUALITY( ddgraph.getColMap()->getGlobalElement(myrow_lcl[0]), mymiddle );
        }
        // also, the row map and column map should be equivalent
        TEST_EQUALITY( ddgraph.getGlobalNumCols(), (global_size_t)(numImages) );
        TEST_EQUALITY( ddgraph.getGlobalNumRows(), (global_size_t)(3*numImages) );
        TEST_EQUALITY( ddgraph.getGlobalNumDiags(), (global_size_t)numImages );
        TEST_EQUALITY_CONST( ddgraph.getNodeNumDiags(), 1 );
        STD_TESTS(ddgraph);
      }
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, NonLocals, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    // create a Map, one row per processor
    const size_t numLocal = 1;
    RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,numLocal,comm);
    GO myrowind = map->getGlobalElement(0);
    for (int T=0; T<4; ++T) {
      ProfileType pftype = ( (T & 1) == 1 ) ? StaticProfile : DynamicProfile;
      OptimizeOption os  = ( (T & 2) == 2 ) ? DoOptimizeStorage : DoNotOptimizeStorage;
      {
        // create a diagonal graph, where the graph entries are contributed by a single off-node contribution, no filtering
        // let node i contribute to row i+1, where node the last node contributes to row 0
        GRAPH diaggraph(map,1,pftype);
        GO grow = myImageID+1;
        if (grow == numImages) {
          grow = 0;
        }
        diaggraph.insertGlobalIndices(grow, tuple<GO>(grow));
        // before globalAssemble(), there should be no local entries if numImages > 1
        ArrayView<const GO> myrow_gbl; 
        diaggraph.getGlobalRowView(myrowind, myrow_gbl);
        TEST_EQUALITY( myrow_gbl.size(), (numImages == 1 ? 1 : 0) );
        diaggraph.globalAssemble();   // after globalAssemble(), there should be one local entry per row, corresponding to the diagonal
        diaggraph.getGlobalRowView(myrowind,myrow_gbl);
        TEST_COMPARE_ARRAYS( myrow_gbl, tuple<GO>(myrowind) );
        if (pftype == StaticProfile) { // no room for more
          TEST_THROW( diaggraph.insertGlobalIndices(myrowind,tuple<GO>(myrowind)), std::runtime_error );
        }
        diaggraph.fillComplete(os);
        // after fillComplete(), there should be a single entry on my row, corresponding to the diagonal
        ArrayView<const LO> myrow_lcl; 
        diaggraph.getLocalRowView(0, myrow_lcl);
        TEST_EQUALITY_CONST( myrow_lcl.size(), 1 );
        if (myrow_lcl.size() == 1) {
          TEST_EQUALITY( diaggraph.getColMap()->getGlobalElement(myrow_lcl[0]), myrowind );
        }
        // also, the row map and column map should be equivalent
        TEST_EQUALITY_CONST( diaggraph.getRowMap()->isSameAs(*diaggraph.getColMap()), true );
        TEST_EQUALITY( diaggraph.getGlobalNumDiags(), (global_size_t)numImages );
        TEST_EQUALITY_CONST( diaggraph.getNodeNumDiags(), 1 );
        TEST_EQUALITY_CONST( diaggraph.isUpperTriangular(), true );
        TEST_EQUALITY_CONST( diaggraph.isLowerTriangular(), true );
        STD_TESTS(diaggraph);
      }
      {
        // create a next-door-neighbor graph (tridiagonal plus corners), where the graph entries are contributed by single off-node contribution, no filtering
        // let node i add the contributions for column i of the graph: (i-1,i), (i,i), (i+1,i)
        // allocate only as much space as we need
        // some hacking here to support this test when numImages == 1 or 2
        GRAPH ngraph(map,3,pftype);
        Array<GO> grows(3);
        grows[0] = (numImages+myImageID-1) % numImages;   // my left neighbor
        grows[1] = (numImages+myImageID  ) % numImages;   // myself
        grows[2] = (numImages+myImageID+1) % numImages;   // my right neighbor
        ngraph.insertGlobalIndices(grows[0],tuple<GO>(myImageID)); // ^^^^^^^^^^^^^^^^^^^^^^^
        ngraph.insertGlobalIndices(grows[1],tuple<GO>(myImageID)); // add me to the graph for my neighbors
        ngraph.insertGlobalIndices(grows[2],tuple<GO>(myImageID)); // vvvvvvvvvvvvvvvvvvvvvvv
        // before globalAssemble(), there should be a single local entry on parallel runs, three on serial runs
        ArrayView<const GO> myrow_gbl; 
        ngraph.getGlobalRowView(myrowind, myrow_gbl);
        TEST_EQUALITY_CONST( myrow_gbl.size(), (numImages == 1 ? 3 : 1) );
        ngraph.globalAssemble();    // after globalAssemble(), storage should be maxed out
        TEST_EQUALITY( ngraph.getNumEntriesInLocalRow(0), ngraph.getNumAllocatedEntriesInLocalRow(0) );
        if (pftype == StaticProfile) {
          TEST_THROW( ngraph.insertGlobalIndices(myImageID,tuple<GO>(myImageID)), std::runtime_error );  // adding an addition entry under static allocation should fail
        }
        ngraph.fillComplete(os);
        // after fillComplete(), there should be entries for me and my neighbors on my row
        ArrayView<const LO> myrow_lcl; 
        ngraph.getLocalRowView(0, myrow_lcl);
        {
          // check indices on my row
          typename Array<GO>::iterator glast;
          sort(grows.begin(),grows.end());
          glast = unique(grows.begin(),grows.end());
          size_t numunique = glast - grows.begin();
          // test the test: numunique == min(numImages,3)
          TEST_EQUALITY( numunique, (size_t)min(numImages,3) );
          TEST_EQUALITY_CONST( (size_t)myrow_lcl.size(), numunique );
          if ((size_t)myrow_lcl.size() == numunique) {
            size_t numinds;
            Array<GO> inds(numunique+1);
            TEST_THROW(   ngraph.getGlobalRowCopy(myrowind,inds(0,numunique-1), numinds), std::runtime_error );
            TEST_NOTHROW( ngraph.getGlobalRowCopy(myrowind,inds(0,numunique), numinds) );
            TEST_NOTHROW( ngraph.getGlobalRowCopy(myrowind,inds(), numinds) );
            sort(inds.begin(), inds.begin()+numinds);
            TEST_COMPARE_ARRAYS( inds(0,numinds), grows(0,numunique) );
          }
        }
        TEST_EQUALITY_CONST( ngraph.getRowMap()->isSameAs(*ngraph.getColMap()), (numImages==1 ? true : false) );
        TEST_EQUALITY( ngraph.getGlobalNumDiags(), (global_size_t)numImages );
        TEST_EQUALITY( ngraph.getNodeNumDiags(), 1 );
        STD_TESTS(ngraph);
      }
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, Describable, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    // create Map
    RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,3,comm);
    {
      GRAPH graph(map,1,StaticProfile);
      // test labeling
      const string lbl("graphA");
      string desc1 = graph.description();
      out << desc1 << endl;
      graph.setObjectLabel(lbl);
      string desc2 = graph.description();
      out << desc2 << endl;
      TEST_EQUALITY( graph.getObjectLabel(), lbl );
    }
    {
      GRAPH graph(map,1,StaticProfile);
      // test describing at different verbosity levels
      if (myImageID==0) out << "Describing with verbosity VERB_DEFAULT..." << endl;
      graph.describe(out);
      comm->barrier();
      comm->barrier();
      if (myImageID==0) out << "Describing with verbosity VERB_NONE..." << endl;
      graph.describe(out,VERB_NONE);
      comm->barrier();
      comm->barrier();
      if (myImageID==0) out << "Describing with verbosity VERB_LOW..." << endl;
      graph.describe(out,VERB_LOW);
      comm->barrier();
      comm->barrier();
      if (myImageID==0) out << "Describing with verbosity VERB_MEDIUM..." << endl;
      graph.describe(out,VERB_MEDIUM);
      comm->barrier();
      comm->barrier();
      if (myImageID==0) out << "Describing with verbosity VERB_HIGH..." << endl;
      graph.describe(out,VERB_HIGH);
      comm->barrier();
      comm->barrier();
      if (myImageID==0) out << "Describing with verbosity VERB_EXTREME..." << endl;
      graph.describe(out,VERB_EXTREME);
      comm->barrier();
      comm->barrier();
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, ActiveFill, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create Map
    RCP<const Map<LO,GO> > map = createContigMap<LO,GO>(INVALID,1,comm);
    {
      GRAPH graph(map,map,0,DynamicProfile);
      TEST_EQUALITY_CONST( graph.isFillActive(),   true );
      TEST_EQUALITY_CONST( graph.isFillComplete(), false );
      graph.insertLocalIndices( 0, tuple<LO>(0) );
      //
      graph.fillComplete(DoNotOptimizeStorage);
      TEST_EQUALITY_CONST( graph.isFillActive(),   false );
      TEST_EQUALITY_CONST( graph.isFillComplete(), true );
      TEST_THROW( graph.insertLocalIndices( 0, tuple<LO>(0) ), std::runtime_error );
      TEST_THROW( graph.removeLocalIndices( 0 ),               std::runtime_error );
      TEST_THROW( graph.globalAssemble(),                      std::runtime_error );
      TEST_THROW( graph.fillComplete(),                        std::runtime_error );
    }
    {
      GRAPH graph(map,map,0,DynamicProfile);
      TEST_EQUALITY_CONST( graph.isFillActive(),   true );
      TEST_EQUALITY_CONST( graph.isFillComplete(), false );
      graph.insertLocalIndices( 0, tuple<LO>(0) );
      //
      graph.fillComplete(DoNotOptimizeStorage);
      TEST_EQUALITY_CONST( graph.isFillActive(),   false );
      TEST_EQUALITY_CONST( graph.isFillComplete(), true );
      //
      graph.resumeFill();
      TEST_EQUALITY_CONST( graph.isFillActive(),   true );
      TEST_EQUALITY_CONST( graph.isFillComplete(), false );
      TEST_NOTHROW( graph.insertLocalIndices( 0, tuple<LO>(0) ) );
      //
      TEST_NOTHROW( graph.fillComplete()                        );
      TEST_EQUALITY_CONST( graph.isFillActive(),   false );
      TEST_EQUALITY_CONST( graph.isFillComplete(), true );
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, Typedefs, LO, GO )
  {
    typedef CrsGraph<LO,GO,Node> GRAPH;
    typedef typename GRAPH::local_ordinal_type  local_ordinal_type;
    typedef typename GRAPH::global_ordinal_type global_ordinal_type;
    typedef typename GRAPH::node_type           node_type;
    TEST_EQUALITY_CONST( (is_same< local_ordinal_type  , LO  >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< global_ordinal_type , GO  >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< node_type           , Node>::value) == true, true );
    typedef RowGraph<LO,GO,Node> RGRAPH;
    typedef typename RGRAPH::local_ordinal_type  rgraph_local_ordinal_type;
    typedef typename RGRAPH::global_ordinal_type rgraph_global_ordinal_type;
    typedef typename RGRAPH::node_type           rgraph_node_type;
    TEST_EQUALITY_CONST( (is_same< rgraph_local_ordinal_type  , LO  >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< rgraph_global_ordinal_type , GO  >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< rgraph_node_type           , Node>::value) == true, true );
  }

  // 
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_LO_GO( LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, EmptyGraphAlloc0, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, EmptyGraphAlloc1, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, BadConst  , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, BadGIDs   , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, BadLIDs   , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, insert_remove_LIDs   , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, NonLocals , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, DottedDiag , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, WithStaticProfile , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, CopiesAndViews, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, WithColMap,     LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, Describable   , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, EmptyFillComplete, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, Typedefs      , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, Bug20100622K  , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, ActiveFill    , LO, GO )

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
