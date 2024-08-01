// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_TestingUtilities.hpp"

namespace { // (anonymous)

  using Tpetra::TestingUtilities::arcp_from_view;
  using Teuchos::arcp;
  using Teuchos::arcpClone;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::Comm;
  using Teuchos::null;
  using Teuchos::outArg;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  using Teuchos::tuple;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using std::endl;
  using std::max;
  using std::min;
  using GST = Tpetra::global_size_t;

  double errorTolSlack = 1e+1;
  std::string filedir;

  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();

#define STD_TESTS(graph) \
  { \
    auto STCOMM = graph.getComm(); \
    auto STMYGIDS = graph.getRowMap()->getLocalElementList(); \
    size_t STMAX = 0; \
    for (size_t STR = 0; STR < graph.getLocalNumRows(); ++STR) { \
      TEST_EQUALITY( graph.getNumEntriesInLocalRow (STR), graph.getNumEntriesInGlobalRow (STMYGIDS[STR]) ); \
      STMAX = std::max (STMAX, graph.getNumEntriesInLocalRow(STR)); \
    } \
    TEST_EQUALITY( graph.getLocalMaxNumRowEntries(), STMAX ); \
    GST STGMAX; \
    reduceAll<int, GST> (*STCOMM, Teuchos::REDUCE_MAX, STMAX, outArg (STGMAX)); \
    TEST_EQUALITY( graph.getGlobalMaxNumRowEntries(), STGMAX ); \
  }


  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption(
        "filedir",&filedir,"Directory of expected input files.");
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, BadConst, LO, GO , Node )
  {
    using Teuchos::REDUCE_MIN;
    using map_type = Tpetra::Map<LO, GO, Node>;
    using graph_type = Tpetra::CrsGraph<LO, GO, Node>;
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();

    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?

    // Create a contiguous Map to use as the graph's row Map.
    const size_t numLocal = 10;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocal, indexBase, comm));

    RCP<graph_type> badgraph;

    // CrsGraph's constructor should throw on invalid allocation hint.
    const size_t allocHint = Teuchos::OrdinalTraits<size_t>::invalid ();
    TEST_THROW( badgraph = rcp (new graph_type (map, allocHint)),
                std::invalid_argument );

    // Make sure that the test passed on all processes.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // CrsGraph's constructor should throw if input allocation hint
    // has the wrong length, as it does in this case.
    ArrayRCP<size_t> hints = arcp<size_t> (numLocal + 1);
    std::fill (hints.begin (), hints.end (), static_cast<size_t> (1));
    TEST_THROW( badgraph = rcp (new graph_type (map, hints (0, numLocal+1))),
                std::invalid_argument ); // too many entries

    // Make sure that the test passed on all processes.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // CrsGraph's constructor should throw if input allocation hint
    // has the wrong length, as it does in this case.
    TEST_THROW( badgraph = rcp (new graph_type (map, hints (0, numLocal-1))),
                std::invalid_argument ); // too few entries

    // Make sure that the test passed on all processes.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // mfh 13 Aug 2014: It's not reasonable to ask the constructor to
    // check all the values.  That might be reasonable for a debug
    // build, but not for a release build.

    // hints[0] = Teuchos::OrdinalTraits<size_t>::invalid ();
    // TEST_THROW( badgraph = rcp (new graph_type (map, hints (0, numLocal))),
    //             std::invalid_argument ); // invalid value

    // // Make sure that the test passed on all processes.
    // lclSuccess = success ? 1 : 0;
    // gblSuccess = 0;
    // reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    // TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, BadGIDs, LO, GO , Node )
  {
    const bool debug = Tpetra::Details::Behavior::debug("CrsGraph");
    if (debug) {
      using GRAPH = Tpetra::CrsGraph<LO, GO, Node>;
      using map_type = Tpetra::Map<LO, GO, Node>;

      // get a comm
      RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
      const int myRank = comm->getRank();
      // create a Map
      const size_t numLocal = 10;
      const GO indexBase = 0;
      RCP<const map_type> map =
        rcp (new map_type (INVALID, numLocal, indexBase, comm));
      {
        Array<GO> gids(1);
        // This GID is off process but still in the domain Map on
        // every process but the last.  On the last process, this GID
        // is not in the domain Map on ANY process.  This makes
        // inserting it an error.
        gids[0] = myRank*numLocal+numLocal;
        // In debug mode, CrsGraph::fillComplete discovers on the last
        // process that this GID is bad, because it's not in the
        // domain Map.
        GRAPH goodgraph(map,1);
        goodgraph.insertGlobalIndices(map->getMinGlobalIndex(),gids());
        TEST_THROW( goodgraph.fillComplete(), std::runtime_error );
      }
      // All procs fail if any process fails
      int globalSuccess_int = -1;
      reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, ExcessAllocation, LO, GO , Node )
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    typedef Tpetra::CrsGraph<LO, GO, Node>  GRPH;
    using map_type = Tpetra::Map<LO, GO, Node>;

    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    // create a Map with numLocal entries per node
    const int numLocal = 10;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocal, indexBase, comm));

    // send in a parameterlist, check the defaults
    RCP<ParameterList> params = parameterList();
    // create static-profile graph, fill-complete without inserting (and therefore, without allocating)
    GRPH graph (map, 3);
    for (GO i = map->getMinGlobalIndex(); i <= map->getMaxGlobalIndex(); ++i) {
      graph.insertGlobalIndices (i, tuple<GO> (i));
    }
    params->set("Optimize Storage",false);
    graph.fillComplete(params);
    TEST_EQUALITY(graph.getLocalNumEntries(), (size_t)numLocal);
    TEST_EQUALITY_CONST(graph.isStorageOptimized(), false);
    //
    graph.resumeFill();
    for (int i=0; i < numLocal; ++i) graph.removeLocalIndices(i);
    params->set("Optimize Storage",true);
    graph.fillComplete(params);
    TEST_EQUALITY_CONST(params->get<bool>("Optimize Storage"), true);
    TEST_EQUALITY(graph.getLocalNumEntries(), 0);
    TEST_EQUALITY_CONST(graph.isStorageOptimized(), true);

    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));

    if (gblSuccess == 1) {
      out << "Succeeded on all processes!" << endl;
    } else {
      out << "FAILED on at least one process!" << endl;
    }
    TEST_EQUALITY_CONST(gblSuccess, 1);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, insert_remove_LIDs, LO, GO , Node )
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using GRAPH = Tpetra::CrsGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int myRank = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocal, indexBase, comm));

    Array<LO> lids(1);
    lids[0] = myRank;
    GRAPH diaggraph(map,map,1);
    TEST_EQUALITY(diaggraph.hasColMap(), true);
    const LO row = 0;
    // insert
    diaggraph.insertLocalIndices(row, lids());
    TEST_EQUALITY( static_cast<size_t> (diaggraph.getNumEntriesInLocalRow (row)),
                   static_cast<size_t> (lids.size ()) );
    // remove the row
    diaggraph.removeLocalIndices(row);
    TEST_EQUALITY(diaggraph.getNumEntriesInLocalRow(row), 0)
    // now inserting the index again, should make the row-length be 1 again...
    diaggraph.insertLocalIndices(row, lids());
      TEST_EQUALITY(static_cast<size_t> (diaggraph.getNumEntriesInLocalRow (row)),
                    static_cast<size_t> (lids.size ()) );

    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));

    if (gblSuccess == 1) {
      out << "Succeeded on all processes!" << endl;
    } else {
      out << "FAILED on at least one process!" << endl;
    }
    TEST_EQUALITY_CONST(gblSuccess, 1);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, CopiesAndViews, LO, GO , Node )
  {
    using GRAPH = Tpetra::CrsGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int myRank = comm->getRank();
    const int numProcs = comm->getSize();
    if (numProcs < 2) return;
    // create a Map, one row per processor
    const size_t numLocal = 1;
    RCP<const map_type> rmap =
      rcp (new map_type (INVALID, numLocal, 0, comm));
    GO myrowind = rmap->getGlobalElement(0);
    // specify the column map to control ordering
    // construct tridiagonal graph
    Array<GO> ginds(0);
    if (myRank==0) {
      ginds.resize(2);
      ginds[0] = myrowind;
      ginds[1] = myrowind+1;
    }
    else if (myRank==numProcs-1) {
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
    RCP<const map_type> cmap =
      rcp (new map_type (INVALID, ginds, 0, comm));

    RCP<ParameterList> params = parameterList();
    for (int T=0; T<4; ++T) {
      if ( (T & 1) != 1 ) continue;
      params->set("Optimize Storage",((T & 2) == 2));
      GRAPH trigraph(rmap,cmap, ginds.size());   // only allocate as much room as necessary
      size_t numindices;
      {

        typename GRAPH::global_inds_host_view_type GView;
        typename GRAPH::local_inds_host_view_type LView;
        typename GRAPH::nonconst_global_inds_host_view_type GCopy("gcopy",4);
        typename GRAPH::nonconst_local_inds_host_view_type LCopy("lcopy",4);
        // at this point, there are no global or local indices, but views and copies should succeed
        TEST_NOTHROW(trigraph.getLocalRowCopy(0,LCopy,numindices));
        TEST_NOTHROW(trigraph.getLocalRowView(0,LView));
        TEST_NOTHROW(trigraph.getGlobalRowCopy(myrowind,GCopy,numindices));
        TEST_NOTHROW(trigraph.getGlobalRowView(myrowind,GView));
      }
        // use multiple inserts: this illustrated an overwrite bug for column-map-specified graphs
      typedef typename Teuchos::ArrayView<const GO>::size_type size_type;
      for (size_type j=0; j < ginds.size(); ++j) {
        trigraph.insertGlobalIndices(myrowind,ginds(j,1));
      }
      TEST_EQUALITY( trigraph.getNumEntriesInLocalRow(0), trigraph.getNumAllocatedEntriesInLocalRow(0) ); // test that we only allocated as much room as necessary
      // Attempt to insert one additional entry
      // in my row that is not already in the row, and verify that it
      // throws an exception.
      TEST_THROW( trigraph.insertGlobalIndices(myrowind,tuple<GO>(myrowind+2)), std::runtime_error );
      trigraph.fillComplete(params);
      // check that inserting global entries throws (inserting local entries is still allowed)
      {
        Array<GO> zero(0);
        TEST_THROW( trigraph.insertGlobalIndices(0,zero()), std::runtime_error );
      }
      {
        // check for throws and no-throws/values
        typename GRAPH::global_inds_host_view_type GView;
        typename GRAPH::local_inds_host_view_type LView;
        typename GRAPH::nonconst_global_inds_host_view_type GCopy_short("gshort",1), GCopy("gcopy",4);
        typename GRAPH::nonconst_local_inds_host_view_type LCopy_short("lshort",1), LCopy("lcopy",4);
        TEST_THROW( trigraph.getGlobalRowView(myrowind,GView), std::runtime_error );
        TEST_THROW( trigraph.getLocalRowCopy(0, LCopy_short,numindices), std::runtime_error );
        TEST_THROW( trigraph.getGlobalRowCopy(myrowind,GCopy_short,numindices), std::runtime_error );
        TEST_NOTHROW( trigraph.getLocalRowView(0,LView) );
        TEST_COMPARE_ARRAYS( LView, linds );
        TEST_NOTHROW( trigraph.getLocalRowCopy(0,LCopy,numindices) );
        TEST_COMPARE_ARRAYS( arcp_from_view(LCopy,numindices), linds );
        TEST_NOTHROW( trigraph.getGlobalRowCopy(myrowind,GCopy,numindices) );
        TEST_COMPARE_ARRAYS( arcp_from_view(GCopy,numindices), ginds );
      }
      STD_TESTS(trigraph);
      
      // All procs fail if any node fails
      int globalSuccess_int = -1;
      reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, WithStaticProfile, LO, GO , Node )
  {
    using GRAPH = Tpetra::CrsGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int myRank = comm->getRank();
    // create a Map, one row per processor
    const size_t numLocal = 1;
    RCP<const map_type> map = rcp (new map_type (INVALID, numLocal, 0, comm));

    // add too many entries to a static graph
    // let node i contribute to row i+1, where node the last node contributes to row 0
    GRAPH diaggraph(map,1);
    GO grow = myRank;
    Array<GO> colinds(1);
    colinds[0] = grow;
    TEST_NOTHROW( diaggraph.insertGlobalIndices(grow,colinds()) );
    TEST_THROW( diaggraph.insertGlobalIndices(grow, Teuchos::tuple<GO> (grow+1)), std::runtime_error );

    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, EmptyGraphAlloc0, LO, GO , Node )
  {
    using crs_graph_type = Tpetra::CrsGraph<LO, GO, Node>;
    using row_graph_type = Tpetra::RowGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    out << "CrsGrap EmptyGraphAlloc0 test" << endl;
    Teuchos::OSTab tab0 (out);

    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int numProcs = size(*comm);
    const int myRank = comm->getRank();
    // create the empty graph
    // this one is empty (because of rows distribution) on only one node (demonstrating a previous bug)
    if (numProcs > 1) {
      // create a Map
      const size_t numLocal = (myRank == 1 ? 0 : 1);
      TEUCHOS_TEST_FOR_EXCEPTION( myRank != 1 && numLocal != 1, std::logic_error, "numLocal assumed below to be 1.");
      RCP<const map_type> map =
        rcp (new map_type (INVALID, (myRank == 1 ? 0 : numLocal), 0, comm));

      Teuchos::OSTab tab1 (out);

      RCP<row_graph_type> test_row;
      // allocate
      RCP<crs_graph_type> test_crs = rcp (new crs_graph_type (map, 1));
      // invalid, because none are allocated yet
      TEST_EQUALITY_CONST( test_crs->getLocalAllocationSize(), STINV );
      if (myRank != 1) {
        test_crs->insertGlobalIndices (map->getMinGlobalIndex (),
                                       tuple<GO> (map->getMinGlobalIndex ()));
      }
      test_crs->fillComplete ();
      TEST_EQUALITY( test_crs->getLocalAllocationSize(), numLocal );
      test_row = test_crs;
      RCP<const map_type> cmap = test_row->getColMap();
      TEST_EQUALITY( cmap->getGlobalNumElements(), (size_t)numProcs-1 );
      TEST_EQUALITY( test_row->getGlobalNumRows(), (size_t)numProcs-1 );
      TEST_EQUALITY( test_row->getLocalNumRows(), numLocal );
      TEST_EQUALITY( test_row->getGlobalNumCols(), (size_t)numProcs-1 );
      TEST_EQUALITY( test_row->getLocalNumCols(), numLocal );
      TEST_EQUALITY( test_row->getIndexBase(), 0 );

      TEST_EQUALITY( test_row->getGlobalNumEntries(), (size_t)numProcs-1 );
      TEST_EQUALITY( test_row->getLocalNumEntries(), numLocal );
      TEST_EQUALITY( test_row->getGlobalMaxNumRowEntries(), 1 );
      TEST_EQUALITY( test_row->getLocalMaxNumRowEntries(), numLocal );
      STD_TESTS((*test_row));
    }

    // this one is empty on all nodes because of zero allocation size
    {
      // create a Map
      const size_t numLocal = 10;
      RCP<const map_type> map =
        rcp (new map_type (INVALID, numLocal, 0, comm));

      {
        RCP<row_graph_type> zero;
        // allocate with no space
        RCP<crs_graph_type> zero_crs = rcp (new crs_graph_type (map, 0));
        // invalid, because none are allocated yet
        TEST_EQUALITY_CONST( zero_crs->getLocalAllocationSize(), STINV );
        zero_crs->fillComplete ();
        // zero, because none were allocated.
        TEST_EQUALITY_CONST( zero_crs->getLocalAllocationSize(), 0 );
        zero = zero_crs;
        RCP<const map_type> cmap = zero->getColMap();
        TEST_EQUALITY( cmap->getGlobalNumElements(), 0 );
        TEST_EQUALITY( zero->getGlobalNumRows(), numProcs*numLocal );
        TEST_EQUALITY( zero->getLocalNumRows(), numLocal );
        TEST_EQUALITY( zero->getGlobalNumCols(), numProcs*numLocal );
        TEST_EQUALITY( zero->getLocalNumCols(), 0 );
        TEST_EQUALITY( zero->getIndexBase(), 0 );

        TEST_EQUALITY( zero->getGlobalNumEntries(), 0 );
        TEST_EQUALITY( zero->getLocalNumEntries(), 0 );
        TEST_EQUALITY( zero->getGlobalMaxNumRowEntries(), 0 );
        TEST_EQUALITY( zero->getLocalMaxNumRowEntries(), 0 );
        STD_TESTS((*zero));
      }
    }
    // this one is empty on all nodes because of zero rows
    {
      // create a Map
      const size_t numLocal = 0;
      RCP<const map_type> map = rcp (new map_type (INVALID, numLocal, 0, comm));

      RCP<row_graph_type> zero;
      // allocate with no space
      RCP<crs_graph_type> zero_crs = rcp (new crs_graph_type (map, 0));
      // invalid, because none are allocated yet
        TEST_EQUALITY_CONST( zero_crs->getLocalAllocationSize(), STINV );
      zero_crs->fillComplete ();
      // zero, because none were allocated.
      TEST_EQUALITY_CONST( zero_crs->getLocalAllocationSize(), 0 );
      zero = zero_crs;
      RCP<const map_type> cmap = zero->getColMap();
      TEST_EQUALITY( cmap->getGlobalNumElements(), 0 );
      TEST_EQUALITY( zero->getGlobalNumRows(), numProcs*numLocal );
      TEST_EQUALITY( zero->getLocalNumRows(), numLocal );
      TEST_EQUALITY( zero->getGlobalNumCols(), numProcs*numLocal );
      TEST_EQUALITY( zero->getLocalNumCols(), 0 );
      TEST_EQUALITY( zero->getIndexBase(), 0 );

      TEST_EQUALITY( zero->getGlobalNumEntries(), 0 );
      TEST_EQUALITY( zero->getLocalNumEntries(), 0 );
      TEST_EQUALITY( zero->getGlobalMaxNumRowEntries(), 0 );
      TEST_EQUALITY( zero->getLocalMaxNumRowEntries(), 0 );
      STD_TESTS((*zero));
    }

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, EmptyGraphAlloc1, LO, GO , Node )
  {
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    using map_type = Tpetra::Map<LO, GO, Node>;

    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    int numProcs = size(*comm);
    // create a Map
    const size_t numLocal = 10;
    RCP<const map_type> map = rcp (new map_type (INVALID, numLocal, 0, comm));
    // create the empty graph
    RCP<Tpetra::RowGraph<LO,GO,Node> > zero;

    // allocated with space for one entry per row
    RCP<graph_type> zero_crs = rcp (new graph_type (map,1));
    TEST_EQUALITY( zero_crs->getLocalAllocationSize(), STINV ); // zero, because none are allocated yet

    zero_crs->fillComplete ();
    zero = zero_crs;
    RCP<const map_type> cmap = zero->getColMap();
    TEST_EQUALITY( cmap->getGlobalNumElements(), 0 );
    TEST_EQUALITY( zero->getGlobalNumRows(), numProcs*numLocal );
    TEST_EQUALITY( zero->getLocalNumRows(), numLocal );
    TEST_EQUALITY( zero->getGlobalNumCols(), numProcs*numLocal );
    TEST_EQUALITY( zero->getLocalNumCols(), 0 );
    TEST_EQUALITY( zero->getIndexBase(), 0 );

    TEST_EQUALITY( zero->getGlobalNumEntries(), 0 );
    TEST_EQUALITY( zero->getLocalNumEntries(), 0 );
    TEST_EQUALITY( zero->getGlobalMaxNumRowEntries(), 0 );
    TEST_EQUALITY( zero->getLocalMaxNumRowEntries(), 0 );
    STD_TESTS((*zero));

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, DottedDiag, LO, GO , Node )
  {
    using GRAPH = Tpetra::CrsGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int numProcs = comm->getSize();
    // create a Map, three rows per processor
    const size_t numLocal = 3;
    RCP<const map_type> map = rcp (new map_type (INVALID, numLocal, 0, comm));
    GO mymiddle = map->getGlobalElement(1);  // get my middle row

    for (int T=0; T<4; ++T) {
      if ( (T & 1) != 1 ) continue;
      RCP<ParameterList> params = parameterList ();
      params->set("Optimize Storage",((T & 2) == 2));

      // create a diagonal graph, but where only my middle row has an entry
      ArrayRCP<size_t> toalloc = arcpClone<size_t>( tuple<size_t>(0,1,0) );
      GRAPH ddgraph(map, toalloc ());
      ddgraph.insertGlobalIndices(mymiddle, tuple<GO>(mymiddle));
      {
        // before globalAssemble(), there should be one local entry on middle, none on the others
        typename GRAPH::global_inds_host_view_type myrow_gbl;
        ddgraph.getGlobalRowView(mymiddle-1,myrow_gbl); TEST_EQUALITY( myrow_gbl.size(), 0 );
        ddgraph.getGlobalRowView(mymiddle  ,myrow_gbl); TEST_COMPARE_ARRAYS( myrow_gbl, tuple<GO>(mymiddle) );
        ddgraph.getGlobalRowView(mymiddle+1,myrow_gbl); TEST_EQUALITY( myrow_gbl.size(), 0 );
      }
      TEST_THROW( ddgraph.insertGlobalIndices(mymiddle-1,tuple<GO>(mymiddle+1)), std::runtime_error );
      TEST_THROW( ddgraph.insertGlobalIndices(mymiddle  ,tuple<GO>(mymiddle+1)), std::runtime_error );
      TEST_THROW( ddgraph.insertGlobalIndices(mymiddle+1,tuple<GO>(mymiddle+1)), std::runtime_error );
      ddgraph.fillComplete(params);
      // after fillComplete(), there should be a single entry on my middle, corresponding to the diagonal, none on the others
      {
        typename GRAPH::local_inds_host_view_type myrow_lcl;
        TEST_EQUALITY_CONST( ddgraph.getNumEntriesInLocalRow(0), 0 );
        TEST_EQUALITY_CONST( ddgraph.getNumEntriesInLocalRow(2), 0 );
        ddgraph.getLocalRowView(1,myrow_lcl);
        TEST_EQUALITY_CONST( myrow_lcl.size(), 1 );
        if (myrow_lcl.size() == 1) {
          TEST_EQUALITY( ddgraph.getColMap()->getGlobalElement(myrow_lcl[0]), mymiddle );
        }
      }
      // also, the row map and column map should be equivalent
      TEST_EQUALITY( ddgraph.getGlobalNumCols(), static_cast<GST> (3*numProcs) );
      TEST_EQUALITY( ddgraph.getGlobalNumRows(), ddgraph.getGlobalNumCols() );

      STD_TESTS(ddgraph);
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, NonLocals, LO, GO , Node )
  {
    using Teuchos::as;
    using std::endl;
    using GRAPH = Tpetra::CrsGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    out << "Tpetra::CrsGraph: Test insert into nonowned rows" << endl;
    Teuchos::OSTab tab0 (out);

    // Get a communicator and Kokkos Node instance
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    // Create a Map with one row per process
    const size_t numLocal = 1;
    RCP<const map_type> map = rcp (new map_type (INVALID, numLocal, 0, comm));

    GO myrowind = map->getGlobalElement (0);

    Teuchos::OSTab tab1 (out);

    {
      Teuchos::OSTab tab2 (out);
      for (bool optimizeStorage : {false, true}) {
        out << "Optimize Storage: " << (optimizeStorage ? "true" : "false")
            << endl;
        Teuchos::OSTab tab3 (out);

        RCP<ParameterList> params = parameterList ();
        params->set ("Optimize Storage", optimizeStorage);

        // Diagonal graph
        {
          out << "Diagonal graph test" << endl;
          Teuchos::OSTab tab4 (out);

          // create a diagonal graph, where the graph entries are
          // contributed by a single off-node contribution, no
          // filtering.  let node i contribute to row i+1, where node
          // the last node contributes to row 0
          GRAPH diaggraph (map, 1);
          GO grow = myRank+1;
          if (as<int> (grow) == numProcs) {
            grow = 0;
          }
          diaggraph.insertGlobalIndices (grow, tuple<GO> (grow));
          // before globalAssemble(), there should be no local entries if 
          // numProcs > 1
          {
            typename GRAPH::global_inds_host_view_type myrow_gbl;
            diaggraph.getGlobalRowView (myrowind, myrow_gbl);
            TEST_EQUALITY( myrow_gbl.size (), (numProcs == 1 ? 1 : 0) );
          }
          diaggraph.globalAssemble ();
          // after globalAssemble(), there should be one local entry per
          // row, corresponding to the diagonal
          {
            typename GRAPH::global_inds_host_view_type myrow_gbl;
            diaggraph.getGlobalRowView (myrowind, myrow_gbl);
            TEST_COMPARE_ARRAYS( myrow_gbl, tuple<GO> (myrowind) );
          }

          { // no room for more
            out << "Attempt to insert global column index " << (myrowind+1) 
                << " into global row " << myrowind 
                << "; it should throw, because the graph"
                << " has an upper bound of one entry "
                << "per row, and already has a different column index " 
                << grow << " in this row."
                << endl;
            TEST_THROW( diaggraph.insertGlobalIndices(myrowind,
                                                      tuple<GO>(myrowind+1)),
                        std::runtime_error );
          }

          diaggraph.fillComplete (params);

          // after fillComplete(), there should be a single entry on my
          // row, corresponding to the diagonal
          {
            typename GRAPH::local_inds_host_view_type myrow_lcl;
            diaggraph.getLocalRowView (0, myrow_lcl);
            TEST_EQUALITY_CONST( myrow_lcl.size (), 1 );
            if (myrow_lcl.size() == 1) {
              TEST_EQUALITY( 
                   diaggraph.getColMap()->getGlobalElement(myrow_lcl[0]),
                   myrowind );
            }
          }
          // also, the row map and column map should be equivalent
          TEST_EQUALITY_CONST(
               diaggraph.getRowMap()->isSameAs(*diaggraph.getColMap()), true );

          STD_TESTS(diaggraph);
        }

        // Next-door-neighbor graph
        {
          out << "Next-door-neighbor graph test" << endl;
          Teuchos::OSTab tab4 (out);

          // create a next-door-neighbor graph (tridiagonal plus
          // corners), where the graph entries are contributed by single
          // off-node contribution, no filtering.  let node i add the
          // contributions for column i of the graph: (i-1,i), (i,i),
          // (i+1,i). allocate only as much space as we need. some
          // hacking here to support this test when numProcs == 1 or 2
          GRAPH ngraph(map,3);
          Array<GO> grows(3);
          grows[0] = (numProcs+myRank-1) % numProcs;   // my left neighbor
          grows[1] = (numProcs+myRank  ) % numProcs;   // myself
          grows[2] = (numProcs+myRank+1) % numProcs;   // my right neighbor

          // Add me to the graph for my neighbors
          ngraph.insertGlobalIndices (grows[0], tuple<GO> (myRank));
          ngraph.insertGlobalIndices (grows[1], tuple<GO> (myRank));
          ngraph.insertGlobalIndices (grows[2], tuple<GO> (myRank));

          // after globalAssemble(), storage should be maxed out
          out << "Calling globalAssemble()" << endl;
          ngraph.globalAssemble();
          TEST_EQUALITY( ngraph.getNumEntriesInLocalRow(0),
                         (numProcs == 1 ? 1 
                                        : ngraph.getNumAllocatedEntriesInLocalRow(0) ));
          out << "Calling fillComplete(params)" << endl;
          ngraph.fillComplete (params);

          {
            // after fillComplete(), there should be entries for me and my
            // neighbors on my row
            typename GRAPH::local_inds_host_view_type myrow_lcl;
            ngraph.getLocalRowView (0, myrow_lcl);
            // check indices on my row
            typename Array<GO>::iterator glast;
            std::sort (grows.begin (), grows.end ());
            glast = std::unique (grows.begin (), grows.end ());
            size_t numunique = glast - grows.begin ();
            // test the test: numunique == min(numProcs,3)
            TEST_EQUALITY( numunique, (size_t)min(numProcs,3) );
            TEST_EQUALITY_CONST( (size_t)myrow_lcl.size(), numunique );
            if ((size_t)myrow_lcl.size() == numunique) {
              size_t numinds;
              typename GRAPH::nonconst_global_inds_host_view_type inds("right",numunique),inds_short("short",numunique-1),inds_long("long",numunique+1);
              TEST_THROW(
                   ngraph.getGlobalRowCopy (myrowind, inds_short,
                                            numinds),
                   std::runtime_error );
              TEST_NOTHROW(
                   ngraph.getGlobalRowCopy (myrowind, inds_long,
                                            numinds) );
              TEST_NOTHROW(
                   ngraph.getGlobalRowCopy (myrowind,inds, numinds) );
              Tpetra::sort(inds, numinds);
              TEST_COMPARE_ARRAYS( arcp_from_view(inds,numinds), grows (0, numunique) );

              out << "On Proc 0:" << endl;
              Teuchos::OSTab tab5 (out);
              out << "numinds: " << numinds << endl
                  << "inds(0,numinds): " << arcp_from_view(inds, numinds) << endl
                  << "numunique: " << numunique << endl
                  << "grows(0,numunique): " << grows (0, numunique) << endl;
            }
          }
          TEST_EQUALITY_CONST( ngraph.getRowMap ()->isSameAs (* (ngraph.getColMap ())),
                               (numProcs == 1 ? true : false) );

          out << "Concluding with standard graph tests" << endl;
          STD_TESTS(ngraph);
        }
      } // optimizeStorage
    }

    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, ConstrFromRowPtrsAndLocalInds, LO, GO , Node )
  {
    using crs_graph_type = Tpetra::CrsGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    using row_ptrs_device_view_type = typename crs_graph_type::row_ptrs_device_view_type;
    using local_inds_dualv_type = typename crs_graph_type::local_inds_dualv_type;
    using local_inds_wdv_type = typename crs_graph_type::local_inds_wdv_type;

    // Get a communicator
    const auto comm = Tpetra::getDefaultComm();

    // Create a Map with one row per process
    constexpr size_t numLocal = 1;
    const auto map = Teuchos::make_rcp<map_type>(INVALID, numLocal, 0, comm);

    // Create a diagonal graph by creating the local row pointer and column indices arrays and passing them
    // to the CrsGraph constructor that takes these arrays as arguments. Note that the CrsGraph class stores
    // the local column indices array as a wrapped dual view. Here, we construct this local column indices array
    // as a wrapped dual view and pass it as such to the CrsGraph constructor.
    const typename row_ptrs_device_view_type::non_const_type rowPointers(Kokkos::view_alloc(Kokkos::WithoutInitializing, "row pointers"), 2);
    {
      const auto rowPointers_h = Kokkos::create_mirror_view(Kokkos::view_alloc(Kokkos::WithoutInitializing), rowPointers);
      rowPointers_h(0) = 0; rowPointers_h(1) = 1;
      Kokkos::deep_copy(rowPointers, rowPointers_h);
    }
    
    local_inds_wdv_type columnIndices;
    {
      local_inds_dualv_type columnIndices_("column indices", 1);
      columnIndices = local_inds_wdv_type(columnIndices_);
    }

    const auto crs_graph = Teuchos::make_rcp<crs_graph_type>(
      rowPointers,
      columnIndices,
      map,
      map,
      Teuchos::null,
      Teuchos::null,
      Teuchos::null,
      Teuchos::null
    );

    // Check that after the construction of the CrsGraph, the use counts of the host and device views of the local
    // column indices wrapped dual view are equal. This property is key to ensure that the sync mechanism of the
    // wrapped dual view works correctly.
    TEST_EQUALITY(columnIndices.host_view_use_count(), columnIndices.device_view_use_count());

    // Check that the right arrays are stored in the constructed CrsGraph instance. Note that these checks also
    // check that we can call getHostView() and getDeviceView() on the local column indices wrapped dual view.
    TEST_ASSERT(crs_graph->getLocalRowPtrsDevice() == rowPointers);
    TEST_ASSERT(crs_graph->getLocalIndicesHost() == columnIndices.getHostView(Tpetra::Access::ReadOnly));
    TEST_ASSERT(crs_graph->getLocalIndicesDevice() == columnIndices.getDeviceView(Tpetra::Access::ReadOnly));
  }

//
// INSTANTIATIONS
//

// Tests to build and run.  We will instantiate them over all enabled
// LocalOrdinal (LO), GlobalOrdinal (GO), and Node (NODE) types.
#define UNIT_TEST_GROUP( LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, EmptyGraphAlloc0,              LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, EmptyGraphAlloc1,              LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, ExcessAllocation,              LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, BadConst,                      LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, insert_remove_LIDs,            LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, NonLocals,                     LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, DottedDiag,                    LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, WithStaticProfile,             LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, CopiesAndViews,                LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, BadGIDs,                       LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, ConstrFromRowPtrsAndLocalInds, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

} // namespace (anonymous)
