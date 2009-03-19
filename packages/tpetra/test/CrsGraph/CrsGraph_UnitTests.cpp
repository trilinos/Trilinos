#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_as.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsGraph.hpp"

#ifdef HAVE_TPETRA_TRIUTILS
#include <iohb.h>
#endif

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Tpetra::Map;
  using Tpetra::DefaultPlatform;
  using Tpetra::Platform;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::arrayView;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;
  using Tpetra::MultiVector;
  using std::endl;
  using std::swap;
  using std::min;
  using std::max;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Tpetra::CrsGraph;
  using Tpetra::RowGraph;
  using Tpetra::INSERT;
  using Tpetra::Import;
  using std::string;
  using Teuchos::tuple;

  bool testMpi = true;
  double errorTolSlack = 1e+1;
  string filedir;

#define STD_TESTS(graph) \
  { \
    RCP<const Comm<int> > STCOMM = graph.getComm(); \
    ArrayView<const GO> STMYGIDS = graph.getRowMap().getMyGlobalEntries(); \
    Teuchos_Ordinal STMAX = 0; \
    for (Teuchos_Ordinal STR=0; STR<graph.numLocalRows(); ++STR) { \
      TEST_EQUALITY( graph.numEntriesForMyRow(STR), graph.numEntriesForGlobalRow( STMYGIDS[STR] ) ); \
      STMAX = std::max( STMAX, graph.numEntriesForMyRow(STR) ); \
    } \
    TEST_EQUALITY( graph.myMaxNumRowEntries(), STMAX ); \
    GO STGMAX; \
    reduceAll( *STCOMM, Teuchos::REDUCE_MAX, STMAX, &STGMAX ); \
    TEST_EQUALITY( graph.globalMaxNumRowEntries(), STGMAX ); \
  }


  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption(
        "filedir",&filedir,"Directory of expected matrix files.");
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, BadConst, LO, GO )
  {
    typedef CrsGraph<LO,GO> GRAPH;
    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const Teuchos_Ordinal indexBase = 0;
    const Teuchos_Ordinal numLocal = 10;
    Map<LO,GO> map(INVALID,numLocal,indexBase,comm);
    {
      // bad constructor
      TEST_THROW( GRAPH badgraph(map,-1), std::runtime_error ); // allocation hint must be >= 0
    }
    {
      // bad constructor
      Array<Teuchos_Ordinal> hints(numLocal+1);
      std::fill(hints.begin(),hints.end(),1);
      hints[0] = -1;    // invalid
      TEST_THROW( GRAPH badgraph(map,hints(0,numLocal+1)), std::runtime_error ); // too many
      TEST_THROW( GRAPH badgraph(map,hints(0,numLocal-1)), std::runtime_error ); // too few
      TEST_THROW( GRAPH badgraph(map,hints(0,numLocal))  , std::runtime_error ); // one is invalid
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, BadGIDs, LO, GO )
  {
    typedef CrsGraph<LO,GO> GRAPH;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    // create a Map
    const Teuchos_Ordinal indexBase = 0;
    const Teuchos_Ordinal numLocal = 10;
    Map<LO,GO> map(INVALID,numLocal,indexBase,comm);
    {
      Array<GO> gids(1);
      gids[0] = myImageID*numLocal+numLocal;    // off this node, except on the last proc, where it is off the map
      // bad gid on the last node (not in domain map), discovered at fillComplete()
      GRAPH goodgraph(map,1);
      goodgraph.insertGlobalIndices(map.getMinGlobalIndex(),gids());
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
      goodgraph.insertGlobalIndices(map.getMinGlobalIndex(),gids());
      goodgraph.fillComplete();
      TEST_EQUALITY( goodgraph.numEntriesForGlobalRow(map.getMinGlobalIndex()),
                     (myImageID == numImages-1) ? 0 : 1 );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, BadLIDs, LO, GO )
  {
    typedef CrsGraph<LO,GO> GRAPH;
    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    // create a Map
    const Teuchos_Ordinal indexBase = 0;
    const Teuchos_Ordinal numLocal = 10;
    Map<LO,GO> map(INVALID,numLocal,indexBase,comm);
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
        TEST_THROW(diaggraph.insertMyIndices(numLocal,lids()), std::runtime_error);
      }
      {
        GRAPH diaggraph(map,map,1);
        TEST_EQUALITY(diaggraph.hasColMap(), true);
        diaggraph.insertMyIndices(0,lids());
        TEST_EQUALITY( diaggraph.numEntriesForMyRow(0), (myImageID == numImages-1) ? 0 : 1 );
      }
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, CopiesAndViews, LO, GO )
  {
    typedef CrsGraph<LO,GO> GRAPH;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    if (numImages < 2) return;
    // create a Map, one row per processor
    const Teuchos_Ordinal indexBase = 0;
    const Teuchos_Ordinal numLocal = 1;
    Map<LO,GO> rmap(INVALID,numLocal,indexBase,comm);
    GO myrowind = rmap.getGlobalIndex(0);
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
    for (unsigned int i=0; i<linds.size(); ++i) linds[i] = i;
    Map<LO,GO> cmap(INVALID,ginds,0,comm);
    GRAPH trigraph(rmap,cmap,3);
    Array<GO> GCopy(4); Array<LO> LCopy(4);
    ArrayView<const GO> GView; ArrayView<const LO> LView;
    Teuchos_Ordinal numindices;
    // at this point, there are no global or local indices, but views and copies should succeed
    trigraph.extractMyRowCopy(0,LCopy,numindices);
    trigraph.extractMyRowConstView(0,LView);
    trigraph.extractGlobalRowCopy(myrowind,GCopy,numindices);
    trigraph.extractGlobalRowConstView(myrowind,GView);
    trigraph.insertGlobalIndices(myrowind,ginds);
    trigraph.fillComplete();
    // check for throws and no-throws/values
    TEST_THROW( trigraph.extractGlobalRowConstView(myrowind,GView    ), std::runtime_error );
    TEST_THROW( trigraph.extractMyRowCopy(    0       ,LCopy(0,1),numindices), std::runtime_error );
    TEST_THROW( trigraph.extractGlobalRowCopy(myrowind,GCopy(0,1),numindices), std::runtime_error );
    TEST_NOTHROW( trigraph.extractMyRowConstView(0,LView) );
    TEST_COMPARE_ARRAYS( LView, linds );
    TEST_NOTHROW( trigraph.extractMyRowCopy(0,LCopy,numindices) );
    TEST_COMPARE_ARRAYS( LCopy(0,numindices), linds );
    TEST_NOTHROW( trigraph.extractGlobalRowCopy(myrowind,GCopy,numindices) );
    TEST_COMPARE_ARRAYS( GCopy(0,numindices), ginds );
    STD_TESTS(trigraph);
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, NonLocals, LO, GO )
  {
    typedef CrsGraph<LO,GO> GRAPH;
    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    if (numImages < 3) return;
    // create a Map, one row per processor
    const Teuchos_Ordinal indexBase = 0;
    const Teuchos_Ordinal numLocal = 1;
    Map<LO,GO> map(INVALID,numLocal,indexBase,comm);
    GO myrowind = map.getGlobalIndex(0);
    {
      // create a diagonal graph, where the graph entries are contributed by a single off-node contribution, no filtering
      // let node i contribute to row i+1, where node the last node contributes to row 0
      GRAPH diaggraph(map,1);
      GO grow = myImageID+1;
      if (grow == numImages) {
        grow = 0;
      }
      Array<GO> colinds(1);
      colinds[0] = grow;
      diaggraph.insertGlobalIndices(grow,colinds());
      // before the fillComplete() (and the associated globalAssemble()), there should be no local entries
      ArrayView<const GO> myrow_gbl;
      diaggraph.extractGlobalRowConstView(myrowind,myrow_gbl);
      TEST_EQUALITY_CONST( myrow_gbl.size(), 0 );
      TEST_NOTHROW( diaggraph.fillComplete() );
      // after fillComplete(), there should be a single entry on my row, corresponding to the diagonal
      ArrayView<const LO> myrow_lcl;
      diaggraph.extractMyRowConstView(0,myrow_lcl);
      TEST_EQUALITY_CONST( myrow_lcl.size(), 1 );
      if (myrow_lcl.size() == 1) {
        TEST_EQUALITY( diaggraph.getColMap().getGlobalIndex(myrow_lcl[0]), myrowind );
      }
      // also, the row map and column map should be equivalent
      TEST_EQUALITY_CONST( diaggraph.getRowMap().isSameAs(diaggraph.getColMap()), true );
      TEST_EQUALITY( diaggraph.numGlobalDiagonals(), numImages );
      TEST_EQUALITY( diaggraph.numMyDiagonals(), 1 );
      STD_TESTS(diaggraph);
    }
    {
      // create a next-door-neighbor graph (tridiagonal plus corners), where the graph entries are contributed by single off-node contribution, no filtering
      // let node i add the contributions for column i of the graph: (i-1,i), (i,i), (i+1,i)
      // require at least two procs
      GRAPH trigraph(map,3);
      Array<GO> colinds(1), grows(3);
      colinds[0] = myImageID;
      grows[0] = (numImages+myImageID-1) % numImages;
      grows[1] = (numImages+myImageID  ) % numImages;
      grows[2] = (numImages+myImageID+1) % numImages;
      trigraph.insertGlobalIndices(grows[0],colinds());
      trigraph.insertGlobalIndices(grows[1],colinds());
      trigraph.insertGlobalIndices(grows[2],colinds());
      // before the fillComplete() (and the associated globalAssemble()), there should be a single local entry
      ArrayView<const GO> myrow_gbl;
      trigraph.extractGlobalRowConstView(myrowind,myrow_gbl);
      TEST_EQUALITY_CONST( myrow_gbl.size(), 1 );
      TEST_NOTHROW( trigraph.fillComplete() );
      // after fillComplete(), there should be three entries on my row
      ArrayView<const LO> myrow_lcl;
      trigraph.extractMyRowConstView(0,myrow_lcl);
      TEST_EQUALITY_CONST( myrow_lcl.size(), 3 );
      if (myrow_lcl.size() == 3) {
        sort(grows.begin(),grows.end());
        Teuchos_Ordinal numinds;
        Array<GO> inds(4);
        TEST_THROW(   trigraph.extractGlobalRowCopy(myrowind,inds(0,2),numinds), std::runtime_error );
        TEST_NOTHROW( trigraph.extractGlobalRowCopy(myrowind,inds(0,3),numinds) );
        TEST_NOTHROW( trigraph.extractGlobalRowCopy(myrowind,inds(0,4),numinds) );
        TEST_EQUALITY(numinds, 3);
        sort(inds.begin(),inds.begin()+numinds);
        TEST_COMPARE_ARRAYS( inds(0,3), grows );
      }
      TEST_EQUALITY_CONST( trigraph.getRowMap().isSameAs(trigraph.getColMap()), false );
      TEST_EQUALITY( trigraph.numGlobalDiagonals(), numImages );
      TEST_EQUALITY( trigraph.numMyDiagonals(), 1 );
      STD_TESTS(trigraph);
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraph, EmptyGraph, LO, GO )
  {
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    int numImages = size(*comm);
    // create a Map
    const Teuchos_Ordinal indexBase = 0;
    const Teuchos_Ordinal numLocal = 10;
    Map<LO,GO> map(INVALID,numLocal,indexBase,comm);
    // create the empty graph
    RCP<RowGraph<LO,GO> > zero;
    {
      RCP<CrsGraph<LO,GO> > zero_crs = rcp(new CrsGraph<LO,GO>(map,0));
      zero_crs->fillComplete();
      zero = zero_crs;
    }
    Map<LO,GO> cmap = zero->getColMap();
    TEST_EQUALITY( cmap.getNumGlobalEntries(), 0 );
    TEST_EQUALITY( zero->numGlobalRows(), numImages*numLocal );
    TEST_EQUALITY( zero->numLocalRows(), numLocal );
    TEST_EQUALITY( zero->numGlobalCols(), numImages*numLocal );
    TEST_EQUALITY( zero->numLocalCols(), 0 );
    TEST_EQUALITY( zero->getIndexBase(), 0 );
    TEST_EQUALITY( zero->upperTriangular(), true );
    TEST_EQUALITY( zero->lowerTriangular(), true );
    TEST_EQUALITY( zero->numGlobalDiagonals(), 0 );
    TEST_EQUALITY( zero->numMyDiagonals(), 0 );
    TEST_EQUALITY( zero->numGlobalEntries(), 0 );
    TEST_EQUALITY( zero->numMyEntries(), 0 );
    TEST_EQUALITY( zero->globalMaxNumRowEntries(), 0 );
    TEST_EQUALITY( zero->myMaxNumRowEntries(), 0 );
    STD_TESTS((*zero));

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
  #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_LO_GO( LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, EmptyGraph, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, BadConst  , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, BadGIDs   , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, BadLIDs   , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, NonLocals , LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraph, CopiesAndViews, LO, GO )

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
     UNIT_TEST_GROUP_LO_GO(int,int)
# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

     UNIT_TEST_GROUP_LO_GO(short,int)
     UNIT_TEST_GROUP_LO_GO(int,int)

     typedef long int LongInt;
     UNIT_TEST_GROUP_LO_GO(int,LongInt)
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
        typedef long long int LongLongInt;
        UNIT_TEST_GROUP_LO_GO(int,LongLongInt)
#    endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
