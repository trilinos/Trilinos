/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_CrsGraph.hpp>
#include "Tpetra_Details_determineLocalTriangularStructure.hpp"
#include <Tpetra_TestingUtilities.hpp>
#include <type_traits> // std::is_same

namespace {

  template<class LO, class GO, class NT>
  Tpetra::Details::LocalTriangularStructureResult<LO>
  getLocalTriangularStructure (const Tpetra::RowGraph<LO, GO, NT>& G)
  {
    using Tpetra::Details::determineLocalTriangularStructure;
    using crs_graph_type = Tpetra::CrsGraph<LO, GO, NT>;

    const crs_graph_type& G_crs = dynamic_cast<const crs_graph_type&> (G);

    auto G_lcl = G_crs.getLocalGraph ();
    auto lclRowMap = G.getRowMap ()->getLocalMap ();
    auto lclColMap = G.getColMap ()->getLocalMap ();
    return determineLocalTriangularStructure (G_lcl, lclRowMap, lclColMap, true);
  }

  using Tpetra::DynamicProfile;
  using Tpetra::ProfileType;
  using Tpetra::StaticProfile;
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

#define STD_TESTS(graph) \
  { \
    auto STCOMM = graph.getComm(); \
    auto STMYGIDS = graph.getRowMap()->getNodeElementList(); \
    size_t STMAX = 0; \
    for (size_t STR = 0; STR < graph.getNodeNumRows(); ++STR) { \
      TEST_EQUALITY( graph.getNumEntriesInLocalRow (STR), graph.getNumEntriesInGlobalRow (STMYGIDS[STR]) ); \
      STMAX = std::max (STMAX, graph.getNumEntriesInLocalRow(STR)); \
    } \
    TEST_EQUALITY( graph.getNodeMaxNumRowEntries(), STMAX ); \
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

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, BadConst, LO, GO , Node )
  {
    using Teuchos::REDUCE_MIN;
    using map_type = Tpetra::Map<LO, GO, Node>;
    using graph_type = Tpetra::CrsGraph<LO, GO, Node>;
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();

    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?

    // Create a contiguous Map to use as the graph's row Map.
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
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
    TEST_THROW( badgraph = rcp (new graph_type (map, hints.persistingView (0, numLocal+1))),
                std::invalid_argument ); // too many entries

    // Make sure that the test passed on all processes.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // CrsGraph's constructor should throw if input allocation hint
    // has the wrong length, as it does in this case.
    TEST_THROW( badgraph = rcp (new graph_type (map, hints.persistingView (0, numLocal-1))),
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
    // TEST_THROW( badgraph = rcp (new graph_type (map, hints.persistingView (0, numLocal))),
    //             std::invalid_argument ); // invalid value

    // // Make sure that the test passed on all processes.
    // lclSuccess = success ? 1 : 0;
    // gblSuccess = 0;
    // reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    // TEST_EQUALITY_CONST( gblSuccess, 1 );
  }


// mfh 05 Apr 2013: CrsGraph only tests for bad nonowned GIDs in a
// debug build.  The BadGIDs test fails in a release build.
#ifdef HAVE_TPETRA_DEBUG
  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, BadGIDs, LO, GO , Node )
  {
    using GRAPH = Tpetra::CrsGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
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
      gids[0] = myRank*numLocal+numLocal;    // off this node, except on the last proc, where it is off the map
      // bad gid on the last node (not in domain map), discovered at fillComplete()
      GRAPH goodgraph(map,1);
      goodgraph.insertGlobalIndices(map->getMinGlobalIndex(),gids());
      TEST_THROW( goodgraph.fillComplete(), std::runtime_error );
    }
    // All procs fail if any process fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }
#endif // HAVE_TPETRA_DEBUG


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, ExcessAllocation, LO, GO , Node )
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    typedef Tpetra::CrsGraph<LO, GO, Node>  GRPH;
    using map_type = Tpetra::Map<LO, GO, Node>;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    // create a Map with numLocal entries per node
    const int numLocal = 10;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocal, indexBase, comm));

    // Test (GitHub Issue) #2565 fix, while we're at it.
    //
    // clts == 0: Don't set "compute local triangular constants"
    // clts == 1: Set "compute local triangular constants" to false
    // clts == 2: Set "compute local triangular constants" to true
    for (int clts : {0, 1, 2}) {
      // send in a parameterlist, check the defaults
      RCP<ParameterList> params = parameterList();
      if (clts == 1) {
        params->set ("compute local triangular constants", false);
      }
      else if (clts == 2) {
        params->set ("compute local triangular constants", true);
      }

      // create static-profile graph, fill-complete without inserting (and therefore, without allocating)
      GRPH graph (map, 3, StaticProfile);
      for (GO i = map->getMinGlobalIndex(); i <= map->getMaxGlobalIndex(); ++i) {
        graph.insertGlobalIndices (i, tuple<GO> (i));
      }
      params->set("Optimize Storage",false);
      graph.fillComplete(params);
      TEST_EQUALITY(graph.getNodeNumEntries(), (size_t)numLocal);
      TEST_EQUALITY_CONST(graph.isStorageOptimized(), false);
      //
      graph.resumeFill();
      for (int i=0; i < numLocal; ++i) graph.removeLocalIndices(i);
      params->set("Optimize Storage",true);
      graph.fillComplete(params);
      TEST_EQUALITY_CONST(params->get<bool>("Optimize Storage"), true);
      TEST_EQUALITY(graph.getNodeNumEntries(), 0);
      TEST_EQUALITY_CONST(graph.getProfileType(), StaticProfile);
      TEST_EQUALITY_CONST(graph.isStorageOptimized(), true);
    }

    // Test (GitHub Issue) #2565 fix, while we're at it.
    //
    // clts == 0: Don't set "compute local triangular constants"
    // clts == 1: Set "compute local triangular constants" to false
    // clts == 2: Set "compute local triangular constants" to true
    for (int clts : {0, 1, 2}) {
      // send in a parameterlist, check the defaults
      RCP<ParameterList> params = parameterList();
      if (clts == 1) {
        params->set ("compute local triangular constants", false);
      }
      else if (clts == 2) {
        params->set ("compute local triangular constants", true);
      }

      // create dynamic-profile graph, fill-complete without inserting (and therefore, without allocating)
      GRPH graph (map, 3, DynamicProfile);
      for (GO i = map->getMinGlobalIndex(); i <= map->getMaxGlobalIndex(); ++i) {
        graph.insertGlobalIndices (i, tuple<GO> (i));
      }
      params->set("Optimize Storage",false);
      graph.fillComplete(params);
      TEST_EQUALITY(graph.getNodeNumEntries(), (size_t)numLocal);
      TEST_EQUALITY_CONST(graph.isStorageOptimized(), false);
      //
      graph.resumeFill();
      for (int i=0; i < numLocal; ++i) graph.removeLocalIndices(i);
      params->set("Optimize Storage",true);
      graph.fillComplete(params);
      TEST_EQUALITY_CONST(params->get<bool>("Optimize Storage"), true);
      TEST_EQUALITY(graph.getNodeNumEntries(), 0);
      TEST_EQUALITY_CONST(graph.getProfileType(), StaticProfile);
      TEST_EQUALITY_CONST(graph.isStorageOptimized(), true);
    }

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


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, insert_remove_LIDs, LO, GO , Node )
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using GRAPH = Tpetra::CrsGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int myRank = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocal, indexBase, comm));
    {
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
      TEST_EQUALITY( static_cast<size_t> (diaggraph.getNumEntriesInLocalRow (row)),
                     static_cast<size_t> (lids.size ()) );
    }

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


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, CopiesAndViews, LO, GO , Node )
  {
    using GRAPH = Tpetra::CrsGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
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

    // Test (GitHub Issue) #2565 fix, while we're at it.
    //
    // clts == 0: Don't set "compute local triangular constants"
    // clts == 1: Set "compute local triangular constants" to false
    // clts == 2: Set "compute local triangular constants" to true
    for (int clts : {0, 1, 2}) {
      RCP<ParameterList> params = parameterList();
      if (clts == 1) {
        params->set ("compute local triangular constants", false);
      }
      else if (clts == 2) {
        params->set ("compute local triangular constants", true);
      }

      for (int T=0; T<4; ++T) {
        ProfileType pftype = ( (T & 1) == 1 ) ? StaticProfile : DynamicProfile;
        params->set("Optimize Storage",((T & 2) == 2));
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
        typedef typename Teuchos::ArrayView<const GO>::size_type size_type;
        for (size_type j=0; j < ginds.size(); ++j) {
          trigraph.insertGlobalIndices(myrowind,ginds(j,1));
        }
        TEST_EQUALITY( trigraph.getNumEntriesInLocalRow(0), trigraph.getNumAllocatedEntriesInLocalRow(0) ); // test that we only allocated as much room as necessary
        // If StaticProfile, then attempt to insert one additional entry
        // in my row that is not already in the row, and verify that it
        // throws an exception.
        if (pftype == StaticProfile) {
          TEST_THROW( trigraph.insertGlobalIndices(myrowind,tuple<GO>(myrowind+2)), std::runtime_error );
        }
        trigraph.fillComplete(params);
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
      reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, WithStaticProfile, LO, GO , Node )
  {
    using GRAPH = Tpetra::CrsGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    // what happens when we call CrsGraph::submitEntry() for a row that isn't on the Map?
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int myRank = comm->getRank();
    // create a Map, one row per processor
    const size_t numLocal = 1;
    RCP<const map_type> map = rcp (new map_type (INVALID, numLocal, 0, comm));
    {
      // add too many entries to a static graph
      // let node i contribute to row i+1, where node the last node contributes to row 0
      GRAPH diaggraph(map,1,StaticProfile);
      GO grow = myRank;
      Array<GO> colinds(1);
      colinds[0] = grow;
      TEST_NOTHROW( diaggraph.insertGlobalIndices(grow,colinds()) );
      TEST_THROW( diaggraph.insertGlobalIndices(grow, Teuchos::tuple<GO> (grow+1)), std::runtime_error );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, EmptyGraphAlloc0, LO, GO , Node )
  {
    using crs_graph_type = Tpetra::CrsGraph<LO, GO, Node>;
    using row_graph_type = Tpetra::RowGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    out << "CrsGrap EmptyGraphAlloc0 test" << endl;
    Teuchos::OSTab tab0 (out);

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
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

      // Test (GitHub Issue) #2565 fix, while we're at it.
      //
      // clts == 0: Don't set "compute local triangular constants"
      // clts == 1: Set "compute local triangular constants" to false
      // clts == 2: Set "compute local triangular constants" to true
      // clts == 3: Don't give a ParameterList to fillComplete at all
      for (int clts : {0, 1, 2, 3}) {

        if (clts == 0) {
          out << "Don't set \"compute local triangular constants\"" << endl;
        }
        else if (clts == 1) {
          out << "Set \"compute local triangular constants to false\"" << endl;
        }
        else if (clts == 2) {
          out << "Set \"compute local triangular constants to true\"" << endl;
        }
        else if (clts == 3) {
          out << "Don't give a ParameterList to fillComplete at all" << endl;
        }
        Teuchos::OSTab tab1 (out);

        RCP<row_graph_type> test_row;
        {
          // allocate with no space
          RCP<crs_graph_type> test_crs = rcp (new crs_graph_type (map, 0));
          // invalid, because none are allocated yet
          TEST_EQUALITY_CONST( test_crs->getNodeAllocationSize(), STINV );
          if (myRank != 1) {
            test_crs->insertGlobalIndices (map->getMinGlobalIndex (),
                                           tuple<GO> (map->getMinGlobalIndex ()));
          }
          if (clts == 3) {
            test_crs->fillComplete ();
          }
          else {
            RCP<ParameterList> params (new ParameterList);
            if (clts == 1) {
              params->set ("compute local triangular constants", false);
            }
            else if (clts == 2) {
              params->set ("compute local triangular constants", true);
            }
            test_crs->fillComplete (params);
          }
          TEST_EQUALITY( test_crs->getNodeAllocationSize(), numLocal );
          test_row = test_crs;
        }
        RCP<const map_type> cmap = test_row->getColMap();
        TEST_EQUALITY( cmap->getGlobalNumElements(), (size_t)numProcs-1 );
        TEST_EQUALITY( test_row->getGlobalNumRows(), (size_t)numProcs-1 );
        TEST_EQUALITY( test_row->getNodeNumRows(), numLocal );
        TEST_EQUALITY( test_row->getGlobalNumCols(), (size_t)numProcs-1 );
        TEST_EQUALITY( test_row->getNodeNumCols(), numLocal );
        TEST_EQUALITY( test_row->getIndexBase(), 0 );

        auto lclTri = getLocalTriangularStructure (*test_row);
        TEST_ASSERT( lclTri.couldBeLowerTriangular );
        TEST_ASSERT( lclTri.couldBeUpperTriangular );
        TEST_EQUALITY( lclTri.diagCount, static_cast<LO> (numLocal) );
        GO gblDiagCount = 0;
        reduceAll<int, GO> (*comm, REDUCE_SUM, static_cast<GO> (lclTri.diagCount), outArg (gblDiagCount));
        TEST_EQUALITY( gblDiagCount, static_cast<GO> (numProcs-1) );

        TEST_EQUALITY( test_row->getGlobalNumEntries(), (size_t)numProcs-1 );
        TEST_EQUALITY( test_row->getNodeNumEntries(), numLocal );
        TEST_EQUALITY( test_row->getGlobalMaxNumRowEntries(), 1 );
        TEST_EQUALITY( test_row->getNodeMaxNumRowEntries(), numLocal );
        STD_TESTS((*test_row));
      }
    }

    // this one is empty on all nodes because of zero allocation size
    {
      // create a Map
      const size_t numLocal = 10;
      RCP<const map_type> map =
        rcp (new map_type (INVALID, numLocal, 0, comm));

      // Test (GitHub Issue) #2565 fix, while we're at it.
      //
      // clts == 0: Don't set "compute local triangular constants"
      // clts == 1: Set "compute local triangular constants" to false
      // clts == 2: Set "compute local triangular constants" to true
      // clts == 3: Don't give a ParameterList to fillComplete at all
      for (int clts : {0, 1, 2, 3}) {
        RCP<row_graph_type> zero;
        {
          // allocate with no space
          RCP<crs_graph_type> zero_crs = rcp (new crs_graph_type (map, 0));
          // invalid, because none are allocated yet
          TEST_EQUALITY_CONST( zero_crs->getNodeAllocationSize(), STINV );
          if (clts == 3) {
            zero_crs->fillComplete ();
          }
          else {
            RCP<ParameterList> params (new ParameterList);
            if (clts == 1) {
              params->set ("compute local triangular constants", false);
            }
            else if (clts == 2) {
              params->set ("compute local triangular constants", true);
            }
            zero_crs->fillComplete (params);
          }
          // zero, because none were allocated.
          TEST_EQUALITY_CONST( zero_crs->getNodeAllocationSize(), 0 );
          zero = zero_crs;
        }
        RCP<const map_type> cmap = zero->getColMap();
        TEST_EQUALITY( cmap->getGlobalNumElements(), 0 );
        TEST_EQUALITY( zero->getGlobalNumRows(), numProcs*numLocal );
        TEST_EQUALITY( zero->getNodeNumRows(), numLocal );
        TEST_EQUALITY( zero->getGlobalNumCols(), numProcs*numLocal );
        TEST_EQUALITY( zero->getNodeNumCols(), 0 );
        TEST_EQUALITY( zero->getIndexBase(), 0 );

        auto lclTri = getLocalTriangularStructure (*zero);
        TEST_ASSERT( lclTri.couldBeLowerTriangular );
        TEST_ASSERT( lclTri.couldBeUpperTriangular );
        TEST_EQUALITY( lclTri.diagCount, static_cast<LO> (0) );
        GO gblDiagCount = 0;
        reduceAll<int, GO> (*comm, REDUCE_SUM, static_cast<GO> (lclTri.diagCount), outArg (gblDiagCount));
        TEST_EQUALITY( gblDiagCount, static_cast<GO> (0) );

        TEST_EQUALITY( zero->getGlobalNumEntries(), 0 );
        TEST_EQUALITY( zero->getNodeNumEntries(), 0 );
        TEST_EQUALITY( zero->getGlobalMaxNumRowEntries(), 0 );
        TEST_EQUALITY( zero->getNodeMaxNumRowEntries(), 0 );
        STD_TESTS((*zero));
      }
    }
    // this one is empty on all nodes because of zero rows
    {
      // create a Map
      const size_t numLocal = 0;
      RCP<const map_type> map = rcp (new map_type (INVALID, numLocal, 0, comm));

      // Test (GitHub Issue) #2565 fix, while we're at it.
      //
      // clts == 0: Don't set "compute local triangular constants"
      // clts == 1: Set "compute local triangular constants" to false
      // clts == 2: Set "compute local triangular constants" to true
      // clts == 3: Don't give a ParameterList to fillComplete at all
      for (int clts : {0, 1, 2, 3}) {
        RCP<row_graph_type> zero;
        {
          // allocate with no space
          RCP<crs_graph_type> zero_crs = rcp (new crs_graph_type (map, 0));
          // invalid, because none are allocated yet
          TEST_EQUALITY_CONST( zero_crs->getNodeAllocationSize(), STINV );
          if (clts == 3) {
            zero_crs->fillComplete ();
          }
          else {
            RCP<ParameterList> params (new ParameterList);
            if (clts == 1) {
              params->set ("compute local triangular constants", false);
            }
            else if (clts == 2) {
              params->set ("compute local triangular constants", true);
            }
            zero_crs->fillComplete (params);
          }
          // zero, because none were allocated.
          TEST_EQUALITY_CONST( zero_crs->getNodeAllocationSize(), 0 );
          zero = zero_crs;
        }
        RCP<const map_type> cmap = zero->getColMap();
        TEST_EQUALITY( cmap->getGlobalNumElements(), 0 );
        TEST_EQUALITY( zero->getGlobalNumRows(), numProcs*numLocal );
        TEST_EQUALITY( zero->getNodeNumRows(), numLocal );
        TEST_EQUALITY( zero->getGlobalNumCols(), numProcs*numLocal );
        TEST_EQUALITY( zero->getNodeNumCols(), 0 );
        TEST_EQUALITY( zero->getIndexBase(), 0 );

        auto lclTri = getLocalTriangularStructure (*zero);
        TEST_ASSERT( lclTri.couldBeLowerTriangular );
        TEST_ASSERT( lclTri.couldBeUpperTriangular );
        TEST_EQUALITY( lclTri.diagCount, static_cast<LO> (0) );
        GO gblDiagCount = 0;
        reduceAll<int, GO> (*comm, REDUCE_SUM, static_cast<GO> (lclTri.diagCount), outArg (gblDiagCount));
        TEST_EQUALITY( gblDiagCount, static_cast<GO> (0) );

        TEST_EQUALITY( zero->getGlobalNumEntries(), 0 );
        TEST_EQUALITY( zero->getNodeNumEntries(), 0 );
        TEST_EQUALITY( zero->getGlobalMaxNumRowEntries(), 0 );
        TEST_EQUALITY( zero->getNodeMaxNumRowEntries(), 0 );
        STD_TESTS((*zero));
      }
    }

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, EmptyGraphAlloc1, LO, GO , Node )
  {
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    using map_type = Tpetra::Map<LO, GO, Node>;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    int numProcs = size(*comm);
    // create a Map
    const size_t numLocal = 10;
    RCP<const map_type> map = rcp (new map_type (INVALID, numLocal, 0, comm));
    // create the empty graph
    RCP<Tpetra::RowGraph<LO,GO,Node> > zero;

    // Test (GitHub Issue) #2565 fix, while we're at it.
    //
    // clts == 0: Don't set "compute local triangular constants"
    // clts == 1: Set "compute local triangular constants" to false
    // clts == 2: Set "compute local triangular constants" to true
    // clts == 3: Don't give a ParameterList to fillComplete at all
    for (int clts : {0, 1, 2, 3}) {
      {
        // allocated with space for one entry per row
        RCP<graph_type> zero_crs = rcp (new graph_type (map,1));
        TEST_EQUALITY( zero_crs->getNodeAllocationSize(), STINV ); // zero, because none are allocated yet

        if (clts == 3) {
          zero_crs->fillComplete ();
        }
        else {
          RCP<ParameterList> params (new ParameterList);
          if (clts == 1) {
            params->set ("compute local triangular constants", false);
          }
          else if (clts == 2) {
            params->set ("compute local triangular constants", true);
          }
          zero_crs->fillComplete (params);
        }
        zero = zero_crs;
      }
      RCP<const map_type> cmap = zero->getColMap();
      TEST_EQUALITY( cmap->getGlobalNumElements(), 0 );
      TEST_EQUALITY( zero->getGlobalNumRows(), numProcs*numLocal );
      TEST_EQUALITY( zero->getNodeNumRows(), numLocal );
      TEST_EQUALITY( zero->getGlobalNumCols(), numProcs*numLocal );
      TEST_EQUALITY( zero->getNodeNumCols(), 0 );
      TEST_EQUALITY( zero->getIndexBase(), 0 );

      auto lclTri = getLocalTriangularStructure (*zero);
      TEST_ASSERT( lclTri.couldBeLowerTriangular );
      TEST_ASSERT( lclTri.couldBeUpperTriangular );
      TEST_EQUALITY( lclTri.diagCount, static_cast<LO> (0) );
      GO gblDiagCount = 0;
      reduceAll<int, GO> (*comm, REDUCE_SUM, static_cast<GO> (lclTri.diagCount), outArg (gblDiagCount));
      TEST_EQUALITY( gblDiagCount, static_cast<GO> (0) );

      TEST_EQUALITY( zero->getGlobalNumEntries(), 0 );
      TEST_EQUALITY( zero->getNodeNumEntries(), 0 );
      TEST_EQUALITY( zero->getGlobalMaxNumRowEntries(), 0 );
      TEST_EQUALITY( zero->getNodeMaxNumRowEntries(), 0 );
      STD_TESTS((*zero));

      // All procs fail if any proc fails
      int globalSuccess_int = -1;
      reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, DottedDiag, LO, GO , Node )
  {
    using GRAPH = Tpetra::CrsGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int numProcs = comm->getSize();
    // create a Map, three rows per processor
    const size_t numLocal = 3;
    RCP<const map_type> map = rcp (new map_type (INVALID, numLocal, 0, comm));
    GO mymiddle = map->getGlobalElement(1);  // get my middle row

    for (int T=0; T<4; ++T) {
      ProfileType pftype = ( (T & 1) == 1 ) ? StaticProfile : DynamicProfile;
      RCP<ParameterList> params = parameterList ();
      params->set("Optimize Storage",((T & 2) == 2));

      // Test (GitHub Issue) #2565 fix, while we're at it.
      //
      // clts == 0: Don't set "compute local triangular constants"
      // clts == 1: Set "compute local triangular constants" to false
      // clts == 2: Set "compute local triangular constants" to true
      for (int clts : {0, 1, 2}) {
        if (clts == 1) {
          params->set ("compute local triangular constants", false);
        }
        else if (clts == 2) {
          params->set ("compute local triangular constants", true);
        }
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
          TEST_THROW( ddgraph.insertGlobalIndices(mymiddle-1,tuple<GO>(mymiddle+1)), std::runtime_error );
          TEST_THROW( ddgraph.insertGlobalIndices(mymiddle  ,tuple<GO>(mymiddle+1)), std::runtime_error );
          TEST_THROW( ddgraph.insertGlobalIndices(mymiddle+1,tuple<GO>(mymiddle+1)), std::runtime_error );
        }
        ddgraph.fillComplete(params);
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
        TEST_EQUALITY( ddgraph.getGlobalNumCols(), static_cast<GST> (3*numProcs) );
        TEST_EQUALITY( ddgraph.getGlobalNumRows(), ddgraph.getGlobalNumCols() );

        auto lclTri = getLocalTriangularStructure (ddgraph);
        TEST_EQUALITY( lclTri.diagCount, static_cast<LO> (1) );
        GO gblDiagCount = 0;
        reduceAll<int, GO> (*comm, REDUCE_SUM, static_cast<GO> (lclTri.diagCount), outArg (gblDiagCount));
        TEST_EQUALITY( gblDiagCount, static_cast<GO> (numProcs) );

        STD_TESTS(ddgraph);
      }
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, NonLocals, LO, GO , Node )
  {
    using Teuchos::as;
    using std::endl;
    using GRAPH = Tpetra::CrsGraph<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

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

    // Test (GitHub Issue) #2565 fix, while we're at it.
    //
    // clts == 0: Don't set "compute local triangular constants"
    // clts == 1: Set "compute local triangular constants" to false
    // clts == 2: Set "compute local triangular constants" to true
    for (int clts : {0, 1, 2}) {
      if (clts == 0) {
        out << "Don't set \"compute local triangular constants\"" << endl;
      }
      else if (clts == 1) {
        out << "Set \"compute local triangular constants\" to false" << endl;
      }
      else if (clts == 2) {
        out << "Set \"compute local triangular constants\" to true" << endl;
      }
      Teuchos::OSTab tab1 (out);

      for (ProfileType pftype : {StaticProfile, DynamicProfile}) {
        out << "ProfileType: "
            << (pftype == StaticProfile ? "StaticProfile" : "DynamicProfile")
            << endl;
        Teuchos::OSTab tab2 (out);
        for (bool optimizeStorage : {false, true}) {
          out << "Optimize Storage: " << (optimizeStorage ? "true" : "false")
              << endl;
          Teuchos::OSTab tab3 (out);

          RCP<ParameterList> params = parameterList ();
          if (clts == 1) {
            params->set ("compute local triangular constants", false);
          }
          else if (clts == 2) {
            params->set ("compute local triangular constants", true);
          }
          params->set ("Optimize Storage", optimizeStorage);

          // Diagonal graph
          {
            out << "Diagonal graph test" << endl;
            Teuchos::OSTab tab4 (out);

            // create a diagonal graph, where the graph entries are
            // contributed by a single off-node contribution, no
            // filtering.  let node i contribute to row i+1, where node
            // the last node contributes to row 0
            GRAPH diaggraph (map, 1, pftype);
            GO grow = myRank+1;
            if (as<int> (grow) == numProcs) {
              grow = 0;
            }
            diaggraph.insertGlobalIndices (grow, tuple<GO> (grow));
            // before globalAssemble(), there should be no local entries if numProcs > 1
            ArrayView<const GO> myrow_gbl;
            diaggraph.getGlobalRowView (myrowind, myrow_gbl);
            TEST_EQUALITY( myrow_gbl.size (), (numProcs == 1 ? 1 : 0) );
            diaggraph.globalAssemble ();
            // after globalAssemble(), there should be one local entry per
            // row, corresponding to the diagonal
            diaggraph.getGlobalRowView (myrowind, myrow_gbl);
            TEST_COMPARE_ARRAYS( myrow_gbl, tuple<GO> (myrowind) );

            if (pftype == StaticProfile) { // no room for more
              out << "Attempt to insert global column index " << (myrowind+1) << " into"
                " global row " << myrowind << "; it should throw, because the graph"
                " is StaticProfile, has an upper bound of one entry per row, and "
                "already has a different column index " << grow << " in this row."
                  << endl;
              TEST_THROW( diaggraph.insertGlobalIndices(myrowind,tuple<GO>(myrowind+1)),
                          std::runtime_error );
            }

            diaggraph.fillComplete (params);

            // after fillComplete(), there should be a single entry on my
            // row, corresponding to the diagonal
            ArrayView<const LO> myrow_lcl;
            diaggraph.getLocalRowView (0, myrow_lcl);
            TEST_EQUALITY_CONST( myrow_lcl.size (), 1 );
            if (myrow_lcl.size() == 1) {
              TEST_EQUALITY( diaggraph.getColMap ()->getGlobalElement (myrow_lcl[0]),
                             myrowind );
            }
            // also, the row map and column map should be equivalent
            TEST_EQUALITY_CONST( diaggraph.getRowMap()->isSameAs(*diaggraph.getColMap()), true );

            auto lclTri = getLocalTriangularStructure (diaggraph);
            TEST_ASSERT( lclTri.couldBeLowerTriangular );
            TEST_ASSERT( lclTri.couldBeUpperTriangular );
            TEST_EQUALITY( lclTri.diagCount, static_cast<LO> (1) );
            GO gblDiagCount = 0;
            reduceAll<int, GO> (*comm, REDUCE_SUM, static_cast<GO> (lclTri.diagCount), outArg (gblDiagCount));
            TEST_EQUALITY( gblDiagCount, static_cast<GO> (numProcs) );

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
            GRAPH ngraph(map,3,pftype);
            Array<GO> grows(3);
            grows[0] = (numProcs+myRank-1) % numProcs;   // my left neighbor
            grows[1] = (numProcs+myRank  ) % numProcs;   // myself
            grows[2] = (numProcs+myRank+1) % numProcs;   // my right neighbor

            // Add me to the graph for my neighbors
            ngraph.insertGlobalIndices (grows[0], tuple<GO> (myRank));
            ngraph.insertGlobalIndices (grows[1], tuple<GO> (myRank));
            ngraph.insertGlobalIndices (grows[2], tuple<GO> (myRank));

            // before globalAssemble(), there should be a single local
            // entry on parallel runs, three on serial runs
            ArrayView<const GO> myrow_gbl;
            ngraph.getGlobalRowView (myrowind, myrow_gbl);
            TEST_EQUALITY_CONST( myrow_gbl.size(), (numProcs == 1 ? 3 : 1) );

            // after globalAssemble(), storage should be maxed out
            out << "Calling globalAssemble()" << endl;
            ngraph.globalAssemble();
            TEST_EQUALITY( ngraph.getNumEntriesInLocalRow(0),
                           ngraph.getNumAllocatedEntriesInLocalRow(0) );
            out << "Calling fillComplete(params)" << endl;
            ngraph.fillComplete (params);

            // after fillComplete(), there should be entries for me and my
            // neighbors on my row
            ArrayView<const LO> myrow_lcl;
            ngraph.getLocalRowView (0, myrow_lcl);
            out << "Returned view of column indices on Proc 0: "
                << Teuchos::toString (myrow_lcl) << endl;

            {
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
                Array<GO> inds(numunique+1);
                TEST_THROW(
                           ngraph.getGlobalRowCopy (myrowind, inds (0, numunique-1), numinds),
                           std::runtime_error );
                TEST_NOTHROW(
                             ngraph.getGlobalRowCopy (myrowind, inds (0, numunique), numinds) );
                TEST_NOTHROW( ngraph.getGlobalRowCopy (myrowind,inds (), numinds) );
                std::sort (inds.begin (), inds.begin () + numinds);
                TEST_COMPARE_ARRAYS( inds (0, numinds), grows (0, numunique) );

                out << "On Proc 0:" << endl;
                Teuchos::OSTab tab5 (out);
                out << "numinds: " << numinds << endl
                    << "inds(0,numinds): " << inds (0, numinds) << endl
                    << "numunique: " << numunique << endl
                    << "grows(0,numunique): " << grows (0, numunique) << endl;
              }
            }
            TEST_EQUALITY_CONST( ngraph.getRowMap ()->isSameAs (* (ngraph.getColMap ())),
                                 (numProcs == 1 ? true : false) );

            auto lclTri = getLocalTriangularStructure (ngraph);
            TEST_EQUALITY( lclTri.diagCount, 1 );
            GO gblDiagCount = 0;
            reduceAll<int, GO> (*comm, REDUCE_SUM, static_cast<GO> (lclTri.diagCount), outArg (gblDiagCount));
            TEST_EQUALITY( gblDiagCount, static_cast<GO> (numProcs) );

            out << "Concluding with standard graph tests" << endl;
            STD_TESTS(ngraph);
          }
        } // optimizeStorage
      } // pftype
    } // clts

    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, NodeConversion, LO, GO, N2 )
  {
    using std::cerr;
    using std::endl;
    using Teuchos::REDUCE_MIN;
    using N1 = Tpetra::Map<>::node_type; // default Tpetra Node type
    using Map1 = Tpetra::Map<LO, GO, N1>;
    using Graph1 = Tpetra::CrsGraph<LO, GO, N1>;
    using Graph2 = Tpetra::CrsGraph<LO, GO, N2>;

    out << "Tpetra test: CrsGraph, NodeConversion" << endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();
    out << "Number of processes: " << numProcs << endl;

    const size_t numLocal = 10;
    const GST numGlobal = numProcs*numLocal;

    // create a contiguous uniform distributed map with numLocal entries per node
    out << "Creating Map" << endl;
    RCP<const Map1> map1 = rcp (new Map1 (numGlobal,0,comm));
    out << "Creating CrsGraph" << endl;
    RCP<Graph1>       A1 = createCrsGraph(map1,3);
    RCP<N2> n2; // this can be null; we keep it to help type deduction

    // empty source, not filled

    {
      out << "Testing clone (1)" << endl;

      RCP<ParameterList> plClone = parameterList();
      // default: plClone->set("fillComplete clone",true);
      RCP<Graph2> A2;
      try {
        A2 = A1->template clone<N2> (n2, plClone);
      } catch (std::exception& e) {
        std::ostringstream err;
        err << "Process " << myRank << ": clone raised exception: "
            << e.what () << endl;
        cerr << err.str ();
      }

      int lclSuccess = (success && ! A2.is_null ()) ? 1 : 0;
      int gblSuccess = 1;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));

      if (gblSuccess == 1) {
        out << "Clone succeeded on all processes!" << endl;
      } else {
        out << "Clone FAILED on at least one process!" << endl;
      }
      TEST_EQUALITY_CONST( gblSuccess , 1)

      if (gblSuccess == 1 && ! A2.is_null ()) {
        out << "Testing status of newly created graph" << endl;
        TEST_EQUALITY_CONST( A2->isFillComplete(), true );
        TEST_EQUALITY_CONST( A2->isStorageOptimized(), true );
        TEST_EQUALITY_CONST( A2->getNodeNumEntries(), (size_t)0 );
        TEST_EQUALITY_CONST( A2->getNodeAllocationSize(), (size_t)0 );
      }
    }

    out << "Filling graph" << endl;
    // one entry per row
    for (GO grow =map1->getMinGlobalIndex();
            grow<=map1->getMaxGlobalIndex();
            ++grow)
    {
      if (grow == map1->getMinGlobalIndex())      A1->insertGlobalIndices(grow, tuple<GO>(grow,grow+1));
      else if (grow == map1->getMaxGlobalIndex()) A1->insertGlobalIndices(grow, tuple<GO>(grow-1,grow));
      else                                        A1->insertGlobalIndices(grow, tuple<GO>(grow-1,grow,grow+1));
    }
    // source has global indices, not filled, dynamic profile
    out << "Testing clone (2)" << endl;
    {
      RCP<ParameterList> plClone = parameterList();
      plClone->set("fillComplete clone",false);
      plClone->set("Static profile clone",false);
      // default: plClone->set("Locally indexed clone",false);
      RCP<Graph2> A2 = A1->template clone<N2>(n2,plClone);
      TEST_EQUALITY_CONST( A2->hasColMap(), false );
      TEST_EQUALITY_CONST( A2->isFillComplete(), false );
      TEST_EQUALITY_CONST( A2->isGloballyIndexed(), true );
      TEST_EQUALITY_CONST( A2->getNodeAllocationSize(), (size_t)(numLocal*3-2) );
      TEST_EQUALITY( A2->getNodeNumEntries(), A1->getNodeNumEntries() );
      TEST_NOTHROW( A2->insertGlobalIndices( map1->getMaxGlobalIndex(), tuple<GO>(map1->getMinGlobalIndex()) ) );
      TEST_NOTHROW( A2->insertGlobalIndices( map1->getMinGlobalIndex(), tuple<GO>(map1->getMaxGlobalIndex()) ) );
      TEST_NOTHROW( A2->fillComplete() );
      TEST_EQUALITY_CONST( A2->getNodeNumEntries(), A1->getNodeNumEntries()+2 );
    }

    // source has local indices
    out << "Calling fillComplete on original" << endl;
    A1->fillComplete();

    out << "Testing clone (3)" << endl;
    {
      RCP<ParameterList> plClone = parameterList();
      plClone->set("Static profile clone", false);
      RCP<ParameterList> plCloneFill = sublist(plClone,"fillComplete");
      plCloneFill->set("Optimize Storage",false);
      RCP<Graph2> A2 = A1->template clone<N2>(n2,plClone);

      out << "Finished clone; testing result" << endl;
      TEST_EQUALITY_CONST( A2->isFillComplete(), true );
      TEST_EQUALITY_CONST( A2->isStorageOptimized(), false );
      TEST_EQUALITY_CONST( A2->getNodeNumEntries(), A1->getNodeNumEntries() );

      out << "Calling resumeFill on result of clone" << endl;
      A2->resumeFill();

      out << "Filling result of clone" << endl;
      for (LO lrow = map1->getMinLocalIndex();
              lrow < map1->getMaxLocalIndex()-1;
              ++lrow)
      {
        TEST_NOTHROW( A2->insertLocalIndices(lrow, tuple<LO>(lrow+2)) );
      }
      out << "Calling fillComplete on result of clone" << endl;
      A2->fillComplete();

      // Call computeGlobal constants multiple times to make sure nothing crashes
      A2->computeGlobalConstants();
      A2->computeGlobalConstants();
      A2->computeGlobalConstants();


      out << "Testing result of clone" << endl;
      TEST_EQUALITY_CONST( A2->isFillComplete(), true );
      TEST_EQUALITY_CONST( A2->isStorageOptimized(), true );
      TEST_EQUALITY_CONST( A2->getNodeNumEntries(), A1->getNodeNumEntries()+numLocal-2 );
    }

    {
      int lclSuccess = success ? 1 : 0;
      int gblSuccess = 1;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST( gblSuccess , 1)
      if (gblSuccess == 1) {
        out << "Test succeeded on all processes!" << endl;
      } else {
        out << "Test FAILED on at least one process!" << endl;
      }
    }
  }

//
// INSTANTIATIONS
//

// Tests to build and run in both debug and release modes.  We will
// instantiate them over all enabled local ordinal (LO), global
// ordinal (GO), and Kokkos Node (NODE) types.
#define UNIT_TEST_GROUP_DEBUG_AND_RELEASE( LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, EmptyGraphAlloc0,  LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, EmptyGraphAlloc1,  LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, ExcessAllocation,  LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, BadConst,          LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, insert_remove_LIDs, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, NonLocals,         LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, DottedDiag,        LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, WithStaticProfile, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, CopiesAndViews,    LO, GO, NODE )

// Test(s) for "Node conversion" (i.e., the clone() template method of
// CrsGraph).  We will instantiate them over all enabled Kokkos Node
// (N2) types.
#define NC_TESTS(N2)
  //    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, NodeConversion, int, int, N2 )

// mfh 05 Apr 2013: CrsGraph only tests for bad nonowned GIDs in a
// debug build.  The BadGIDs test fails in a release build.
#ifdef HAVE_TPETRA_DEBUG

#define UNIT_TEST_GROUP_DEBUG_ONLY( LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, BadGIDs, LO, GO, NODE )

#else // NOT HAVE_TPETRA_DEBUG

#define UNIT_TEST_GROUP_DEBUG_ONLY( LO, GO, NODE )

#endif // HAVE_TPETRA_DEBUG


    TPETRA_ETI_MANGLING_TYPEDEFS()

    TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP_DEBUG_AND_RELEASE )

    TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP_DEBUG_ONLY )

    TPETRA_INSTANTIATE_N(NC_TESTS)
}


