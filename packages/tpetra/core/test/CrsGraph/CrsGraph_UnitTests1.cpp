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
#include <Tpetra_TestingUtilities.hpp>
#include <type_traits> // std::is_same

namespace { // (anonymous)
  using Tpetra::TestingUtilities::getDefaultComm;
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
  using Teuchos::tuple;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using std::endl;
  using std::max;
  using std::min;
  typedef Tpetra::global_size_t GST;

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
    Teuchos::reduceAll<int, GST> (*STCOMM, Teuchos::REDUCE_MAX, STMAX, Teuchos::outArg (STGMAX)); \
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

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, EmptyFillComplete, LO, GO , Node )
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef Tpetra::CrsGraph<LO, GO, Node> GRAPH;
    typedef Tpetra::Map<LO, GO, Node> map_type;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map with numLocal entries per node
    const size_t numLocal = 10;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocal, indexBase, comm));
    {
      // create static-profile graph, fill-complete without inserting
      // (and therefore, without allocating)
      GRAPH graph(map,1,StaticProfile);
      graph.fillComplete();
    }
    {
      // create dynamic-profile graph, fill-complete without inserting
      // (and therefore, without allocating)
      GRAPH graph(map,1,DynamicProfile);
      graph.fillComplete();
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
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, SortingTests, LO, GO , Node )
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef Tpetra::CrsGraph<LO, GO, Node> GRAPH;
    typedef Tpetra::Map<LO, GO, Node> map_type;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 10;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocal, indexBase, comm));

    GRAPH graph(map,map,4);
    TEST_EQUALITY_CONST(graph.isSorted(), true);
    // insert entries; shouldn't be sorted anymore
    for (GO i = map->getMinGlobalIndex (); i <= map->getMaxGlobalIndex (); ++i) {
      Array<GO> curInds;
      const GO firstInd = (i+5) % map->getMaxAllGlobalIndex ();
      if (map->isNodeGlobalElement (firstInd)) {
        curInds.push_back (firstInd);
      }
      curInds.push_back (i);
      const GO lastInd = (i-5) % map->getMaxAllGlobalIndex ();
      if (map->isNodeGlobalElement (lastInd)) {
        curInds.push_back (lastInd);
      }
      graph.insertGlobalIndices (i, curInds ());
    }
    TEST_EQUALITY_CONST(graph.isSorted(), false);
    // fill complete; should be sorted now
    RCP<ParameterList> params = parameterList();
    params->set("Optimize Storage",false);
    graph.fillComplete(params);
    {
      bool sortingCheck = true;
      for (LO i=map->getMinLocalIndex(); i <= map->getMaxLocalIndex(); ++i) {
        ArrayView<const LO> inds;
        graph.getLocalRowView(i,inds);
        for (int j=1; j < (int)inds.size(); ++j) {
          if (inds[j-1] > inds[j]) {sortingCheck = false; break;}
        }
      }
      TEST_EQUALITY_CONST(sortingCheck, graph.isSorted() );
    }
    // resume fill; should still be sorted
    graph.resumeFill();
    TEST_EQUALITY_CONST(graph.isSorted(), true);
    {
      bool sortingCheck = true;
      for (LO i=map->getMinLocalIndex(); i <= map->getMaxLocalIndex(); ++i) {
        ArrayView<const LO> inds;
        graph.getLocalRowView(i,inds);
        for (int j=1; j < (int)inds.size(); ++j) {
          if (inds[j-1] > inds[j]) {sortingCheck = false; break;}
        }
      }
      TEST_EQUALITY_CONST(sortingCheck, graph.isSorted() );
    }
    // insert a column-index; currently, this invalidates sorting, though it may change in the future
    graph.insertLocalIndices(0, tuple<LO>(0));
    TEST_EQUALITY_CONST(graph.isSorted(), false);
    // fill complete, check one more time
    params->set("Optimize Storage",true);
    graph.fillComplete(params);
    {
      bool sortingCheck = true;
      for (LO i=map->getMinLocalIndex(); i <= map->getMaxLocalIndex(); ++i) {
        ArrayView<const LO> inds;
        graph.getLocalRowView(i,inds);
        for (int j=1; j < (int)inds.size(); ++j) {
          if (inds[j-1] > inds[j]) {sortingCheck = false; break;}
        }
      }
      TEST_EQUALITY_CONST(sortingCheck, graph.isSorted() );
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
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, Bug20100622K, LO, GO , Node )
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef Tpetra::CrsGraph<LO, GO, Node> GRAPH;
    typedef Tpetra::Map<LO, GO, Node> map_type;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numProcs = comm->getSize();
    const int myRank = comm->getRank();

    // test filtering
    if (numProcs > 1) {
      // we only need two procs to demonstrate this bug. ignore the others.
      Array<GO> mygids, domaingids;
      if (myRank == 0) {
        mygids.push_back(0);
        mygids.push_back(1);
        mygids.push_back(2);

        domaingids.push_back(0);
        domaingids.push_back(1);
      }
      else if (myRank == 1) {
        mygids.push_back(2);
        mygids.push_back(3);

        domaingids.push_back(2);
        domaingids.push_back(3);
      }
      RCP<const map_type> rmap =
        rcp (new map_type (INVALID, mygids (), 0, comm));

      RCP<const map_type> dmap =
        rcp (new map_type (INVALID, domaingids (), 0, comm));


      RCP<GRAPH> G = rcp(new GRAPH(rmap,2) );

      if (myRank == 0) {
        G->insertGlobalIndices(0, tuple<GO>(0));
        G->insertGlobalIndices(1, tuple<GO>(0,1));
        G->insertGlobalIndices(2, tuple<GO>(2,1));
      }
      else if (myRank == 1) {
        G->insertGlobalIndices(2, tuple<GO>(2, 1));
        G->insertGlobalIndices(3, tuple<GO>(3));
      }

      try {
        G->fillComplete(dmap,rmap);
      }
      catch (std::exception& e) {
        std::cerr << "G->fillComplete() raised an exception! "
                  << e.what () << std::endl;
        success = false;
      }
    }
    // All procs fail if any node fails
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
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, WithColMap, LO, GO , Node )
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef Tpetra::CrsGraph<LO, GO, Node> GRAPH;
    typedef Tpetra::Map<LO, GO, Node> map_type;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numProcs = comm->getSize();
    if (numProcs > 1) {
      const size_t numLocal = 1;
      RCP<const map_type> rmap =
        rcp (new map_type (INVALID, numLocal, 0, comm));
      RCP<const map_type> cmap =
        rcp (new map_type (INVALID, numLocal, 0, comm));
      // must allocate enough for all submitted indices.
      RCP<GRAPH> G = rcp(new GRAPH(rmap,cmap,2,StaticProfile) );
      TEST_EQUALITY_CONST( G->hasColMap(), true );
      const GO myrowind = rmap->getGlobalElement(0);

      // mfh 16 May 2013: insertGlobalIndices doesn't do column Map
      // filtering anymore, so we have to test whether each of the
      // column indices to insert is invalid.
      if (cmap->isNodeGlobalElement (myrowind)) {
        if (cmap->isNodeGlobalElement (myrowind+1)) {
          TEST_NOTHROW( G->insertGlobalIndices( myrowind, tuple<GO>(myrowind,myrowind+1) ) );
        }
        else {
          TEST_NOTHROW( G->insertGlobalIndices( myrowind, tuple<GO>(myrowind) ) );
        }
      }
      TEST_NOTHROW( G->fillComplete() );
      TEST_EQUALITY( G->getRowMap(), rmap );
      TEST_EQUALITY( G->getColMap(), cmap );
      TEST_EQUALITY( G->getNumEntriesInGlobalRow(myrowind), 1 );
    }
    // All procs fail if any node fails
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
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, Describable, LO, GO , Node )
  {
    typedef Tpetra::CrsGraph<LO, GO, Node> GRAPH;
    typedef Tpetra::Map<LO, GO, Node> map_type;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myRank = comm->getRank();
    // create Map
    RCP<const map_type> map = rcp (new map_type (INVALID, 3, 0, comm));
    {
      GRAPH graph(map,1,StaticProfile);
      // test labeling
      const std::string lbl("graphA");
      std::string desc1 = graph.description();
      out << desc1 << endl;
      graph.setObjectLabel(lbl);
      std::string desc2 = graph.description();
      out << desc2 << endl;
      TEST_EQUALITY( graph.getObjectLabel(), lbl );
    }
    {
      GRAPH graph(map,1,StaticProfile);
      // test describing at different verbosity levels
      if (myRank==0) out << "Describing with verbosity VERB_DEFAULT..." << endl;
      graph.describe(out);
      comm->barrier();
      comm->barrier();
      if (myRank==0) out << "Describing with verbosity VERB_NONE..." << endl;
      graph.describe(out,VERB_NONE);
      comm->barrier();
      comm->barrier();
      if (myRank==0) out << "Describing with verbosity VERB_LOW..." << endl;
      graph.describe(out,VERB_LOW);
      comm->barrier();
      comm->barrier();
      if (myRank==0) out << "Describing with verbosity VERB_MEDIUM..." << endl;
      graph.describe(out,VERB_MEDIUM);
      comm->barrier();
      comm->barrier();
      if (myRank==0) out << "Describing with verbosity VERB_HIGH..." << endl;
      graph.describe(out,VERB_HIGH);
      comm->barrier();
      comm->barrier();
      if (myRank==0) out << "Describing with verbosity VERB_EXTREME..." << endl;
      graph.describe(out,VERB_EXTREME);
      comm->barrier();
      comm->barrier();
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, ActiveFill, LO, GO , Node )
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef Tpetra::CrsGraph<LO, GO, Node> GRAPH;
    typedef Tpetra::Map<LO, GO, Node> map_type;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create Map
    RCP<const map_type> map = rcp (new map_type (INVALID, 1, 0, comm));
    RCP<ParameterList> params = parameterList();
    {
      GRAPH graph(map,map,0,DynamicProfile);
      TEST_EQUALITY_CONST( graph.isFillActive(),   true );
      TEST_EQUALITY_CONST( graph.isFillComplete(), false );
      graph.insertLocalIndices( 0, tuple<LO>(0) );
      //
      params->set("Optimize Storage",false);
      graph.fillComplete(params);
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
      params->set("Optimize Storage",false);
      graph.fillComplete(params);
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
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, Typedefs, LO, GO , Node )
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;

    typedef Tpetra::CrsGraph<LO, GO, Node> GRAPH;
    typedef typename GRAPH::local_ordinal_type  local_ordinal_type;
    typedef typename GRAPH::global_ordinal_type global_ordinal_type;
    typedef typename GRAPH::node_type           node_type;
    // Commas in macro invocations may cause trouble.
    {
      const bool LO_same = std::is_same<local_ordinal_type, LO>::value;
      TEST_ASSERT( LO_same );
      const bool GO_same = std::is_same<global_ordinal_type, GO>::value;
      TEST_ASSERT( GO_same );
      const bool Node_same = std::is_same<node_type, Node>::value;
      TEST_ASSERT( Node_same );
    }

    typedef Tpetra::RowGraph<LO,GO,Node> RGRAPH;
    typedef typename RGRAPH::local_ordinal_type  rgraph_local_ordinal_type;
    typedef typename RGRAPH::global_ordinal_type rgraph_global_ordinal_type;
    typedef typename RGRAPH::node_type           rgraph_node_type;
    {
      const bool LO_same = std::is_same<rgraph_local_ordinal_type, LO>::value;
      TEST_ASSERT( LO_same );
      const bool GO_same = std::is_same<rgraph_global_ordinal_type, GO>::value;
      TEST_ASSERT( GO_same );
      const bool Node_same = std::is_same<rgraph_node_type, Node>::value;
      TEST_ASSERT( Node_same );
    }

    RCP<const Comm<int> > comm = getDefaultComm ();
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));

    if (gblSuccess == 1) {
      out << "Succeeded on all processes!" << endl;
    } else {
      out << "FAILED on at least one process!" << endl;
    }
    TEST_EQUALITY_CONST(gblSuccess , 1);
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, TwoArraysESFC, LO, GO , Node )
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef Tpetra::CrsGraph<LO, GO, Node> GRAPH;
    typedef Tpetra::Map<LO, GO, Node> map_type;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numProcs = comm->getSize();
    // test filtering
    if (numProcs > 1) {
      const size_t numLocal = 2;
      const RCP<const map_type> rmap =
        rcp (new map_type (INVALID, numLocal, 0, comm));
      ArrayRCP<size_t> rowptr(numLocal+1);
      ArrayRCP<LO>     colind(numLocal); // one unknown per row
      rowptr[0] = 0; rowptr[1] = 1; rowptr[2] = 2;
      colind[0] = Teuchos::as<LO>(0);
      colind[1] = Teuchos::as<LO>(1);

      RCP<GRAPH> G = rcp(new GRAPH(rmap,rmap,rowptr,colind) );
      TEST_EQUALITY_CONST( G->hasColMap(), true );

      TEST_NOTHROW( G->expertStaticFillComplete(rmap,rmap) );
      TEST_EQUALITY( G->getRowMap(), rmap );
      TEST_EQUALITY( G->getColMap(), rmap );
    }

    // All procs fail if any node fails
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
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, SetAllIndices, LO, GO , Node )
  {
    typedef Tpetra::CrsGraph<LO, GO, Node> GRAPH;
    typedef Tpetra::Map<LO, GO, Node> map_type;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numProcs = comm->getSize();
    // test filtering
    if (numProcs > 1) {
      const size_t numLocal = 2;
      RCP<const map_type> rmap =
        rcp (new map_type (INVALID, numLocal, 0, comm));
      ArrayRCP<size_t> rowptr(numLocal+1);
      ArrayRCP<LO>     colind(numLocal); // one unknown per row
      rowptr[0] = 0; rowptr[1] = 1; rowptr[2] = 2;
      colind[0] = Teuchos::as<LO>(0);
      colind[1] = Teuchos::as<LO>(1);

      RCP<GRAPH> G = rcp(new GRAPH(rmap,rmap,0,StaticProfile) );
      TEST_NOTHROW( G->setAllIndices(rowptr,colind) );
      TEST_EQUALITY_CONST( G->hasColMap(), true );

      TEST_NOTHROW( G->expertStaticFillComplete(rmap,rmap) );
      TEST_EQUALITY( G->getRowMap(), rmap );
      TEST_EQUALITY( G->getColMap(), rmap );
    }

    // All procs fail if any node fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

//
// INSTANTIATIONS
//

// Tests to build and run in both debug and release modes.  We will
// instantiate them over all enabled local ordinal (LO), global
// ordinal (GO), and Kokkos Node (NODE) types.
#define UNIT_TEST_GROUP_DEBUG_AND_RELEASE( LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, WithColMap,        LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, Describable,       LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, EmptyFillComplete, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, Typedefs,          LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, Bug20100622K,      LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, ActiveFill,        LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, SortingTests,      LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, TwoArraysESFC,     LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, SetAllIndices,     LO, GO, NODE )

    TPETRA_ETI_MANGLING_TYPEDEFS()

    TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP_DEBUG_AND_RELEASE )

} // namespace (anonymous)


