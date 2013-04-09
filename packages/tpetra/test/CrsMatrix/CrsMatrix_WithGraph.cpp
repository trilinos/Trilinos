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

#include <Tpetra_TestingUtilities.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
// mfh 08 Mar 2013: This include isn't being used here, so I'm
// commenting it out to speed up compilation time.
//#include <Tpetra_CrsMatrixMultiplyOp.hpp>

// TODO: add test where some nodes have zero rows
// TODO: add test where non-"zero" graph is used to build matrix; if no values are added to matrix, the operator effect should be zero. This tests that matrix values are initialized properly.
// TODO: add test where dynamic profile initially has no allocation, then entries are added. this will test new view functionality.

namespace Teuchos {
  template <>
    ScalarTraits<int>::magnitudeType
    relErr( const int &s1, const int &s2 )
    {
      typedef ScalarTraits<int> ST;
      return ST::magnitude(s1-s2);
    }

  template <>
    ScalarTraits<char>::magnitudeType
    relErr( const char &s1, const char &s2 )
    {
      typedef ScalarTraits<char> ST;
      return ST::magnitude(s1-s2);
    }
}

namespace {

  // no ScalarTraits<>::eps() for integer types

  template <class Scalar, bool hasMachineParameters> struct TestingTolGuts {};

  template <class Scalar>
  struct TestingTolGuts<Scalar, true> {
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType testingTol()
      { return Teuchos::ScalarTraits<Scalar>::eps(); }
  };

  template <class Scalar>
  struct TestingTolGuts<Scalar, false> {
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType testingTol()
      { return 0; }
  };

  template <class Scalar>
  static typename Teuchos::ScalarTraits<Scalar>::magnitudeType testingTol()
  {
    return TestingTolGuts<Scalar, Teuchos::ScalarTraits<Scalar>::hasMachineParameters>::
      testingTol();
  }

  using Tpetra::TestingUtilities::getNode;
  using Tpetra::TestingUtilities::getDefaultComm;

  using std::endl;
  using std::swap;

  using std::string;

  using Teuchos::TypeTraits::is_same;
  using Teuchos::as;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::arcp;
  using Teuchos::outArg;
  using Teuchos::arcpClone;
  using Teuchos::arrayView;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using Teuchos::ETransp;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Teuchos::EDiag;
  using Teuchos::UNIT_DIAG;
  using Teuchos::NON_UNIT_DIAG;
  using Teuchos::EUplo;
  using Teuchos::UPPER_TRI;
  using Teuchos::LOWER_TRI;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;

  using Tpetra::Map;
  using Tpetra::MultiVector;
  using Tpetra::Vector;
  using Tpetra::Operator;
  using Tpetra::CrsMatrix;
  using Tpetra::CrsGraph;
  using Tpetra::RowMatrix;
  using Tpetra::Import;
  using Tpetra::global_size_t;
  using Tpetra::createNonContigMapWithNode;
  using Tpetra::createUniformContigMapWithNode;
  using Tpetra::createContigMapWithNode;
  using Tpetra::createLocalMapWithNode;
  // mfh 08 Mar 2013: This isn't being used here, so I'm commenting it
  // out to save compilation time.
  //using Tpetra::createCrsMatrixMultiplyOp;
  using Tpetra::createVector;
  using Tpetra::createCrsMatrix;
  using Tpetra::DefaultPlatform;
  using Tpetra::ProfileType;
  using Tpetra::StaticProfile;
  using Tpetra::DynamicProfile;
  using Tpetra::OptimizeOption;
  using Tpetra::DoOptimizeStorage;
  using Tpetra::DoNotOptimizeStorage;
  using Tpetra::GloballyDistributed;
  using Tpetra::INSERT;


  double errorTolSlack = 1e+1;
  string filedir;

template <class tuple, class T>
inline void tupleToArray(Array<T> &arr, const tuple &tup)
{
  arr.assign(tup.begin(), tup.end());
}

#define STD_TESTS(matrix) \
  { \
    using Teuchos::outArg; \
    RCP<const Comm<int> > STCOMM = matrix.getComm(); \
    ArrayView<const GO> STMYGIDS = matrix.getRowMap()->getNodeElementList(); \
    ArrayView<const LO> loview; \
    ArrayView<const Scalar> sview; \
    size_t STMAX = 0; \
    for (size_t STR=0; STR < matrix.getNodeNumRows(); ++STR) { \
      const size_t numEntries = matrix.getNumEntriesInLocalRow(STR); \
      TEST_EQUALITY( numEntries, matrix.getNumEntriesInGlobalRow( STMYGIDS[STR] ) ); \
      matrix.getLocalRowView(STR,loview,sview); \
      TEST_EQUALITY( static_cast<size_t>(loview.size()), numEntries ); \
      TEST_EQUALITY( static_cast<size_t>( sview.size()), numEntries ); \
      STMAX = std::max( STMAX, numEntries ); \
    } \
    TEST_EQUALITY( matrix.getNodeMaxNumRowEntries(), STMAX ); \
    global_size_t STGMAX; \
    Teuchos::reduceAll<int,global_size_t>( *STCOMM, Teuchos::REDUCE_MAX, STMAX, outArg(STGMAX) ); \
    TEST_EQUALITY( matrix.getGlobalMaxNumRowEntries(), STGMAX ); \
  }


  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption(
        "filedir",&filedir,"Directory of expected matrix files.");
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
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, AdvancedGraphUsage, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // generate a tridiagonal matrix
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map with numLocal entries per node
    const size_t numLocal = 10;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    {
      CrsGraph<LO,GO,Node> diaggraph(map,1,StaticProfile);
      // A pre-constructed graph must be fill complete before being used to construct a CrsMatrix
      TEST_THROW( MAT matrix(rcpFromRef(diaggraph)), std::runtime_error );
    }
    {
      // create a simple diagonal graph
      CrsGraph<LO,GO,Node> diaggraph(map,1,StaticProfile);
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        diaggraph.insertGlobalIndices(r,tuple(r));
      }
      // fill-complete the graph, but do not optimize the storage
      RCP<ParameterList> params = parameterList();
      params->set("Optimize Storage",false);
      diaggraph.fillComplete(params);
      TEST_EQUALITY_CONST( diaggraph.isFillComplete(), true );
      TEST_EQUALITY_CONST( diaggraph.isStorageOptimized(), false );
      // matrix constructed with non-storage-optimized graph
      MAT mat1(rcpFromRef(diaggraph));
      // fill complete the matrix and ask it to optimize storage.
      // this is not allowed on a static graph, and will either throw an exception or ignore the request to optimize storage.
      params->set("Optimize Storage",false);
#ifdef HAVE_TPETRA_THROW_ABUSE_WARNINGS
      TEST_THROW( mat1.fillComplete(params), std::runtime_error );
      TEST_EQUALITY_CONST( mat1.isFillComplete(), false );
      TEST_EQUALITY_CONST( mat1.isStorageOptimized(), false );
#else
      mat1.fillComplete(params);
      TEST_EQUALITY_CONST( mat1.isFillComplete(), true );
      TEST_EQUALITY_CONST( mat1.isStorageOptimized(), false );
#endif
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, WithGraph, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // generate a tridiagonal matrix
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node>  MAT;
    typedef CrsGraph<LO,GO,Node>         GRPH;
    typedef Vector<Scalar,LO,GO,Node>       V;
    typedef typename ST::magnitudeType    Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Scalar SONE  = ST::one();
    const Scalar SZERO = ST::zero();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = size(*comm);
    // create a Map with numLocal entries per node
    const size_t numLocal = 10;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    RCP<ParameterList> params = parameterList();
    params->set("Optimize Storage",true);
    RCP<ParameterList> fillparams = sublist(params,"Local Sparse Ops");
    fillparams->set("Prepare Solve", true);
    {
      //////////////////////////////////
      // create a simple tridiagonal graph
      GRPH trigraph(map,3,StaticProfile);
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        if (r == map->getMinAllGlobalIndex()) {
          trigraph.insertGlobalIndices(r,tuple(r,r+1));
        }
        else if (r == map->getMaxAllGlobalIndex()) {
          trigraph.insertGlobalIndices(r,tuple(r-1,r));
        }
        else {
          trigraph.insertGlobalIndices(r,tuple(r-1,r,r+1));
        }
      }
      trigraph.fillComplete(params);
      // create a matrix using the tri-diagonal graph and test allowed functionality
      MAT matrix(rcpFromRef(trigraph));
      TEST_EQUALITY_CONST( matrix.getProfileType() == StaticProfile, true );
      // insert throws exception: not allowed with static graph
      TEST_THROW( matrix.insertGlobalValues(map->getMinGlobalIndex(),tuple<GO>(map->getMinGlobalIndex()),tuple(ST::one())), std::runtime_error );
      // suminto and replace are allowed
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        if (r == map->getMinAllGlobalIndex()) {
          matrix.replaceGlobalValues(r, tuple(r,r+1), tuple(ST::one(),ST::one()) );
        }
        else if (r == map->getMaxAllGlobalIndex()) {
          matrix.replaceGlobalValues(r, tuple(r-1,r), tuple(ST::one(),ST::one()) );
        }
        else {
          matrix.replaceGlobalValues(r, tuple(r-1,r,r+1), tuple(ST::one(),ST::one(),ST::one()) );
        }
      }
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        // increment the diagonals
        matrix.sumIntoGlobalValues(r, tuple(r), tuple(ST::one()) );
      }
      matrix.fillComplete();
      TEST_EQUALITY( matrix.getNodeNumDiags(), numLocal );
      TEST_EQUALITY( matrix.getGlobalNumDiags(), numImages*numLocal );
      TEST_EQUALITY( matrix.getGlobalNumEntries(), 3*numImages*numLocal - 2 );
      V dvec(map,false);
      dvec.randomize();
      matrix.getLocalDiagCopy(dvec);
      Array<Scalar> expectedDiags(numLocal, static_cast<Scalar>(2));
      ArrayRCP<const Scalar> dvec_view = dvec.get1dView();
      if (ST::isOrdinal) {
        TEST_COMPARE_ARRAYS(expectedDiags(), dvec_view);
      } else {
        TEST_COMPARE_FLOATING_ARRAYS( expectedDiags(), dvec_view, MT::zero() );
      }
    }
    {
      // create a simple diagonal graph
      GRPH diaggraph(map,1,StaticProfile);
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        diaggraph.insertGlobalIndices(r,tuple(r));
      }
      diaggraph.fillComplete(params);
      // Bug verification:
      //  Tpetra::CrsMatrix constructed with a graph was experiencing a seg-fault if setAllToScalar is called before
      //  some other call allocates memory. This was because setAllToScalar has an incorrect if-statement that is
      //  not allocating memory.
      // This bug has been fixed. Furthermore, CrsMatrix no longer utilizes lazy allocation when constructed with a graph.
      // However, we will leave this test in place, because it still demonstrates valid behavior.
      MAT matrix(rcpFromRef(diaggraph));
      TEST_NOTHROW( matrix.setAllToScalar( ST::one() ) );
    }
    {
      // create a simple diagonal graph
      GRPH diaggraph(map,1,StaticProfile);
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        diaggraph.insertGlobalIndices(r,tuple(r));
      }
      diaggraph.fillComplete(params);
      TEST_EQUALITY_CONST( diaggraph.isFillComplete(), true );
      TEST_EQUALITY_CONST( diaggraph.isStorageOptimized(), true );
      TEST_EQUALITY_CONST( diaggraph.isUpperTriangular(), true );
      TEST_EQUALITY_CONST( diaggraph.isLowerTriangular(), true );
      // Bug verification:
      // Tpetra::CrsMatrix constructed with a Optimized, Fill-Complete graph will not call fillLocalMatrix()
      // in optimizeStorage(), because it returns early due to picking up the storage optimized bool from the graph.
      // As a result, the local mat-vec and mat-solve operations are never initialized, and localMultiply() and localSolve()
      // fail with a complaint regarding the initialization of these objects.
      MAT matrix(rcpFromRef(diaggraph));
      TEST_NOTHROW( matrix.setAllToScalar( ST::one() ) );
      matrix.fillComplete(params);
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
      TEST_EQUALITY_CONST( matrix.isStorageOptimized(), true );
      TEST_EQUALITY_CONST( matrix.isUpperTriangular(), true );
      TEST_EQUALITY_CONST( matrix.isLowerTriangular(), true );
      // init x to ones(); multiply into y, solve in-situ in y, check result
      V x(map,false), y(map,false);
      x.putScalar(SONE);
      TEST_NOTHROW( matrix.localMultiply(x,y,NO_TRANS,SONE,SZERO) );
      TEST_NOTHROW( matrix.localSolve(y,y,NO_TRANS) );
      ArrayRCP<const Scalar> x_view = x.get1dView();
      ArrayRCP<const Scalar> y_view = y.get1dView();
      if (ST::isOrdinal) {
        TEST_COMPARE_ARRAYS( y_view, x_view );
      } else {
        TEST_COMPARE_FLOATING_ARRAYS( y_view, x_view, MT::zero() );
      }
    }
    {
      // create a simple diagonal graph
      RCP<GRPH> diaggraph = rcp( new GRPH(map,1,StaticProfile) );
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        diaggraph->insertGlobalIndices(r,tuple(r));
      }
      diaggraph->fillComplete(params);
      TEST_EQUALITY_CONST( diaggraph->isFillComplete(), true );
      TEST_EQUALITY_CONST( diaggraph->isStorageOptimized(), true );
      TEST_EQUALITY_CONST( diaggraph->isUpperTriangular(), true );
      TEST_EQUALITY_CONST( diaggraph->isLowerTriangular(), true );
      // construct a matrix with the graph from another matrix
      MAT matrix1(diaggraph);
      TEST_EQUALITY( matrix1.getCrsGraph(), diaggraph );
      MAT matrix2( matrix1.getCrsGraph() );
      TEST_EQUALITY( matrix2.getCrsGraph(), matrix1.getCrsGraph() );
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, WithColMap, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // generate a tridiagonal matrix
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Scalar SONE  = ST::one();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map with numLocal entries per node using a pre-existing column map.
    // ensure:
    // * that the matrix uses this col map
    // * that it performs filtering during insertions
    // * that we can perform local or global insertions
    const size_t numLocal = 10; TEUCHOS_TEST_FOR_EXCEPTION( numLocal < 2, std::logic_error, "Test assumes that numLocal be greater than 1.");
    // these maps are equalivalent, but we should keep two distinct maps just to verify the general use case.
    RCP<const Map<LO,GO,Node> > rmap = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    RCP<const Map<LO,GO,Node> > cmap = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    //////////////////////////////////
    // add tridiagonal entries, but use a diagonal column map.
    // result should be block diagonal matrix, with no importer/exporter.
    //
    // run this test twice; once where we insert global indices and once where we insert local indices
    // both are allowed with a specified column map; however, we can only test one at a time.
    //
    // the first time, use const NNZ
    // the second, use NNZ array
    {
      MAT bdmat(rmap,cmap,3,StaticProfile);
      TEST_EQUALITY(bdmat.getRowMap(), rmap);
      TEST_EQUALITY_CONST(bdmat.hasColMap(), true);
      TEST_EQUALITY(bdmat.getColMap(), cmap);
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        // use global for the first one to verify that the matrix allows it
        // r-1 might be invalid, but the column map filtering should address that.
        bdmat.insertGlobalValues(r,tuple<GO>(r-1,r,r+1),tuple<Scalar>(SONE,SONE,SONE));
      }
      TEST_NOTHROW(bdmat.fillComplete());
      // nothing should have changed with regard to the row and column maps of the matrix
      TEST_EQUALITY(bdmat.getRowMap(), rmap);
      TEST_EQUALITY_CONST(bdmat.hasColMap(), true);
      TEST_EQUALITY(bdmat.getColMap(), cmap);
      // check that filtering happened
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        if (r == rmap->getMinGlobalIndex() || r == rmap->getMaxGlobalIndex()) {
          TEST_EQUALITY_CONST(bdmat.getNumEntriesInGlobalRow(r), 2);
        }
        else {
          TEST_EQUALITY_CONST(bdmat.getNumEntriesInGlobalRow(r), 3);
        }
      }
    }
    {
      ArrayRCP<size_t> nnzperrow = arcp<size_t>(numLocal);
      std::fill(nnzperrow.begin(), nnzperrow.end(), 3);
      MAT bdmat(rmap,cmap,nnzperrow,StaticProfile);
      TEST_EQUALITY(bdmat.getRowMap(), rmap);
      TEST_EQUALITY_CONST(bdmat.hasColMap(), true);
      TEST_EQUALITY(bdmat.getColMap(), cmap);
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        // use local for the rest. need the column map
        // column map and row map are the same, so we only have to do one translation
        LO lid = cmap->getLocalElement(r);
        // as above, filtering via column map (required to happen for local and global) will save us for the invalid r-1
        bdmat.insertLocalValues(lid,tuple<LO>(lid-1,lid,lid+1),tuple<Scalar>(SONE,SONE,SONE));
      }
      TEST_NOTHROW(bdmat.fillComplete());
      // nothing should have changed with regard to the row and column maps of the matrix
      TEST_EQUALITY(bdmat.getRowMap(), rmap);
      TEST_EQUALITY_CONST(bdmat.hasColMap(), true);
      TEST_EQUALITY(bdmat.getColMap(), cmap);
      // check that filtering happened
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        if (r == rmap->getMinGlobalIndex() || r == rmap->getMaxGlobalIndex()) {
          TEST_EQUALITY_CONST(bdmat.getNumEntriesInGlobalRow(r), 2);
        }
        else {
          TEST_EQUALITY_CONST(bdmat.getNumEntriesInGlobalRow(r), 3);
        }
      }
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, WithGraph_replaceLocal, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // generate a tridiagonal matrix
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Vector<Scalar,LO,GO,Node> V;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = size(*comm);
    // create a Map
    const size_t numLocal = 10;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    CrsGraph<LO,GO,Node> graph(map,3,StaticProfile);
    for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
      if (r == map->getMinAllGlobalIndex()) {
        graph.insertGlobalIndices(r,tuple(r,r+1));
      }
      else if (r == map->getMaxAllGlobalIndex()) {
        graph.insertGlobalIndices(r,tuple(r-1,r));
      }
      else {
        graph.insertGlobalIndices(r,tuple(r-1,r,r+1));
      }
    }
    graph.fillComplete();
    // create a matrix using the graph
    MAT matrix(rcpFromRef(graph));
    TEST_EQUALITY_CONST( matrix.getProfileType() == StaticProfile, true );
    // insert throws exception: not allowed with static graph
    TEST_THROW( matrix.insertGlobalValues(map->getMinGlobalIndex(),tuple<GO>(map->getMinGlobalIndex()),tuple(ST::one())), std::runtime_error );
    // suminto and replace are allowed
    for (LO r=map->getMinLocalIndex(); r <= map->getMaxLocalIndex(); ++r) {
      if (r == map->getMinLocalIndex()) {
        matrix.replaceLocalValues(r, tuple(r,r+1), tuple(ST::one(),ST::one()) );
      }
      else if (r == map->getMaxLocalIndex()) {
        matrix.replaceLocalValues(r, tuple(r-1,r), tuple(ST::one(),ST::one()) );
      }
      else {
        matrix.replaceLocalValues(r, tuple(r-1,r,r+1), tuple(ST::one(),ST::one(),ST::one()) );
      }
    }
    for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
      // increment the diagonals
      matrix.sumIntoGlobalValues(r, tuple(r), tuple(ST::one()) );
    }
    matrix.fillComplete();
    TEST_EQUALITY( matrix.getNodeNumDiags(), numLocal );
    TEST_EQUALITY( matrix.getGlobalNumDiags(), numImages*numLocal );
    TEST_EQUALITY( matrix.getGlobalNumEntries(), 3*numImages*numLocal - 2 );
    V dvec(map,false);
    dvec.randomize();
    matrix.getLocalDiagCopy(dvec);
    Array<Scalar> expectedDiags(numLocal, static_cast<Scalar>(2));
    ArrayRCP<const Scalar> dvec_view = dvec.get1dView();
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS( expectedDiags(), dvec_view );
    } else {
      TEST_COMPARE_FLOATING_ARRAYS( expectedDiags(), dvec_view, MT::zero() );
    }

    // Test the precomputed offsets version of getLocalDiagCopy().
    V dvec2 (map, false);
    dvec2.randomize ();
    ArrayRCP<size_t> offsets;
    matrix.getLocalDiagOffsets (offsets);
    TEST_EQUALITY( matrix.getNodeNumRows(), Teuchos::as<size_t>(offsets.size()) );
    matrix.getLocalDiagCopy (dvec2, offsets ());
    ArrayRCP<const Scalar> dvec2_view = dvec2.get1dView ();
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS( expectedDiags(), dvec2_view );
    } else {
      TEST_COMPARE_FLOATING_ARRAYS( expectedDiags(), dvec2_view, MT::zero() );
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ExceedStaticAlloc, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // test that an exception is thrown when we exceed statically allocated memory
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = size(*comm);
    // create a Map
    const size_t numLocal = 10;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    {
      MAT matrix(map,1,StaticProfile);
      // room for one on each row
      for (GO r=map->getMinGlobalIndex(); r<=map->getMaxGlobalIndex(); ++r)
      {
        matrix.insertGlobalValues(r,tuple(r),tuple(ST::one()));
      }
      // no room for any more
      GO r = map->getMinGlobalIndex();
      TEST_THROW( matrix.insertGlobalValues( r, tuple(r+1), tuple(ST::one()) ), std::runtime_error );
    }
    if (numImages > 1) {
      // add too many entries globally
      MAT matrix(map,1,StaticProfile);
      // room for one on each row
      for (GO r=map->getMinGlobalIndex(); r<=map->getMaxGlobalIndex(); ++r)
      {
        matrix.insertGlobalValues(r,tuple(r),tuple(ST::one()));
      }
      // always room for non-locals
      GO r = map->getMaxGlobalIndex() + 1;
      if (r > map->getMaxAllGlobalIndex()) r = map->getMinAllGlobalIndex();
      TEST_NOTHROW( matrix.insertGlobalValues( r, tuple(r), tuple(ST::one()) ) );
      // after communicating non-locals, failure trying to add them
      TEST_THROW( matrix.globalAssemble(), std::runtime_error );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, WithGraph, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, WithColMap, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, AdvancedGraphUsage, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, WithGraph_replaceLocal, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ExceedStaticAlloc, LO, GO, SCALAR, NODE ) \

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}
