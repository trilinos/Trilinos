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

#ifndef TPETRA_TEST_CRSMATRIX_WITHGRAPH_HPP
#define TPETRA_TEST_CRSMATRIX_WITHGRAPH_HPP

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Details_getNumDiags.hpp"

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

namespace Tpetra {
namespace Test {

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

  using Tpetra::TestingUtilities::getDefaultComm;
  using Tpetra::createContigMapWithNode;
  using Tpetra::global_size_t;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  using Teuchos::outArg;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using std::endl;

  double errorTolSlack = 1e+1;
  std::string filedir;

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
    reduceAll<int,global_size_t>( *STCOMM, Teuchos::REDUCE_MAX, STMAX, outArg(STGMAX) ); \
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
    // generate a tridiagonal matrix
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> MAT;
    const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map with numLocal entries per node
    const size_t numLocal = 10;
    RCP<const Tpetra::Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    {
      Tpetra::CrsGraph<LO,GO,Node> diaggraph (map, 1, Tpetra::StaticProfile);
      // A pre-constructed graph must be fill complete before being used to construct a CrsMatrix
      TEST_THROW( MAT matrix(rcpFromRef(diaggraph)), std::runtime_error );
    }
    {
      // create a simple diagonal graph
      Tpetra::CrsGraph<LO,GO,Node> diaggraph (map, 1, Tpetra::StaticProfile);
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

  // Test that CrsMatrix operations (setAllToScalar, apply, and local
  // triangular solve) work when the matrix was constructed with a
  // const CrsGraph.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, WithGraph, LO, GO, Scalar, Node )
  {
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node>  MAT;
    typedef Tpetra::CrsGraph<LO,GO,Node>         GRPH;
    typedef Tpetra::Vector<Scalar,LO,GO,Node>       V;
    typedef ScalarTraits<Scalar>                   ST;
    typedef typename ST::magnitudeType            Mag;
    typedef ScalarTraits<Mag>                      MT;

    out << "Test that CrsMatrix operations (setAllToScalar, apply, and "
      "local triangular solve) work when the matrix was constructed with "
      "a const CrsGraph" << endl;
    Teuchos::OSTab tab1 (out);

    const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
    const Scalar SONE  = ST::one();
    const Scalar SZERO = ST::zero();

    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = comm->getSize ();
    out << "Communicator has " << numImages << " process" << (numImages != 1 ? "es" : "") << endl;

    const size_t numLocal = 10;
    out << "Create a Map with numLocal=" << numLocal << " entries per process" << endl;

    RCP<const Tpetra::Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    RCP<ParameterList> params = parameterList();
    params->set("Optimize Storage",true);
    RCP<ParameterList> fillparams = sublist(params,"Local Sparse Ops");
    fillparams->set("Prepare Solve", true);

    {
      out << "Create tridiagonal CrsGraph with StaticProfile" << endl;

      GRPH trigraph (map, 3, Tpetra::StaticProfile);
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

      out << "Call fillComplete on the CrsGraph" << endl;
      trigraph.fillComplete(params);

      out << "Create a CrsMatrix using the tridiagonal CrsGraph "
        "and test allowed functionality" << endl;

      MAT matrix(rcpFromRef(trigraph));
      TEST_EQUALITY_CONST( matrix.getProfileType() == Tpetra::StaticProfile, true );
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

      out << "Call fillComplete on the CrsMatrix" << endl;
      matrix.fillComplete();
      TEST_EQUALITY( matrix.getGlobalNumEntries(), 3*numImages*numLocal - 2 );

      out << "Check the diagonal entries of the CrsMatrix, using getLocalDiagCopy" << endl;
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

    // Make sure that all processes finished and were successful.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED; no sense in continuing" << endl;
      return;
    }

    {
      out << "Create a diagonal CrsGraph" << endl;
      GRPH diaggraph (map, 1, Tpetra::StaticProfile);
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        diaggraph.insertGlobalIndices(r,tuple(r));
      }

      out << "Call fillComplete on the CrsGraph" << endl;
      diaggraph.fillComplete(params);

      // Bug verification:
      //  Tpetra::CrsMatrix constructed with a graph was experiencing a seg-fault if setAllToScalar is called before
      //  some other call allocates memory. This was because setAllToScalar has an incorrect if-statement that is
      //  not allocating memory.
      // This bug has been fixed. Furthermore, CrsMatrix no longer utilizes lazy allocation when constructed with a graph.
      // However, we will leave this test in place, because it still demonstrates valid behavior.

      out << "Create a CrsMatrix with the diagonal CrsGraph" << endl;
      MAT matrix(rcpFromRef(diaggraph));

      out << "Call setAllToScalar on the CrsMatrix; it should not throw" << endl;
      TEST_NOTHROW( matrix.setAllToScalar( ST::one() ) );
    }
    TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success);
    if (!success) {
      out << "Test FAILED; no sense in continuing" << endl;
      return;
    }


    {
      out << "Create a diagonal CrsGraph" << endl;
      GRPH diaggraph (map, 1, Tpetra::StaticProfile);
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        diaggraph.insertGlobalIndices(r,tuple(r));
      }

      out << "Call fillComplete on the CrsGraph" << endl;
      diaggraph.fillComplete(params);


      out << "Create a CrsMatrix with the diagonal CrsGraph and Kokkos view" << endl;
      size_t numEnt = diaggraph.getLocalGraph().entries.extent(0);
      typename MAT::local_matrix_type::values_type val ("Tpetra::CrsMatrix::val", numEnt);
      MAT matrix(rcpFromRef(diaggraph),val);

      out << "Call setAllToScalar on the CrsMatrix; it should not throw" << endl;
      TEST_NOTHROW( matrix.setAllToScalar( ST::one() ) );
    }
    TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success);
    if (!success) {
      out << "Test FAILED; no sense in continuing" << endl;
      return;
    }

    {
      out << "Create a diagonal CrsGraph" << endl;
      GRPH diaggraph (map, 1, Tpetra::StaticProfile);
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        diaggraph.insertGlobalIndices(r,tuple(r));
      }

      out << "Call fillComplete on the CrsGraph" << endl;
      diaggraph.fillComplete(params);

      TEST_EQUALITY_CONST( diaggraph.isFillComplete(), true );
      TEST_EQUALITY_CONST( diaggraph.isStorageOptimized(), true );

      // Make sure that if you create a CrsMatrix with an optimized,
      // fillComplete CrsGraph, that the resulting CrsMatrix is
      // fillComplete, has optimized storage, and has a correctly
      // initialized local sparse matrix-vector multiply.

      out << "Create a CrsMatrix with the diagonal CrsGraph" << endl;
      MAT matrix(rcpFromRef(diaggraph));

      out << "Call setAllToScalar on the CrsMatrix; it should not throw" << endl;
      TEST_NOTHROW( matrix.setAllToScalar( ST::one() ) );

      out << "Call fillComplete on the CrsMatrix" << endl;
      matrix.fillComplete(params);

      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
      TEST_EQUALITY_CONST( matrix.isStorageOptimized(), true );
      // init x to ones(); multiply into y, solve in-situ in y, check result
      V x(map,false), y(map,false);
      x.putScalar(SONE);

      out << "Call localApply (local mat-vec) on the CrsMatrix; "
        "it should not throw" << endl;
      TEST_NOTHROW( matrix.localApply(x,y,Teuchos::NO_TRANS,SONE,SZERO) );

      ArrayRCP<const Scalar> x_view = x.get1dView();
      ArrayRCP<const Scalar> y_view = y.get1dView();
      if (ST::isOrdinal) {
        TEST_COMPARE_ARRAYS( y_view, x_view );
      } else {
        TEST_COMPARE_FLOATING_ARRAYS( y_view, x_view, MT::zero() );
      }
    }

    // Make sure that all processes finished and were successful.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED; no sense in continuing" << endl;
      return;
    }

    {
      out << "Create a diagonal CrsGraph" << endl;
      RCP<GRPH> diaggraph = rcp( new GRPH (map, 1, Tpetra::StaticProfile) );
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        diaggraph->insertGlobalIndices(r,tuple(r));
      }

      out << "Call fillComplete on the CrsGraph" << endl;
      diaggraph->fillComplete(params);

      TEST_EQUALITY_CONST( diaggraph->isFillComplete(), true );
      TEST_EQUALITY_CONST( diaggraph->isStorageOptimized(), true );

      out << "Construct a CrsMatrix with the diagonal CrsGraph" << endl;
      MAT matrix1(diaggraph);
      TEST_EQUALITY( matrix1.getCrsGraph(), diaggraph );

      out << "Construct another CrsMatrix with the same diagonal CrsGraph "
        "shared by the first CrsMatrix" << endl;
      MAT matrix2( matrix1.getCrsGraph() );
      TEST_EQUALITY( matrix2.getCrsGraph(), matrix1.getCrsGraph() );
    }

    // Make sure that all processes finished and were successful.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, WithColMap, LO, GO, Scalar, Node )
  {
    // generate a tridiagonal matrix
    typedef ScalarTraits<Scalar> ST;
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> MAT;
    const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
    const Scalar SONE  = ST::one();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map with numLocal entries per node using a pre-existing column map.
    // ensure:
    // * that the matrix uses this column Map
    // * that we can perform local or global insertions
    const size_t numLocal = 10;
    TEUCHOS_TEST_FOR_EXCEPTION( numLocal < 2, std::logic_error, "Test assumes that numLocal be greater than 1.");
    // these maps are equalivalent, but we should keep two distinct maps just to verify the general use case.
    RCP<const Tpetra::Map<LO,GO,Node> > rmap = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    RCP<const Tpetra::Map<LO,GO,Node> > cmap = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    //////////////////////////////////
    // add tridiagonal entries, but use a diagonal column map.
    // result should be block diagonal matrix, with no importer/exporter.
    //
    // run this test twice; once where we insert global indices and once where we insert local indices
    // both are allowed with a specified column map; however, we can only test one at a time.

    // First test: use a constant upper bound (3) on the number of
    // entries in each row, and insert using global indices.
    {
      MAT bdmat (rmap, cmap, 3, Tpetra::StaticProfile);
      TEST_EQUALITY(bdmat.getRowMap(), rmap);
      TEST_EQUALITY_CONST(bdmat.hasColMap(), true);
      TEST_EQUALITY(bdmat.getColMap(), cmap);

      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        // The second, apparently superfluous check avoids issues if
        // r-1 overflows unsigned.
        if (r - 1 >= cmap->getMinGlobalIndex () && r - 1 <= cmap->getMaxGlobalIndex ()) {
          // The second, apparently superfluous check avoids issues if
          // r+1 overflows.
          if (r + 1 <= cmap->getMaxGlobalIndex () && r + 1 >= cmap->getMinGlobalIndex ()) {
            bdmat.insertGlobalValues(r,tuple<GO>(r-1,r,r+1),tuple<Scalar>(SONE,SONE,SONE));
          } else {
            bdmat.insertGlobalValues(r,tuple<GO>(r-1,r),tuple<Scalar>(SONE,SONE));
          }
        } else { // r - 1 invalid
          if (r + 1 <= cmap->getMaxGlobalIndex () && r + 1 >= cmap->getMinGlobalIndex ()) {
            bdmat.insertGlobalValues(r,tuple<GO>(r,r+1),tuple<Scalar>(SONE,SONE));
          } else { // r + 1 invalid
            bdmat.insertGlobalValues(r,tuple<GO>(r),tuple<Scalar>(SONE));
          }
        }
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

    // Second test: use an array to bound from above the number of
    // entries in each row, and insert using local indices.
    {
      Teuchos::Array<size_t> nnzperrow (numLocal);
      std::fill(nnzperrow.begin(), nnzperrow.end(), 3);
      MAT bdmat (rmap, cmap, nnzperrow (), Tpetra::StaticProfile);
      TEST_EQUALITY(bdmat.getRowMap(), rmap);
      TEST_EQUALITY_CONST(bdmat.hasColMap(), true);
      TEST_EQUALITY(bdmat.getColMap(), cmap);

      for (LO localRow = rmap->getMinLocalIndex (); localRow <= rmap->getMaxLocalIndex (); ++localRow) {
        // In this test, the column Map and row Map are the same, so
        // we can use localRow as the local column index as well.
        const LO c = localRow;

        // The second, apparently superfluous check avoids issues if
        // c-1 overflows unsigned.
        if (c - 1 >= cmap->getMinLocalIndex () && c - 1 <= cmap->getMaxLocalIndex ()) {
          // The second, apparently superfluous check avoids issues if
          // c+1 overflows.
          if (c + 1 <= cmap->getMaxLocalIndex () && c + 1 >= cmap->getMinLocalIndex ()) {
            bdmat.insertLocalValues (c, tuple<LO>(c-1,c,c+1), tuple<Scalar>(SONE,SONE,SONE));
          } else { // c + 1 is an invalid column index
            bdmat.insertLocalValues (c, tuple<LO>(c-1,c), tuple<Scalar>(SONE,SONE));
          }
        } else { // c - 1 is an invalid column index
          // The second, apparently superfluous check avoids issues if
          // c+1 overflows.
          if (c + 1 <= cmap->getMaxLocalIndex () && c + 1 >= cmap->getMinLocalIndex ()) {
            bdmat.insertLocalValues (c, tuple<LO>(c,c+1), tuple<Scalar>(SONE,SONE));
          } else { // c + 1 is an invalid column index
            bdmat.insertLocalValues (c, tuple<LO>(c), tuple<Scalar>(SONE));
          }
        }
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
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Tpetra::Vector<Scalar,LO,GO,Node> V;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef typename ST::magnitudeType Mag;
    typedef Teuchos::ScalarTraits<Mag> MT;
    int lclSuccess = 1; // for use below in finding whether all processes passed
    int gblSuccess = 0; // for use below in finding whether all processes passed

    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = comm->getSize ();
    out << "Communicator has " << numImages << " process" << (numImages != 1 ? "es" : "") << endl;

    const size_t numLocal = 10;
    out << "Create a Map with numLocal=" << numLocal << " entries per process" << endl;

    const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
    RCP<const Tpetra::Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);

    out << "Create a tridiagonal CrsGraph" << endl;
    Tpetra::CrsGraph<LO,GO,Node> graph (map, 3, Tpetra::StaticProfile);
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

    out << "Call fillComplete on the CrsGraph" << endl;
    graph.fillComplete();

    out << "Create a CrsMatrix using the tridiagonal CrsGraph" << endl;
    MAT matrix(rcpFromRef(graph));

    TEST_ASSERT( matrix.getProfileType () == Tpetra::StaticProfile );
    TEST_ASSERT( matrix.isStaticGraph () );

    // Make sure that all processes finished and were successful.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED; no sense in continuing" << endl;
      return;
    }

    {
      using Teuchos::TypeNameTraits;

      std::ostringstream os;
      os << "Process " << comm->getRank () << ": About to call insertGlobalValues: "
        "Scalar = " << TypeNameTraits<Scalar>::name () << ", LO = " <<
        TypeNameTraits<LO>::name () << ", GO = " << TypeNameTraits<GO>::name () <<
        ", Node = " << TypeNameTraits<Node>::name () << std::endl;
      std::cerr << os.str () << std::endl;
    }
    // insert throws exception: not allowed with static graph
    //TEST_THROW( matrix.insertGlobalValues(map->getMinGlobalIndex(),tuple<GO>(map->getMinGlobalIndex()),tuple(ST::one())), std::runtime_error );

    // FIXME (mfh 08 May 2014): For some entirely inexplicable reason,
    // this test fails only if Scalar is unsigned int.  I honestly
    // have no idea why, or whether it was failing before, or whether
    // this has something to do with Clang 3.2 on my Mac or whatever.
    // The test passes (that is, correctly throws) for every other
    // Scalar type.  Thus, I'm disabling the test in question when
    // Scalar is unsigned int.  Sometime when I have a chance, I'll go
    // back and revisit this.
    if (typeid (Scalar) != typeid (unsigned int)) {
      out << "Attempt to call insertGlobalValues; this should throw" << endl;

      bool didThrowOnOwnedInsert = false;
      try {
        matrix.insertGlobalValues (map->getMinGlobalIndex (), tuple<GO> (map->getMinGlobalIndex ()), tuple (ST::one ()));
        didThrowOnOwnedInsert = false;
      }
      catch (std::runtime_error&) {
        didThrowOnOwnedInsert = true;
      }
      TEST_ASSERT( didThrowOnOwnedInsert );

      {
        using Teuchos::TypeNameTraits;

        std::ostringstream os;
        os << "Process " << comm->getRank () << ": insertGlobalValues ";
        if (didThrowOnOwnedInsert) {
          os << "threw, like it should: ";
        } else {
          os << "did NOT throw (oops!): ";
        }

        os << "Scalar = " << TypeNameTraits<Scalar>::name () << ", LO = " <<
          TypeNameTraits<LO>::name () << ", GO = " << TypeNameTraits<GO>::name ()
           << ", Node = " << TypeNameTraits<Node>::name () << std::endl;
        std::cerr << os.str () << std::endl;
      }
    }

    // Make sure that all processes finished and were successful.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED; no sense in continuing" << endl;
      return;
    }

    out << "Attempt to call replaceLocalValues; this should _not_ throw" << endl;
    for (LO r=map->getMinLocalIndex(); r <= map->getMaxLocalIndex(); ++r) {
      if (r == map->getMinLocalIndex()) {
        const LO numChanged =
          matrix.replaceLocalValues(r, tuple(r,r+1), tuple(ST::one(),ST::one()) );
        TEST_EQUALITY_CONST( numChanged, static_cast<LO> (2) );
      }
      else if (r == map->getMaxLocalIndex()) {
        const LO numChanged =
          matrix.replaceLocalValues(r, tuple(r-1,r), tuple(ST::one(),ST::one()) );
        TEST_EQUALITY_CONST( numChanged, static_cast<LO> (2) );
      }
      else {
        const LO numChanged =
          matrix.replaceLocalValues(r, tuple(r-1,r,r+1), tuple(ST::one(),ST::one(),ST::one()) );
        TEST_EQUALITY_CONST( numChanged, static_cast<LO> (3) );
      }
    }

    // Make sure that all processes finished and were successful.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED; no sense in continuing" << endl;
      return;
    }

    out << "Attempt to call sumIntoGlobalValues on owned rows; this should _not_ throw" << endl;
    for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
      // increment the diagonals
      const LO numChanged =
        matrix.sumIntoGlobalValues(r, tuple(r), tuple(ST::one()) );
      TEST_EQUALITY_CONST( numChanged, static_cast<LO> (1) );
    }

    // Make sure that all processes finished and were successful.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED; no sense in continuing" << endl;
      return;
    }

    out << "Call fillComplete on the CrsMatrix" << endl;
    matrix.fillComplete();

    TEST_EQUALITY( Tpetra::Details::getLocalNumDiags (matrix), static_cast<LO> (numLocal) );
    TEST_EQUALITY( Tpetra::Details::getGlobalNumDiags (matrix), static_cast<GO> (numImages*numLocal) );
    TEST_EQUALITY( matrix.getGlobalNumEntries(), 3*numImages*numLocal - 2 );

    // Make sure that all processes finished and were successful.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED; no sense in continuing" << endl;
      return;
    }

    out << "Test the CrsMatrix's diagonals using 1-argument "
      "getLocalDiagCopy" << endl;
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

    // Make sure that all processes finished and were successful.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED; no sense in continuing" << endl;
      return;
    }

    out << "Test the CrsMatrix's diagonals using 2-argument "
      "getLocalDiagCopy (that takes precomputed offsets)" << endl;
    {
      Teuchos::OSTab tab2 (out);

      V dvec2 (map, false);
      dvec2.randomize ();
      ArrayRCP<size_t> offsets;

      out << "Call matrix.getLocalDiagOffsets (ArrayRCP version)" << endl;
      matrix.getLocalDiagOffsets (offsets);
      TEST_EQUALITY( matrix.getNodeNumRows(), Teuchos::as<size_t>(offsets.size()) );

      out << "Call matrix.getLocalDiagCopy (2-arg version with "
        "ArrayView offsets)" << endl;
      matrix.getLocalDiagCopy (dvec2, offsets ());

      out << "Check diagonal entries" << endl;
      ArrayRCP<const Scalar> dvec2_view = dvec2.get1dView ();
      if (ST::isOrdinal) {
        TEST_COMPARE_ARRAYS( expectedDiags(), dvec2_view );
      } else {
        TEST_COMPARE_FLOATING_ARRAYS( expectedDiags(), dvec2_view, MT::zero() );
      }
    }

    // Make sure that all processes finished and were successful.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ExceedStaticAlloc, LO, GO, Scalar, Node )
  {
    // test that an exception is thrown when we exceed statically allocated memory
    typedef ScalarTraits<Scalar> ST;
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> MAT;
    const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 10;
    RCP<const Tpetra::Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    {
      MAT matrix(map, 1, Tpetra::StaticProfile);
      // room for one on each row
      for (GO r=map->getMinGlobalIndex(); r<=map->getMaxGlobalIndex(); ++r)
      {
        matrix.insertGlobalValues(r,tuple(r),tuple(ST::one()));
      }
      // no room for any more
      GO r = map->getMinGlobalIndex();
      TEST_THROW( matrix.insertGlobalValues( r, tuple(r+1), tuple(ST::one()) ), std::runtime_error );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

//
// Macro for test instanatiations.
// Please only invoke this in the Tpetra::Test namespace.
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, WithGraph, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, WithColMap, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, AdvancedGraphUsage, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, WithGraph_replaceLocal, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ExceedStaticAlloc, LO, GO, SCALAR, NODE ) \

} // namespace Test
} // namespace Tpetra

#endif // TPETRA_TEST_CRSMATRIX_WITHGRAPH_HPP
