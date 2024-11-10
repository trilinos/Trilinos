// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"

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

  using std::endl;

  using Teuchos::as;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::arrayView;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  //using Teuchos::VERB_MEDIUM;
  //using Teuchos::VERB_HIGH;
  //using Teuchos::VERB_EXTREME;
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
  using Tpetra::createContigMapWithNode;
  using Tpetra::createLocalMapWithNode;
  using Tpetra::createVector;
  using Tpetra::OptimizeOption;
  using Tpetra::DoOptimizeStorage;
  using Tpetra::DoNotOptimizeStorage;
  using Tpetra::GloballyDistributed;
  using Tpetra::INSERT;
  using namespace Tpetra::TestingUtilities;


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
    ArrayView<const GO> STMYGIDS = matrix.getRowMap()->getLocalElementList(); \
    typename MAT::local_inds_host_view_type loview; \
    typename MAT::values_host_view_type sview; \
    size_t STMAX = 0; \
    for (size_t STR=0; STR < matrix.getLocalNumRows(); ++STR) { \
      const size_t numEntries = matrix.getNumEntriesInLocalRow(STR); \
      TEST_EQUALITY( numEntries, matrix.getNumEntriesInGlobalRow( STMYGIDS[STR] ) ); \
      matrix.getLocalRowView(STR,loview,sview); \
      TEST_EQUALITY( static_cast<size_t>(loview.size()), numEntries ); \
      TEST_EQUALITY( static_cast<size_t>( sview.size()), numEntries ); \
      STMAX = std::max( STMAX, numEntries ); \
    } \
    TEST_EQUALITY( matrix.getLocalMaxNumRowEntries(), STMAX ); \
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
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, EmptyFillComplete, LO, GO, Scalar, Node )
  {
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef CrsGraph<LO,GO,Node>  GRPH;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    // create a Map with numLocal entries per process
    const size_t numLocal = 10;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    {
      // send in a parameterlist, check the defaults
      RCP<ParameterList> defparams = parameterList();
      // create static-profile matrix, fill-complete without inserting (and therefore, without allocating)
      MAT matrix(map,1);
      matrix.fillComplete(defparams);
      TEST_EQUALITY_CONST(defparams->get<bool>("Optimize Storage"), true);
    }

    {
      // send in a parameterlist, check the defaults
      RCP<ParameterList> defparams = parameterList();
      // create static-profile graph, fill-complete without inserting (and therefore, without allocating)
      GRPH graph(map,1);
      graph.fillComplete(defparams);
      TEST_EQUALITY_CONST(defparams->get<bool>("Optimize Storage"), true);
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, CopiesAndViews, LO, GO, Scalar, Node )
  {
    // test that an exception is thrown when we exceed statically allocated memory
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t numImages = size(*comm);
    const size_t myImageID = rank(*comm);
    if (numImages < 2) return;
    // create a Map, one row per processor
    const size_t numLocal = 1;
    RCP<const Map<LO,GO,Node> > rmap = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    GO myrowind = rmap->getGlobalElement(0);
    // specify the column map to control ordering
    // construct tridiagonal graph
    Array<GO> ginds;
    Array<LO> linds;
    if (myImageID==0) {
      tupleToArray( ginds, tuple<GO>(myrowind,myrowind+1) );
      tupleToArray( linds, tuple<LO>(0,1) );
    }
    else if (myImageID==numImages-1) {
      tupleToArray( ginds , tuple<GO>(myrowind-1,myrowind) );
      tupleToArray( linds , tuple<LO>(0,1) );
    }
    else {
      tupleToArray( ginds , tuple<GO>(myrowind-1,myrowind,myrowind+1) );
      tupleToArray( linds , tuple<LO>(0,1,2) );
    }
    Array<Scalar> vals(ginds.size(),ST::one());
    RCP<Map<LO,GO,Node> > cmap = rcp( new Map<LO,GO,Node>(INVALID,ginds(),0,comm) );
    RCP<ParameterList> params = parameterList();
    for (int T=0; T<4; ++T) {
      if ( (T & 1) != 1 ) continue;
      params->set("Optimize Storage",((T & 2) == 2));
      MAT matrix(rmap,cmap, ginds.size());   // only allocate as much room as necessary
      RowMatrix<Scalar,LO,GO,Node> &rowmatrix = matrix;

      typename MAT::nonconst_global_inds_host_view_type GCopy("gids",4);
      typename MAT::nonconst_local_inds_host_view_type LCopy("lids",4);
      typename MAT::nonconst_values_host_view_type SCopy("vals",4);

      typename MAT::nonconst_global_inds_host_view_type GCopy_toshort("gids",1);
      typename MAT::nonconst_local_inds_host_view_type LCopy_toshort("lids",1);
      typename MAT::nonconst_values_host_view_type SCopy_toshort("vals",1);

      typename MAT::global_inds_host_view_type CGView; 
      typename MAT::local_inds_host_view_type CLView; 
      typename MAT::values_host_view_type CSView;

      size_t numentries;
      // at this point, the graph has not allocated data as global or local, so we can do views/copies for either local or global
      matrix.getLocalRowCopy(0,LCopy,SCopy,numentries);
      matrix.getLocalRowView(0,CLView,CSView);

      matrix.getGlobalRowCopy(myrowind,GCopy,SCopy,numentries);
      matrix.getGlobalRowView(myrowind,CGView,CSView);

      // use multiple inserts: this illustrated an overwrite bug for column-map-specified graphs
      for (size_t j=0; j<(size_t)ginds.size(); ++j) {
        matrix.insertGlobalValues(myrowind,ginds(j,1),tuple(ST::one()));
      }
      TEST_EQUALITY( matrix.getNumEntriesInLocalRow(0), matrix.getCrsGraph()->getNumAllocatedEntriesInLocalRow(0) ); // test that we only allocated as much room as necessary
      // Before Mar 2019, insertion in to a matrix (graph) would append new indices/values
      // to the appropriate row. The result was that insertion would result in duplicate
      // indices (which would later be compressed out at fillComplete). A side-effect was
      // that insertion would throw if the total number of inserted indices exceeded the
      // space allocated - not just the number of unique indices. The current behavior now
      // checks if an index exists the graph and only inserts if it doesn't. Thus, the
      // following test needs to be modified to insert another *unique* index and not a
      // repeat.
      TEST_THROW( matrix.insertGlobalValues(myrowind, tuple(myrowind+5), tuple(ST::one())), std::runtime_error );
      matrix.fillComplete(params);
      // check for throws and no-throws/values
      TEST_THROW( matrix.getGlobalRowView(myrowind,CGView,CSView), std::runtime_error );
      TEST_THROW( matrix.getLocalRowCopy(0,LCopy_toshort,SCopy_toshort,numentries), std::runtime_error );
      TEST_THROW( matrix.getGlobalRowCopy(myrowind,GCopy_toshort,SCopy_toshort,numentries), std::runtime_error );

      //
      TEST_NOTHROW( matrix.getLocalRowView(0,CLView,CSView) );
      TEST_COMPARE_ARRAYS( CLView, linds );
      TEST_COMPARE_ARRAYS( CSView, vals  );
      //
      TEST_NOTHROW( matrix.getLocalRowCopy(0,LCopy,SCopy,numentries) );
      TEST_COMPARE_ARRAYS( arcp_from_view(LCopy,numentries), linds );
      TEST_COMPARE_ARRAYS( arcp_from_view(SCopy,numentries), vals );
      //
      TEST_NOTHROW( matrix.getGlobalRowCopy(myrowind,GCopy,SCopy,numentries) );
      TEST_COMPARE_ARRAYS( arcp_from_view(GCopy,numentries), ginds );
      TEST_COMPARE_ARRAYS( arcp_from_view(SCopy,numentries), vals  );
      //
      STD_TESTS(rowmatrix);
    }
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, Transpose, LO, GO, Scalar, Node )
  {
    using std::endl;
    out << "CrsMatrix Transpose test" << endl;
    Teuchos::OSTab tab0 (out);

    // this is the same matrix as in test NonSquare, but we will apply the transpose
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int M = 3;
    const int P = 5;
    const int N = comm->getSize();
    const int myImageID = comm->getRank();

    // create Maps
    // matrix is M*N-by-P
    //                  col
    //            0        1                  P-1
    //    0  [0        MN              ... (P-1)MN     ]
    //    .  [...      ...                 ...         ]
    //    0  [M-1      MN+M-1              (P-1)MN+M-1 ]
    //p   1  [M        MN+M                            ]
    //r   .  [...      ...                             ] = [A_ij], where A_ij = i+jMN
    //o   1  [2M-1     MN+2M-1                         ]
    //c   .  [...                                      ]
    //   N-1 [(N-1)M   MN+(N-1)(M-1)                   ]
    //    .  [...      ...                             ]
    //   N-1 [MN-1     MN+MN-1                         ]
    //
    // row map, range map is [0,M-1] [M,2M-1] [2M,3M-1] ... [MN-M,MN-1]
    // domain map will be a non-distributed map for a vector of length P
    //
    // input multivector is
    //  col      0            1         ...        numVecs-1
    //     0 [0         MN                    (numVecs-1)MN     ]
    // p   . [...       ...                   ...               ]
    // r   0 [M-1       MN+M-1                (numVecs-1)MN+M-1 ]
    // o   1 [M         MN+M                                    ]
    // c   . [...       ...                                     ] = [X_ij], where X_ij = i+jMN
    //     1 [2M-1      MN+2M-1                                 ]
    //     . [...       ...                                     ]
    //    N-1[(N-1)M    MN+(N-1)(M-1)                           ]
    //     . [...       ...                                     ]
    //    N-1[MN-1      MN+MN-1                                 ]
    //
    // output multivector is not-distributed
    // the result of the transpose multiplication should be
    //              MN-1              MN-1
    // (A^T*X)_ij = sum_k A_ki X_kj = sum (k+iMN)(k+jMN)
    //              k=0               k=0
    //   MN-1
    // = sum k(i+j)MN + ij(MN)(MN) + k^2 = (i+j)(MN)^2(MN-1)/2 + ij(MN)^3 + (MN)(MN-1)(2MN-1)/6
    //   k=0
    //

    out << "Construct Maps" << endl;
    const int numVecs = 3;
    RCP<const Map<LO,GO,Node> > rowmap = createContigMapWithNode<LO,GO,Node>(INVALID,M,comm);
    RCP<const Map<LO,GO,Node> > lclmap = createLocalMapWithNode<LO,GO,Node>(P,comm);

    // create the matrix
    out << "Create matrix" << endl;
    MAT A (rowmap, P);
    for (int i = 0; i < M; ++i) {
      for (int j=0; j < P; ++j) {
        A.insertGlobalValues (static_cast<GO> (M*myImageID+i), tuple<GO> (j), tuple<Scalar> (M*myImageID+i + j*M*N));
      }
    }

    // call fillComplete()
    out << "Call fillComplete on the matrix" << endl;
    A.fillComplete (lclmap, rowmap);
    A.describe (out, VERB_LOW);

    out << "Build input and output multivectors" << endl;

    // build the input multivector X
    MV X (rowmap, numVecs);
    for (int i = myImageID*M; i < myImageID*M+M; ++i) {
      for (int j = 0; j < numVecs; ++j) {
        X.replaceGlobalValue (i, j, static_cast<Scalar> (i + j*M*N));
      }
    }

    // build the expected output multivector B
    MV Bexp (lclmap, numVecs);
    MV Bout (lclmap, numVecs);
    for (int i = 0; i < P; ++i) {
      for (int j = 0; j < numVecs; ++j) {
        Bexp.replaceGlobalValue (i, j, static_cast<Scalar> ((i+j)*(M*N)*(M*N)*(M*N-1)/2 + i*j*(M*N)*(M*N)*(M*N) + (M*N)*(M*N-1)*(2*M*N-1)/6));
      }
    }

    // test the action
    out << "Test applying the (conjugate) transpose of the matrix" << endl;
    Bout.randomize ();
    X.setObjectLabel ("Input");
    //X.describe(out, VERB_EXTREME);
    A.apply (X, Bout, CONJ_TRANS);
    Bout.setObjectLabel("Actual output");
    Bexp.setObjectLabel("Expected output");
    //Bout.describe(out, VERB_EXTREME);
    //Bexp.describe(out, VERB_EXTREME);
    Bout.update (-ST::one (), Bexp, ST::one());
    Array<Mag> norms (numVecs);
    Array<Mag> zeros (numVecs, MT::zero ());
    Bout.norm1 (norms ());
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms,zeros);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
    }

    out << "Done with test" << endl;
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, Transpose,      LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, EmptyFillComplete, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, CopiesAndViews, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}
