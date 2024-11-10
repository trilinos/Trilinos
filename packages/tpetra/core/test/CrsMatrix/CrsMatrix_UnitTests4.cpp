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
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Tpetra_Details_getNumDiags.hpp"
#include "Tpetra_Details_extractBlockDiagonal.hpp"
#include "Tpetra_Details_scaleBlockDiagonal.hpp"
#include "Tpetra_Apply_Helpers.hpp"
#include "TpetraUtils_MatrixGenerator.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include <type_traits> // std::is_same


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
  using Teuchos::arcp;
  using Teuchos::outArg;
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
  using Tpetra::createContigMapWithNode;

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
    ArrayView<const LO> loview; \
    ArrayView<const Scalar> sview; \
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

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, TheEyeOfTruthDistAlloc, LO, GO, Scalar, Node )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();

    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.randomize();

    // create the identity matrix
    RCP<RowMatrix<Scalar,LO,GO,Node> > eye;
    {
      RCP<MAT> eye_crs = rcp(new MAT(map,1) );
      if (myImageID == 0) {
        for (int i=0; i<map->getGlobalNumEntries(); ++i) {
          eye_crs->insertGlobalValues(i,tuple<GO>(i),tuple<Scalar>(ST::one()));
        }
      }
      eye_crs->fillComplete();
      eye = eye_crs;
    }
    // test the properties
    TEST_EQUALITY(eye->getGlobalNumEntries()  , numImages*numLocal);
    TEST_EQUALITY(eye->getLocalNumEntries()      , numLocal);
    TEST_EQUALITY(eye->getGlobalNumRows()      , numImages*numLocal);
    TEST_EQUALITY(eye->getLocalNumRows()          , numLocal);
    TEST_EQUALITY(eye->getLocalNumCols()          , numLocal);
    TEST_EQUALITY( Tpetra::Details::getGlobalNumDiags (*eye), static_cast<GO> (numImages*numLocal) );
    TEST_EQUALITY( Tpetra::Details::getLocalNumDiags (*eye), static_cast<LO> (numLocal) );
    TEST_EQUALITY(eye->getGlobalMaxNumRowEntries(), 1);
    TEST_EQUALITY(eye->getLocalMaxNumRowEntries()    , 1);
    TEST_EQUALITY(eye->getIndexBase()          , 0);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getColMap())   , true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getDomainMap()), true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getRangeMap()) , true);
    // test the action
    mvres.randomize();
    eye->apply(mvrand,mvres);
    mvres.update(-ST::one(),mvrand,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm1(norms());
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms,zeros);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, AlphaBetaMultiply, LO, GO, Scalar, Node )
  {
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef  Operator<Scalar,LO,GO,Node> OP;
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    const size_t THREE = 3;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t myImageID = comm->getRank();
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,THREE,comm);

    /* Create the identity matrix, three rows per proc */
    RCP<OP> AOp;

    RCP<MAT> A = rcp(new MAT(map,1));
    A->insertGlobalValues(3*myImageID,  tuple<GO>(3*myImageID  ), tuple<Scalar>(ST::one()) );
    A->insertGlobalValues(3*myImageID+1,tuple<GO>(3*myImageID+1), tuple<Scalar>(ST::one()) );
    A->insertGlobalValues(3*myImageID+2,tuple<GO>(3*myImageID+2), tuple<Scalar>(ST::one()) );
    A->fillComplete();
    AOp = A;

    MV X(map,1), Y(map,1), Z(map,1);
    const Scalar alpha = ST::random(),
                  beta = ST::random();
    X.randomize();
    Y.randomize();
    // Keep copies for later testing of CrsMatrixMultiplyOp
    MV X_copy (X, Teuchos::Copy);
    MV Y_copy (Y, Teuchos::Copy);

    // Z = alpha*X + beta*Y
    Z.update(alpha,X,beta,Y,ST::zero());
    // test the action: Y = alpha*I*X + beta*Y = alpha*X + beta*Y = Z
    AOp->apply(X,Y,NO_TRANS,alpha,beta);

    // mfh 07 Dec 2018: Little test for CrsMatrixMultiplyOp; it
    // doesn't get tested much elsewhere.  (It used to be part of
    // CrsMatrix's implementation, so it got more exercise before.)
    Tpetra::CrsMatrixMultiplyOp<Scalar, Scalar, LO, GO, Node> multOp (A);
    multOp.apply (X_copy, Y_copy, NO_TRANS, alpha, beta);

    Array<Mag> normY(1), normZ(1);
    Z.norm1(normZ());
    Y.norm1(normY());
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(normY,normZ);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(normY,normZ,2.0*testingTol<Mag>());
    }

    Array<Mag> normYcopy(1);
    Y_copy.norm1 (normYcopy ());
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(normYcopy,normZ);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(normYcopy,normZ,2.0*testingTol<Mag>());
    }
  }

  // Make sure that CrsMatrix and RowMatrix have the correct typedefs,
  // and that the typedefs match up with their corresponding template
  // parameters.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, Typedefs, LO, GO, Scalar, Node )
  {
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename MAT::scalar_type         scalar_type;
    typedef typename MAT::local_ordinal_type  local_ordinal_type;
    typedef typename MAT::global_ordinal_type global_ordinal_type;
    typedef typename MAT::node_type           node_type;
    static_assert (std::is_same<scalar_type, Scalar>::value,
                   "CrsMatrix<Scalar, ...>::scalar_type != Scalar");
    static_assert (std::is_same<local_ordinal_type, LO>::value,
                   "CrsMatrix<Scalar, LO, ...>::local_ordinal_type != LO");
    static_assert (std::is_same<global_ordinal_type, GO>::value,
                   "CrsMatrix<Scalar, LO, GO, ...>::global_ordinal_type != GO");
    static_assert (std::is_same<node_type, Node>::value,
                   "CrsMatrix<Scalar, LO, GO, Node>::node_type != Node");

    typedef RowMatrix<Scalar,LO,GO,Node> RMAT;
    typedef typename RMAT::scalar_type         rmat_scalar_type;
    typedef typename RMAT::local_ordinal_type  rmat_local_ordinal_type;
    typedef typename RMAT::global_ordinal_type rmat_global_ordinal_type;
    typedef typename RMAT::node_type           rmat_node_type;
    static_assert (std::is_same<rmat_scalar_type, Scalar>::value,
                   "RowMatrix<Scalar, ...>::scalar_type != Scalar");
    static_assert (std::is_same<rmat_local_ordinal_type, LO>::value,
                   "RowMatrix<Scalar, LO, ...>::local_ordinal_type != LO");
    static_assert (std::is_same<rmat_global_ordinal_type, GO>::value,
                   "RowMatrix<Scalar, LO, GO, ...>::global_ordinal_type != GO");
    static_assert (std::is_same<rmat_node_type, Node>::value,
                   "RowMatrix<Scalar, LO, GO, Node>::node_type != Node");
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ActiveFill, LO, GO, Scalar, Node )
  {
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,1,comm);
    const Scalar SZERO = ScalarTraits<Scalar>::zero();
    {
      MAT matrix(map,map,1);
      TEST_EQUALITY_CONST( matrix.isFillActive(),   true );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), false );
      matrix.insertLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) );
      //
      RCP<ParameterList> params = parameterList();
      params->set("Optimize Storage",false);
      matrix.fillComplete(params);
      TEST_EQUALITY_CONST( matrix.isFillActive(),   false );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );

      LO numValid = 0;

      // It's forbidden to call any of the *LocalValues methods if the
      // matrix is fill complete (not fill active).

      TEST_THROW( matrix.insertLocalValues ( 0, tuple<LO>(0), tuple<Scalar>(0) ), std::runtime_error );
      // numValid = matrix.insertLocalValues (0, tuple<LO> (0), tuple<Scalar> (0));
      // TEST_EQUALITY( numValid, 0 );

      numValid = matrix.replaceLocalValues (0, tuple<LO> (0), tuple<Scalar> (0));
      TEST_EQUALITY( numValid, Teuchos::OrdinalTraits<LO>::invalid () );

      numValid = matrix.replaceLocalValues (0, tuple<LO> (0), tuple<Scalar> (0));
      TEST_EQUALITY( numValid, Teuchos::OrdinalTraits<LO>::invalid () );

      numValid = matrix.sumIntoLocalValues (0, tuple<LO> (0), tuple<Scalar> (0));
      TEST_EQUALITY( numValid, Teuchos::OrdinalTraits<LO>::invalid () );

      TEST_THROW( matrix.globalAssemble(),      std::runtime_error );
      TEST_THROW( matrix.fillComplete(),        std::runtime_error );
    }
    {
      MAT matrix(map,map,1);
      TEST_EQUALITY_CONST( matrix.isFillActive(),   true );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), false );
      matrix.insertLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) );

      RCP<ParameterList> params = parameterList();
      params->set("Optimize Storage",false);
      matrix.fillComplete(params);
      TEST_EQUALITY_CONST( matrix.isFillActive(),   false );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );

      matrix.resumeFill();
      TEST_EQUALITY_CONST( matrix.isFillActive(),   true );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), false );
      TEST_NOTHROW( matrix.insertLocalValues ( 0, tuple<LO>(0), tuple<Scalar>(0) ) );
      TEST_NOTHROW( matrix.replaceLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) ) );
      TEST_NOTHROW( matrix.sumIntoLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) ) );
      TEST_NOTHROW( matrix.setAllToScalar(SZERO)                                   );
      TEST_NOTHROW( matrix.scale(SZERO)                                            );
      TEST_NOTHROW( matrix.globalAssemble()                                        );

      TEST_NOTHROW( matrix.fillComplete()                        );
      TEST_EQUALITY_CONST( matrix.isFillActive(),   false );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ThreeArraysESFC, LO, GO, Scalar, Node )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.randomize();

    // create the identity matrix, via three arrays constructor
    ArrayRCP<size_t> rowptr(numLocal+1);
    ArrayRCP<LO>     colind(numLocal); // one unknown per row
    ArrayRCP<Scalar> values(numLocal); // one unknown per row

    for(size_t i=0; i<numLocal; i++){
      rowptr[i] = i;
      colind[i] = Teuchos::as<LO>(i);
      values[i] = ScalarTraits<Scalar>::one();
    }
    rowptr[numLocal]=numLocal;

    RCP<CrsMatrix<Scalar,LO,GO,Node> > eye = rcp(new MAT(map,map,rowptr,colind,values));
    TEST_NOTHROW( eye->expertStaticFillComplete(map,map) );

    // test the properties
    TEST_EQUALITY(eye->getGlobalNumEntries()  , numImages*numLocal);
    TEST_EQUALITY(eye->getLocalNumEntries()      , numLocal);
    TEST_EQUALITY(eye->getGlobalNumRows()      , numImages*numLocal);
    TEST_EQUALITY(eye->getLocalNumRows()          , numLocal);
    TEST_EQUALITY(eye->getLocalNumCols()          , numLocal);
    TEST_EQUALITY( Tpetra::Details::getGlobalNumDiags (*eye), static_cast<GO> (numImages*numLocal) );
    TEST_EQUALITY( Tpetra::Details::getLocalNumDiags (*eye), static_cast<LO> (numLocal) );
    TEST_EQUALITY(eye->getGlobalMaxNumRowEntries(), 1);
    TEST_EQUALITY(eye->getLocalMaxNumRowEntries()    , 1);
    TEST_EQUALITY(eye->getIndexBase()          , 0);
    TEST_EQUALITY_CONST(eye->getRowMap()!=Teuchos::null, true);
    TEST_EQUALITY_CONST(eye->getColMap()!=Teuchos::null, true);
    TEST_EQUALITY_CONST(eye->getDomainMap()!=Teuchos::null, true);
    TEST_EQUALITY_CONST(eye->getRangeMap()!=Teuchos::null, true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getColMap())   , true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getDomainMap()), true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getRangeMap()) , true);

    // test the action
    mvres.randomize();
    eye->apply(mvrand,mvres);
    mvres.update(-ST::one(),mvrand,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm1(norms());
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms,zeros);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, SetAllValues, LO, GO, Scalar, Node )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t numImages = comm->getSize();

    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.randomize();

    // create the identity matrix, via three arrays constructor
    ArrayRCP<size_t> rowptr(numLocal+1);
    ArrayRCP<LO>     colind(numLocal); // one unknown per row
    ArrayRCP<Scalar> values(numLocal); // one unknown per row

    for(size_t i=0; i<numLocal; i++){
      rowptr[i] = i;
      colind[i] = Teuchos::as<LO>(i);
      values[i] = ScalarTraits<Scalar>::one();
    }
    rowptr[numLocal]=numLocal;

    RCP<CrsMatrix<Scalar,LO,GO,Node> > eye = rcp(new MAT(map,map,0));
    TEST_NOTHROW( eye->setAllValues(rowptr,colind,values) );
    TEST_NOTHROW( eye->expertStaticFillComplete(map,map) );

    // test that it's possible to call setAllValues on a fillComplete'd matrix
    TEST_NOTHROW( eye->resumeFill() );
    TEST_NOTHROW( eye->setAllValues(rowptr,colind,values) );
    TEST_NOTHROW( eye->expertStaticFillComplete(map,map) );

    // test that it's possible to call setAllValues using a KokkosSparse::CrsMatrix
    RCP<CrsMatrix<Scalar,LO,GO,Node> > eye2 = rcp(new MAT(map,map,0));
    TEST_NOTHROW( eye2->fillComplete(map,map) );

    TEST_EQUALITY(eye2->getGlobalNumEntries()  , 0);
    TEST_EQUALITY(eye2->getLocalNumEntries()      , 0);
    auto pupil = eye->getLocalMatrixDevice();
    TEST_NOTHROW( eye2->resumeFill() );
    TEST_NOTHROW( eye2->setAllValues(pupil) );
    TEST_NOTHROW( eye2->expertStaticFillComplete(map,map) );
    TEST_EQUALITY(eye2->getGlobalNumEntries()  , numImages*numLocal);
    TEST_EQUALITY(eye2->getLocalNumEntries()      , numLocal);
    
    

    // test the properties
    TEST_EQUALITY(eye->getGlobalNumEntries()  , numImages*numLocal);
    TEST_EQUALITY(eye->getLocalNumEntries()      , numLocal);
    TEST_EQUALITY(eye->getGlobalNumRows()      , numImages*numLocal);
    TEST_EQUALITY(eye->getLocalNumRows()          , numLocal);
    TEST_EQUALITY(eye->getLocalNumCols()          , numLocal);
    TEST_EQUALITY( Tpetra::Details::getGlobalNumDiags (*eye), static_cast<GO> (numImages*numLocal) );
    TEST_EQUALITY( Tpetra::Details::getLocalNumDiags (*eye), static_cast<LO> (numLocal) );
    TEST_EQUALITY(eye->getGlobalMaxNumRowEntries(), 1);
    TEST_EQUALITY(eye->getLocalMaxNumRowEntries()    , 1);
    TEST_EQUALITY(eye->getIndexBase()          , 0);
    TEST_EQUALITY_CONST(eye->getRowMap()!=Teuchos::null, true);
    TEST_EQUALITY_CONST(eye->getColMap()!=Teuchos::null, true);
    TEST_EQUALITY_CONST(eye->getDomainMap()!=Teuchos::null, true);
    TEST_EQUALITY_CONST(eye->getRangeMap()!=Teuchos::null, true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getColMap())   , true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getDomainMap()), true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getRangeMap()) , true);
    // test the action
    mvres.randomize();
    eye->apply(mvrand,mvres);
    mvres.update(-ST::one(),mvrand,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm1(norms());
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms,zeros);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, RemoveZeros, LO, GO, Scalar, Node )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t numImages = comm->getSize();

    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.randomize();

    // create the identity matrix, via three arrays constructor
    ArrayRCP<size_t> rowptr(numLocal+1);
    ArrayRCP<LO>     colind(numLocal); // one unknown per row
    ArrayRCP<Scalar> values(numLocal); // one unknown per row

    for(size_t i=0; i<numLocal; i++){
      rowptr[i] = i;
      colind[i] = Teuchos::as<LO>(i);
      values[i] = ScalarTraits<Scalar>::one();
    }
    values[numLocal-1] = MT::eps();
    rowptr[numLocal]=numLocal;

    RCP<CrsMatrix<Scalar,LO,GO,Node> > eye = rcp(new MAT(map,map,0));
    TEST_NOTHROW( eye->setAllValues(rowptr,colind,values) );
    TEST_NOTHROW( eye->expertStaticFillComplete(map,map) );

    TEST_EQUALITY(eye->getGlobalNumEntries()  , numImages*numLocal);

    TEST_NOTHROW(removeCrsMatrixZeros(*eye));
    TEST_EQUALITY(eye->getGlobalNumEntries()  , numImages*numLocal);
    Mag threshold = 10*MT::eps();
    TEST_NOTHROW(removeCrsMatrixZeros(*eye,threshold));
    TEST_EQUALITY(eye->getGlobalNumEntries()  , numImages*(numLocal-1));

  }


  // mfh 17 Mar 2014: In CrsMatrix, in both the classic and Kokkos
  // refactor versions, the undocumented parameter "Preserve Local
  // Graph" now defaults to true.  This makes the following scenario
  // work by default:
  //
  //   1. Create a CrsMatrix A that creates and owns its graph (i.e., don't
  //      use the constructor that takes an RCP<const Tpetra::CrsGraph> or
  //      a local graph)
  //   2. Set an entry in the matrix A, and call fillComplete on it
  //   3. Create a CrsMatrix B using A's graph (obtained via
  //      A.getCrsGraph()), so that B has a const (a.k.a. "static") graph
  //   4. Change a value in B (you can't change its structure), and call
  //      fillComplete on B
  //
  // Before this commit, the above scenario didn't work by default.
  // This is because A's first fillComplete call would call
  // fillLocalGraphAndMatrix, which by default sets the local graph to
  // null.  As a result, from that point,
  // A.getCrsGraph()->getLocalGraph() returns null, which makes B's
  // fillComplete throw an exception.  The only way to make this
  // scenario work was to set A's "Preserve Local Graph" parameter to
  // true.  (It defaulted to false.)
  //
  // The idea behind this nonintuitive behavior was for the local
  // sparse ops object to own all the data.  This might make sense if
  // it is a third-party library that takes CSR's three arrays and
  // copies them into its own storage format.  In that case, it might
  // be a good idea to free the original three CSR arrays, in order to
  // avoid duplicate storage.  However, resumeFill never had a way to
  // get that data back out of the local sparse ops object.  Rather
  // than try to implement that, it's easier just to make "Preserve
  // Local Graph" default to true.
  //
  // The possible data duplication mentioned in the previous paragraph
  // can never happen with the Kokkos refactor version of CrsMatrix,
  // since it insists on controlling the matrix representation itself.
  // This makes the code shorter and easier to read, and also ensures
  // efficient fill.  That will in turn make the option unnecessary.
  //
  // Thanks to Andrey for pointing this out and giving a test case
  // (given below, with only minor modifications from his patch).
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, GraphOwnedByFirstMatrixSharedBySecond, LO, GO, Scalar, Node )
  {
    typedef Tpetra::global_size_t GST;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    auto comm = Tpetra::getDefaultComm ();
    auto map = createContigMapWithNode<LO, GO, Node> (INVALID, 1, comm);

    // construct matrix
    CrsMatrix<Scalar,LO,GO,Node> A (map, map, 1);
    A.insertLocalValues (0, tuple<LO> (0), tuple<Scalar> (STS::zero ()));
    A.fillComplete (map, map);

    // construct second matrix using the first matrix's graph
    CrsMatrix<Scalar,LO,GO,Node> B (A.getCrsGraph ());
    B.resumeFill ();
    B.replaceLocalValues (0, tuple<LO> (0), tuple<Scalar> (STS::one ()));
    B.fillComplete (map, map);
  }

 ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ExtractBlockDiagonal, LO, GO, Scalar, Node )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    Scalar SC_one = ST::one();

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t numImages = comm->getSize();

    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.randomize();

    // create the identity matrix, via three arrays constructor
    ArrayRCP<size_t> rowptr(numLocal+1);
    ArrayRCP<LO>     colind(numLocal); // one unknown per row
    ArrayRCP<Scalar> values(numLocal); // one unknown per row

    for(size_t i=0; i<numLocal; i++){
      rowptr[i] = i;
      colind[i] = Teuchos::as<LO>(i);
      values[i] = ScalarTraits<Scalar>::one();
    }
    rowptr[numLocal]=numLocal;

    RCP<CrsMatrix<Scalar,LO,GO,Node> > eye = rcp(new MAT(map,map,0));
    TEST_NOTHROW( eye->setAllValues(rowptr,colind,values) );
    TEST_NOTHROW( eye->expertStaticFillComplete(map,map) );

    // Now, rip out some diagonals.
    RCP<MV> diag1 = rcp(new MV(map,1));  diag1->putScalar(SC_one);
    RCP<MV> diag2 = rcp(new MV(map,2));  diag2->putScalar(SC_one);
    RCP<MV> diag5 = rcp(new MV(map,5));  diag5->putScalar(SC_one);

    Tpetra::Details::extractBlockDiagonal(*eye,*diag1);
    Tpetra::Details::extractBlockDiagonal(*eye,*diag2);
    Tpetra::Details::extractBlockDiagonal(*eye,*diag5);

    // Check norms
    Array<Mag> norms1(1), norms2(2), norms5(5);
    diag1->norm1(norms1());
    diag2->norm1(norms2());
    diag5->norm1(norms5());

    Array<Mag> cmp1(1), cmp2(2), cmp5(5);
    cmp1[0] = numLocal*numImages;
    cmp2[0] = numLocal*numImages/2;
    cmp2[1] = numLocal*numImages/2;

    cmp5[0] = numLocal*numImages/5;
    cmp5[1] = numLocal*numImages/5;
    cmp5[2] = numLocal*numImages/5;
    cmp5[3] = numLocal*numImages/5;
    cmp5[4] = numLocal*numImages/5;

    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms1,cmp1);
      TEST_COMPARE_ARRAYS(norms2,cmp2);
      TEST_COMPARE_ARRAYS(norms5,cmp5);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms1,cmp1,MT::zero());
      TEST_COMPARE_FLOATING_ARRAYS(norms2,cmp2,MT::zero());
      TEST_COMPARE_FLOATING_ARRAYS(norms5,cmp5,MT::zero());
    }
  }

 ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ScaleBlockDiagonal, LO, GO, Scalar, Node )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    Scalar SC_one = ST::one();

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t numImages = comm->getSize();

    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.randomize();

    // create the identity matrix, via three arrays constructor
    ArrayRCP<size_t> rowptr(numLocal+1);
    ArrayRCP<LO>     colind(numLocal); // one unknown per row
    ArrayRCP<Scalar> values(numLocal); // one unknown per row

    for(size_t i=0; i<numLocal; i++){
      rowptr[i] = i;
      colind[i] = Teuchos::as<LO>(i);
      values[i] = SC_one + SC_one;
    }
    rowptr[numLocal]=numLocal;

    RCP<CrsMatrix<Scalar,LO,GO,Node> > eye2 = rcp(new MAT(map,map,0));
    TEST_NOTHROW( eye2->setAllValues(rowptr,colind,values) );
    TEST_NOTHROW( eye2->expertStaticFillComplete(map,map) );

    // Now, rip out some diagonals.
    RCP<MV> diag1 = rcp(new MV(map,1));  diag1->putScalar(SC_one);
    RCP<MV> diag2 = rcp(new MV(map,2));  diag2->putScalar(SC_one);
    RCP<MV> diag5 = rcp(new MV(map,5));  diag5->putScalar(SC_one);

    Tpetra::Details::extractBlockDiagonal(*eye2,*diag1);
    Tpetra::Details::extractBlockDiagonal(*eye2,*diag2);
    Tpetra::Details::extractBlockDiagonal(*eye2,*diag5);

    // Now, let's rescale some vectors
    RCP<MV> toScale1 = rcp(new MV(map,2)); toScale1->putScalar(SC_one);
    RCP<MV> toScale2 = rcp(new MV(map,2)); toScale2->putScalar(SC_one);
    RCP<MV> toScale5 = rcp(new MV(map,2)); toScale5->putScalar(SC_one);

    Tpetra::Details::inverseScaleBlockDiagonal(*diag1,true,*toScale1);
    Tpetra::Details::inverseScaleBlockDiagonal(*diag2,false,*toScale2);
    Tpetra::Details::inverseScaleBlockDiagonal(*diag5,false,*toScale5);

    // Check norms
    Array<Mag> norms1(2), norms2(2), norms5(2);
    toScale1->norm1(norms1());
    toScale2->norm1(norms2());
    toScale5->norm1(norms5());

    Array<Mag> cmp1(2), cmp2(2), cmp5(2);
    cmp1[0] = numLocal*numImages / 2;
    cmp1[1] = numLocal*numImages / 2;

    cmp2[0] = numLocal*numImages / 2;
    cmp2[1] = numLocal*numImages / 2;

    cmp5[0] = numLocal*numImages / 2;
    cmp5[1] = numLocal*numImages / 2;

    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms1,cmp1);
      TEST_COMPARE_ARRAYS(norms2,cmp2);
      TEST_COMPARE_ARRAYS(norms5,cmp5);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms1,cmp1,MT::zero());
      TEST_COMPARE_FLOATING_ARRAYS(norms2,cmp2,MT::zero());
      TEST_COMPARE_FLOATING_ARRAYS(norms5,cmp5,MT::zero());
    }

  }

 ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ScaleBlockDiagonal_Forward, LO, GO, Scalar, Node )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    Scalar SC_one = ST::one();

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t numImages = comm->getSize();

    const size_t numLocal = 8;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);

    // create block lower-triangular, via three arrays constructor
    ArrayRCP<size_t> rowptr(numLocal+1);
    ArrayRCP<LO>     colind((int)(1.5*numLocal)); // 1.5 unknowns per row
    ArrayRCP<Scalar> values((int)(1.5*numLocal)); // 1.5 two unknowns per row

    int nnz=0;
    for(size_t i=0; i<numLocal; i++){
      rowptr[i] = nnz;
      if(i%2 == 0) {
        // Diagonal entry only (2)
        colind[nnz] = Teuchos::as<LO>(i);
        values[nnz] = SC_one + SC_one;
        nnz++;
      }
      else {
        // Sub-diagonal entry (-1) and diagonal entry (2)
        colind[nnz] = Teuchos::as<LO>(i-1);
        values[nnz] = - SC_one;
        nnz++;
        colind[nnz] = Teuchos::as<LO>(i);
        values[nnz] = SC_one + SC_one;
        nnz++;
      }

    }
    rowptr[numLocal]=nnz;

    RCP<CrsMatrix<Scalar,LO,GO,Node> > blockTri = rcp(new MAT(map,map,0));
    TEST_NOTHROW( blockTri->setAllValues(rowptr,colind,values) );
    TEST_NOTHROW( blockTri->expertStaticFillComplete(map,map) );

    // Now, rip out some diagonals.
    RCP<MV> diag2 = rcp(new MV(map,2));  diag2->putScalar(SC_one);
    RCP<MV> diag4 = rcp(new MV(map,4));  diag4->putScalar(SC_one);

    Tpetra::Details::extractBlockDiagonal(*blockTri,*diag2);
    Tpetra::Details::extractBlockDiagonal(*blockTri,*diag4);

    // Make some vectors
    RCP<MV> toScale2 = rcp(new MV(map,1));
    RCP<MV> toScale4 = rcp(new MV(map,1)); 
    {
      auto v2 = toScale2->getDataNonConst(0);
      auto v4 = toScale4->getDataNonConst(0);
      for(size_t i=0; i<numLocal; i++){
        if(i%2 == 0) {
          v2[i] = SC_one;
          v4[i] = SC_one;
        }
        else {        
          v2[i] = SC_one+SC_one;
          v4[i] = SC_one+SC_one;
        }
      }
    }    
    
    // Scale some vectors
    Tpetra::Details::inverseScaleBlockDiagonal(*diag2,false,*toScale2);
    Tpetra::Details::inverseScaleBlockDiagonal(*diag4,false,*toScale4);

    Kokkos::fence();

    // Check norms
    Array<Mag> norms2(1), norms4(1);
    toScale2->norm1(norms2());
    toScale4->norm1(norms4());

    Array<Mag> cmp2(1), cmp4(1);
    cmp2[0] = 0.875*numLocal*numImages;
    cmp4[0] = 0.875*numLocal*numImages;

    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms2,cmp2);
      TEST_COMPARE_ARRAYS(norms4,cmp4);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms2,cmp2,100*MT::eps());
      TEST_COMPARE_FLOATING_ARRAYS(norms4,cmp4,100*MT::eps());
    }

  }


 ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ScaleBlockDiagonal_Transpose, LO, GO, Scalar, Node )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    Scalar SC_one = ST::one();

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const size_t numImages = comm->getSize();

    const size_t numLocal = 8;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);

    // create block lower-triangular, via three arrays constructor
    ArrayRCP<size_t> rowptr(numLocal+1);
    ArrayRCP<LO>     colind((int)(1.5*numLocal)); // 1.5 unknowns per row
    ArrayRCP<Scalar> values((int)(1.5*numLocal)); // 1.5 two unknowns per row

    int nnz=0;
    for(size_t i=0; i<numLocal; i++){
      rowptr[i] = nnz;
      if(i%2 == 0) {
        // Diagonal entry only (2)
        colind[nnz] = Teuchos::as<LO>(i);
        values[nnz] = SC_one + SC_one;
        nnz++;
      }
      else {
        // Sub-diagonal entry (-1) and diagonal entry (2)
        colind[nnz] = Teuchos::as<LO>(i-1);
        values[nnz] = - SC_one;
        nnz++;
        colind[nnz] = Teuchos::as<LO>(i);
        values[nnz] = SC_one + SC_one;
        nnz++;
      }

    }
    rowptr[numLocal]=nnz;

    RCP<CrsMatrix<Scalar,LO,GO,Node> > blockTri = rcp(new MAT(map,map,0));
    TEST_NOTHROW( blockTri->setAllValues(rowptr,colind,values) );
    TEST_NOTHROW( blockTri->expertStaticFillComplete(map,map) );

    // Now, rip out some diagonals.
    RCP<MV> diag2 = rcp(new MV(map,2));  diag2->putScalar(SC_one);
    RCP<MV> diag4 = rcp(new MV(map,4));  diag4->putScalar(SC_one);

    Tpetra::Details::extractBlockDiagonal(*blockTri,*diag2);
    Tpetra::Details::extractBlockDiagonal(*blockTri,*diag4);

    // Make some vectors
    RCP<MV> toScale2 = rcp(new MV(map,1));
    RCP<MV> toScale4 = rcp(new MV(map,1)); 
    {
      auto v2 = toScale2->getDataNonConst(0);
      auto v4 = toScale4->getDataNonConst(0);
      for(size_t i=0; i<numLocal; i++){
        if(i%2 == 0) {
          v2[i] = SC_one;
          v4[i] = SC_one;
        }
        else {        
          v2[i] = SC_one+SC_one;
          v4[i] = SC_one+SC_one;
        }
      }    
    }

    // Now, let's rescale some vectors
    Tpetra::Details::inverseScaleBlockDiagonal(*diag2,true,*toScale2);
    Tpetra::Details::inverseScaleBlockDiagonal(*diag4,true,*toScale4);

    Kokkos::fence();

    // Check norms
    Array<Mag> norms2(1), norms4(1);
    toScale2->norm1(norms2());
    toScale4->norm1(norms4());

    Array<Mag> cmp2(1), cmp4(1);
    cmp2[0] = numLocal*numImages;
    cmp4[0] = numLocal*numImages;

    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms2,cmp2);
      TEST_COMPARE_ARRAYS(norms4,cmp4);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms2,cmp2,100*MT::eps());
      TEST_COMPARE_FLOATING_ARRAYS(norms4,cmp4,100*MT::eps());
    }

  }


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ApplyHelpers, LO, GO, Scalar, Node )
  {
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    //    typedef  Operator<Scalar,LO,GO,Node> OP;
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef Map<LO,GO,Node> MAP;
    typedef typename ST::magnitudeType Mag;
    using values_type = typename MAT::local_matrix_device_type::values_type;
    using range_policy = Kokkos::RangePolicy<typename Node::device_type::execution_space>;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    int numRanks = comm->getSize();

    Scalar ZERO = ST::zero(), ONE = ST::one();
    Mag MT_ZERO = Teuchos::ScalarTraits<Mag>::zero();
    int nsize=100;

    /* Create the identity matrix, three rows per proc */
    RCP<MAT> A1 = Tpetra::Utils::MatrixGenerator<MAT>::generate_miniFE_matrix(nsize, comm);
    if(!A1->isFillComplete()) A1->fillComplete();
    RCP<const MAP> map = A1->getRowMap();

    /* Now take the same thing, but multiply the entries by 2*/
    auto lclMtx1 = A1->getLocalMatrixDevice();
    values_type A1_values = lclMtx1.values;
    values_type A2_values("values",A1_values.extent(0));
    values_type A2_abs_values("absolute values",A1_values.extent(0));
    Kokkos::parallel_for( range_policy(0,A1_values.extent(0)),KOKKOS_LAMBDA(const int i){
        A2_values[i] = A1_values[i] + A1_values[i];
        A2_abs_values[i] = Kokkos::abs(A2_values[i]);
      });

    RCP<MAT> A2 = rcp(new MAT(map,A1->getColMap(),lclMtx1.graph.row_map,lclMtx1.graph.entries,A2_values));
    A2->expertStaticFillComplete(A1->getDomainMap(),A1->getRangeMap(),A1->getGraph()->getImporter(),A1->getGraph()->getExporter());
    // Entrywise absolute value of A2, used for choosing a tolerance
    RCP<MAT> A2_abs = rcp(new MAT(map,A1->getColMap(),lclMtx1.graph.row_map,lclMtx1.graph.entries,A2_abs_values));
    A2_abs->expertStaticFillComplete(A1->getDomainMap(),A1->getRangeMap(),A1->getGraph()->getImporter(),A1->getGraph()->getExporter());

    /* Allocate multivectors */
    MV X(map,2), Y1a(map,2), Y1b(map,2), Y2a(map,2), Y2b(map,2), compare(map,2);
    MV Y_abs(map,2);
    Array<Mag> norm(2), exact(2,MT_ZERO);
    X.putScalar(ONE);

    A2_abs->apply(X, Y_abs);
    // both columns of X are filled with 1, so just use the first column here
    Mag absNorm = Y_abs.getVector(0)->norm1();

    // Do a std::vector version 
    {
      Y1a.putScalar(ZERO);Y1b.putScalar(ZERO);
      Y2a.putScalar(ZERO);Y2b.putScalar(ZERO);

      std::vector<MAT*> matrices = {A1.get(),A2.get()};
      std::vector<MV*>  outvec   = {&Y1b,&Y2b};

      // Comparison guy
      A1->apply(X,Y1a);
      A2->apply(X,Y2a);

      // Apply Helpers
      Tpetra::batchedApply(matrices,X,outvec);
     
      // Compare
      compare.update(ONE,Y1a,-ONE,Y1b,ZERO);
      compare.norm1(norm());

      if (ST::isOrdinal) {TEST_COMPARE_ARRAYS(exact,norm);}
      else{TEST_ABSOLUTE_COMPARE_FLOATING_ARRAYS(exact,norm,absNorm*testingTol<Mag>());}

      compare.update(ONE,Y2a,-ONE,Y2b,ZERO);
      compare.norm1(norm());
      if (ST::isOrdinal) {TEST_COMPARE_ARRAYS(exact,norm);}
      else {TEST_ABSOLUTE_COMPARE_FLOATING_ARRAYS(exact,norm,absNorm*testingTol<Mag>());}
    }

    // Do a std::vector version w/ a nice combination of options
    std::vector<Scalar> alpha = {ONE+ONE, ONE, ZERO};
    std::vector<Scalar> beta  = {ONE+ONE, ONE, ZERO};
    for(int i=0; i<(int)alpha.size(); i++) {
      for(int j=0; j<(int)beta.size(); j++) {
        Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList());
        Y1a.putScalar(ONE);Y1b.putScalar(ONE);
        Y2a.putScalar(ZERO);Y2b.putScalar(ZERO);
        
        std::vector<MAT*> matrices = {A1.get(),A2.get()};
        std::vector<MV*>  outvec   = {&Y1b,&Y2b};
        
        // Comparison guy
        A1->apply(X,Y1a,Teuchos::NO_TRANS,alpha[i],beta[j]);
        A2->apply(X,Y2a,Teuchos::NO_TRANS,alpha[i],beta[j]);
      
        // Apply Helpers
        Tpetra::batchedApply(matrices,X,outvec,alpha[i],beta[j],params);
      
        // Check to make sure this guy decided we can batch appropriately
        bool did_we_batch = params->get<bool>("can batch");
        bool should_we_batch = (numRanks>1)?true:false;
        TEST_EQUALITY(did_we_batch,should_we_batch);
        
        // Compare
        compare.update(ONE,Y1a,-ONE,Y1b,ZERO);
        compare.norm1(norm());
        if (ST::isOrdinal) {TEST_COMPARE_ARRAYS(exact,norm);}
        else{TEST_ABSOLUTE_COMPARE_FLOATING_ARRAYS(exact,norm,absNorm*testingTol<Mag>());}
        
        compare.update(ONE,Y2a,-ONE,Y2b,ZERO);
        compare.norm1(norm());
        if (ST::isOrdinal) {TEST_COMPARE_ARRAYS(exact,norm);}
        else {TEST_ABSOLUTE_COMPARE_FLOATING_ARRAYS(exact,norm,absNorm*testingTol<Mag>());}
      }
    }

    // Teuchos::Array version
    // Do a std::vector version 
    {
      Y1a.putScalar(ZERO);Y2a.putScalar(ZERO);
      Y1b.putScalar(ZERO);Y2b.putScalar(ZERO);

      Teuchos::Array<MAT*> matrices(2); matrices[0]=A1.get(); matrices[1]=A2.get();
      Teuchos::Array<MV*>  outvec(2);   outvec[0]=&Y1b; outvec[1]=&Y2b;

      // Comparison guy
      A1->apply(X,Y1a);
      A2->apply(X,Y2a);

      // Apply Helpers
      Tpetra::batchedApply(matrices,X,outvec);

      // Compare
      compare.update(ONE,Y1a,-ONE,Y1b,ZERO);
      compare.norm1(norm());
      if (ST::isOrdinal) {TEST_COMPARE_ARRAYS(exact,norm);}
      else{TEST_ABSOLUTE_COMPARE_FLOATING_ARRAYS(exact,norm,absNorm*testingTol<Mag>());}

      compare.update(ONE,Y2a,-ONE,Y2b,ZERO);
      compare.norm1(norm());
      if (ST::isOrdinal) {TEST_COMPARE_ARRAYS(exact,norm);}
      else {TEST_ABSOLUTE_COMPARE_FLOATING_ARRAYS(exact,norm,absNorm*testingTol<Mag>());}
    }

  }


//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, AlphaBetaMultiply, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ActiveFill,        LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, Typedefs,          LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ThreeArraysESFC,   LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, SetAllValues,      LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, GraphOwnedByFirstMatrixSharedBySecond, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ExtractBlockDiagonal,      LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ApplyHelpers,      LO, GO, SCALAR, NODE )  

#define UNIT_TEST_GROUP_NO_ORDINAL_SCALAR( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ScaleBlockDiagonal,      LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ScaleBlockDiagonal_Forward,     LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ScaleBlockDiagonal_Transpose,     LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, RemoveZeros,       LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( UNIT_TEST_GROUP_NO_ORDINAL_SCALAR )

}
