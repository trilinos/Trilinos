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

// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_ConfigDefs.hpp>

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
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, TheEyeOfTruthDistAlloc, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
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
    TEST_EQUALITY(eye->getNodeNumEntries()      , numLocal);
    TEST_EQUALITY(eye->getGlobalNumRows()      , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumRows()          , numLocal);
    TEST_EQUALITY(eye->getNodeNumCols()          , numLocal);
    TEST_EQUALITY(eye->getGlobalNumDiags() , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumDiags()     , numLocal);
    TEST_EQUALITY(eye->getGlobalMaxNumRowEntries(), 1);
    TEST_EQUALITY(eye->getNodeMaxNumRowEntries()    , 1);
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


  // mfh 08 Mar 2013: This test wasn't being instantiated, so I
  // disabled it to save compilation time.
#if 0
  ////
  TEUCHOS_UNIT_TEST( CrsMatrix, Convert )
  {
    typedef KokkosClassic::DefaultNode::DefaultNodeType Node;
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<double> ST;
    typedef OrdinalTraits<int> LOT;
    typedef CrsMatrix<double,int,int,Node> DMat;
    typedef CrsMatrix<   int,int,int,Node> IMat;
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numProcs = comm->getSize();

    global_size_t numGlobal = 2*numProcs;
    RCP<const Map<int,int,Node> > map = createUniformContigMapWithNode<int,int>(numGlobal,comm,node);

    RCP<DMat> dmatrix = createCrsMatrix<double>(map, 3);

    const int maxglobalrow = map->getMaxAllGlobalIndex();
    for(int r = map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r)
    {
      if (r==0) {
        dmatrix->insertGlobalValues( r, tuple<int>(0,1), tuple<double>(1,1) );
      }
      else if( r == maxglobalrow ) {
        dmatrix->insertGlobalValues( r, tuple<int>(maxglobalrow-1,maxglobalrow), tuple<double>(1,1) );
      }
      else{
        dmatrix->insertGlobalValues( r, tuple<int>(r-1,r,r+1), tuple<double>(1,1,1) );
      }
    }

    RCP<ParameterList> params = parameterList();
    params->set("Preserve Local Graph",true);
    dmatrix->fillComplete(params);

    RCP<IMat> imatrix = dmatrix->convert<int>();

    // check graphs
    TEST_EQUALITY( dmatrix->getGraph(),    imatrix->getGraph() );
    TEST_EQUALITY( dmatrix->getCrsGraph(), imatrix->getCrsGraph() );
    // check entries
    for (int i=map->getMinLocalIndex(); i <= map->getMaxLocalIndex(); ++i) {
      ArrayView<const int> indsd, indsi, valsi;
      ArrayView<const double> valsd;
      dmatrix->getLocalRowView(i, indsd, valsd);
      imatrix->getLocalRowView(i, indsi, valsi);
      TEST_COMPARE_ARRAYS(indsd,indsi);
      TEST_COMPARE_ARRAYS(valsd,valsi);
    }
  }
#endif // 0

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, AlphaBetaMultiply, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef  Operator<Scalar,LO,GO,Node> OP;
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    const size_t THREE = 3;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t myImageID = comm->getRank();
    // create a Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,THREE,comm,node);

    /* Create the identity matrix, three rows per proc */
    RCP<OP> AOp;
    {
      RCP<MAT> A = rcp(new MAT(map,1));
      A->insertGlobalValues(3*myImageID,  tuple<GO>(3*myImageID  ), tuple<Scalar>(ST::one()) );
      A->insertGlobalValues(3*myImageID+1,tuple<GO>(3*myImageID+1), tuple<Scalar>(ST::one()) );
      A->insertGlobalValues(3*myImageID+2,tuple<GO>(3*myImageID+2), tuple<Scalar>(ST::one()) );
      A->fillComplete();
      AOp = A;
    }
    MV X(map,1), Y(map,1), Z(map,1);
    const Scalar alpha = ST::random(),
                  beta = ST::random();
    X.randomize();
    Y.randomize();
    // Z = alpha*X + beta*Y
    Z.update(alpha,X,beta,Y,ST::zero());
    // test the action: Y = alpha*I*X + beta*Y = alpha*X + beta*Y = Z
    AOp->apply(X,Y,NO_TRANS,alpha,beta);
    //
    Array<Mag> normY(1), normZ(1);
    Z.norm1(normZ());
    Y.norm1(normY());
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(normY,normZ);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(normY,normZ,2.0*testingTol<Mag>());
    }
  }

  // mfh 08 Mar 2013: The MixedMultiplyOp test wasn't being
  // instantiated (at the end of this file) anyway, so I'm commenting
  // it out for now to speed up compilation of this file.
#if 0
  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, MixedMultiplyOp, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef CrsMatrix<int,LO,GO,Node> IntMAT;
    typedef  Operator<Scalar,LO,GO,Node> OP;
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const size_t THREE = 3;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t myImageID = comm->getRank();
    // create a Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,THREE,comm,node);

    /* Create the integer identity matrix, three rows per proc, wrapped in a Op<Scalar>  */
    RCP<OP> AOp;
    {
      RCP<IntMAT> A = rcp(new IntMAT(map,1));
      A->insertGlobalValues(3*myImageID,  tuple<GO>(3*myImageID  ), tuple<int>(1) );
      A->insertGlobalValues(3*myImageID+1,tuple<GO>(3*myImageID+1), tuple<int>(1) );
      A->insertGlobalValues(3*myImageID+2,tuple<GO>(3*myImageID+2), tuple<int>(1) );
      A->fillComplete();
      AOp = createCrsMatrixMultiplyOp<Scalar>(A.getConst());
    }
    MV X(map,1), Y(map,1), Z(map,1);
    X.randomize();
    Y.randomize();
    // Z = X + Y
    Z.update(ST::one(),X,ST::one(),Y,ST::zero());
    // test the action: Y = I*X + Y = X + Y == Z
    AOp->apply(X,Y,NO_TRANS,ST::one(),ST::one());
    // Z -= Y  -> zero
    Z.update(-ST::one(),Y,ST::one());
    Array<Mag> norms(1), zeros(1,MT::zero());
    Z.norm1(norms());
    if (ST::isOrdinal) {
      TEST_COMPARE_ARRAYS(norms,zeros);
    } else {
      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
    }
  }
#endif // 0

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, Typedefs, LO, GO, Scalar, Node )
  {
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename MAT::scalar_type         scalar_type;
    typedef typename MAT::local_ordinal_type  local_ordinal_type;
    typedef typename MAT::global_ordinal_type global_ordinal_type;
    typedef typename MAT::node_type           node_type;
    TEST_EQUALITY_CONST( (is_same< scalar_type         , Scalar >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< local_ordinal_type  , LO     >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< global_ordinal_type , GO     >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< node_type           , Node   >::value) == true, true );
    typedef RowMatrix<Scalar,LO,GO,Node> RMAT;
    typedef typename RMAT::scalar_type         rmat_scalar_type;
    typedef typename RMAT::local_ordinal_type  rmat_local_ordinal_type;
    typedef typename RMAT::global_ordinal_type rmat_global_ordinal_type;
    typedef typename RMAT::node_type           rmat_node_type;
    TEST_EQUALITY_CONST( (is_same< rmat_scalar_type         , Scalar >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< rmat_local_ordinal_type  , LO     >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< rmat_global_ordinal_type , GO     >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< rmat_node_type           , Node   >::value) == true, true );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ActiveFill, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,1,comm,node);
    const Scalar SZERO = ScalarTraits<Scalar>::zero();
    {
      MAT matrix(map,map,0,DynamicProfile);
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

      //TEST_THROW( matrix.replaceLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) ), std::runtime_error );
      numValid = matrix.replaceLocalValues (0, tuple<LO> (0), tuple<Scalar> (0));
      TEST_EQUALITY( numValid, Teuchos::OrdinalTraits<LO>::invalid () );

      //TEST_THROW( matrix.replaceLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) ), std::runtime_error );
      numValid = matrix.replaceLocalValues (0, tuple<LO> (0), tuple<Scalar> (0));
      TEST_EQUALITY( numValid, Teuchos::OrdinalTraits<LO>::invalid () );

      //TEST_THROW( matrix.sumIntoLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) ), std::runtime_error );
      numValid = matrix.sumIntoLocalValues (0, tuple<LO> (0), tuple<Scalar> (0));
      TEST_EQUALITY( numValid, Teuchos::OrdinalTraits<LO>::invalid () );

      TEST_THROW( matrix.setAllToScalar(SZERO), std::runtime_error );
      TEST_THROW( matrix.scale(SZERO),          std::runtime_error );
      TEST_THROW( matrix.globalAssemble(),      std::runtime_error );
      TEST_THROW( matrix.fillComplete(),        std::runtime_error );
    }
    {
      MAT matrix(map,map,0,DynamicProfile);
      TEST_EQUALITY_CONST( matrix.isFillActive(),   true );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), false );
      matrix.insertLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) );
      //
      RCP<ParameterList> params = parameterList();
      params->set("Optimize Storage",false);
      matrix.fillComplete(params);
      TEST_EQUALITY_CONST( matrix.isFillActive(),   false );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
      //
      matrix.resumeFill();
      TEST_EQUALITY_CONST( matrix.isFillActive(),   true );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), false );
      TEST_NOTHROW( matrix.insertLocalValues ( 0, tuple<LO>(0), tuple<Scalar>(0) ) );
      TEST_NOTHROW( matrix.replaceLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) ) );
      TEST_NOTHROW( matrix.sumIntoLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) ) );
      TEST_NOTHROW( matrix.setAllToScalar(SZERO)                                   );
      TEST_NOTHROW( matrix.scale(SZERO)                                            );
      TEST_NOTHROW( matrix.globalAssemble()                                        );
      //
      TEST_NOTHROW( matrix.fillComplete()                        );
      TEST_EQUALITY_CONST( matrix.isFillActive(),   false );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ThreeArraysESFC, LO, GO, Scalar, Node )
  {

    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = comm->getSize();
    // create a Map
    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
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
    TEST_EQUALITY(eye->getNodeNumEntries()      , numLocal);
    TEST_EQUALITY(eye->getGlobalNumRows()      , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumRows()          , numLocal);
    TEST_EQUALITY(eye->getNodeNumCols()          , numLocal);
    TEST_EQUALITY(eye->getGlobalNumDiags() , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumDiags()     , numLocal);
    TEST_EQUALITY(eye->getGlobalMaxNumRowEntries(), 1);
    TEST_EQUALITY(eye->getNodeMaxNumRowEntries()    , 1);
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

    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = comm->getSize();
    // create a Map
    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
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

    // test the properties
    TEST_EQUALITY(eye->getGlobalNumEntries()  , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumEntries()      , numLocal);
    TEST_EQUALITY(eye->getGlobalNumRows()      , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumRows()          , numLocal);
    TEST_EQUALITY(eye->getNodeNumCols()          , numLocal);
    TEST_EQUALITY(eye->getGlobalNumDiags() , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumDiags()     , numLocal);
    TEST_EQUALITY(eye->getGlobalMaxNumRowEntries(), 1);
    TEST_EQUALITY(eye->getNodeMaxNumRowEntries()    , 1);
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

    // get a Kokkos node and a comm
    RCP<Node> node = getNode<Node> ();
    RCP<const Comm<int> > comm = getDefaultComm ();

    // create Map
    RCP<const Map<LO,GO,Node> > map =
      createContigMapWithNode<LO,GO> (INVALID, 1, comm, node);

    // construct matrix
    CrsMatrix<Scalar,LO,GO,Node> A (map, map, 0, DynamicProfile);
    A.insertLocalValues (0, tuple<LO> (0), tuple<Scalar> (STS::zero ()));
    A.fillComplete (map, map);

    // construct second matrix using the first matrix's graph
    CrsMatrix<Scalar,LO,GO,Node> B (A.getCrsGraph ());
    B.resumeFill ();
    B.replaceLocalValues (0, tuple<LO> (0), tuple<Scalar> (STS::one ()));
    B.fillComplete (map, map);
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
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, GraphOwnedByFirstMatrixSharedBySecond, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}


