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
#include "TpetraExt_MatrixMatrix.hpp"

namespace {
  using std::endl;

  using Teuchos::as;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::arcp;
  using Teuchos::outArg;
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
  using Tpetra::createUniformContigMapWithNode;
  using Tpetra::createVector;
  using Tpetra::createCrsMatrix;
  using Tpetra::OptimizeOption;
  using Tpetra::GloballyDistributed;
  using Tpetra::INSERT;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
  }

  // Compute the Frobenius norm of the matrix.
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  getNorm (const RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& matrix)
  {
    typedef LocalOrdinal LO;
    typedef Scalar ST;
    typedef Teuchos::ScalarTraits<ST> STS;
    typedef typename STS::magnitudeType MT;
    typedef Teuchos::ScalarTraits<MT> STM;
    typedef CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> MAT;

    MT mySum = STM::zero ();
    Array<LO> inds (matrix->getLocalMaxNumRowEntries ());
    Array<Scalar> vals (matrix->getLocalMaxNumRowEntries ());

    const size_t myNumRows = matrix->getLocalNumRows ();
    for (size_t i = 0; i < myNumRows; ++i) {
      const LO myRow = as<LO> (i);
      const size_t numRowEnts = matrix->getNumEntriesInLocalRow (myRow);
      typename MAT::local_inds_host_view_type indsView;
      typename MAT::values_host_view_type valsView;
      matrix->getLocalRowView (myRow, indsView, valsView);
      for (size_t j = 0; j < numRowEnts; ++j) {
        const ST curVal = valsView[j];
        mySum += STS::real (curVal) * STS::real (curVal) +
          STS::imag (curVal) * STS::imag (curVal);
      }
    }
    MT totalSum = 0;
    Teuchos::reduceAll (* (matrix->getComm ()), Teuchos::REDUCE_SUM, 1, &mySum, &totalSum);
    return STM::squareroot (totalSum);
  }

  //
  // UNIT TEST(S)
  //

  // Construct two tridiagonal matrices and scale them by a vector,
  // one on the left and the other on the right.  Then compare the
  // result to the known correct matrix.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, LeftRightScale, LO, GO, Scalar, Node )
  {
    typedef Vector<Scalar,LO,GO,Node> VEC;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    int numProcs = comm->getSize();
    int myRank = comm->getRank();

    global_size_t numGlobal = 4*numProcs;
    RCP<const Map<LO,GO,Node> > map = createUniformContigMapWithNode<LO,GO,Node>(numGlobal,comm);
    RCP<VEC> vector = createVector<Scalar, LO, GO, Node>(map);
    vector->putScalar(2);

    RCP<MAT> matrix = createCrsMatrix<Scalar, LO, GO, Node>(map, 3);
    RCP<MAT> matrix2= createCrsMatrix<Scalar, LO, GO, Node>(map, 3);
    RCP<MAT> answerMatrix = createCrsMatrix<Scalar, LO, GO, Node>(map, 3);

    Array<Scalar> vals = tuple<Scalar>(1,2,3);
    Array<Scalar> answerVals = tuple<Scalar>(2,4,6);
    Array<GO> cols(3,0);
    for(
      GO i = Teuchos::as<GO>(myRank)*4;
      i<(Teuchos::as<GO>(myRank)*4)+4;
      ++i)
    {
      if(i==0){
        cols = tuple<GO>(0,1);
        matrix->insertGlobalValues(i, cols(), vals(1,2));
        matrix2->insertGlobalValues(i, cols(), vals(1,2));
        answerMatrix->insertGlobalValues(i, cols(), answerVals(1,2));
      }
      else if(i==(Teuchos::as<GO>(numProcs-1)*4)+3){
        cols = tuple<GO>(numGlobal-2, numGlobal-1);
        matrix->insertGlobalValues(i, cols(), vals(0,2));
        matrix2->insertGlobalValues(i, cols(), vals(0,2));
        answerMatrix->insertGlobalValues(i, cols(), answerVals(0,2));
      }
      else{
        cols = tuple<GO>(i-1,i,i+1);
        matrix->insertGlobalValues(i, cols(), vals());
        matrix2->insertGlobalValues(i, cols(), vals());
        answerMatrix->insertGlobalValues(i, cols(), answerVals());
      }
    }

    matrix->fillComplete();
    matrix2->fillComplete();
    answerMatrix->fillComplete();
    matrix->leftScale(*vector);
    matrix2->rightScale(*vector);

    Scalar sOne = ScalarTraits<Scalar>::one();

    RCP<MAT> diffMat1 = Tpetra::MatrixMatrix::add(sOne, false, *matrix, -sOne, false, *answerMatrix);
    RCP<MAT> diffMat2 = Tpetra::MatrixMatrix::add(sOne, false, *matrix2, -sOne, false, *answerMatrix);

    Scalar epsilon1 = getNorm(diffMat1)/getNorm(answerMatrix);
    Scalar epsilon2 = getNorm(diffMat2)/getNorm(answerMatrix);
    TEST_COMPARE(ScalarTraits<Scalar>::real(epsilon1), <, 1e-10)
    TEST_COMPARE(ScalarTraits<Scalar>::imag(epsilon1), <, 1e-10)
    TEST_COMPARE(ScalarTraits<Scalar>::real(epsilon2), <, 1e-10)
    TEST_COMPARE(ScalarTraits<Scalar>::imag(epsilon2), <, 1e-10)
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, LeftRightScale, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}
