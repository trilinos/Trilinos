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
#include <TpetraExt_MatrixMatrix.hpp>

namespace {
  using Tpetra::TestingUtilities::getNode;
  using Tpetra::TestingUtilities::getDefaultComm;

  using std::endl;
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
  using Tpetra::createVector;
  using Tpetra::createCrsMatrix;
  using Tpetra::DefaultPlatform;
  using Tpetra::ProfileType;
  using Tpetra::StaticProfile;
  using Tpetra::DynamicProfile;
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
    typedef typename ArrayView<const LO>::size_type size_type;

    MT mySum = STM::zero ();
    Array<LO> inds (matrix->getNodeMaxNumRowEntries ());
    Array<Scalar> vals (matrix->getNodeMaxNumRowEntries ());

    const size_t myNumRows = matrix->getNodeNumRows ();
    for (size_t i = 0; i < myNumRows; ++i) {
      const LO myRow = as<LO> (i);
      const size_t numRowEnts = matrix->getNumEntriesInLocalRow (myRow);
      ArrayView<const LO> indsView = inds.view (0, as<size_type> (numRowEnts));
      ArrayView<const ST> valsView = vals.view (0, as<size_type> (numRowEnts));
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
    RCP<Node> node = getNode<Node>();
    typedef Vector<Scalar,LO,GO,Node> VEC;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    int numProcs = comm->getSize();
    int myRank = comm->getRank();

    global_size_t numGlobal = 4*numProcs;
    RCP<const Map<LO,GO,Node> > map = createUniformContigMapWithNode<LO,GO>(numGlobal,comm,node);
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

    RCP<MAT> diffMat1 = createCrsMatrix<Scalar, LO, GO, Node>(map,3);
    RCP<MAT> diffMat2 = createCrsMatrix<Scalar, LO, GO, Node>(map,3);
    Scalar sOne = ScalarTraits<Scalar>::one();
    Tpetra::MatrixMatrix::Add(*matrix, false, -sOne, *answerMatrix, false, sOne, diffMat1);
    diffMat1->fillComplete();
    Tpetra::MatrixMatrix::Add(*matrix2, false, -sOne, *answerMatrix, false, sOne, diffMat2);
    diffMat2->fillComplete();
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
