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

#include <Teuchos_UnitTestHarness.hpp>
#include <TpetraCore_ETIHelperMacros.h>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#include <MatrixMarket_Tpetra.hpp>

namespace { // (anonymous)

  // mfh 30 Jul 2014: This is a Tpetra port of Xpetra's CrsMatrix,
  // Tpetra_ReplaceLocalValues test, which lives in
  // xpetra/test/CrsMatrix/CrsMatrix_UnitTests.cpp.
  //
  // The Xpetra version of this test currently fails (as of 30 Jul
  // 2014) with the Kokkos refactor version of Tpetra, and has for a
  // couple weeks at least.  Porting this test to Tpetra should help
  // us figure out if this is a Tpetra or an Xpetra issue.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ReplaceLocalValues, Scalar, LO, GO, Node )
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    using std::cerr;
    using std::endl;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
    typedef Tpetra::Vector<Scalar, LO, GO, Node> vec_type;
    typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
    typedef typename Teuchos::Array<LO>::size_type size_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType MT;
    typedef Teuchos::ScalarTraits<MT> STM;
    const Scalar ZERO = STS::zero ();
    const Scalar ONE = STS::one ();
    const Scalar FIVE = static_cast<Scalar> (5.0);
    const bool extraDebug = false; // for extra-verbose output
    int lclSuccess = 0;
    int gblSuccess = 0;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();

    if (myRank == 0) {
      cerr << "Tpetra replaceLocalValues test" << endl
           << "Create Map and matrix" << endl;
    }

    // generate problem
    LO nEle = 63;
    RCP<const map_type> map = rcp (new map_type (nEle, 0, comm));

      TEUCHOS_TEST_FOR_EXCEPTION(
        ! Kokkos::is_initialized (), std::logic_error,
        "Kokkos is not initialized!" );


    RCP<crs_matrix_type> matrix = rcp (new crs_matrix_type (map, 10));
    const LO NumMyElements = map->getNodeNumElements ();
    Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList ();

    // Make the matrix the identity matrix.
    if (myRank == 0) {
      cerr << "Fill matrix by calling insertGlobalValues" << endl;
    }
    for (LO i = 0; i < NumMyElements; ++i) {
      matrix->insertGlobalValues (MyGlobalElements[i],
                                  Teuchos::tuple<GO> (MyGlobalElements[i]),
                                  Teuchos::tuple<Scalar> (ONE));
    }

    if (myRank == 0) {
      cerr << "Call fillComplete on the matrix" << endl;
    }

    matrix->fillComplete ();
    TEST_ASSERT( matrix->isFillComplete () );
    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (myRank == 0) {
      cerr << "The matrix is now the identity matrix.  Make sure that mat-vec works" << endl;
    }
    comm->barrier ();
    {
      vec_type x (matrix->getDomainMap ());
      x.putScalar (ONE);
      vec_type y (matrix->getDomainMap ());
      y.putScalar (ZERO);
      matrix->apply (x, y);

      if (myRank == 0) {
        cerr << "Test that matrix->apply(x,y) with x all ones fills y "
          "with all ones, by computing the 1-norm of y" << endl;
      }
      comm->barrier ();
      {
        const MT expectedNorm1 = static_cast<MT> (y.getGlobalLength ());
        const MT actualNorm1 = y.norm1 ();
        TEST_EQUALITY( actualNorm1, expectedNorm1 );

        std::ostringstream os;
        os << "  Proc " << myRank << ": ";
        if (actualNorm1 != expectedNorm1) {
          os << "INCORRECT: {expectedNorm1: " << expectedNorm1
             << ", actualNorm1: " << actualNorm1 << "}" << endl;
        } else {
          os << "Correct" << endl;
        }
        cerr << os.str ();
      }
    }

    if (myRank == 0) {
      cerr << "Call resumeFill on the matrix" << endl;
    }
    matrix->resumeFill ();
    comm->barrier ();
    if (myRank == 0) {
      cerr << "Got through resumeFill on the matrix" << endl;
    }
    TEST_ASSERT( ! matrix->isFillComplete () );

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // Change the 0,0 local entry, on each process, to be FIVE.
    if (myRank == 0) {
      cerr << "Modify entries of the matrix using replaceLocalValues, "
        "and test the result before calling fillComplete" << endl;
    }
    Teuchos::Array<LO> indout (1, 0);
    Teuchos::Array<Scalar> valout (1, FIVE);

    // Every process should have a local row index 0.
    TEST_ASSERT( map->isNodeLocalElement (0) );
    TEST_ASSERT( ! matrix->getColMap ().is_null () );

    LO actualLclNumReplaced = 0;
    LO expectedLclNumReplaced = 0;
    if (map->isNodeLocalElement (0) && ! matrix->getColMap ().is_null ()) {
      expectedLclNumReplaced = 1;

      bool validLocalColumnIndices = true;
      for (size_type k = 0; k < indout.size (); ++k) {
        if (! matrix->getColMap ()->isNodeLocalElement (indout[k])) {
          validLocalColumnIndices = false;
          break;
        }
      }
      // Every process should have a local column index 0.
      TEST_ASSERT( validLocalColumnIndices );
      if (validLocalColumnIndices) {
        // Make sure that we are changing the first diagonal entry on
        // this process.  We determine whether a matrix is diagonal
        // using global indices.
        TEST_ASSERT( matrix->getColMap ()->getGlobalElement (indout[0]) ==
                     map->getGlobalElement (0) );
        // Replace the local (0,0) entry with valout[0].  We know from
        // the above test that the local (0,0) entry is the first
        // diagonal entry on the calling process.
        actualLclNumReplaced =
          matrix->replaceLocalValues (0, indout.view (0, indout.size ()),
                                      valout.view (0, valout.size ()));
      }

      // Make sure that replaceLocalValues worked, by getting the
      // values in the local row 0.
      const size_t numEnt = matrix->getNumEntriesInLocalRow (0);
      TEST_EQUALITY_CONST( numEnt, static_cast<size_t> (1) );

      if (numEnt == static_cast<size_t> (1)) {
        Teuchos::Array<LO> ind (numEnt);
        Teuchos::Array<Scalar> val (numEnt);
        size_t numEntOut = 0;
        matrix->getLocalRowCopy (0, ind (), val (), numEntOut);
        TEST_EQUALITY( numEnt, numEntOut );

        if (numEntOut == static_cast<size_t> (1)) {
          TEST_EQUALITY( ind[0], 0 );
          TEST_EQUALITY( val[0], FIVE );
        }
      }
    }

    TEST_EQUALITY( actualLclNumReplaced, expectedLclNumReplaced );

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (myRank == 0) {
      cerr << "Call fillComplete on matrix for the second time" << endl;
    }
    matrix->fillComplete ();
    comm->barrier ();
    if (myRank == 0) {
      cerr << "Got through fillComplete on matrix for the second time" << endl;
    }
    TEST_ASSERT( matrix->isFillComplete () );

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (myRank == 0) {
      cerr << "Test the result of replaceLocalValues after fillComplete" << endl;
    }
    if (map->isNodeLocalElement (0)) {
      // Make sure that replaceLocalValues worked, by getting the
      // values in the local row 0.
      const size_t numEnt = matrix->getNumEntriesInLocalRow (0);
      TEST_EQUALITY_CONST( numEnt, static_cast<size_t> (1) );

      if (numEnt == static_cast<size_t> (1)) {
        Teuchos::Array<LO> ind (numEnt);
        Teuchos::Array<Scalar> val (numEnt);
        size_t numEntOut = 0;
        matrix->getLocalRowCopy (0, ind (), val (), numEntOut);
        TEST_EQUALITY( numEnt, numEntOut );

        if (numEntOut == static_cast<size_t> (1)) {
          TEST_EQUALITY( ind[0], 0 );
          TEST_EQUALITY( val[0], FIVE );
        }
      }
    }

    // Extract and print the diagonal entries of the matrix.
    if (extraDebug) {
      if (myRank == 0) {
        cerr << "Diagonal entries of the matrix:" << endl;
      }
      vec_type D (matrix->getRowMap ());
      D.putScalar (ZERO);
      matrix->getLocalDiagCopy (D);
      writer_type::writeDense (cerr, D);
      cerr << endl;
    }

    // Print the sparse matrix itself.
    if (extraDebug) {
      if (myRank == 0) {
        cerr << "The sparse matrix itself:" << endl;
      }
      writer_type::writeSparse (cerr, matrix);
      cerr << endl;
    }

    RCP<vec_type> vec = rcp (new vec_type (map));
    vec->putScalar (ONE);
    if (myRank == 0) {
      cerr << "Test that vec->putScalar(ONE) filled vec with ones:" << endl;
    }
    comm->barrier ();
    {
      const MT N = static_cast<MT> (vec->getGlobalLength ());
      const MT expectedNorm2 = STM::squareroot (N);
      const MT actualNorm2 = vec->norm2 ();
      TEST_EQUALITY( actualNorm2, expectedNorm2 );

      std::ostringstream os;
      os << "  Proc " << myRank << ": ";
      if (actualNorm2 != expectedNorm2) {
        os << "2-norm INCORRECT: {expectedNorm2: " << expectedNorm2
           << ", actualNorm2: " << actualNorm2 << "}" << endl;
      } else {
        os << "2-norm Correct" << endl;
      }
      cerr << os.str ();
    }

    comm->barrier ();
    {
      const MT N = static_cast<MT> (vec->getGlobalLength ());
      const MT expectedNorm1 = N;
      const MT actualNorm1 = vec->norm1 ();
      TEST_EQUALITY( actualNorm1, expectedNorm1 );

      std::ostringstream os;
      os << "  Proc " << myRank << ": ";
      if (actualNorm1 != expectedNorm1) {
        os << "1-norm INCORRECT: {expectedNorm1: " << expectedNorm1
             << ", actualNorm1: " << actualNorm1 << "}" << endl;
      } else {
        os << "1-norm Correct" << endl;
      }
      cerr << os.str ();
    }

    RCP<const map_type> rangeMap = matrix->getRangeMap ();
    TEST_ASSERT( ! rangeMap.is_null () );
    RCP<vec_type> vec_sol = rcp (new vec_type (rangeMap));
    vec_sol->putScalar (ZERO);

    if (myRank == 0) {
      cerr << "Test that vec_sol->putScalar(ZERO) filled vec with zeros:"
           << endl;
    }
    comm->barrier ();
    {
      const MT expectedNorm2 = STM::zero ();
      const MT actualNorm2 = vec_sol->norm2 ();
      TEST_EQUALITY( actualNorm2, expectedNorm2 );

      std::ostringstream os;
      os << "  Proc " << myRank << ": ";
      if (actualNorm2 != expectedNorm2) {
        os << "2-norm INCORRECT: {expectedNorm2: " << expectedNorm2
           << ", actualNorm2: " << actualNorm2 << "}" << endl;
      } else {
        os << "2-norm Correct" << endl;
      }
      cerr << os.str ();
    }

    comm->barrier ();
    {
      const MT expectedNorm1 = STM::zero ();
      const MT actualNorm1 = vec_sol->norm1 ();
      TEST_EQUALITY( actualNorm1, expectedNorm1 );

      std::ostringstream os;
      os << "  Proc " << myRank << ": ";
      if (actualNorm1 != expectedNorm1) {
        os << "1-norm INCORRECT: {expectedNorm1: " << expectedNorm1
             << ", actualNorm1: " << actualNorm1 << "}" << endl;
      } else {
        os << "1-norm Correct" << endl;
      }
      cerr << os.str ();
    }
    comm->barrier ();

    // Compute vec_sol := matrix*vec.  The result _should_ be a vector
    // of ones everywhere, except for the entry at local index zero
    // (on every process), which should be FIVE.
    if (myRank == 0) {
      cerr << "Apply the matrix to vec, writing the result in vec_sol" << endl;
    }
    //matrix->apply (*vec, *vec_sol, Teuchos::NO_TRANS, ONE, ZERO);
    matrix->apply (*vec, *vec_sol);

    // All entries of vec_sol should be 1, except for the first local
    // entry on each process, which should be 5.  Test this using the
    // 1-norm (to avoid rounding error issues with square root).
    if (myRank == 0) {
      cerr << "Test the 1-norm of vec_sol:" << endl;
    }
    comm->barrier ();
    {
      MT lclExpectedNorm1 = STS::magnitude (FIVE);
      const size_type numLcl = vec_sol->getLocalLength ();
      for (size_type k = 1; k < numLcl; ++k) {
        lclExpectedNorm1 += STS::magnitude (ONE);
      }
      MT expectedNorm1 = STM::zero ();
      reduceAll<int, MT> (*comm, REDUCE_SUM, lclExpectedNorm1,
                          outArg (expectedNorm1));

      const MT actualNorm1 = vec_sol->norm1 ();
      TEST_EQUALITY( actualNorm1, expectedNorm1 );

      std::ostringstream os;
      os << "  Proc " << myRank << ": ";
      if (actualNorm1 != expectedNorm1) {
        os << "INCORRECT: {expectedNorm1: " << expectedNorm1
           << ", actualNorm1: " << actualNorm1 << "}" << endl;
      } else {
        os << "  Correct" << endl;
      }
      cerr << os.str ();
    }

    if (myRank == 0) {
      cerr << "Test the 2-norm of vec_sol:" << endl;
    }
    comm->barrier ();
    {
      MT lclExpectedNorm2 = STS::magnitude (FIVE) * STS::magnitude (FIVE);
      const size_type numLcl = vec_sol->getLocalLength ();
      for (size_type k = 1; k < numLcl; ++k) {
        lclExpectedNorm2 += STS::magnitude (ONE) * STS::magnitude (ONE);
      }
      MT expectedNorm2 = STM::zero ();
      reduceAll<int, MT> (*comm, REDUCE_SUM, lclExpectedNorm2,
                          outArg (expectedNorm2));
      expectedNorm2 = STM::squareroot (expectedNorm2);

      const MT actualNorm2 = vec_sol->norm2 ();
      TEST_EQUALITY( actualNorm2, expectedNorm2 );

      std::ostringstream os;
      os << "  Proc " << myRank << ": ";
      if (actualNorm2 != expectedNorm2) {
        os << "INCORRECT: {expectedNorm2: " << expectedNorm2
           << ", actualNorm2: " << actualNorm2 << "}" << endl;
      } else {
        os << "Correct" << endl;
      }
      cerr << os.str ();
    }

    if (myRank == 0) {
      cerr << "Test the entries of vec_sol" << endl;
    }
    comm->barrier ();
    if (rangeMap->getNodeNumElements () > 0) {
      // Test this both for a const view and for a nonconst view.
      // This may also be a test for {T,X}petra::MultiVector::getData
      // and {T,X}petra::MultiVector::getDataNonConst.

      // Create the const view.
      Teuchos::ArrayRCP<const Scalar> outData = vec_sol->getData (0);
      TEST_ASSERT( static_cast<size_t> (outData.size ()) == rangeMap->getNodeNumElements () );
      if (static_cast<size_t> (outData.size ()) == rangeMap->getNodeNumElements () &&
          outData.size () > static_cast<size_type> (0)) {
        TEST_EQUALITY( outData[0], FIVE );
        if (outData[0] != FIVE) {
          std::ostringstream os;
          os << "Proc " << myRank << ": outData[0] = " << outData[0]
             << " != FIVE" << endl;
          cerr << os.str ();
        }
      }
      if (rangeMap->getNodeNumElements () > static_cast<size_t> (1)) {
        bool allOnes = true;
        for (size_t k = 1; k < rangeMap->getNodeNumElements (); ++k) {
          if (outData[k] != ONE) {
            allOnes = false;
          }
        }
        TEST_ASSERT( allOnes );
      }

      // Invalidate the const view, before creating a nonconst view.
      outData = Teuchos::null;
      // Create the nonconst view.
      Teuchos::ArrayRCP<Scalar> outDataNonConst = vec_sol->getDataNonConst (0);
      TEST_ASSERT( static_cast<size_t> (outDataNonConst.size ()) == rangeMap->getNodeNumElements () );
      if (static_cast<size_t> (outDataNonConst.size ()) == rangeMap->getNodeNumElements () &&
          outDataNonConst.size () > static_cast<size_type> (0)) {
        TEST_EQUALITY( outDataNonConst[0], FIVE );
      }
      if (rangeMap->getNodeNumElements () > static_cast<size_t> (1)) {
        bool allOnes = true;
        for (size_t k = 1; k < rangeMap->getNodeNumElements (); ++k) {
          if (outDataNonConst[k] != ONE) {
            allOnes = false;
          }
        }
        TEST_ASSERT( allOnes );
      }
    }

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (myRank == 0) {
      cerr << "Construct the expected answer vector" << endl;
    }

    RCP<vec_type> vectest = rcp (new vec_type (map));
    {
      vectest->putScalar (ONE);
      Teuchos::ArrayRCP<Scalar> vectestData = vectest->getDataNonConst (0);
      vectestData[0] = FIVE;
      vectestData = Teuchos::null; // as the semantics require

      if (myRank == 0) {
        cerr << "Test the expected answer vector's 1-norm:" << endl;
      }
      comm->barrier ();
      {
        MT lclExpectedNorm1 = STS::magnitude (FIVE);
        const size_type numLcl = vec_sol->getLocalLength ();
        for (size_type k = 1; k < numLcl; ++k) {
          lclExpectedNorm1 += STS::magnitude (ONE);
        }
        MT expectedNorm1 = STM::zero ();
        reduceAll<int, MT> (*comm, REDUCE_SUM, lclExpectedNorm1,
                            outArg (expectedNorm1));
        const MT actualNorm1 = vectest->norm1 ();
        TEST_EQUALITY( actualNorm1, expectedNorm1 );

        std::ostringstream os;
        os << "  Proc " << myRank << ": ";
        if (actualNorm1 != expectedNorm1) {
          os << "INCORRECT: {expectedNorm1: " << expectedNorm1
             << ", actualNorm1: " << actualNorm1 << "}" << endl;
        } else {
          os << "Correct" << endl;
        }
        cerr << os.str ();
      }
    }

    // vec_sol := vec_sol - vectest.
    vec_sol->update (-ONE, *vectest, ONE);

    if (myRank == 0) {
      cerr << "Test the solution error's 1-norm:" << endl;
    }
    comm->barrier ();
    {
      const MT expectedNorm1 = STM::zero ();
      const MT actualNorm1 = vec_sol->norm1 ();
      TEST_EQUALITY( actualNorm1, expectedNorm1 );

      std::ostringstream os;
      os << "  Proc " << myRank << ": ";
      if (actualNorm1 != expectedNorm1) {
        os << "INCORRECT: {expectedNorm1: " << expectedNorm1
           << ", actualNorm1: " << actualNorm1 << "}" << endl;
      } else {
        os << "Correct" << endl;
      }
      cerr << os.str ();
    }

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (myRank == 0) {
      cerr << "All done!  Test " << (gblSuccess ? "succeeded" : "FAILED") << "."
           << endl;
    }
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ReplaceLocalValues, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

} // namespace (anonymous)

