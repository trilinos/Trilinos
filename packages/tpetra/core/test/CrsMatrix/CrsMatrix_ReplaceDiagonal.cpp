// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <TpetraCore_ETIHelperMacros.h>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

#include <Tpetra_replaceDiagonalCrsMatrix.hpp>

#include <MatrixMarket_Tpetra.hpp>

namespace { // (anonymous)


  // Unit test of replacing the diagonal of a Tpetra::CrsMatrix
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ReplaceDiagonal, Scalar, LO, GO, Node )
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
    using std::size_t;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
    typedef Tpetra::Vector<Scalar, LO, GO, Node> vec_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType MT;
    const Scalar SC_ONE = STS::one();
    typedef Teuchos::OrdinalTraits<LO> LOT;
    const LO LO_INVALID = LOT::invalid();
    const LO LO_ONE = LOT::one();
    const GO GO_ONE = Teuchos::OrdinalTraits<GO>::one();
    int lclSuccess = 0;
    int gblSuccess = 0;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const size_t numProc = comm->getSize();
    const size_t myProc = comm->getRank();

    if (myProc == 0) {
      cerr << "Tpetra replaceDiagonal test" << endl
           << "Create Map and matrix" << endl;
    }

    // create a Map
    RCP<const map_type> map = Tpetra::createContigMapWithNode<LO,GO,Node> (LO_INVALID,
                                                                           LO_ONE + LO_ONE,
                                                                           comm);

    // Create a matrix with at most 3 entries per row
    RCP<crs_matrix_type> matrix = rcp (new crs_matrix_type (map, 3));
    const Scalar rankAsScalar = static_cast<Scalar>(static_cast<MT>(comm->getRank()));

    Teuchos::Array<Scalar> vals = {{SC_ONE, rankAsScalar + SC_ONE, SC_ONE}};
    for(size_t lclRowIdx = 0; lclRowIdx < 2; ++lclRowIdx) {
      const GO gblRowIdx = Teuchos::as<GO>(2*myProc + lclRowIdx);
      Teuchos::Array<GO> cols = {{gblRowIdx - GO_ONE, gblRowIdx, gblRowIdx + GO_ONE}};

      if((myProc == 0) && (lclRowIdx == 0)) { // First row of the matrix
        matrix->insertGlobalValues(gblRowIdx, cols(1, 2), vals(1, 2));
      } else if((myProc == numProc - 1) && (lclRowIdx == 1)) { // Last row of the matrix
        matrix->insertGlobalValues(gblRowIdx, cols(0, 2), vals(0, 2));
      } else {
        matrix->insertGlobalValues(gblRowIdx, cols(), vals());
      }
    }

    matrix->fillComplete();
    TEST_ASSERT(matrix->isFillComplete());

    // Make sure that all processes got this far.
    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      success = success && (gblSuccess == 1);
      TEST_EQUALITY_CONST( gblSuccess, 1 );
    }

    if (myProc == 0) {
      cerr << "The matrix is now a FEM-type tri-diagonal mass matrix. "
           << "Replace all diagonal entries by the rank of the owning proc." 
           << endl;
    }
    comm->barrier ();
    {
      /* Replace the diagonal of the matrix by the ID of the owning MPI rank
       *
       * 1. Create map
       * 2. Create vector with new diagonal values
       * 3. Replace the diagonal
       * 4. Test for
       *    - successful replacement of diagonal values
       *    - unchanged off-diagonal values
       */

      // Create vector with new diagonal values
      RCP<vec_type> newDiag = rcp(new vec_type(matrix->getRowMap()));
      newDiag->putScalar(rankAsScalar);

      // Replace the diagonal
      LO numReplacedDiagEntries = 
         Tpetra::replaceDiagonalCrsMatrix<Scalar,LO,GO,Node>(*matrix, *newDiag);

      // Tests
      {
        /* Test if every row has been touched.
         *
         * Every row has just one diagonal element, so we expect
         * the number of replaced diagonal entries to match the
         * local number of rows.
         */
	const LO lclNumRows = static_cast<LO> (matrix->getLocalNumRows ());
        TEST_EQUALITY(numReplacedDiagEntries, lclNumRows);

        /* Test for successful replacement
         *
         * 1. Extract diagonal copy
         * 2. Test if diagonal element matches rank ID the we intended to set
         */

        vec_type diagCopy (matrix->getRowMap ());
        matrix->getLocalDiagCopy (diagCopy);
	auto diagCopyData = diagCopy.getLocalViewHost(Tpetra::Access::ReadOnly);

	using impl_scalar_type = typename vec_type::impl_scalar_type;
	// If Scalar is std::complex<T>, impl_scalar_type is
	// Kokkos::complex<T>.  Otherwise, Scalar and impl_scalar_type
	// are the same.
	const impl_scalar_type rankAsImplScalarType =
	  static_cast<impl_scalar_type> (rankAsScalar);
        for (size_t i = 0; i < diagCopyData.size(); ++i) {
          TEST_EQUALITY_CONST(diagCopyData(i,0), rankAsImplScalarType);
	}

        /* Test that no other matrix values were changed
         */
        
        for (size_t i = 0; i < matrix->getRowMap()->getLocalNumElements(); i++) {
          typename crs_matrix_type::local_inds_host_view_type lcols;
          typename crs_matrix_type::values_host_view_type lvals;
          matrix->getLocalRowView(static_cast<LO>(i), lcols, lvals);
          GO gI = matrix->getRowMap()->getGlobalElement(i);
          auto j = lcols.size();
          for (j = 0; j < lcols.size(); j++) {
            GO gJ = matrix->getColMap()->getGlobalElement(lcols[j]);
            if (gI == gJ) {
              // Diagonal entry; value should be row GO
              TEST_EQUALITY_CONST(lvals[j], rankAsImplScalarType);
            }
            else {
              // Nondiagonal entry; value should be SC_MONE
              TEST_EQUALITY_CONST(lvals[j], SC_ONE);
            }
          }
        }
      }
    }

    // Make sure that all processes got this far.
    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      success = success && (gblSuccess == 1);
      TEST_EQUALITY_CONST( gblSuccess, 1 );
    }

    if (myProc == 0) {
      cerr << "All done!  Test " << (gblSuccess ? "succeeded" : "FAILED") << "."
           << endl;
    }
  }

  // Unit test of replacing the diagonal of a Tpetra::CrsMatrix with 
  // overlapped row map
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ReplaceDiagonalOverlapRowMap, Scalar, LO, GO, Node )
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
    using std::size_t;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
    typedef Tpetra::Vector<Scalar, LO, GO, Node> vec_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    const Scalar SC_MONE = -STS::one();
    int lclSuccess = 0;
    int gblSuccess = 0;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int nProc = comm->getSize();
    const int myProc = comm->getRank();

    if (myProc == 0) {
      cerr << "Tpetra replaceDiagonalOverlapRowMap test" << endl
           << "Create Map and matrix" << endl;
    }
    
    RCP<crs_matrix_type> matrix;
    LO nMyDiags = 0;
    { // Create a tri-diagonal matrix with a 2D blockwise parallel distribution

      // create 2D processor layout
      int nProcRow = std::sqrt(nProc);
      int nProcCol = nProc / nProcRow;
      while (nProcRow * nProcCol != nProc) {
        nProcRow--;
        nProcCol = nProc / nProcRow;
      }
      if (myProc == 0) {
        cerr << "Tpetra replaceDiagonalOverlapRowMap test" << endl
             << "Processor layout:  " << nProcRow << " x " << nProcCol << endl;
      }
      int myProcRow = myProc % nProcRow;
      int myProcCol = myProc / nProcRow;
  
      // create a Map with overlapped entries 
      const int nRowPerProc = 10;
      const int nRows = nRowPerProc * nProcRow;
      const int nCols = nRows; // square matrix
      GO myFirstRow = myProcRow * nRowPerProc;
      GO myLastRow = myFirstRow + nRowPerProc - 1;
      GO myFirstCol = myProcCol * (nCols / nProcCol);
      GO myLastCol = (myProcCol < nProcCol-1 ? (myProcCol+1)*(nCols/nProcCol) 
                                             : nRows) - 1;

      const int maxNZPerRow = 3;
      Teuchos::Array<int> nNZPerRow(nRowPerProc, 0);
      Teuchos::Array<GO> myRows;
      Teuchos::Array<GO> myJ;

      // Build tridiagonal matrix using 2D processor layout
      for (int i = 0; i < nRows; i++) {
        if (i >= myFirstRow && i <= myLastRow) {
          // Diagonal entry
          if (i >= myFirstCol && i <= myLastCol) {
            nNZPerRow[i-myFirstRow]++;
            myJ.push_back(i);
            nMyDiags++;
          }
          // i+1 entry
          if (i+1 < nRows && i+1 >= myFirstCol && i+1 <= myLastCol) {
            nNZPerRow[i-myFirstRow]++;
            myJ.push_back(i+1);
          }
          // i-1 entry
          if (i-1 >= 0 && i-1 >= myFirstCol && i-1 <= myLastCol) {
            nNZPerRow[i-myFirstRow]++;
            myJ.push_back(i-1);
          }
          if (nNZPerRow[i-myFirstRow]) myRows.push_back(i);
        }
      }
    
      Tpetra::global_size_t dummy = 
              Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
      RCP<const map_type> map = rcp(new map_type(dummy, myRows(), 0, comm));

      matrix = rcp (new crs_matrix_type (map, 3));

      int curJ = 0;
      Teuchos::Array<Scalar> vals(maxNZPerRow, SC_MONE);
      for (int i = 0; i < nRowPerProc; i++) {
        int nnz = nNZPerRow[i];
        if (nnz > 0) 
          matrix->insertGlobalValues(i+myFirstRow, myJ(curJ,nnz), vals(0,nnz));
        curJ += nnz;
      }

      matrix->fillComplete();
      TEST_ASSERT(matrix->isFillComplete());
    }

    // Make sure that all processes got this far.
    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      success = success && (gblSuccess == 1);
      TEST_EQUALITY_CONST( gblSuccess, 1 );
    }

    if (myProc == 0) {
      cerr << "The matrix is now a FEM-type tri-diagonal mass matrix. "
           << "Replace diagonal entries by the global ID of row." << endl;
    }
    comm->barrier ();
    {
      /* Replace the diagonal of the matrix by the row's GNO
       *
       * 2. Create vector with new diagonal values using row map
       * 3. Replace the diagonal
       * 4. Test for
       *    - successful replacement of diagonal values
       *    - unchanged off-diagonal values 
       */

      // Create vector with new diagonal values (row GO)
      RCP<vec_type> newDiag = rcp(new vec_type(matrix->getRowMap()));
      {
        auto newDiagData = newDiag->getLocalViewHost(Tpetra::Access::OverwriteAll);
        for (size_t i = 0; i < newDiag->getLocalLength(); i++) 
          newDiagData(i,0) = newDiag->getMap()->getGlobalElement(i);
      }

      // Replace the diagonal
      LO numReplacedDiagEntries = 
         Tpetra::replaceDiagonalCrsMatrix<Scalar,LO,GO,Node>(*matrix, *newDiag);

      // Tests
      {
        /* Test whether number of replaced diagonal entries is correct. */
        TEST_EQUALITY_CONST(numReplacedDiagEntries, nMyDiags);

        /* Test for successful replacement
         * Diagonal entries should be global ID of row
         * Non-diagonal entries should stil be SC_MONE.
         */
        
	using impl_scalar_type = typename crs_matrix_type::impl_scalar_type;
        for (size_t i = 0; i < matrix->getRowMap()->getLocalNumElements(); i++) {
          typename crs_matrix_type::local_inds_host_view_type lcols;
          typename crs_matrix_type::values_host_view_type lvals;
          matrix->getLocalRowView(static_cast<LO>(i), lcols, lvals);
          GO gI = matrix->getRowMap()->getGlobalElement(i);
          auto j = lcols.size();
          for (j = 0; j < lcols.size(); j++) {
            GO gJ = matrix->getColMap()->getGlobalElement(lcols[j]);
            if (gI == gJ) {
              // Diagonal entry; value should be row GO
	      const impl_scalar_type expected = 
    	            static_cast<impl_scalar_type> (gI);
              TEST_EQUALITY_CONST(lvals[j], expected);
            }
            else {
              // Nondiagonal entry; value should be SC_MONE
              TEST_EQUALITY_CONST(lvals[j], SC_MONE);
            }
          }
        }
      }
    }

    // Make sure that all processes got this far.
    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      success = success && (gblSuccess == 1);
      TEST_EQUALITY_CONST( gblSuccess, 1 );
    }

    if (myProc == 0) {
      cerr << "All done!  Test " << (gblSuccess ? "succeeded" : "FAILED") << "."
           << endl;
    }
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ReplaceDiagonal, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ReplaceDiagonalOverlapRowMap, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

} // namespace (anonymous)

