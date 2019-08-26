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
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();

    if (myImageID == 0) {
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
      const GO gblRowIdx = Teuchos::as<GO>(2*myImageID + lclRowIdx);
      Teuchos::Array<GO> cols = {{gblRowIdx - GO_ONE, gblRowIdx, gblRowIdx + GO_ONE}};

      if((myImageID == 0) && (lclRowIdx == 0)) { // First row of the matrix
        matrix->insertGlobalValues(gblRowIdx, cols(1, 2), vals(1, 2));
      } else if((myImageID == numImages - 1) && (lclRowIdx == 1)) { // Last row of the matrix
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

    if (myImageID == 0) {
      cerr << "The matrix is now a FEM-type tri-diagonal mass matrix. "
           << "Let's replace all diagonal entries by the global ID of the owning proc." << endl;
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
       *    - unchanged off-diagonal values (not implemented yet)
       */

      // Create vector with new diagonal values
      RCP<vec_type> newDiag = rcp(new vec_type(matrix->getRowMap()));
      newDiag->putScalar(rankAsScalar);

      // Replace the diagonal
      LO numReplacedDiagEntries = ::Tpetra::replaceDiagonalCrsMatrix<Scalar,LO,GO,Node>(*matrix, *newDiag);

      // Tests
      {
        /* Test if every row has been touched.
         *
         * Every row has just one diagonal element, so we expect
         * the number of replaced diagonal entries to match the
         * local number of rows.
         */
	const LO lclNumRows = static_cast<LO> (matrix->getNodeNumRows ());
        TEST_EQUALITY(numReplacedDiagEntries, lclNumRows);

        /* Test for successful replacement
         *
         * 1. Extract diagonal copy
         * 2. Test if diagonal element matches rank ID the we intended to set
         */

        vec_type diagCopy (matrix->getRowMap ());
        matrix->getLocalDiagCopy (diagCopy);
	diagCopy.sync_host ();
	auto diagCopyData = diagCopy.getLocalViewHost ();

	using impl_scalar_type = typename vec_type::impl_scalar_type;
	// If Scalar is std::complex<T>, impl_scalar_type is
	// Kokkos::complex<T>.  Otherwise, Scalar and impl_scalar_type
	// are the same.
	const impl_scalar_type rankAsImplScalarType =
	  static_cast<impl_scalar_type> (rankAsScalar);
        for (size_t i = 0; i < diagCopyData.size(); ++i) {
          TEST_EQUALITY_CONST(diagCopyData(i,0), rankAsImplScalarType);
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

    if (myImageID == 0) {
      cerr << "All done!  Test " << (gblSuccess ? "succeeded" : "FAILED") << "."
           << endl;
    }
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ReplaceDiagonal, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

} // namespace (anonymous)

