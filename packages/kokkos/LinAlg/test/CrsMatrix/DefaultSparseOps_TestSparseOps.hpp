/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER
*/

#ifndef __CrsMatrix_TestSparseOps_hpp
#define __CrsMatrix_TestSparseOps_hpp

#include <Teuchos_MatrixMarket_Raw_Reader.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_BLAS.hpp>
#include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_DefaultArithmetic.hpp>

/// \file CrsMatrix_TestSparseOps.hpp
/// \brief Header file with helper functions for testing Kokkos sparse kernels.
///

/// \class TestSparseOps
/// \brief Helper functions for tests of Kokkos sparse kernels.
///
/// \tparam SparseOpsType Any class which implements the Kokkos sparse
///   kernels interface, for which DefaultHostSparseOps is a working
///   example and EmptySparseKernel (in
///   kokkos/LinAlg/examples/KokkosExamples_EmptySparseKernelClass.hpp)
///   is a stub.
template<class SparseOpsType>
class TestSparseOps {
public:
  typedef SparseOpsType sparse_ops_type;
  typedef typename sparse_ops_type::scalar_type scalar_type;
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;
  typedef typename sparse_ops_type::ordinal_type ordinal_type;
  typedef typename sparse_ops_type::node_type node_type;

  typedef typename sparse_ops_type::template graph<ordinal_type, node_type>::graph_type graph_type;
  typedef typename sparse_ops_type::template matrix<scalar_type, ordinal_type, node_type>::matrix_type matrix_type;

  /// \brief Read a CSR-format sparse matrix from a Matrix Market file.
  ///
  /// CSR stands for "compressed sparse row."  It's the desired input
  /// format for Kokkos' sparse kernels interface.
  ///
  /// \param numRows [out] The number of rows in the sparse matrix.
  /// \param numCols [out] The number of columns in the sparse matrix.
  /// \param rowptr [out] The first of the three CSR arrays.  For row
  ///   r (zero-based row indices), rowptr[r] .. rowptr[r+1]-1 give
  ///   the index range of colind and values for the entries of that
  ///   row.
  /// \param colind [out] The second of the three CSR arrays; the column
  ///   indices.
  /// \param values [out] The third of the three CSR arrays; the values of
  ///   the matrix.
  /// \param filename [in] The name of the Matrix Market - format file to read.
  void
  readFile (size_t& numRows,
            size_t& numCols,
            Teuchos::ArrayRCP<const size_t>& rowptr,
            Teuchos::ArrayRCP<const ordinal_type>& colind,
            Teuchos::ArrayRCP<const scalar_type>& values,
            const std::string& filename) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::as;
    using Teuchos::arcp;
    using Teuchos::arcp_const_cast;
    using Teuchos::null;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // The reader wants ptr to have the same type of entries as ind.
    // We'll copy right before leaving the routine.
    ArrayRCP<ordinal_type> ptr;
    ArrayRCP<ordinal_type> ind;
    ArrayRCP<scalar_type>  val;
    ordinal_type nrows = 0, ncols = 0;

    Teuchos::MatrixMarket::Raw::Reader<scalar_type, ordinal_type> reader;

    // In "intolerant" mode, this will throw an exception if there is
    // a syntax error in the file.
    (void) reader.readFile (ptr, ind, val, nrows, ncols, filename);

    typedef ArrayRCP<size_t>::size_type size_type;
    ArrayRCP<size_t> ptrout (static_cast<size_type> (nrows + 1));
    for (size_type k = 0; k <= nrows; ++k) {
      ptrout[k] = as<size_t> (ptr[k]);
    }
    // Now we're done with ptr.
    ptr = null;

    // Assign the output arguments.
    numRows = as<size_t> (nrows);
    numCols = as<size_t> (ncols);
    rowptr = arcp_const_cast<const size_t> (ptrout);
    colind = arcp_const_cast<const ordinal_type> (ind);
    values = arcp_const_cast<const scalar_type> (val);
  }

  /// \brief Initialize and return a sparse kernels object.
  ///
  /// A Kokkos sparse kernels object turns a sparse matrix
  /// (represented in compressed sparse row format, more or less) into
  /// an opaque implementation of sparse matrix-(multi)vector multiply
  /// and sparse triangular solve.
  ///
  /// \param node [in/out] Kokkos Node instance, which the
  ///   constructors of graph_type and SparseOpsType require.
  /// \param ptr [in] The first of the three CSR arrays.  For row r
  ///   (zero-based row indices), ptr[r] .. ptr[r+1]-1 give the index
  ///   range of ind and val for the entries of that row.
  /// \param ind [in] The second of the three CSR arrays; the column
  ///   indices.
  /// \param val [in] The third of the three CSR arrays; the values of
  ///   the matrix.
  ///
  /// After calling this method, you can set ptr, ind, and val to
  /// null.  This may free memory if the SparseOpsType copies into its
  /// own internal format instead of just using the original arrays.
  Teuchos::RCP<SparseOpsType>
  makeSparseOps (const Teuchos::RCP<node_type>& node,
                 const Teuchos::ArrayRCP<const size_t>& ptr,
                 const Teuchos::ArrayRCP<const ordinal_type>& ind,
                 const Teuchos::ArrayRCP<const scalar_type>& val) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::null;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    const size_t numRows = static_cast<size_t> (ptr.size() == 0 ? 0 : ptr.size() - 1);
    RCP<ParameterList> graphParams = parameterList ("Graph");
    RCP<graph_type> graph = rcp (new graph_type (numRows, node, graphParams));

    graph->setStructure (ptr, ind);
    RCP<ParameterList> matParams = parameterList ("Matrix");
    RCP<matrix_type> matrix = rcp (new matrix_type (graph, matParams));
    matrix->setValues (val);

    Teuchos::EUplo uplo = Teuchos::UNDEF_TRI;
    Teuchos::EDiag diag = Teuchos::NON_UNIT_DIAG;
    RCP<ParameterList> finParams = parameterList ("Finalize");
    SparseOpsType::finalizeGraphAndMatrix (uplo, diag, *graph, *matrix, finParams);

    RCP<SparseOpsType> ops = rcp (new SparseOpsType (node));
    ops->setGraphAndMatrix (graph, matrix);
    return ops;
  }

  /// \brief Convert the given dense triangular matrix to a sparse kernels object.
  ///
  /// \param A [in] The dense (lower or upper) triangular matrix.
  /// \param uplo [in] Whether A is lower (LOWER_TRI) or upper
  ///   (UPPER_TRI) triangular.
  /// \param diag [in] Whether A has an implicit unit diagonal
  ///   (UNIT_DIAG) or not (NON_UNIT_DIAG).
  Teuchos::RCP<SparseOpsType>
  denseTriToSparseOps (const Teuchos::SerialDenseMatrix<int, scalar_type>& A,
                       const Teuchos::RCP<node_type>& node,
                       const Teuchos::EUplo uplo,
                       const Teuchos::EDiag diag) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::arcp;
    using Teuchos::arcp_const_cast;
    using Teuchos::as;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;

    const ordinal_type N = A.numRows ();
    // If we're not storing the diagonal entries, subtract N off the
    // total number of entries to store.
    const ordinal_type NNZ = diag == Teuchos::UNIT_DIAG ?
      (N*(N-1)) / 2 : // UNIT_DIAG
      (N*(N+1)) / 2;  // NON_UNIT_DIAG

    ArrayRCP<size_t> ptr (N+1);
    ArrayRCP<ordinal_type> ind (NNZ);
    ArrayRCP<scalar_type> val (NNZ);

    ordinal_type counter = 0;
    for (ordinal_type i = 0; i < N; ++i) {
      ptr[i] = counter;
      if (uplo == Teuchos::UPPER_TRI) {
        const ordinal_type lowerBound = (diag == Teuchos::UNIT_DIAG) ? i+1 : i;
        for (ordinal_type j = lowerBound; j < N; ++j) {
          ind[counter] = j;
          val[counter] = A(i, j);
          ++counter;
        }
      }
      else { // uplo == Teuchos::LOWER_TRI
        const ordinal_type upperBound = (diag == Teuchos::UNIT_DIAG) ? i : i+1;
        for (ordinal_type j = 0; j < upperBound; ++j) {
          ind[counter] = j;
          val[counter] = A(i, j);
          ++counter;
        }
      }
    }
    ptr[N] = counter;

    TEUCHOS_TEST_FOR_EXCEPTION(counter != NNZ, std::logic_error,
      "TestSparseOps::denseTriToSparseOps: Expected " << NNZ << " entries in "
      "the sparse matrix, but got " << counter << " instead.  Please report "
      "this bug (in tests) to the Kokkos developers.");

    RCP<graph_type> graph =
      rcp (new graph_type (as<size_t> (N), node, null));
    graph->setStructure (arcp_const_cast<const size_t> (ptr),
                         arcp_const_cast<const ordinal_type> (ind));
    RCP<matrix_type> matrix =
      rcp (new matrix_type (rcp_const_cast<const graph_type> (graph), null));
    matrix->setValues (arcp_const_cast<const scalar_type> (val));

    SparseOpsType::finalizeGraphAndMatrix (uplo, diag, *graph, *matrix, null);
    RCP<SparseOpsType> ops = rcp (new SparseOpsType (node));
    ops->setGraphAndMatrix (graph, matrix);
    return ops;
  }

  /// Return an initialized Kokkos::MultiVector, filled with zeros.
  ///
  /// \param node [in] The Kokkos Node instance.
  /// \param numRows [in] The number of rows in the MultiVector.
  /// \param numCols [in] The number of columns in the MultiVector.
  Teuchos::RCP<Kokkos::MultiVector<scalar_type, node_type> >
  makeMultiVector (const Teuchos::RCP<node_type>& node,
                   const ordinal_type numRows,
                   const ordinal_type numCols) const
  {
    using Teuchos::arcp;
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Kokkos::MultiVector<scalar_type, node_type> MV;
    typedef Kokkos::DefaultArithmetic<MV> MVT;

    RCP<MV> X = rcp (new MV (node));
    X->initializeValues (as<size_t> (numRows),
                         as<size_t> (numCols),
                         arcp<scalar_type> (numRows*numCols),
                         as<size_t> (numRows));
    MVT::Init (*X, STS::zero ());
    return X;
  }

  /// \brief Return the maximum relative error between the two multivectors.
  ///
  /// We define "maximum relative error" between X and Y as
  /// \f\[
  ///   \max_i \| X_i - Y_i \|-2 / \| X_i \|_2,
  /// \f\]
  /// where \f$X_i\f$ indicates the i-th column of X.
  ///
  /// \param X [in] The "correct" multivector, against which to compare Y.
  /// \param Y [in] The "test" multivector, which we hope is equal or
  ///   close to X.
  /// \param Z [in] A "scratch" multivector to use for storing X - Y.
  magnitude_type
  maxRelativeError (const Teuchos::RCP<const Kokkos::MultiVector<scalar_type, node_type> >& X,
                    const Teuchos::RCP<const Kokkos::MultiVector<scalar_type, node_type> >& Y,
                    const Teuchos::RCP<Kokkos::MultiVector<scalar_type, node_type> >& Z) const
  {
    using Teuchos::Array;
    using Teuchos::as;
    typedef typename Array<magnitude_type>::size_type size_type;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;
    typedef Kokkos::MultiVector<scalar_type, node_type> MV;
    typedef Kokkos::DefaultArithmetic<MV> MVT;

    const ordinal_type numCols = as<ordinal_type> (X->getNumCols ());
    Array<magnitude_type> numerNorms (numCols);
    Array<magnitude_type> denomNorms (numCols);

    MVT::Assign (*Z, (const MV) *Y); // Z := Y
    MVT::GESUM (*Z, STS::one(), (const MV) *X, -STS::one()); // Z := X - Z
    MVT::NormInf ((const MV) *Z, numerNorms ());
    MVT::NormInf ((const MV) *X, denomNorms ());

    magnitude_type maxRelNorm = STM::zero();
    for (size_type k = 0; k < numCols; ++k) {
      // If the norm of the current column of X is zero, use the absolute error.
      const magnitude_type relNorm = (denomNorms[k] == STM::zero ()) ?
        numerNorms[k] :
        numerNorms[k] / denomNorms[k];
      if (relNorm > maxRelNorm) {
        maxRelNorm = relNorm;
      }
    }
    return maxRelNorm;
  }

  /// \brief Test sparse matrix-(multi)vector multiply and sparse triangular solve.
  ///
  /// \param node [in/out] The Kokkos Node instance.
  /// \param N [in] The number of rows (and columns) in the sparse
  ///   matrices to test.
  /// \param tol [in] Tolerance for relative errors.
  ///
  /// Test methodology
  /// ================
  ///
  /// Summary: Imitate the L and U factors produced by an LU
  /// factorization (no pivoting).  Use the factors to test sparse
  /// matrix-vector multiply and sparse triangular solve.
  ///
  /// 1. Make random dense $\hat{L}$ and $\hat{U}$ matrices.
  ///    $\hat{L}$ is upper triangular with unit diagonal, and
  ///    $\hat{U}$ is lower triangular.
  /// 2. Compute a random multivector X.  Give X > 1 column, to make
  ///    sure that the routine can handle multivectors correctly.
  /// 3. Compute $\hat{Y} = \hat{L} (\hat{U} X)$ using the BLAS.
  /// 4. Convert $\hat{L}$ to L and $\hat{U}$ to U, where L and U are
  ///    "sparse" matrices.  We put "sparse" in quotes because they
  ///    are actually dense, but stored sparsely.
  /// 5. Test sparse matrix-(multi)vector multiply by computing $Y =
  ///    L (U X)$ and making sure that $\|Y - \hat{Y}\| / \|\hat{Y}\|
  ///    \geq \tau$ for some reasonable tolerance $\tau$.
  /// 6. Test sparse triangular solve:
  ///    a. Compute $\hat{Z} = \hat{L}^{-1} \hat{Y}$ and $Z = L^{-1}
  ///       \hat{Y}$, and compare the results.
  ///    b. Compute $\hat{W} = \hat{U}^{-1} \hat{Z}$ and $W = L^{-1}
  ///       \hat{Z}$, and compare the results.
  ///
  /// Possible issues
  /// ===============
  ///
  /// Random square matrices tend to be well conditioned on average
  /// (this statement sounds vague, but can be made rigorous).
  /// However, that doesn't say anything about the product of random
  /// triangular matrices L and U.  That means the test is not
  /// guaranteed to succeed.
  void
  testSparseOps (const Teuchos::RCP<node_type>& node,
                 const ordinal_type N,
                 const magnitude_type tol) const
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::LEFT_SIDE;
    using Teuchos::UPPER_TRI;
    using Teuchos::NO_TRANS;
    using Teuchos::NON_UNIT_DIAG;
    using Teuchos::LOWER_TRI;
    using Teuchos::UNIT_DIAG;
    typedef Array<size_t>::size_type size_type;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;
    typedef Kokkos::MultiVector<scalar_type, node_type> MV;
    typedef Kokkos::DefaultArithmetic<MV> MVT;
    // Teuchos' BLAS and LAPACK wrappers do accept an "OrdinalType"
    // template parameter, but only work for int (unless you've build
    // your BLAS and LAPACK with 64-bit integer indices).  Thus, in
    // turn, we use int for SDM's ordinal_type for compatibility.
    typedef Teuchos::SerialDenseMatrix<int, scalar_type> dense_matrix_type;
    typedef Teuchos::BLAS<int, scalar_type> blas_type;
    typedef Teuchos::LAPACK<int, scalar_type> lapack_type;

    // A record of error messages reported by any failed tests.
    // We run _all_ the tests first, then report any failures.
    std::ostringstream err;
    // Whether any tests have failed thus far.
    bool success = true;

    // Relative error of the current operation.
    magnitude_type relErr = STM::zero ();

    dense_matrix_type L_dense (N, N), U_dense (N, N);
    L_dense.random ();
    U_dense.random ();

    // Make the diagonal of U larger and positive, to ensure that U is
    // nonsingular and well conditioned.
    for (ordinal_type i = 0; i < N; ++i) {
      const scalar_type ten = as<scalar_type> (10);
      U_dense(i, i) = STS::magnitude(U_dense(i,i)) * ten;
    }

    // Convert L and U into separate sparse matrices.
    RCP<SparseOpsType> L_sparse =
      denseTriToSparseOps (L_dense, node, LOWER_TRI, UNIT_DIAG);
    RCP<SparseOpsType> U_sparse =
      denseTriToSparseOps (U_dense, node, UPPER_TRI, NON_UNIT_DIAG);

    // Compute a random input multivector.  Don't give it more columns
    // than N, but do give it more than one column if the size of the
    // matrices allows it.
    const ordinal_type numColsX = std::min (N, 5);
    RCP<MV> X = makeMultiVector (node, N, numColsX);
    MVT::Random (*X);

    // We make MVs both for results and for scratch space (since the
    // BLAS' dense triangular solve overwrites its input).
    RCP<MV> Y     = makeMultiVector (node, N, numColsX);
    RCP<MV> Y_hat = makeMultiVector (node, N, numColsX);
    RCP<MV> Z     = makeMultiVector (node, N, numColsX);
    RCP<MV> Z_hat = makeMultiVector (node, N, numColsX);
    RCP<MV> W     = makeMultiVector (node, N, numColsX);
    RCP<MV> W_hat = makeMultiVector (node, N, numColsX);

    // Compute Z_hat := U_dense * X.  First copy X into Z_hat, since TRMM
    // overwrites its input.
    MVT::Assign (*Z_hat, (const MV) *X);
    blas_type blas;
    blas.TRMM (LEFT_SIDE, UPPER_TRI, NO_TRANS, NON_UNIT_DIAG, N, numColsX,
               STS::one (), U_dense.values (), U_dense.stride (),
               Z_hat->getValuesNonConst ().getRawPtr (),
               as<int> (Z_hat->getStride ()));
    // Compute Z := U_sparse * X.
    U_sparse->template multiply<scalar_type, scalar_type> (NO_TRANS, STS::one (),
                                                           (const MV) *X, *Z);
    // Compare Z and Z_hat.  Use W as scratch space.
    relErr = maxRelativeError (Z_hat, Z, W);
    if (relErr > tol) {
      err << "Sparse matrix-(multi)vector multiply test failed for upper "
        "triangular matrix U.  Maximum relative error " << relErr << " exceeds "
        "the given tolerance " << tol << ".\n";
      success = false;
    }

    // Compute Y_hat := L_dense * Z_hat.  First copy Z_hat into Y_hat,
    // since TRMM overwrites its input.
    MVT::Assign (*Y_hat, (const MV) *Z_hat);
    blas.TRMM (LEFT_SIDE, LOWER_TRI, NO_TRANS, UNIT_DIAG, N, numColsX,
               STS::one (), L_dense.values (), L_dense.stride (),
               Y_hat->getValuesNonConst ().getRawPtr (),
               as<int> (Y_hat->getStride ()));

    // FIXME (mfh 19 June 2012) There's a bit of a problem with
    // DefaultHostSparseOps: its multiply() method doesn't appear to
    // respect the implicit unit diagonal indication.  Thus, for now,
    // we don't run this test.
    if (false) {
      // Compute Y = L_sparse * Z.
      L_sparse->template multiply<scalar_type, scalar_type> (NO_TRANS, STS::one (),
                                                             (const MV) *Z, *Y);
      // Compare Y and Y_hat.  Use W as scratch space.
      relErr = maxRelativeError (Y_hat, Y, W);
      if (relErr > tol) {
        err << "Sparse matrix-(multi)vector multiply test failed for lower "
          "triangular matrix L with implicit unit diagonal.  Maximum relative "
          "error " << relErr << " exceeds the given tolerance " << tol << ".\n";
        success = false;
      }
    }

    //
    // Test sparse triangular solve.
    //

    // Compute Z_hat = L_dense \ Y_hat.
    MVT::Assign (*Z_hat, (const MV) *Y_hat);
    blas.TRSM (LEFT_SIDE, LOWER_TRI, NO_TRANS, UNIT_DIAG, N, numColsX,
               STS::one (), L_dense.values (), L_dense.stride (),
               Z_hat->getValuesNonConst ().getRawPtr (),
               as<int> (Z_hat->getStride ()));
    // Compute Z = L_sparse \ Y_hat.
    L_sparse->template solve<scalar_type, scalar_type> (NO_TRANS,
                                                        (const MV) *Y_hat, *Z);
    // Compare Z and Z_hat.  Use W as scratch space.
    relErr = maxRelativeError (Z_hat, Z, W);
    if (relErr > tol) {
      err << "Sparse triangular solve test failed for lower triangular matrix "
        "L with implicit unit diagonal.  Maximum relative error " << relErr
          << " exceeds the given tolerance " << tol << ".\n";
      success = false;
    }

    // Compute W_hat = U_dense \ Z_hat.
    MVT::Assign (*W_hat, (const MV) *Z_hat);
    blas.TRSM (LEFT_SIDE, UPPER_TRI, NO_TRANS, NON_UNIT_DIAG, N, numColsX,
               STS::one (), U_dense.values (), U_dense.stride (),
               W_hat->getValuesNonConst ().getRawPtr (),
               as<int> (W_hat->getStride ()));
    // Compute W = U_sparse \ Z_hat.
    U_sparse->template solve<scalar_type, scalar_type> (NO_TRANS,
                                                        (const MV) *Z_hat, *W);
    // Compare W and W_hat.  Use Z as scratch space.
    relErr = maxRelativeError (W_hat, W, Z);
    if (relErr > tol) {
      err << "Sparse triangular solve test failed for upper triangular matrix "
        "U.  Maximum relative error " << relErr << " exceeds the given "
        "tolerance " << tol << ".\n";
      success = false;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(! success, std::runtime_error,
      "One or more sparse ops tests failed.  Here is the full report:\n"
      << err.str());

    const bool alternateTest = false;
    if (alternateTest) {
      dense_matrix_type A (N, N);
      A.random ();

      lapack_type lapack;
      int info = 0;
      Array<int> ipiv (N);
      lapack.GETRF (N, N, A.values (), A.stride (), ipiv.getRawPtr (), &info);
      TEUCHOS_TEST_FOR_EXCEPTION(info < 0, std::logic_error, "LAPACK's _GETRF "
        "routine reported that its " << (-info) << "-th argument had an illegal "
        "value.  This probably indicates a bug in the way Kokkos is calling the "
        "routine.  Please report this bug to the Kokkos developers.");
      TEUCHOS_TEST_FOR_EXCEPTION(info > 0, std::runtime_error, "LAPACK's _GETRF "
        "routine reported that the " << info << "-th diagonal element of the U "
        "factor is exactly zero.  This indicates that the matrix A is singular.  "
        "A is pseudorandom, so it is possible but unlikely that it actually is "
        "singular.  More likely is that the pseudorandom number generator isn't "
        "working correctly.  This is not a Kokkos bug, but it could be a Teuchos "
        "bug, since Teuchos is invoking the generator.");
    }
  }
};

#endif // __CrsMatrix_TestSparseOps_hpp
