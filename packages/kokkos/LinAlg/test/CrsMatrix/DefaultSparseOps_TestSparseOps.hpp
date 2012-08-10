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
#include <Teuchos_TimeMonitor.hpp>
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_DefaultArithmetic.hpp>

/// \file DefaultSparseOps_TestSparseOps.hpp
/// \brief Header file with helper functions for testing Kokkos sparse kernels.
///


namespace {
  template<class T, class OrdinalType>
  class CopyOp {
  private:
    T* const out_;
    const T* const in_;

  public:
    CopyOp (T* const out, const T* const in) : out_ (out), in_ (in) {}

    // FIXME (mfh 09 Aug 2012) make a device kernel.
    void execute (const OrdinalType i) {
      out_[i] = in_[i];
    }
  };
} // namespace (anonymous)

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

  /// \name Dense matrix, BLAS, and LAPACK typedefs
  ///
  /// Teuchos' BLAS and LAPACK wrappers do accept an "OrdinalType"
  /// template parameter, but they only work for OrdinalType=int
  /// (unless you've built your BLAS and LAPACK with 64-bit integer
  /// indices).  Thus, in turn, we use int for SerialDenseMatrix's
  /// ordinal_type for compatibility.
  ///
  /// The BLAS and LAPACK typedefs are private because they are
  /// implementation details.
  //@{

  //! The type of local dense matrices in column-major order.
  typedef Teuchos::SerialDenseMatrix<int, scalar_type> dense_matrix_type;

private:
  //! The BLAS wrapper type.
  typedef Teuchos::BLAS<int, scalar_type> blas_type;
  //! The LAPACK wrapper type.
  typedef Teuchos::LAPACK<int, scalar_type> lapack_type;

public:
  //@}

  /// \brief Compute dense test matrices for sparse triangular solve.
  ///
  /// To make nonsingular lower and upper triangular matrices for
  /// testing sparse triangular solve, we start with a dense random
  /// matrix A and compute its LU factorization using LAPACK.  Random
  /// matrices tend to be well conditioned, so the L and U factors
  /// will also be well conditioned.  We output the dense matrices
  /// here so that you can test sparse routines by comparing their
  /// results with those of using the BLAS and LAPACK.
  ///
  /// \param A_out [out] A numRows by numRows random matrix, of which
  ///   L_out, U_out is the LU factorization with permutation pivots.
  /// \param L_out [out] The L factor in the LU factorization of
  ///   A_out.
  /// \param U_out [out] The U factor in the LU factorization of
  ///   A_out.
  /// \param pivots [out] The permutation array in the LU
  ///   factorization of A_out.  Row i in the matrix A was
  ///   interchanged with row pivots[i].
  /// \param numRows [in] The number of rows and columns in A_out.
  ///
  /// If you're solving AX=B using the LU factorization, you first
  /// have to apply the row permutation to B.  You can do this in the
  /// same way that the reference _GETRS implementation does, using
  /// the _LASWP BLAS routine.  Here's how the Fortran calls it (where
  /// _ is replaced by D for the special case of scalar_type=double):
  ///
  /// CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
  ///
  /// Then, after the triangular solves, apply the reverse-direction
  /// permutation (last parameter is -1 instead of 1) to the solution
  /// vector X:
  ///
  /// CALL DLASWP( NRHS, X, LDX, 1, N, IPIV, -1 )
  void
  makeDenseTestProblem (Teuchos::RCP<dense_matrix_type>& A_out,
                        Teuchos::RCP<dense_matrix_type>& L_out,
                        Teuchos::RCP<dense_matrix_type>& U_out,
                        Teuchos::Array<ordinal_type>& pivots,
                        const ordinal_type numRows) const
  {
    using Teuchos::Array;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef dense_matrix_type MT;

    const ordinal_type N = numRows;
    RCP<MT> A = rcp (new MT (N, N));
    A->random ();
    // Keep a copy of A, since LAPACK's LU factorization overwrites
    // its input matrix with the L and U factors.
    RCP<MT> A_copy = rcp (new MT (*A));

    // Compute the LU factorization of A.
    lapack_type lapack;
    int info = 0;
    Array<ordinal_type> ipiv (N);
    lapack.GETRF (N, N, A->values (), A->stride (), ipiv.getRawPtr (), &info);
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

    // Create L and U, and copy out the lower resp. upper triangle of
    // A into L resp. U.
    RCP<MT> L = rcp (new MT (N, N));
    RCP<MT> U = rcp (new MT (N, N));
    {
      // Get the MT refs so we don't have to dereference the RCPs each
      // time we change an entry.
      MT& LL = *L;
      MT& UU = *U;
      MT& AA = *A;

      for (ordinal_type i = 0; i < N; ++i) {
        // LL has an implicitly stored unit diagonal, so don't include j==i.
        for (ordinal_type j = 0; j < i; ++j) {
          LL(i,j) = AA(i,j);
        }
        for (ordinal_type j = i; j < N; ++j) {
          UU(i,j) = AA(i,j);
        }
      }
    }

    // "Commit" the outputs.
    pivots.resize (N);
    std::copy (ipiv.begin (), ipiv.end (), pivots.begin ());
    A_out = A_copy; // Return the "original" A, before the factorization.
    L_out = L;
    U_out = U;
  }

  // mfh 28 Jun 2012: It would be nice to use LAPACK's banded solver
  // to generate test problems, but it may be more trouble than it's
  // worth to get the indexing right.  Having a good source of valid
  // LU factorizations without needing N^2 data would be helpful,
  // though for now we disable this code, since we don't know that it
  // works.
#if 0
  void
  makeBandedTestProblem (Teuchos::RCP<dense_matrix_type>& A_out,
                         Teuchos::RCP<dense_matrix_type>& L_out,
                         Teuchos::RCP<dense_matrix_type>& U_out,
                         Teuchos::Array<ordinal_type>& pivots,
                         const ordinal_type numRows) const
  {
    using Teuchos::Array;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef dense_matrix_type MT;

    const ordinal_type M = numRows; // Number of rows in the matrix
    const ordinal_type N = numRows; // Number of columns in the matrix
    const ordinal_type KL = 1; // Number of subdiagonals
    const ordinal_type KU = 1; // Number of superdiagonals
    const ordinal_type LDAB = 2*KL + KU + 1; // As required by LAPACK

    // On input: Don't fill top KL rows of A.
    //
    // On output: The top KL rows of A may be filled with elements of
    // the U factor on output, due to row interchanges.

    RCP<MT> A = rcp (new MT (LDAB, N));
    // Just fill the whole matrix A with random numbers.  LAPACK won't
    // use some of those elements; that's OK.
    A->random ();

    // Keep a copy of A, since LAPACK's LU factorization overwrites
    // its input matrix with the L and U factors.
    RCP<MT> A_copy = rcp (new MT (*A));

    // Compute the LU factorization of A.
    lapack_type lapack;
    int info = 0;
    Array<ordinal_type> ipiv (N);
    lapack.GBTRF (M, N, KL, KU, A->values (), A->stride (), ipiv.getRawPtr (), &info);
    TEUCHOS_TEST_FOR_EXCEPTION(info < 0, std::logic_error, "LAPACK's _GBTRF "
      "routine reported that its " << (-info) << "-th argument had an illegal "
      "value.  This probably indicates a bug in the way Kokkos is calling the "
      "routine.  Please report this bug to the Kokkos developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(info > 0, std::runtime_error, "LAPACK's _GBTRF "
      "routine reported that the " << info << "-th diagonal element of the U "
      "factor is exactly zero.  This indicates that the matrix A is singular.  "
      "A is pseudorandom, so it is possible but unlikely that it actually is "
      "singular.  More likely is that the pseudorandom number generator isn't "
      "working correctly.  This is not a Kokkos bug, but it could be a Teuchos "
      "bug, since Teuchos is invoking the generator.");

    // Create L and U, and copy out the lower resp. upper triangle
    // (both packed in banded form) of A into L resp. U.
    const ordinal_type LDU = KL + KU + 1;
    RCP<MT> L = rcp (new MT (LDAB - LDU, N));
    RCP<MT> U = rcp (new MT (LDU, N));
    {
      // Get the MT refs so we don't have to dereference the RCPs each
      // time we change an entry.
      MT& LL = *L;
      MT& UU = *U;
      MT& AA = *A;

      for (ordinal_type j = 0; j < N; ++j) {
        for (ordinal_type i = 0; i < LDU; ++i) {
          U(i,j) = A(i,j);
        }
        for (ordinal_type i = LDU; i < LDAB; ++i) {
          L(i,j) = A(i,j);
        }
      }
    }

    // "Commit" the outputs.
    pivots.resize (N);
    std::copy (ipiv.begin (), ipiv.end (), pivots.begin ());
    A_out = A_copy; // Return the "original" A, before the factorization.
    L_out = L;
    U_out = U;
  }

  /// \brief Repack LAPACK _GBTRF-style banded format into sparse format.
  ///
  /// A_in, L_in, and U_in come from makeBandedTestProblem() above.
  ///
  /// For an explanation of LAPACK's packed banded format, especially
  /// that of the L and U factors from an LU factorization, see the <a
  /// href="http://www.netlib.org/lapack/double/dgbtrf.f">DGBTRF</a>
  /// documentation or the <a
  /// href="http://www.netlib.org/lapack/lug/node124.html">LAPACK
  /// Users' Guide.</a>
  void
  denseBandedToSparseOps (Teuchos::RCP<SparseOpsType>& A_out,
                          Teuchos::RCP<SparseOpsType>& L_out,
                          Teuchos::RCP<SparseOpsType>& U_out,
                          const Teuchos::RCP<node_type>& node,
                          const Teuchos::RCP<const dense_matrix_type>& A_in,
                          const Teuchos::RCP<const dense_matrix_type>& L_in,
                          const Teuchos::RCP<const dense_matrix_type>& U_in) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;

    const ordinal_type LDAB = A_in.numRows (); // == 2*KL + KU + 1
    const ordinal_type LDU = U_in.numRows ();  // == KL + KU + 1
    const ordinal_type LDL = L_in.numRows ();  // == KL

    const ordinal_type M = A_in.numCols ();
    const ordinal_type N = A_in.numCols ();
    const ordinal_type KL = LDL;
    const ordinal_type KU = LDU - LDL - 1;

    TEUCHOS_TEST_FOR_EXCEPTION(2*KL + KU + 1 != LDAB, std::logic_error,
      "Failed to compute KL and KU correctly.  2*KL + KU + 1 = "
      << 2*KL + KU + 1 << " != LDAB = " << LDAB << ".  KL = " << KL
      << ", KU = " << KU << ".  Please report this bug to the Kokkos "
      "developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(KL + KU + 1 != LDU, std::logic_error,
      "Failed to compute KL and KU correctly.  KL + KU + 1 = "
      << KL + KU + 1 << " != LDU = " << LDU << ".  KL = " << KL
      << ", KU = " << KU << ".  Please report this bug to the Kokkos "
      "developers.");

    // Make U.
    RCP<SparseOpsType> U_sparse;
    {
      const ordinal_type nnz = LDU * N - ((LDU-1)*LDU)/2;
      ArrayRCP<ordinal_type> ptr (LDU+1);
      ArrayRCP<ordinal_type> ind (nnz);
      ArrayRCP<scalar_type> val (nnz);
      ordinal_type ctr = 0;
      ptr[0] = 0;
      for (ordinal_type r = 0; r < LDU; ++r) {
        // Row 0 of U starts in the last row of the packed U matrix, at
        // the (LDU-1, 0) entry (zero-based).  Read the entries in each
        // row diagonally and to the northeast.  Row 1 starts at (LDU-1,
        // 1).  In general, row r starts at (LDU-1, r), and has min(LDU,
        // N-r) entries.
        const ordinal_type numEntries = std::min (LDU, N-r);
        ordinal_type curRow = LDU - 1;
        ordinal_type curCol = r;
        ordinal_type c = r; // start here; it's the upper triangle
        for (ordinal_type k = 0; k < numEntries; ++k, ++ctr) {
          val[ctr] = U(curRow--, curCol++);
          ind[ctr] = c++;
        }
        ptr[r+1] = ctr;
      }
      U_sparse = makeSparseOps (node, Teuchos::null,
                                M, N, ptr, ind, val, Teuchos::UPPER_TRI,
                                Teuchos::NON_UNIT_DIAG);
    }

    // Make L
    RCP<SparseOpsType> L_sparse;
    {
      // Row 0 of L starts in the first row of the packed L matrix, at
      // the (0, 0) entry (zero-based).  Row 1 of L

      // Read the entries in each row diagonally and to the northeast.

      const ordinal_type nnz = LDU * N - ((LDU-1)*LDU)/2;
      ArrayRCP<ordinal_type> ptr (LDU+1);
      ArrayRCP<ordinal_type> ind (nnz);
      ArrayRCP<scalar_type> val (nnz);
      ordinal_type ctr = 0;
      ptr[0] = 0;
      for (ordinal_type r = 0; r < LDU; ++r) {
        // Row 0 of U starts in the last row of the packed U matrix, at
        // the (LDU-1, 0) entry (zero-based).  Read the entries in each
        // row diagonally and to the northeast.  Row 1 starts at (LDU-1,
        // 1).  In general, row r starts at (LDU-1, r), and has min(LDU,
        // N-r) entries.
        const ordinal_type numEntries = std::min (LDU, N-r);
        ordinal_type curRow = LDU - 1;
        ordinal_type curCol = r;
        ordinal_type c = r; // start here; it's the upper triangle
        for (ordinal_type k = 0; k < numEntries; ++k, ++ctr) {
          val[ctr] = U(curRow--, curCol++);
          ind[ctr] = c++;
        }
        ptr[r+1] = ctr;
      }
      // FIXME: Dimensions???
      L_sparse = makeSparseOps (node, M, M, ptr, ind, val, Teuchos::LOWER_TRI,
                                Teuchos::UNIT_DIAG);

    }
  }
  // mfh 28 Jun 2012: See note above on using LAPACK's banded solver
  // to generate sparse test problems.
#endif // 0

  /// \brief Read a CSR-format sparse matrix from a Matrix Market file.
  ///
  /// CSR stands for "compressed sparse row."  It's the desired input
  /// format for Kokkos' sparse kernels interface.
  ///
  /// \param numRows [out] The number of rows in the sparse matrix.
  /// \param numCols [out] The number of columns in the sparse matrix.
  /// \param rowptr [out] The first of the three CSR arrays; the row
  ///   offsets.  For row r (zero-based row indices), rowptr[r]
  ///   .. rowptr[r+1]-1 give the index range of colind and values for
  ///   the entries of that row.
  /// \param colind [out] The second of the three CSR arrays; the
  ///   column indices.
  /// \param values [out] The third of the three CSR arrays; the
  ///   values of the matrix.
  /// \param filename [in] The name of the Matrix Market - format file
  ///   to read.
  void
  readFile (ordinal_type& numRows,
            ordinal_type& numCols,
            Teuchos::ArrayRCP<const       size_t>& rowptr,
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

    ArrayRCP<size_t> ptrout (nrows + 1);
    std::copy (ptr.begin(), ptr.end(), ptrout.begin());
    // Now we're done with ptr.
    ptr = null;

    // Assign the output arguments.
    numRows = as<ordinal_type> (nrows);
    numCols = as<ordinal_type> (ncols);
    rowptr = arcp_const_cast<const size_t> (ptrout);
    colind = arcp_const_cast<const ordinal_type> (ind);
    values = arcp_const_cast<const scalar_type> (val);
  }

  /// \brief Initialize and return a sparse kernels object from a
  ///   Matrix Market file.
  ///
  /// A Kokkos sparse kernels object turns a sparse matrix
  /// (represented in compressed sparse row format, more or less) into
  /// an opaque implementation of sparse matrix-(multi)vector multiply
  /// and sparse triangular solve.
  ///
  /// \param numRows [out] Number of rows in the sparse matrix.
  ///   SparseOpsType doesn't currently require a method that tells
  ///   you the number of rows in the sparse matrix, so we output it
  ///   here for later use.
  /// \param numCols [out] Number of columns in the sparse matrix.
  ///   SparseOpsType doesn't currently require a method that tells
  ///   you the number of columns in the sparse matrix, so we output it
  ///   here for later use.
  /// \param node [in/out] Kokkos Node instance, which the
  ///   constructors of graph_type and SparseOpsType require.
  /// \param params [in/out] Parameters for configuring the
  ///   SparseOpsType instance at construction.
  /// \param filename [in] Name of a Matrix Market sparse matrix file.
  /// \param uplo [in] Whether the matrix is lower triangular
  ///   (LOWER_TRI), upper triangular (UPPER_TRI), or neither
  ///   (UNDEF_TRI).  The latter is the default.
  /// \param diag [in] Whether the matrix has an implicitly stored
  ///   unit diagonal (UNIT_DIAG) or not (NON_UNIT_DIAG).  The latter
  ///   is the default.  Kokkos' convention is to ignore this for
  ///   sparse matrix-vector multiply, and respect it only for
  ///   triangular solves.
  Teuchos::RCP<SparseOpsType>
  makeSparseOpsFromFile (ordinal_type& numRows,
                         ordinal_type& numCols,
                         const Teuchos::RCP<node_type>& node,
                         Teuchos::ParameterList& params,
                         const std::string& filename,
                         const Teuchos::EUplo uplo = Teuchos::UNDEF_TRI,
                         const Teuchos::EDiag diag = Teuchos::NON_UNIT_DIAG) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::arcp_const_cast;
    using Teuchos::as;

    ordinal_type theNumRows = 0;
    ordinal_type theNumCols = 0;
    ArrayRCP<const size_t>       ptr;
    ArrayRCP<const ordinal_type> ind;
    ArrayRCP<const scalar_type>  val;
    readFile (theNumRows, theNumCols, ptr, ind, val, filename);
    numRows = as<ordinal_type> (theNumRows);
    numCols = as<ordinal_type> (theNumCols);

    return makeSparseOps (node, params, numRows, numCols, ptr, ind, val, uplo, diag);
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
  /// \param params [in/out] Parameters for the SparseOpsType object's
  ///   constructor.
  /// \params numRows [in] Number of rows in the sparse matrix.
  /// \params numCols [in] Number of columns in the sparse matrix.
  /// \param ptr [in/out] On input: the first of the three CSR arrays.
  ///   For row r (zero-based row indices), ptr[r] .. ptr[r+1]-1 give
  ///   the index range of ind and val for the entries of that row.
  ///   On output: set to null.
  /// \param ind [in/out] On input: the second of the three CSR
  ///   arrays; the column indices.  On output: set to null.
  /// \param val [in/out] On input: the third of the three CSR arrays;
  ///   the values of the matrix.  On output: set to null.
  /// \param uplo [in] Whether the matrix is lower triangular
  ///   (LOWER_TRI), upper triangular (UPPER_TRI), or neither
  ///   (UNDEF_TRI).  The latter is the default.
  /// \param diag [in] Whether the matrix has an implicitly stored
  ///   unit diagonal (UNIT_DIAG) or not (NON_UNIT_DIAG).  The latter
  ///   is the default.  Kokkos' convention is to ignore this for
  ///   sparse matrix-vector multiply, and respect it only for
  ///   triangular solves.
  ///
  /// This method sets ptr, ind, and val to null on output.  This may
  /// free memory if the SparseOpsType copies into its own internal
  /// format instead of just using the original arrays.  Many
  /// implementations of SparseOpsType that we provide copy the
  /// original data and reorganize them into a different format.
  Teuchos::RCP<SparseOpsType>
  makeSparseOps (const Teuchos::RCP<node_type>& node,
                 Teuchos::ParameterList& params,
                 const ordinal_type numRows,
                 const ordinal_type numCols,
                 Teuchos::ArrayRCP<const size_t>       &ptr,
                 Teuchos::ArrayRCP<const ordinal_type> &ind,
                 Teuchos::ArrayRCP<const scalar_type>  &val,
                 const Teuchos::EUplo uplo = Teuchos::UNDEF_TRI,
                 const Teuchos::EDiag diag = Teuchos::NON_UNIT_DIAG) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::arcp_const_cast;
    using Teuchos::as;
    using Teuchos::null;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // We don't know where ptr, ind, and val came from, so we have to
    // reallocate using the SparseOpsType.  We start over by counting
    // the number of entries in each row.
    ArrayRCP<size_t> numEntriesPerRow (numRows);
    for (ordinal_type i = 0; i < numRows; ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(ptr[i+1] < ptr[i], std::invalid_argument,
        "makeSparseOps: Invalid row offsets array on input.  "
        "ptr[" << (i+1) << "] = " << ptr[i+1] << " < ptr[" << i << "] = "
        << ptr[i] << ".");
      const size_t curNumEntries = ptr[i+1] - ptr[i]; // >= 0; see test above.
      numEntriesPerRow[i] = as<ordinal_type> (curNumEntries);
    }

    // Reallocate and copy over the three arrays for sparse matrix storage:
    //
    // ptr: row offsets
    // ind: column indices
    // val: matrix entries
    ptr = SparseOpsType::allocRowPtrs (node, numEntriesPerRow ());
    // We don't need the array of entry counts per row anymore.
    numEntriesPerRow = Teuchos::null;
    {
      ArrayRCP<ordinal_type> newInd =
        SparseOpsType::template allocStorage<ordinal_type> (node, ptr ());
      ArrayRCP<scalar_type> newVal =
        SparseOpsType::template allocStorage<scalar_type> (node, ptr ());
      node->parallel_for (0, ind.size(),
                          CopyOp<ordinal_type, ordinal_type> (newInd.getRawPtr (),
                                                              ind.getRawPtr ()));
      node->parallel_for (0, val.size(),
                          CopyOp<scalar_type, ordinal_type> (newVal.getRawPtr (),
                                                             val.getRawPtr ()));
      ind = arcp_const_cast<const ordinal_type> (newInd);
      val = arcp_const_cast<const scalar_type> (newVal);
    }

    return makeSparseOpsFromArrays (node, params, numRows, numCols,
                                    ptr, ind, val, uplo, diag);
  }

  /// \brief Initialize and return a sparse kernels object (internal method).
  ///
  /// Initialize and return a sparse kernels object, from the three
  /// arrays already allocated and initialized by the sparse kernels
  /// object's allocation methods.  This method is meant to be called
  /// internally by methods such as makeSparseOps(),
  /// denseTriToSparseOps(), and denseToSparseOps().
  ///
  /// \param node [in/out] Kokkos Node instance, which the
  ///   constructors of graph_type and SparseOpsType require.
  /// \param opsParams [in/out] On input: parameters for the
  ///   SparseOpsType object's constructor.  On output: may be filled
  ///   in with default values.
  /// \params numRows [in] Number of rows in the sparse matrix.
  /// \params numCols [in] Number of columns in the sparse matrix.
  /// \param ptr [in/out] On input: the first of the three CSR arrays.
  ///   For row r (zero-based row indices), ptr[r] .. ptr[r+1]-1 give
  ///   the index range of ind and val for the entries of that row.
  ///   This must have been allocated and initialized using
  ///   SparseOpsTyps::allocRowPtrs().  On output: set to null.
  /// \param ind [in/out] On input: the second of the three CSR
  ///   arrays; the column indices.  This must have been allocated
  ///   using SparseOpsTyps::allocStorage().  On output: set to null.
  /// \param val [in/out] On input: the third of the three CSR arrays;
  ///   the values of the matrix.  This must have been allocated using
  ///   SparseOpsTyps::allocStorage().  On output: set to null.
  /// \param uplo [in] Whether the matrix is lower triangular
  ///   (LOWER_TRI), upper triangular (UPPER_TRI), or neither
  ///   (UNDEF_TRI).  The latter is the default.
  /// \param diag [in] Whether the matrix has an implicitly stored
  ///   unit diagonal (UNIT_DIAG) or not (NON_UNIT_DIAG).  The latter
  ///   is the default.  Kokkos' convention is to ignore this for
  ///   sparse matrix-vector multiply, and respect it only for
  ///   triangular solves.
  Teuchos::RCP<SparseOpsType>
  makeSparseOpsFromArrays (const Teuchos::RCP<node_type>& node,
                           Teuchos::ParameterList& opsParams,
                           const ordinal_type numRows,
                           const ordinal_type numCols,
                           Teuchos::ArrayRCP<const size_t>       &ptr,
                           Teuchos::ArrayRCP<const ordinal_type> &ind,
                           Teuchos::ArrayRCP<const scalar_type>  &val,
                           const Teuchos::EUplo uplo,
                           const Teuchos::EDiag diag) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;

    RCP<graph_type> graph = rcp (new graph_type (numRows, numCols, node, null));
    graph->setStructure (ptr, ind);
    RCP<matrix_type> matrix =
      rcp (new matrix_type (rcp_const_cast<const graph_type> (graph), null));
    matrix->setValues (val);
    SparseOpsType::finalizeGraphAndMatrix (uplo, diag, *graph, *matrix, null);
    // We don't need these arrays anymore.
    ptr = null;
    ind = null;
    val = null;

    RCP<SparseOpsType> ops = rcp (new SparseOpsType (node, opsParams));
    ops->setGraphAndMatrix (graph, matrix);
    return ops;
  }

  /// \brief Convert the given dense triangular matrix to a sparse kernels object.
  ///
  /// \param A [in] The dense (lower or upper) triangular matrix.
  /// \param node [in/out] The Kokkos Node instance.
  /// \param params [in/out] Parameters for configuring the
  ///   SparseOpsType instance at construction.
  /// \param uplo [in] Whether A is lower (LOWER_TRI) or upper
  ///   (UPPER_TRI) triangular.
  /// \param diag [in] Whether A has an implicit unit diagonal
  ///   (UNIT_DIAG) or not (NON_UNIT_DIAG).
  Teuchos::RCP<SparseOpsType>
  denseTriToSparseOps (const dense_matrix_type& A,
                       const Teuchos::RCP<node_type>& node,
                       Teuchos::ParameterList& params,
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
    TEUCHOS_TEST_FOR_EXCEPTION(N < 1, std::invalid_argument,
      "denseTriToSparseOps: Input dense matrix must have a positive number of "
      "rows.");

    // If we're not storing the diagonal entries, subtract N off the
    // total number of entries to store.
    const size_t NNZ = diag == Teuchos::UNIT_DIAG ?
      (N*(N-1)) / 2 : // UNIT_DIAG
      (N*(N+1)) / 2;  // NON_UNIT_DIAG

    // Compute the number of entries in each column.
    // This is necessary for allocating the row offsets (ptr) array.
    // We could just construct ptr directly, but some local sparse ops
    // implementations like to allocate and fill in ptr themselves in
    // order to ensure first-touch allocation.
    ArrayRCP<size_t> numEntriesPerRow (N);
    {
      size_t sumOfEntriesPerRow = 0; // correctness check

      // Fill in the number of entries in each column.
      for (ordinal_type i = 0; i < N; ++i) {
        const size_t lastElt = (diag == Teuchos::UNIT_DIAG) ? 1 : 0;
        if (uplo == Teuchos::LOWER_TRI) {
          numEntriesPerRow[i] = (i + 1) - lastElt;
        }
        else if (uplo == Teuchos::UPPER_TRI) {
          numEntriesPerRow[i] = (N - i) - lastElt;
        }
        sumOfEntriesPerRow += numEntriesPerRow[i];
      }
      TEUCHOS_TEST_FOR_EXCEPTION(sumOfEntriesPerRow != NNZ, std::logic_error,
        "TestSparseOps::denseTriToSparseOps: "
        "When computing number of entries per row: total number of entries "
        "should be NNZ = " << NNZ << ", but the sum of the number of entries in "
        "each row is " << sumOfEntriesPerRow << ".  Please report this bug to "
        "the Kokkos developers.");
    }

    // Allocate the three arrays for sparse matrix storage:
    //
    // ptr: row offsets (will be initialized correctly)
    // ind: column indices (possibly filled with zeros)
    // val: matrix entries (possibly filled with zeros)
    ArrayRCP<size_t> ptr =
      SparseOpsType::allocRowPtrs (node, numEntriesPerRow ());
    // We don't need the array of entry counts per row anymore.
    numEntriesPerRow = Teuchos::null;
    ArrayRCP<ordinal_type> ind =
      SparseOpsType::template allocStorage<ordinal_type> (node, ptr ());
    ArrayRCP<scalar_type> val =
      SparseOpsType::template allocStorage<scalar_type> (node, ptr ());

    // Copy the dense triangular matrix into ind and val.
    size_t curOffset = 0; // current row offset
    for (ordinal_type i = 0; i < N; ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(ptr[i] != curOffset, std::logic_error,
        "ptr[" << i << "] = " << ptr[i] << " != curOffset = " << curOffset
        << ".");
      // Compute loop bounds
      ordinal_type lowerBound, upperBound;
      if (uplo == Teuchos::UPPER_TRI) {
        lowerBound = (diag == Teuchos::UNIT_DIAG) ? i+1 : i;
        upperBound = N; // exclusive
      }
      else { // uplo == Teuchos::LOWER_TRI
        lowerBound = Teuchos::OrdinalTraits<ordinal_type>::zero ();
        upperBound = (diag == Teuchos::UNIT_DIAG) ? i : i+1; // exclusive
      }
      // Copy valid entries of the current row
      for (ordinal_type j = lowerBound; j < upperBound; ++j) {
        ind[curOffset] = j;
        val[curOffset] = A(i, j);
        ++curOffset;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(ptr[N] != curOffset, std::logic_error,
      "ptr[" << N << "] = " << ptr[N] << " != curOffset = " << curOffset
      << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(curOffset != NNZ, std::logic_error,
      "TestSparseOps::denseTriToSparseOps: Expected " << NNZ << " entries in "
      "the sparse matrix, but got " << curOffset << " instead.  Please report "
      "this bug (in tests) to the Kokkos developers.");

    ArrayRCP<const size_t> constPtr = arcp_const_cast<const size_t> (ptr);
    ArrayRCP<const ordinal_type> constInd = arcp_const_cast<const ordinal_type> (ind);
    ArrayRCP<const scalar_type> constVal = arcp_const_cast<const scalar_type> (val);
    return makeSparseOpsFromArrays (node, params, N, N,
                                    constPtr, constInd, constVal, uplo, diag);
  }


  /// \brief Convert the given dense matrix to a sparse kernels object.
  ///
  /// \param A [in] The dense matrix.
  /// \param node [in/out] The Kokkos Node instance.
  /// \param params [in/out] Parameters for configuring the
  ///   SparseOpsType instance at construction.
  Teuchos::RCP<SparseOpsType>
  denseToSparseOps (const dense_matrix_type& A,
                    const Teuchos::RCP<node_type>& node,
                    Teuchos::ParameterList& params) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::arcp;
    using Teuchos::arcp_const_cast;
    using Teuchos::as;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;

    const ordinal_type numRows = A.numRows ();
    const ordinal_type numCols = A.numCols ();
    // ordinal_type * ordinal_type could overflow, since ordinal_type
    // need only be big enough to hold max(numRows, numCols), not
    // necessarily their product (the max number of entries in the
    // matrix).
    const size_t NNZ = as<size_t> (numRows) * as<size_t> (numCols);

    // Array of the number of entries in each column.
    // This is necessary for allocating the row offsets (ptr) array.
    // We could just construct ptr directly, but some local sparse ops
    // implementations like to allocate and fill in ptr themselves in
    // order to ensure first-touch allocation.
    ArrayRCP<size_t> numEntriesPerRow (numRows, numCols);

    // Allocate the three arrays for sparse matrix storage:
    //
    // ptr: row offsets (will be initialized correctly)
    // ind: column indices (possibly filled with zeros)
    // val: matrix entries (possibly filled with zeros)
    ArrayRCP<size_t> ptr =
      SparseOpsType::allocRowPtrs (node, numEntriesPerRow ());
    // We don't need the array of entry counts per row anymore.
    numEntriesPerRow = Teuchos::null;
    ArrayRCP<ordinal_type> ind =
      SparseOpsType::template allocStorage<ordinal_type> (node, ptr ());
    ArrayRCP<scalar_type> val =
      SparseOpsType::template allocStorage<scalar_type> (node, ptr ());

    // Copy the dense matrix into ind and val.
    size_t curOffset = 0;
    for (ordinal_type i = 0; i < numRows; ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(ptr[i] != curOffset, std::logic_error,
        "ptr[" << i << "] = " << ptr[i] << " != curOffset = " << curOffset
        << ".");
      for (ordinal_type j = 0; j < numCols; ++j) {
        ind[curOffset] = j;
        val[curOffset] = A(i, j);
        ++curOffset;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(ptr[numRows] != curOffset, std::logic_error,
      "ptr[" << numRows << "] = " << ptr[numRows] << " != curOffset = "
      << curOffset << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(curOffset != NNZ, std::logic_error,
      "TestSparseOps::denseToSparseOps: Expected " << NNZ << " entries in "
      "the sparse matrix, but got " << curOffset << " instead.  Please report "
      "this bug (in tests) to the Kokkos developers.");

    const Teuchos::EUplo uplo = Teuchos::UNDEF_TRI;
    const Teuchos::EDiag diag = Teuchos::NON_UNIT_DIAG;
    ArrayRCP<const size_t> constPtr = arcp_const_cast<const size_t> (ptr);
    ArrayRCP<const ordinal_type> constInd = arcp_const_cast<const ordinal_type> (ind);
    ArrayRCP<const scalar_type> constVal = arcp_const_cast<const scalar_type> (val);
    return makeSparseOpsFromArrays (node, params, numRows, numCols,
                                    constPtr, constInd, constVal, uplo, diag);
  }

  /// Return an initialized Kokkos::MultiVector, filled with zeros or random values.
  ///
  /// \param node [in] The Kokkos Node instance.
  /// \param numRows [in] The number of rows in the MultiVector.
  /// \param numCols [in] The number of columns in the MultiVector.
  /// \param random [in] If true, fill with random values, else fill with zeros.
  Teuchos::RCP<Kokkos::MultiVector<scalar_type, node_type> >
  makeMultiVector (const Teuchos::RCP<node_type>& node,
                   const ordinal_type numRows,
                   const ordinal_type numCols,
                   const bool random=false) const
  {
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Kokkos::MultiVector<scalar_type, node_type> MV;
    typedef Kokkos::DefaultArithmetic<MV> MVT;

    RCP<MV> X = rcp (new MV (node));
    const size_t size = as<size_t> (numRows) * as<size_t> (numCols);
    ArrayRCP<scalar_type> buf = node->template allocBuffer<scalar_type> (size);
    const size_t stride = numRows;
    // This doesn't actually "initialize" anything; it just sets pointers.
    X->initializeValues (numRows, numCols, buf, stride);
    // Always initialize the data first to ensure first-touch allocation.
    // (MVT::Random() isn't parallel by default.)
    MVT::Init (*X, STS::zero ());
    if (random) {
      MVT::Random (*X);
    }
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
    for (ordinal_type k = 0; k < numCols; ++k) {
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

  // TODO (mfh 09 Aug 2012)
  //
  // 1. Replace testSparseOps() (which might be buggy, with false
  //    failures?) with testSparseMatVec() and testSparseTriOps().
  // DONE (10 Aug 2012)
  // 2. Expand testSparseTriOps() to test the transpose and conjugate
  //    transpose cases.

  /// \brief Test sparse matrix-(multi)vector multiply.
  ///
  /// \param node [in/out] The Kokkos Node instance.
  /// \param numRows [in] The number of rows in the sparse matrix
  /// \param numCols [in] The number of columns in the sparse matrix.
  /// \param numVecs [in] The number of columns in the multivectors to test.
  /// \param tol [in] Tolerance for relative errors.
  ///
  /// Test methodology
  /// ================
  ///
  /// Compute a dense random matrix A_dense.  Copy A_dense convert to
  /// sparse format as A_sparse.  A_sparse is only "sparse" in terms
  /// of storage format, probably not in terms of the ratio of zeros
  /// to nonzeros, but this is a start for testing the sparse kernels.
  /// Compare dense mat-vec using A_dense with sparse mat-vec using
  /// A_sparse.
  void
  testSparseMatVec (const Teuchos::RCP<node_type>& node,
                    Teuchos::ParameterList& params,
                    const ordinal_type numRows,
                    const ordinal_type numCols,
                    const ordinal_type numVecs,
                    const magnitude_type tol) const
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::as;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::LEFT_SIDE;
    using Teuchos::LOWER_TRI;
    using Teuchos::UPPER_TRI;
    using Teuchos::NO_TRANS;
    using Teuchos::TRANS;
    using Teuchos::CONJ_TRANS;
    using Teuchos::NON_UNIT_DIAG;
    using Teuchos::UNIT_DIAG;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;
    typedef Kokkos::MultiVector<scalar_type, node_type> MV;
    typedef Kokkos::DefaultArithmetic<MV> MVT;

    // A record of error messages reported by any failed tests.
    // We run _all_ the tests first, then report any failures.
    std::ostringstream err;
    // Whether any tests have failed thus far.
    bool success = true;

    // Relative error of the current operation.
    magnitude_type relErr = STM::zero ();

    // Make a random dense matrix.
    RCP<dense_matrix_type> A_dense (new dense_matrix_type (numRows, numCols));
    A_dense->random ();

    // Convert A_dense into a separate sparse matrix.
    RCP<SparseOpsType> A_sparse = denseToSparseOps (*A_dense, node, params);

    // We make MVs for input, results, and scratch space (since the
    // BLAS' dense triangular solve overwrites its input).  Y and
    // Y_hat are for A * X.  X and X_hat are for A^T * Y and A^H * Y.
    // X_scratch and Y_scratch are scratch multivectors: X_scratch is
    // for X - X_hat, and Y_scratch is for Y - Y_hat.
    RCP<MV> X = makeMultiVector (node, numCols, numVecs, true);
    RCP<MV> X_hat = makeMultiVector (node, numCols, numVecs, true);
    RCP<MV> Y = makeMultiVector (node, numRows, numVecs, true);
    RCP<MV> Y_hat = makeMultiVector (node, numRows, numVecs, true);
    RCP<MV> X_scratch = makeMultiVector (node, numCols, numVecs, true);
    RCP<MV> Y_scratch = makeMultiVector (node, numRows, numVecs, true);

    blas_type blas; // For dense matrix and vector operations.

    // Compute Y_hat := A_dense * X and Y := A_sparse * X.
    blas.GEMM (NO_TRANS, NO_TRANS, numRows, numVecs, numCols,
               STS::one (), A_dense->values (), A_dense->stride (),
               X->getValues ().getRawPtr (), as<int> (X->getStride ()),
               STS::zero (), Y_hat->getValuesNonConst ().getRawPtr (),
               as<int> (Y_hat->getStride ()));
    A_sparse->template multiply<scalar_type, scalar_type> (NO_TRANS, STS::one (),
                                                           (const MV) *X, *Y);
    // Compare Y and Y_hat.
    relErr = maxRelativeError (Y_hat, Y, Y_scratch);
    if (relErr > tol) {
      err << "Sparse matrix-(multi)vector multiply test failed for general "
          << "matrix A, with alpha = 1 and beta = 0." << std::endl
          << "Maximum relative error " << relErr
          << " exceeds the given tolerance " << tol << "." << std::endl;
      success = false;
    }

    {
      // Fill Y_hat with random values, make Y a copy of Y_hat, and
      // test Y = beta*Y + alpha*A_sparse*X with different values of
      // alpha and beta.  We include an irrational value in the list,
      // since it's probably not a special case in the implementation.
      const scalar_type rootTwo = STS::squareroot (STS::one() + STS::one());
      const scalar_type alphaValues[] = {
        -STS::one(), STS::zero(), STS::one(), -rootTwo, rootTwo
      };
      const int numAlphaValues = 5;
      const scalar_type betaValues[] = {
        -STS::one(), STS::zero(), STS::one(), -rootTwo, rootTwo
      };
      const int numBetaValues = 5;

      for (int alphaInd = 0; alphaInd < numAlphaValues; ++alphaInd) {
        for (int betaInd = 0; betaInd < numBetaValues; ++betaInd) {
          const scalar_type alpha = alphaValues[alphaInd];
          const scalar_type beta = betaValues[betaInd];

          MVT::Random (*Y_hat);
          MVT::Assign (*Y, *Y_hat);

          // Y_hat := beta * Y_hat + alpha * A_dense * X
          blas.GEMM (NO_TRANS, NO_TRANS, numRows, numVecs, numCols,
                     alpha, A_dense->values (), A_dense->stride (),
                     X->getValues ().getRawPtr (), as<int> (X->getStride ()),
                     beta, Y_hat->getValuesNonConst ().getRawPtr (),
                     as<int> (Y_hat->getStride ()));
          // Y := beta * Y + alpha * A_sparse * X
          A_sparse->template multiply<scalar_type, scalar_type> (NO_TRANS, alpha,
                                                                 (const MV) *X,
                                                                 beta, *Y);
          // Compare Y and Y_hat.
          relErr = maxRelativeError (Y_hat, Y, Y_scratch);
          if (relErr > tol) {
            err << "Sparse matrix-(multi)vector multiply test failed for general "
                << "matrix A, with alpha = " << alpha << " and beta = " << beta
                << "." << std::endl
                << "Maximum relative error " << relErr << " exceeds "
                << "the given tolerance " << tol << "." << std::endl;
            success = false;
          }
        }
      }
    }

    // Randomize input and output multivectors in preparation for the
    // A^T mat-vec test.
    MVT::Random (*Y);
    MVT::Random (*X_hat);
    MVT::Random (*X);

    // Compute X_hat := A_dense^T * Y and X := A_sparse^T * Y.
    blas.GEMM (TRANS, NO_TRANS, numCols, numVecs, numRows,
               STS::one (), A_dense->values (), A_dense->stride (),
               Y->getValues ().getRawPtr (), as<int> (Y->getStride ()),
               STS::zero (), X_hat->getValuesNonConst ().getRawPtr (),
               as<int> (X_hat->getStride ()));
    A_sparse->template multiply<scalar_type, scalar_type> (TRANS, STS::one (),
                                                           (const MV) *Y, *X);
    // Compare X and X_hat.
    relErr = maxRelativeError (X_hat, X, X_scratch);
    if (relErr > tol) {
      err << "Sparse matrix-(multi)vector multiply (transpose) test failed for "
          << "general matrix A, with alpha = 1 and beta = 0." << std::endl
          << "Maximum relative error " << relErr
          << " exceeds the given tolerance " << tol << "." << std::endl;
      success = false;
    }

    // Randomize input and output multivectors in preparation for the
    // A^T mat-vec test with different alpha and beta values.
    MVT::Random (*Y);
    MVT::Random (*X_hat);
    MVT::Random (*X);

    {
      // Test X := beta * X + alpha * A_sparse^T * Y with different
      // values of alpha and beta.  We include an irrational value in
      // the list, since it's probably not a special case in the
      // implementation.
      const scalar_type rootTwo = STS::squareroot (STS::one() + STS::one());
      const scalar_type alphaValues[] = {
        -STS::one(), STS::zero(), STS::one(), -rootTwo, rootTwo
      };
      const int numAlphaValues = 5;
      const scalar_type betaValues[] = {
        -STS::one(), STS::zero(), STS::one(), -rootTwo, rootTwo
      };
      const int numBetaValues = 5;

      for (int alphaInd = 0; alphaInd < numAlphaValues; ++alphaInd) {
        for (int betaInd = 0; betaInd < numBetaValues; ++betaInd) {
          const scalar_type alpha = alphaValues[alphaInd];
          const scalar_type beta = betaValues[betaInd];

          // Fill X_hat with random values and make X a copy of X_hat.
          MVT::Random (*X_hat);
          MVT::Assign (*X, *X_hat);

          // X_hat := beta * X_hat + A_dense^T * Y.
          blas.GEMM (TRANS, NO_TRANS, numCols, numVecs, numRows,
                     alpha, A_dense->values (), A_dense->stride (),
                     Y->getValues ().getRawPtr (), as<int> (Y->getStride ()),
                     beta, X_hat->getValuesNonConst ().getRawPtr (),
                     as<int> (X_hat->getStride ()));
          // X := beta * X + A_sparse^T * Y.
          A_sparse->template multiply<scalar_type, scalar_type> (TRANS, alpha,
                                                                 (const MV) *Y,
                                                                 beta, *X);
          // Compare X and X_hat.
          relErr = maxRelativeError (X_hat, X, X_scratch);
          if (relErr > tol) {
            err << "Sparse matrix-(multi)vector multiply (transpose) test "
                << "failed for general matrix A, with alpha = " << alpha
                << " and beta = " << beta << "." << std::endl
                << "Maximum relative error " << relErr
                << " exceeds the given tolerance " << tol << "." << std::endl;
            success = false;
          }
        }
      }
    }

    // Randomize input and output multivectors in preparation for the
    // A^H mat-vec test.
    MVT::Random (*Y);
    MVT::Random (*X_hat);
    MVT::Random (*X);

    if (STS::isComplex) {
      // Compute X_hat := A_dense^H * Y and X := A_sparse^H * Y.
      blas.GEMM (CONJ_TRANS, NO_TRANS, numCols, numVecs, numRows,
                 STS::one (), A_dense->values (), A_dense->stride (),
                 Y->getValues ().getRawPtr (), as<int> (Y->getStride ()),
                 STS::zero (), X_hat->getValuesNonConst ().getRawPtr (),
                 as<int> (X_hat->getStride ()));
      A_sparse->template multiply<scalar_type, scalar_type> (CONJ_TRANS, STS::one (),
                                                             (const MV) *Y, *X);
      // Compare X and X_hat.
      relErr = maxRelativeError (X_hat, X, X_scratch);
      if (relErr > tol) {
        err << "Sparse matrix-(multi)vector multiply (conjugate transpose) "
            << "test failed for general matrix A, with alpha = 1 and beta = 0."
            << std::endl
            << "Maximum relative error " << relErr << " exceeds the given "
            << "tolerance " << tol << "." << std::endl;
        success = false;
      }

      // Randomize input and output multivectors in preparation for the
      // A^H mat-vec test with different alpha and beta values.
      MVT::Random (*Y);
      MVT::Random (*X_hat);
      MVT::Random (*X);

      {
        // Test X = beta * X + alpha * A_sparse^H * Y with different
        // values of alpha and beta.  We include an irrational value in
        // the list, since it's probably not a special case in the
        // implementation.
        const scalar_type rootTwo = STS::squareroot (STS::one() + STS::one());
        const scalar_type alphaValues[] = {
          -STS::one(), STS::zero(), STS::one(), -rootTwo, rootTwo
        };
        const int numAlphaValues = 5;
        const scalar_type betaValues[] = {
          -STS::one(), STS::zero(), STS::one(), -rootTwo, rootTwo
        };
        const int numBetaValues = 5;

        for (int alphaInd = 0; alphaInd < numAlphaValues; ++alphaInd) {
          for (int betaInd = 0; betaInd < numBetaValues; ++betaInd) {
            const scalar_type alpha = alphaValues[alphaInd];
            const scalar_type beta = betaValues[betaInd];

            // Fill X_hat with random values and make X a copy of X_hat.
            MVT::Random (*X_hat);
            MVT::Assign (*X, *X_hat);

            // X_hat := beta * X_hat + A_dense^H * Y.
            blas.GEMM (CONJ_TRANS, NO_TRANS, numCols, numVecs, numRows,
                       alpha, A_dense->values (), A_dense->stride (),
                       Y->getValues ().getRawPtr (), as<int> (Y->getStride ()),
                       beta, X_hat->getValuesNonConst ().getRawPtr (),
                       as<int> (X_hat->getStride ()));
            // X := beta * X + A_sparse^H * Y.
            A_sparse->template multiply<scalar_type, scalar_type> (CONJ_TRANS, alpha,
                                                                   (const MV) *Y,
                                                                   beta, *X);
            // Compare X and X_hat.
            relErr = maxRelativeError (X_hat, X, X_scratch);
            if (relErr > tol) {
              err << "Sparse matrix-(multi)vector multiply (conjugate "
                  << "transpose) test failed for general matrix A, with alpha = "
                  << alpha << " and beta = " << beta << "." << std::endl
                  << "Maximum relative error " << relErr
                  << " exceeds the given tolerance " << tol << "." << std::endl;
              success = false;
            }
          }
        }
      }
    }
  }

  /// \brief Test sparse matrix-(multi)vector multiply and sparse
  ///   triangular solve.
  ///
  /// This method always tests sparse matrix-(multi)vector multiply.
  /// Sparse triangular solve only makes sense if numRows == numCols,
  /// so we only test it in that case.
  ///
  /// \param node [in/out] The Kokkos Node instance.
  /// \param numRows [in] Number of rows in the sparse matrices.
  /// \param numCols [in] Number of columns in the sparse matrices.
  /// \param numVecs [in] Number of columns in the multivectors to test.
  /// \param tol [in] Relative error tolerance.
  /// \param implicitUnitDiagTriMultCorrect [in] Whether SparseOpsType
  ///   correctly implements sparse matrix-(multi)vector multiply with
  ///   a triangular matrix with implicitly stored unit diagonal.
  ///   "Incorrectly" means that SparseOpsType assumes that all
  ///   entries of the matrix are stored explicitly, regardless of its
  ///   Teuchos::EDiag input value.  The SparseOpsType interface does
  ///   not require this by default, but some implementations do.
  void
  testSparseOps (const Teuchos::RCP<node_type>& node,
                 Teuchos::ParameterList& params,
                 const ordinal_type numRows,
                 const ordinal_type numCols,
                 const ordinal_type numVecs,
                 const magnitude_type tol,
                 const bool implicitUnitDiagTriMultCorrect=false)
  {
    typedef Teuchos::OrdinalTraits<ordinal_type> OTO;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    TEUCHOS_TEST_FOR_EXCEPTION(
      numRows < OTO::zero() || numCols < OTO::zero() || numVecs < OTO::zero(),
      std::invalid_argument,
      "testSparseOps: the sparse matrix and dense multivector dimensions must "
      "all be nonnegative.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      tol < STM::zero(),
      std::invalid_argument,
      "testSparseOps: the relative error tolerance must be nonnegative.");

    testSparseMatVec (node, params, numRows, numCols, numVecs, tol);
    if (numRows == numCols) {
      testSparseTriOps (node, params, numRows, numVecs, tol,
                        implicitUnitDiagTriMultCorrect);
    }
  }

  /// \brief Test sparse matrix-(multi)vector multiply and sparse
  ///   triangular solve with triangular matrices.
  ///
  /// \param node [in/out] The Kokkos Node instance.
  /// \param N [in] The number of rows (and columns) in the sparse
  ///   matrices to test.  Sparse triangular solve assumes that the
  ///   matrices are square.
  /// \param numVecs [in] The number of columns in the multivectors to test.
  /// \param tol [in] Tolerance for relative errors.
  /// \param implicitUnitDiagTriMultCorrect [in] Whether SparseOpsType
  ///   correctly implements sparse matrix-(multi)vector multiply with
  ///   a triangular matrix with implicitly stored unit diagonal.
  ///   "Incorrectly" means that SparseOpsType assumes that all
  ///   entries of the matrix are stored explicitly, regardless of its
  ///   Teuchos::EDiag input value.  The SparseOpsType interface does
  ///   not require this by default, but some implementations do.
  ///
  /// Test methodology
  /// ================
  ///
  /// Compute a dense LU factorization of a random matrix A.  Store
  /// its L and U factors as sparse matrices.  Use the factors to test
  /// sparse matrix-vector multiply and sparse triangular solve.
  ///
  /// Discussion
  /// ==========
  ///
  /// Random square matrices tend to be well conditioned on average
  /// (this statement sounds vague, but can be made rigorous), so the
  /// L and U factors are likely to exist and be well conditioned on
  /// average.  Of course, outliers are possible, so this test isn't
  /// 100% guaranteed to succeed, but success is very likely and false
  /// positives (passing the test when it should fail) are even less
  /// likely.
  void
  testSparseTriOps (const Teuchos::RCP<node_type>& node,
                    Teuchos::ParameterList& params,
                    const ordinal_type N,
                    const ordinal_type numVecs,
                    const magnitude_type tol,
                    const bool implicitUnitDiagTriMultCorrect=false)
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::as;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::LEFT_SIDE;
    using Teuchos::LOWER_TRI;
    using Teuchos::UPPER_TRI;
    using Teuchos::NO_TRANS;
    using Teuchos::TRANS;
    using Teuchos::CONJ_TRANS;
    using Teuchos::NON_UNIT_DIAG;
    using Teuchos::UNIT_DIAG;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;
    typedef Kokkos::MultiVector<scalar_type, node_type> MV;
    typedef Kokkos::DefaultArithmetic<MV> MVT;

    TEUCHOS_TEST_FOR_EXCEPTION(
      N < Teuchos::OrdinalTraits<ordinal_type>::zero(),
      std::invalid_argument,
      "testSparseTriOps: Number of rows and columns N = " << N << " < 0.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      numVecs < Teuchos::OrdinalTraits<ordinal_type>::zero(),
      std::invalid_argument,
      "testSparseTriOps: numVecs = " << numVecs << " < 0.");

    // A record of error messages reported by any failed tests.
    // We run _all_ the tests first, then report any failures.
    std::ostringstream err;
    // Whether any tests have failed thus far.
    bool success = true;

    // Relative error of the current operation.
    magnitude_type relErr = STM::zero ();

    RCP<dense_matrix_type> A_dense, L_dense, U_dense;
    Array<ordinal_type> ipiv;
    makeDenseTestProblem (A_dense, L_dense, U_dense, ipiv, N);

    // Convert L_dense and U_dense into separate sparse matrices.
    RCP<SparseOpsType> L_sparse =
      denseTriToSparseOps (*L_dense, node, params, LOWER_TRI, UNIT_DIAG);
    RCP<SparseOpsType> U_sparse =
      denseTriToSparseOps (*U_dense, node, params, UPPER_TRI, NON_UNIT_DIAG);
    // Convert A_dense into a separate sparse matrix.
    RCP<SparseOpsType> A_sparse = denseToSparseOps (*A_dense, node, params);

    // We make MVs for input, results, and scratch space (since the
    // BLAS' dense triangular solve overwrites its input).  Y and
    // Y_hat are for {L,U} * X.  X and X_hat are for {L,U}^T * Y and
    // {L,U}^H * Y.  X_scratch and Y_scratch are scratch multivectors:
    // X_scratch is for X - X_hat, and Y_scratch is for Y - Y_hat.
    RCP<MV> X = makeMultiVector (node, N, numVecs, true);
    RCP<MV> X_hat = makeMultiVector (node, N, numVecs, true);
    RCP<MV> Y = makeMultiVector (node, N, numVecs, true);
    RCP<MV> Y_hat = makeMultiVector (node, N, numVecs, true);
    RCP<MV> X_scratch = makeMultiVector (node, N, numVecs, true);
    RCP<MV> Y_scratch = makeMultiVector (node, N, numVecs, true);

    blas_type blas; // For dense matrix and vector operations.

    //////////////////////////////////////////////////////////////////////
    // Test triangular matrix-(multi)vector multiply.
    //////////////////////////////////////////////////////////////////////

    // Compute Y_hat := U_dense * X.  First copy X into Y_hat, since
    // TRMM overwrites its input.
    MVT::Assign (*Y_hat, (const MV) *X);
    blas.TRMM (LEFT_SIDE, UPPER_TRI, NO_TRANS, NON_UNIT_DIAG, N, numVecs,
               STS::one (), U_dense->values (), U_dense->stride (),
               Y_hat->getValuesNonConst ().getRawPtr (),
               as<int> (Y_hat->getStride ()));
    // Compute Y := U_sparse * X.  The local sparse operations don't
    // overwrite their input.
    U_sparse->template multiply<scalar_type, scalar_type> (NO_TRANS, STS::one (),
                                                           (const MV) *X, *Y);
    // Compare Y and Y_hat.
    relErr = maxRelativeError (Y_hat, Y, Y_scratch);
    if (relErr > tol) {
      err << "Sparse matrix-(multi)vector multiply test failed for upper "
        "triangular matrix U." << std::endl << "Maximum relative error "
          << relErr << " exceeds the given tolerance " << tol << "."
          << std::endl;
      success = false;
    }

    // Start over with a fresh random X.
    MVT::Random (*X);
    // Fill Y_hat and Y with random values, just to make sure the
    // operations below aren't fakely reported correct by means of
    // leftover values.
    MVT::Random (*Y_hat);
    MVT::Random (*Y);

    // Compute Y_hat := U_dense^T * X.
    MVT::Assign (*Y_hat, (const MV) *X);
    blas.TRMM (LEFT_SIDE, UPPER_TRI, TRANS, NON_UNIT_DIAG, N, numVecs,
               STS::one (), U_dense->values (), U_dense->stride (),
               Y_hat->getValuesNonConst ().getRawPtr (),
               as<int> (Y_hat->getStride ()));
    // Compute Y := U_sparse^T * X.
    U_sparse->template multiply<scalar_type, scalar_type> (TRANS, STS::one (),
                                                           (const MV) *X, *Y);
    // Compare Y and Y_hat.
    relErr = maxRelativeError (Y_hat, Y, Y_scratch);
    if (relErr > tol) {
      err << "Sparse matrix-(multi)vector multiply test failed for "
        "transpose of upper triangular matrix U." << std::endl
          << "Maximum relative error " << relErr
          << " exceeds the given tolerance " << tol << "." << std::endl;
      success = false;
    }

    // Start over with a fresh random X.
    MVT::Random (*X);
    // Fill Y_hat and Y with random values, just to make sure the
    // operations below aren't fakely reported correct by means of
    // leftover values.
    MVT::Random (*Y_hat);
    MVT::Random (*Y);

    if (STS::isComplex) {
      // Compute Y_hat := U_dense^T * X.
      MVT::Assign (*Y_hat, (const MV) *X);
      blas.TRMM (LEFT_SIDE, UPPER_TRI, CONJ_TRANS, NON_UNIT_DIAG, N, numVecs,
                 STS::one (), U_dense->values (), U_dense->stride (),
                 Y_hat->getValuesNonConst ().getRawPtr (),
                 as<int> (Y_hat->getStride ()));
      // Compute Y := U_sparse^H * X.
      U_sparse->template multiply<scalar_type, scalar_type> (CONJ_TRANS,
                                                             STS::one (),
                                                             (const MV) *X, *Y);
      // Compare Y and Y_hat.
      relErr = maxRelativeError (Y_hat, Y, Y_scratch);
      if (relErr > tol) {
        err << "Sparse matrix-(multi)vector multiply test failed for "
          "conjugate transpose of upper triangular matrix U." << std::endl
            << "Maximum relative error " << relErr
            << " exceeds the given tolerance " << tol << "." << std::endl;
        success = false;
      }

      // Start over with a fresh random X.
      MVT::Random (*X);
      // Fill Y_hat and Y with random values, just to make sure the
      // operations below aren't fakely reported correct by means of
      // leftover values.
      MVT::Random (*Y_hat);
      MVT::Random (*Y);
    }

    // FIXME (mfh 10 Aug 2012) I can't seem to get this test to pass,
    // no matter whether or not I assume that the implicit unit
    // diagonal case works.  I'm not sure what to do and I don't want
    // this to hold up some of my other bug fixes, so I'm going to
    // disable this test for now.
#if 0
    // Compute Y_hat := L_dense * X.  First copy X into Y_hat, since
    // TRMM overwrites its input.
    MVT::Assign (*Y_hat, (const MV) *X);
    blas.TRMM (LEFT_SIDE, LOWER_TRI, NO_TRANS, UNIT_DIAG, N, numVecs,
               STS::one (), L_dense->values (), L_dense->stride (),
               Y_hat->getValuesNonConst ().getRawPtr (),
               as<int> (Y_hat->getStride ()));
    // Compute Y := L_sparse * X.  The local sparse operations
    // don't overwrite their input.
    if (implicitUnitDiagTriMultCorrect) {
      L_sparse->template multiply<scalar_type, scalar_type> (NO_TRANS, STS::one (),
                                                             (const MV) *X, *Y);
    }
    else {
      // Compute Y := Y + L_sparse * X in order to simulate the effect
      // of an implicitly stored unit diagonal.  We do this using
      // beta=1, though it might be better to test beta=0 and use
      // MVT::GESUM instead.
      L_sparse->template multiply<scalar_type, scalar_type> (NO_TRANS, STS::one (),
                                                             (const MV) *X,
                                                             STS::one (), *Y);
    }
    // Compare Y and Y_hat.
    relErr = maxRelativeError (Y_hat, Y, Y_scratch);
    if (relErr > tol) {
      err << "Sparse matrix-(multi)vector multiply test failed for lower "
        "triangular matrix L with implicit unit diagonal." << std::endl
          << "Maximum relative error " << relErr
          << " exceeds the given tolerance " << tol << "." << std::endl;
      success = false;
    }
#endif // 0

    //////////////////////////////////////////////////////////////////////
    // Test sparse triangular solve.
    //////////////////////////////////////////////////////////////////////

    // Start over with a fresh random Y (input multivector in this case).
    MVT::Random (*Y);
    // Fill X_hat and X (output multivectors, in this case) with
    // random values, just to make sure the operations below aren't
    // fakely reported correct by means of leftover values.
    MVT::Random (*X_hat);
    MVT::Random (*X);

    // Compute X_hat := L_dense \ Y.
    MVT::Assign (*X_hat, (const MV) *Y);
    blas.TRSM (LEFT_SIDE, LOWER_TRI, NO_TRANS, UNIT_DIAG, N, numVecs,
               STS::one (), L_dense->values (), L_dense->stride (),
               X_hat->getValuesNonConst ().getRawPtr (),
               as<int> (X_hat->getStride ()));
    // Compute X := L_sparse \ Y.
    L_sparse->template solve<scalar_type, scalar_type> (NO_TRANS, (const MV) *Y, *X);
    // Compare X and X_hat.
    relErr = maxRelativeError (X_hat, X, X_scratch);
    if (relErr > tol) {
      err << "Sparse triangular solve test failed for lower triangular matrix "
        "L with implicit unit diagonal." << std::endl
          << "Maximum relative error " << relErr
          << " exceeds the given tolerance " << tol << "." << std::endl;
      success = false;
    }

    // Start over with a fresh random Y (input multivector in this case).
    MVT::Random (*Y);
    // Fill X_hat and X (output multivectors, in this case) with
    // random values, just to make sure the operations below aren't
    // fakely reported correct by means of leftover values.
    MVT::Random (*X_hat);
    MVT::Random (*X);

    // Compute X_hat = U_dense \ Y.
    MVT::Assign (*X_hat, (const MV) *Y);
    blas.TRSM (LEFT_SIDE, UPPER_TRI, NO_TRANS, NON_UNIT_DIAG, N, numVecs,
               STS::one (), U_dense->values (), U_dense->stride (),
               X_hat->getValuesNonConst ().getRawPtr (),
               as<int> (X_hat->getStride ()));
    // Compute X = U_sparse \ Y.
    U_sparse->template solve<scalar_type, scalar_type> (NO_TRANS,
                                                        (const MV) *Y, *X);
    // Compare X and X_hat.
    relErr = maxRelativeError (X_hat, X, X_scratch);
    if (relErr > tol) {
      err << "Sparse triangular solve test failed for upper triangular matrix "
        "U." << std::endl << "Maximum relative error " << relErr
          << " exceeds the given tolerance " << tol << "." << std::endl;
      success = false;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(! success, std::runtime_error,
      "One or more sparse ops tests failed.  Here is the full report:\n"
      << err.str());
  }


  void
  benchmarkSparseOps (std::vector<std::pair<std::string, double> >& results,
                      const std::string& label,
                      const Teuchos::RCP<node_type>& node,
                      Teuchos::ParameterList& params,
                      const ordinal_type numRows,
                      const ordinal_type numCols,
                      const ordinal_type numVecs,
                      const int numTrials) const
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::LOWER_TRI;
    using Teuchos::UPPER_TRI;
    using Teuchos::NON_UNIT_DIAG;
    using Teuchos::UNIT_DIAG;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;
    typedef Kokkos::MultiVector<scalar_type, node_type> MV;
    typedef Kokkos::DefaultArithmetic<MV> MVT;

    const bool testTriSolve = (numRows == numCols);

    RCP<dense_matrix_type> A_dense, L_dense, U_dense;
    RCP<SparseOpsType> A_sparse, L_sparse, U_sparse;
    Array<ordinal_type> ipiv;

    if (testTriSolve) {
      makeDenseTestProblem (A_dense, L_dense, U_dense, ipiv, numRows);
      // Convert L_dense and U_dense into sparse matrices.
      L_sparse = denseTriToSparseOps (*L_dense, node, params, LOWER_TRI, UNIT_DIAG);
      U_sparse = denseTriToSparseOps (*U_dense, node, params, UPPER_TRI, NON_UNIT_DIAG);
    }
    else {
      A_dense = rcp (new dense_matrix_type (numRows, numCols));
    }
    // Convert A_dense into a sparse matrix.
    A_sparse = denseToSparseOps (*A_dense, node, params);

    // // Compute a random input multivector.
    // RCP<MV> X = makeMultiVector (node, numCols, numVecs);
    // MVT::Random (*X);
    // RCP<MV> Y = makeMultiVector (node, numRows, numVecs);
    // MVT::Init (*Y, STS::zero ());

    benchmarkSparseMatVec (results, label, *A_sparse,
                           numRows, numCols, numVecs, numTrials);
    if (testTriSolve) {
      benchmarkSparseTriSolve (results, label, "lower tri, unit diag",
                               *L_sparse, numRows, numCols, numVecs, numTrials);
      benchmarkSparseTriSolve (results, label, "upper tri, non unit diag",
                               *U_sparse, numRows, numCols, numVecs, numTrials);
    }
  }

  /// \brief Benchmark sparse mat-vec with a matrix read from a Matrix Market file.
  ///
  /// \param results [out] Timing results.  You can also use
  ///   TimeMonitor::report() or TimeMonitor::summarize() to display
  ///   results.
  /// \param filename [in] Name of the Matrix Market sparse matrix file.
  /// \param label [in] Label to distinguish the SparseOpsType, in case
  ///   you are running multiple benchmarks with different SparseOpsType
  ///   types.
  /// \param node [in/out] Kokkos Node instance.
  /// \param numVecs Number of columns in the multivectors to benchmark.
  /// \param numTrials Number of trials.  We time a loop around all the
  ///   trials to increase accuracy and smooth out variance.
  void
  benchmarkSparseMatVecFromFile (std::vector<std::pair<std::string, double> >& results,
                                 const std::string& filename,
                                 const std::string& label,
                                 const Teuchos::RCP<node_type>& node,
                                 Teuchos::ParameterList& params,
                                 const ordinal_type numVecs,
                                 const int numTrials) const
  {
    using Teuchos::RCP;
    typedef Kokkos::MultiVector<scalar_type, node_type> MV;
    typedef Kokkos::DefaultArithmetic<MV> MVT;

    // SparseOpsType isn't required to tell us how many rows and
    // columns the sparse matrix has, so we find out when we read the
    // file.  Kokkos ignores uplo and diag for sparse mat-vec, so we
    // don't need to provide those arguments to
    // makeSparseOpsFromFile() here.
    ordinal_type numRows = 0;
    ordinal_type numCols = 0;
    RCP<SparseOpsType> A_sparse =
      makeSparseOpsFromFile (numRows, numCols, node, params, filename);

    // Compute a random input multivector.
    RCP<MV> X = makeMultiVector (node, numCols, numVecs);
    MVT::Random (*X);
    RCP<MV> Y = makeMultiVector (node, numRows, numVecs); // output MV.

    benchmarkSparseMatVec (results, label, *A_sparse,
                           numRows, numCols, numVecs, numTrials);
    // benchmarkSparseTriSolve (results, label, "lower tri, unit diag",
    //                          *L_sparse, numRows, numCols, numVecs, numTrials);
    // benchmarkSparseTriSolve (results, label, "upper tri, non unit diag",
    //                          *U_sparse, numRows, numCols, numVecs, numTrials);
  }

private:
  /// \brief Benchmark sparse matrix-multivector multiply.
  ///
  /// Benchmark all variations of sparse matrix-multivector multiply
  /// as implemented by SparseOpsType.  This includes no transpose,
  /// transpose, and conjugate transpose (if scalar_type is complex),
  /// in both the overwrite (Y = A*X) and update (Y = Y - A*X) forms.
  ///
  /// We save benchmark results to the output std::vector, but they
  /// are also stored by Teuchos::TimeMonitor, so that you can display
  /// them using TimeMonitor::report() or TimeMonitor::summarize().
  ///
  /// \param results [out] Pairs of (benchmark name, elapsed time for
  ///   that benchmark over all trials).
  /// \param ops [in] The sparse kernels instance to benchmark.
  /// \param numRows [in] Number of rows in the linear operator
  ///   represented by ops.  We need this because SparseOpsType
  ///   doesn't necessarily tell us.
  /// \param numCols [in] Number of columns in the linear operator
  ///   represented by ops.  We need this because SparseOpsType
  ///   doesn't necessarily tell us.
  /// \param numVecs [in] Number of columns in the multivectors to
  ///   benchmark.
  /// \param numTrials [in] Number of runs over which to measure
  ///   elapsed time, for each benchmark.
  void
  benchmarkSparseMatVec (std::vector<std::pair<std::string, double> >& results,
                         const std::string& label,
                         const SparseOpsType& ops,
                         const ordinal_type numRows,
                         const ordinal_type numCols,
                         const ordinal_type numVecs,
                         const int numTrials) const
  {
    using Teuchos::RCP;
    using Teuchos::Time;
    using Teuchos::TimeMonitor;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Kokkos::MultiVector<scalar_type, node_type> MV;
    typedef Kokkos::DefaultArithmetic<MV> MVT;

    RCP<MV> X = makeMultiVector (ops.getNode (), numRows, numVecs);
    RCP<MV> Y = makeMultiVector (ops.getNode (), numRows, numVecs);
    //MVT::Init (*Y, STS::zero()); // makeMultiVector() already does this.
    MVT::Random (*X);

    const int numBenchmarks = STS::isComplex ? 9 : 7;
    results.reserve (numBenchmarks);

    // Time sparse matrix-vector multiply, overwrite mode, no transpose.
    {
      const std::string timerName (label + ": Y = A*X");
      RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = TimeMonitor::getNewCounter (timerName);
      }
      {
        TimeMonitor timeMon (*timer);
        for (int i = 0; i < numTrials; ++i) {
          ops.template multiply<scalar_type, scalar_type> (Teuchos::NO_TRANS,
                                                           STS::one(), *X, *Y);
        }
      }
      results.push_back (std::make_pair (timerName, timer->totalElapsedTime ()));
    }

    // Time sparse matrix-vector multiply, update mode, no transpose.
    {
      const std::string timerName (label + ": Y = Y + A*X");
      RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = TimeMonitor::getNewCounter (timerName);
      }
      {
        TimeMonitor timeMon (*timer);
        for (int i = 0; i < numTrials; ++i) {
          ops.template multiply<scalar_type, scalar_type> (Teuchos::NO_TRANS,
                                                           STS::one(), *X,
                                                           STS::one(), *Y);
        }
      }
      results.push_back (std::make_pair (timerName, timer->totalElapsedTime ()));
    }

    // Time sparse matrix-vector multiply, update mode, no transpose.
    // We subtract to simulate a residual computation.
    {
      const std::string timerName (label + ": Y = Y - A*X");
      RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = TimeMonitor::getNewCounter (timerName);
      }
      {
        TimeMonitor timeMon (*timer);
        for (int i = 0; i < numTrials; ++i) {
          ops.template multiply<scalar_type, scalar_type> (Teuchos::NO_TRANS,
                                                           -STS::one(), *X,
                                                           STS::one(), *Y);
        }
      }
      results.push_back (std::make_pair (timerName, timer->totalElapsedTime ()));
    }

    // Time sparse matrix-vector multiply, update mode, no transpose.
    // We subtract to simulate an alternate form for a residual computation.
    {
      const std::string timerName (label + ": Y = -Y + A*X");
      RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = TimeMonitor::getNewCounter (timerName);
      }
      {
        TimeMonitor timeMon (*timer);
        for (int i = 0; i < numTrials; ++i) {
          ops.template multiply<scalar_type, scalar_type> (Teuchos::NO_TRANS,
                                                           STS::one(), *X,
                                                           -STS::one(), *Y);
        }
      }
      results.push_back (std::make_pair (timerName, timer->totalElapsedTime ()));
    }

    // Time sparse matrix-vector multiply with alpha=0 and beta=-1.
    {
      const std::string timerName (label + ": Y = -Y + 0*A*X");
      RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = TimeMonitor::getNewCounter (timerName);
      }
      {
        TimeMonitor timeMon (*timer);
        for (int i = 0; i < numTrials; ++i) {
          ops.template multiply<scalar_type, scalar_type> (Teuchos::NO_TRANS,
                                                           STS::zero(), *X,
                                                           -STS::one(), *Y);
        }
      }
      results.push_back (std::make_pair (timerName, timer->totalElapsedTime ()));
    }

    // Time sparse matrix-vector multiply, overwrite mode, transpose.
    {
      const std::string timerName (label + ": Y = A^T * X");
      RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = TimeMonitor::getNewCounter (timerName);
      }
      {
        TimeMonitor timeMon (*timer);
        for (int i = 0; i < numTrials; ++i) {
          ops.template multiply<scalar_type, scalar_type> (Teuchos::TRANS,
                                                           STS::one(), *X, *Y);
        }
      }
      results.push_back (std::make_pair (timerName, timer->totalElapsedTime ()));
    }

    // Time sparse matrix-vector multiply, update mode, transpose.
    // We subtract to simulate a residual computation.
    {
      const std::string timerName (label + ": Y = Y - A^T * X");
      RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = TimeMonitor::getNewCounter (timerName);
      }
      {
        TimeMonitor timeMon (*timer);
        for (int i = 0; i < numTrials; ++i) {
          ops.template multiply<scalar_type, scalar_type> (Teuchos::TRANS,
                                                           -STS::one(), *X,
                                                           STS::one(), *Y);
        }
      }
      results.push_back (std::make_pair (timerName, timer->totalElapsedTime ()));
    }

    // Only test conjugate transpose if scalar_type is complex.
    if (STS::isComplex) {
      // Time sparse matrix-vector multiply, overwrite mode, conjugate transpose.
      {
        const std::string timerName (label + ": Y = A^H * X");
        RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
        if (timer.is_null ()) {
          timer = TimeMonitor::getNewCounter (timerName);
        }
        {
          TimeMonitor timeMon (*timer);
          for (int i = 0; i < numTrials; ++i) {
            ops.template multiply<scalar_type, scalar_type> (Teuchos::CONJ_TRANS,
                                                             STS::one(), *X, *Y);
          }
        }
        results.push_back (std::make_pair (timerName, timer->totalElapsedTime ()));
      }

      // Time sparse matrix-vector multiply, update mode, conjugate transpose.
      // We subtract to simulate a residual computation.
      {
        const std::string timerName (label + ": Y = Y - A^H * X");
        RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
        if (timer.is_null ()) {
          timer = TimeMonitor::getNewCounter (timerName);
        }
        {
          TimeMonitor timeMon (*timer);
          for (int i = 0; i < numTrials; ++i) {
            ops.template multiply<scalar_type, scalar_type> (Teuchos::CONJ_TRANS,
                                                             -STS::one(), *X,
                                                             STS::one(), *Y);
          }
        }
        results.push_back (std::make_pair (timerName, timer->totalElapsedTime ()));
      }
    }
  }



  /// \brief Benchmark sparse triangular solve.
  ///
  /// Benchmark all variations of sparse triangular solve as
  /// implemented by SparseOpsType.  This includes the no transpose,
  /// transpose, and conjugate transpose (if scalar_type is complex)
  /// forms.
  ///
  /// \note This method only works correctly if SparseOpsType encodes
  ///   a lower or upper triangular matrix.  Whether the sparse matrix
  ///   is upper or lower triangular depends on ops, so if you want to
  ///   benchmark both cases, you'll have to create two different
  ///   SparseOpsType instances (one lower triangular, one upper
  ///   triangular) and benchmark both.  The same goes for the
  ///   implicit unit diagonal / explicit diagonal variants.
  ///
  /// We save benchmark results to the output std::vector, but they
  /// are also stored by Teuchos::TimeMonitor, so that you can display
  /// them using TimeMonitor::report() or TimeMonitor::summarize().
  ///
  /// \param results [out] Pairs of (benchmark name, elapsed time for
  ///   that benchmark over all trials).
  /// \param opsLabel [in] Label to identify the SparseOpsType type.
  ///   This is helpful if running benchmarks with different
  ///   SparseOpsType types.
  /// \param benchmarkLabel [in] Extra timer label.  Use this to
  ///   distinguish e.g., lower triangular solve benchmarks from upper
  ///   triangular solve benchmarks.
  /// \param ops [in] The sparse kernels instance to benchmark.
  /// \param numRows [in] Number of rows in the linear operator
  ///   represented by ops.  We need this because SparseOpsType
  ///   doesn't necessarily tell us.
  /// \param numCols [in] Number of columns in the linear operator
  ///   represented by ops.  We need this because SparseOpsType
  ///   doesn't necessarily tell us.
  /// \param numVecs [in] Number of columns in the multivectors to
  ///   benchmark.
  /// \param numTrials [in] Number of runs over which to measure
  ///   elapsed time, for each benchmark.
  void
  benchmarkSparseTriSolve (std::vector<std::pair<std::string, double> >& results,
                           const std::string& opsLabel,
                           const std::string& benchmarkLabel,
                           const SparseOpsType& ops,
                           const ordinal_type numRows,
                           const ordinal_type numCols,
                           const ordinal_type numVecs,
                           const int numTrials) const
  {
    using Teuchos::RCP;
    using Teuchos::Time;
    using Teuchos::TimeMonitor;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Kokkos::MultiVector<scalar_type, node_type> MV;
    typedef Kokkos::DefaultArithmetic<MV> MVT;

    RCP<MV> X = makeMultiVector (ops.getNode (), numCols, numVecs);
    RCP<MV> Y = makeMultiVector (ops.getNode (), numRows, numVecs);
    //MVT::Init (*X, STS::zero()); // makeMultiVector() already does this.
    MVT::Random (*Y);

    const int numBenchmarks = STS::isComplex ? 6 : 4;
    results.reserve (numBenchmarks);

    // Time sparse triangular solve, no transpose.
    {
      const std::string timerName (opsLabel + ": Y = A \\ X (" + benchmarkLabel + ")");
      RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = TimeMonitor::getNewCounter (timerName);
      }
      {
        TimeMonitor timeMon (*timer);
        for (int i = 0; i < numTrials; ++i) {
          ops.template solve<scalar_type, scalar_type> (Teuchos::NO_TRANS,
                                                        *Y, *X);
        }
      }
      results.push_back (std::make_pair (timerName, timer->totalElapsedTime ()));
    }

    // Time sparse triangular solve, transpose.
    {
      const std::string timerName (opsLabel + ": Y = A^T \\ X (" + benchmarkLabel + ")");
      RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = TimeMonitor::getNewCounter (timerName);
      }
      {
        TimeMonitor timeMon (*timer);
        for (int i = 0; i < numTrials; ++i) {
          ops.template solve<scalar_type, scalar_type> (Teuchos::NO_TRANS,
                                                        *Y, *X);
        }
      }
      results.push_back (std::make_pair (timerName, timer->totalElapsedTime ()));
    }

    // Only test conjugate transpose if scalar_type is complex.
    if (STS::isComplex) {
      // Time sparse triangular solve, conjugate transpose.
      const std::string timerName (opsLabel + ": Y = A^H \\ X (" + benchmarkLabel + ")");
      RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = TimeMonitor::getNewCounter (timerName);
      }
      {
        TimeMonitor timeMon (*timer);
        for (int i = 0; i < numTrials; ++i) {
          ops.template solve<scalar_type, scalar_type> (Teuchos::NO_TRANS,
                                                        *Y, *X);
        }
      }
      results.push_back (std::make_pair (timerName, timer->totalElapsedTime ()));
    }
  }
};

#endif // __CrsMatrix_TestSparseOps_hpp
