// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tsqr.hpp
/// \brief Parallel Tall Skinny QR (TSQR) implementation

#ifndef TSQR_TSQR_HPP
#define TSQR_TSQR_HPP

#include "Tsqr_ApplyType.hpp"
#include "Tsqr_Matrix.hpp"
#include "Tsqr_MessengerBase.hpp"
#include "Tsqr_DistTsqr.hpp"
#include "Tsqr_SequentialTsqr.hpp"
#include "Tsqr_Util.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace TSQR {
/// \class Tsqr
/// \brief Parallel Tall Skinny QR (TSQR) factorization
/// \author Mark Hoemmen
///
/// This class computes the parallel Tall Skinny QR (TSQR)
/// factorization of a matrix distributed in block rows across one
/// or more MPI processes.  The parallel critical path length for
/// TSQR is independent of the number of columns in the matrix,
/// unlike ScaLAPACK's comparable QR factorization (P_GEQR2),
/// Modified Gram-Schmidt, or Classical Gram-Schmidt.
///
/// \tparam LocalOrdinal Index type that can address all elements of
///   a matrix, when treated as a 1-D array.  That is, for A[i +
///   LDA*j], the index i + LDA*j must fit in a LocalOrdinal.
///
/// \tparam Scalar The type of the matrix entries.
///
/// \tparam NodeTsqrType The intranode (single-MPI-process) part of
///   TSQR.  Defaults to \c SequentialTsqr, which provides a
///   sequential cache-blocked implementation.  Any class
///   implementing the same compile-time interface is valid.  We
///   provide \c NodeTsqr as an archetype of the "NodeTsqrType"
///   concept, but it is not necessary that NodeTsqrType derive from
///   that abstract base class.  Inheriting from \c NodeTsqr is
///   useful, though, because it provides default implementations of
///   some routines that are not performance-critical.
///
/// \note TSQR only needs to know about the local ordinal type (used
///   to index matrix entries on a single process), not about the
///   global ordinal type (used to index matrix entries globally,
///   i.e., over all processes).  For some distributed linear
///   algebra libraries, such as Epetra, the local and global
///   ordinal types are the same (int, in the case of Epetra).  For
///   other distributed linear algebra libraries, such as Tpetra,
///   the local and global ordinal types may be different.
template <class LocalOrdinal,
          class Scalar>
class Tsqr {
 public:
  typedef MatView<LocalOrdinal, Scalar> mat_view_type;
  typedef MatView<LocalOrdinal, const Scalar> const_mat_view_type;
  typedef Matrix<LocalOrdinal, Scalar> matrix_type;

  typedef Scalar scalar_type;
  typedef LocalOrdinal ordinal_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType magnitude_type;

  using node_tsqr_type = NodeTsqr<LocalOrdinal, Scalar>;
  using dist_tsqr_type = DistTsqr<LocalOrdinal, Scalar>;
  typedef typename Teuchos::RCP<node_tsqr_type> node_tsqr_ptr;
  typedef typename Teuchos::RCP<dist_tsqr_type> dist_tsqr_ptr;
  /// \typedef rank_type
  /// \brief "Rank" here means MPI rank, not linear algebra rank.
  typedef typename dist_tsqr_type::rank_type rank_type;

  using NodeOutput = typename node_tsqr_type::factor_output_type;
  using DistOutput = typename dist_tsqr_type::FactorOutput;

  /// \typedef FactorOutput
  /// \brief Return value of \c factor().
  ///
  /// Part of the implicit representation of the Q factor returned
  /// by \c factor().  The other part of that representation is
  /// stored in the A matrix on output.
  using FactorOutput =
      std::pair<Teuchos::RCP<NodeOutput>, DistOutput>;

  /// \brief Constructor
  ///
  /// \param nodeTsqr [in/out] Previously initialized NodeTsqrType
  ///   object.  This takes care of the intraprocess part of TSQR.
  ///
  /// \param distTsqr [in/out] Previously initialized DistTsqrType
  ///   object.  This takes care of the interprocess part of TSQR.
  Tsqr(const node_tsqr_ptr& nodeTsqr,
       const dist_tsqr_ptr& distTsqr)
    : nodeTsqr_(nodeTsqr)
    , distTsqr_(distTsqr) {
    TEUCHOS_ASSERT(!nodeTsqr_.is_null());
    TEUCHOS_ASSERT(!distTsqr_.is_null());
  }

  /// \brief Cache size hint in bytes used by the intraprocess part
  ///   of TSQR.
  ///
  /// This value may differ from the cache size hint given to the
  /// constructor of the NodeTsqrType object, since that constructor
  /// input is merely a suggestion.
  size_t cache_size_hint() const {
    return nodeTsqr_->cache_size_hint();
  }

  /// \brief Does the R factor have a nonnegative diagonal?
  ///
  /// Tsqr implements a QR factorization (of a distributed matrix).
  /// Some, but not all, QR factorizations produce an R factor whose
  /// diagonal may include negative entries.  This Boolean tells you
  /// whether Tsqr promises to compute an R factor whose diagonal
  /// entries are all nonnegative.
  ///
  bool QR_produces_R_factor_with_nonnegative_diagonal() const {
    // Tsqr computes an R factor with nonnegative diagonal, if and
    // only if all QR factorization steps (both intraprocess and
    // interprocess) produce an R factor with a nonnegative
    // diagonal.
    return nodeTsqr_->QR_produces_R_factor_with_nonnegative_diagonal() &&
           distTsqr_->QR_produces_R_factor_with_nonnegative_diagonal();
  }

  /// \brief Whether the implementation wants device memory for
  ///   "large" arrays, like the input matrix, and the output Q
  ///   factor or C apply result.
  bool wants_device_memory() const {
    return nodeTsqr_->wants_device_memory();
  }

  /// \brief Compute QR factorization with explicit Q factor: "raw"
  ///   arrays interface, for column-major data.
  ///
  /// This method computes the "thin" QR factorization, like
  /// Matlab's [Q,R] = qr(A,0), of a matrix A with more rows than
  /// columns.  Like Matlab, it computes the Q factor "explicitly,"
  /// that is, as a matrix represented in the same format as in the
  /// input matrix A (rather than the implicit representation
  /// returned by factor()).  You may use this to orthogonalize a
  /// block of vectors which are the columns of A on input.
  ///
  /// Calling this method may be faster than calling factor() and
  /// explicit_Q() in sequence, if you know that you only want the
  /// explicit version of the Q factor.  This method is especially
  /// intended for orthogonalizing the columns of a distributed
  /// block of vectors, stored on each process as a column-major
  /// matrix.  Examples in Trilinos include Tpetra::MultiVector and
  /// Epetra_MultiVector.
  ///
  /// \warning No output array may alias any other input or output
  ///   array.  In particular, Q may <i>not</i> alias A.  If you try
  ///   this, you will certainly give you the wrong answer.
  ///
  /// \param numRows [in] Number of rows in my process' part of A.
  ///   Q must have the same number of rows as A on each process.
  ///
  /// \param numCols [in] Number of columns in A, Q, and R.
  ///
  /// \param A [in/out] On input: my process' part of the matrix to
  ///   factor; the matrix is distributed over the participating
  ///   processors.  My process' part of the matrix is stored in
  ///   column-major order with column stride LDA.  The columns of A
  ///   on input are the vectors to orthogonalize.  On output:
  ///   overwritten with garbage.
  ///
  /// \param LDA [in] Leading dimension (column stride) of my
  ///   process' part of A.
  ///
  /// \param Q [out] On output: my process' part of the explicit Q
  ///   factor.  My process' part of the matrix is stored in
  ///   column-major order with column stride LDQ (which may differ
  ///   from LDA).  Q may <i>not</i> alias A.
  ///
  /// \param LDQ [in] Leading dimension (column stride) of my
  ///   process' part of Q.
  ///
  /// \param R [out] On output: the R factor, which is square with
  ///   the same number of rows and columns as the number of columns
  ///   in A.  The R factor is replicated on all processes in the
  ///   matrix's communicator.  It is stored in column-major order
  ///   with column stride LDR.
  ///
  /// \param LDR [in] Leading dimension (column stride) of my
  ///   process' part of R.
  ///
  /// \param forceNonnegativeDiagonal [in] If true, then (if
  ///   necessary) do extra work (modifying both the Q and R
  ///   factors) in order to force the R factor to have a
  ///   nonnegative diagonal.
  void
  factorExplicit(const LocalOrdinal numRows,
                 const LocalOrdinal numCols,
                 Scalar A[],
                 const LocalOrdinal LDA,
                 Scalar Q[],
                 const LocalOrdinal LDQ,
                 Scalar R[],
                 const LocalOrdinal LDR,
                 const bool forceNonnegativeDiagonal = false) {
    constexpr bool contiguousCacheBlocks = false;
    this->factorExplicitRaw(numRows, numCols,
                            A, LDA, Q, LDQ, R, LDR,
                            contiguousCacheBlocks,
                            forceNonnegativeDiagonal);
  }

  /// \brief Generalization of factorExplicit (see above) that
  ///   permits input(/output) A and output Q to have contiguous
  ///   cache blocks.
  ///
  /// If contiguousCacheBlocks is true, then both A and Q are stored
  /// with contiguous cache blocks.  This method ignores LDA and LDQ
  /// in that case.  R still has the usual column-major storage
  /// format.
  void
  factorExplicitRaw(const LocalOrdinal numRows,
                    const LocalOrdinal numCols,
                    Scalar A[],
                    const LocalOrdinal LDA,
                    Scalar Q[],
                    const LocalOrdinal LDQ,
                    Scalar R[],
                    const LocalOrdinal LDR,
                    const bool contiguousCacheBlocks,
                    const bool forceNonnegativeDiagonal = false) {
    const char prefix[] = "TSQR::Tsqr::factorExplicitRaw: ";

    // Sanity checks for matrix dimensions.
    if (numRows < numCols) {
      std::ostringstream os;
      os << "In Tsqr::factorExplicitRaw: input matrix A has " << numRows
         << " local rows, and " << numCols << " columns.  The input "
                                              "matrix must have at least as many rows on each processor as "
                                              "there are columns.";
      throw std::invalid_argument(os.str());
    }

    // Check for quick exit, based on matrix dimensions.
    if (numCols == 0) {
      return;
    }

    // Fill R initially with zeros.
    mat_view_type R_view(numCols, numCols, R, LDR);
    try {
      deep_copy(R_view, Scalar{});
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, prefix << "deep_copy(R_view, 0.0) threw: " << e.what());
    }

    // Compute the local QR factorization, in place in A, with the R
    // factor written to R.
    Teuchos::RCP<NodeOutput> nodeResults;
    try {
      nodeResults =
          nodeTsqr_->factor(numRows, numCols, A, LDA, R, LDR,
                            contiguousCacheBlocks);
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, prefix << "nodeTsqr_->factor(...) threw: " << e.what());
    }

    // Prepare the output matrix Q by filling with zeros.
    try {
      nodeTsqr_->fill_with_zeros(numRows, numCols, Q, LDQ,
                                 contiguousCacheBlocks);
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, prefix << "nodeTsqr_->fill_with_zeros(...) threw: " << e.what());
    }

    // Wrap the output matrix Q in a "view."
    mat_view_type Q_rawView(numRows, numCols, Q, LDQ);
    // Wrap the uppermost cache block of Q.  We will need to extract
    // its numCols x numCols uppermost block below.  We can't just
    // extract the numCols x numCols top block from all of Q, in
    // case Q is arranged using contiguous cache blocks.
    mat_view_type Q_top_block =
        nodeTsqr_->top_block(Q_rawView, contiguousCacheBlocks);
    if (Q_top_block.extent(0) < numCols) {
      std::ostringstream os;
      os << "The top block of Q has too few rows.  This means that the "
         << "the intranode TSQR implementation has a bug in its top_block"
         << "() method.  The top block should have at least " << numCols
         << " rows, but instead has only " << Q_top_block.extent(1)
         << " rows.";
      throw std::logic_error(os.str());
    }
    // Use the numCols x numCols top block of Q and the local R
    // factor (computed above) to compute the distributed-memory
    // part of the QR factorization.
    {
      mat_view_type Q_top(numCols, numCols, Q_top_block.data(),
                          Q_top_block.stride(1));

      if (nodeTsqr_->wants_device_memory()) {
        // DistTsqr doesn't know what to do with device memory, so
        // if Q_top is device memory, we need to work in a host copy
        // and copy back to Q_top.  Q_top is an output argument
        // here, so we can just fill Q_top_copy with zeros.
        matrix_type Q_top_copy(Q_top.extent(0), Q_top.extent(1),
                               Scalar{});
        try {
          distTsqr_->factorExplicit(R_view, Q_top_copy.view(),
                                    forceNonnegativeDiagonal);
        } catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, prefix << "distTsqr_->"
                                                                         "factorExplicit (wants_device_memory()=true case) "
                                                                         "threw: "
                                                                      << e.what());
        }
        try {
          nodeTsqr_->copy_from_host(Q_top, Q_top_copy.view());
        } catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, prefix << "nodeTsqr_->"
                                                                         "copy_from_host threw: "
                                                                      << e.what());
        }
      } else {
        try {
          distTsqr_->factorExplicit(R_view, Q_top,
                                    forceNonnegativeDiagonal);
        } catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, prefix << "distTsqr_->"
                                                                         "factorExplicit (wants_device_memory()=false case) "
                                                                         "threw: "
                                                                      << e.what());
        }
      }
    }
    // Apply the local part of the Q factor to the result of the
    // distributed-memory QR factorization, to get the explicit Q
    // factor.
    try {
      nodeTsqr_->apply(ApplyType::NoTranspose,
                       numRows, numCols, A, LDA,
                       *nodeResults, numCols, Q, LDQ,
                       contiguousCacheBlocks);
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, prefix << "nodeTsqr_->"
                                                                     "apply threw: "
                                                                  << e.what());
    }

    // If necessary, and if the user asked, force the R factor to
    // have a nonnegative diagonal.
    if (forceNonnegativeDiagonal &&
        !QR_produces_R_factor_with_nonnegative_diagonal()) {
      // We ignore contiguousCacheBlocks here, since we're only
      // looking at the top block of Q.
      try {
        nodeTsqr_->force_nonnegative_diagonal(numRows, numCols,
                                              Q, LDQ, R, LDR);
      } catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, prefix << "nodeTsqr_->"
                                                                       "force_nonnegative_diagonal threw: "
                                                                    << e.what());
      }
    }
  }

  /// \brief Compute QR factorization of the global dense matrix A
  ///
  /// Compute the QR factorization of the tall and skinny dense
  /// matrix A.  The matrix A is distributed in a row block layout
  /// over all the MPI processes.  A_local contains the matrix data
  /// for this process.
  ///
  /// \param nrows_local [in] Number of rows of this node's local
  ///   component (A_local) of the matrix.  May differ on different
  ///   nodes.  Precondition: nrows_local >= ncols.
  ///
  /// \param ncols [in] Number of columns in the matrix to factor.
  ///   Should be the same on all nodes.
  ///   Precondition: nrows_local >= ncols.
  ///
  /// \param A_local [in,out] On input, this node's local component of
  ///   the matrix, stored as a general dense matrix in column-major
  ///   order.  On output, overwritten with an implicit representation
  ///   of the Q factor.
  ///
  /// \param lda_local [in] Leading dimension of A_local.
  ///   Precondition: lda_local >= nrows_local.
  ///
  /// \param R [out] The final R factor of the QR factorization of the
  ///   global matrix A.  An ncols by ncols upper triangular matrix with
  ///   leading dimension ldr.
  ///
  /// \param ldr [in] Leading dimension of the matrix R.
  ///
  /// \param contiguousCacheBlocks [in] Whether or not cache blocks
  ///   of A_local are stored contiguously.  The default value of
  ///   false means that A_local uses ordinary column-major
  ///   (Fortran-style) order.  Otherwise, the details of the format
  ///   depend on the specific NodeTsqrType.  Tsqr's cache_block()
  ///   and un_cache_block() methods may be used to convert between
  ///   cache-blocked and non-cache-blocked (column-major order)
  ///   formats.
  ///
  /// \return Part of the representation of the implicitly stored Q
  ///   factor.  It should be passed into apply() or explicit_Q() as
  ///   the "factorOutput" parameter.  The other part of the
  ///   implicitly stored Q factor is stored in A_local (the input
  ///   is overwritten).  Both parts go together.
  FactorOutput
  factor(const LocalOrdinal nrows_local,
         const LocalOrdinal ncols,
         Scalar A_local[],
         const LocalOrdinal lda_local,
         Scalar R[],
         const LocalOrdinal ldr,
         const bool contiguousCacheBlocks = false) {
    mat_view_type R_view(ncols, ncols, R, ldr);
    deep_copy(R_view, Scalar{});
    auto nodeResults =
        nodeTsqr_->factor(nrows_local, ncols, A_local, lda_local,
                          R_view.data(), R_view.stride(1),
                          contiguousCacheBlocks);
    DistOutput distResults = distTsqr_->factor(R_view);
    return {nodeResults, distResults};
  }

  /// \brief Apply Q factor to the global dense matrix C
  ///
  /// Apply the Q factor (computed by factor() and represented
  /// implicitly) to the global dense matrix C, consisting of all
  /// nodes' C_local matrices stacked on top of each other.
  ///
  /// \param [in] If "N", compute Q*C.  If "T", compute Q^T * C.
  ///   If "H" or "C", compute Q^H * C.  (The last option may not
  ///   be implemented in all cases.)
  ///
  /// \param nrows_local [in] Number of rows of this node's local
  ///   component (C_local) of the matrix C.  Should be the same on
  ///   this node as the nrows_local argument with which factor() was
  ///   called  Precondition: nrows_local >= ncols.
  ///
  /// \param ncols_Q [in] Number of columns in Q.  Should be the same
  ///   on all nodes.  Precondition: nrows_local >= ncols_Q.
  ///
  /// \param Q_local [in] Same as A_local output of factor()
  ///
  /// \param ldq_local [in] Same as lda_local of factor()
  ///
  /// \param factor_output [in] Return value of factor()
  ///
  /// \param ncols_C [in] Number of columns in C.  Should be the same
  ///   on all nodes.  Precondition: nrows_local >= ncols_C.
  ///
  /// \param C_local [in,out] On input, this node's local component of
  ///   the matrix C, stored as a general dense matrix in column-major
  ///   order.  On output, overwritten with this node's component of
  ///   op(Q)*C, where op(Q) = Q, Q^T, or Q^H.
  ///
  /// \param ldc_local [in] Leading dimension of C_local.
  ///   Precondition: ldc_local >= nrows_local.
  ///
  /// \param contiguousCacheBlocks [in] Whether or not the cache
  ///   blocks of Q and C are stored contiguously.
  void
  apply(const std::string& op,
        const LocalOrdinal nrows_local,
        const LocalOrdinal ncols_Q,
        const Scalar Q_local[],
        const LocalOrdinal ldq_local,
        const FactorOutput& factor_output,
        const LocalOrdinal ncols_C,
        Scalar C_local[],
        const LocalOrdinal ldc_local,
        const bool contiguousCacheBlocks = false) {
    ApplyType applyType(op);

    // This determines the order in which we apply the intranode
    // part of the Q factor vs. the internode part of the Q factor.
    const bool transposed = applyType.transposed();

    // View of this node's local part of the matrix C.
    mat_view_type C_view(nrows_local, ncols_C, C_local, ldc_local);

    // Identify top ncols_C by ncols_C block of C_local.
    // top_block() will set C_top_view to have the correct leading
    // dimension, whether or not cache blocks are stored
    // contiguously.
    //
    // C_top_view is the topmost cache block of C_local.  It has at
    // least as many rows as columns, but it likely has more rows
    // than columns.
    mat_view_type C_view_top_block =
        nodeTsqr_->top_block(C_view, contiguousCacheBlocks);

    // View of the topmost ncols_C by ncols_C block of C.
    mat_view_type C_top_view(ncols_C, ncols_C, C_view_top_block.data(),
                             C_view_top_block.stride(1));

    // DistTsqr doesn't know what to do with device memory, so we
    // need to copy the top block of C if applicable.  The NodeTsqr
    // implementation can decide if that's necessary.
    //
    // That "matrix_type C_top" is the temporary copy of C_top_view.
    // C_top_view here is the "top block of C" that might live in
    // device memory.

    if (!transposed) {
      // C_top (small compact storage) gets a deep copy of the top
      // ncols_C by ncols_C block of C_local.
      matrix_type C_top = nodeTsqr_->copy_to_host(C_top_view);
      // Compute in place on all processes' C_top blocks.
      distTsqr_->apply(applyType, C_top.extent(1), ncols_Q,
                       C_top.data(), C_top.stride(1),
                       factor_output.second);
      // Copy result back to the top block of C_local.
      nodeTsqr_->copy_from_host(C_top_view, C_top.view());
      // Apply the local Q factor to C_local.
      nodeTsqr_->apply(applyType, nrows_local, ncols_Q,
                       Q_local, ldq_local, *(factor_output.first),
                       ncols_C, C_local, ldc_local,
                       contiguousCacheBlocks);
    } else {
      // Apply the (transpose of the) local Q factor to C_local.
      nodeTsqr_->apply(applyType, nrows_local, ncols_Q,
                       Q_local, ldq_local, *(factor_output.first),
                       ncols_C, C_local, ldc_local,
                       contiguousCacheBlocks);

      // C_top (small compact storage) gets a deep copy of the top
      // ncols_C by ncols_C block of C_local.
      matrix_type C_top = nodeTsqr_->copy_to_host(C_top_view);

      // Compute in place on all processors' C_top blocks.
      distTsqr_->apply(applyType, ncols_C, ncols_Q, C_top.data(),
                       C_top.stride(1), factor_output.second);
      // Copy result back to the top block of C_local.
      nodeTsqr_->copy_from_host(C_top_view, C_top.view());
    }
  }

  /// \brief Compute the explicit Q factor from result of factor().
  ///
  /// Compute the explicit version of the Q factor computed by
  /// factor() and represented implicitly (via Q_local_in and
  /// factor_output).
  ///
  /// \param nrows_local [in] Number of rows of this node's local
  ///   component (Q_local_in) of the matrix Q_local_in.  Also, the
  ///   number of rows of this node's local component (Q_local_out) of
  ///   the output matrix.  Should be the same on this node as the
  ///   nrows_local argument with which factor() was called.
  ///   Precondition: nrows_local >= ncols_Q_in.
  ///
  /// \param ncols_Q_in [in] Number of columns in the original matrix
  ///   A, whose explicit Q factor we are computing.  Should be the
  ///   same on all nodes.  Precondition: nrows_local >= ncols_Q_in.
  ///
  /// \param Q_local_in [in] Same as A_local output of factor().
  ///
  /// \param ldq_local_in [in] Same as lda_local of factor()
  ///
  /// \param factorOutput [in] Return value of factor().
  ///
  /// \param ncols_Q_out [in] Number of columns of the explicit Q
  ///   factor to compute.  Should be the same on all nodes.
  ///
  /// \param Q_local_out [out] This node's component of the Q factor
  ///   (in explicit form).
  ///
  /// \param ldq_local_out [in] Leading dimension of Q_local_out.
  ///
  /// \param contiguousCacheBlocks [in] Whether or not cache blocks
  ///   in Q_local_in and Q_local_out are stored contiguously.
  void
  explicit_Q(const LocalOrdinal nrows_local,
             const LocalOrdinal ncols_Q_in,
             const Scalar Q_local_in[],
             const LocalOrdinal ldq_local_in,
             const FactorOutput& factorOutput,
             const LocalOrdinal ncols_Q_out,
             Scalar Q_local_out[],
             const LocalOrdinal ldq_local_out,
             const bool contiguousCacheBlocks = false) {
    nodeTsqr_->fill_with_zeros(nrows_local, ncols_Q_out, Q_local_out,
                               ldq_local_out, contiguousCacheBlocks);
    // "Rank" here means MPI rank, not linear algebra rank.
    const rank_type myRank = distTsqr_->rank();
    if (myRank == 0) {
      // View of this node's local part of the matrix Q_out.
      mat_view_type Q_out_view(nrows_local, ncols_Q_out,
                               Q_local_out, ldq_local_out);

      // View of the topmost cache block of Q_out.  It is
      // guaranteed to have at least as many rows as columns.
      mat_view_type Q_out_top =
          nodeTsqr_->top_block(Q_out_view, contiguousCacheBlocks);

      // Q_out_top is device memory, so we shouldn't write directly
      // to it.  Instead, let NodeTsqr fill it with the first
      // ncols_Q_out columns of the identity matrix.  Note that
      // we've already filled Q_out with zeros above.
      nodeTsqr_->set_diagonal_entries_to_one(Q_out_top);
    }
    apply("N", nrows_local,
          ncols_Q_in, Q_local_in, ldq_local_in, factorOutput,
          ncols_Q_out, Q_local_out, ldq_local_out,
          contiguousCacheBlocks);
  }

  /// \brief Compute Q*B
  ///
  /// Compute matrix-matrix product Q*B, where Q is nrows by ncols
  /// and B is ncols by ncols.  Respect cache blocks of Q.
  void
  Q_times_B(const LocalOrdinal nrows,
            const LocalOrdinal ncols,
            Scalar Q[],
            const LocalOrdinal ldq,
            const Scalar B[],
            const LocalOrdinal ldb,
            const bool contiguousCacheBlocks = false) const {
    // This requires no internode communication.  However, the work
    // is not redundant, since each MPI process has a different Q.
    nodeTsqr_->Q_times_B(nrows, ncols, Q, ldq, B, ldb,
                         contiguousCacheBlocks);
    // We don't need a barrier after this method, unless users are
    // doing mean MPI_Get() things.
  }

  /// \brief Reveal the rank of the R factor, using the SVD.
  ///
  /// Compute the singular value decomposition (SVD) of the R
  /// factor: \f$R = U \Sigma V^*\f$, not in place.  Use the
  /// resulting singular values to compute the numerical rank of R,
  /// with respect to the relative tolerance tol.  If R is full
  /// rank, return without modifying R.  If R is not full rank,
  /// overwrite R with \f$\Sigma \cdot V^*\f$.
  ///
  /// \return Numerical rank of R: 0 <= rank <= ncols.
  LocalOrdinal
  reveal_R_rank(const LocalOrdinal ncols,
                Scalar R[],
                const LocalOrdinal ldr,
                Scalar U[],
                const LocalOrdinal ldu,
                const magnitude_type& tol) const {
    // Forward the request to the intranode TSQR implementation.
    // Currently, this work is performed redundantly on all MPI
    // processes, without communication or agreement.
    //
    // FIXME (mfh 26 Aug 2010) This be a problem if your cluster is
    // heterogeneous, because then you might obtain different
    // integer rank results.  This is because heterogeneous nodes
    // might each compute the rank-revealing decomposition with
    // slightly different rounding error.
    return nodeTsqr_->reveal_R_rank(ncols, R, ldr, U, ldu, tol);
  }

  /// \brief Rank-revealing decomposition
  ///
  /// Using the R factor and explicit Q factor from
  /// factorExplicit(), compute the singular value decomposition
  /// (SVD) of R (\f$R = U \Sigma V^*\f$).  If R is full rank (with
  /// respect to the given relative tolerance tol), don't change Q
  /// or R.  Otherwise, compute \f$Q := Q \cdot U\f$ and \f$R :=
  /// \Sigma V^*\f$ in place (the latter may be no longer upper
  /// triangular).
  ///
  /// \param nrows [in] Number of rows in Q (on the calling process)
  ///
  /// \param ncols [in] Number of columns in Q (on the calling
  ///   process)
  ///
  /// \param Q [in/out] On input: explicit Q factor computed by
  ///   factorExplicit().  (Must be an orthogonal resp. unitary
  ///   matrix.)  On output: If R is of full numerical rank with
  ///   respect to the tolerance tol, Q is unmodified.  Otherwise, Q
  ///   is updated so that the first rank columns of Q are a basis
  ///   for the column space of A (the original matrix whose QR
  ///   factorization was computed by factorExplicit()).  The
  ///   remaining columns of Q are a basis for the null space of A.
  ///
  /// \param ldq [in] Leading dimension / column stride of Q
  ///
  /// \param R [in/out] On input: ncols by ncols upper triangular
  ///   matrix with leading dimension ldr >= ncols.  On output: if
  ///   input is full rank, R is unchanged on output.  Otherwise, if
  ///   \f$R = U \Sigma V^*\f$ is the SVD of R, on output R is
  ///   overwritten with $\Sigma \cdot V^*$.  This is also an ncols by
  ///   ncols matrix, but may not necessarily be upper triangular.
  ///
  /// \param ldr [in] Leading dimension / column stride of R
  ///
  /// \param tol [in] Relative tolerance for computing the numerical
  ///   rank of the matrix R.
  ///
  /// \param contiguousCacheBlocks [in] Whether or not the cache
  ///   blocks of Q are stored contiguously.  Defaults to false,
  ///   which means that Q uses the ordinary column-major layout on
  ///   each MPI process.
  ///
  /// \return Rank \f$r\f$ of R: \f$ 0 \leq r \leq ncols\f$.
  LocalOrdinal
  revealRankRaw(const LocalOrdinal nrows,
                const LocalOrdinal ncols,
                Scalar Q[],
                const LocalOrdinal ldq,
                Scalar R[],
                const LocalOrdinal ldr,
                const magnitude_type& tol,
                const bool contiguousCacheBlocks = false) const {
    // Take the easy exit if available.
    if (ncols == 0) {
      return 0;
    }
    // FIXME (mfh 16 Jul 2010) We _should_ compute the SVD of R (as
    // the copy B) on Proc 0 only.  This would ensure that all
    // processors get the same SVD and rank (esp. in a heterogeneous
    // computing environment).  For now, we just do this computation
    // redundantly, and hope that all the returned rank values are
    // the same.
    matrix_type U(ncols, ncols, Scalar{});
    const ordinal_type rank =
        reveal_R_rank(ncols, R, ldr, U.data(), U.stride(1), tol);
    if (rank < ncols) {
      // If R is not full rank: reveal_R_rank() already computed
      // the SVD \f$R = U \Sigma V^*\f$ of (the input) R, and
      // overwrote R with \f$\Sigma V^*\f$.  Now, we compute \f$Q
      // := Q \cdot U\f$, respecting cache blocks of Q.
      Q_times_B(nrows, ncols, Q, ldq, U.data(), U.stride(1),
                contiguousCacheBlocks);
    }
    return rank;
  }

  /// \brief Cache-block A_in into A_out.
  ///
  /// Cache-block the given A_in matrix, writing the results to
  /// A_out.
  void
  cache_block(const LocalOrdinal nrows_local,
              const LocalOrdinal ncols,
              Scalar A_local_out[],
              const Scalar A_local_in[],
              const LocalOrdinal lda_local_in) const {
    nodeTsqr_->cache_block(nrows_local, ncols,
                           A_local_out,
                           A_local_in, lda_local_in);
  }

  /// \brief Un-cache-block A_in into A_out.
  ///
  /// "Un"-cache-block the given A_in matrix, writing the results to
  /// A_out.
  void
  un_cache_block(const LocalOrdinal nrows_local,
                 const LocalOrdinal ncols,
                 Scalar A_local_out[],
                 const LocalOrdinal lda_local_out,
                 const Scalar A_local_in[]) const {
    nodeTsqr_->un_cache_block(nrows_local, ncols,
                              A_local_out, lda_local_out,
                              A_local_in);
  }

 private:
  node_tsqr_ptr nodeTsqr_;
  dist_tsqr_ptr distTsqr_;
};  // class Tsqr

}  // namespace TSQR

#endif  // TSQR_TSQR_HPP
