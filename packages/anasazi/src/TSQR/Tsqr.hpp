// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef __TSQR_Tsqr_hpp
#define __TSQR_Tsqr_hpp

#include <Tsqr_ApplyType.hpp>
#include <Tsqr_Matrix.hpp>
#include <Tsqr_MessengerBase.hpp>
#include <Tsqr_DistTsqr.hpp>
#include <Tsqr_SequentialTsqr.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <Tsqr_Util.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class TSQR factorization, in parallel for a tall skinny matrix
  /// distributed in block rows across one or more MPI processes.
  ///
  /// \note TSQR only needs to know about the local ordinal type, not
  /// about the global ordinal type.
  ///
  /// LocalOrdinal: index type that can address all elements of a
  /// matrix (when treated as a 1-D array, so for A[i + LDA*j], the
  /// number i + LDA*j must fit in a LocalOrdinal).
  ///
  /// Scalar: the type of the matrix entries.
  ///
  /// NodeTsqr: how Tsqr does computations on a single node.  Defaults
  /// to sequential (cache-blocked) TSQR.  TbbTsqr< LocalOrdinal,
  /// Scalar > is also valid.  We provide NodeTsqr.hpp as a model of
  /// the "node TSQR" concept, but it is not necessary that NodeTsqr
  /// here derive from that abstract base class.  (This is
  /// compile-time polymorphism.)
  template< class LocalOrdinal, 
	    class Scalar, 
	    class NodeTsqr = SequentialTsqr< LocalOrdinal, Scalar > >
  class Tsqr {
  private:
    typedef MatView< LocalOrdinal, Scalar > mat_view;
    typedef ConstMatView< LocalOrdinal, Scalar > const_mat_view;

  public:
    typedef Scalar scalar_type;
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
    typedef LocalOrdinal ordinal_type;
    typedef typename NodeTsqr::FactorOutput NodeOutput;
    typedef typename DistTsqr< LocalOrdinal, Scalar >::FactorOutput DistOutput;
    typedef std::pair< NodeOutput, DistOutput > FactorOutput;

    /// TSQR constructor; sets up tuning parameters
    ///
    /// \param node_tsqr [in/out] Previously initialized NodeTsqr
    ///   object.
    ///
    /// \param messenger [in] object handling internode communication
    ///   for all nodes participating in the factorization.
    Tsqr (NodeTsqr& node_tsqr,
	  MessengerBase< Scalar >* const messenger) :
      node_tsqr_ (node_tsqr),
      dist_ (messenger),
      messenger_ (messenger) {}

    /// Cache block size (in bytes) used by the underlying per-node
    /// TSQR implementation.
    size_t cache_block_size() const { return node_tsqr_.cache_block_size(); }

    /// Whether or not all diagonal entries of the R factor computed
    /// by the QR factorization are guaranteed to be nonnegative.
    ///
    /// \note This property holds if all QR factorization steps (both
    ///   intranode and internode) produce an R factor with a
    ///   nonnegative diagonal.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      return node_tsqr_.QR_produces_R_factor_with_nonnegative_diagonal() &&
	dist_.QR_produces_R_factor_with_nonnegative_diagonal();
    }

    /// \brief Compute QR factorization of the global dense matrix A
    ///
    /// Compute the QR factorization of the dense matrix A, consisting
    /// of consisting of all nodes' A_local matrices stacked on top of
    /// each other.
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
    /// \return Three arrays with the same number of elements: 
    ///   1. The results of sequential TSQR on this node
    ///   2. An array of "local Q factors" from parallel TSQR
    ///   3. An array of "local tau arrays" from parallel TSQR
    ///   These together form an implicit representation of
    ///   the Q factor.  They should be passed into the apply() and
    ///   explicit_Q() functions as the "factor_output" parameter.
    FactorOutput
    factor (const LocalOrdinal nrows_local,
	    const LocalOrdinal ncols, 
	    Scalar A_local[],
	    const LocalOrdinal lda_local,
	    Scalar R[],
	    const LocalOrdinal ldr,
	    const bool contiguous_cache_blocks = false)
    {
      MatView< LocalOrdinal, Scalar > R_view (ncols, ncols, R, ldr);
      R_view.fill (Scalar(0));
      NodeOutput node_results = 
	node_tsqr_.factor (nrows_local, ncols, A_local, lda_local, 
			   R_view.get(), R_view.lda(), 
			   contiguous_cache_blocks);
      DistOutput dist_results = dist_.factor (R_view);
      return std::make_pair (node_results, dist_results);
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
    /// \param contiguous_cache_blocks [in] Whether or not the cache
    ///   blocks of Q and C are stored contiguously.
    ///
    void
    apply (const std::string& op,
	   const LocalOrdinal nrows_local,
	   const LocalOrdinal ncols_Q,
	   const Scalar Q_local[],
	   const LocalOrdinal ldq_local,
	   const FactorOutput& factor_output,
	   const LocalOrdinal ncols_C,
	   Scalar C_local[],
	   const LocalOrdinal ldc_local,
	   const bool contiguous_cache_blocks = false)
    {
      ApplyType apply_type (op);

      // "Conjugate transpose" means the same thing as "Transpose"
      // when Scalar is real.
      if (apply_type == ApplyType::ConjugateTranspose && 
	  ! ScalarTraits< Scalar >::is_complex)
	apply_type = ApplyType::Transpose;
      // This determines the order in which we apply the intranode
      // part of the Q factor vs. the internode part of the Q factor.
      const bool transposed = apply_type.transposed();

      // View of this node's local part of the matrix C.
      mat_view C_view (nrows_local, ncols_C, C_local, ldc_local);

      // Identify top ncols_C by ncols_C block of C_local.
      // top_block() will set C_top_view to have the correct leading
      // dimension, whether or not cache blocks are stored
      // contiguously.
      //
      // C_top_view is the topmost cache block of C_local.  It has at
      // least as many rows as columns, but it likely has more rows
      // than columns.
      mat_view C_view_top_block = 
	node_tsqr_.top_block (C_view, contiguous_cache_blocks);

      // View of the topmost ncols_C by ncols_C block of C.
      mat_view C_top_view (ncols_C, ncols_C, C_view_top_block.get(), 
			   C_view_top_block.lda());

      if (! transposed) 
	{
	  // C_top (small compact storage) gets a deep copy of the top
	  // ncols_C by ncols_C block of C_local.
      	  Matrix< LocalOrdinal, Scalar > C_top (C_top_view);

	  // Compute in place on all processors' C_top blocks.
	  dist_.apply (apply_type, C_top.ncols(), ncols_Q, C_top.get(), 
		       C_top.lda(), factor_output.second);

	  // Copy the result from C_top back into the top ncols_C by
	  // ncols_C block of C_local.
	  C_top_view.copy (C_top);

	  // Apply the local Q factor (in Q_local and
	  // factor_output.first) to C_local.
	  node_tsqr_.apply (apply_type, nrows_local, ncols_Q, 
			    Q_local, ldq_local, factor_output.first, 
			    ncols_C, C_local, ldc_local,
			    contiguous_cache_blocks);
	}
      else
	{
	  // Apply the (transpose of the) local Q factor (in Q_local
	  // and factor_output.first) to C_local.
	  node_tsqr_.apply (apply_type, nrows_local, ncols_Q, 
			    Q_local, ldq_local, factor_output.first, 
			    ncols_C, C_local, ldc_local,
			    contiguous_cache_blocks);

	  // C_top (small compact storage) gets a deep copy of the top
	  // ncols_C by ncols_C block of C_local.
      	  Matrix< LocalOrdinal, Scalar > C_top (C_top_view);

	  // Compute in place on all processors' C_top blocks.
	  dist_.apply (apply_type, ncols_C, ncols_Q, C_top.get(), 
		       C_top.lda(), factor_output.second);

	  // Copy the result from C_top back into the top ncols_C by
	  // ncols_C block of C_local.
	  C_top_view.copy (C_top);
	}
    }

    /// \brief Compute the explicit Q factor from factor()
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
    /// \param ncols_Q_out [in] Number of columns of the explicit Q
    ///   factor to compute.  Should be the same on all nodes.
    ///
    /// \param Q_local_out [out] This node's component of the Q factor
    ///   (in explicit form).
    ///
    /// \param ldq_local_out [in] Leading dimension of Q_local_out.
    ///
    /// \param factor_output [in] Return value of factor().
    ///
    /// \param comm [in] MPI communicator object for all nodes
    ///   participating in the operation.
    void
    explicit_Q (const LocalOrdinal nrows_local,
		const LocalOrdinal ncols_Q_in,
		const Scalar Q_local_in[],
		const LocalOrdinal ldq_local_in,
		const FactorOutput& factor_output,
		const LocalOrdinal ncols_Q_out,
		Scalar Q_local_out[],
		const LocalOrdinal ldq_local_out,
		const bool contiguous_cache_blocks = false)
    {
      node_tsqr_.fill_with_zeros (nrows_local, ncols_Q_out, Q_local_out,
				  ldq_local_out, contiguous_cache_blocks);
      const int my_rank = messenger_->rank();
      if (my_rank == 0)
	{
	  // View of this node's local part of the matrix Q_out.
	  mat_view Q_out_view (nrows_local, ncols_Q_out, 
			       Q_local_out, ldq_local_out);

	  // View of the topmost cache block of Q_out.  It is
	  // guaranteed to have at least as many rows as columns.
	  mat_view Q_out_top = 
	    node_tsqr_.top_block (Q_out_view, contiguous_cache_blocks);

	  // Fill (topmost cache block of) Q_out with the first
	  // ncols_Q_out columns of the identity matrix.
	  for (LocalOrdinal j = 0; j < ncols_Q_out; ++j)
	    Q_out_top(j, j) = Scalar(1);
	}
      apply ("N", nrows_local, 
	     ncols_Q_in, Q_local_in, ldq_local_in, factor_output,
	     ncols_Q_out, Q_local_out, ldq_local_out,
	     contiguous_cache_blocks);
    }

    /// \brief Compute Q*B
    ///
    /// Compute matrix-matrix product Q*B, where Q is nrows by ncols
    /// and B is ncols by ncols.  Respect cache blocks of Q.
    void
    Q_times_B (const LocalOrdinal nrows,
	       const LocalOrdinal ncols,
	       Scalar Q[],
	       const LocalOrdinal ldq,
	       const Scalar B[],
	       const LocalOrdinal ldb,
	       const bool contiguous_cache_blocks = false) const
    {
      // This requires no internode communication.
      node_tsqr_.Q_times_B (nrows, ncols, Q, ldq, B, ldb, 
			    contiguous_cache_blocks);
      // We don't need a barrier after this method, unless users are
      // doing mean MPI_Get() things.
    }

    /// Compute SVD \f$R = U \Sigma V^*\f$, not in place.  Use the
    /// resulting singular values to compute the numerical rank of R,
    /// with respect to the relative tolerance tol.  If R is full
    /// rank, return without modifying R.  If R is not full rank,
    /// overwrite R with \f$\Sigma \cdot V^*\f$.
    ///
    /// \return Numerical rank of R: 0 <= rank <= ncols.
    LocalOrdinal
    reveal_R_rank (const LocalOrdinal ncols,
		   Scalar R[],
		   const LocalOrdinal ldr,
		   Scalar U[],
		   const LocalOrdinal ldu,
		   const magnitude_type tol) const 
    {
      return node_tsqr_.reveal_R_rank (ncols, R, ldr, U, ldu, tol);
    }

    /// \brief Rank-revealing decomposition
    ///
    /// Using the R factor from factor() and the explicit Q factor
    /// from explicit_Q(), compute the SVD of R (\f$R = U \Sigma
    /// V^*\f$).  R.  If R is full rank (with respect to the given
    /// relative tolerance tol), don't change Q or R.  Otherwise,
    /// compute \f$Q := Q \cdot U\f$ and \f$R := \Sigma V^*\f$ in
    /// place (the latter may be no longer upper triangular).
    ///
    /// \param R [in/out] On input: ncols by ncols upper triangular
    /// matrix with leading dimension ldr >= ncols.  On output: if
    /// input is full rank, R is unchanged on output.  Otherwise, if
    /// \f$R = U \Sigma V^*\f$ is the SVD of R, on output R is
    /// overwritten with $\Sigma \cdot V^*$.  This is also an ncols by
    /// ncols matrix, but may not necessarily be upper triangular.
    ///
    /// \return Rank \f$r\f$ of R: \f$ 0 \leq r \leq ncols\f$.
    ///
    LocalOrdinal
    reveal_rank (const LocalOrdinal nrows,
		 const LocalOrdinal ncols,
		 Scalar Q[],
		 const LocalOrdinal ldq,
		 Scalar R[],
		 const LocalOrdinal ldr,
		 const magnitude_type tol,
		 const bool contiguous_cache_blocks = false) const
    {
      // Take the easy exit if available.
      if (ncols == 0)
	return 0;

      //
      // FIXME (mfh 16 Jul 2010) We _should_ compute the SVD of R (as
      // the copy B) on Proc 0 only.  This would ensure that all
      // processors get the same SVD and rank (esp. in a heterogeneous
      // computing environment).  For now, we just do this computation
      // redundantly, and hope that all the returned rank values are
      // the same.
      //
      Matrix< LocalOrdinal, Scalar > U (ncols, ncols, Scalar(0));
      const LocalOrdinal rank = 
	reveal_R_rank (ncols, R, ldr, U.get(), U.lda(), tol);
      
      if (rank < ncols)
	{
	  // cerr << ">>> Rank of R: " << rank << " < ncols=" << ncols << endl;
	  // cerr << ">>> Resulting U:" << endl;
	  // print_local_matrix (cerr, ncols, ncols, R, ldr);
	  // cerr << endl;

	  // If R is not full rank: reveal_R_rank() already computed
	  // the SVD \f$R = U \Sigma V^*\f$ of (the input) R, and
	  // overwrote R with \f$\Sigma V^*\f$.  Now, we compute \f$Q
	  // := Q \cdot U\f$, respecting cache blocks of Q.
	  Q_times_B (nrows, ncols, Q, ldq, U.get(), U.lda(), 
		     contiguous_cache_blocks);
	}
      return rank;
    }

    /// \brief Cache-block A_in into A_out
    ///
    /// Cache-block the given A_in matrix, writing the results to
    /// A_out.
    void
    cache_block (const LocalOrdinal nrows_local,
		 const LocalOrdinal ncols,
		 Scalar A_local_out[],
		 const Scalar A_local_in[],
		 const LocalOrdinal lda_local_in) const
    {
      node_tsqr_.cache_block (nrows_local, ncols, A_local_out, A_local_in, lda_local_in);
    }

    /// \brief Un-cache-block A_in into A_out
    ///
    /// "Un"-cache-block the given A_in matrix, writing the results to
    /// A_out.
    void
    un_cache_block (const LocalOrdinal nrows_local,
		    const LocalOrdinal ncols,
		    Scalar A_local_out[],
		    const LocalOrdinal lda_local_out,		    
		    const Scalar A_local_in[]) const
    {
      node_tsqr_.un_cache_block (nrows_local, ncols, A_local_out, lda_local_out, A_local_in);
    }

  private:
    NodeTsqr& node_tsqr_;
    DistTsqr< LocalOrdinal, Scalar > dist_;
    MessengerBase< Scalar >* messenger_;
  }; // class Tsqr

} // namespace TSQR

#endif // __TSQR_Tsqr_hpp
