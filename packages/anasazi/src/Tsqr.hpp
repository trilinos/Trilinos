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
    typedef LocalOrdinal ordinal_type;
    typedef typename NodeTsqr::FactorOutput NodeOutput;
    typedef typename DistTsqr< LocalOrdinal, Scalar >::FactorOutput DistOutput;
    typedef std::pair< NodeOutput, DistOutput > FactorOutput;

    /// Identifying string for this particular type of TSQR
    /// implementation.  Different for each NodeTsqr parameter.
    // std::string identify() const;

    /// TSQR constructor; sets up tuning parameters
    ///
    /// \param cache_block_size [in] Size (in bytes) of cache block to
    ///   use in the node part of TSQR.  If zero or not specified, a
    ///   reasonable default is used.
    ///
    /// \param messenger [in] object handling internode communication
    ///   for all nodes participating in the factorization.
    Tsqr (NodeTsqr& node_tsqr,
	  MessengerBase< Scalar >* const messenger,
	  const size_t cache_block_size = 0) :
      node_tsqr_ (node_tsqr),
      dist_ (messenger),
      messenger_ (messenger) {}

    /// Cache block size (in bytes) used by the underlying per-node
    /// TSQR implementation.
    size_t cache_block_size() const { return node_tsqr_.cache_block_size(); }

    /// Whether or not all diagonal entries of the R factor computed
    /// by the QR factorization are guaranteed to be nonnegative.
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
    ///
    /// \param nrows_local [in] Number of rows of this node's local
    ///   component (C_local) of the matrix C.  Should be the same on
    ///   this node as the nrows_local argument with which factor() was
    ///   called  Precondition: nrows_local >= ncols.
    ///
    /// \param ncols_C [in] Number of columns in C.  Should be the same
    ///   on all nodes.  Precondition: nrows_local >= ncols_C.
    ///
    /// \param C_local [in,out] On input, this node's local component of
    ///   the matrix C, stored as a general dense matrix in column-major
    ///   order.  On output, overwritten with this node's component of 
    ///   op(Q)*C, where op(Q) = Q or Q^T.
    ///
    /// \param ldc_local [in] Leading dimension of C_local.  
    ///   Precondition: ldc_local >= nrows_local.
    ///
    /// \param Q_local [in] Same as A_local output of factor()
    ///
    /// \param ldq_local [in] Same as lda_local of factor()
    ///
    /// \param factor_output [in] Return value of factor()
    ///
    /// \param comm [in] Communicator object for all nodes
    ///   participating in the operation.
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
      const bool transposed = apply_type.transposed();

      if (apply_type == ApplyType::ConjugateTranspose)
	{
	  if (ScalarTraits< Scalar >::is_complex)
	    throw std::logic_error("TSQR::apply: applying Q^H for complex "
				   "scalar types not yet implemented");
	}

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
	  dist_.apply (op, C_top.ncols(), ncols_Q, C_top.get(), C_top.lda(),
		       factor_output.second);

	  // Copy the result from C_top back into the top ncols_C by
	  // ncols_C block of C_local.
	  C_top_view.copy (C_top);

	  // Apply the local Q factor (in Q_local and
	  // factor_output.first) to C_local.
	  node_tsqr_.apply (op, nrows_local, ncols_Q, Q_local, ldq_local, 
			    factor_output.first, ncols_C, C_local, ldc_local,
			    contiguous_cache_blocks);
	}
      else
	{
	  // Apply the (transpose of the) local Q factor (in Q_local
	  // and factor_output.first) to C_local.
	  node_tsqr_.apply (op, nrows_local, ncols_Q, Q_local, ldq_local, 
			    factor_output.first, ncols_C, C_local, ldc_local,
			    contiguous_cache_blocks);

	  // C_top (small compact storage) gets a deep copy of the top
	  // ncols_C by ncols_C block of C_local.
      	  Matrix< LocalOrdinal, Scalar > C_top (C_top_view);

	  // Compute in place on all processors' C_top blocks.
	  dist_.apply (op, ncols_C, ncols_Q, C_top.get(), C_top.lda(),
		       factor_output.second);

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
		Scalar Q_local_in[],
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

    /// Cache-block the given A_in matrix, writing the results to A_out.
    void
    cache_block (const LocalOrdinal nrows_local,
		 const LocalOrdinal ncols,
		 Scalar A_local_out[],
		 const Scalar A_local_in[],
		 const LocalOrdinal lda_local_in) const
    {
      node_tsqr_.cache_block (nrows_local, ncols, A_local_out, A_local_in, lda_local_in);
    }

    /// "Un"-cache-block the given A_in matrix, writing the results to A_out.
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
