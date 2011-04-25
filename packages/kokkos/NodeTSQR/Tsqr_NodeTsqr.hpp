//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef __TSQR_Tsqr_NodeTsqr_hpp
#define __TSQR_Tsqr_NodeTsqr_hpp

#include <Tsqr_ApplyType.hpp>
#include <Tsqr_Lapack.hpp>
#include <Tsqr_Matrix.hpp>

#include <Teuchos_Describable.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class NodeTsqr
  /// \brief Common interface and functionality for intranode TSQR.
  ///
  /// NodeTsqr provides a generic interface for TSQR operations within
  /// a node ("intranode").  It also implements rank-revealing
  /// functionality used by all intranode TSQR implementations.
  ///
  /// NodeTsqr has three template parameters:  
  /// - the (local) ordinal type (how we index into a matrix)
  /// - the scalar type (the type of elements stored in the matrix)
  /// - the type returned by factor() (FactorOutputType)
  ///
  /// We template on FactorOutputType for compile-time polymorphism.
  /// This lets us define the factor method() without constraining
  /// subclasses to set up a class hierarchy for that type.  All
  /// NodeTsqr subclasses' factor() methods return some composition of
  /// std::pair, std::vector, and Scalar.  It would be silly to wrap
  /// that up in an abstract base class, since we don't intend that
  /// returned object to be polymorphic. 
  ///
  /// Templating on FactorOutputType means that we can't use run-time
  /// polymorphism to swap between NodeTsqr subclasses, since the
  /// latter are really subclasses of different NodeTsqr
  /// instantiations (i.e., different FactorOutputType types).  The
  /// point is to share common functionality, not a common interface;
  /// we only define the interface here as an easy syntactic
  /// enforcement of compile-time polymorphism (to ensure that
  /// subclasses implemented the right methods).  Run-time
  /// polymorphism of different NodeTsqr subclasses would not be
  /// useful.  This is because ultimately each subclass is bound to a
  /// Kokkos Node type, and those only use compile-time polymorphism.
  ///
  template<class Ordinal, class Scalar, class FactorOutputType>
  class NodeTsqr : public Teuchos::Describable {
  public:
    typedef Ordinal ordinal_type;
    typedef Scalar scalar_type;
    typedef FactorOutputType factor_output_type;

    //! Constructor
    NodeTsqr() {}

    //! Virtual destructor ensures safe polymorphic destruction.
    virtual ~NodeTsqr() {}

    /// \brief One-line description of this object.
    ///
    /// This implements Teuchos::Describable::description().
    /// Subclasses should override this to provide a more specific
    /// description of their implementation.  Subclasses may also
    /// implement Teuchos::Describable::describe(), which for this
    /// class has a simple default implementation that calls
    /// description() with appropriate indenting.
    virtual std::string description () const {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os << "NodeTsqr<Ordinal=" << TypeNameTraits<Ordinal>::name()
	 << ", Scalar=" << TypeNameTraits<Scalar>::name()
	 << ", ...>: Intranode Tall Skinny QR (TSQR)";
      return os.str();
    }

    //! Compute QR factorization (implicitly stored Q factor) of A.
    factor_output_type
    factor (const Ordinal nrows,
	    const Ordinal ncols, 
	    Scalar A[],
	    const Ordinal lda,
	    Scalar R[],
	    const Ordinal ldr,
	    const bool contiguousCacheBlocks=false)
    {
      return factorImpl (nrows, ncols, A, lda, R, ldr, contiguousCacheBlocks);
    }

    //! Apply implicit Q factor stored in Q and factorOutput to C.
    void
    apply (const ApplyType& applyType,
	   const Ordinal nrows,
	   const Ordinal ncols_Q,
	   const Scalar Q[],
	   const Ordinal ldq,
	   const FactorOutput& factorOutput,
	   const Ordinal ncols_C,
	   Scalar C[],
	   const Ordinal ldc,
	   const bool contiguousCacheBlocks=false) 
    {
      applyImpl (applyType, nrows, ncols_Q, Q, ldq, factorOutput, 
		 ncols_C, C, ldc, contiguousCacheBlocks);
    }

    //! Compute explicit Q factor from result of factor().
    void
    explicit_Q (const Ordinal nrows,
		const Ordinal ncols_Q,
		const Scalar Q[],
		const Ordinal ldq,
		const factor_output_type& factorOutput,
		const Ordinal ncols_C,
		Scalar C[],
		const Ordinal ldc,
		const bool contiguousCacheBlocks=false)
    {
      explicit_Q_impl (nrows, ncols_Q, Q, ldq, factorOutput, 
		       ncols_C, C, ldc, contiguousCacheBlocks);
    }

    /// \brief Cache block A_in into A_out.
    ///
    /// \param nrows [in] Number of rows in A_in and A_out.
    /// \param ncols [in] Number of columns in A_in and A_out.
    /// \param A_out [out] Result of cache-blocking A_in.
    /// \param A_in [in] Matrix to cache block, stored in column-major
    ///   order with leading dimension lda_in.
    /// \param lda_in [in] Leading dimension of A_in.  (See the LAPACK
    ///   documentation for a definition of "leading dimension.")
    ///   lda_in >= nrows.
    virtual void
    cache_block (const Ordinal nrows,
		 const Ordinal ncols, 
		 Scalar A_out[],
		 const Scalar A_in[],
		 const Ordinal lda_in) const = 0;

    /// \brief Un - cache block A_in into A_out.
    ///
    /// A_in is a matrix produced by \c cache_block().  It is
    /// organized as contiguously stored cache blocks.  This method
    /// reorganizes A_in into A_out as an ordinary matrix stored in
    /// column-major order with leading dimension lda_out.
    ///
    /// \param nrows [in] Number of rows in A_in and A_out.
    /// \param ncols [in] Number of columns in A_in and A_out.
    /// \param A_out [out] Result of un-cache-blocking A_in.
    ///   Matrix stored in column-major order with leading
    ///   dimension lda_out.
    /// \param lda_out [in] Leading dimension of A_out.  (See the
    ///   LAPACK documentation for a definition of "leading
    ///   dimension.")  lda_out >= nrows.
    /// \param A_in [in] Matrix to un-cache-block.
    virtual void
    un_cache_block (const Ordinal nrows,
		    const Ordinal ncols,
		    Scalar A_out[],
		    const Ordinal lda_out,		    
		    const Scalar A_in[]) const = 0;

    /// \brief Compute Q*B
    ///
    /// Compute matrix-matrix product Q*B, where Q is nrows by ncols
    /// and B is ncols by ncols.  Respect cache blocks of Q.
    void
    Q_times_B (const Ordinal nrows,
	       const Ordinal ncols,
	       Scalar Q[],
	       const Ordinal ldq,
	       const Scalar B[],
	       const Ordinal ldb,
	       const bool contiguousCacheBlocks=false) const
    {
      Q_times_B_impl (nrows, ncols, Q, ldq, B, ldb, contiguousCacheBlocks);
    }

    /// \brief Fill the nrows by ncols matrix A with zeros.
    /// 
    /// Fill the matrix A with zeros, in a way that respects the cache
    /// blocking scheme.
    ///
    /// \param nrows [in] Number of rows in A
    /// \param ncols [in] Number of columns in A 
    /// \param A [out] nrows by ncols column-major-order dense matrix 
    ///   with leading dimension lda
    /// \param lda [in] Leading dimension of A: lda >= nrows
    /// \param contiguousCacheBlocks [in] Whether the cache blocks
    ///   in A are stored contiguously.
    void
    fill_with_zeros (const Ordinal nrows,
		     const Ordinal ncols,
		     Scalar A[],
		     const Ordinal lda, 
		     const bool contiguousCacheBlocks=false)
    {
      return fill_with_zeros_impl (nrows, ncols, A, lda, contiguousCacheBlocks);
    }

    /// \brief Return view of topmost cache block of C
    ///
    /// \param C [in] Matrix (view), supporting the usual nrows(),
    ///   ncols(), get(), lda() interface.
    /// \param contiguousCacheBlocks [in] Whether the cache blocks
    ///   in C are stored contiguously.
    ///
    /// Return a view of the topmost cache block (on this node) of the
    /// given matrix C.  This is not necessarily square, though it
    /// must have at least as many rows as columns.  For a square
    /// ncols by ncols block, as needed by Tsqr::apply(), do as 
    /// follows:
    /// \code 
    /// MatrixViewType top = this->top_block (C, contig);
    /// MatView< Ordinal, Scalar > square (ncols, ncols, top.get(), top.lda());
    /// \endcode
    MatView<Ordinal, Scalar>
    top_block (const MatView<Ordinal, Scalar>& C,
	       const bool contiguousCacheBlocks=false) const
    {
      return top_block_impl (C, contiguousCacheBlocks);
    }

    /// \brief Does factor() compute R with nonnegative diagonal?
    ///
    /// When using a QR factorization to orthogonalize a block of
    /// vectors, computing an R factor with nonnegative diagonal
    /// ensures that in exact arithmetic, the result of the
    /// orthogonalization (orthogonalized vectors Q and their
    /// coefficients R) are the same as would be produced by
    /// Gram-Schmidt orthogonalization.
    virtual bool 
    QR_produces_R_factor_with_nonnegative_diagonal () const = 0;

    /// \brief Reveal rank of TSQR's R factor.
    ///
    /// Compute SVD \f$R = U \Sigma V^*\f$, not in place.  Use the
    /// resulting singular values to compute the numerical rank of R,
    /// with respect to the relative tolerance tol.  If R is full
    /// rank, return without modifying R.  If R is not full rank,
    /// overwrite R with \f$\Sigma \cdot V^*\f$.
    ///
    /// \return Numerical rank of R: 0 <= rank <= ncols.
    Ordinal
    reveal_R_rank (const Ordinal ncols,
		   Scalar R[],
		   const Ordinal ldr,
		   Scalar U[],
		   const Ordinal ldu,
		   const typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol) const;

    /// \brief Compute rank-revealing decomposition.
    ///
    /// Using the R factor from factor() and the explicit Q factor
    /// from explicit_Q(), compute the SVD of R (\f$R = U \Sigma
    /// V^*\f$).  R.  If R is full rank (with respect to the given
    /// relative tolerance tol), don't change Q or R.  Otherwise,
    /// compute \f$Q := Q \cdot U\f$ and \f$R := \Sigma V^*\f$ in
    /// place (the latter may be no longer upper triangular).
    ///
    /// \return Rank \f$r\f$ of R: \f$ 0 \leq r \leq ncols\f$.
    ///
    Ordinal
    reveal_rank (const Ordinal nrows,
		 const Ordinal ncols,
		 Scalar Q[],
		 const Ordinal ldq,
		 Scalar R[],
		 const Ordinal ldr,
		 const typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol,
		 const bool contiguousCacheBlocks=false) const;

  protected:

    /// \brief Implementation of factor().
    ///
    /// Compute QR factorization (implicitly stored Q factor) of A.
    /// This is a separate method so that factor() can have a default
    /// parameter value.
    virtual factor_output_type
    factorImpl (const Ordinal nrows,
		const Ordinal ncols, 
		Scalar A[],
		const Ordinal lda,
		Scalar R[],
		const Ordinal ldr,
		const bool contiguousCacheBlocks) = 0;

    /// \brief Implementation of apply().
    ///
    /// Apply implicit Q factor stored in Q and factorOutput to C.
    /// This is a separate method so that apply() can have a default
    /// parameter value.
    virtual void
    applyImpl (const ApplyType& applyType,
	       const Ordinal nrows,
	       const Ordinal ncols_Q,
	       const Scalar Q[],
	       const Ordinal ldq,
	       const FactorOutput& factorOutput,
	       const Ordinal ncols_C,
	       Scalar C[],
	       const Ordinal ldc,
	       const bool contiguousCacheBlocks) = 0;

    /// \brief Implementation of explicit_Q().
    ///
    /// This is a separate method so that explicit_Q() can have a
    /// default parameter value.
    virtual void
    explicit_Q_impl (const Ordinal nrows,
		     const Ordinal ncols_Q,
		     const Scalar Q[],
		     const Ordinal ldq,
		     const factor_output_type& factorOutput,
		     const Ordinal ncols_C,
		     Scalar C[],
		     const Ordinal ldc,
		     const bool contiguousCacheBlocks) = 0;

    /// \brief Implementation of Q_times_B().
    ///
    /// This is a separate method so that Q_times_B() can have a
    /// default parameter value.
    virtual void
    Q_times_B_impl (const Ordinal nrows,
		    const Ordinal ncols,
		    Scalar Q[],
		    const Ordinal ldq,
		    const Scalar B[],
		    const Ordinal ldb,
		    const bool contiguousCacheBlocks) const = 0;

    /// \brief Implementation of fill_with_zeros().
    ///
    /// This is a separate method so that fill_with_zeros() can have a
    /// default parameter value.
    virtual void
    fill_with_zeros_impl (const Ordinal nrows,
			  const Ordinal ncols,
			  Scalar A[],
			  const Ordinal lda, 
			  const bool contiguousCacheBlocks) = 0;

    /// \brief Implementation of top_block().
    ///
    /// This is a separate method so that top_block() can have a
    /// default parameter value.
    virtual MatView<Ordinal, Scalar>
    top_block_impl (const MatView<Ordinal, Scalar>& C,
		    const bool contiguousCacheBlocks) const = 0;
  };


  template<class Ordinal, class Scalar, class FactorOutputType>
  Ordinal
  NodeTsqr<Ordinal, Scalar, FactorOutputType>::
  reveal_R_rank (const Ordinal ncols,
		 Scalar R[],
		 const Ordinal ldr,
		 Scalar U[],
		 const Ordinal ldu,
		 const typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol) const
  {
    using Teuchos::TypeNameTraits;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef typename STS::magnitudeType magnitude_type;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    TEST_FOR_EXCEPTION(tol < 0, std::invalid_argument,
		       "In NodeTsqr::reveal_R_rank: numerical rank tolerance "
		       "(tol = " << tol << ") is negative.");
    TEST_FOR_EXCEPTION(ncols < 0, std::invalid_argument,
		       "In NodeTsqr::reveal_R_rank: number of columns "
		       "(ncols = " << ncols << ") is negative.");
    TEST_FOR_EXCEPTION(ldr < ncols, std::invalid_argument,
		       "In NodeTsqr::reveal_R_ank: stride of R (ldr = " 
		       << ldr << ") is less than the number of columns "
		       "(ncols = " << ncols << ").");
    TEST_FOR_EXCEPTION(ldu < ncols, std::invalid_argument,
		       "In NodeTsqr::reveal_R_rank: stride of U (ldu = " 
		       << ldu << ") is less than the number of columns "
		       "(ncols = " << ncols << ")");

    // Zero columns always means rank zero.
    if (ncols == 0)
      return 0;
    //
    // Compute the SVD (singular value decomposition) of the R
    // factor, using LAPACK's GESVD routine.  We do so in a deep
    // copy (B) because LAPACK overwrites the input.  If the R
    // factor is full rank (expected to be the common case), we need
    // to leave it alone (so that it stays upper triangular).
    //
    LAPACK<Ordinal, Scalar> lapack;
    MatView<Ordinal, Scalar> R_view (ncols, ncols, R, ldr);
    Matrix<Ordinal, Scalar> B (R_view); // B := R (deep copy)
    MatView<Ordinal, Scalar> U_view (ncols, ncols, U, ldu);
    Matrix<Ordinal, Scalar> VT (ncols, ncols, Scalar(0));

    // Set up workspace for the SVD.
    std::vector<magnitude_type> svd_rwork (5*ncols);
    std::vector<magnitude_type> singular_values (ncols);
    Ordinal svd_lwork = -1; // -1 for LWORK query; will be changed
    int svd_info = 0;

    // LAPACK workspace ("LWORK") query for SVD.  The workspace
    // ("WORK") array is always of Scalar type, even in the complex
    // case.
    {
      // Exception messages in this scope all start with this.
      const char prefix[] = "In NodeTsqr::reveal_R_rank: LAPACK SVD (_GESVD) "
	"workspace query returned ";
      // std::logic_error messages in this scope all end with this.
      const char postfix[] = ".  Please report this bug to the Kokkos "
	"developers.";

      Scalar svd_lwork_scalar = Teuchos::ScalarTraits<Scalar>::zero();
      lapack.GESVD ("A", "A", ncols, ncols, B.get(), B.lda(), 
		    &singular_values[0], U_view.get(), U_view.lda(), 
		    VT.get(), VT.lda(), &svd_lwork_scalar, svd_lwork, 
		    &svd_rwork[0], &svd_info);
      // Failure of the LAPACK workspace query is a logic error (a
      // bug) because we have already validated the matrix
      // dimensions above.
      TEST_FOR_EXCEPTION(svd_info != 0, std::logic_error,
			 prefix << "a nonzero INFO = " << svd_info 
			 << postfix);
      // LAPACK returns the workspace array length as a Scalar.  We
      // have to convert it back to an Ordinal in order to allocate
      // the workspace array and pass it in to LAPACK as the LWORK
      // argument.  Ordinal definitely must be a signed type, since
      // LWORK = -1 indicates a workspace query.
      svd_lwork = static_cast<Ordinal> (svd_lwork_scalar);

      // LAPACK should always return an LWORK that fits in Ordinal,
      // but it's a good idea to check anyway.  We do so by checking
      // whether casting back from Ordinal to Scalar gives the same
      // original Scalar result.  This should work unless Scalar and
      // Ordinal are user-defined types with weird definitions of
      // the type casts.
      TEST_FOR_EXCEPTION(static_cast<Scalar> (svd_lwork) != svd_lwork_scalar,
			 std::logic_error,
			 prefix << "a workspace array length (LWORK) of type "
			 "Scalar=" << TypeNameTraits<Scalar>::name() 
			 << " that does not fit in an Ordinal=" 
			 << TypeNameTraits<Ordinal>::name() << " type.  "
			 "As a Scalar, LWORK=" << svd_lwork_scalar 
			 << ", but cast to Ordinal, LWORK=" << svd_lwork 
			 << postfix);
      // Make sure svd_lwork is nonnegative.  (Ordinal must be a
      // signed type, as we explain above, so this test should never
      // signal any unsigned-to-signed conversions from the compiler.
      // If it does, you're probably using the wrong Ordinal type.
      TEST_FOR_EXCEPTION(svd_lwork < 0, std::logic_error,
			 prefix << "a negative workspace array length (LWORK)"
			 " = " << svd_lwork << postfix);
    }
    // Allocate workspace for LAPACK's SVD routine.
    std::vector<Scalar> svd_work (svd_lwork);

    // Compute SVD $B := U \Sigma V^*$.  B is overwritten, which is
    // why we copied R into B (so that we don't overwrite R if R is
    // full rank).
    lapack.GESVD ("A", "A", ncols, ncols, B.get(), B.lda(), 
		  &singular_values[0], U_view.get(), U_view.lda(), 
		  VT.get(), VT.lda(), &svd_work[0], svd_lwork, 
		  &svd_rwork[0], &svd_info);
    //
    // Compute the numerical rank of B, using the given relative
    // tolerance and the computed singular values.  GESVD computes
    // singular values in decreasing order and they are all
    // nonnegative.  We know by now that ncols > 0.
    //
    // The tolerance "tol" is relative to the largest singular
    // value, which is the 2-norm of the matrix.
    const magnitude_type absolute_tolerance = tol * singular_values[0];
    Ordinal rank = 0;
    for (Ordinal k = 0; k < ncols; ++k)
      // First branch of the IF ensures correctness even if all the
      // singular values are zero and the absolute tolerance is
      // zero.  Recall that LAPACK sorts the singular values in
      // decreasing order.
      if (singular_values[k] == STM::zero() || 
	  singular_values[k] < absolute_tolerance)
	{
	  rank = k;
	  break;
	}
    // Don't modify Q or R, if R is full rank.
    if (rank < ncols)
      { //
	// R is not full rank.  
	//
	// 1. Compute \f$R := \Sigma V^*\f$.
	// 2. Return rank (0 <= rank < ncols).
	//
	// Compute R := \Sigma VT.  \Sigma is diagonal so we apply it
	// column by column (normally one would think of this as row by
	// row, but this "Hadamard product" formulation iterates more
	// efficiently over VT).  
	//
	// After this computation, R may no longer be upper triangular.
	// R may be zero if all the singular values are zero, but we
	// don't need to check for this case; it's rare in practice, and
	// the computations below will be correct regardless.
	for (Ordinal j = 0; j < ncols; ++j)
	  {
	    const Scalar* const VT_j = &VT(0,j);
	    Scalar* const R_j = &R_view(0,j);

	    for (Ordinal i = 0; i < ncols; ++i)
	      R_j[i] = singular_values[i] * VT_j[i];
	  }
      }
    return rank;
  }

  template<class Ordinal, class Scalar, class FactorOutputType>
  Ordinal
  NodeTsqr<Ordinal, Scalar, FactorOutputType>::
  reveal_rank (const Ordinal nrows,
	       const Ordinal ncols,
	       Scalar Q[],
	       const Ordinal ldq,
	       Scalar R[],
	       const Ordinal ldr,
	       const typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol,
	       const bool contiguousCacheBlocks) const
  {
    // Take the easy exit if available.
    if (ncols == 0)
      return 0;
    // Matrix to hold the left singular vectors of the R factor.
    Matrix<Ordinal, Scalar> U (ncols, ncols, Scalar(0));
    // Compute numerical rank of the R factor using the SVD.
    // Store the left singular vectors in U.
    const Ordinal rank = 
      reveal_R_rank (ncols, R, ldr, U.get(), U.ldu(), tol);

    // If R is full rank, we're done.  Otherwise, reveal_R_rank()
    // already computed the SVD \f$R = U \Sigma V^*\f$ of (the
    // input) R, and overwrote R with \f$\Sigma V^*\f$.  Now, we
    // compute \f$Q := Q \cdot U\f$, respecting cache blocks of Q.
    if (rank < ncols)
      Q_times_B (nrows, ncols, Q, ldq, U.get(), U.lda(), 
		 contiguousCacheBlocks);
    return rank;
  }

} // namespace TSQR


#endif __TSQR_Tsqr_NodeTsqr_hpp
