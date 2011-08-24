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

/// \file Tsqr_NodeTsqr.hpp
/// \brief Common interface and functionality for intranode TSQR.
///
#ifndef __TSQR_Tsqr_NodeTsqr_hpp
#define __TSQR_Tsqr_NodeTsqr_hpp

#include <Tsqr_ApplyType.hpp>
#include <Tsqr_Lapack.hpp>
#include <Tsqr_Matrix.hpp>

#include <Teuchos_Describable.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include <vector>

namespace TSQR {

  /// \class NodeTsqr
  /// \brief Common interface and functionality for intranode TSQR.
  ///
  /// NodeTsqr provides a generic interface for TSQR operations within
  /// a node ("intranode").  It also implements rank-revealing
  /// functionality used by all intranode TSQR implementations.
  ///
  /// \tparam Ordinal The (local) Ordinal type; the type of indices
  ///   into a matrix on a node
  /// \tparam Scalar Tthe type of elements stored in the matrix
  /// \tparam FactorOutputType The type returned by factor().
  ///
  /// We template on FactorOutputType for compile-time polymorphism.
  /// This lets subclasses define the \c factor() method, without
  /// constraining them to inherit their particular FactorOutputType
  /// from a common abstract base class.  FactorOutputType is meant to
  /// be either just a simple composition of std::pair and
  /// std::vector, or a simple struct.  Its contents are specific to
  /// each intranode TSQR implementation.  and are not intended to be
  /// polymorphic, so it would not make sense for all the different
  /// FactorOutputType types to inherit from a common base class.
  ///
  /// Templating on FactorOutputType means that we can't use run-time
  /// polymorphism to swap between NodeTsqr subclasses, since the
  /// latter are really subclasses of different NodeTsqr
  /// instantiations (i.e., different FactorOutputType types).
  /// However, inheriting from different specializations of NodeTsqr
  /// does enforce correct compile-time polymorphism in a syntactic
  /// way.  It also avoids repeated code for common functionality.
  /// Full run-time polymorphism of different NodeTsqr subclasses
  /// would not be useful.  This is because ultimately each subclass
  /// is bound to a Kokkos Node type, and those only use compile-time
  /// polymorphism.
  ///
  template<class Ordinal, class Scalar, class FactorOutputType>
  class NodeTsqr : public Teuchos::Describable {
  public:
    typedef Ordinal ordinal_type;
    typedef Scalar scalar_type;
    typedef FactorOutputType factor_output_type;

    //! Constructor
    NodeTsqr() {}

    //! Virtual destructor, for memory safety of derived classes.
    virtual ~NodeTsqr() {}

    /// \brief Cache size hint (in bytes) used for the factorization.
    ///
    /// This method is deprecated, because the name is misleading.
    /// Please call \c cache_size_hint() instead.
    virtual size_t TEUCHOS_DEPRECATED cache_block_size() const = 0;

    //! Cache size hint (in bytes) used for the factorization.
    virtual size_t cache_size_hint() const = 0;

    /// \brief One-line description of this object.
    ///
    /// This implements \c Teuchos::Describable::description().
    /// Subclasses should override this to provide a more specific
    /// description of their implementation.  Subclasses may also
    /// implement \c Teuchos::Describable::describe(), which for this
    /// class has a simple default implementation that calls
    /// description() with appropriate indenting.
    virtual std::string description () const {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os << "NodeTsqr<Ordinal=" << TypeNameTraits<Ordinal>::name()
	 << ", Scalar=" << TypeNameTraits<Scalar>::name()
	 << ", ...>: Intranode Tall Skinny QR (TSQR), with cache size hint " 
	 << cache_size_hint();
      return os.str();
    }

    /// \brief Compute the QR factorization of A.
    ///
    /// The resulting Q factor is stored implicitly in two parts.  The
    /// first part is stored in place in the A matrix, and thus
    /// overwrites the input matrix.  The second part is stored in the
    /// returned factor_output_type object.  Both parts must be passed
    /// into \c apply() or \c explicit_Q().
    ///
    /// \param nrows [in] Number of rows in the matrix A to factor.
    /// \param ncols [in] Number of columns in the matrix A to factor.
    /// \param A [in/out] On input: the matrix to factor.  It is
    ///   stored either in column-major order with leading dimension
    ///   (a.k.a. stride) lda, or with contiguous cache blocks (if
    ///   contiguousCacheBlocks is true) according to the prevailing
    ///   cache blocking strategy.  Use the \c cache_block() method to
    ///   convert a matrix in column-major order to the latter format,
    ///   and the \c un_cache_block() method to convert it back.  On
    ///   output: part of the implicit representation of the Q factor.
    ///   (The returned object is the other part of that
    ///   representation.)
    /// \param lda [in] Leading dimension (a.k.a. stride) of the
    ///   matrix A to factor.
    /// \param R [out] The ncols x ncols R factor.
    /// \param ldr [in] leading dimension (a.k.a. stride) of the R 
    ///   factor.
    /// \param contiguousCacheBlocks [in] Whether the cache blocks of
    ///   A are stored contiguously.  If you don't know what this
    ///   means, put "false" here.
    ///
    /// \return Part of the implicit representation of the Q factor.
    ///   The other part is the A matrix on output.
    virtual factor_output_type
    factor (const Ordinal nrows,
	    const Ordinal ncols, 
	    Scalar A[],
	    const Ordinal lda,
	    Scalar R[],
	    const Ordinal ldr,
	    const bool contiguousCacheBlocks) const = 0;

    /// \brief Apply the implicit Q factor from \c factor() to C.
    ///
    /// \param applyType [in] Whether to apply Q, Q^T, or Q^H to C.
    /// \param nrows [in] Number of rows in Q and C.
    /// \param ncols [in] Number of columns in in Q.
    /// \param Q [in] Part of the implicit representation of the Q
    ///   factor; the A matrix output of \c factor().  See the \c
    ///   factor() documentation for details.
    /// \param ldq [in] Leading dimension (a.k.a. stride) of Q, if Q
    ///   is stored in column-major order (not contiguously cache
    ///   blocked).
    /// \param factorOutput [in] Return value of factor(),
    ///   corresponding to Q.
    /// \param ncols_C [in] Number of columns in the matrix C.  This
    ///   may be different than the number of columns in Q.  There is
    ///   no restriction on this value, but we optimize performance
    ///   for the case ncols_C == ncols_Q.
    /// \param C [in/out] On input: Matrix to which to apply the Q
    ///   factor.  On output: Result of applying the Q factor (or Q^T,
    ///   or Q^H, depending on applyType) to C.
    /// \param ldc [in] leading dimension (a.k.a. stride) of C, if C
    ///   is stored in column-major order (not contiguously cache
    ///   blocked).
    /// \param contiguousCacheBlocks [in] Whether the cache blocks of
    ///   Q and C are stored contiguously.  If you don't know what
    ///   this means, put "false" here.
    virtual void
    apply (const ApplyType& applyType,
	   const Ordinal nrows,
	   const Ordinal ncols_Q,
	   const Scalar Q[],
	   const Ordinal ldq,
	   const FactorOutputType& factorOutput,
	   const Ordinal ncols_C,
	   Scalar C[],
	   const Ordinal ldc,
	   const bool contiguousCacheBlocks) const = 0;

    /// \brief Compute the explicit Q factor from the result of \c factor().
    ///
    /// This is equivalent to calling \c apply() on the first ncols_C
    /// columns of the identity matrix (suitably cache-blocked, if
    /// applicable).
    ///
    /// \param nrows [in] Number of rows in Q and C.
    /// \param ncols [in] Number of columns in in Q.
    /// \param Q [in] Part of the implicit representation of the Q
    ///   factor; the A matrix output of \c factor().  See the \c
    ///   factor() documentation for details.
    /// \param ldq [in] Leading dimension (a.k.a. stride) of Q, if Q
    ///   is stored in column-major order (not contiguously cache
    ///   blocked).
    /// \param factorOutput [in] Return value of factor(),
    ///   corresponding to Q.
    /// \param ncols_C [in] Number of columns in the matrix C.  This
    ///   may be different than the number of columns in Q, in which
    ///   case that number of columns of the Q factor will be
    ///   computed.  There is no restriction on this value, but we
    ///   optimize performance for the case ncols_C == ncols_Q.
    /// \param C [out] The first ncols_C columns of the Q factor.
    /// \param ldc [in] leading dimension (a.k.a. stride) of C, if C
    ///   is stored in column-major order (not contiguously cache
    ///   blocked).
    /// \param contiguousCacheBlocks [in] Whether the cache blocks of
    ///   Q and C are stored contiguously.  If you don't know what
    ///   this means, put "false" here.
    virtual void
    explicit_Q (const Ordinal nrows,
		const Ordinal ncols_Q,
		const Scalar Q[],
		const Ordinal ldq,
		const factor_output_type& factorOutput,
		const Ordinal ncols_C,
		Scalar C[],
		const Ordinal ldc,
		const bool contiguousCacheBlocks) const = 0;

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
    virtual void
    Q_times_B (const Ordinal nrows,
	       const Ordinal ncols,
	       Scalar Q[],
	       const Ordinal ldq,
	       const Scalar B[],
	       const Ordinal ldb,
	       const bool contiguousCacheBlocks) const = 0;

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
    virtual void
    fill_with_zeros (const Ordinal nrows,
		     const Ordinal ncols,
		     Scalar A[],
		     const Ordinal lda, 
		     const bool contiguousCacheBlocks) const = 0;
    
  protected:

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
    virtual ConstMatView<Ordinal, Scalar>
    const_top_block (const ConstMatView<Ordinal, Scalar>& C,
		     const bool contiguousCacheBlocks) const = 0;

  public:

    /// \brief Return view of topmost cache block of C.
    ///
    /// \param C [in] View of a matrix C.
    /// \param contiguousCacheBlocks [in] Whether the cache blocks
    ///   in C are stored contiguously.
    ///
    /// Return a view of the topmost cache block (on this node) of the
    /// given matrix C.  This is not necessarily square, though it
    /// must have at least as many rows as columns.  For a view of the
    /// first C.ncols() rows of that block, which methods like
    /// Tsqr::apply() need, do the following:
    /// \code 
    /// MatrixViewType top = this->top_block (C, contig);
    /// MatView<Ordinal, Scalar> square (ncols, ncols, top.get(), top.lda());
    /// \endcode
    ///
    /// Models for MatrixViewType are \c MatView and \c ConstMatView.
    /// MatrixViewType must have member functions nrows(), ncols(),
    /// get(), and lda(), and its constructor must take the same four
    /// arguments as the constructor of \c ConstMatView.
    template<class MatrixViewType>
    MatrixViewType
    top_block (const MatrixViewType& C, 
	       const bool contiguous_cache_blocks) const 
    {
      // The *_top_block() methods don't actually modify the data, so
      // it's safe to handle the matrix's data as const within this
      // method.  The only cast from const to nonconst may be in the
      // return value, but there it's legitimate since we're just
      // using the same constness as C has.
      ConstMatView<Ordinal, Scalar> C_view (C.nrows(), C.ncols(), 
					    C.get(), C.lda());
      ConstMatView<Ordinal, Scalar> C_top = 
	const_top_block (C_view, contiguous_cache_blocks);
      TEST_FOR_EXCEPTION(C_top.nrows() < C_top.ncols(), std::logic_error,
			 "The subclass of NodeTsqr has a bug in const_top_block"
			 "(); it returned a block with fewer rows than columns "
			 "(" << C_top.nrows() << " rows and " << C_top.ncols() 
			 << " columns).  Please report this bug to the Kokkos "
			 "developers.");
      typedef typename MatrixViewType::pointer_type ptr_type;
      return MatrixViewType (C_top.nrows(), C_top.ncols(), 
			     const_cast<ptr_type> (C_top.get()), 
			     C_top.lda());
    }

    /// \brief Does factor() compute R with nonnegative diagonal?
    ///
    /// When using a QR factorization to orthogonalize a block of
    /// vectors, computing an R factor with nonnegative diagonal
    /// ensures that in exact arithmetic, the result of the
    /// orthogonalization (orthogonalized vectors Q and their
    /// coefficients R) are the same as would be produced by
    /// Gram-Schmidt orthogonalization.
    ///
    /// This distinction is important because LAPACK's QR
    /// factorization (_GEQRF) may (and does, in practice) compute an
    /// R factor with negative diagonal entries.
    virtual bool 
    QR_produces_R_factor_with_nonnegative_diagonal () const = 0;

    /// \brief Reveal rank of TSQR's R factor.
    ///
    /// Compute the singular value decomposition (SVD) \f$R = U \Sigma
    /// V^*\f$.  This is done not in place, so that the original R is
    /// not affected.  Use the resulting singular values to compute
    /// the numerical rank of R, with respect to the relative
    /// tolerance tol.  If R is full rank, return without modifying R.
    /// If R is not full rank, overwrite R with \f$\Sigma \cdot
    /// V^*\f$.
    ///
    /// \param ncols [in] Number of (rows and) columns in R.
    /// \param R [in/out] ncols x ncols upper triangular matrix,
    ///   stored in column-major order with leading dimension ldr.
    /// \param ldr [in] Leading dimension of the matrix R.
    /// \param U [out] Left singular vectors of the matrix R; 
    ///   an ncols x ncols matrix with leading dimension ldu.
    /// \param ldu [in] Leading dimension of the matrix U.
    /// \param tol [in] Numerical rank tolerance; relative to 
    ///   the largest nonzero singular value of R.
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
		 const bool contiguousCacheBlocks) const;
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
      // LWORK = -1 indicates a workspace query.  If Scalar is
      // complex, LAPACK had better return something with a zero
      // imaginary part, since I can't allocate imaginary amounts of
      // memory!  (Take the real part to avoid rounding error, since
      // magnitude() may be implemented using a square root always...)
      svd_lwork = static_cast<Ordinal> (Teuchos::ScalarTraits<Scalar>::real (svd_lwork_scalar));

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
    Ordinal rank = ncols; // "innocent unless proven guilty"
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


#endif // __TSQR_Tsqr_NodeTsqr_hpp
