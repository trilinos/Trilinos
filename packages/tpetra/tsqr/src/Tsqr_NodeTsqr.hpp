// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tsqr_NodeTsqr.hpp
/// \brief Common interface and functionality for intranode TSQR.

#ifndef __TSQR_Tsqr_NodeTsqr_hpp
#define __TSQR_Tsqr_NodeTsqr_hpp

#include "Tsqr_ApplyType.hpp"
#include "Tsqr_Matrix.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_Describable.hpp"
#include "Tsqr_Impl_Lapack.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include <vector>

namespace TSQR {
  namespace Impl {
    template<class Ordinal, class Scalar>
    class NodeFactorOutput {
    public:
      virtual ~NodeFactorOutput() = default;
    };
  } // namespace Impl

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
  template<class Ordinal, class Scalar>
  class NodeTsqr : public Teuchos::Describable {
  public:
    using ordinal_type = Ordinal;
    using scalar_type = Scalar;
    using magnitude_type =
      typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
    using factor_output_type = Impl::NodeFactorOutput<Ordinal, Scalar>;
    using mat_view_type = MatView<Ordinal, Scalar>;
    using const_mat_view_type = MatView<Ordinal, const Scalar>;

    //! Constructor
    NodeTsqr() = default;

    //! Virtual destructor, for memory safety of derived classes.
    virtual ~NodeTsqr() = default;

    //! List of valid parameters for the NodeTsqr subclass.
    virtual Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters () const = 0;

    //! Validate and read in parameters.
    virtual void
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& paramList) = 0;

    //! Whether the subclass wants large arrays as GPU device memory.
    virtual bool wants_device_memory () const { return false; }

    /// \brief Whether this object is ready to perform computations.
    ///
    /// Some NodeTsqr subclasses require additional initialization
    /// after construction, before they can perform computations.
    /// Call this method to make sure that the subclass instance is
    /// fully initialized, before calling any of its computational
    /// methods.
    virtual bool ready() const = 0;

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
    virtual Teuchos::RCP<factor_output_type>
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
           const factor_output_type& factorOutput,
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

    /// \brief Force the diagonal entries of the R factor to be
    ///   nonnegative, and change the columns of Q (result of
    ///   explicit_Q) to match (if needed).
    virtual void
    force_nonnegative_diagonal (const Ordinal nrows,
                                const Ordinal ncols,
                                Scalar Q[],
                                const Ordinal ldq,
                                Scalar R[],
                                const Ordinal ldr) const
    {
      mat_view_type Q_view (nrows, ncols, Q, ldq);
      mat_view_type R_view (ncols, ncols, R, ldr);

      // The complex-arithmetic specialization does nothing, since
      // _GEQR{2,F} for complex arithmetic returns an R factor with
      // nonnegative diagonal already.  However, we need the code to
      // compile regardless.
      using STS = Teuchos::ScalarTraits<Scalar>;
      if (! STS::isComplex) {
        using mag_type = typename STS::magnitudeType;
        constexpr mag_type ZERO {};

        for (Ordinal k = 0; k < ncols; ++k) {
          if (STS::real (R_view(k,k)) < ZERO) {
            // Scale column k of Q_view.
            Scalar* const Q_k = &Q_view(0,k);
            for (Ordinal i = 0; i < nrows; ++i) {
              Q_k[i] = -Q_k[i];
            }
            // Scale row k of R_view.  R_view is upper triangular,
            // so we only have to scale right of (and including) the
            // diagonal entry.
            for (int j = k; j < ncols; ++j) {
              R_view(k,j) = -R_view(k,j);
            }
          }
        }
      }
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
    /// \param C [in] Matrix (view), supporting the usual extent(0),
    ///   extent(1), data(), stride(1) interface.
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
    /// mat_view_type square (ncols, ncols, top.data(), top.stride(1));
    /// \endcode
    virtual const_mat_view_type
    const_top_block (const const_mat_view_type& C,
                     const bool /* contiguousCacheBlocks */) const {
      return C;
    }

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
    /// first C.extent(1) rows of that block, which methods like
    /// Tsqr::apply() need, do the following:
    /// \code
    /// MatrixViewType top = this->top_block (C, contig);
    /// mat_view_type square (ncols, ncols, top.data(), top.stride(1));
    /// \endcode
    ///
    /// A model for MatrixViewType is MatView.  MatrixViewType must
    /// have member functions extent(0), extent(1), data(), and
    /// stride(1), and its constructor must take the same four
    /// arguments as the constructor of MatView.
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
      const_mat_view_type C_view (C.extent(0), C.extent(1), C.data(), C.stride(1));
      const_mat_view_type C_top =
        const_top_block (C_view, contiguous_cache_blocks);
      TEUCHOS_TEST_FOR_EXCEPTION(C_top.extent(0) < C_top.extent(1), std::logic_error,
                         "The subclass of NodeTsqr has a bug in const_top_block"
                         "(); it returned a block with fewer rows than columns "
                         "(" << C_top.extent(0) << " rows and " << C_top.extent(1)
                         << " columns).  Please report this bug to the Kokkos "
                         "developers.");
      using pointer = typename MatrixViewType::pointer;
      return MatrixViewType (C_top.extent(0), C_top.extent(1),
                             const_cast<pointer> (C_top.data()),
                             C_top.stride(1));
    }

    /// \brief Copy from "native" NodeTsqr device storage, to a packed
    ///   host matrix.
    virtual Matrix<Ordinal, Scalar>
    copy_to_host (const MatView<Ordinal, Scalar>& C) const
    {
      // FIXME (mfh 17 Dec 2019) Need to reimplement in
      // CuSolverNodeTsqr, since C is device memory there.
      //
      // The same concerns as in CuSolverNodeTsqr::extract_R, about
      // Kokkos::deep_copy not wanting to copy from noncontiguous
      // device memory to contiguous host memory, apply here.
      return Matrix<Ordinal, Scalar> (C);
    }

    /// \brief Copy from a host matrix, to "native" NodeTsqr device
    ///   storage.
    virtual void
    copy_from_host (const MatView<Ordinal, Scalar>& C_device,
                    const MatView<Ordinal, const Scalar>& C_host) const
    {
      // FIXME (mfh 17 Dec 2019) Need to reimplement in
      // CuSolverNodeTsqr, since C_device is device memory there.
      //
      // The same concerns as in CuSolverNodeTsqr::extract_R, about
      // Kokkos::deep_copy not wanting to copy between noncontiguous
      // device memory and contiguous host memory, apply here.
      deep_copy (C_device, C_host);
    }

    //! Set the first C.extent(1) diagonal entries of C to 1.0.
    virtual void
    set_diagonal_entries_to_one
      (const MatView<Ordinal, Scalar>& C) const
    {
      // NOTE (mfh 17 Dec 2019) Downstream classes must reimplement
      // this if C is device memory for those classes.  See
      // wants_device_memory above.
      const Ordinal ncols = C.extent (1);
      for (Ordinal j = 0; j < ncols; ++j) {
        C(j,j) = Scalar (1.0);
      }
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


  template<class Ordinal, class Scalar>
  Ordinal
  NodeTsqr<Ordinal, Scalar>::
  reveal_R_rank (const Ordinal ncols,
                 Scalar R[],
                 const Ordinal ldr,
                 Scalar U[],
                 const Ordinal ldu,
                 const typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol) const
  {
    using Teuchos::as;
    using Teuchos::TypeNameTraits;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    TEUCHOS_TEST_FOR_EXCEPTION(tol < 0, std::invalid_argument,
                       "In NodeTsqr::reveal_R_rank: numerical rank tolerance "
                       "(tol = " << tol << ") is negative.");
    TEUCHOS_TEST_FOR_EXCEPTION(ncols < 0, std::invalid_argument,
                       "In NodeTsqr::reveal_R_rank: number of columns "
                       "(ncols = " << ncols << ") is negative.");
    TEUCHOS_TEST_FOR_EXCEPTION(ldr < ncols, std::invalid_argument,
                       "In NodeTsqr::reveal_R_ank: stride of R (ldr = "
                       << ldr << ") is less than the number of columns "
                       "(ncols = " << ncols << ").");
    TEUCHOS_TEST_FOR_EXCEPTION(ldu < ncols, std::invalid_argument,
                       "In NodeTsqr::reveal_R_rank: stride of U (ldu = "
                       << ldu << ") is less than the number of columns "
                       "(ncols = " << ncols << ")");

    // Zero columns always means rank zero.
    if (ncols == 0) {
      return 0;
    }
    //
    // Compute the SVD (singular value decomposition) of the R
    // factor, using LAPACK's GESVD routine.  We do so in a deep
    // copy (B) because LAPACK overwrites the input.  If the R
    // factor is full rank (expected to be the common case), we need
    // to leave it alone (so that it stays upper triangular).
    //
    Impl::Lapack<Scalar> lapack;
    mat_view_type R_view (ncols, ncols, R, ldr);
    Matrix<Ordinal, Scalar> B (R_view); // B := R (deep copy)
    mat_view_type U_view (ncols, ncols, U, ldu);
    Matrix<Ordinal, Scalar> VT (ncols, ncols, Scalar(0));

    // Set up workspace for the SVD.
    std::vector<magnitude_type> svd_rwork (5*ncols);
    std::vector<magnitude_type> singular_values (ncols);
    Ordinal svd_lwork = -1; // -1 for LWORK query; will be changed

    // LAPACK workspace ("LWORK") query for SVD.  The workspace
    // ("WORK") array is always of Scalar type, even in the complex
    // case.
    {
      // Exception messages in this scope all start with this.
      const char prefix[] = "In NodeTsqr::reveal_R_rank: LAPACK SVD "
        "(_GESVD) workspace query returned ";
      // std::logic_error messages in this scope all end with this.
      const char postfix[] = ".  Please report this bug to the Kokkos "
        "developers.";

      Scalar svd_lwork_scalar {};
      lapack.GESVD ('A', 'A', ncols, ncols, B.data(), B.stride(1),
                    singular_values.data(), U_view.data(), U_view.stride(1),
                    VT.data(), VT.stride(1), &svd_lwork_scalar, svd_lwork,
                    svd_rwork.data());
      // LAPACK returns the workspace array length as a Scalar.  We
      // have to convert it back to an Ordinal in order to allocate
      // the workspace array and pass it in to LAPACK as the LWORK
      // argument.  Ordinal definitely must be a signed type, since
      // LWORK = -1 indicates a workspace query.  If Scalar is
      // complex, LAPACK had better return something with a zero
      // imaginary part, since I can't allocate imaginary amounts of
      // memory!  (Take the real part to avoid rounding error, since
      // magnitude() may be implemented using a square root always...)
      svd_lwork = as<Ordinal> (STS::real (svd_lwork_scalar));

      // LAPACK should always return an LWORK that fits in Ordinal,
      // but it's a good idea to check anyway.  We do so by checking
      // whether casting back from Ordinal to Scalar gives the same
      // original Scalar result.  This should work unless Scalar and
      // Ordinal are user-defined types with weird definitions of
      // the type casts.
      TEUCHOS_TEST_FOR_EXCEPTION
        (as<Scalar> (svd_lwork) != svd_lwork_scalar, std::logic_error,
         prefix << "a workspace array length (LWORK) of type Scalar="
         << TypeNameTraits<Scalar>::name() << " that does not fit in "
         << "Ordinal=" << TypeNameTraits<Ordinal>::name() << " type."
         "  As a Scalar, LWORK=" << svd_lwork_scalar << ", but cast "
         << "to Ordinal, LWORK=" << svd_lwork << postfix);
      // Make sure svd_lwork is nonnegative.  (Ordinal must be a
      // signed type, as we explain above, so this test should never
      // signal any unsigned-to-signed conversions from the compiler.
      // If it does, you're probably using the wrong Ordinal type.
      TEUCHOS_TEST_FOR_EXCEPTION
        (svd_lwork < 0, std::logic_error, prefix << "a negative "
         "workspace array length (LWORK) = " << svd_lwork << postfix);
    }
    // Allocate workspace for LAPACK's SVD routine.
    std::vector<Scalar> svd_work (svd_lwork);

    // Compute SVD $B := U \Sigma V^*$.  B is overwritten, which is
    // why we copied R into B (so that we don't overwrite R if R is
    // full rank).
    lapack.GESVD ('A', 'A', ncols, ncols, B.data(), B.stride(1),
                  singular_values.data(), U_view.data(), U_view.stride(1),
                  VT.data(), VT.stride(1), svd_work.data(), svd_lwork,
                  svd_rwork.data());
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
          singular_values[k] < absolute_tolerance) {
        rank = k;
        break;
      }
    // Don't modify Q or R, if R is full rank.
    if (rank < ncols) { // R is not full rank.
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
      for (Ordinal j = 0; j < ncols; ++j) {
        const Scalar* const VT_j = &VT(0,j);
        Scalar* const R_j = &R_view(0,j);

        for (Ordinal i = 0; i < ncols; ++i) {
          R_j[i] = singular_values[i] * VT_j[i];
        }
      }
    }
    return rank;
  }

  template<class Ordinal, class Scalar>
  Ordinal
  NodeTsqr<Ordinal, Scalar>::
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
    if (ncols == 0) {
      return 0;
    }
    // Matrix to hold the left singular vectors of the R factor.
    Matrix<Ordinal, Scalar> U (ncols, ncols, Scalar {});
    // Compute numerical rank of the R factor using the SVD.
    // Store the left singular vectors in U.
    const Ordinal rank =
      reveal_R_rank (ncols, R, ldr, U.data(), U.stride(1), tol);

    // If R is full rank, we're done.  Otherwise, reveal_R_rank()
    // already computed the SVD \f$R = U \Sigma V^*\f$ of (the
    // input) R, and overwrote R with \f$\Sigma V^*\f$.  Now, we
    // compute \f$Q := Q \cdot U\f$, respecting cache blocks of Q.
    if (rank < ncols) {
      Q_times_B (nrows, ncols, Q, ldq, U.data(), U.stride(1),
                 contiguousCacheBlocks);
    }
    return rank;
  }

} // namespace TSQR


#endif // __TSQR_Tsqr_NodeTsqr_hpp
