// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tsqr_SequentialTsqr.hpp
/// \brief Implementation of the sequential cache-blocked part of TSQR.

#ifndef TSQR_SEQUENTIALTSQR_HPP
#define TSQR_SEQUENTIALTSQR_HPP

#include "Tsqr_ApplyType.hpp"
#include "Tsqr_Matrix.hpp"
#include "Tsqr_CacheBlockingStrategy.hpp"
#include "Tsqr_CacheBlocker.hpp"
#include "Tsqr_Impl_CombineUser.hpp"
#include "Tsqr_NodeTsqr.hpp"
#include "Tsqr_Util.hpp"
#include "Tsqr_Impl_SystemBlas.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <algorithm>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace TSQR {
  namespace Impl {
    template<class LocalOrdinal, class Scalar>
    class SequentialTsqrFactorOutput :
      public NodeFactorOutput<LocalOrdinal, Scalar>
    {
    private:
      using my_data_type = std::vector<std::vector<Scalar>>;
    public:
      SequentialTsqrFactorOutput () = default;
      ~SequentialTsqrFactorOutput () override = default;

      void add_and_consume (std::vector<Scalar>&& tau) {
        data_.emplace_back (tau);
      }
      typename my_data_type::const_iterator begin() const {
        return data_.begin();
      }
      typename my_data_type::const_reverse_iterator rbegin() const {
        return data_.rbegin();
      }
    private:
      my_data_type data_;
    };
  } // namespace Impl

  /// \class SequentialTsqr
  /// \brief Sequential cache-blocked TSQR factorization.
  /// \author Mark Hoemmen
  ///
  /// TSQR (Tall Skinny QR) is a collection of different algorithms
  /// for computing the QR factorization of a "tall and skinny" matrix
  /// (with many more rows than columns).  We use it in Trilinos as an
  /// orthogonalization method for Epetra_MultiVector and
  /// Tpetra::MultiVector.  (In this context, TSQR is provided as an
  /// "OrthoManager" in Anasazi and Belos; you do not have to use it
  /// directly.)  For details, see e.g., our 2008 University of
  /// California Berkeley technical report (Demmel, Grigori, Hoemmen,
  /// and Langou), or our Supercomputing 2009 paper (Demmel, Hoemmen,
  /// Mohiyuddin, and Yelick).
  ///
  /// SequentialTsqr implements the "sequential TSQR" algorithm of the
  /// aforementioned 2008 technical report.  It breaks up the matrix
  /// by rows into "cache blocks," and iterates over consecutive cache
  /// blocks.  The input matrix may be in either the conventional
  /// LAPACK-style column-major layout, or in a "cache-blocked"
  /// layout.  We provide conversion routines between these two
  /// formats.  Users should not attempt to construct a matrix in the
  /// latter format themselves.  In our experience, the performance
  /// difference between the two formats is not significant, but this
  /// may be different on different architectures.
  ///
  /// SequentialTsqr is designed to be used as the "intranode TSQR"
  /// part of the full TSQR implementation in Tsqr.  The Tsqr class
  /// can use any of various intranode TSQR implementations.
  /// SequentialTsqr is an appropriate choice when running in MPI-only
  /// mode.  Other intranode TSQR implementations, such as TbbTsqr
  /// (which has been removed temporarily) are appropriate for hybrid
  /// parallelism (MPI + threads).
  ///
  /// SequentialTsqr is unlikely to benefit from a multithreaded BLAS
  /// implementation.  In fact, implementations of LAPACK's QR
  /// factorization generally do not show performance benefits from
  /// multithreading when factoring tall skinny matrices.  (See our
  /// Supercomputing 2009 paper and my IPDPS 2011 paper.)  This is why
  /// we built other intranode TSQR factorizations that do effectively
  /// exploit thread-level parallelism, such as TbbTsqr.
  ///
  /// \note To implementers: SequentialTsqr cannot currently be a
  ///   Teuchos::ParameterListAcceptorDefaultBase, because the latter
  ///   uses RCP, and RCPs (more specifically, their reference counts)
  ///   are not currently thread safe.  TbbTsqr uses SequentialTsqr in
  ///   parallel to implement each thread's cache-blocked TSQR.  This
  ///   can be fixed as soon as RCPs are made thread safe.
  template<class LocalOrdinal, class Scalar>
  class SequentialTsqr :
    public NodeTsqr<LocalOrdinal, Scalar>,
    private Impl::CombineUser<LocalOrdinal, Scalar>
  {
  private:
    using base_type = NodeTsqr<LocalOrdinal, Scalar>;
    using my_factor_output_type =
      Impl::SequentialTsqrFactorOutput<LocalOrdinal, Scalar>;

  public:
    using ordinal_type = typename base_type::ordinal_type;
    using scalar_type = typename base_type::scalar_type;
    using mat_view_type = typename base_type::mat_view_type;
    using const_mat_view_type =
      typename base_type::const_mat_view_type;
    using magnitude_type = typename base_type::magnitude_type;
    using factor_output_type = typename base_type::factor_output_type;

  private:
    Combine<ordinal_type, scalar_type>&
    getMyCombine (const ordinal_type /* maxNumCols */) const
    {
      // FIXME (mfh 20 Dec 2019) If SequentialTsqr has more than one
      // cache block, it only passes tests if you use CombineNative.
      // This likely explains why it fails with complex Scalar types,
      // since CombineNative just uses CombineDefault in that case.  I
      // tried making SequentialTsqr's implementation of
      // QR_produces_R_factor_with_nonnegative_diagonal always return
      // false, but that didn't help, so the issue likely is
      // CombineDefault.
      return this->getCombine ("CombineNative");
    }

    /// \brief Factor the first cache block of the matrix.
    ///
    /// Compute the QR factorization of the first cache block A_top.
    /// Overwrite the upper triangle of A_top with the R factor, and
    /// return a view of the R factor (stored in place in A_top).
    /// Overwrite the (strict) lower triangle of A_top, and the
    /// A_top.extent(1) entries of tau, with an implicit representation
    /// of the Q factor.
    ///
    /// \param combine [in/out] Implementation of linear algebra
    ///   primitives.  This is non-const because some Combine
    ///   implementations use scratch space.
    ///
    /// \param A_top [in/out] On input: the first (topmost) cache
    ///   block of the matrix.  Prerequisite: A_top.extent(0) >=
    ///   A.top.extent(1).  On output, the upper triangle of A_top is
    ///   overwritten with the R factor, and the lower trapezoid of
    ///   A_top is overwritten with part of the implicit
    ///   representation of the Q factor.
    ///
    /// \param tau [out] Array of length >= A_top.extent(1).  On output:
    ///   the TAU array (see the LAPACK documentation for _GEQRF).
    ///
    /// \param work [out] Workspace array of length >= A_top.extent(1).
    ///
    /// \return A view of the upper triangle of A_top, containing the
    ///   R factor.
    mat_view_type
    factor_first_block (Combine<LocalOrdinal, Scalar>& combine,
                        const mat_view_type& A_top,
                        std::vector<Scalar>& tau,
                        Scalar work[],
                        const LocalOrdinal lwork) const
    {
      combine.factor_first (A_top, tau.data (), work, lwork);
      const LocalOrdinal ncols = A_top.extent (1);
      return partition_2x1 (A_top, ncols).first;
    }

  public:
    /// \brief The standard constructor.
    ///
    /// \param cacheSizeHint [in] Cache size hint in bytes to use in
    ///   the sequential TSQR factorization.  If 0, the implementation
    ///   will pick a reasonable size.  Good nondefault choices are
    ///   the amount of per-CPU highest-level private cache, or the
    ///   amount of lowest-level shared cache divided by the number of
    ///   CPU cores sharing it.  We recommend experimenting to find
    ///   the best value.  Too large a value is worse than too small a
    ///   value, though an excessively small value will result in
    ///   extra computation and may also cause a slow down.
    ///
    /// \param sizeOfScalar [in] The number of bytes required to store
    ///   a Scalar value.  This is used to compute the dimensions of
    ///   cache blocks.  If sizeof(Scalar) correctly reports the size
    ///   of the representation of Scalar in memory, you can use the
    ///   default.  The default is correct for float, double, and any
    ///   of various fixed-length structs (like double-double and
    ///   quad-double).  It should also work for std::complex<T> where
    ///   T is anything in the previous sentence's list.  It does
    ///   <it>not</it> work for arbitrary-precision types whose
    ///   storage is dynamically allocated, even if the amount of
    ///   storage is a constant.  In the latter case, you should
    ///   specify a nondefault value.
    ///
    /// \note sizeOfScalar affects performance, not correctness (more
    ///   or less -- it should never be zero, for example).  It's OK
    ///   for it to be a slight overestimate.  Being much too big may
    ///   affect performance by underutilizing the cache.  Being too
    ///   small may also affect performance by thrashing the cache.
    ///
    /// \note If Scalar is an arbitrary-precision type whose
    ///   representation length can change at runtime, you should
    ///   construct a new SequentialTsqr object whenever the
    ///   representation length changes.
    SequentialTsqr (const size_t cacheSizeHint = 0,
                    const size_t sizeOfScalar = sizeof(Scalar)) :
      strategy_ (cacheSizeHint, sizeOfScalar)
    {}

    /// \brief Alternate constructor for a given cache blocking strategy.
    ///
    /// The cache blocking strategy stores the same information as
    /// would be passed into the standard constructor: the cache block
    /// size, and the size of the Scalar type.
    ///
    /// \param strategy [in] Cache blocking strategy to use (copied).
    ///
    SequentialTsqr (const CacheBlockingStrategy<LocalOrdinal, Scalar>& strategy) :
      strategy_ (strategy)
    {}

    /// \brief Alternate constructor that takes a list of parameters.
    ///
    /// See the documentation of \c setParameterList() for the list of
    /// currently understood parameters.  The constructor ignores
    /// parameters that it doesn't understand.
    ///
    /// \param plist [in/out] On input: List of parameters.  On
    ///   output: Missing parameters are filled in with default
    ///   values.
    SequentialTsqr (const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      setParameterList (params);
    }

    /// \brief List of valid parameters for SequentialTsqr.
    ///
    /// \note This object has to create a new parameter list each
    ///   time, since it cannot cache an RCP (due to thread safety --
    ///   TbbTsqr invokes multiple instances of SequentialTsqr in
    ///   parallel).
    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters () const override
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      const size_t cacheSizeHint = 0;
      const size_t sizeOfScalar = sizeof(Scalar);

      RCP<ParameterList> plist = parameterList ("NodeTsqr");
      plist->set ("Cache Size Hint", cacheSizeHint,
                  "Cache size hint in bytes (as a size_t) to use for intranode"
                  "TSQR.  If zero, TSQR will pick a reasonable default.  "
                  "The size should correspond to that of the largest cache that "
                  "is private to each CPU core, if such a private cache exists; "
                  "otherwise, it should correspond to the amount of shared "
                  "cache, divided by the number of cores sharing that cache.");
      plist->set ("Size of Scalar", sizeOfScalar, "Size of the Scalar type.  "
                  "Default is sizeof(Scalar).  Only set if sizeof(Scalar) does "
                  "not describe how much memory a Scalar type takes.");
      return plist;
    }

    /// \brief Set parameters.
    ///
    /// \param plist [in/out] On input: List of parameters.  On
    ///   output: Missing parameters are filled in with default
    ///   values.
    ///
    /// For a list of currently understood parameters, see the
    /// parameter list returned by \c getValidParameters().
    void
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist) override
    {
      using Teuchos::Exceptions::InvalidParameter;
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      RCP<ParameterList> params = plist.is_null() ?
        parameterList (*getValidParameters()) : plist;

      const std::string cacheSizeHintName ("Cache Size Hint");
      const std::string sizeOfScalarName ("Size of Scalar");
      // In order to avoid calling getValidParameters() and
      // constructing a default list, we set missing values here to
      // their defaults.  This duplicates default values set in
      // getValidParameters(), so if you change those, be careful to
      // change them here.
      size_t cacheSizeHint = 0;
      size_t sizeOfScalar = sizeof(Scalar);

      if (params->isType<size_t> (cacheSizeHintName)) {
        cacheSizeHint = params->get<size_t> (cacheSizeHintName);
      }
      else {
        params->set (cacheSizeHintName, cacheSizeHint);
      }

      if (params->isType<size_t> (sizeOfScalarName)) {
        sizeOfScalar = params->get<size_t> (sizeOfScalarName);
      }
      else {
        params->set (sizeOfScalarName, sizeOfScalar);
      }

      // Reconstruct the cache blocking strategy, since we may have
      // changed parameters.
      strategy_ = CacheBlockingStrategy<LocalOrdinal, Scalar> (cacheSizeHint,
                                                               sizeOfScalar);
    }

    /// \brief One-line description of this object.
    ///
    /// This implements Teuchos::Describable::description().  For now,
    /// SequentialTsqr uses the default implementation of
    /// Teuchos::Describable::describe().
    std::string description () const override {
      std::ostringstream os;
      os << "Intranode Tall Skinny QR (TSQR): sequential cache-blocked "
        "implementation with cache size hint " << this->cache_size_hint()
         << " bytes.";
      return os.str();
    }

    //! Whether this object is ready to perform computations.
    bool ready() const override {
      return true;
    }

    //! Whether factor() promises to compute R with a nonnegative diagonal.
    bool
    QR_produces_R_factor_with_nonnegative_diagonal () const override
    {
      // FIXME (19 Dec 2019) If the combine type is dynamic, we can't
      // answer this question without knowing the number of columns.
      // Just guess for now.
      constexpr LocalOrdinal fakeNumCols = 10;
      auto& c = this->getMyCombine (fakeNumCols);
      return c.QR_produces_R_factor_with_nonnegative_diagonal ();
    }

    /// \brief Cache size hint (in bytes) used for the factorization.
    ///
    /// This may be different than the cache size hint argument
    /// specified in the constructor.  SequentialTsqr treats that as a
    /// hint, not a command.
    size_t cache_size_hint () const override {
      return strategy_.cache_size_hint();
    }

    /// \brief Extract R factor from \c factor() results.
    ///
    /// The five-argument version of \c factor() leaves the R factor
    /// in place in the matrix A.  This method copies the R factor out
    /// of A into a separate matrix R in column-major order
    /// (regardless of whether A was stored with contiguous cache
    /// blocks).
    void
    extract_R (const LocalOrdinal nrows,
               const LocalOrdinal ncols,
               const Scalar A[],
               const LocalOrdinal lda,
               Scalar R[],
               const LocalOrdinal ldr,
               const bool contiguous_cache_blocks) const
    {
      const_mat_view_type A_view (nrows, ncols, A, lda);

      // Identify top cache block of A
      const_mat_view_type A_top = this->top_block (A_view, contiguous_cache_blocks);

      // Fill R (including lower triangle) with zeros.
      mat_view_type R_view (ncols, ncols, R, ldr);
      deep_copy (R_view, Scalar {});

      // Copy out the upper triangle of the R factor from A into R.
      copy_upper_triangle (R, A_top);
    }

    /// \brief Compute the QR factorization of the matrix A.
    ///
    /// See the \c NodeTsqr documentation for details.  This version
    /// of factor() is more useful than the five-argument version,
    /// when using SequentialTsqr as the intranode TSQR implementation
    /// in \c Tsqr.  The five-argument version is more useful when
    /// using SequentialTsqr inside of another intranode TSQR
    /// implementation, such as TbbTsqr.
    Teuchos::RCP<factor_output_type>
    factor (const LocalOrdinal nrows,
            const LocalOrdinal ncols,
            Scalar A[],
            const LocalOrdinal lda,
            Scalar R[],
            const LocalOrdinal ldr,
            const bool contigCacheBlocks) const override
    {
      using LO = LocalOrdinal;
      CacheBlocker<LO, Scalar> blocker (nrows, ncols, strategy_);
      auto& combine = this->getMyCombine (ncols);
      const LO lwork = combine.work_size (nrows, ncols, ncols);
      std::vector<Scalar> work (lwork);
      Teuchos::RCP<my_factor_output_type> tau_arrays
        (new my_factor_output_type);

      // We say "A_rest" because it points to the remaining part of
      // the matrix left to factor; at the beginning, the "remaining"
      // part is the whole matrix, but that will change as the
      // algorithm progresses.
      //
      // Note: if the cache blocks are stored contiguously, lda won't
      // be the correct leading dimension of A, but it won't matter:
      // we only ever operate on A_cur here, and A_cur's leading
      // dimension is set correctly by split_top_block.
      mat_view_type A_rest (nrows, ncols, A, lda);
      // This call modifies A_rest.
      mat_view_type A_cur =
        blocker.split_top_block (A_rest, contigCacheBlocks);

      // Factor the topmost block of A.
      std::vector<Scalar> tau_first (ncols);
      mat_view_type R_view =
        factor_first_block (combine, A_cur, tau_first,
                            work.data (), lwork);
      tau_arrays->add_and_consume (std::move (tau_first));

      while (! empty (A_rest)) {
        A_cur = blocker.split_top_block (A_rest, contigCacheBlocks);
        std::vector<Scalar> tau (ncols);
        combine.factor_inner (R_view, A_cur, tau.data (),
                              work.data (), lwork);
        tau_arrays->add_and_consume (std::move (tau));
      }

      // Copy the R factor resulting from the factorization out of
      // R_view (a view of the topmost cache block of A) into the R
      // output argument.
      mat_view_type R_out (ncols, ncols, R, ldr);
      deep_copy (R_out, Scalar {});
      copy_upper_triangle (R_out, R_view);
      return tau_arrays;
    }

    /// \brief The number of cache blocks that factor() would use.
    ///
    /// The \c factor() method breaks the input matrix A into one or
    /// more cache blocks.  This method reports how many cache blocks
    /// \c factor() would use, without actually factoring the matrix.
    ///
    /// \param nrows [in] Number of rows in the matrix A.
    /// \param ncols [in] Number of columns in the matrix A.
    /// \param A [in] The matrix A.  If contiguous_cache_blocks is
    ///   false, A is stored in column-major order; otherwise, A is
    ///   stored with contiguous cache blocks (as the \c cache_block()
    ///   method would do).
    /// \param lda [in] If the matrix A is stored in column-major
    ///   order: the leading dimension (a.k.a. stride) of A.
    ///   Otherwise, the value of this parameter doesn't matter.
    /// \param contigCacheBlocks [in] Whether the cache blocks
    ///   in the matrix A are stored contiguously.
    ///
    /// \return Number of cache blocks in the matrix A: a positive integer.
    LocalOrdinal
    factor_num_cache_blocks (const LocalOrdinal nrows,
                             const LocalOrdinal ncols,
                             const Scalar A[],
                             const LocalOrdinal lda,
                             const bool contigCacheBlocks) const
    {
      using LO = LocalOrdinal;
      CacheBlocker<LO, Scalar> blocker (nrows, ncols, strategy_);
      LO count = 0;

      const_mat_view_type A_rest (nrows, ncols, A, lda);
      if (empty (A_rest)) {
        return count;
      }
      const_mat_view_type A_cur =
        blocker.split_top_block (A_rest, contigCacheBlocks);
      ++count; // first factor step

      while (! empty (A_rest)) {
        A_cur = blocker.split_top_block (A_rest, contigCacheBlocks);
        ++count; // next factor step
      }
      return count;
    }

    /// \brief Apply the implicit Q factor to the matrix C.
    ///
    /// See the \c NodeTsqr documentation for details.
    void
    apply (const ApplyType& apply_type,
           const LocalOrdinal nrows,
           const LocalOrdinal ncols_Q,
           const Scalar Q[],
           const LocalOrdinal ldq,
           const factor_output_type& factor_output,
           const LocalOrdinal ncols_C,
           Scalar C[],
           const LocalOrdinal ldc,
           const bool contigCacheBlocks) const override
    {
      using LO = LocalOrdinal;
      const char prefix[] = "TSQR::SequentialTsqr::apply: ";

      // Quick exit and error tests
      if (ncols_Q == 0 || ncols_C == 0 || nrows == 0) {
        return;
      }
      else if (ldc < nrows) {
        std::ostringstream os;
        os << prefix << "ldc (= " << ldc << ") < nrows (= "
           << nrows << ")";
        throw std::invalid_argument (os.str());
      }
      else if (ldq < nrows) {
        std::ostringstream os;
        os << prefix << "ldq (= " << ldq << ") < nrows (= "
           << nrows << ")";
        throw std::invalid_argument (os.str());
      }

      const my_factor_output_type& tau_arrays = [&] () {
        const my_factor_output_type* tau_arrays_ptr =
          dynamic_cast<const my_factor_output_type*> (&factor_output);
        if (tau_arrays_ptr == nullptr) {
          using Teuchos::demangleName;
          using Teuchos::TypeNameTraits;
          using Teuchos::typeName;
          std::ostringstream os;
          os << prefix << "Input factor_output_type object was not "
            "created by the same type of SequentialTsqr object as "
            "this one.  This object has type " << typeName (*this) <<
            " and its subclass of factor_output_type has type " <<
            TypeNameTraits<my_factor_output_type>::name () << ", but "
            "the input factor_output_type object has dynamic type "
            << demangleName (typeid (factor_output).name ());
          throw std::invalid_argument (os.str ());
        }
        return *tau_arrays_ptr;
      } ();

      // If contiguous cache blocks are used, then we have to use the
      // same convention as we did for factor().  Otherwise, we are
      // free to choose the cache block dimensions as we wish in
      // apply(), independently of what we did in factor().
      CacheBlocker<LO, Scalar> blocker (nrows, ncols_Q, strategy_);
      auto& combine =
        this->getMyCombine (std::max (ncols_Q, ncols_C));
      const LO lwork = combine.work_size (nrows, ncols_Q, ncols_C);
      std::vector<Scalar> work (lwork);

      const bool transposed = apply_type.transposed ();

      // We say "*_rest" because it points to the remaining part of
      // the matrix left to factor; at the beginning, the "remaining"
      // part is the whole matrix, but that will change as the
      // algorithm progresses.
      //
      // Note: if the cache blocks are stored contiguously, ldq
      // resp. ldc won't be the correct leading dimension, but it
      // won't matter, since we only read the leading dimension of
      // return values of split_top_block() / split_bottom_block(),
      // which are set correctly (based e.g., on whether or not we are
      // using contiguous cache blocks).
      const_mat_view_type Q_rest (nrows, ncols_Q, Q, ldq);
      mat_view_type C_rest (nrows, ncols_C, C, ldc);

      // Identify the top ncols_C by ncols_C block of C.  C_rest is
      // not modified.
      mat_view_type C_top =
        blocker.top_block (C_rest, contigCacheBlocks);

      if (transposed) {
        const_mat_view_type Q_cur =
          blocker.split_top_block (Q_rest, contigCacheBlocks);
        mat_view_type C_cur =
          blocker.split_top_block (C_rest, contigCacheBlocks);

        // Apply the topmost block of Q.
        auto tau_iter = tau_arrays.begin();
        const std::vector<Scalar>& tau_first = *tau_iter++;
        combine.apply_first (apply_type, Q_cur, tau_first.data (),
                             C_cur, work.data (), lwork);
        while (! empty (Q_rest)) {
          Q_cur = blocker.split_top_block (Q_rest, contigCacheBlocks);
          C_cur = blocker.split_top_block (C_rest, contigCacheBlocks);
          const Scalar* tau = tau_iter->data ();
          combine.apply_inner (apply_type, Q_cur, tau, C_top, C_cur,
                               work.data (), lwork);
          tau_iter++;
        }
      }
      else {
        // Start with the last local Q factor and work backwards up
        // the matrix.
        auto tau_iter = tau_arrays.rbegin ();
        const_mat_view_type Q_cur =
          blocker.split_bottom_block (Q_rest, contigCacheBlocks);
        mat_view_type C_cur =
          blocker.split_bottom_block (C_rest, contigCacheBlocks);
        while (! empty (Q_rest)) {
          const Scalar* tau = tau_iter->data ();
          combine.apply_inner (apply_type, Q_cur, tau, C_top, C_cur,
                               work.data (), lwork);
          tau_iter++;
          Q_cur =
            blocker.split_bottom_block (Q_rest, contigCacheBlocks);
          C_cur =
            blocker.split_bottom_block (C_rest, contigCacheBlocks);
        }
        // Apply to last (topmost) cache block.
        const std::vector<Scalar>& tau_first = *tau_iter++;
        combine.apply_first (apply_type, Q_cur, tau_first.data (),
                             C_cur, work.data (), lwork);
      }
    }

    /// \brief Compute the explicit Q factor from the result of factor().
    ///
    /// See the \c NodeTsqr documentation for details.
    void
    explicit_Q (const LocalOrdinal nrows,
                const LocalOrdinal ncols_Q,
                const Scalar Q[],
                const LocalOrdinal ldq,
                const factor_output_type& factor_output,
                const LocalOrdinal ncols_C,
                Scalar C[],
                const LocalOrdinal ldc,
                const bool contigCacheBlocks) const override
    {
      mat_view_type C_view (nrows, ncols_C, C, ldc);
      deep_copy (C_view, Scalar {});
      // Don't just call set_diagonal_entries_to_one(C_view), because
      // that doesn't respect contigCacheBlocks.
      auto C_top = this->top_block (C_view, contigCacheBlocks);
      deep_copy (C_top, Scalar {});
      this->set_diagonal_entries_to_one (C_top);
      apply (ApplyType::NoTranspose,
             nrows, ncols_Q, Q, ldq, factor_output,
             ncols_C, C, ldc, contigCacheBlocks);
    }

    /// \brief Compute Q := Q*B.
    ///
    /// See the NodeTsqr documentation for details.
    void
    Q_times_B (const LocalOrdinal nrows,
               const LocalOrdinal ncols,
               Scalar Q[],
               const LocalOrdinal ldq,
               const Scalar B[],
               const LocalOrdinal ldb,
               const bool contigCacheBlocks) const override
    {
      using Teuchos::NO_TRANS;
      using LO = LocalOrdinal;

      // Take the easy exit if available.
      if (ncols == 0 || nrows == 0) {
        return;
      }

      // Compute Q := Q*B by iterating through cache blocks of Q.
      // This iteration works much like iteration through cache blocks
      // of A in factor() (which see).  Here, though, each cache block
      // computation is completely independent of the others; a slight
      // restructuring of this code would parallelize nicely using
      // OpenMP.
      CacheBlocker<LO, Scalar> blocker (nrows, ncols, strategy_);
      Impl::SystemBlas<Scalar> blas;
      mat_view_type Q_rest (nrows, ncols, Q, ldq);
      Matrix<LO, Scalar> Q_cur_copy (0, 0); // will be resized
      while (! empty (Q_rest)) {
        mat_view_type Q_cur =
          blocker.split_top_block (Q_rest, contigCacheBlocks);

        // GEMM doesn't like aliased arguments, so we use a copy.
        // We only copy the current cache block, rather than all of
        // Q; this saves memory.
        Q_cur_copy.reshape (Q_cur.extent (0), ncols);
        deep_copy (Q_cur_copy, Q_cur);
        // Q_cur := Q_cur_copy * B.
        constexpr Scalar ZERO {};
        constexpr Scalar ONE (1.0);
        blas.GEMM (NO_TRANS, NO_TRANS,
                   Q_cur.extent (0), ncols, ncols,
                   ONE, Q_cur_copy.data (), Q_cur_copy.stride (1),
                   B, ldb,
                   ZERO, Q_cur.data (), Q_cur.stride (1));
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
    void
    cache_block (const LocalOrdinal nrows,
                 const LocalOrdinal ncols,
                 Scalar A_out[],
                 const Scalar A_in[],
                 const LocalOrdinal lda_in) const override
    {
      CacheBlocker<LocalOrdinal, Scalar> blocker
        (nrows, ncols, strategy_);
      blocker.cache_block (nrows, ncols, A_out, A_in, lda_in);
    }

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
    void
    un_cache_block (const LocalOrdinal nrows,
                    const LocalOrdinal ncols,
                    Scalar A_out[],
                    const LocalOrdinal lda_out,
                    const Scalar A_in[]) const override
    {
      CacheBlocker<LocalOrdinal, Scalar> blocker
        (nrows, ncols, strategy_);
      blocker.un_cache_block (nrows, ncols, A_out, lda_out, A_in);
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
    /// \param contigCacheBlocks [in] Whether the cache blocks
    ///   in A are stored contiguously.
    void
    fill_with_zeros (const LocalOrdinal nrows,
                     const LocalOrdinal ncols,
                     Scalar A[],
                     const LocalOrdinal lda,
                     const bool contigCacheBlocks) const override
    {
      CacheBlocker<LocalOrdinal, Scalar> blocker
        (nrows, ncols, strategy_);
      blocker.fill_with_zeros (nrows, ncols, A, lda,
                               contigCacheBlocks);
    }

  protected:
    /// \brief Return the topmost cache block of the matrix C.
    ///
    /// NodeTsqr's top_block() method must be implemented using
    /// subclasses' const_top_block() method, since top_block() is a
    /// template method and template methods cannot be virtual.
    ///
    /// \param C [in] View of a matrix, with at least as many rows as
    ///   columns.
    /// \param contigCacheBlocks [in] Whether the cache blocks of C
    ///   are stored contiguously.
    ///
    /// \return View of the topmost cache block of the matrix C.
    const_mat_view_type
    const_top_block (const const_mat_view_type& C,
                     const bool contigCacheBlocks) const override
    {
      // The CacheBlocker object knows how to construct a view of the
      // top cache block of C.  This is complicated because cache
      // blocks (in C) may or may not be stored contiguously.  If they
      // are stored contiguously, the CacheBlocker knows the right
      // layout, based on the cache blocking strategy.
      using blocker_type = CacheBlocker<LocalOrdinal, Scalar>;
      blocker_type blocker (C.extent (0), C.extent (1), strategy_);

      // This is a view of the topmost cache block of C.  C_top_block
      // should have >= ncols rows, otherwise either cache blocking is
      // broken or the input matrix C itself had fewer rows than
      // columns.
      return blocker.top_block (C, contigCacheBlocks);
    }

  private:
    //! Strategy object that helps us cache block matrices.
    CacheBlockingStrategy<LocalOrdinal, Scalar> strategy_;
  };

} // namespace TSQR

#endif // TSQR_SEQUENTIALTSQR_HPP
