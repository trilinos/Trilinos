// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tsqr_DistTsqr.hpp
/// \brief Internode part of TSQR.
///
#ifndef TSQR_DISTTSQR_HPP
#define TSQR_DISTTSQR_HPP

#include "Tsqr_DistTsqrHelper.hpp"
#include "Tsqr_DistTsqrRB.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

#include <utility> // std::pair

namespace TSQR {
  /// \class DistTsqr
  /// \brief Internode part of TSQR.
  /// \author Mark Hoemmen
  ///
  /// \tparam LocalOrdinal Index type for dense matrices of Scalar.
  /// \tparam Scalar Value type for matrices to factor.
  ///
  /// This class combines the square R factors computed by the
  /// intranode TSQR factorization (NodeTsqr subclass) on individual
  /// MPI processes.
  template<class LocalOrdinal, class Scalar>
  class DistTsqr : public Teuchos::ParameterListAcceptorDefaultBase {
  public:
    using scalar_type = Scalar;
    using ordinal_type = LocalOrdinal;

  private:
    using VecVec = std::vector<std::vector<scalar_type>>;

  public:
    using mat_view_type = MatView<ordinal_type, scalar_type>;
    using FactorOutput = std::pair<VecVec, VecVec>;
    using rank_type = int;

    /// \brief Constructor (that accepts a parameter list).
    ///
    /// \param plist [in/out] List of parameters for configuring TSQR.
    ///   The specific parameter keys that are read depend on the TSQR
    ///   implementation.  For details, call \c getValidParameters()
    ///   and examine the documentation embedded therein.
    DistTsqr (const Teuchos::RCP<Teuchos::ParameterList>& plist)
    {
      setParameterList (plist);
    }

    //! Constructor (that uses default parameters).
    DistTsqr ()
    {
      setParameterList (Teuchos::null);
    }

    void
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      RCP<ParameterList> params = plist.is_null() ?
        parameterList (*getValidParameters()) : plist;

      // Do nothing for now, other than store the list.
      this->setMyParamList (params);
    }

    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters() const
    {
      return Teuchos::parameterList ("DistTsqr"); // Empty list for now.
    }

    /// \brief Finish initialization using the messenger object.
    ///
    /// \param messenger [in/out] An object handling communication
    ///   between (MPI) process(es).
    void init (const Teuchos::RCP<MessengerBase<scalar_type> >& messenger)
    {
      messenger_ = messenger;
      reduceBroadcastImpl_ =
        Teuchos::rcp (new DistTsqrRB<ordinal_type, scalar_type> (messenger_));
    }

    /// \brief Rank of this (MPI) process.
    ///
    /// Rank is computed via MPI_Comm_rank() on the underlying
    /// communicator, if the latter is an MPI communicator.  If it's a
    /// serial "communicator," the rank is always zero.
    rank_type rank() const {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! ready (), std::logic_error, "Before using DistTsqr "
         "computational methods, you must first call init() with a "
         "valid MessengerBase instance.");
      return messenger_->rank();
    }

    /// \brief Total number of MPI processes in this communicator.
    ///
    /// The size is communicated via MPI_Comm_size() on the underlying
    /// communicator, if the latter is an MPI communicator.  If it's a
    /// serial "communicator," the size is always one.
    rank_type size() const {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! ready (), std::logic_error, "Before using DistTsqr "
         "computational methods, you must first call init() with a "
         "valid MessengerBase instance.");
      return messenger_->size();
    }

    virtual ~DistTsqr () = default;

    /// \brief Does the R factor have a nonnegative diagonal?
    ///
    /// DistTsqr implements a QR factorization (of a distributed
    /// matrix with a special structure).  Some, but not all, QR
    /// factorizations produce an R factor whose diagonal may include
    /// negative entries.  This Boolean tells you whether DistTsqr
    /// promises to compute an R factor whose diagonal entries are all
    /// nonnegative.
    bool
    QR_produces_R_factor_with_nonnegative_diagonal () const
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! ready (), std::logic_error, "Before using DistTsqr "
         "computational methods, you must first call init() with a "
         "valid MessengerBase instance.");
      TEUCHOS_ASSERT( reduceBroadcastImpl_.getRawPtr () != nullptr );
      return reduceBroadcastImpl_->
        QR_produces_R_factor_with_nonnegative_diagonal ();
    }

    /// \brief Internode TSQR with explicit Q factor.
    ///
    /// Call this routine, instead of \c factor() and \c explicit_Q(),
    /// if you want to compute the QR factorization and only want the
    /// Q factor in explicit form (i.e., as a matrix).
    ///
    /// \param R_mine [in/out] View of a matrix with at least as many
    ///   rows as columns.  On input: upper triangular matrix (R
    ///   factor from intranode TSQR); different on each process..  On
    ///   output: R factor from intranode QR factorization; bitwise
    ///   identical on all processes, since it is effectively
    ///   broadcast from Proc 0.
    ///
    /// \param Q_mine [out] View of a matrix with the same number of
    ///   rows as R_mine has columns.  On output: this process'
    ///   component of the internode Q factor.  (Write into the top
    ///   block of this process' entire Q factor, fill the rest of Q
    ///   with zeros, and call intranode TSQR's apply() on it, to get
    ///   the final explicit Q factor.)
    ///
    /// \param forceNonnegativeDiagonal [in] If true, then (if
    ///   necessary) do extra work (modifying both the Q and R
    ///   factors) in order to force the R factor to have a
    ///   nonnegative diagonal.
    void
    factorExplicit (mat_view_type R_mine,
                    mat_view_type Q_mine,
                    const bool forceNonnegativeDiagonal=false)
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! ready (), std::logic_error, "Before using DistTsqr "
         "computational methods, you must first call init() with a "
         "valid MessengerBase instance.");
      reduceBroadcastImpl_->factorExplicit (R_mine, Q_mine,
                                            forceNonnegativeDiagonal);
    }

    /// \brief Get cumulative timings for \c factorExplicit().
    ///
    /// Fill in the timings vector with cumulative timings from
    /// factorExplicit().  The vector gets resized to fit all the
    /// timings.
    void
    getFactorExplicitTimings (std::vector<TimeStats>& stats) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! ready (), std::logic_error, "Before using DistTsqr "
         "computational methods, you must first call init() with a "
         "valid MessengerBase instance.");
      reduceBroadcastImpl_->getStats (stats);
    }

    /// \brief Get labels for timings for \c factorExplicit().
    ///
    /// Fill in the labels vector with the string labels for the
    /// timings from factorExplicit().  The vector gets resized to fit
    /// all the labels.
    void
    getFactorExplicitTimingLabels (std::vector<std::string>& labels) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! ready (), std::logic_error, "Before using DistTsqr "
         "computational methods, you must first call init() with a "
         "valid MessengerBase instance.");
      reduceBroadcastImpl_->getStatsLabels (labels);
    }

    /// \brief Compute QR factorization of R factors, one per MPI process.
    ///
    /// Compute the QR factorization of the P*ncols by ncols matrix
    /// consisting of all P nodes' R_mine upper triangular matrices
    /// stacked on top of each other.  Generally these upper triangular
    /// matrices should come from the QR factorization (perhaps computed
    /// by sequential or node-parallel TSQR) of a general matrix on each
    /// node.
    ///
    /// \note "ncols" below is the number of columns in the matrix to
    ///   factor.  Should be the same on all nodes.
    ///
    /// \param R_mine [in,out] On input, an ncols by ncols upper triangular
    ///   matrix with leading dimension ncols, stored unpacked (as a general
    ///   matrix).  Elements below the diagonal are ignored.  On output, the
    ///   final R factor of the QR factorization of all nodes' different
    ///   R_mine inputs.  The final R factor is replicated over all nodes.
    ///
    /// \return Two arrays with the same number of elements: first, an
    ///   array of "local Q factors," and second, an array of "local tau
    ///   arrays."  These together form an implicit representation of
    ///   the Q factor.  They should be passed into the apply() and
    ///   explicit_Q() functions as the "factorOutput" parameter.
    FactorOutput
    factor (mat_view_type R_mine)
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! ready (), std::logic_error, "Before using DistTsqr "
         "computational methods, you must first call init() with a "
         "valid MessengerBase instance.");
      VecVec Q_factors, tau_arrays;
      DistTsqrHelper<ordinal_type, scalar_type> helper;
      const ordinal_type ncols = R_mine.extent(1);

      std::vector<scalar_type> R_local (ncols * ncols);
      MatView<ordinal_type, scalar_type> R_local_view
        (ncols, ncols, R_local.data(), ncols);
      deep_copy (R_local_view, R_mine);

      const int P = messenger_->size();
      const int my_rank = messenger_->rank();
      const int first_tag = 0;

      const ordinal_type lwork = helper.work_size (ncols);
      std::vector<scalar_type> work (lwork);
      helper.factor_helper (ncols, R_local, my_rank, 0, P-1,
                            first_tag, messenger_.get (),
                            Q_factors, tau_arrays,
                            work.data (), lwork);
      deep_copy (R_mine, R_local_view);
      return std::make_pair (Q_factors, tau_arrays);
    }

    //! Apply the result of \c factor() to the distributed matrix C.
    void
    apply (const ApplyType& apply_type,
           const ordinal_type ncols_C,
           const ordinal_type ncols_Q,
           scalar_type C_mine[],
           const ordinal_type ldc_mine,
           const FactorOutput& factor_output)
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! ready (), std::logic_error, "Before using DistTsqr "
         "computational methods, you must first call init() with a "
         "valid MessengerBase instance.");
      const bool transposed = apply_type.transposed();
      TEUCHOS_TEST_FOR_EXCEPTION(transposed, std::logic_error,
                                 "DistTsqr: Applying Q^T or Q^H has not yet "
                                 "been implemented.");
      const int P = messenger_->size();
      const int my_rank = messenger_->rank();
      const int first_tag = 0;
      std::vector<scalar_type> C_other (ncols_C * ncols_C);
      DistTsqrHelper<ordinal_type, scalar_type> helper;
      const ordinal_type lwork = helper.work_size (ncols_C);
      std::vector<scalar_type> work (lwork);

      const VecVec& Q_factors = factor_output.first;
      const VecVec& tau_arrays = factor_output.second;

      // assert (Q_factors.size() == tau_arrays.size());
      const int cur_pos = Q_factors.size() - 1;

      helper.apply_helper (apply_type, ncols_C, ncols_Q, C_mine,
                           ldc_mine, C_other.data (), my_rank, 0, P-1,
                           first_tag, messenger_.get (), Q_factors,
                           tau_arrays, cur_pos, work.data (), lwork);
    }

    //! Apply the result of \c factor() to compute the explicit Q factor.
    void
    explicit_Q (const ordinal_type ncols_Q,
                scalar_type Q_mine[],
                const ordinal_type ldq_mine,
                const FactorOutput& factor_output)
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! ready (), std::logic_error, "Before using DistTsqr "
         "computational methods, you must first call init() with a "
         "valid MessengerBase instance.");
      MatView<ordinal_type, scalar_type> Q_mine_view
        (ncols_Q, ncols_Q, Q_mine, ldq_mine);
      deep_copy (Q_mine_view, scalar_type {});

      const int myRank = messenger_->rank ();
      if (myRank == 0) {
        for (ordinal_type j = 0; j < ncols_Q; ++j) {
          Q_mine_view(j,j) = scalar_type (1.0);
        }
      }
      apply (ApplyType::NoTranspose, ncols_Q, ncols_Q,
             Q_mine, ldq_mine, factor_output);
    }

  private:
    Teuchos::RCP<MessengerBase<scalar_type>> messenger_;
    Teuchos::RCP<DistTsqrRB<ordinal_type, scalar_type>> reduceBroadcastImpl_;

    /// \brief Whether this object is ready to perform computations.
    ///
    /// It is <i>not</i> ready until after \c init() has been called.
    bool ready() const {
      return ! messenger_.is_null () &&
        ! reduceBroadcastImpl_.is_null ();
    }
  };

} // namespace TSQR

#endif // TSQR_DISTTSQR_HPP
