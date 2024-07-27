// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_DISTTSQRRB_HPP
#define TSQR_DISTTSQRRB_HPP

#include "Tsqr_ApplyType.hpp"
#include "Tsqr_Impl_CombineUser.hpp"
#include "Tsqr_Matrix.hpp"
#include "Tsqr_StatTimeMonitor.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace TSQR {

  /// \namespace details
  /// \brief TSQR implementation details.
  /// \author Mark Hoemmen
  ///
  /// \warning TSQR users should not use anything in this namespace.
  ///   They should not even assume that the namespace will continue
  ///   to exist between releases.  The namespace's name itself or
  ///   anything it contains may change at any time.
  namespace details {

    // Force the diagonal of R_mine to be nonnegative, where
    // Q_mine*R_mine is a QR factorization.
    //
    // We only made this a class because C++ (pre-C++11) does not
    // allow partial specialization of template functions.
    template<class LocalOrdinal, class Scalar, bool isComplex>
    class NonnegDiagForcer {
    public:
      typedef MatView<LocalOrdinal, Scalar> mat_view_type;

      // Force the diagonal of R_mine to be nonnegative, where
      // Q_mine*R_mine is a QR factorization.
      void force (mat_view_type Q_mine, mat_view_type R_mine);
    };

    // The complex-arithmetic specialization does nothing, since
    // _GEQR{2,F} for complex arithmetic returns an R factor with
    // nonnegative diagonal already.
    template<class LocalOrdinal, class Scalar>
    class NonnegDiagForcer<LocalOrdinal, Scalar, true> {
    public:
      typedef MatView<LocalOrdinal, Scalar> mat_view_type;

      void force (mat_view_type Q_mine, mat_view_type R_mine) {
        (void) Q_mine;
        (void) R_mine;
      }
    };

    // Real-arithmetic specialization.
    template<class LocalOrdinal, class Scalar>
    class NonnegDiagForcer<LocalOrdinal, Scalar, false> {
    public:
      typedef MatView<LocalOrdinal, Scalar> mat_view_type;

      void force (mat_view_type Q_mine, mat_view_type R_mine) {
        typedef Teuchos::ScalarTraits<Scalar> STS;

        if (Q_mine.extent(0) > 0 && Q_mine.extent(1) > 0) {
          for (int k = 0; k < R_mine.extent(1); ++k) {
            if (R_mine(k,k) < STS::zero()) {
              // Scale column k of Q_mine.  We use a raw pointer since
              // typically there are many rows in Q_mine, so this
              // operation should be fast.
              Scalar* const Q_k = &Q_mine(0,k);
              for (int i = 0; i < Q_mine.extent(0); ++i) {
                Q_k[i] = -Q_k[i];
              }
              // Scale row k of R_mine.  R_mine is upper triangular,
              // so we only have to scale right of (and including) the
              // diagonal entry.
              for (int j = k; j < R_mine.extent(1); ++j) {
                R_mine(k,j) = -R_mine(k,j);
              }
            }
          }
        }
      }
    };
  } // namespace details

  /// \class DistTsqrRB
  /// \brief Reduce-and-Broadcast (RB) version of DistTsqr.
  /// \author Mark Hoemmen
  ///
  /// \tparam LocalOrdinal Corresponds to the "local ordinal" template
  ///   parameter of Tpetra objects (though TSQR is not Tpetra-specific).
  ///
  /// \tparam Scalar Corresponds to the "scalar" template parameter of
  ///   Tpetra objects (though TSQR is not Tpetra-specific).
  ///
  /// This class implements the Reduce-and-Broadcast (RB) version of
  /// DistTsqr.  DistTsqr factors a vertical stack of n by n R
  /// factors, one per MPI process.  Only the final R factor is
  /// broadcast.  The implicit Q factor data stay on the MPI process
  /// where they were computed.
  template<class LocalOrdinal, class Scalar>
  class DistTsqrRB : private Impl::CombineUser<LocalOrdinal, Scalar> {
  public:
    using ordinal_type = LocalOrdinal;
    using scalar_type = Scalar;
    using magnitude_type =
      typename Teuchos::ScalarTraits<scalar_type>::magnitudeType;
    using mat_view_type = MatView<ordinal_type, scalar_type>;
    using matrix_type = Matrix<ordinal_type, scalar_type>;
    using rank_type = int;

    /// \brief Constructor
    ///
    /// \param messenger [in/out] Smart pointer to a wrapper handling
    ///   communication between MPI process(es).
    DistTsqrRB (const Teuchos::RCP< MessengerBase< scalar_type > >& messenger) :
      messenger_ (messenger),
      totalTime_ (Teuchos::TimeMonitor::getNewTimer ("DistTsqrRB::factorExplicit() total time")),
      reduceCommTime_ (Teuchos::TimeMonitor::getNewTimer ("DistTsqrRB::factorReduce() communication time")),
      reduceTime_ (Teuchos::TimeMonitor::getNewTimer ("DistTsqrRB::factorReduce() total time")),
      bcastCommTime_ (Teuchos::TimeMonitor::getNewTimer ("DistTsqrRB::explicitQBroadcast() communication time")),
      bcastTime_ (Teuchos::TimeMonitor::getNewTimer ("DistTsqrRB::explicitQBroadcast() total time"))
    {}

    /// \brief Fill stats with cumulative timings from \c factorExplicit().
    ///
    /// Fill in the timings vector with cumulative timings from
    /// factorExplicit().  The vector gets resized if necessary to fit
    /// all the timings.
    void
    getStats (std::vector< TimeStats >& stats) const
    {
      const int numTimers = 5;
      stats.resize (std::max (stats.size(), static_cast<size_t>(numTimers)));

      stats[0] = totalStats_;
      stats[1] = reduceCommStats_;
      stats[2] = reduceStats_;
      stats[3] = bcastCommStats_;
      stats[4] = bcastStats_;
    }

    /// \brief Fill labels with timer labels from \c factorExplicit().
    ///
    /// Fill in the labels vector with the string labels for the
    /// timings from factorExplicit().  The vector gets resized if
    /// necessary to fit all the labels.
    void
    getStatsLabels (std::vector<std::string>& labels) const
    {
      const int numTimers = 5;
      labels.resize (std::max (labels.size (), size_t (numTimers)));

      labels[0] = totalTime_->name();
      labels[1] = reduceCommTime_->name();
      labels[2] = reduceTime_->name();
      labels[3] = bcastCommTime_->name();
      labels[4] = bcastTime_->name();
    }

    /// Whether or not all diagonal entries of the R factor computed
    /// by the QR factorization are guaranteed to be nonnegative.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      // FIXME (20 Dec 2019) If the combine type is dynamic, we can't
      // answer this question without knowing the number of columns.
      // Just guess for now.
      constexpr LocalOrdinal fakeNumCols = 10;
      auto& c = this->getCombine (fakeNumCols);
      return c.QR_produces_R_factor_with_nonnegative_diagonal ();
    }

    /// \brief Internode TSQR with explicit Q factor
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
      StatTimeMonitor totalMonitor (*totalTime_, totalStats_);

      // Dimension sanity checks.  R_mine should have at least as many
      // rows as columns (since we will be working on the upper
      // triangle).  Q_mine should have the same number of rows as
      // R_mine has columns, but Q_mine may have any number of
      // columns.  (It depends on how many columns of the explicit Q
      // factor we want to compute.)
      if (R_mine.extent(0) < R_mine.extent(1)) {
        std::ostringstream os;
        os << "R factor input has fewer rows (" << R_mine.extent(0)
           << ") than columns (" << R_mine.extent(1) << ")";
        // This is a logic error because TSQR users should not be
        // calling this method directly.
        throw std::logic_error (os.str());
      }
      else if (Q_mine.extent(0) != R_mine.extent(1)) {
        std::ostringstream os;
        os << "Q factor input must have the same number of rows as the R "
          "factor input has columns.  Q has " << Q_mine.extent(0)
           << " rows, but R has " << R_mine.extent(1) << " columns.";
        // This is a logic error because TSQR users should not be
        // calling this method directly.
        throw std::logic_error (os.str());
      }

      // The factorization is a recursion over processors [P_first, P_last].
      const rank_type P_mine = messenger_->rank();
      const rank_type P_first = 0;
      const rank_type P_last = messenger_->size() - 1;

      // Intermediate Q factors are stored implicitly.  QFactors[k] is
      // an upper triangular matrix of Householder reflectors, and
      // tauArrays[k] contains its corresponding scaling factors (TAU,
      // in LAPACK notation).  These two arrays will be filled in by
      // factorReduce().  Different MPI processes will have different
      // numbers of elements in these arrays.  In fact, on some
      // processes these arrays may be empty on output.  This is a
      // feature, not a bug!
      //
      // Even though QFactors and tauArrays have the same type has the
      // first resp. second elements of DistTsqr::FactorOutput, they
      // are not compatible with the output of DistTsqr::factor() and
      // cannot be used as the input to DistTsqr::apply() or
      // DistTsqr::explicit_Q().  This is because factor() computes a
      // general factorization suitable for applying Q (or Q^T or Q^*)
      // to any compatible matrix, whereas factorExplicit() computes a
      // factorization specifically for the purpose of forming the
      // explicit Q factor.  The latter lets us use a broadcast to
      // compute Q, rather than a more message-intensive all-to-all
      // (butterfly).
      std::vector< matrix_type > QFactors;
      std::vector< std::vector< scalar_type > > tauArrays;

      {
        StatTimeMonitor reduceMonitor (*reduceTime_, reduceStats_);
        factorReduce (R_mine, P_mine, P_first, P_last, QFactors, tauArrays);
      }

      if (QFactors.size() != tauArrays.size()) {
        std::ostringstream os;
        os << "QFactors and tauArrays should have the same number of element"
          "s after factorReduce() returns, but they do not.  QFactors has "
           << QFactors.size() << " elements, but tauArrays has "
           << tauArrays.size() << " elements.";
        throw std::logic_error (os.str());
      }

      deep_copy (Q_mine, scalar_type {});
      if (messenger_->rank() == 0) {
        for (ordinal_type j = 0; j < Q_mine.extent(1); ++j) {
          Q_mine(j, j) = scalar_type (1);
        }
      }
      // Scratch space for computing results to send to other processors.
      matrix_type Q_other (Q_mine.extent(0), Q_mine.extent(1), scalar_type {});
      const rank_type numSteps = QFactors.size() - 1;

      {
        StatTimeMonitor bcastMonitor (*bcastTime_, bcastStats_);
        explicitQBroadcast (R_mine, Q_mine, Q_other.view(),
                            P_mine, P_first, P_last,
                            numSteps, QFactors, tauArrays);
      }

      if (forceNonnegativeDiagonal &&
          ! QR_produces_R_factor_with_nonnegative_diagonal()) {
        using STS = Teuchos::ScalarTraits<Scalar>;
        details::NonnegDiagForcer<LocalOrdinal, Scalar, STS::isComplex> forcer;
        forcer.force (Q_mine, R_mine);
      }
    }

  private:

    void
    factorReduce (mat_view_type R_mine,
                  const rank_type P_mine,
                  const rank_type P_first,
                  const rank_type P_last,
                  std::vector<matrix_type>& QFactors,
                  std::vector<std::vector<scalar_type>>& tauArrays)
    {
      if (P_last < P_first) {
        std::ostringstream os;
        os << "factorReduce: Interval [P_first=" << P_first
           << ", P_last=" << P_last << "] is invalid.";
        throw std::logic_error (os.str());
      }
      else if (P_mine < P_first || P_mine > P_last) {
        std::ostringstream os;
        os << "factorReduce: P_mine=" << P_mine << " is not in "
           << "current process rank interval [P_first=" << P_first
           << ", P_last=" << P_last << "]";
        throw std::logic_error (os.str());
      }
      else if (P_last == P_first) {
        return; // skip singleton intervals (see explanation below)
      }
      else {
        // Recurse on two intervals: [P_first, P_mid-1] and [P_mid,
        // P_last].  For example, if [P_first, P_last] = [0, 9], P_mid
        // = floor( (0+9+1)/2 ) = 5 and the intervals are [0,4] and
        // [5,9].
        //
        // If [P_first, P_last] = [4,6], P_mid = floor( (4+6+1)/2 ) =
        // 5 and the intervals are [4,4] (a singleton) and [5,6].  The
        // latter case shows that singleton intervals may arise.  We
        // treat them as a base case in the recursion.  Process 4
        // won't be skipped completely, though; it will get combined
        // with the result from [5,6].

        // Adding 1 and doing integer division works like "ceiling."
        const rank_type P_mid = (P_first + P_last + 1) / 2;

        if (P_mine < P_mid) { // Interval [P_first, P_mid-1]
          factorReduce (R_mine, P_mine, P_first, P_mid - 1,
                        QFactors, tauArrays);
        }
        else { // Interval [P_mid, P_last]
          factorReduce (R_mine, P_mine, P_mid, P_last,
                        QFactors, tauArrays);
        }

        // This only does anything if P_mine is either P_first or P_mid.
        if (P_mine == P_first) {
          const ordinal_type numCols = R_mine.extent(1);
          matrix_type R_other (numCols, numCols);
          recv_R (R_other, P_mid);

          std::vector<scalar_type> tau (numCols);

          auto& combine = this->getCombine (numCols);
          const ordinal_type lwork =
            combine.work_size (2 * numCols, numCols, numCols);
          work_.resize (lwork);
          combine.factor_pair (R_mine, R_other.view (),
                               tau.data (), work_.data (), lwork);
          QFactors.push_back (R_other);
          tauArrays.push_back (tau);
        }
        else if (P_mine == P_mid) {
          send_R (R_mine, P_first);
        }
      }
    }

    void
    explicitQBroadcast (mat_view_type R_mine,
                        mat_view_type Q_mine,
                        mat_view_type Q_other, // workspace
                        const rank_type P_mine,
                        const rank_type P_first,
                        const rank_type P_last,
                        const rank_type curpos,
                        std::vector<matrix_type>& QFactors,
                        std::vector<std::vector<scalar_type>>& tauArrays)
    {
      using LO = LocalOrdinal;

      if (P_last < P_first) {
        std::ostringstream os;
        os << "explicitQBroadcast: interval [P_first=" << P_first
           << ", P_last=" << P_last << "] is invalid.";
        throw std::logic_error (os.str());
      }
      else if (P_mine < P_first || P_mine > P_last) {
        std::ostringstream os;
        os << "explicitQBroadcast: P_mine=" << P_mine << " is not "
          "in current process rank interval [P_first = " << P_first
           << ", P_last = " << P_last << "]";
        throw std::logic_error (os.str());
      }
      else if (P_last == P_first) {
        return; // skip singleton intervals
      }
      else {
        // Adding 1 and integer division works like "ceiling."
        const rank_type P_mid = (P_first + P_last + 1) / 2;
        rank_type newpos = curpos;
        if (P_mine == P_first) {
          if (curpos < 0) {
            std::ostringstream os;
            os << "Programming error: On the current P_first (= "
               << P_first << ") proc: curpos (= " << curpos << ") < 0";
            throw std::logic_error (os.str());
          }
          // Q_impl, tau: implicitly stored local Q factor.
          auto Q_bot = QFactors[curpos].view ();
          const scalar_type* tau = tauArrays[curpos].data ();

          // Apply implicitly stored local Q factor to
          //   [Q_mine;
          //    Q_other]
          // where Q_other = zeros(Q_mine.extent(0), Q_mine.extent(1)).
          // Overwrite both Q_mine and Q_other with the result.
          deep_copy (Q_other, scalar_type {});

          const LO pair_nrows
            (Q_mine.extent (0) + Q_other.extent (0));
          const LO pair_ncols (Q_mine.extent (1));
          auto& combine = this->getCombine (pair_ncols);
          const LO lwork =
            combine.work_size (pair_nrows, pair_ncols, pair_ncols);
          if (lwork > LO (work_.size ())) {
            work_.resize (lwork);
          }
          combine.apply_pair (ApplyType::NoTranspose, Q_bot, tau,
                              Q_mine, Q_other, work_.data (), lwork);
          // Send the resulting Q_other, and the final R factor, to P_mid.
          send_Q_R (Q_other, R_mine, P_mid);
          newpos = curpos - 1;
        }
        else if (P_mine == P_mid) {
          // P_first computed my explicit Q factor component.
          // Receive it, and the final R factor, from P_first.
          recv_Q_R (Q_mine, R_mine, P_first);
        }

        if (P_mine < P_mid) { // Interval [P_first, P_mid-1]
          explicitQBroadcast (R_mine, Q_mine, Q_other,
                              P_mine, P_first, P_mid - 1,
                              newpos, QFactors, tauArrays);
        }
        else { // Interval [P_mid, P_last]
          explicitQBroadcast (R_mine, Q_mine, Q_other,
                              P_mine, P_mid, P_last,
                              newpos, QFactors, tauArrays);
        }
      }
    }

    template< class ConstMatrixType1, class ConstMatrixType2 >
    void
    send_Q_R (const ConstMatrixType1& Q,
              const ConstMatrixType2& R,
              const rank_type destProc)
    {
      StatTimeMonitor bcastCommMonitor (*bcastCommTime_, bcastCommStats_);

      const ordinal_type R_numCols = R.extent(1);
      const ordinal_type Q_size = Q.extent(0) * Q.extent(1);
      const ordinal_type R_size = (R_numCols * (R_numCols + 1)) / 2;
      const ordinal_type numElts = Q_size + R_size;

      // Don't shrink the workspace array; doing so would still be
      // correct, but may require reallocation of data when it needs
      // to grow again.
      work_.resize (numElts);

      // Pack the Q data into the workspace array.
      mat_view_type Q_contig (Q.extent (0), Q.extent (1),
                              work_.data (), Q.extent (0));
      deep_copy (Q_contig, Q);
      // Pack the R data into the workspace array.
      pack_R (R, &work_[Q_size]);
      messenger_->send (work_.data (), numElts, destProc, 0);
    }

    template< class MatrixType1, class MatrixType2 >
    void
    recv_Q_R (MatrixType1& Q,
              MatrixType2& R,
              const rank_type srcProc)
    {
      StatTimeMonitor bcastCommMonitor (*bcastCommTime_, bcastCommStats_);

      const ordinal_type R_numCols = R.extent(1);
      const ordinal_type Q_size = Q.extent(0) * Q.extent(1);
      const ordinal_type R_size = (R_numCols * (R_numCols + 1)) / 2;
      const ordinal_type numElts = Q_size + R_size;

      // Don't shrink the workspace array; doing so would still be
      // correct, but may require reallocation of data when it needs
      // to grow again.
      work_.resize (numElts);

      messenger_->recv (work_.data (), numElts, srcProc, 0);

      // Unpack the C data from the workspace array.
      deep_copy (Q, mat_view_type (Q.extent (0), Q.extent (1),
                                   work_.data (), Q.extent (0)));
      // Unpack the R data from the workspace array.
      unpack_R (R, &work_[Q_size]);
    }

    template< class ConstMatrixType >
    void
    send_R (const ConstMatrixType& R, const rank_type destProc)
    {
      StatTimeMonitor reduceCommMonitor (*reduceCommTime_, reduceCommStats_);

      const ordinal_type numCols = R.extent(1);
      const ordinal_type numElts = (numCols * (numCols+1)) / 2;

      // Don't shrink the workspace array; doing so would still be
      // correct, but may require reallocation of data when it needs
      // to grow again.
      work_.resize (numElts);
      // Pack the R data into the workspace array.
      pack_R (R, work_.data ());
      messenger_->send (work_.data (), numElts, destProc, 0);
    }

    template< class MatrixType >
    void
    recv_R (MatrixType& R, const rank_type srcProc)
    {
      StatTimeMonitor reduceCommMonitor (*reduceCommTime_, reduceCommStats_);

      const ordinal_type numCols = R.extent(1);
      const ordinal_type numElts = (numCols * (numCols+1)) / 2;

      // Don't shrink the workspace array; doing so would still be
      // correct, but may require reallocation of data when it needs
      // to grow again.
      work_.resize (numElts);
      messenger_->recv (work_.data (), numElts, srcProc, 0);
      // Unpack the R data from the workspace array.
      unpack_R (R, work_.data ());
    }

    template< class MatrixType >
    static void
    unpack_R (MatrixType& R, const scalar_type buf[])
    {
      // FIXME (mfh 08 Dec 2019) Rewrite to use deep_copy; we don't
      // want to access Matrix or MatView entries on host directly any
      // more.
      ordinal_type curpos = 0;
      for (ordinal_type j = 0; j < R.extent(1); ++j) {
        scalar_type* const R_j = &R(0, j);
        for (ordinal_type i = 0; i <= j; ++i) {
          R_j[i] = buf[curpos++];
        }
      }
    }

    template< class ConstMatrixType >
    static void
    pack_R (const ConstMatrixType& R, scalar_type buf[])
    {
      ordinal_type curpos = 0;
      for (ordinal_type j = 0; j < R.extent(1); ++j) {
        const scalar_type* const R_j = &R(0, j);
        for (ordinal_type i = 0; i <= j; ++i) {
          buf[curpos++] = R_j[i];
        }
      }
    }

  private:
    Teuchos::RCP<MessengerBase<scalar_type>> messenger_;
    std::vector<scalar_type> work_;

    // Timers for various phases of the factorization.  Time is
    // cumulative over all calls of factorExplicit().
    Teuchos::RCP<Teuchos::Time> totalTime_;
    Teuchos::RCP<Teuchos::Time> reduceCommTime_;
    Teuchos::RCP<Teuchos::Time> reduceTime_;
    Teuchos::RCP<Teuchos::Time> bcastCommTime_;
    Teuchos::RCP<Teuchos::Time> bcastTime_;

    TimeStats totalStats_;
    TimeStats reduceCommStats_;
    TimeStats reduceStats_;
    TimeStats bcastCommStats_;
    TimeStats bcastStats_;
  };

} // namespace TSQR

#endif // TSQR_DISTTSQRRB_HPP
