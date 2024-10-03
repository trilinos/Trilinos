// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_DISTTSQRHELPER_HPP
#define TSQR_DISTTSQRHELPER_HPP

#include "Tsqr_MatView.hpp"
#include "Tsqr_MessengerBase.hpp"
#include "Tsqr_Impl_CombineUser.hpp"
#include "Tsqr_Util.hpp"

#include <algorithm> // std::min, std::max
#include <sstream>
#include <stdexcept>
#include <vector>

namespace TSQR {

  /// \class DistTsqrHelper
  /// \brief Implementation of the internode part of TSQR.
  ///
  /// Implementation of the internode part of TSQR (used by DistTsqr).
  /// The only reason to mess with this class is if you want to change
  /// how the internode part of TSQR is implemented.
  template<class LocalOrdinal, class Scalar>
  class DistTsqrHelper :
    private Impl::CombineUser<LocalOrdinal, Scalar> {
  public:
    using ordinal_type = LocalOrdinal;
    using scalar_type = Scalar;

    ordinal_type work_size (const ordinal_type ncols) {
      auto& combine = this->getCombine (ncols);
      return combine.work_size (2*ncols, ncols, ncols);
    }

    void
    factor_pair (const ordinal_type ncols,
                 std::vector<scalar_type>& R_mine,
                 const ordinal_type P_mine,
                 const ordinal_type P_other,
                 const ordinal_type tag,
                 MessengerBase<scalar_type>* const messenger,
                 std::vector<std::vector<scalar_type>>& Q_factors,
                 std::vector<std::vector<scalar_type>>& tau_arrays,
                 scalar_type work[],
                 const ordinal_type lwork)
    {
      using std::endl;
      using std::ostringstream;
      using std::vector;
      using LO = ordinal_type;
      if (P_mine == P_other) {
        return; // nothing to do
      }
      const int P_top = std::min (P_mine, P_other);
      const int P_bot = std::max (P_mine, P_other);
      const LO nelts = ncols * ncols;
      const LO ldr = ncols;
      MatView<LO, scalar_type> R_mine_view
        (ncols, ncols, R_mine.data (), ldr);
      vector<scalar_type> R_other (nelts);
      MatView<LO, scalar_type> R_other_view
        (ncols, ncols, R_other.data (), ldr);
      vector<scalar_type> tau (ncols);

      // Send and receive R factor.
      messenger->swapData (R_mine.data (), R_other.data (),
                           nelts, P_other, tag);

      auto& combine = this->getCombine (ncols);
      if (P_mine == P_top) {
        combine.factor_pair (R_mine_view, R_other_view,
                             tau.data(), work, lwork);
        Q_factors.push_back (R_other);
        tau_arrays.push_back (tau);
      }
      else if (P_mine == P_bot) {
        combine.factor_pair (R_other_view, R_mine_view,
                             tau.data (), work, lwork);
        Q_factors.push_back (R_mine);
        // Make sure that the "bottom" processor gets the current R
        // factor, which is returned in R_mine.
        deep_copy (R_mine_view, R_other_view);
        tau_arrays.push_back (tau);
      }
      else {
        ostringstream os;
        os << "Should never get here: P_mine (= " << P_mine
           << ") not one of P_top, P_bot = " << P_top << ", " << P_bot;
        throw std::logic_error (os.str ());
      }
    }

    void
    factor_helper (const ordinal_type ncols,
                   std::vector<scalar_type>& R_mine,
                   const ordinal_type my_rank,
                   const ordinal_type P_first,
                   const ordinal_type P_last,
                   const ordinal_type tag,
                   MessengerBase<scalar_type>* const messenger,
                   std::vector<std::vector<scalar_type>>& Q_factors,
                   std::vector<std::vector<scalar_type>>& tau_arrays,
                   scalar_type work[],
                   const ordinal_type lwork)
    {
      using std::endl;
      using std::ostringstream;
      using std::vector;

      if (P_last <= P_first) {
        return;
      }
      else {
        const int P = P_last - P_first + 1;
        // Whether the interval [P_first, P_last] has an even number
        // of elements.  Our interval splitting scheme ensures that
        // the interval [P_first, P_mid - 1] always has an even number
        // of elements.
        const bool b_even = (P % 2 == 0);
        // We split the interval [P_first, P_last] into 2 intervals:
        // [P_first, P_mid-1], and [P_mid, P_last].  We bias the
        // splitting procedure so that the lower interval always has
        // an even number of processor ranks, and never has fewer
        // processor ranks than the higher interval.
        const int P_mid = b_even ? (P_first + P/2) : (P_first + P/2 + 1);

        if (my_rank < P_mid) { // Interval [P_first, P_mid-1]
          factor_helper (ncols, R_mine, my_rank, P_first, P_mid - 1,
                         tag + 1, messenger, Q_factors, tau_arrays,
                         work, lwork);

          // If there aren't an even number of processors in the
          // original interval, then the last processor in the lower
          // interval has to skip this round.
          if (b_even || my_rank < P_mid - 1) {
            const int my_offset = my_rank - P_first;
            const int P_other = P_mid + my_offset;
            if (P_other < P_mid || P_other > P_last) {
              throw std::logic_error ("P_other not in [P_mid,P_last] range");
            }
            factor_pair (ncols, R_mine, my_rank, P_other, tag,
                         messenger, Q_factors, tau_arrays,
                         work, lwork);
          }
          // If I'm skipping this round, get the "current" R factor
          // from P_mid.
          if (! b_even && my_rank == P_mid - 1) {
            const int theTag = 142; // magic constant
            messenger->recv (R_mine.data (), ncols*ncols, P_mid,
                             theTag);
          }
        }
        else { // Interval [P_mid, P_last]
          factor_helper (ncols, R_mine, my_rank, P_mid, P_last,
                         tag + 1, messenger, Q_factors, tau_arrays,
                         work, lwork);
          const int my_offset = my_rank - P_mid;
          const int P_other = P_first + my_offset;

          if (P_other < P_first || P_other >= P_mid) {
            throw std::logic_error ("P_other not in [P_first,"
                                    "P_mid-1] range");
          }
          factor_pair (ncols, R_mine, my_rank, P_other, tag,
                       messenger, Q_factors, tau_arrays, work, lwork);

          // If Proc P_mid-1 is skipping this round, Proc P_mid will
          // send it the "current" R factor.
          if (! b_even) {
            const int theTag = 142; // magic constant
            messenger->send (R_mine.data(), ncols*ncols, P_mid-1, theTag);
          }
        }
      }
    }

    void
    apply_pair (const ApplyType& apply_type,
                const ordinal_type ncols_C,
                const ordinal_type ncols_Q,
                scalar_type C_mine[],
                const ordinal_type ldc_mine,
                scalar_type C_other[], // contiguous ncols_C x ncols_C scratch
                const ordinal_type P_mine,
                const ordinal_type P_other,
                const ordinal_type tag,
                MessengerBase<scalar_type>* const messenger,
                const std::vector<scalar_type>& Q_cur,
                const std::vector<scalar_type>& tau_cur,
                scalar_type work[],
                const ordinal_type lwork)
    {
      using std::endl;
      using std::ostringstream;
      using std::vector;
      using LO = ordinal_type;
      using const_mat_view_type = MatView<LO, const scalar_type>;
      using mat_view_type = MatView<LO, scalar_type>;

      if (P_mine == P_other) {
        return; // nothing to do
      }
      const int P_top = std::min (P_mine, P_other);
      const int P_bot = std::max (P_mine, P_other);
      const LO nelts = ncols_C * ncols_C;
      const LO ldq = ncols_Q;
      const LO ldc_other = ncols_C;

      // Send and receive C_mine resp. C_other to the other processor of
      // the pair.
      messenger->swapData (C_mine, C_other, nelts, P_other, tag);

      const_mat_view_type Q_bot
        (ncols_Q, ncols_Q, Q_cur.data (), ldq);
      auto& combine = this->getCombine (std::max (ncols_Q, ncols_C));
      if (P_mine == P_top) {
        mat_view_type C_top (ncols_Q, ncols_C, C_mine, ldc_mine);
        mat_view_type C_bot (ncols_Q, ncols_C, C_other, ldc_other);
        combine.apply_pair (apply_type, Q_bot, tau_cur.data (),
                            C_top, C_bot, work, lwork);
      }
      else if (P_mine == P_bot) {
        mat_view_type C_top (ncols_Q, ncols_C, C_other, ldc_other);
        mat_view_type C_bot (ncols_Q, ncols_C, C_mine, ldc_mine);
        combine.apply_pair (apply_type, Q_bot, tau_cur.data (),
                            C_top, C_bot, work, lwork);
      }
      else {
        ostringstream os;
        os << "Should never get here: P_mine (= " << P_mine
           << ") not one of P_top, P_bot = " << P_top << ", "
           << P_bot;
        throw std::logic_error (os.str ());
      }
    }

    void
    apply_helper (const ApplyType& apply_type,
                  const ordinal_type ncols_C,
                  const ordinal_type ncols_Q,
                  scalar_type C_mine[],
                  const ordinal_type ldc_mine,
                  scalar_type C_other[], // contiguous ncols_C x ncols_C scratch
                  const ordinal_type my_rank,
                  const ordinal_type P_first,
                  const ordinal_type P_last,
                  const ordinal_type tag,
                  MessengerBase<scalar_type>* const messenger,
                  const std::vector<std::vector<scalar_type>>& Q_factors,
                  const std::vector<std::vector<scalar_type>>& tau_arrays,
                  const ordinal_type cur_pos,
                  scalar_type work[],
                  const ordinal_type lwork)
    {
      using std::endl;
      using std::ostringstream;
      using std::vector;

      if (P_last <= P_first) {
        return;
      }
      else {
        const int P = P_last - P_first + 1;
        // Whether the interval [P_first, P_last] has an even number
        // of elements.  Our interval splitting scheme ensures that
        // the interval [P_first, P_mid - 1] always has an even number
        // of elements.
        const bool b_even = (P % 2 == 0);
        // We split the interval [P_first, P_last] into 2 intervals:
        // [P_first, P_mid-1], and [P_mid, P_last].  We bias the
        // splitting procedure so that the lower interval always has
        // an even number of processor ranks, and never has fewer
        // processor ranks than the higher interval.
        const int P_mid = b_even ? (P_first + P/2) : (P_first + P/2 + 1);

        if (my_rank < P_mid) { // Interval [P_first, P_mid - 1]
          const bool b_participating = b_even || my_rank < P_mid - 1;

          if (cur_pos < 0) {
            ostringstream os;
            os << "On Proc " << my_rank << ": cur_pos (= " << cur_pos
               << ") < 0; lower interval [" << P_first << "," << (P_mid-1)
               << "]; original interval [" << P_first << "," << P_last
               << "]" << endl;
            throw std::logic_error (os.str());
          }

          // If there aren't an even number of processors in the
          // original interval, then the last processor in the lower
          // interval has to skip this round.  Since we skip this
          // round, don't decrement cur_pos (else we'll skip an entry
          // and eventually fall off the front of the array.
          int new_cur_pos;
          if (b_even || my_rank < P_mid - 1) {
            if (! b_participating) {
              throw std::logic_error("Should never get here");
            }

            const int my_offset = my_rank - P_first;
            const int P_other = P_mid + my_offset;
            // assert (P_mid <= P_other && P_other <= P_last);
            if (P_other < P_mid || P_other > P_last) {
              throw std::logic_error("Should never get here");
            }
            apply_pair (apply_type, ncols_C, ncols_Q, C_mine,
                        ldc_mine, C_other, my_rank, P_other,
                        tag, messenger, Q_factors[cur_pos],
                        tau_arrays[cur_pos], work, lwork);
            new_cur_pos = cur_pos - 1;
          }
          else {
            if (b_participating) {
              throw std::logic_error("Should never get here");
            }
            new_cur_pos = cur_pos;
          }
          apply_helper (apply_type, ncols_C, ncols_Q, C_mine,
                        ldc_mine, C_other, my_rank, P_first,
                        P_mid - 1, tag + 1, messenger, Q_factors,
                        tau_arrays, new_cur_pos, work, lwork);
        }
        else {
          if (cur_pos < 0) {
            ostringstream os;
            os << "On Proc " << my_rank << ": cur_pos (= " << cur_pos
               << ") < 0; upper interval [" << P_mid << "," << P_last
               << "]; original interval [" << P_first << "," << P_last
               << "]" << endl;
            throw std::logic_error (os.str ());
          }

          const int my_offset = my_rank - P_mid;
          const int P_other = P_first + my_offset;
          // assert (0 <= P_other && P_other < P_mid);
          apply_pair (apply_type, ncols_C, ncols_Q, C_mine, ldc_mine,
                      C_other, my_rank, P_other, tag, messenger,
                      Q_factors[cur_pos], tau_arrays[cur_pos],
                      work, lwork);
          apply_helper (apply_type, ncols_C, ncols_Q, C_mine, ldc_mine,
                        C_other, my_rank, P_mid, P_last, tag + 1,
                        messenger, Q_factors, tau_arrays, cur_pos - 1,
                        work, lwork);
        }
      }
    }
  };

} // namespace TSQR

#endif // TSQR_DISTTSQRHELPER_HPP
