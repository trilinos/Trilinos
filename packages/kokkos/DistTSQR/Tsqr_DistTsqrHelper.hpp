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

#ifndef __TSQR_Tsqr_DistTsqrHelper_hpp
#define __TSQR_Tsqr_DistTsqrHelper_hpp

#include <Tsqr_MatView.hpp>
#include <Tsqr_MessengerBase.hpp>
#include <Tsqr_Combine.hpp>
#include <Tsqr_Util.hpp>

#include <algorithm> // std::min, std::max
#include <sstream>
#include <stdexcept>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class DistTsqrHelper
  /// \brief Implementation of the internode part of TSQR
  ///
  /// Implementation of the internode part of TSQR (used by DistTsqr).
  /// The only reason to mess with this class is if you want to change
  /// how the internode part of TSQR is implemented.
  template< class LocalOrdinal, class Scalar >
  class DistTsqrHelper {
  public:
    DistTsqrHelper () {}

    void
    factor_pair (const LocalOrdinal ncols,
		 std::vector< Scalar >& R_mine,
		 const LocalOrdinal P_mine,
		 const LocalOrdinal P_other,
		 const LocalOrdinal tag, 
		 MessengerBase< Scalar >* const messenger,
		 std::vector< std::vector< Scalar > >& Q_factors,
		 std::vector< std::vector< Scalar > >& tau_arrays,
		 std::vector< Scalar >& work) 
    {
      using std::endl;
      using std::ostringstream;
      using std::vector;

      if (P_mine == P_other)
	return; // nothing to do

      const int P_top = std::min (P_mine, P_other);
      const int P_bot = std::max (P_mine, P_other);
      const LocalOrdinal nelts = ncols * ncols;
      const LocalOrdinal ldr = ncols;
      vector< Scalar > R_other (nelts);
      vector< Scalar > tau (ncols);

      // Send and receive R factor.
      messenger->swapData (&R_mine[0], &R_other[0], nelts, P_other, tag);

      Combine< LocalOrdinal, Scalar > combine;
      if (P_mine == P_top)
	{
	  combine.factor_pair (ncols, &R_mine[0], ldr, &R_other[0], ldr, &tau[0], &work[0]);
	  Q_factors.push_back (R_other);
	  tau_arrays.push_back (tau);
	}
      else if (P_mine == P_bot)
	{
	  combine.factor_pair (ncols, &R_other[0], ldr, &R_mine[0], ldr, &tau[0], &work[0]);
	  Q_factors.push_back (R_mine);
	  // Make sure that the "bottom" processor gets the current R
	  // factor, which is returned in R_mine.
	  copy_matrix (ncols, ncols, &R_mine[0], ldr, &R_other[0], ldr);
	  tau_arrays.push_back (tau);
	}
      else
	{
	  // mfh 16 Apr 2010: the troubles with assert statements are as follows:
	  //
	  // 1. They go away in a release build.
	  // 2. They don't often print out useful diagnostic information.
	  // 3. If you mistype the assert, like "assert(errcode = 1);" instead of 
	  //    "assert(errcode == 1)", you'll get false positives.
	  ostringstream os;
	  os << "Should never get here: P_mine (= " << P_mine 
	     << ") not one of P_top, P_bot = " << P_top << ", " << P_bot;
	  throw std::logic_error (os.str());
	}
    }

    void
    factor_helper (const LocalOrdinal ncols,
		   std::vector< Scalar >& R_mine,
		   const LocalOrdinal my_rank,
		   const LocalOrdinal P_first,
		   const LocalOrdinal P_last,
		   const LocalOrdinal tag,
		   MessengerBase< Scalar >* const messenger,
		   std::vector< std::vector< Scalar > >& Q_factors,
		   std::vector< std::vector< Scalar > >& tau_arrays,
		   std::vector< Scalar >& work)
    {
      using std::endl;
      using std::ostringstream;
      using std::vector;

      if (P_last <= P_first)
	return;
      else
	{
	  const int P = P_last - P_first + 1;
	  // Whether the interval [P_first, P_last] has an even number of
	  // elements.  Our interval splitting scheme ensures that the
	  // interval [P_first, P_mid - 1] always has an even number of
	  // elements.
	  const bool b_even = (P % 2 == 0);
	  // We split the interval [P_first, P_last] into 2 intervals:
	  // [P_first, P_mid-1], and [P_mid, P_last].  We bias the
	  // splitting procedure so that the lower interval always has an
	  // even number of processor ranks, and never has fewer processor
	  // ranks than the higher interval.
	  const int P_mid = b_even ? (P_first + P/2) : (P_first + P/2 + 1); 

	  if (my_rank < P_mid) // Interval [P_first, P_mid-1]
	    {
	      factor_helper (ncols, R_mine, my_rank, P_first, P_mid - 1, 
			     tag + 1, messenger, Q_factors, tau_arrays, work);

	      // If there aren't an even number of processors in the
	      // original interval, then the last processor in the lower
	      // interval has to skip this round.
	      if (b_even || my_rank < P_mid - 1)
		{
		  const int my_offset = my_rank - P_first;
		  const int P_other = P_mid + my_offset;
		  if (P_other < P_mid || P_other > P_last)
		    throw std::logic_error ("P_other not in [P_mid,P_last] range");

		  factor_pair (ncols, R_mine, my_rank, P_other, tag, 
			       messenger, Q_factors, tau_arrays, work);
		}

	      // If I'm skipping this round, get the "current" R factor
	      // from P_mid.
	      if (! b_even && my_rank == P_mid - 1)
		{
		  const int theTag = 142; // magic constant
		  messenger->recv (&R_mine[0], ncols*ncols, P_mid, theTag);
		}
	    }
	  else // Interval [P_mid, P_last]
	    {
	      factor_helper (ncols, R_mine, my_rank, P_mid, P_last, 
			     tag + 1, messenger, Q_factors, tau_arrays, work);

	      const int my_offset = my_rank - P_mid;
	      const int P_other = P_first + my_offset;

	      if (P_other < P_first || P_other >= P_mid)
		throw std::logic_error ("P_other not in [P_first,P_mid-1] range");
	      factor_pair (ncols, R_mine, my_rank, P_other, tag, 
			   messenger, Q_factors, tau_arrays, work);

	      // If Proc P_mid-1 is skipping this round, Proc P_mid will
	      // send it the "current" R factor.
	      if (! b_even)
		{
		  const int theTag = 142; // magic constant
		  messenger->send (&R_mine[0], ncols*ncols, P_mid-1, theTag);
		}
	    }
	}
    }

    void
    apply_pair (const ApplyType& apply_type,
		const LocalOrdinal ncols_C,
		const LocalOrdinal ncols_Q,
		Scalar C_mine[],
		const LocalOrdinal ldc_mine,
		Scalar C_other[], // contiguous ncols_C x ncols_C scratch
		const LocalOrdinal P_mine,
		const LocalOrdinal P_other,
		const LocalOrdinal tag, 
		MessengerBase< Scalar >* const messenger,
		const std::vector< Scalar >& Q_cur,
		const std::vector< Scalar >& tau_cur,
		std::vector< Scalar >& work)
    {
      using std::endl;
      using std::ostringstream;
      using std::vector;

      if (P_mine == P_other)
	return; // nothing to do
    
      const int P_top = std::min (P_mine, P_other);
      const int P_bot = std::max (P_mine, P_other);
    
      const LocalOrdinal nelts = ncols_C * ncols_C;
      const LocalOrdinal ldq = ncols_Q;
      const LocalOrdinal ldc_other = ncols_C;

      // Send and receive C_mine resp. C_other to the other processor of
      // the pair.
      messenger->swapData (&C_mine[0], &C_other[0], nelts, P_other, tag);

      Combine< LocalOrdinal, Scalar > combine;
      if (P_mine == P_top)
	combine.apply_pair (apply_type, ncols_C, ncols_Q, &Q_cur[0], ldq, 
			    &tau_cur[0], C_mine, ldc_mine, C_other, ldc_other, 
			    &work[0]);
      else if (P_mine == P_bot)
	combine.apply_pair (apply_type, ncols_C, ncols_Q, &Q_cur[0], ldq, 
			    &tau_cur[0], C_other, ldc_other, C_mine, ldc_mine, 
			    &work[0]);
      else
	{
	  ostringstream os;
	  os << "Should never get here: P_mine (= " << P_mine 
	     << ") not one of P_top, P_bot = " << P_top << ", " << P_bot;
	  throw std::logic_error (os.str());
	}
    }

    void
    apply_helper (const ApplyType& apply_type,
		  const LocalOrdinal ncols_C,
		  const LocalOrdinal ncols_Q,
		  Scalar C_mine[],
		  const LocalOrdinal ldc_mine,
		  Scalar C_other[], // contiguous ncols_C x ncols_C scratch
		  const LocalOrdinal my_rank,
		  const LocalOrdinal P_first,
		  const LocalOrdinal P_last,
		  const LocalOrdinal tag,
		  MessengerBase< Scalar >* const messenger,
		  const std::vector< std::vector< Scalar > >& Q_factors,
		  const std::vector< std::vector< Scalar > >& tau_arrays,
		  const LocalOrdinal cur_pos,
		  std::vector< Scalar >& work)
    {
      using std::endl;
      using std::ostringstream;
      using std::vector;

      if (P_last <= P_first)
	return;
      else
	{
	  const int P = P_last - P_first + 1;
	  // Whether the interval [P_first, P_last] has an even number of
	  // elements.  Our interval splitting scheme ensures that the
	  // interval [P_first, P_mid - 1] always has an even number of
	  // elements.
	  const bool b_even = (P % 2 == 0);
	  // We split the interval [P_first, P_last] into 2 intervals:
	  // [P_first, P_mid-1], and [P_mid, P_last].  We bias the
	  // splitting procedure so that the lower interval always has an
	  // even number of processor ranks, and never has fewer processor
	  // ranks than the higher interval.
	  const int P_mid = b_even ? (P_first + P/2) : (P_first + P/2 + 1); 

	  if (my_rank < P_mid) // Interval [P_first, P_mid - 1]
	    {
	      const bool b_participating = b_even || my_rank < P_mid - 1;

	      if (cur_pos < 0)
		{
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
	      if (b_even || my_rank < P_mid - 1)
		{
		  if (! b_participating)
		    throw std::logic_error("Should never get here");

		  const int my_offset = my_rank - P_first;
		  const int P_other = P_mid + my_offset;
		  // assert (P_mid <= P_other && P_other <= P_last);
		  if (P_other < P_mid || P_other > P_last)
		    throw std::logic_error("Should never get here");

		  apply_pair (apply_type, ncols_C, ncols_Q, C_mine, ldc_mine, 
			      C_other, my_rank, P_other, tag, messenger, 
			      Q_factors[cur_pos], tau_arrays[cur_pos], work);
		  new_cur_pos = cur_pos - 1;
		}
	      else
		{
		  if (b_participating)
		    throw std::logic_error("Should never get here");

		  new_cur_pos = cur_pos;
		}
	      apply_helper (apply_type, ncols_C, ncols_Q, C_mine, ldc_mine, 
			    C_other, my_rank, P_first, P_mid - 1, tag + 1, 
			    messenger, Q_factors, tau_arrays, new_cur_pos, 
			    work);
	    }
	  else
	    {
	      if (cur_pos < 0)
		{
		  ostringstream os;
		  os << "On Proc " << my_rank << ": cur_pos (= " << cur_pos 
		     << ") < 0; upper interval [" << P_mid << "," << P_last
		     << "]; original interval [" << P_first << "," << P_last 
		     << "]" << endl;
		  throw std::logic_error (os.str());
		}

	      const int my_offset = my_rank - P_mid;
	      const int P_other = P_first + my_offset;
	      // assert (0 <= P_other && P_other < P_mid);
	      apply_pair (apply_type, ncols_C, ncols_Q, C_mine, ldc_mine, 
			  C_other, my_rank, P_other, tag, messenger, 
			  Q_factors[cur_pos], tau_arrays[cur_pos], work);
	      apply_helper (apply_type, ncols_C, ncols_Q, C_mine, ldc_mine, 
			    C_other, my_rank, P_mid, P_last, tag + 1, 
			    messenger, Q_factors, tau_arrays, cur_pos - 1, 
			    work);
	    }
	}
    }
  };

} // namespace TSQR

#endif // __TSQR_Tsqr_DistTsqrHelper_hpp
