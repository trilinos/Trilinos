#ifndef __TSQR_DistTsqrRB_hpp
#define __TSQR_DistTsqrRB_hpp

#include <Tsqr_ApplyType.hpp>
#include <Tsqr_Combine.hpp>
#include <Tsqr_Matrix.hpp>
#include <Tsqr_ScalarTraits.hpp>

#include <Teuchos_RCP.hpp>

#include <algorithm>
#include <sstream>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class DistTsqrRB
  /// \brief Reduce-and-Broadcast (RB) version of DistTsqr
  ///
  /// Reduce-and-Broadcast (RB) version of DistTsqr, which factors a
  /// vertical stack of n by n R factors, one per MPI process.  Only
  /// the final R factor is broadcast; the implicit Q factor data stay
  /// on the MPI process where they are computed.
  template< class LocalOrdinal, class Scalar >
  class DistTsqrRB {
  public:
    typedef LocalOrdinal ordinal_type;
    typedef Scalar scalar_type;
    typedef ScalarTraits< Scalar >::magnitude_type magnitude_type;
    typedef MatView< LocalOrdinal, Scalar > matview_type;
    typedef Matrix< LocalOrdinal, Scalar > matrix_type;
    typedef int rank_type;
    typedef Combine< LocalOrdinal, Scalar > combine_type;

    /// Constructor
    ///
    /// \param messenger [in/out] Smart pointer to a wrapper handling
    /// communication between MPI process(es).
    DistTsqrRB (const Teuchos::RCP< MessengerBase< Scalar > >& messenger) :
      messenger_ (messenger)
    {}

    /// Whether or not all diagonal entries of the R factor computed
    /// by the QR factorization are guaranteed to be nonnegative.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      return Combine< LocalOrdinal, Scalar >::QR_produces_R_factor_with_nonnegative_diagonal();
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
    void
    factorExplicit (matview_type& R_mine, matview_type& Q_mine)
    {
      // Dimension sanity checks.  R_mine should have at least as many
      // rows as columns (since we will be working on the upper
      // triangle).  Q_mine should have the same number of rows as
      // R_mine has columns, but Q_mine may have any number of
      // columns.  (It depends on how many columns of the explicit Q
      // factor we want to compute.)
      if (R_mine.nrows() < R_mine.ncols())
	{
	  std::ostringstream os;
	  os << "R factor input has fewer rows (" << R_mine.nrows() 
	     << ") than columns (" << R_mine.ncols() << ")";
	  // This is a logic error because TSQR users should not be
	  // calling this method directly.
	  throw std::logic_error (os.str());
	}
      else if (Q_mine.nrows() != R_mine.ncols())
	{
	  std::ostringstream os;
	  os << "Q factor input must have the same number of rows as the R "
	    "factor input has columns.  Q has " << Q_mine.nrows() 
	     << " rows, but R has " << R_mine.ncols() << " columns.";
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
      factorReduce (R_mine, P_mine, P_first, P_last, QFactors, tauArrays);
      if (QFactors.size() != tauArrays.size())
	{
	  std::ostringstream os;
	  os << "QFactors and tauArrays should have the same number of element"
	    "s after factorReduce() returns, but they do not.  QFactors has " 
	     << QFactors.size() << " elements, but tauArrays has " 
	     << tauArrays.size() << " elements.";
	  throw std::logic_error (os.str());
	}

      Q_mine.fill (scalar_type (0));
      if (messenger_->rank() == 0)
	{
	  for (ordinal_type j = 0; j < Q_mine.ncols(); ++j)
	    Q_mine(j, j) = scalar_type (1);
	}
      // Scratch space for computing results to send to other processors.
      matrix_type Q_other (Q_mine.nrows(), Q_mine.ncols(), scalar_type (0));
      const rank_type numSteps = QFactors.size() - 1;
      explicitQBroadcast (R_mine, Q_mine, Q_other.view(), 
			  P_mine, P_first, P_last,
			  numSteps, QFactors, tauArrays);
    }

  private:

    void
    factorReduce (matview_type& R_mine,
		  const rank_type P_mine, 
		  const rank_type P_first,
		  const rank_type P_last,
		  std::vector< matrix_type >& QFactors,
		  std::vector< std::vector< scalar_type > >& tauArrays)
    {
      if (P_last < P_first)
	{
	  std::ostringstream os;
	  os << "Programming error in factorReduce() recursion: interval "
	    "[P_first, P_last] is invalid: P_first = " << P_first 
	     << ", P_last = " << P_last << ".";
	  throw std::logic_error (os.str());
	}
      else if (P_mine < P_first || P_mine > P_last)
	{
	  std::ostringstream os;
	  os << "Programming error in factorReduce() recursion: P_mine (= " 
	     << P_mine << ") is not in current process rank interval " 
	     << "[P_first = " << P_first << ", P_last = " << P_last << "]";
	  throw std::logic_error (os.str());
	}
      else if (P_last == P_first)
	return; // skip singleton intervals (see explanation below)
      else
	{
	  // Recurse on two intervals: [P_first, P_mid-1] and [P_mid,
	  // P_last].  For example, if [P_first, P_last] = [0, 9],
	  // P_mid = floor( (0+9+1)/2 ) = 5 and the intervals are
	  // [0,4] and [5,9].  
	  // 
	  // If [P_first, P_last] = [4,6], P_mid = floor( (4+6+1)/2 )
	  // = 5 and the intervals are [4,4] (a singleton) and [5,6].
	  // The latter case shows that singleton intervals may arise.
	  // We treat them as a base case in the recursion.  Process 4
	  // won't be skipped completely, though; it will get combined
	  // with the result from [5,6].

	  // Adding 1 and doing integer division works like "ceiling."
	  const rank_type P_mid = (P_first + P_last + 1) / 2;

	  if (P_mine < P_mid) // Interval [P_first, P_mid-1]
	    factorReduce (R_mine, P_mine, P_first, P_mid - 1,
			  QFactors, tauArrays);
	  else // Interval [P_mid, P_last]
	    factorReduce (R_mine, P_mine, P_mid, P_last,
			  QFactors, tauArrays);

	  // This only does anything if P_mine is either P_first or P_mid.
	  if (P_mine == P_first)
	    {
	      const ordinal_type numCols = R_mine.ncols();
	      matrix_type R_other (numCols, numCols);
	      recv_R (R_other.view(), P_mid);

	      std::vector< scalar_type > tau (numCols);
	      // Don't shrink the workspace array; doing so may
	      // require expensive reallocation every time we send /
	      // receive data.
	      work_.resize (std::max (work_.size(), numCols));
	      combine_.factor_pair (numCols, R_mine.get(), R_mine.lda(), 
				    R_other.get(), R_other.lda(), 
				    &tau[0], &work_[0]);
	      QFactors.push_back (R_other);
	      tauArrays.push_back (tau);
	    }
	  else if (P_mine == P_mid)
	    send_R (R_mine, P_first);
	}
    }

    void
    explicitQBroadcast (matview_type& R_mine,
			matview_type& Q_mine,
			matview_type& Q_other, // workspace
			const rank_type P_mine, 
			const rank_type P_first,
			const rank_type P_last,
			const rank_type curpos,
			std::vector< matrix_type >& QFactors,
			std::vector< std::vector< scalar_type > >& tauArrays)
    {
      if (P_last < P_first)
	{
	  std::ostringstream os;
	  os << "Programming error in explicitQBroadcast() recursion: interval"
	    " [P_first, P_last] is invalid: P_first = " << P_first 
	     << ", P_last = " << P_last << ".";
	  throw std::logic_error (os.str());
	}
      else if (P_mine < P_first || P_mine > P_last)
	{
	  std::ostringstream os;
	  os << "Programming error in explicitQBroadcast() recursion: P_mine "
	    "(= " << P_mine << ") is not in current process rank interval " 
	     << "[P_first = " << P_first << ", P_last = " << P_last << "]";
	  throw std::logic_error (os.str());
	}
      else if (P_last == P_first)
	return; // skip singleton intervals
      else
	{
	  // Adding 1 and integer division works like "ceiling."
	  const rank_type P_mid = (P_first + P_last + 1) / 2;

	  if (P_mine == P_first)
	    {
	      if (curpos < 0)
		{
		  std::ostringstream os;
		  os << "Programming error: On the current P_first (= " 
		     << P_first << ") proc: curpos (= " << curpos << ") < 0";
		  throw std::logic_error (os.str());
		}
	      // Q_impl, tau: implicitly stored local Q factor.
	      matrix_type& Q_impl = QFactors[curpos];
	      std::vector< scalar_type >& tau = tauArrays[curpos];
	      
	      // Apply implicitly stored local Q factor to 
	      //   [Q_mine; 
	      //    Q_other]
	      // where Q_other = zeros(Q_mine.nrows(), Q_mine.ncols()).
	      // Overwrite both Q_mine and Q_other with the result.
	      Q_other.fill (scalar_type (0));
	      combine_.apply_pair (applyType, Q_mine.ncols(), Q_impl.ncols(),
				   Q_impl.get(), Q_impl.lda(), &tau[0],
				   Q_mine.get(), Q_mine.lda(),
				   Q_other.get(), Q_other.lda(), &work_[0]);
	      // Send the resulting Q_other, and the final R factor, to P_mid.
	      send_Q_R (Q_other, R_mine, P_mid)
	    }
	  else if (P_mine == P_mid)
	    // P_first computed my explicit Q factor component.
	    // Receive it, and the final R factor, from P_first.
	    receive_Q_R (Q_mine, R_mine, P_first);

	  if (P_mine < P_mid) // Interval [P_first, P_mid-1]
	    explicitQBroadcast (R_mine, Q_mine, Q_other, 
				P_mine, P_first, P_mid - 1,
				curpos - 1, QFactors, tauArrays);
	  else // Interval [P_mid, P_last]
	    explicitQBroadcast (R_mine, Q_mine, Q_other, 
				P_mine, P_mid, P_last,
				curpos - 1, QFactors, tauArrays);
	}
    }

    void
    send_C_R (const matview_type& Q, const matview_type& R, const rank_type destProc) 
    {
      const ordinal_type R_numCols = R.ncols();
      const ordinal_type Q_size = Q.nrows() * Q.ncols();
      const ordinal_type R_size = (R_numCols * (R_numCols + 1)) / 2;
      const ordinal_type numElts = Q_size + R_size;

      // Don't shrink the workspace array; doing so would still be
      // correct, but may require reallocation of data when it needs
      // to grow again.
      work_.resize (std::max (work_.size(), numElts));

      // Pack the Q data into the workspace array.
      matview_type Q_contig (Q.nrows(), Q.ncols(), &work_[0], Q.nrows());
      Q_contig.copy (Q);
      // Pack the R data into the workspace array.
      pack_R (R, &work_[Q_size]);
      messenger_->send (&work_[0], numElts, destProc, 0);
    }

    void
    recv_Q_R (matview_type& Q, matview_type& R, const rank_type srcProc)
    {
      const ordinal_type R_numCols = R.ncols();
      const ordinal_type Q_size = Q.nrows() * Q.ncols();
      const ordinal_type R_size = (R_numCols * (R_numCols + 1)) / 2;
      const ordinal_type numElts = Q_size + R_size;

      // Don't shrink the workspace array; doing so would still be
      // correct, but may require reallocation of data when it needs
      // to grow again.
      work_.resize (std::max (work_.size(), numElts));
      messenger_->recv (&work_[0], numElts, destProc, 0);

      // Unpack the C data from the workspace array.
      Q.copy (matview_type (Q.nrows(), Q.ncols(), &work_[0], Q.nrows()));
      // Unpack the R data from the workspace array.
      unpack_R (R, &work_[Q_size]);
    }

    void
    send_R (const matrix_type& R, const rank_type destProc)
    {
      const ordinal_type numCols = R.ncols();
      const ordinal_type numElts = (numCols * (numCols+1)) / 2;

      // Don't shrink the workspace array; doing so would still be
      // correct, but may require reallocation of data when it needs
      // to grow again.
      work_.resize (std::max (work_.size(), numElts));
      // Pack the R data into the workspace array.
      pack_R (R, &work_[0]);
      messenger_->send (&work_[0], numElts, destProc, 0);
    }

    void
    recv_R (matrix_type& R, const rank_type srcProc)
    {
      const ordinal_type numCols = R.ncols();
      const ordinal_type numElts = (numCols * (numCols+1)) / 2;

      // Don't shrink the workspace array; doing so would still be
      // correct, but may require reallocation of data when it needs
      // to grow again.
      work_.resize (std::max (work_.size(), numElts));
      messenger_->recv (&work_[0], numElts, srcProc, 0);
      // Unpack the R data from the workspace array.
      unpack_R (R, &work_[0]);
    }

    static void 
    unpack_R (matrix_type& R, const scalar_type buf[])
    {
      ordinal_type curpos = 0;
      for (ordinal_type j = 0; j < R.ncols(); ++j)
	{
	  const scalar_type* const R_j = &R(0, j);
	  for (ordinal_type i = 0; i <= j; ++i)
	    R_j[i] = buf[curpos++];
	}
    }

    static void 
    pack_R (const matrix_type& R, scalar_type buf[])
    {
      ordinal_type curpos = 0;
      for (ordinal_type j = 0; j < R.ncols(); ++j)
	{
	  scalar_type* const R_j = &R(0, j);
	  for (ordinal_type i = 0; i <= j; ++i)
	    buf[curpos++] = R_j[i];
	}
    }

  private:
    combine_type combine_;
    Teuchos::RCP< MessengerBase< scalar_type > > messenger_;
    std::vector< scalar_type > work_;
  };

} // namespace TSQR

#endif // __TSQR_DistTsqrRB_hpp
