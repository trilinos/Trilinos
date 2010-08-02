// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#ifndef __TSQR_Test_SeqTest_hpp
#define __TSQR_Test_SeqTest_hpp

#include <Tsqr_Random_NormalGenerator.hpp>
#include <Tsqr_nodeTestProblem.hpp>
#include <Tsqr_verifyTimerConcept.hpp>

#include <Tsqr_Blas.hpp>
#include <Tsqr_Lapack.hpp>
#include <Tsqr_LocalVerify.hpp>
#include <Tsqr_Matrix.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <Tsqr_SequentialTsqr.hpp>
#include <Tsqr_Util.hpp>

#include <algorithm>
#include <cstring> // size_t definition
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Test {

    template< class Ordinal, class Scalar >
    static Ordinal
    lworkQueryLapackQr (LAPACK< Ordinal, Scalar >& lapack,
			const Ordinal nrows,
			const Ordinal ncols,
			const Ordinal lda)
    {
      typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
      using std::ostringstream;
      using std::endl;

      Scalar d_lwork_geqrf = Scalar (0);
      int INFO = 0;
      lapack.GEQRF (nrows, ncols, NULL, lda, NULL, &d_lwork_geqrf, -1, &INFO);
      if (INFO != 0)
	{
	  ostringstream os;
	  os << "LAPACK _GEQRF workspace size query failed: INFO = " << INFO;
	  // It's a logic error and not a runtime error, because the
	  // LWORK query should only fail if the input parameters have
	  // invalid (e.g., out of range) values.
	  throw std::logic_error (os.str());
	}

      Scalar d_lwork_orgqr = Scalar (0);
      // A workspace query appropriate for computing the explicit Q
      // factor (nrows x ncols) in place, from the QR factorization of
      // an nrows x ncols matrix with leading dimension lda.
      lapack.ORGQR (nrows, ncols, ncols, NULL, lda, NULL, &d_lwork_orgqr, -1, &INFO);
      if (INFO != 0)
	{
	  ostringstream os;
	  os << "LAPACK _ORGQR workspace size query failed: INFO = " << INFO;
	  // It's a logic error and not a runtime error, because the
	  // LWORK query should only fail if the input parameters have
	  // invalid (e.g., out of range) values.
	  throw std::logic_error (os.str());
	}

      // LAPACK workspace queries do return their results as a
      // double-precision floating-point value, but LAPACK promises
      // that that value will fit in an int.  Thus, we don't need to
      // check for valid casts to int below.  I include the checks
      // just to be "bulletproof" and also to show how to do the
      // checks for later reference.
      const magnitude_type lwork_geqrf_test = 
	static_cast< magnitude_type > (static_cast< Ordinal > (ScalarTraits< Scalar >::abs (d_lwork_geqrf)));
      if (lwork_geqrf_test != ScalarTraits< Scalar >::abs (d_lwork_geqrf))
	{
	  ostringstream os;
	  os << "LAPACK _GEQRF workspace query returned a result, " 
	     << d_lwork_geqrf << ", bigger than the max Ordinal value, " 
	     << std::numeric_limits< Ordinal >::max();
	  throw std::range_error (os.str());
	}
      const Scalar lwork_orgqr_test = 
	static_cast< magnitude_type > (static_cast< Ordinal > (ScalarTraits< Scalar >::abs ((d_lwork_orgqr))));
      if (lwork_orgqr_test != ScalarTraits< Scalar >::abs (d_lwork_orgqr))
	{
	  ostringstream os;
	  os << "LAPACK _ORGQR workspace query returned a result, " 
	     << d_lwork_orgqr << ", bigger than the max Ordinal value, " 
	     << std::numeric_limits< Ordinal >::max();
	  throw std::range_error (os.str());
	}
      return std::max (static_cast< Ordinal > (ScalarTraits< Scalar >::abs (d_lwork_geqrf)),
		       static_cast< Ordinal > (ScalarTraits< Scalar >::abs (d_lwork_orgqr)));
    }


    /// Test the accuracy of sequential TSQR on an nrows by ncols
    /// matrix (using the given cache block size (in bytes)), and
    /// print the results to stdout.
    void
    verifySeqTsqr (const int nrows, 
		   const int ncols, 
		   const size_t cache_block_size,
		   const bool test_complex_arithmetic,
		   const bool save_matrices,
		   const bool contiguous_cache_blocks,
		   const bool human_readable,
		   const bool b_debug = false);


    /// Test the accuracy of LAPACK's QR factorization on an nrows by
    /// ncols matrix, and print the results to stdout.
    void
    verifyLapack (const int nrows, 
		  const int ncols, 
		  const bool test_complex_arithmetic,
		  const bool human_readable,
		  const bool b_debug = false);
        

    /// Test the runtime (over ntrials trials) of sequential TSQR, on
    /// an nrows by ncols matrix (using the given cache block size (in
    /// bytes)), and print the results to stdout.
    ///
    /// \param human_readable [in] If true, print the benchmark
    /// results to stdout in human-readable format.  Otherwise, print
    /// them as two rows of comma-delimited ASCII, in an abbreviated
    /// format suitable for automatic processing.
    template< class Ordinal, class Scalar, class Generator, class TimerType >
    void
    benchmarkSeqTsqr (Generator& generator,
		      const Ordinal ntrials,
		      const Ordinal nrows, 
		      const Ordinal ncols, 
		      const size_t cache_block_size,
		      const bool contiguous_cache_blocks,
		      const bool human_readable)
    {
      using std::cerr;
      using std::cout;
      using std::endl;

      TSQR::Test::verifyTimerConcept< TimerType >();

      SequentialTsqr< Ordinal, Scalar > actor (cache_block_size);

      Matrix< Ordinal, Scalar > A (nrows, ncols);
      Matrix< Ordinal, Scalar > A_copy (nrows, ncols);
      Matrix< Ordinal, Scalar > Q (nrows, ncols);
      Matrix< Ordinal, Scalar > R (ncols, ncols);
      const Ordinal lda = nrows;
      const Ordinal ldq = nrows;

      // Create a test problem
      nodeTestProblem (generator, nrows, ncols, A.get(), lda, false);

      // Copy A into A_copy, since TSQR overwrites the input
      A_copy.copy (A);

      // Benchmark sequential TSQR for ntrials trials.
      //
      // Name of timer doesn't matter here; we only need the timing.
      TimerType timer("SeqTSQR");
      timer.start();
      for (int trial_num = 0; trial_num < ntrials; ++trial_num)
	{
	  // Factor the matrix and extract the resulting R factor
	  typedef typename SequentialTsqr< Ordinal, Scalar >::FactorOutput factor_output_type;
	  factor_output_type factor_output = 
	    actor.factor (nrows, ncols, A_copy.get(), lda, R.get(), R.lda(),
			  contiguous_cache_blocks);
	  // Compute the explicit Q factor.  Unlike with LAPACK QR,
	  // this doesn't happen in place: the implicit Q factor is
	  // stored in A_copy, and the explicit Q factor is written to
	  // Q.
	  actor.explicit_Q (nrows, ncols, A_copy.get(), lda, factor_output, 
			    ncols, Q.get(), ldq, contiguous_cache_blocks);
	}
      const double seq_tsqr_timing = timer.stop();

      // Print the results  
      if (human_readable)
	cout << "Sequential (cache-blocked) TSQR:" << endl
	     << "# rows = " << nrows << endl
	     << "# columns = " << ncols << endl
	     << "cache block # bytes = " << actor.cache_block_size() << endl
	     << "contiguous cache blocks? " << contiguous_cache_blocks << endl
	     << "# trials = " << ntrials << endl
	     << "Total time (s) = " << seq_tsqr_timing << endl 
	     << endl;
      else
	cout << "SeqTSQR" 
	     << "," << nrows 
	     << "," << ncols 
	     << "," << actor.cache_block_size()
	     << "," << contiguous_cache_blocks
	     << "," << ntrials << "," 
	     << seq_tsqr_timing << endl;
    }


    /// Test the runtime (over ntrials trials) of LAPACK QR, on an
    /// nrows by ncols matrix, and print the results to stdout.
    ///
    /// \param human_readable [in] If true, print the benchmark
    /// results to stdout in human-readable format.  Otherwise, print
    /// them as two rows of comma-delimited ASCII, in an abbreviated
    /// format suitable for automatic processing.
    template< class Ordinal, class Scalar, class Generator, class TimerType >
    void
    benchmarkLapack (Generator& generator,
		     const int ntrials,
		     const Ordinal nrows, 
		     const Ordinal ncols, 
		     const bool human_readable)
    {
      using std::ostringstream;
      using std::cerr;
      using std::cout;
      using std::endl;

      TSQR::Test::verifyTimerConcept< TimerType >();

      LAPACK< Ordinal, Scalar > lapack;
      Matrix< Ordinal, Scalar > A (nrows, ncols);
      Matrix< Ordinal, Scalar > Q (nrows, ncols);
      Matrix< Ordinal, Scalar > R (ncols, ncols);
      const Ordinal lda = nrows;
      const Ordinal ldq = nrows;
      const Ordinal ldr = ncols;

      // Create a test problem
      nodeTestProblem (generator, nrows, ncols, A.get(), lda, false);

      // Copy A into Q, since LAPACK QR overwrites the input.  We only
      // need Q because LAPACK's computation of the explicit Q factor
      // occurs in place.  This doesn't work with TSQR.  To give
      // LAPACK QR the fullest possible advantage over TSQR, we don't
      // allocate an A_copy here (as we would when benchmarking TSQR).
      Q.copy (A);

      // Determine the required workspace for the factorization
      const Ordinal lwork = lworkQueryLapackQr (lapack, nrows, ncols, lda);
      std::vector< Scalar > work (lwork);
      std::vector< Scalar > tau (ncols);

      // Benchmark LAPACK's QR factorization for ntrials trials
      TimerType timer("LapackQR");
      timer.start();
      for (int trial_num = 0; trial_num < ntrials; ++trial_num)
	{
	  // Compute the QR factorization
	  int info = 0; // INFO is always an int
	  lapack.GEQRF (nrows, ncols, Q.get(), ldq, &tau[0], &work[0], lwork, &info);
	  if (info != 0)
	    {
	      ostringstream os;
	      os << "LAPACK QR factorization (_GEQRF) failed: INFO = " << info;
	      throw std::runtime_error (os.str());
	    }

	  // Extract the upper triangular factor R from Q (where it
	  // was computed in place by GEQRF), since ORGQR will
	  // overwrite all of Q with the explicit Q factor.
	  copy_upper_triangle (nrows, ncols, R.get(), ldr, Q.get(), ldq);

	  // Compute the explicit Q factor
	  lapack.ORGQR (nrows, ncols, ncols, Q.get(), ldq,
			&tau[0], &work[0], lwork, &info);
	  if (info != 0)
	    {
	      ostringstream os;
	      os << "LAPACK explicit Q computation (_ORGQR) failed: INFO = " << info;
	      throw std::runtime_error (os.str());
	    }
	}
      const double lapack_timing = timer.stop();

      // Print the results  
      if (human_readable)
	cout << "LAPACK\'s QR factorization (DGEQRF + DORGQR):" << endl
	     << "nrows = " << nrows << endl
	     << "ncols = " << ncols << endl
	     << "ntrials = " << ntrials << endl
	     << "Total time (s) = " << lapack_timing << endl << endl;
      else
	// "0" refers to the cache block size, which is not applicable
	// in this case.
	cout << "LAPACK" 
	     << "," << nrows 
	     << "," << ncols
	     << "," << 0 
	     << "," << ntrials 
	     << "," << lapack_timing
	     << endl;
    }

  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_SeqTest_hpp
