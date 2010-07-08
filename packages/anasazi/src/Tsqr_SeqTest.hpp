#ifndef __TSQR_Test_SeqTest_hpp
#define __TSQR_Test_SeqTest_hpp

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
      else if (d_lwork_geqrf < 0)
	{
	  ostringstream os;
	  os << "LAPACK _GEQRF workspace size query returned negative LWORK = " 
	     << d_lwork_geqrf;
	  // It's a logic error and not a runtime error, because no
	  // LAPACK LWORK query should ever return a negative LWORK
	  // value.
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
      else if (d_lwork_orgqr < 0)
	{
	  ostringstream os;
	  os << "LAPACK _ORGQR workspace size query returned negative LWORK = " 
	     << d_lwork_orgqr;
	  // It's a logic error and not a runtime error, because no
	  // LAPACK LWORK query should ever return a negative LWORK
	  // value.
	  throw std::logic_error (os.str());
	}

      // LAPACK workspace queries do return their results as a
      // double-precision floating-point value, but LAPACK promises
      // that that value will fit in an int.  Thus, we don't need to
      // check for valid casts to int below.  I include the checks
      // just to be "bulletproof" and also to show how to do the
      // checks for later reference.
      const Scalar lwork_geqrf_test = static_cast<Scalar >(static_cast<Ordinal >(d_lwork_geqrf));
      if (lwork_geqrf_test != d_lwork_geqrf)
	{
	  ostringstream os;
	  os << "LAPACK _GEQRF workspace query returned a result, " 
	     << d_lwork_geqrf << ", bigger than the max Ordinal value, " 
	     << std::numeric_limits< Ordinal >::max();
	  throw std::range_error (os.str());
	}
      const Scalar lwork_orgqr_test = static_cast<Scalar >(static_cast<Ordinal >(d_lwork_orgqr));
      if (lwork_orgqr_test != d_lwork_orgqr)
	{
	  ostringstream os;
	  os << "LAPACK _ORGQR workspace query returned a result, " 
	     << d_lwork_orgqr << ", bigger than the max Ordinal value, " 
	     << std::numeric_limits< Ordinal >::max();
	  throw std::range_error (os.str());
	}
      return std::max (static_cast< Ordinal > (d_lwork_geqrf), 
		       static_cast< Ordinal > (d_lwork_orgqr));
    }

    /// Test the accuracy of sequential TSQR on an nrows by ncols
    /// matrix (using the given cache block size (in bytes)), and
    /// print the results to stdout.
    template< class Ordinal, class Scalar, class Generator >
    void
    verifySeqTsqr (Generator& generator,
		   const Ordinal nrows, 
		   const Ordinal ncols, 
		   const size_t cache_block_size,
		   const bool contiguous_cache_blocks,
		   const bool human_readable,
		   const bool b_debug = false)
    {
      typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
      using std::cerr;
      using std::cout;
      using std::endl;

      SequentialTsqr< Ordinal, Scalar > actor (cache_block_size);

      if (b_debug)
	{
	  cerr << "Sequential TSQR test problem:" << endl
	       << "* " << nrows << " x " << ncols << endl
	       << "* Cache block of " << actor.cache_block_size() << " bytes" << endl;
	  if (contiguous_cache_blocks)
	    cerr << "* Contiguous cache blocks" << endl;
	}

      Matrix< Ordinal, Scalar > A (nrows, ncols);
      Matrix< Ordinal, Scalar > A_copy (nrows, ncols);
      Matrix< Ordinal, Scalar > Q (nrows, ncols);
      Matrix< Ordinal, Scalar > R (ncols, ncols);
      if (std::numeric_limits< Scalar >::has_quiet_NaN)
	{
	  A.fill (std::numeric_limits< Scalar>::quiet_NaN());
	  A_copy.fill (std::numeric_limits< Scalar >::quiet_NaN());
	  Q.fill (std::numeric_limits< Scalar >::quiet_NaN());
	  R.fill (std::numeric_limits< Scalar >::quiet_NaN());
	}
      const Ordinal lda = nrows;
      const Ordinal ldq = nrows;
      const Ordinal ldr = ncols;

      // Create a test problem
      nodeTestProblem (generator, nrows, ncols, A.get(), A.lda(), true);

      if (b_debug)
	cerr << "-- Generated test problem" << endl;

      // Copy A into A_copy, since TSQR overwrites the input.  If
      // specified, rearrange the data in A_copy so that the data in
      // each cache block is contiguously stored.	  
      if (! contiguous_cache_blocks)
	{
	  A_copy.copy (A);
	  if (b_debug)
	    cerr << "-- Copied test problem from A into A_copy" << endl;
	}
      else
	{
	  actor.cache_block (nrows, ncols, A_copy.get(), A.get(), A.lda());
	  if (b_debug)
	    cerr << "-- Reorganized test matrix to have contiguous "
	      "cache blocks" << endl;

	  // Verify cache blocking, when in debug mode.
	  if (b_debug)
	    {
	      Matrix< Ordinal, Scalar > A2 (nrows, ncols);
	      if (std::numeric_limits< Scalar >::has_quiet_NaN)
		A2.fill (std::numeric_limits< Scalar >::quiet_NaN());

	      actor.un_cache_block (nrows, ncols, A2.get(), A2.lda(), A_copy.get());
	      if (A == A2)
		{
		  if (b_debug)
		    cerr << "-- Cache blocking test succeeded!" << endl;
		}
	      else
		throw std::logic_error ("Cache blocking failed");
	    }
	}

      // Fill R with zeros, since the factorization may not overwrite
      // the strict lower triangle of R.
      R.fill (Scalar(0));

      // Factor the matrix and compute the explicit Q factor
      typedef typename SequentialTsqr< Ordinal, Scalar >::FactorOutput factor_output_type;
      factor_output_type factor_output = 
	actor.factor (nrows, ncols, A_copy.get(), A_copy.lda(), 
		      R.get(), R.lda(), contiguous_cache_blocks);
      if (b_debug)
	cerr << "-- Finished SequentialTsqr::factor" << endl;

      actor.explicit_Q (nrows, ncols, A_copy.get(), lda, factor_output,
			ncols, Q.get(), Q.lda(), contiguous_cache_blocks);
      if (b_debug)
	cerr << "-- Finished SequentialTsqr::explicit_Q" << endl;

      // "Un"-cache-block the output, if contiguous cache blocks were
      // used.  This is only necessary because local_verify() doesn't
      // currently support contiguous cache blocks.
      if (contiguous_cache_blocks)
	{
	  // Use A_copy as temporary storage for un-cache-blocking Q.
	  actor.un_cache_block (nrows, ncols, A_copy.get(), A_copy.lda(), Q.get());
	  Q.copy (A_copy);
	  if (b_debug)
	    cerr << "-- Un-cache-blocked output Q factor" << endl;
	}

      // Print out the R factor
      if (b_debug)
	{
	  cerr << endl << "-- R factor:" << endl;
	  print_local_matrix (cerr, ncols, ncols, R.get(), R.lda());
	  cerr << endl;
	}

      // Validate the factorization
      std::pair< magnitude_type, magnitude_type > results =
	local_verify (nrows, ncols, A.get(), lda, Q.get(), ldq, R.get(), ldr);
      if (b_debug)
	cerr << "-- Finished local_verify" << endl;

      // Print the results
      if (human_readable)
	cout << "Sequential (cache-blocked) TSQR:" << endl
	     << "Relative residual: " << results.first << endl
	     << "Relative orthogonality: " << results.second 
	     << endl << endl;
      else
	cout << "SeqTSQR"
	     << "," << nrows
	     << "," << ncols
	     << "," << actor.cache_block_size()
	     << "," << contiguous_cache_blocks 
	     << "," << results.first 
	     << "," << results.second
	     << endl;
    }

    /// Test the accuracy of LAPACK's QR factorization on an nrows by
    /// ncols matrix, and print the results to stdout.
    template< class Ordinal, class Scalar, class Generator >
    void
    verifyLapack (Generator& generator,
		  const Ordinal nrows, 
		  const Ordinal ncols, 
		  const bool human_readable,
		  const bool b_debug = false)
    {
      typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
      using std::ostringstream;
      using std::cerr;
      using std::cout;
      using std::endl;

      // Initialize LAPACK.
      LAPACK< Ordinal, Scalar > lapack;

      if (b_debug)
	cerr << "LAPACK test problem:" << endl
	     << "* " << nrows << " x " << ncols << endl;

      Matrix< Ordinal, Scalar > A (nrows, ncols);
      Matrix< Ordinal, Scalar > A_copy (nrows, ncols);
      Matrix< Ordinal, Scalar > Q (nrows, ncols);
      Matrix< Ordinal, Scalar > R (ncols, ncols);
      if (std::numeric_limits< Scalar >::has_quiet_NaN)
	{
	  A.fill (std::numeric_limits< Scalar>::quiet_NaN());
	  A_copy.fill (std::numeric_limits< Scalar >::quiet_NaN());
	  Q.fill (std::numeric_limits< Scalar >::quiet_NaN());
	  R.fill (std::numeric_limits< Scalar >::quiet_NaN());
	}
      const Ordinal lda = nrows;
      const Ordinal ldq = nrows;
      const Ordinal ldr = ncols;

      // Create a test problem
      nodeTestProblem (generator, nrows, ncols, A.get(), A.lda(), true);

      if (b_debug)
	cerr << "-- Generated test problem" << endl;

      // Copy A into A_copy, since LAPACK QR overwrites the input.
      A_copy.copy (A);
      if (b_debug)
	cerr << "-- Copied test problem from A into A_copy" << endl;

      // Now determine the required workspace for the factorization.
      const Ordinal lwork = lworkQueryLapackQr (lapack, nrows, ncols, A_copy.lda());
      std::vector< Scalar > work (lwork);
      std::vector< Scalar > tau (ncols);

      // Fill R with zeros, since the factorization may not overwrite
      // the strict lower triangle of R.
      R.fill (Scalar(0));

      // Compute the QR factorization
      int info = 0; // INFO is always an int
      lapack.GEQRF (nrows, ncols, A_copy.get(), A_copy.lda(), 
		    &tau[0], &work[0], lwork, &info);
      if (info != 0)
	{
	  ostringstream os;
	  os << "LAPACK QR factorization (_GEQRF) failed: INFO = " << info;
	  throw std::runtime_error (os.str());
	}

      // Copy out the R factor from A_copy (where we computed the QR
      // factorization in place) into R.
      copy_upper_triangle (ncols, ncols, R.get(), ldr, A_copy.get(), lda);

      // The explicit Q factor will be computed in place, so copy the
      // result of the factorization into Q.
      Q.copy (A_copy);

      // Compute the explicit Q factor
      lapack.ORGQR (nrows, ncols, ncols, Q.get(), ldq, &tau[0], &work[0], lwork, &info);
      if (info != 0)
	{
	  ostringstream os;
	  os << "LAPACK explicit Q computation (_ORGQR) failed: INFO = " << info;
	  throw std::runtime_error (os.str());
	}
  
      // Validate the factorization
      std::pair< magnitude_type, magnitude_type > results = 
	local_verify (nrows, ncols, A.get(), lda, Q.get(), ldq, R.get(), ldr);

      // Print the results
      if (human_readable)
	cout << "LAPACK QR (DGEQRF and DORGQR):" << endl
	     << "Relative residual: " << results.first << endl
	     << "Relative orthogonality: " << results.second << endl;
      else
	cout << "LAPACK"
	     << "," << nrows
	     << "," << ncols
	     << "," << size_t(0) // cache_block_size
	     << "," << false     // contiguous_cache_blocks 
	     << "," << results.first 
	     << "," << results.second
	     << endl;
    }

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
      const bool transposed = false;

      Matrix< Ordinal, Scalar > A (nrows, ncols);
      Matrix< Ordinal, Scalar > A_copy (nrows, ncols);
      Matrix< Ordinal, Scalar > Q (nrows, ncols);
      Matrix< Ordinal, Scalar > R (ncols, ncols);
      const Ordinal lda = nrows;
      const Ordinal ldq = nrows;
      const Ordinal ldr = ncols;

      // Create a test problem
      nodeTestProblem (generator, nrows, ncols, A.get(), lda, false);

      // Copy A into A_copy, since TSQR overwrites the input
      A_copy.copy (A);

      // Benchmark sequential TSQR for ntrials trials
      TimerType timer;
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
      const double seq_tsqr_timing = timer.finish();

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
      TimerType timer;
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
      const double lapack_timing = timer.finish();

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
