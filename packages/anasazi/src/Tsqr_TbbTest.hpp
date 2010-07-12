#ifndef __TSQR_Test_TbbTest_hpp
#define __TSQR_Test_TbbTest_hpp

#include <Tsqr_nodeTestProblem.hpp>
#include <Tsqr_verifyTimerConcept.hpp>

#include <Tsqr_Blas.hpp>
#include <Tsqr_Lapack.hpp>
#include <Tsqr_LocalVerify.hpp>
#include <Tsqr_Matrix.hpp>
#include <Tsqr_Util.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <TbbTsqr.hpp>

#include <algorithm>
#include <cstring> // size_t definition
//#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using std::make_pair;
using std::pair;
using std::vector;

using std::cerr;
using std::cout;
using std::endl;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Test {

    /// Test the accuracy of Intel TBB TSQR on an nrows by ncols
    /// matrix (using the given number of cores and the given cache
    /// block size (in bytes)), and print the results to stdout.
    template< class Ordinal, class Scalar, class Generator >
    void
    verifyTbbTsqr (Generator& generator,
		   const Ordinal nrows, 
		   const Ordinal ncols, 
		   const int num_cores,
		   const size_t cache_block_size,
		   const bool contiguous_cache_blocks,
		   const bool human_readable,
		   const bool b_debug = false)
    {
      // Need c++0x to have a default template parameter argument for
      // a template function, otherwise we would have templated this
      // function on TimerType and made TrivialTimer the default.
      // TimerType is only used instead of TbbTsqr.
      typedef TSQR::TBB::TrivialTimer TimerType; 
      typedef TSQR::TBB::TbbTsqr< Ordinal, Scalar, TimerType > node_tsqr_type;
      typedef typename node_tsqr_type::FactorOutput factor_output_type;
      typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
      using std::cerr;
      using std::cout;
      using std::endl;

      node_tsqr_type actor (num_cores, cache_block_size);

      if (b_debug)
	{
	  cerr << "Intel TBB TSQR test problem:" << endl
	       << "* " << nrows << " x " << ncols << endl
	       << "* # cores: " << num_cores << endl
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
      factor_output_type factor_output = 
	actor.factor (nrows, ncols, A_copy.get(), A_copy.lda(), R.get(), 
		      R.lda(), contiguous_cache_blocks);
      if (b_debug)
	cerr << "-- Finished TbbTsqr::factor" << endl;
      actor.explicit_Q (nrows, ncols, A_copy.get(), A_copy.lda(), factor_output,
			ncols, Q.get(), Q.lda(), contiguous_cache_blocks);
      if (b_debug)
	cerr << "-- Finished TbbTsqr::explicit_Q" << endl;

      // "Un"-cache-block the output Q (the explicit Q factor), if
      // contiguous cache blocks were used.  This is only necessary
      // because local_verify() doesn't currently support contiguous
      // cache blocks.
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
	cout << "Parallel (via Intel\'s Threading Building Blocks) / cache-blocked) TSQR:" << endl
	     << "# rows = " << nrows << endl
	     << "# columns = " << ncols << endl
	     << "# cores: " << num_cores << endl
	     << "cache block # bytes = " << actor.cache_block_size() << endl
	     << "contiguous cache blocks? " << contiguous_cache_blocks << endl
	     << "Relative residual $\\|A - Q*R\\|_2 / \\|A\\|_2$ = " 
	     << results.first << endl
	     << "Relative orthogonality $\\|I - Q^T*Q\\|_2$ = " 
	     << results.second << endl
	     << endl;
      else
	cout << "TbbTSQR"
	     << "," << nrows
	     << "," << ncols
	     << "," << num_cores
	     << "," << actor.cache_block_size()
	     << "," << contiguous_cache_blocks 
	     << "," << results.first 
	     << "," << results.second
	     << endl;
    }

    /// Benchmark Intel TBB TSQR vs. LAPACK's QR, and print the
    /// results to stdout.
    template< class Ordinal, class Scalar, class Generator, class TimerType >
    void
    benchmarkTbbTsqr (Generator& generator,
		      const int ntrials,
		      const Ordinal nrows, 
		      const Ordinal ncols, 
		      const int num_cores,
		      const size_t cache_block_size,
		      const bool contiguous_cache_blocks,
		      const bool human_readable)
    {
      typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
      using TSQR::TBB::TbbTsqr;
      using std::cerr;
      using std::cout;
      using std::endl;

      TSQR::Test::verifyTimerConcept< TimerType >();

      TbbTsqr< Ordinal, Scalar, TimerType > actor (num_cores, cache_block_size);

      Matrix< Ordinal, Scalar > A (nrows, ncols);
      Matrix< Ordinal, Scalar > A_copy (nrows, ncols);
      Matrix< Ordinal, Scalar > Q (nrows, ncols);
      Matrix< Ordinal, Scalar > R (ncols, ncols, Scalar(0));

      // Fill R with zeros, since the factorization may not overwrite
      // the strict lower triangle of R.
      R.fill (Scalar(0));

      // Create a test problem
      nodeTestProblem (generator, nrows, ncols, A.get(), A.lda(), false);

      // Copy A into A_copy, since TSQR overwrites the input.  If
      // specified, rearrange the data in A_copy so that the data in
      // each cache block is contiguously stored.	  
      if (contiguous_cache_blocks)
	actor.cache_block (nrows, ncols, A_copy.get(), A.get(), A.lda());
      else
	A_copy.copy (A);

      // Benchmark TBB-based TSQR for ntrials trials.
      //
      // Name of timer doesn't matter here; we only need the timing.
      TimerType timer("TbbTSQR");
      timer.start();
      for (int trial_num = 0; trial_num < ntrials; ++trial_num)
	{
	  // Factor the matrix in-place in A_copy, and extract the
	  // resulting R factor into R.
	  typedef typename TbbTsqr< Ordinal, Scalar, TimerType >::FactorOutput factor_output_type;
	  factor_output_type factor_output = 
	    actor.factor (nrows, ncols, A_copy.get(), A_copy.lda(), 
			  R.get(), R.lda(), contiguous_cache_blocks);
	  // Compute the explicit Q factor (which was stored
	  // implicitly in A_copy and factor_output) and store in Q.
	  // We don't need to un-cache-block the output, because we
	  // aren't verifying it here.
	  actor.explicit_Q (nrows, ncols, A_copy.get(), A_copy.lda(), 
			    factor_output, ncols, Q.get(), Q.lda(), 
			    contiguous_cache_blocks);
	}
      const double tbb_tsqr_timing = timer.stop();

      // Print the results
      if (human_readable)
	{
	  cout << "(Intel TBB / cache-blocked) TSQR:" << endl
	       << "# rows = " << nrows << endl
	       << "# columns = " << ncols << endl
	       << "# cores: " << num_cores << endl
	       << "cache block # bytes = " << actor.cache_block_size() << endl
	       << "contiguous cache blocks? " << contiguous_cache_blocks << endl
	       << "# trials = " << ntrials << endl
	       << "Total time (s) = " << tbb_tsqr_timing << endl
	       << "Total time (s) in factor() (min over all tasks): " 
	       << (ntrials * actor.min_seq_factor_timing()) << endl
	       << "Total time (s) in factor() (max over all tasks): " 
	       << (ntrials * actor.max_seq_factor_timing()) << endl
	       << "Total time (s) in apply() (min over all tasks): " 
	       << (ntrials * actor.min_seq_apply_timing()) << endl
	       << "Total time (s) in apply() (max over all tasks): " 
	       << (ntrials * actor.max_seq_apply_timing()) << endl
	       << endl << endl;
	}
      else
	// We don't include {min,max}_seq_apply_timing() here, because
	// those times don't benefit from the accuracy of benchmarking
	// for ntrials > 1.  Thus, it's misleading to include them
	// with tbb_tsqr_timing, the total time over ntrials trials.
	cout << "TbbTSQR"
	     << "," << nrows
	     << "," << ncols
	     << "," << num_cores
	     << "," << actor.cache_block_size()
	     << "," << contiguous_cache_blocks 
	     << "," << ntrials
	     << "," << tbb_tsqr_timing 
	     << endl;
    }
  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_TbbTest_hpp
