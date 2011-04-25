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

#ifndef __TSQR_Test_TbbTest_hpp
#define __TSQR_Test_TbbTest_hpp

#include <Tsqr_nodeTestProblem.hpp>
#include <Tsqr_verifyTimerConcept.hpp>
#include <Tsqr_Random_NormalGenerator.hpp>

#include <Tsqr_Blas.hpp>
#include <Tsqr_Lapack.hpp>
#include <Tsqr_LocalVerify.hpp>
#include <Tsqr_Matrix.hpp>
#include <Tsqr_Util.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <TbbTsqr.hpp>

#include <Teuchos_Time.hpp>

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
    template< class Ordinal, class Scalar >
    void
    verifyTbbTsqr (const std::string& scalarTypeName,
		   TSQR::Random::NormalGenerator< Ordinal, Scalar >& generator,
		   const Ordinal nrows, 
		   const Ordinal ncols, 
		   const int num_cores,
		   const size_t cache_block_size,
		   const bool contiguous_cache_blocks,
		   const bool printFieldNames,
		   const bool human_readable,
		   const bool b_debug = false)
    {
      typedef Teuchos::Time timer_type;
      typedef TSQR::TBB::TbbTsqr< Ordinal, Scalar, timer_type > node_tsqr_type;
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
      std::vector< magnitude_type > results =
	local_verify (nrows, ncols, A.get(), lda, Q.get(), ldq, R.get(), ldr);
      if (b_debug)
	cerr << "-- Finished local_verify" << endl;

      // Print the results
      if (human_readable)
	cout << "Parallel (via Intel\'s Threading Building Blocks) / cache-blocked) TSQR:" << endl
	     << "Scalar type: " << scalarTypeName << endl
	     << "# rows = " << nrows << endl
	     << "# columns = " << ncols << endl
	     << "# cores: " << num_cores << endl
	     << "cache block # bytes = " << actor.cache_block_size() << endl
	     << "contiguous cache blocks? " << contiguous_cache_blocks << endl
	     << "Absolute residual $\\|A - Q*R\\|_2$: "
	     << results[0] << endl
	     << "Absolute orthogonality $\\|I - Q^T*Q\\|_2$: " 
	     << results[1] << endl
	     << "Test matrix norm $\\| A \\|_F$: "
	     << results[2] << endl
	     << endl;
      else
	{
	  if (printFieldNames)
	    {
	      const char prefix[] = "%";
	      cout << prefix
		   << "method"
		   << ",scalarType"
		   << ",numRows"
		   << ",numCols"
		   << ",numThreads"
		   << ",cacheBlockSize"
		   << ",contiguousCacheBlocks"
		   << ",absFrobResid"
		   << ",absFrobOrthog"
		   << ",frobA"
		   << endl;
	    }
	  cout << "TbbTsqr"
	       << "," << scalarTypeName
	       << "," << nrows
	       << "," << ncols
	       << "," << num_cores
	       << "," << actor.cache_block_size()
	       << "," << contiguous_cache_blocks 
	       << "," << results[0]
	       << "," << results[1]
	       << "," << results[2]
	       << endl;
	}
    }

    /// Benchmark Intel TBB TSQR vs. LAPACK's QR, and print the
    /// results to stdout.
    ///
    /// \note c++0x support is need in order to have a default
    /// template parameter argument for a template function, otherwise
    /// we would have templated this function on TimerType and made
    /// Teuchos::Time the default.
    template< class Ordinal, class Scalar >
    void
    benchmarkTbbTsqr (const std::string& scalarTypeName,
		      const int ntrials,
		      const Ordinal nrows, 
		      const Ordinal ncols, 
		      const int num_cores,
		      const size_t cache_block_size,
		      const bool contiguous_cache_blocks,
		      const bool printFieldNames,
		      const bool human_readable)
    {
      using TSQR::TBB::TbbTsqr;
      using std::cerr;
      using std::cout;
      using std::endl;

      typedef Teuchos::Time timer_type;
      typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
      typedef Ordinal ordinal_type;
      typedef Scalar scalar_type;
      typedef Matrix< ordinal_type, scalar_type > matrix_type;
      typedef TbbTsqr< ordinal_type, scalar_type, timer_type > node_tsqr_type;

      // Pseudorandom normal(0,1) generator.  Default seed is OK,
      // because this is a benchmark, not an accuracy test.
      TSQR::Random::NormalGenerator< ordinal_type, scalar_type > generator;

      // Set up TSQR implementation.
      node_tsqr_type actor (num_cores, cache_block_size);

      matrix_type A (nrows, ncols);
      matrix_type A_copy (nrows, ncols);
      matrix_type Q (nrows, ncols);
      matrix_type R (ncols, ncols, scalar_type(0));

      // Fill R with zeros, since the factorization may not overwrite
      // the strict lower triangle of R.
      R.fill (scalar_type(0));

      // Create a test problem
      nodeTestProblem (generator, nrows, ncols, A.get(), A.lda(), false);

      // Copy A into A_copy, since TSQR overwrites the input.  If
      // specified, rearrange the data in A_copy so that the data in
      // each cache block is contiguously stored.	  
      if (contiguous_cache_blocks)
	actor.cache_block (nrows, ncols, A_copy.get(), A.get(), A.lda());
      else
	A_copy.copy (A);

      // Do a few timing runs and throw away the results, just to warm
      // up any libraries that do autotuning.
      const int numWarmupRuns = 5;
      for (int warmupRun = 0; warmupRun < numWarmupRuns; ++warmupRun)
	{
	  // Factor the matrix in-place in A_copy, and extract the
	  // resulting R factor into R.
	  typedef typename node_tsqr_type::FactorOutput factor_output_type;
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

      // Benchmark TBB-based TSQR for ntrials trials.
      //
      // Name of timer doesn't matter here; we only need the timing.
      timer_type timer("TbbTsqr");
      timer.start();
      for (int trial_num = 0; trial_num < ntrials; ++trial_num)
	{
	  // Factor the matrix in-place in A_copy, and extract the
	  // resulting R factor into R.
	  typedef typename node_tsqr_type::FactorOutput factor_output_type;
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
	  cout << "(Intel TBB / cache-blocked) TSQR cumulative timings:" << endl
	       << "Scalar type: " << scalarTypeName << endl
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
	  cout << "(Intel TBB / cache-blocked) TSQR per-invocation timings:" << endl;
	  
	  std::vector< TimeStats > stats;
	  actor.getStats (stats);
	  std::vector< std::string > labels;
	  actor.getStatsLabels (labels);

	  const std::string labelLabel ("label");
	  for (std::vector< std::string >::size_type k = 0; k < labels.size(); ++k)
	    {
	      const bool printHeaders = (k == 0);
	      if (stats[k].count() > 0)
		stats[k].print (cout, human_readable, labels[k], labelLabel, printHeaders);
	    }
	}
      else
	{
	  if (printFieldNames)
	    {
	      const char prefix[] = "%";
	      cout << prefix 
		   << "method"
		   << ",scalarType"
		   << ",numRows"
		   << ",numCols"
		   << ",numThreads"
		   << ",cacheBlockSize"
		   << ",contiguousCacheBlocks"
		   << ",numTrials"
		   << ",timing"
		   << endl;
	    }

	  // We don't include {min,max}_seq_apply_timing() here, because
	  // those times don't benefit from the accuracy of benchmarking
	  // for ntrials > 1.  Thus, it's misleading to include them
	  // with tbb_tsqr_timing, the total time over ntrials trials.
	  cout << "TbbTsqr"
	       << "," << scalarTypeName
	       << "," << nrows
	       << "," << ncols
	       << "," << num_cores
	       << "," << actor.cache_block_size()
	       << "," << contiguous_cache_blocks 
	       << "," << ntrials
	       << "," << tbb_tsqr_timing 
	       << endl;
	}
    }
  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_TbbTest_hpp
