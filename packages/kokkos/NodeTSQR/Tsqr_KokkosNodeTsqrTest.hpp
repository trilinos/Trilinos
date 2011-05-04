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

#ifndef __TSQR_Test_KokkosNodeTsqrTest_hpp
#define __TSQR_Test_KokkosNodeTsqrTest_hpp

#include <Tsqr_nodeTestProblem.hpp>
#include <Tsqr_verifyTimerConcept.hpp>
#include <Tsqr_Random_NormalGenerator.hpp>
#include <Tsqr_LocalVerify.hpp>
#include <Tsqr_Matrix.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <Tsqr_KokkosNodeTsqr.hpp>

#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace TSQR {
  namespace Test {

    /// \fn verifyKokkosNodeTsqr
    /// \brief Test accuracy of KokkosNodeTsqr's QR factorization.
    ///
    /// Test the accuracy of KokkosNodeTsqr's QR factorization on a
    /// numRows by numCols matrix, and print results to stdout.
    ///
    /// \param node [in] The Kokkos Node instance on which to execute
    ///   in parallel.
    /// \param gen [in/out] Pseudorandom number generator for the
    ///   normal(0,1) distribution.
    /// \param numRows [in] Number of rows in the test matrix.
    /// \param numCols [in] Number of columns in the test matrix.
    /// \param numPartitions [in] Number of parallel partitions (must
    ///   be a positive integer).
    /// \param cacheSizeHint [in] Cache size hint, in bytes.  Zero
    ///   means pick a reasonable default.
    /// \param contiguousCacheBlocks [in] Whether cache blocks in the
    ///   matrix to factor should be stored contiguously.
    /// \param printFieldNames [in] If humanReadable is true, this is
    ///   ignored; otherwise, whether to print a line of field names
    ///   before the line of output.
    /// \param humanReadable [in] Whether to print output that is easy
    ///   for humans to read, or instead to print output that is easy
    ///   for a script to parse.
    /// \param debug [in] Whether to print extra debugging output to
    ///   stderr.
    ///
    template<class Ordinal, class Scalar, class NodeType>
    void
    verifyKokkosNodeTsqr (const Teuchos::RCP<const NodeType>& node,
			  TSQR::Random::NormalGenerator<Ordinal, Scalar>& gen,
			  const Ordinal numRows,
			  const Ordinal numCols,
			  const int numPartitions,
			  const size_t cacheSizeHint,
			  const bool contiguousCacheBlocks,
			  const bool printFieldNames,
			  const bool humanReadable,
			  const bool debug)
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;
      using Teuchos::TypeNameTraits;
      using std::cerr;
      using std::cout;
      using std::endl;
      typedef TSQR::KokkosNodeTsqr<Ordinal, Scalar, NodeType> node_tsqr_type;
      typedef typename node_tsqr_type::FactorOutput factor_output_type;
      typedef Teuchos::ScalarTraits<Scalar> STS;
      typedef typename STS::magnitudeType magnitude_type;
      typedef Teuchos::Time timer_type;
      typedef Matrix<Ordinal, Scalar> matrix_type;

      const std::string scalarTypeName = TypeNameTraits<Scalar>::name();

      // Set up TSQR implementation.
      RCP<ParameterList> params = parameterList ("Intranode TSQR");
      params->set ("Cache Size Hint", cacheSizeHint);
      params->set ("Num Partitions", numPartitions);
      node_tsqr_type actor (node, params);
      if (debug)
	{
	  cerr << actor.description() << endl;
	  if (contiguousCacheBlocks)
	    cerr << "-- Test with contiguous cache blocks" << endl;
	}

      // Allocate space for test problem.
      matrix_type A (numRows, numCols);
      matrix_type A_copy (numRows, numCols);
      matrix_type Q (numRows, numCols);
      matrix_type R (numCols, numCols);
      if (std::numeric_limits<Scalar>::has_quiet_NaN)
	{
	  A.fill (std::numeric_limits<Scalar>::quiet_NaN());
	  A_copy.fill (std::numeric_limits<Scalar>::quiet_NaN());
	  Q.fill (std::numeric_limits<Scalar>::quiet_NaN());
	  R.fill (std::numeric_limits<Scalar>::quiet_NaN());
	}
      const Ordinal lda = numRows;
      const Ordinal ldq = numRows;
      const Ordinal ldr = numCols;

      // Create a test problem
      nodeTestProblem (gen, numRows, numCols, A.get(), A.lda(), true);

      if (debug)
	cerr << "-- Generated test problem" << endl;

      // Copy A into A_copy, since TSQR overwrites the input.  If
      // specified, rearrange the data in A_copy so that the data in
      // each cache block is contiguously stored.	  
      if (! contiguousCacheBlocks)
	{
	  A_copy.copy (A);
	  if (debug)
	    cerr << "-- Copied test problem from A into A_copy" << endl;
	}
      else
	{
	  actor.cache_block (numRows, numCols, A_copy.get(), A.get(), A.lda());
	  if (debug)
	    cerr << "-- Reorganized test matrix to have contiguous "
	      "cache blocks" << endl;

	  // Verify cache blocking, when in debug mode.
	  if (debug)
	    {
	      matrix_type A2 (numRows, numCols);
	      if (std::numeric_limits< Scalar >::has_quiet_NaN)
		A2.fill (std::numeric_limits<Scalar>::quiet_NaN());

	      actor.un_cache_block (numRows, numCols, A2.get(), A2.lda(), A_copy.get());
	      if (A == A2)
		{
		  if (debug)
		    cerr << "-- Cache blocking test succeeded!" << endl;
		}
	      else
		throw std::logic_error ("Cache blocking failed");
	    }
	}

      // Fill R with zeros, since the factorization may not overwrite
      // the strict lower triangle of R.
      R.fill (STS::zero());

      // Factor the matrix and compute the explicit Q factor
      factor_output_type factor_output = 
	actor.factor (numRows, numCols, A_copy.get(), A_copy.lda(), R.get(), 
		      R.lda(), contiguousCacheBlocks);
      if (debug)
	cerr << "-- Finished factor()" << endl;
      actor.explicit_Q (numRows, numCols, A_copy.get(), A_copy.lda(), 
			factor_output, numCols, Q.get(), Q.lda(), 
			contiguousCacheBlocks);
      if (debug)
	cerr << "-- Finished explicit_Q()" << endl;

      // "Un"-cache-block the output Q (the explicit Q factor), if
      // contiguous cache blocks were used.  This is only necessary
      // because local_verify() doesn't currently support contiguous
      // cache blocks.
      if (contiguousCacheBlocks)
	{
	  // Use A_copy as temporary storage for un-cache-blocking Q.
	  actor.un_cache_block (numRows, numCols, A_copy.get(), 
				A_copy.lda(), Q.get());
	  Q.copy (A_copy);
	  if (debug)
	    cerr << "-- Un-cache-blocked output Q factor" << endl;
	}

      // Print out the R factor
      if (debug)
	{
	  cerr << endl << "-- R factor:" << endl;
	  print_local_matrix (cerr, numCols, numCols, R.get(), R.lda());
	  cerr << endl;
	}

      // Validate the factorization
      std::vector<magnitude_type> results =
	local_verify (numRows, numCols, A.get(), lda, 
		      Q.get(), ldq, R.get(), ldr);
      if (debug)
	cerr << "-- Finished local_verify" << endl;

      // Print the results
      if (humanReadable)
	cout << "KokkosNodeTsqr:" << endl
	     << "Scalar type: " << scalarTypeName << endl
	     << "# rows = " << numRows << endl
	     << "# columns = " << numCols << endl
	     << "# partitions: " << numPartitions << endl
	     << "cache size hint (revised) = " << actor.cache_block_size() << endl
	     << "contiguous cache blocks? " << contiguousCacheBlocks << endl
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
		   << ",numPartitions"
		   << ",cacheSizeHint"
		   << ",contiguousCacheBlocks"
		   << ",absFrobResid"
		   << ",absFrobOrthog"
		   << ",frobA"
		   << endl;
	    }
	  cout << "KokkosNodeTsqr"
	       << "," << scalarTypeName
	       << "," << numRows
	       << "," << numCols
	       << "," << numPartitions
	       << "," << actor.cache_block_size()
	       << "," << contiguousCacheBlocks 
	       << "," << results[0]
	       << "," << results[1]
	       << "," << results[2]
	       << endl;
	}
    }


    /// \fn benchmarkKokkosNodeTsqr
    /// \brief Test performance of KokkosNodeTsqr's QR factorization.
    ///
    /// Compare the performance of KokkosNodeTsqr's QR factorization
    /// to that of LAPACK's QR factorization.  Print results to
    /// stdout.
    ///
    /// \param node [in] The Kokkos Node instance on which to execute
    ///   in parallel.
    /// \param numTrials [in] Number of times to run the benchmark;
    ///   the timing result is cumulative over all trials.  Timing
    ///   over larger numbers of trials improves certainty of the 
    ///   result.
    /// \param numRows [in] Number of rows in the test matrix.
    /// \param numCols [in] Number of columns in the test matrix.
    /// \param numPartitions [in] Number of parallel partitions (must
    ///   be a positive integer).
    /// \param cacheSizeHint [in] Cache size hint, in bytes.  Zero
    ///   means pick a reasonable default.
    /// \param contiguousCacheBlocks [in] Whether cache blocks in the
    ///   matrix to factor should be stored contiguously.
    /// \param printFieldNames [in] If humanReadable is true, this is
    ///   ignored; otherwise, whether to print a line of field names
    ///   before the line of output.
    /// \param humanReadable [in] Whether to print output that is easy
    ///   for humans to read, or instead to print output that is easy
    ///   for a script to parse.
    ///
    template<class Ordinal, class Scalar, class NodeType>
    void
    benchmarkKokkosNodeTsqr (const Teuchos::RCP<const NodeType>& node,
			     const int numTrials,
			     const Ordinal numRows, 
			     const Ordinal numCols, 
			     const int numPartitions,
			     const size_t cacheSizeHint,
			     const bool contiguousCacheBlocks,
			     const bool printFieldNames,
			     const bool humanReadable)
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;
      using Teuchos::TypeNameTraits;
      using std::cerr;
      using std::cout;
      using std::endl;
      typedef TSQR::KokkosNodeTsqr<Ordinal, Scalar, NodeType> node_tsqr_type;
      typedef typename node_tsqr_type::FactorOutput factor_output_type;
      typedef Teuchos::ScalarTraits<Scalar> STS;
      typedef typename STS::magnitudeType magnitude_type;
      typedef Teuchos::Time timer_type;
      typedef Matrix<Ordinal, Scalar> matrix_type;

      const std::string scalarTypeName = TypeNameTraits<Scalar>::name();

      // Pseudorandom normal(0,1) generator.  Default seed is OK,
      // because this is a benchmark, not an accuracy test.
      TSQR::Random::NormalGenerator<Ordinal, Scalar> gen;

      // Set up TSQR implementation.
      RCP<ParameterList> params = parameterList ("Intranode TSQR");
      params->set ("Cache Size Hint", cacheSizeHint);
      params->set ("Num Partitions", numPartitions);
      node_tsqr_type actor (node, params);

      // Allocate space for test problem.
      matrix_type A (numRows, numCols);
      matrix_type A_copy (numRows, numCols);
      matrix_type Q (numRows, numCols);
      matrix_type R (numCols, numCols);

      // Fill R with zeros, since the factorization may not overwrite
      // the strict lower triangle of R.
      R.fill (STS::zero());

      // Create a test problem
      nodeTestProblem (gen, numRows, numCols, A.get(), A.lda(), false);

      // Copy A into A_copy, since TSQR overwrites the input.  If
      // specified, rearrange the data in A_copy so that the data in
      // each cache block is contiguously stored.	  
      if (contiguousCacheBlocks)
	actor.cache_block (numRows, numCols, A_copy.get(), A.get(), A.lda());
      else
	A_copy.copy (A);

      // Do a few timing runs and throw away the results, just to warm
      // up any libraries that do autotuning.
      const int numWarmupRuns = 5;
      for (int warmupRun = 0; warmupRun < numWarmupRuns; ++warmupRun)
	{
	  // Factor the matrix in-place in A_copy, and extract the
	  // resulting R factor into R.
	  factor_output_type factor_output = 
	    actor.factor (numRows, numCols, A_copy.get(), A_copy.lda(), 
			  R.get(), R.lda(), contiguousCacheBlocks);
	  // Compute the explicit Q factor (which was stored
	  // implicitly in A_copy and factor_output) and store in Q.
	  // We don't need to un-cache-block the output, because we
	  // aren't verifying it here.
	  actor.explicit_Q (numRows, numCols, A_copy.get(), A_copy.lda(), 
			    factor_output, numCols, Q.get(), Q.lda(), 
			    contiguousCacheBlocks);
	}

      // Benchmark intranode TSQR for numTrials trials.
      //
      // Name of timer doesn't matter here; we only need the timing.
      timer_type timer("KokkosNodeTsqr");
      timer.start();
      for (int trialNum = 0; trialNum < numTrials; ++trialNum)
	{
	  // Factor the matrix in-place in A_copy, and extract the
	  // resulting R factor into R.
	  factor_output_type factor_output = 
	    actor.factor (numRows, numCols, A_copy.get(), A_copy.lda(), 
			  R.get(), R.lda(), contiguousCacheBlocks);
	  // Compute the explicit Q factor (which was stored
	  // implicitly in A_copy and factor_output) and store in Q.
	  // We don't need to un-cache-block the output, because we
	  // aren't verifying it here.
	  actor.explicit_Q (numRows, numCols, A_copy.get(), A_copy.lda(), 
			    factor_output, numCols, Q.get(), Q.lda(), 
			    contiguousCacheBlocks);
	}
      const double timing = timer.stop();

      // Print the results
      if (humanReadable)
	{
	  cout << "KokkosNodeTsqr cumulative timings:" << endl
	       << "Scalar type: " << scalarTypeName << endl
	       << "# rows = " << numRows << endl
	       << "# columns = " << numCols << endl
	       << "# partitions: " << numPartitions << endl
	       << "Cache size hint (in bytes) = " << actor.cache_block_size() << endl
	       << "Contiguous cache blocks? " << contiguousCacheBlocks << endl
	       << "# trials = " << numTrials << endl
	       << "Total time (s) = " << timing << endl;
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
		   << ",numPartitions"
		   << ",cacheSizeHint"
		   << ",contiguousCacheBlocks"
		   << ",numTrials"
		   << ",timing"
		   << endl;
	    }

	  // We don't include {min,max}_seq_apply_timing() here, because
	  // those times don't benefit from the accuracy of benchmarking
	  // for numTrials > 1.  Thus, it's misleading to include them
	  // with tbb_tsqr_timing, the total time over numTrials trials.
	  cout << "KokkosNodeTsqr"
	       << "," << scalarTypeName
	       << "," << numRows
	       << "," << numCols
	       << "," << numPartitions
	       << "," << actor.cache_block_size()
	       << "," << contiguousCacheBlocks 
	       << "," << numTrials
	       << "," << timing
	       << endl;
	}
    }
  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_KokkosNodeTsqrTest_hpp
