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

#ifndef __TSQR_Test_TsqrTest_hpp
#define __TSQR_Test_TsqrTest_hpp

#include <Tsqr_Config.hpp>

#include <Tsqr.hpp>
#ifdef HAVE_TSQR_INTEL_TBB
#  include <TbbTsqr.hpp>
#endif // HAVE_TSQR_INTEL_TBB
#include <Tsqr_TestSetup.hpp>
#include <Tsqr_GlobalVerify.hpp>
#include <Tsqr_printGlobalMatrix.hpp>
#include <Tsqr_verifyTimerConcept.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <cstring> // size_t
#include <iostream>
#include <stdexcept>
#include <string>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Test {

    template< class TsqrType >
    class TsqrVerifier {
    public:
      typedef TsqrType tsqr_type;
      typedef typename tsqr_type::scalar_type scalar_type;
      typedef typename tsqr_type::ordinal_type ordinal_type;
      typedef Matrix< ordinal_type, scalar_type > matrix_type;
      typedef typename tsqr_type::FactorOutput factor_output_type;
      typedef MessengerBase< scalar_type > messenger_type;
      typedef Teuchos::RCP< messenger_type > messenger_ptr;

      static void
      verify (tsqr_type& tsqr,
	      const messenger_ptr& scalarComm,
	      const matrix_type& A_local,
	      matrix_type& A_copy,
	      matrix_type& Q_local,
	      matrix_type& R,
	      const bool contiguousCacheBlocks,
	      const bool b_debug = false)
      {
	using std::cerr;
	using std::endl;

	const ordinal_type nrows_local = A_local.nrows();
	const ordinal_type ncols = A_local.ncols();

	// If specified, rearrange cache blocks in the copy.
	if (contiguousCacheBlocks)
	  {
	    tsqr.cache_block (nrows_local, ncols, A_copy.get(), 
			      A_local.get(), A_local.lda());
	    if (b_debug)
	      {
		scalarComm->barrier();
		if (scalarComm->rank() == 0)
		  cerr << "-- Cache-blocked input matrix to factor." << endl;
	      }
	  }
	else
	  A_copy.copy (A_local);

	const bool testFactorExplicit = true;
	if (testFactorExplicit)
	  {
	    tsqr.factorExplicit (A_copy.view(), Q_local.view(), R.view(), 
				 contiguousCacheBlocks);
	    if (b_debug)
	      {
		scalarComm->barrier();
		if (scalarComm->rank() == 0)
		  cerr << "-- Finished Tsqr::factorExplicit" << endl;
	      }
	  }
	else
	  {
	    // Factor the (copy of the) matrix.
	    factor_output_type factorOutput = 
	      tsqr.factor (nrows_local, ncols, A_copy.get(), A_copy.lda(), 
			   R.get(), R.lda(), contiguousCacheBlocks);
	    if (b_debug)
	      {
		scalarComm->barrier();
		if (scalarComm->rank() == 0)
		  cerr << "-- Finished Tsqr::factor" << endl;
	      }

	    // Compute the explicit Q factor in Q_local
	    tsqr.explicit_Q (nrows_local, 
			     ncols, A_copy.get(), A_copy.lda(), factorOutput, 
			     ncols, Q_local.get(), Q_local.lda(), 
			     contiguousCacheBlocks);
	    if (b_debug)
	      {
		scalarComm->barrier();
		if (scalarComm->rank() == 0)
		  cerr << "-- Finished Tsqr::explicit_Q" << endl;
	      }
	  }

	// "Un"-cache-block the output, if contiguous cache blocks were
	// used.  This is only necessary because global_verify() doesn't
	// currently support contiguous cache blocks.
	if (contiguousCacheBlocks)
	  {
	    // We can use A_copy as scratch space for un-cache-blocking
	    // Q_local, since we're done using A_copy for other things.
	    tsqr.un_cache_block (nrows_local, ncols, A_copy.get(), 
				 A_copy.lda(), Q_local.get());
	    // Overwrite Q_local with the un-cache-blocked Q factor.
	    Q_local.copy (A_copy);

	    if (b_debug)
	      {
		scalarComm->barrier();
		if (scalarComm->rank() == 0)
		  cerr << "-- Un-cache-blocked output Q factor" << endl;
	      }
	  }
      }
    };

    /// \function verifyTsqr
    /// \brief Test and print to stdout the accuracy of parallel TSQR
    ///
    /// \param which [in] Valid values: "MpiTbbTSQR" (for TBB-parallel
    ///   node-level TSQR underneath MPI-parallel TSQR), "MpiSeqTSQR"
    ///   (for cache-blocked sequential node-level TSQR underneath
    ///   MPI-parallel TSQR)
    ///
    /// \param scalarTypeName [in] Name of the Scalar type
    ///
    /// \param generator [in/out] Normal(0,1) (pseudo)random number
    ///   generator.  Only touched on MPI process 0.  Used to generate
    ///   random test matrices for the factorization.
    ///
    /// \param nrows_global [in] Number of rows in the entire test
    ///   matrix (over all processes) to generate.  The matrix will be
    ///   divided up in blocks of contiguous rows among the processes.
    ///
    /// \param ncols [in] Number of columns in the test matrix to
    ///   generate.
    ///
    /// \param ordinalComm [in/out] Object for communicating Ordinal
    ///   (integer index) objects among the processes
    ///
    /// \param scalarComm [in/out] Object for communicating Scalar
    ///   (matrix data) objects among the processes
    ///
    /// \param num_cores [in] Number of cores to use per MPI process
    ///   for Intel TBB parallelism within that process
    ///
    /// \param cache_size_hint [in] Cache size hint (per core) in
    ///   bytes.  If zero, a sensible default is used.
    ///
    /// \param contiguousCacheBlocks [in] Whether cache blocks
    ///   should be stored contiguously
    ///
    /// \param printFieldNames [in] Whether to print field names (only
    ///   appliable if not human_readable)
    ///
    /// \param human_readable [in] Whether output should be human
    ///   readable, or machine parseable
    ///
    /// \param b_debug [in] Whether to print debug output
    ///
    template< class Ordinal, class Scalar, class Generator >
    void
    verifyTsqr (const std::string& which,
		const std::string& scalarTypeName,
		Generator& generator,
		const Ordinal nrows_global,
		const Ordinal ncols,
		const Teuchos::RCP< MessengerBase< Ordinal > >& ordinalComm,
		const Teuchos::RCP< MessengerBase< Scalar > >& scalarComm,
		const int num_cores = 1,
		const size_t cache_size_hint = 0,
		const bool contiguousCacheBlocks,
		const bool printFieldNames,
		const bool human_readable = false,
		const bool b_debug = false)
    {
      typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
      using std::cerr;
      using std::cout;
      using std::endl;

      const bool b_extra_debug = false;
      const int nprocs = scalarComm->size();
      const int my_rank = scalarComm->rank();
      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "tsqr_verify:" << endl;
	  scalarComm->barrier();
	}
      const Ordinal nrows_local = numLocalRows (nrows_global, my_rank, nprocs);

      // Set up storage for the test problem.
      Matrix< Ordinal, Scalar > A_local (nrows_local, ncols);
      Matrix< Ordinal, Scalar > Q_local (nrows_local, ncols);
      if (std::numeric_limits< Scalar >::has_quiet_NaN)
	{
	  A_local.fill (std::numeric_limits< Scalar >::quiet_NaN());
	  Q_local.fill (std::numeric_limits< Scalar >::quiet_NaN());
	}
      Matrix< Ordinal, Scalar > R (ncols, ncols, Scalar(0));

      // Generate the test problem.
      distributedTestProblem (generator, A_local, ordinalComm.get(), scalarComm.get());
      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Generated test problem." << endl;
	}

      // Make sure that the test problem (the matrix to factor) was
      // distributed correctly.
      if (b_extra_debug && b_debug)
	{
	  if (my_rank == 0)
	    cerr << "Test matrix A:" << endl;
	  scalarComm->barrier ();
	  printGlobalMatrix (cerr, A_local, scalarComm.get(), ordinalComm.get());
	  scalarComm->barrier ();
	}

      // Factoring the matrix stored in A_local overwrites it, so we
      // make a copy of A_local.  Initialize with NaNs to make sure
      // that cache blocking works correctly (if applicable).
      Matrix< Ordinal, Scalar > A_copy (nrows_local, ncols);
      if (std::numeric_limits< Scalar >::has_quiet_NaN)
	A_copy.fill (std::numeric_limits< Scalar >::quiet_NaN());

      // actual_cache_size_hint: "cache_size_hint" is just a
      // suggestion.  TSQR determines the cache size hint itself;
      // this remembers it so we can print it out later.
      size_t actual_cache_size_hint;

      if (which == "MpiTbbTSQR")
	{
#ifdef HAVE_TSQR_INTEL_TBB
	  using Teuchos::RCP;
	  typedef TSQR::TBB::TbbTsqr< Ordinal, Scalar > node_tsqr_type;
	  typedef TSQR::DistTsqr< Ordinal, Scalar > dist_tsqr_type;
	  typedef Tsqr< Ordinal, Scalar, node_tsqr_type, dist_tsqr_type > tsqr_type;

	  RCP< node_tsqr_type > node_tsqr (new node_tsqr_type (num_cores, cache_size_hint));
	  RCP< dist_tsqr_type > dist_tsqr (new dist_tsqr_type (scalarComm));
	  tsqr_type tsqr (node_tsqr, dist_tsqr);
	  
	  // Compute the factorization and explicit Q factor.
	  TsqrVerifier< tsqr_type >::verify (tsqr, scalarComm, A_local, A_copy, 
					     Q_local, R, contiguousCacheBlocks, 
					     b_debug);
	  // Save the "actual" cache block size
	  actual_cache_size_hint = tsqr.cache_size_hint();
#else
	  throw std::logic_error("TSQR not built with Intel TBB support");
#endif // HAVE_TSQR_INTEL_TBB
	}
      else if (which == "MpiSeqTSQR")
	{
	  using Teuchos::RCP;
	  typedef SequentialTsqr< Ordinal, Scalar > node_tsqr_type;
	  typedef TSQR::DistTsqr< Ordinal, Scalar > dist_tsqr_type;
	  typedef Tsqr< Ordinal, Scalar, node_tsqr_type, dist_tsqr_type > tsqr_type;

	  RCP< node_tsqr_type > node_tsqr (new node_tsqr_type (cache_size_hint));
	  RCP< dist_tsqr_type > dist_tsqr (new dist_tsqr_type (scalarComm));
	  tsqr_type tsqr (node_tsqr, dist_tsqr);
	  
	  // Compute the factorization and explicit Q factor.
	  TsqrVerifier< tsqr_type >::verify (tsqr, scalarComm, A_local, A_copy, 
					     Q_local, R, contiguousCacheBlocks, 
					     b_debug);
	  // Save the "actual" cache block size
	  actual_cache_size_hint = tsqr.cache_size_hint();
	}
      else 
	throw std::logic_error("Unknown TSQR implementation type \"" + which + "\"");

      // Print out the Q and R factors
      if (b_extra_debug && b_debug)
	{
	  if (my_rank == 0)
	    cerr << endl << "Q factor:" << endl;
	  scalarComm->barrier();
	  printGlobalMatrix (cerr, Q_local, scalarComm.get(), ordinalComm.get());
	  scalarComm->barrier ();
	  if (my_rank == 0)
	    {
	      cerr << endl << "R factor:" << endl;
	      print_local_matrix (cerr, ncols, ncols, R.get(), R.lda());
	      cerr << endl;
	    }
	  scalarComm->barrier ();
	}

      // Test accuracy of the resulting factorization
      std::vector< magnitude_type > results = 
	global_verify (nrows_local, ncols, A_local.get(), A_local.lda(),
		       Q_local.get(), Q_local.lda(), R.get(), R.lda(), 
		       scalarComm.get());
      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Finished global_verify" << endl;
	}

      // Print the results on Proc 0.
      if (my_rank == 0)
	{
	  if (human_readable)
	    {
	      std::string human_readable_name;

	      if (which == "MpiSeqTSQR")
		human_readable_name = "MPI parallel / cache-blocked TSQR";
	      else if (which == "MpiTbbTSQR")
		{
#ifdef HAVE_TSQR_INTEL_TBB
		  human_readable_name = "MPI parallel / TBB parallel / cache-blocked TSQR";
#else
		  throw std::logic_error("TSQR not built with Intel TBB support");
#endif // HAVE_TSQR_INTEL_TBB
		}
	      else 
		throw std::logic_error("Unknown TSQR implementation type \"" + which + "\"");

	      cout << human_readable_name << ":" << endl
		   << "Scalar type: " << scalarTypeName << endl
		   << "# rows: " << nrows_global << endl
		   << "# columns: " << ncols << endl
		   << "# MPI processes: " << nprocs << endl;
#ifdef HAVE_TSQR_INTEL_TBB
	      if (which == "MpiTbbTSQR")
		cout << "# cores per process = " << num_cores << endl;
#endif // HAVE_TSQR_INTEL_TBB
	      cout << "Cache size hint in bytes: " << actual_cache_size_hint << endl
		   << "Contiguous cache blocks? " << contiguousCacheBlocks << endl
		   << "Absolute residual $\\| A - Q R \\|_2: "
		   << results[0] << endl
		   << "Absolute orthogonality $\\| I - Q^* Q \\|_2$: " 
		   << results[1] << endl
		   << "Test matrix norm $\\| A \\|_F$: "
		   << results[2] << endl
		   << endl;
	    }
	  else
	    {
	      if (printFieldNames)
		{
		  cout << "%"
		       << "method"
		       << ",scalarType"
		       << ",globalNumRows"
		       << ",numCols"
		       << ",numProcs"
		       << ",numCores"
		       << ",cacheSizeHint"
		       << ",contiguousCacheBlocks"
		       << ",absFrobResid"
		       << ",absFrobOrthog"
		       << ",frobA" << endl;
		}

	      cout << which
		   << "," << scalarTypeName
		   << "," << nrows_global
		   << "," << ncols
		   << "," << nprocs;
#ifdef HAVE_TSQR_INTEL_TBB
	      if (which == "MpiTbbTSQR")
		cout << "," << num_cores;
	      else
		cout << ",1";
#else
	      cout << ",1" << endl;
#endif // HAVE_TSQR_INTEL_TBB
	      cout << "," << actual_cache_size_hint
		   << "," << contiguousCacheBlocks 
		   << "," << results[0] 
		   << "," << results[1]
		   << "," << results[2]
		   << endl;
	    }
	}
    }


    template< class TsqrBase, class TimerType >
    double
    do_tsqr_benchmark (const std::string& which,
		       TsqrBase& tsqr, 
		       const Teuchos::RCP< MessengerBase< typename TsqrBase::scalar_type > >& messenger,
		       const Matrix< typename TsqrBase::ordinal_type, typename TsqrBase::scalar_type >& A_local,
		       Matrix< typename TsqrBase::ordinal_type, typename TsqrBase::scalar_type >& A_copy,
		       Matrix< typename TsqrBase::ordinal_type, typename TsqrBase::scalar_type >& Q_local,
		       Matrix< typename TsqrBase::ordinal_type, typename TsqrBase::scalar_type >& R,
		       const int ntrials,
		       const bool contiguousCacheBlocks,
		       const bool human_readable,
		       const bool b_debug = false)
    {
      typedef typename TsqrBase::FactorOutput factor_output_type;
      typedef typename TsqrBase::ordinal_type ordinal_type;
      using std::cerr;
      using std::cout;
      using std::endl;

      const ordinal_type nrows_local = A_local.nrows();
      const ordinal_type ncols = A_local.ncols();

      if (contiguousCacheBlocks)
	{
	  tsqr.cache_block (nrows_local, ncols, A_copy.get(), 
			    A_local.get(), A_local.lda());
	  if (b_debug)
	    {
	      messenger->barrier();
	      if (messenger->rank() == 0)
		cerr << "-- Cache-blocked input matrix to factor." << endl;
	    }
	}
      else
	A_copy.copy (A_local);

      if (b_debug)
	{
	  messenger->barrier();
	  if (messenger->rank() == 0)
	    cerr << "-- Starting timing loop" << endl;
	}

      // Benchmark TSQR for ntrials trials.  The answer (the numerical
      // results of the factorization) is only valid if ntrials == 1,
      // but this is a benchmark and not a verification routine.  Call
      // tsqr_verify() if you want to determine whether TSQR computes
      // the right answer.
      //
      // Name of timer doesn't matter here; we only need the timing.
      TSQR::Test::verifyTimerConcept< TimerType >();
      TimerType timer (which);


      const bool testFactorExplicit = true;
      double tsqr_timing;
      if (testFactorExplicit)
	{
	  timer.start();
	  for (int trial_num = 0; trial_num < ntrials; ++trial_num)
	    tsqr.factorExplicit (A_copy.view(), Q_local.view(), R.view(), 
				 contiguousCacheBlocks);
	  tsqr_timing = timer.stop();
	}
      else
	{
	  timer.start();
	  for (int trial_num = 0; trial_num < ntrials; ++trial_num)
	    {
	      // Factor the matrix and compute the explicit Q factor.
	      // Don't worry about the fact that we're overwriting the
	      // input; this is a benchmark, not a numerical verification
	      // test.  (We have the latter implemented as tsqr_verify()
	      // in this file.)  For the same reason, don't worry about
	      // un-cache-blocking the output (when cache blocks are
	      // stored contiguously).
	      factor_output_type factor_output = 
		tsqr.factor (nrows_local, ncols, A_copy.get(), A_copy.lda(), 
			     R.get(), R.lda(), contiguousCacheBlocks);
	      tsqr.explicit_Q (nrows_local, 
			       ncols, A_copy.get(), A_copy.lda(), factor_output, 
			       ncols, Q_local.get(), Q_local.lda(), 
			       contiguousCacheBlocks);
	      // Timings in debug mode likely won't make sense, because
	      // Proc 0 is outputting the debug messages to cerr.
	      // Nevertheless, we don't put any "if(b_debug)" calls in the
	      // timing loop.
	    }
	  // Compute the resulting total time (in seconds) to execute
	  // ntrials runs of Tsqr::factor() and Tsqr::explicit_Q().  The
	  // time may differ on different MPI processes.
	  tsqr_timing = timer.stop();
	}

      if (b_debug)
	{
	  messenger->barrier();
	  if (messenger->rank() == 0)
	    cerr << "-- Finished timing loop" << endl;
	}
      return tsqr_timing;
    }

    /// \function benchmarkTsqr
    /// \brief Benchmark parallel TSQR and report timings to stdout
    ///
    /// Benchmark the MPI-parallel TSQR implementation specified by
    /// the "which" parameter (either with cache-blocked TSQR or
    /// TBB-parallel cache-blocked TSQR as the node-level
    /// implementation), for "ntrials" trials.  Print the stdout the
    /// cumulative run time (in seconds) for all ntrials trials.
    ///
    /// \param which [in] Valid values: "MpiTbbTSQR" (for TBB-parallel
    ///   node-level TSQR underneath MPI-parallel TSQR), "MpiSeqTSQR"
    ///   (for cache-blocked sequential node-level TSQR underneath
    ///   MPI-parallel TSQR)
    ///
    /// \param scalarTypeName [in] Name of the Scalar type
    ///
    /// \param generator [in/out] Normal(0,1) (pseudo)random number
    ///   generator.  Only touched on MPI process 0.  Used to generate
    ///   random test matrices for the factorization.
    ///
    /// \param ntrials [in] Number of trials to use in the benchmark.
    ///   Reported timings are cumulative over all trials.
    ///
    /// \param nrows_global [in] Number of rows in the entire test
    ///   matrix (over all processes) to generate.  The matrix will be
    ///   divided up in blocks of contiguous rows among the processes.
    ///
    /// \param ncols [in] Number of columns in the test matrix to
    ///   generate.
    ///
    /// \param ordinalComm [in/out] Object for communicating Ordinal
    ///   (integer index) objects among the processes
    ///
    /// \param scalarComm [in/out] Object for communicating Scalar
    ///   (matrix data) objects among the processes
    ///
    /// \param num_cores [in] Number of cores to use per MPI process
    ///   for Intel TBB parallelism within that process
    ///
    /// \param cache_size_hint [in] Cache block size (per core) in
    ///   bytes.  If zero, a sensible default is used.
    ///
    /// \param contiguousCacheBlocks [in] Whether cache blocks
    ///   should be stored contiguously
    ///
    /// \param printFieldNames [in] Whether to print field names (only
    ///   appliable if not human_readable)
    ///
    /// \param human_readable [in] Whether output should be human
    ///   readable, or machine parseable
    ///
    /// \param b_debug [in] Whether to print debug output
    ///
    template< class Ordinal, class Scalar, class Generator, class TimerType >
    void
    benchmarkTsqr (const std::string& which,
		   const std::string& scalarTypeName,
		   Generator& generator,
		   const int ntrials,
		   const Ordinal nrows_global,
		   const Ordinal ncols,
		   const Teuchos::RCP< MessengerBase< Ordinal > >& ordinalComm,
		   const Teuchos::RCP< MessengerBase< Scalar > >& scalarComm,
		   const Ordinal num_cores,
		   const size_t cache_size_hint,
		   const bool contiguousCacheBlocks,
		   const bool printFieldNames,
		   const bool human_readable,
		   const bool b_debug)
    {
      using std::cerr;
      using std::cout;
      using std::endl;

      TSQR::Test::verifyTimerConcept< TimerType >();
      const bool b_extra_debug = false;
      const int nprocs = scalarComm->size();
      const int my_rank = scalarComm->rank();
      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "tsqr_benchmark:" << endl;
	  scalarComm->barrier();
	}
      const Ordinal nrows_local = numLocalRows (nrows_global, my_rank, nprocs);

      // Set up storage for the test problem.
      Matrix< Ordinal, Scalar > A_local (nrows_local, ncols);
      Matrix< Ordinal, Scalar > Q_local (nrows_local, ncols);
      if (std::numeric_limits< Scalar >::has_quiet_NaN)
	{
	  A_local.fill (std::numeric_limits< Scalar >::quiet_NaN());
	  Q_local.fill (std::numeric_limits< Scalar >::quiet_NaN());
	}
      Matrix< Ordinal, Scalar > R (ncols, ncols, Scalar(0));

      // Generate the test problem.
      distributedTestProblem (generator, A_local, ordinalComm.get(), scalarComm.get());
      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Generated test problem." << endl;
	}

      // Make sure that the test problem (the matrix to factor) was
      // distributed correctly.
      if (b_extra_debug && b_debug)
	{
	  if (my_rank == 0)
	    cerr << "Test matrix A:" << endl;
	  scalarComm->barrier ();
	  printGlobalMatrix (cerr, A_local, scalarComm.get(), ordinalComm.get());
	  scalarComm->barrier ();
	}

      // Factoring the matrix stored in A_local overwrites it, so we
      // make a copy of A_local.  If specified, rearrange cache blocks
      // in the copy.  Initialize with NaNs to make sure that cache
      // blocking worked correctly.
      Matrix< Ordinal, Scalar > A_copy (nrows_local, ncols);
      if (std::numeric_limits< Scalar >::has_quiet_NaN)
	A_copy.fill (std::numeric_limits< Scalar >::quiet_NaN());

      // actual_cache_size_hint: "cache_size_hint" is just a
      // suggestion.  TSQR determines the cache block size itself;
      // this remembers it so we can print it out later.
      size_t actual_cache_size_hint;
      // Run time (in seconds, as a double-precision floating-point
      // value) for TSQR on this MPI node.
      double tsqr_timing;

      if (which == "MpiTbbTSQR")
	{
#ifdef HAVE_TSQR_INTEL_TBB
	  using Teuchos::RCP;
	  typedef TSQR::TBB::TbbTsqr< Ordinal, Scalar > node_tsqr_type;
	  typedef TSQR::DistTsqr< Ordinal, Scalar > dist_tsqr_type;
	  typedef Tsqr< Ordinal, Scalar, node_tsqr_type, dist_tsqr_type > tsqr_type;

	  RCP< node_tsqr_type > nodeTsqr (new node_tsqr_type (num_cores, cache_size_hint));
	  RCP< dist_tsqr_type > distTsqr (new dist_tsqr_type (scalarComm));
	  tsqr_type tsqr (nodeTsqr, distTsqr);

	  // Run the benchmark.
	  tsqr_timing = 
	    do_tsqr_benchmark< tsqr_type, TimerType > (which, tsqr, scalarComm, A_local,
						       A_copy, Q_local, R, ntrials, 
						       contiguousCacheBlocks, 
						       human_readable, b_debug);

	  // Save the "actual" cache block size
	  actual_cache_size_hint = tsqr.cache_size_hint();
#else
	  throw std::logic_error("TSQR not built with Intel TBB support");
#endif // HAVE_TSQR_INTEL_TBB
	}
      else if (which == "MpiSeqTSQR")
	{
	  using Teuchos::RCP;
	  typedef SequentialTsqr< Ordinal, Scalar > node_tsqr_type;
	  typedef TSQR::DistTsqr< Ordinal, Scalar > dist_tsqr_type;
	  typedef Tsqr< Ordinal, Scalar, node_tsqr_type, dist_tsqr_type > tsqr_type;

	  // Set up TSQR.
	  RCP< node_tsqr_type > nodeTsqr (new node_tsqr_type (cache_size_hint));
	  RCP< dist_tsqr_type > distTsqr (new dist_tsqr_type (scalarComm));
	  tsqr_type tsqr (nodeTsqr, distTsqr);
	  
	  // Run the benchmark.
	  tsqr_timing = 
	    do_tsqr_benchmark< tsqr_type, TimerType > (which, tsqr, scalarComm, A_local,
						       A_copy, Q_local, R, ntrials, 
						       contiguousCacheBlocks, 
						       human_readable, b_debug);
	  // Save the "actual" cache block size
	  actual_cache_size_hint = tsqr.cache_size_hint();
	}
      else
	throw std::logic_error("Unknown TSQR implementation type \"" + which + "\"");

      // Find the min and max TSQR timing on all processors.
      const double min_tsqr_timing = scalarComm->globalMin (tsqr_timing);
      const double max_tsqr_timing = scalarComm->globalMax (tsqr_timing);

      // Print the results on Proc 0.
      if (my_rank == 0)
	{
	  if (human_readable)
	    {
	      std::string human_readable_name;

	      if (which == "MpiSeqTSQR")
		human_readable_name = "MPI parallel / cache-blocked TSQR";
	      else if (which == "MpiTbbTSQR")
		{
#ifdef HAVE_TSQR_INTEL_TBB
		  human_readable_name = "MPI parallel / TBB parallel / cache-blocked TSQR";
#else
		  throw std::logic_error("TSQR not built with Intel TBB support");
#endif // HAVE_TSQR_INTEL_TBB
		}
	      else 
		throw std::logic_error("Unknown TSQR implementation type \"" + which + "\"");

	      cout << human_readable_name << ":" << endl
		   << "Scalar type: " << scalarTypeName << endl
		   << "# rows: " << nrows_global << endl
		   << "# columns: " << ncols << endl
		   << "# MPI processes: " << nprocs << endl;

#ifdef HAVE_TSQR_INTEL_TBB
	      if (which == "MpiTbbTSQR")
		cout << "# cores per process: " << num_cores << endl;
#endif // HAVE_TSQR_INTEL_TBB

	      cout << "Cache size hint in bytes: " << actual_cache_size_hint << endl
		   << "contiguous cache blocks? " << contiguousCacheBlocks << endl
		   << "# trials: " << ntrials << endl
		   << "Min total time (s) over all MPI processes: " 
		   << min_tsqr_timing << endl
		   << "Max total time (s) over all MPI processes: " 
		   << max_tsqr_timing << endl
		   << endl;
	    }
	  else
	    {
	      if (printFieldNames)
		{
		  cout << "%"
		       << "method"
		       << ",scalarType"
		       << ",globalNumRows"
		       << ",numCols"
		       << ",numProcs"
		       << ",numCores"
		       << ",cacheSizeHint"
		       << ",contiguousCacheBlocks"
		       << ",numTrials"
		       << ",minTiming"
		       << ",maxTiming" 
		       << endl;
		}
	      cout << which
		   << "," << scalarTypeName
		   << "," << nrows_global
		   << "," << ncols 
		   << "," << nprocs;
#ifdef HAVE_TSQR_INTEL_TBB
	      if (which == "MpiTbbTSQR")
		cout << "," << num_cores;
	      else 
		cout << ",1";
#else
	      cout << ",1";
#endif // HAVE_TSQR_INTEL_TBB
	      cout << "," << actual_cache_size_hint
		   << "," << contiguousCacheBlocks
		   << "," << ntrials 
		   << "," << min_tsqr_timing 
		   << "," << max_tsqr_timing 
		   << endl;
	    }
	}
    }


  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_TsqrTest_hpp
