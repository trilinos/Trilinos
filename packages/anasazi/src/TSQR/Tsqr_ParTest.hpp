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

#ifndef __TSQR_Test_DistTest_hpp
#define __TSQR_Test_DistTest_hpp

#include <Tsqr_Config.hpp>
#include <Tsqr_Random_NormalGenerator.hpp>
#include <Tsqr_verifyTimerConcept.hpp>

#include <Tsqr_generateStack.hpp>
#include <Tsqr_DistTsqr.hpp>
#include <Tsqr_GlobalVerify.hpp>
#include <Tsqr_printGlobalMatrix.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Test {

    /// \class DistTsqrVerifier
    /// \brief Template version of DistTsqr accuracy test
    ///
    template< class Ordinal, class Scalar >
    class DistTsqrVerifier {
      TSQR::Random::NormalGenerator< Ordinal, Scalar > gen_;
      Teuchos::RCP< MessengerBase< Scalar > > const scalarComm_;
      std::string scalarTypeName_;
      const bool humanReadable_, debug_;
      std::ostream& out_;
      std::ostream& err_;
      bool testFactorExplicit_, testFactorImplicit_;

    public:
      typedef Ordinal ordinal_type;
      typedef Scalar scalar_type;
      typedef typename ScalarTraits< scalar_type >::magnitude_type magnitude_type;
      typedef typename std::pair< magnitude_type, magnitude_type > result_type;

      /// \brief Constructor, with custom seed value
      ///
      /// \param scalarComm [in/out] Communicator object over which to
      ///   test.
      /// \param seed [in] 4-element vector; the random seed input of
      ///   TSQR::Random::NormalGenerator (which see, since there are
      ///   restrictions on the set of valid seeds)
      /// \param scalarTypeName [in] Human-readable name of the Scalar
      ///   template type parameter
      /// \param out [out] Output stream to which to write results
      /// \param err [out] Output stream to which to write any
      ///   debugging outputs (if applicable) or errors
      /// \param testFactorExplicit [in] Whether to test 
      ///   DistTsqr::factorExplicit()
      /// \param testFactorImplicit [in] Whether to test 
      ///   DistTsqr::factor() and DistTsqr::explicit_Q()
      /// \param humanReadable [in] Whether printed results should be
      ///   easy for humans to read (vs. easy for parsers to parse)
      /// \param debug [in] Whether to write verbose debug output to
      ///   err
      DistTsqrVerifier (const Teuchos::RCP< MessengerBase< Scalar > >& scalarComm,
			const std::vector<int>& seed,
			const std::string& scalarTypeName,
			std::ostream& out,
			std::ostream& err,
			const bool testFactorExplicit,
			const bool testFactorImplicit,
			const bool humanReadable,
			const bool debug) :
	gen_ (seed), 
	scalarComm_ (scalarComm),
	scalarTypeName_ (scalarTypeName), 
	out_ (out), 
	err_ (err),
	testFactorExplicit_ (testFactorExplicit),
	testFactorImplicit_ (testFactorImplicit),
	humanReadable_ (humanReadable), 
	debug_ (debug)
      {}

      /// \brief Constructor, with default seed value
      ///
      /// This constructor sets a default seed (for the pseudorandom
      /// number generator), which is the same seed (0,0,0,1) each
      /// time.
      ///
      /// \param scalarComm [in/out] Communicator object over which to
      ///   test.
      /// \param scalarTypeName [in] Human-readable name of the Scalar
      ///   template type parameter
      /// \param out [out] Output stream to which to write results
      /// \param err [out] Output stream to which to write any
      ///   debugging outputs (if applicable) or errors
      /// \param testFactorExplicit [in] Whether to test 
      ///   DistTsqr::factorExplicit()
      /// \param testFactorImplicit [in] Whether to test 
      ///   DistTsqr::factor() and DistTsqr::explicit_Q()
      /// \param humanReadable [in] Whether printed results should be
      ///   easy for humans to read (vs. easy for parsers to parse)
      /// \param debug [in] Whether to write verbose debug output to
      ///   err
      DistTsqrVerifier (const Teuchos::RCP< MessengerBase< Scalar > >& scalarComm,
			const std::string& scalarTypeName,
			std::ostream& out,
			std::ostream& err,
			const bool testFactorExplicit,
			const bool testFactorImplicit,
			const bool humanReadable,
			const bool debug) :
	scalarComm_ (scalarComm),
	scalarTypeName_ (scalarTypeName), 
	out_ (out), 
	err_ (err),
	testFactorExplicit_ (testFactorExplicit),
	testFactorImplicit_ (testFactorImplicit),
	humanReadable_ (humanReadable), 
	debug_ (debug)
      {}

      void 
      getSeed (std::vector<int>& seed) const
      {
	gen_.getSeed (seed);
      }

      /// \brief Run the DistTsqr accuracy test
      ///
      /// \param numCols [in] Number of columns in the matrix to test.
      ///   Number of rows := (# MPI processors) * ncols.
      void 
      verify (const Ordinal numCols)
      {
	using std::endl;

	const bool extraDebug = false;
	const int myRank = scalarComm_->rank();
	if (debug_)
	  {
	    scalarComm_->barrier();
	    if (myRank == 0)
	      err_ << "Verifying DistTsqr:" << endl;
	    scalarComm_->barrier();
	  }

	// Generate test problem.
	Matrix< Ordinal, Scalar > A_local, Q_local, R;
	testProblem (A_local, Q_local, R, numCols);
	if (debug_)
	  {
	    scalarComm_->barrier();
	    if (myRank == 0)
	      err_ << "-- Generated test problem." << endl;
	    scalarComm_->barrier();
	  }

	// Set up TSQR implementation.
	DistTsqr< Ordinal, Scalar > par (scalarComm_);
	if (extraDebug && debug_)
	  {
	    scalarComm_->barrier();
	    if (myRank == 0)
	      err_ << "-- All MPI process(es) have initialized their "
		"DistTsqr object." << endl << endl;
	  }

	// Test DistTsqr::factor() and DistTsqr::explicit_Q().
	if (testFactorImplicit_)
	  {
	    // Factor the matrix A (copied into R, which will be
	    // overwritten on output)
	    typedef typename DistTsqr< Ordinal, Scalar >::FactorOutput 
	      factor_output_type;
	    factor_output_type factorOutput = par.factor (R.view());
	    if (debug_)
	      {
		scalarComm_->barrier();
		if (myRank == 0)
		  err_ << "-- Finished DistTsqr::factor" << endl;
	      }
	    // Compute the explicit Q factor
	    par.explicit_Q (numCols, Q_local.get(), Q_local.lda(), factorOutput);
	    if (debug_)
	      {
		scalarComm_->barrier();
		if (myRank == 0)
		  err_ << "-- Finished DistTsqr::explicit_Q" << endl;
	      }
	    // Verify the factorization
	    result_type result = 
	      global_verify (numCols, numCols, A_local.get(), A_local.lda(),
			     Q_local.get(), Q_local.lda(), R.get(), R.lda(), 
			     scalarComm_.get());
	    if (debug_)
	      {
		scalarComm_->barrier();
		if (myRank == 0)
		  err_ << "-- Finished global_verify" << endl;
	      }
	    reportResults ("DistTsqr", numCols, result);
	  }

	// Test DistTsqr::factorExplicit()
	if (testFactorExplicit_)
	  {
	    // Factor the matrix and compute the explicit Q factor, both
	    // in a single operation.
	    par.factorExplicit (R.view(), Q_local.view());
	    if (debug_)
	      {
		scalarComm_->barrier();
		if (myRank == 0)
		  err_ << "-- Finished DistTsqr::factorExplicit" << endl;
	      }
	    // Verify the factorization
	    result_type result = 
	      global_verify (numCols, numCols, A_local.get(), A_local.lda(),
			     Q_local.get(), Q_local.lda(), R.get(), R.lda(), 
			     scalarComm_.get());
	    if (debug_)
	      {
		scalarComm_->barrier();
		if (myRank == 0)
		  err_ << "-- Finished global_verify" << endl;
	      }
	    reportResults ("DistTsqrRB", numCols, result);
	  }
      }

    private:
      /// Report verification results.  Call on ALL MPI processes, not
      /// just Rank 0.
      ///
      /// \param method [in] String to print before reporting results
      /// \param numCols [in] Number of columns in the matrix tested.
      /// \param result [in] (relative residual, orthogonality)
      void 
      reportResults (const std::string& method, 
		     const Ordinal numCols, 
		     const result_type& result)
      {
	using std::endl;

	const int numProcs = scalarComm_->size();
	const bool printResults = (scalarComm_->rank() == 0);

	if (printResults)
	  {
	    if (humanReadable_)
	      {
		out_ << method << " accuracy results:" << endl
		     << "Scalar type = " << scalarTypeName_ << endl
		     << "Number of columns = " << numCols << endl
		     << "Number of (MPI) processes = " << numProcs << endl
		     << "Relative residual $\\|A - Q*R\\|_2 / \\|A\\|_2$ = " 
		     << result.first << endl
		     << "Relative orthogonality $\\|I - Q^T*Q\\|_2$ = " 
		     << result.second << endl;
	      }
	    else
	      {
		out_ << method
		     << "," << scalarTypeName_
		     << "," << numCols
		     << "," << numProcs
		     << "," << result.first
		     << "," << result.second
		     << endl;
	      }
	  }
      }

      void 
      testProblem (Matrix< Ordinal, Scalar >& A_local,
		   Matrix< Ordinal, Scalar >& Q_local,
		   Matrix< Ordinal, Scalar >& R,
		   const Ordinal numCols)
      {
	const Ordinal numRowsLocal = numCols;

	// A_local: Space for the matrix A to factor -- local to each
	//   processor.
	//
	// A_global: Global matrix (only nonempty on Proc 0); only
	//   used temporarily.
	Matrix< Ordinal, Scalar > A_global;

	// This modifies A_local on all procs, and A_global on Proc 0.
	par_tsqr_test_problem (gen_, A_local, A_global, numCols, scalarComm_);

	// Copy the test problem input into R, since the factorization
	// will overwrite it in place with the final R factor.
	R.reshape (numCols, numCols);
	R.copy (A_local);

	// Prepare space in which to construct the explicit Q factor
	// (local component on this processor)
	Q_local.reshape (numRowsLocal, numCols);
	Q_local.fill (Scalar(0));
      }
    };


    /// \class DistTsqrBenchmarker
    /// \brief Template version of DistTsqr performance test
    ///
    template< class Ordinal, class Scalar, class TimerType >
    class DistTsqrBenchmarker {
      TSQR::Random::NormalGenerator< Ordinal, Scalar > gen_;
      Teuchos::RCP< MessengerBase< Scalar > > scalarComm_;
      Teuchos::RCP< MessengerBase< double > > doubleComm_; 
      std::string scalarTypeName_;
      const bool humanReadable_, debug_;
      std::ostream& out_;
      std::ostream& err_;
      bool testFactorExplicit_,	testFactorImplicit_;

    public:
      typedef Ordinal ordinal_type;
      typedef Scalar scalar_type;
      typedef typename ScalarTraits< scalar_type >::magnitude_type magnitude_type;
      typedef TimerType timer_type;

      /// \brief Constructor, with custom seed value
      ///
      /// \param scalarComm [in/out] Communicator object over which 
      ///   to test.
      /// \param doubleComm [in/out] Communicator object for doubles,
      ///   used for finding the min and max of timing results over
      ///   all the MPI processes.
      /// \param seed [in] 4-element vector; the random seed input of
      ///   TSQR::Random::NormalGenerator (which see, since there are
      ///   restrictions on the set of valid seeds)
      /// \param scalarTypeName [in] Human-readable name of the Scalar
      ///   template type parameter
      /// \param out [out] Output stream to which to write results
      /// \param err [out] Output stream to which to write any
      ///   debugging outputs (if applicable) or errors
      /// \param testFactorExplicit [in] Whether to test 
      ///   DistTsqr::factorExplicit()
      /// \param testFactorImplicit [in] Whether to test 
      ///   DistTsqr::factor() and DistTsqr::explicit_Q()
      /// \param humanReadable [in] Whether printed results should be
      ///   easy for humans to read (vs. easy for parsers to parse)
      /// \param debug [in] Whether to write verbose debug output to
      ///   err
      DistTsqrBenchmarker (const Teuchos::RCP< MessengerBase< Scalar > >& scalarComm,
			   const Teuchos::RCP< MessengerBase< double > >& doubleComm,
			   const std::vector<int>& seed,
			   const std::string& scalarTypeName,
			   std::ostream& out,
			   std::ostream& err,
			   const bool testFactorExplicit,
			   const bool testFactorImplicit,
			   const bool humanReadable,
			   const bool debug) :
	gen_ (seed), 
	scalarComm_ (scalarComm),
	doubleComm_ (doubleComm),
	scalarTypeName_ (scalarTypeName), 
	out_ (out), 
	err_ (err),
	testFactorExplicit_ (testFactorExplicit),
	testFactorImplicit_ (testFactorImplicit),
	humanReadable_ (humanReadable), 
	debug_ (debug)
      {
      }

      /// \brief Constructor, with default seed value
      ///
      /// This constructor sets a default seed (for the pseudorandom
      /// number generator), which is the same seed (0,0,0,1) each
      /// time.
      ///
      /// \param scalarComm [in/out] Communicator object over which 
      ///   to test.
      /// \param doubleComm [in/out] Communicator object for doubles,
      ///   used for finding the min and max of timing results over
      ///   all the MPI processes.
      /// \param scalarTypeName [in] Human-readable name of the Scalar
      ///   template type parameter
      /// \param out [out] Output stream to which to write results
      /// \param err [out] Output stream to which to write any
      ///   debugging outputs (if applicable) or errors
      /// \param testFactorExplicit [in] Whether to test 
      ///   DistTsqr::factorExplicit()
      /// \param testFactorImplicit [in] Whether to test 
      ///   DistTsqr::factor() and DistTsqr::explicit_Q()
      /// \param humanReadable [in] Whether printed results should be
      ///   easy for humans to read (vs. easy for parsers to parse)
      /// \param debug [in] Whether to write verbose debug output to
      ///   err
      DistTsqrBenchmarker (const Teuchos::RCP< MessengerBase< Scalar > >& scalarComm,
			   const Teuchos::RCP< MessengerBase< double > >& doubleComm,
			   const std::string& scalarTypeName,
			   std::ostream& out,
			   std::ostream& err,
			   const bool testFactorExplicit,
			   const bool testFactorImplicit,
			   const bool humanReadable,
			   const bool debug) :
	scalarComm_ (scalarComm),
	doubleComm_ (doubleComm),
	scalarTypeName_ (scalarTypeName), 
	out_ (out), 
	err_ (err),
	testFactorExplicit_ (testFactorExplicit),
	testFactorImplicit_ (testFactorImplicit),
	humanReadable_ (humanReadable), 
	debug_ (debug)
      {
      }

      void 
      getSeed (std::vector<int>& seed) const
      {
	gen_.getSeed (seed);
      }

      /// \brief Run the DistTsqr benchmark
      ///
      /// \param numTrials [in] Number of times to repeat the computation
      ///   in a single timing run
      /// \param numCols [in] Number of columns in the matrix to test.
      ///   Number of rows := (# MPI processors) * ncols
      void 
      benchmark (const int numTrials, const Ordinal numCols)
      {
	using std::endl;

	// Set up test problem.
	Matrix< Ordinal, Scalar > A_local, Q_local, R;
	testProblem (A_local, Q_local, R, numCols);

	// Set up TSQR implementation.
	DistTsqr< Ordinal, Scalar > par (scalarComm_);

	if (testFactorImplicit_)
	  {
	    std::string timerName ("DistTsqr");

	    // Benchmark DistTsqr (factor() and explicit_Q()) for
	    // numTrials trials.
	    timer_type timer (timerName);
	    timer.start();
	    for (int trialNum = 0; trialNum < numTrials; ++trialNum)
	      {
		// Factor the matrix A (copied into R, which will be
		// overwritten on output)
		typedef typename DistTsqr< Ordinal, Scalar >::FactorOutput 
		  factor_output_type;
		factor_output_type factorOutput = par.factor (R.view());

		// Compute the explicit Q factor
		par.explicit_Q (numCols, Q_local.get(), Q_local.lda(), factorOutput);
	      }
	    // Cumulative timing on this MPI process.
	    // "Cumulative" means the elapsed time of numTrials executions.
	    const double localCumulativeTiming = timer.stop();

	    // reportResults() must be called on all processes, since this
	    // figures out the min and max timings over all processes.
	    reportResults (timerName, numTrials, numCols, localCumulativeTiming);
	  }

	if (testFactorExplicit_)
	  {
	    std::string timerName ("DistTsqrRB");

	    // Benchmark DistTsqr::factorExplicit() for numTrials trials.
	    timer_type timer (timerName);
	    timer.start();
	    for (int trialNum = 0; trialNum < numTrials; ++trialNum)
	      {
		par.factorExplicit (R.view(), Q_local.view());
	      }
	    // Cumulative timing on this MPI process.
	    // "Cumulative" means the elapsed time of numTrials executions.
	    const double localCumulativeTiming = timer.stop();

	    // Per-invocation timings.  localTimings were computed on
	    // this MPI process; globalTimings are statistical
	    // summaries of those over all MPI processes.  We only
	    // collect that data for factorExplicit().
	    std::vector< TimeStats > localTimings;
	    std::vector< TimeStats > globalTimings;
	    par.getFactorExplicitTimings (localTimings);
	    for (std::vector< TimeStats >::size_type k = 0; k < localTimings.size(); ++k)
	      globalTimings.push_back (TimeStats::globalTimeStats (doubleComm_.get(), localTimings[k]));
	    std::vector< std::string > timingLabels;
	    par.getFactorExplicitTimingLabels (timingLabels);

	    if (humanReadable_)
	      out_ << timerName << " per-invocation benchmark results:" << endl;

	    for (std::vector< std::string >::size_type k = 0; k < timingLabels.size(); ++k)
	      {
		out_ << "  " << timingLabels[k] << endl;
		globalTimings[k].print (out_, humanReadable_);
	      }

	    reportResults (timerName, numTrials, numCols, localCumulativeTiming);
	  }


      }

    private:
      /// Report results to the output stream.
      ///
      /// \param method [in] String to print before reporting results
      /// \param numTrials [in] Number of times to repeat the computation
      ///   in a single timing run
      /// \param numCols [in] Number of columns in the matrix to test.
      ///   Number of rows := (# MPI processors) * ncols
      /// \param timing [in] Total benchmark time, as measured on this
      ///   MPI process.  This may differ on each process; we report
      ///   the min and the max.
      ///
      /// \warning Call on ALL MPI processes, not just Rank 0!
      void 
      reportResults (const std::string& method,
		     const int numTrials,
		     const ordinal_type numCols,
		     const double localTiming)
      {
	using std::endl;

	// Find min and max timing over all MPI processes
	TimeStats localStats;
	localStats.update (localTiming);
	TimeStats globalStats = TimeStats::globalTimeStats (doubleComm_.get(), localStats);

	// Only Rank 0 prints the final results.
	const bool printResults = (doubleComm_->rank() == 0);
	if (printResults)
	  {
	    if (humanReadable_)
	      {
		out_ << method << " cumulative benchmark results (total time over all trials):" << endl
		     << "Scalar type = " << scalarTypeName_ << endl
		     << "Number of columns = " << numCols << endl
		     << "Number of (MPI) processes = " << doubleComm_->size() << endl
		     << "Number of trials = " << numTrials << endl
		     << "Min timing (in seconds) = " << globalStats.min() << endl
		     << "Mean timing (in seconds) = " << globalStats.mean() << endl
		     << "Max timing (in seconds) = " << globalStats.max() << endl;
	      }
	    else
	      {
		out_ << method << endl
		     << "," << scalarTypeName_ << endl
		     << "," << numCols << endl
		     << "," << doubleComm_->size() << endl
		     << "," << numTrials << endl
		     << "," << globalStats.min() << endl
		     << "," << globalStats.mean() << endl
		     << "," << globalStats.max() << endl;
	      }
	  }
      }

      void 
      testProblem (Matrix< Ordinal, Scalar >& A_local,
		   Matrix< Ordinal, Scalar >& Q_local,
		   Matrix< Ordinal, Scalar >& R,
		   const Ordinal numCols)
      {
	const Ordinal numRowsLocal = numCols;

	// A_local: Space for the matrix A to factor -- local to each
	//   processor.
	//
	// A_global: Global matrix (only nonempty on Proc 0); only
	//   used temporarily.
	Matrix< Ordinal, Scalar > A_global;

	// This modifies A_local on all procs, and A_global on Proc 0.
	par_tsqr_test_problem (gen_, A_local, A_global, numCols, scalarComm_);

	// Copy the test problem input into R, since the factorization
	// will overwrite it in place with the final R factor.
	R.reshape (numCols, numCols);
	R.copy (A_local);

	// Prepare space in which to construct the explicit Q factor
	// (local component on this processor)
	Q_local.reshape (numRowsLocal, numCols);
	Q_local.fill (Scalar(0));
      }

      /// Make sure that timer_type satisfies the TimerType concept.
      ///
      static void
      conceptChecks () 
      {
	verifyTimerConcept< timer_type >();
      }
    };


  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_DistTest_hpp
