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

#ifndef __TSQR_Test_ParTest_hpp
#define __TSQR_Test_ParTest_hpp

#include <Tsqr_Config.hpp>
#include <Tsqr_Random_NormalGenerator.hpp>

#include <Tsqr_generateStack.hpp>
#include <Tsqr_MessengerBase.hpp>
#include <Tsqr_DistTsqr.hpp>
#include <Tsqr_GlobalVerify.hpp>
#include <Tsqr_printGlobalMatrix.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <Tsqr_verifyTimerConcept.hpp>

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
      MessengerBase< Scalar >* const scalarComm_;
      std::string scalarTypeName_;
      const bool humanReadable_, debug_;
      std::ostream& out_;
      std::ostream& err_;

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
      /// \param humanReadable [in] Whether printed results should be
      ///   easy for humans to read (vs. easy for parsers to parse)
      /// \param debug [in] Whether to write verbose debug output to
      ///   err
      /// \param out [out] Output stream to which to write results
      /// \param err [out] Output stream to which to write any
      ///   debugging outputs (if applicable) or errors
      DistTsqrVerifier (MessengerBase< Scalar >* const scalarComm,
			const std::vector<int>& seed,
			const std::string& scalarTypeName,
			const bool humanReadable = true,
			const bool debug = false,
			std::ostream& out = std::cout,
			std::ostream& err = std::cerr) :
	gen_ (seed), 
	scalarComm_ (scalarComm),
	scalarTypeName_ (scalarTypeName), 
	humanReadable_ (humanReadable), 
	debug_ (debug), 
	out_ (out), 
	err_ (err)
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
      /// \param humanReadable [in] Whether printed results should be
      ///   easy for humans to read (vs. easy for parsers to parse)
      /// \param debug [in] Whether to write verbose debug output to
      ///   err
      /// \param out [out] Output stream to which to write results
      /// \param err [out] Output stream to which to write any
      ///   debugging outputs (if applicable) or errors
      DistTsqrVerifier (MessengerBase< Scalar >* const scalarComm,
			const std::string& scalarTypeName,
			const bool humanReadable = true,
			const bool debug = false,
			std::ostream& out = std::cout,
			std::ostream& err = std::cerr) :
	scalarComm_ (scalarComm),
	scalarTypeName_ (scalarTypeName), 
	humanReadable_ (humanReadable), 
	debug_ (debug), 
	out_ (out), 
	err_ (err)
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
	const Ordinal numRowsLocal = numCols;

	// A_local: Space for the matrix A to factor -- local to each
	//   processor.
	// A_global: Global matrix (only nonempty on Proc 0)
	Matrix< Ordinal, Scalar > A_local, A_global;
	// This modifies A_local on all procs, and A_global as well on
	// Proc 0.
	par_tsqr_test_problem (gen_, A_local, A_global, numCols, scalarComm_);

	if (debug_)
	  {
	    scalarComm_->barrier();
	    if (myRank == 0)
	      err_ << "-- Generated test problem." << endl;
	    scalarComm_->barrier();
	  }

	// Copy the test problem input into R, since the factorization will
	// overwrite it place with the final R factor.
	Matrix< Ordinal, Scalar > R (numCols, numCols);
	R.copy (A_local);

	if (debug_)
	  {
	    scalarComm_->barrier();
	    if (myRank == 0)
	      err_ << "-- Finished copying test problem input into (local) R." << endl;
	  }

	// Set up TSQR.
	DistTsqr< Ordinal, Scalar > par (scalarComm_);
	if (extraDebug && debug_)
	  {
	    scalarComm_->barrier();
	    if (myRank == 0)
	      err_ << "-- All MPI process(es) have initialized their "
		"DistTsqr object." << endl << endl;
	  }

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

	// Prepare space in which to construct the explicit Q factor
	// (local component on this processor)
	Matrix< Ordinal, Scalar > Q_local (numRowsLocal, numCols);

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
	  global_verify (numRowsLocal, numCols, A_local.get(), A_local.lda(),
			 Q_local.get(), Q_local.lda(), R.get(), R.lda(), 
			 scalarComm_);
	if (debug_)
	  {
	    scalarComm_->barrier();
	    if (myRank == 0)
	      err_ << "-- Finished global_verify" << endl;
	  }
	if (myRank == 0)
	  reportResults (numCols, result);
      }

    private:
      /// Report verification results.  Call on ALL MPI processes, not
      /// just Rank 0.
      ///
      /// \param numCols [in] Number of columns in the matrix tested.
      /// \param result [in] (relative residual, orthogonality)
      void 
      reportResults (const Ordinal numCols, const result_type& result)
      {
	using std::endl;

	const int numProcs = scalarComm_->size();
	const bool printResults = (scalarComm_->rank() == 0);

	if (printResults)
	  {
	    if (humanReadable_)
	      {
		out_ << "DistTsqr accuracy results:" << endl
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
		out_ << "DistTsqr" << endl
		     << "," << scalarTypeName_ << endl
		     << "," << numCols << endl
		     << "," << numProcs << endl
		     << "," << result.first << endl
		     << "," << result.second << endl;
	      }
	  }
      }
    };


    /// \class DistTsqrBenchmarker
    /// \brief Template version of DistTsqr performance test
    ///
    template< class Ordinal, class Scalar, class TimerType >
    class DistTsqrBenchmarker {
      TSQR::Random::NormalGenerator< Ordinal, Scalar > gen_;
      MessengerBase< Scalar >* const scalarComm_;
      MessengerBase< double >* const doubleComm_; 
      std::string scalarTypeName_;
      const bool humanReadable_, debug_;
      std::ostream& out_;
      std::ostream& err_;

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
      /// \param humanReadable [in] Whether printed results should be
      ///   easy for humans to read (vs. easy for parsers to parse)
      /// \param debug [in] Whether to write verbose debug output to
      ///   err
      /// \param out [out] Output stream to which to write results
      /// \param err [out] Output stream to which to write any
      ///   debugging outputs (if applicable) or errors
      DistTsqrBenchmarker (MessengerBase< Scalar >* const scalarComm,
			   MessengerBase< double >* const doubleComm,
			   const std::vector<int>& seed,
			   const std::string& scalarTypeName,
			   const bool humanReadable = true,
			   const bool debug = false,
			   std::ostream& out = std::cout,
			   std::ostream& err = std::cerr) :
	gen_ (seed), 
	scalarComm_ (scalarComm),
	doubleComm_ (doubleComm),
	scalarTypeName_ (scalarTypeName), 
	humanReadable_ (humanReadable), 
	debug_ (debug), 
	out_ (out), err_ (err)
      {
	verifyTimerConcept< timer_type >();
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
      /// \param humanReadable [in] Whether printed results should be
      ///   easy for humans to read (vs. easy for parsers to parse)
      /// \param debug [in] Whether to write verbose debug output to
      ///   err
      /// \param out [out] Output stream to which to write results
      /// \param err [out] Output stream to which to write any
      ///   debugging outputs (if applicable) or errors
      DistTsqrBenchmarker (MessengerBase< Scalar >* const scalarComm,
			   MessengerBase< double >* const doubleComm,
			   const std::string& scalarTypeName,
			   const bool humanReadable = true,
			   const bool debug = false,
			   std::ostream& out = std::cout,
			   std::ostream& err = std::cerr) :
	scalarComm_ (scalarComm),
	doubleComm_ (doubleComm),
	scalarTypeName_ (scalarTypeName), 
	humanReadable_ (humanReadable), 
	debug_ (debug), 
	out_ (out), 
	err_ (err)
      {
	verifyTimerConcept< timer_type >();
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

	const Ordinal numRowsLocal = numCols;

	// A_local: Space for the matrix A to factor -- local to each
	//   processor.
	// A_global: Global matrix (only nonempty on Proc 0)
	Matrix< Ordinal, Scalar > A_local, A_global;
	// This modifies A_local on all procs, and A_global as well on
	// Proc 0.
	par_tsqr_test_problem (gen_, A_local, A_global, numCols, scalarComm_);
	// Copy the test problem input into R, since the factorization will
	// overwrite it place with the final R factor.
	Matrix< Ordinal, Scalar > R (numCols, numCols);
	R.copy (A_local);

	// Prepare space in which to construct the explicit Q factor
	// (local component on this processor)
	Matrix< Ordinal, Scalar > Q_local (numRowsLocal, numCols);

	// Set up TSQR.
	DistTsqr< Ordinal, Scalar > par (scalarComm_);

	// Benchmark DistTsqr for numTrials trials.
	//
	// Name of timer doesn't matter here; we only need the timing.
	timer_type timer("DistTsqr");
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
	const double timing = timer.stop();
	// reportResults() must be called on all processes, since this
	// figures out the min and max timings over all processes.
	reportResults (numTrials, numCols, timing);
      }

    private:
      /// Report results to the output stream.
      ///
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
      reportResults (const int numTrials,
		     const ordinal_type numCols,
		     const double localTiming)
      {
	using std::endl;

	// Find min and max timing over all MPI processes
	const double minTiming = doubleComm_->globalMin (localTiming);
	const double maxTiming = doubleComm_->globalMax (localTiming);

	// Only Rank 0 prints the final results.
	const bool printResults = (doubleComm_->rank() == 0);
	if (printResults)
	  {
	    if (humanReadable_)
	      {
		out_ << "DistTsqr benchmark results:" << endl
		     << "Scalar type = " << scalarTypeName_ << endl
		     << "Number of columns = " << numCols << endl
		     << "Number of (MPI) processes = " << doubleComm_->size() << endl
		     << "Number of trials = " << numTrials << endl
		     << "Min timing (in seconds) = " << minTiming << endl
		     << "Max timing (in seconds) = " << maxTiming << endl;
	      }
	    else
	      {
		out_ << "DistTsqr" << endl
		     << "," << scalarTypeName_ << endl
		     << "," << numCols << endl
		     << "," << doubleComm_->size() << endl
		     << "," << numTrials << endl
		     << "," << minTiming << endl
		     << "," << maxTiming << endl;
	      }
	  }
      }
    };


  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_ParTest_hpp
