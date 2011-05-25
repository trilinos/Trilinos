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

#ifndef __Tsqr_CombineBenchmarker_hpp
#define __Tsqr_CombineBenchmarker_hpp

#include <Tsqr_Config.hpp>
#include <Tsqr_Random_NormalGenerator.hpp>
#include <Tsqr_Random_MatrixGenerator.hpp>
#include <Tsqr_verifyTimerConcept.hpp>

#include <Tsqr_ApplyType.hpp>
#include <Tsqr_Matrix.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <Tsqr_Util.hpp>

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>


namespace TSQR {
  namespace Test {

    /// \fn computeTimerResolution
    /// \brief Compute resolution in seconds of the TimerType timer.
    ///
    /// The timer resolution is the smallest time interval for which
    /// TimerType reports a nonzero time.  We measure the timer's
    /// resolution in seconds with a timing loop that does some fake
    /// sequential floating-point work.  We keep increasing the
    /// number of iterations in the timing loop until the timer says
    /// that the loop took nonzero time, and return that nonzero
    /// time.
    template<class TimerType>
    double
    computeTimerResolution () 
    {
      typedef TimerType timer_type;
      timer_type timer ("Timer resolution");

      // Warmup run for the timer.
      for (int warmup = 0; warmup < 5; ++warmup)
	{
	  timer.start();
	  (void) timer.stop();
	}

      // Keep a count of the total number of times timer.stop() is
      // called (once per outer loop iteration).  If bigger than
      // maxCount, assume that the timer is broken and bail out.
      // This prevents an infinite loop.  Since numTrials is
      // multiplied by 10 each time around, maxCount should be more
      // than the base-10 log of the maximum acceptable timer
      // resolution, relative to the time taken by a single
      // floating-point increment operation.
      //
      // For example, suppose that the time per flop is 10^{-9} s,
      // and the timer resolution is 1 s (which is particularly
      // coarse).  We start with numTrials=10 and multiply by 10
      // each time.  Then, it may take 8 or 9 iterations (may have
      // to exceed the timer resolution, or just meet it) so the
      // timing loop exceeds the timer resolution.  Thus, maxCount=9
      // = log10(timer resolution / time per flop) is just enough.
      // Adding a bit more for wiggle room is a good idea.
      //
      // If your machine is exceptionally fast or your compiler
      // autoparallelizes the timing loop, you may have to increase
      // maxCount.  We set it to 20 as a generous upper bound.
      const size_t maxCount = 20;
      size_t count = 0;

      // Don't let numTrials loop around.  (We're multiplying it by
      // 10 each time, so it won't take long for this to happen,
      // especially if the timer granularity is too coarse.)
      const size_t maxNumTrials = std::numeric_limits<size_t>::max() / 10;

      // Do some fake work.  Multiply the number of trials by 10
      // each time, so that resolution is expressed as an order of
      // decimal magnitude.
      size_t numTrials = 1;
      double theTime;
      do {
	const double eps = Teuchos::ScalarTraits<double>::eps();
	double fake = Teuchos::ScalarTraits<double>::one();
	numTrials *= 10;
	++count;

	// The timing loop.
	timer.start();
	for (size_t trial = 0; trial < numTrials; ++trial)
	  fake += eps;
	theTime = timer.stop();

      } while (theTime == 0 && count < maxCount && numTrials < maxNumTrials);

      if (theTime == 0)
	{
	  std::ostringstream os;
	  os << "Maximum number of loops " << maxCount << " exceeded when "
	    "computing timer resolution.  Largest timing loop length tried: " 
	     << numTrials << ".";
	  throw std::logic_error (os.str());
	}
      else if (numTrials >= maxNumTrials)
	{
	  std::ostringstream os;
	  os << "Maximum number of timing loop iterations " << maxNumTrials 
	     << " exceeded when computing timer resolution.  Largest timing "
	    "loop length tried: " << numTrials << ".";
	  throw std::logic_error (os.str());
	}
      else
	return theTime;
    }

    /// \class CombineBenchmarker
    /// \brief Benchmark TSQR::Combine, and also calibrate number of trials.
    /// \author Mark Hoemmen
    ///
    /// This class both calibrates the number of trials for
    /// benchmarking the given Combine implementation for a given
    /// problem size, and does the actual benchmarking given a number
    /// of trials.
    ///
    /// Calibration works as follows:
    /// 1. For the given timer type, approximate the timer's resolution.
    /// 2. For a given set of matrix dimensions, determine how many
    ///    times the Combine benchmarks need to be run in order to get
    ///    the desired accuracy of the timing results.
    ///
    /// Template parameters:
    /// - Ordinal: Combine's index type
    /// - Scalar: The type of matrix entries
    /// - CombineType: The TSQR::Combine implementation
    /// - TimerType: Type of the timer
    ///
    /// TimerType must pass the test given by 
    /// \c TSQR::Test::verifyTimerConcept<TimerType>().
    template<class Ordinal, class Scalar, class CombineType, class TimerType>
    class CombineBenchmarker {
    public:
      typedef Ordinal ordinal_type;
      typedef Scalar scalar_type;
      typedef CombineType combine_type;
      typedef TimerType timer_type;

    private:
      typedef Teuchos::ScalarTraits<scalar_type> STS;
      typedef typename STS::magnitudeType magnitude_type;
      typedef Teuchos::ScalarTraits<magnitude_type> STM;
      typedef TSQR::Random::NormalGenerator<ordinal_type, scalar_type> normgen_type;
      typedef TSQR::Random::MatrixGenerator<ordinal_type, scalar_type, normgen_type> matgen_type;
      typedef Matrix<ordinal_type, scalar_type> matrix_type;

    public:

      /// \brief Constructor with user-specified seed.
      ///
      /// \param timerRes [in] Resolution in seconds of the TimerType
      ///   timer.
      ///
      /// \param iseed [in] Seed for the pseudorandom number
      ///   generator.  See the LAPACK documentation (for the _LARNV
      ///   routines) for requirements on the seed values.
      ///
      CombineBenchmarker (const double timerRes,
			  const std::vector<int>& iseed) : 
	normGenS_ (iseed), 
	normGenM_ (iseed),
	timerResolution_ (timerRes)
      {
	TSQR::Test::verifyTimerConcept<timer_type> ();
      }

      /// \brief Constructor with user-specified seed; computes timer resolution.
      ///
      /// \param iseed [in] Seed for the pseudorandom number
      ///   generator.  See the LAPACK documentation (for the _LARNV
      ///   routines) for requirements on the seed values.
      ///
      CombineBenchmarker (const std::vector<int>& iseed) : 
	normGenS_ (iseed), 
	normGenM_ (iseed),
	timerResolution_ (computeTimerResolution<timer_type> ())
      {}

      /// \brief Constructor with default seed.
      ///
      /// \param timerRes [in] Resolution in seconds of the TimerType
      ///   timer.
      CombineBenchmarker (const double timerRes) :
	timerResolution_ (timerRes)
      {
	TSQR::Test::verifyTimerConcept<timer_type> ();
      }


      //! Constructor with default seed that computes timer resolution.
      CombineBenchmarker () :
	timerResolution_ (computeTimerResolution<timer_type> ())
      {}

      //! Get the current pseudorandom number generator seed.  
      void 
      getSeed (std::vector<int>& iseed) const
      {
	normGenS_.getSeed (iseed);
      }

      //! Smallest time interval (in seconds) which TimerType can measure.
      double 
      timerResolution() const {
	return timerResolution_;
      }

      /// \brief Estimate number of trials for TSQR::Combine on first cache block.
      ///
      /// TSQR::Combine implementations use factor_first() to factor
      /// the "first" numRows by numCols cache block A.  They then
      /// compute the explicit Q factor of the cache block A by
      /// calling apply_first() on a numRows by numCols matrix Q
      /// consisting of the first numCols columns of the identity
      /// matrix.  This benchmark measures the time taken by
      /// repeatedly invoking factor_first() and then apply_first() in
      /// this way.  It returns the number of trials for such
      /// invocations which takes at least as much time as the given
      /// accuracy factor, times the timer resolution.
      /// 
      /// \param numRows [in] Number of rows in the cache block A.
      /// \param numCols [in] Number of columns in the cache block A.
      /// \param accuracyFactor [in] The computed number of trials
      ///   takes at least as much time as this accuracy factor times
      ///   the timer resolution. 
      ///
      /// \return Number of trials, and cumulative time of the
      ///   benchmark over that many trials.  (The second value lets
      ///   you recycle that timing, so you don't have to run that
      ///   benchmark again.)
      std::pair<int, double>
      calibrateFirst (const Ordinal numRows,
		      const Ordinal numCols,
		      const double accuracyFactor)
      {
	if (numRows == 0 || numCols == 0)
	  throw std::invalid_argument("Calibrating timings is impossible for "
				      "a matrix with either zero rows or zero "
				      "columns.");
	else if (accuracyFactor < 0)
	  throw std::invalid_argument("Accuracy factor for Combine numTrials "
				      "calibration must be nonnegative.");
	// Random matrix generator.
	matgen_type matGen (normGenS_);

	// Generate a random cache block A.
	matrix_type A (numRows, numCols);
	std::vector<magnitude_type> sigmas (numCols);
	randomSingularValues (sigmas, numCols);
	matGen.fill_random_svd (numRows, numCols, A.get(), A.lda(), &sigmas[0]);

	// A place to put the Q factor.
	matrix_type Q (numRows, numCols);
	Q.fill (STS::zero());
	for (Ordinal j = 0; j < numCols; ++j)
	  Q(j,j) = STS::one();

	// TAU array (Householder reflector scaling factors).
	std::vector<Scalar> tau (numCols);
	// Work space array for factorization and applying the Q factor.
	std::vector<Scalar> work (numCols);

	// The Combine instance to benchmark.
	combine_type combiner; 

	// A few warmup runs just to avoid timing anomalies.
	const int numWarmupRuns = 3;
	for (int warmupRun = 0; warmupRun < numWarmupRuns; ++warmupRun)
	  {
	    combiner.factor_first (numRows, numCols, A.get(), A.lda(),
				   &tau[0], &work[0]);
	    combiner.apply_first (ApplyType("N"), numRows, numCols, numCols,
				  A.get(), A.lda(), &tau[0],
				  Q.get(), Q.lda(), &work[0]);
	  }

	// How much time numTrials runs must take in order for
	// numTrials to be considered sufficiently large.
	const double minAcceptableTime = accuracyFactor * timerResolution();

	timer_type timer ("Combine first");

	// The actual timing runs.  Repeat, doubling numTrials each
	// time, until the resulting timing is at least timerRes *
	// accuracyFactor.  Also, we mandate somewhat arbitrarily that
	// numTrials >= 4; this gives us some buffer against timer
	// variability.  Finally, don't let numTrials loop around.
	// (We're doubling it each time, so it won't take long for
	// this to happen, especially if something is wrong with the
	// benchmark and it's taking zero time, or if accuracyFactor
	// is too large, or if the timer resolution is too large.)
	const int maxNumTrials = std::numeric_limits<int>::max() / 2;
	double theTime;
	int numTrials = 2;
	do {
	  numTrials *= 2; // First value of numTrials is 4.
	  timer.start();
	  for (int trial = 0; trial < numTrials; ++trial)
	    {
	      combiner.factor_first (numRows, numCols, A.get(), A.lda(),
				     &tau[0], &work[0]);
	      combiner.apply_first (ApplyType("N"), numRows, numCols, numCols,
				    A.get(), A.lda(), &tau[0],
				    Q.get(), Q.lda(), &work[0]);
	    }
	  theTime = timer.stop();
	} while (theTime < minAcceptableTime && numTrials < maxNumTrials);

	return std::make_pair (numTrials, theTime);
      }

      /// \brief Benchmark TSQR::Combine on first cache block.
      ///
      /// TSQR::Combine implementations use factor_first() to factor
      /// the "first" numRows by numCols cache block A.  They then
      /// compute the explicit Q factor of the cache block A by
      /// calling apply_first() on a numRows by numCols matrix Q
      /// consisting of the first numCols columns of the identity
      /// matrix.  This benchmark measures the time taken by
      /// repeatedly invoking factor_first() and then apply_first() in
      /// this way, for numTrials trials.  
      /// 
      /// \param numRows [in] Number of rows in the cache block A.
      /// \param numCols [in] Number of columns in the cache block A.
      /// \param numTrials [in] Number of timing loops.
      ///
      /// \return Cumulative time over numTrials trials.
      double
      benchmarkFirst (const Ordinal numRows,
		      const Ordinal numCols,
		      const int numTrials)
      {
	if (numRows == 0 || numCols == 0)
	  throw std::invalid_argument("Benchmarking does not make sense for "
				      "a matrix with either zero rows or zero "
				      "columns.");
	TEST_FOR_EXCEPTION(numTrials < 1, std::invalid_argument,
			   "The number of trials must be positive, but "
			   "numTrials = " << numTrials << ".");

	// Random matrix generator.
	matgen_type matGen (normGenS_);

	// Generate a random cache block A.
	matrix_type A (numRows, numCols);
	std::vector<magnitude_type> sigmas (numCols);
	randomSingularValues (sigmas, numCols);
	matGen.fill_random_svd (numRows, numCols, A.get(), A.lda(), &sigmas[0]);

	// A place to put the Q factor.
	matrix_type Q (numRows, numCols);
	Q.fill (STS::zero());
	for (Ordinal j = 0; j < numCols; ++j)
	  Q(j,j) = STS::one();

	// TAU array (Householder reflector scaling factors).
	std::vector<Scalar> tau (numCols);
	// Work space array for factorization and applying the Q factor.
	std::vector<Scalar> work (numCols);

	// The Combine instance to benchmark.
	combine_type combiner; 

	// A few warmup runs just to avoid timing anomalies.
	const int numWarmupRuns = 3;
	for (int warmupRun = 0; warmupRun < numWarmupRuns; ++warmupRun)
	  {
	    combiner.factor_first (numRows, numCols, A.get(), A.lda(),
				   &tau[0], &work[0]);
	    combiner.apply_first (ApplyType("N"), numRows, numCols, numCols,
				  A.get(), A.lda(), &tau[0],
				  Q.get(), Q.lda(), &work[0]);
	  }
	//
	// The actual timing runs.
	//
	timer_type timer ("Combine first");
	timer.start();
	for (int trial = 0; trial < numTrials; ++trial)
	  {
	    combiner.factor_first (numRows, numCols, A.get(), A.lda(),
				   &tau[0], &work[0]);
	    combiner.apply_first (ApplyType("N"), numRows, numCols, numCols,
				  A.get(), A.lda(), &tau[0],
				  Q.get(), Q.lda(), &work[0]);
	  }
	return timer.stop();
      }

      /// \brief Estimate number of trials for TSQR::Combine on [R; A];
      ///
      /// TSQR::Combine implementations use factor_inner() to factor a
      /// vertical stack [R; A] of a square R factor and a numRows by
      /// numCols cache block A.  They then compute the explicit Q
      /// factor of [R; A] by calling apply_inner() on a ((numRows +
      /// numCols) by numCols) matrix Q consisting of the first
      /// numCols columns of the identity matrix.  This benchmark
      /// measures the time taken by repeatedly invoking
      /// factor_inner() and then apply_inner() in this way.  It
      /// returns the number of trials for such invocations which
      /// takes at least as much time as the given accuracy factor,
      /// times the timer resolution.
      /// 
      /// \param numRows [in] Number of rows in the cache block A.
      /// \param numCols [in] Number of columns in the cache block A, 
      ///   and number of rows and columns in R.
      /// \param accuracyFactor [in] The computed number of trials
      ///   takes at least as much time as this accuracy factor times
      ///   the timer resolution. 
      ///
      /// \return Number of trials, and cumulative time of the
      ///   benchmark over that many trials.  (The second value lets
      ///   you recycle that timing, so you don't have to run that
      ///   benchmark again.)
      std::pair<int, double>
      calibrateCacheBlock (const Ordinal numRows,
			   const Ordinal numCols,
			   const double accuracyFactor)
      {
	if (numRows == 0 || numCols == 0)
	  throw std::invalid_argument("Calibrating timings is impossible for "
				      "a matrix with either zero rows or zero "
				      "columns.");
	else if (accuracyFactor < 0)
	  throw std::invalid_argument("Accuracy factor for Combine numTrials "
				      "calibration must be nonnegative.");
	// Random matrix generator.
	matgen_type matGen (normGenS_);

	// Generate a random R factor first.
	matrix_type R (numCols, numCols);
	std::vector<magnitude_type> sigmas (numCols);
	randomSingularValues (sigmas, numCols);
	matGen.fill_random_R (numCols, R.get(), R.lda(), &sigmas[0]);

	// Now generate a random cache block.
	matrix_type A (numRows, numCols);
	randomSingularValues (sigmas, numCols);
	matGen.fill_random_svd (numRows, numCols, A.get(), A.lda(), &sigmas[0]);

	// A place to put the Q factor.
	matrix_type Q (numRows + numCols, numCols);
	Q.fill (STS::zero());
	for (Ordinal j = 0; j < numCols; ++j)
	  Q(j,j) = STS::one();

	// TAU array (Householder reflector scaling factors).
	std::vector<Scalar> tau (numCols);
	// Work space array for factorization and applying the Q factor.
	std::vector<Scalar> work (numCols);

	// The Combine instance to benchmark.
	combine_type combiner; 

	// A few warmup runs just to avoid timing anomalies.
	const int numWarmupRuns = 3;
	for (int warmupRun = 0; warmupRun < numWarmupRuns; ++warmupRun)
	  {
	    combiner.factor_inner (numRows, numCols, R.get(), R.lda(),
				   A.get(), A.lda(), &tau[0], &work[0]);
	    combiner.apply_inner (ApplyType("N"), numRows, numCols, numCols,
				  A.get(), A.lda(), &tau[0], 
				  &Q(0, 0), Q.lda(),
				  &Q(numCols, 0), Q.lda(), 
				  &work[0]);
	  }

	// How much time numTrials runs must take in order for
	// numTrials to be considered sufficiently large.
	const double minAcceptableTime = accuracyFactor * timerResolution();

	timer_type timer ("Combine cache block");

	// The actual timing runs.  Repeat, doubling numTrials each
	// time, until the resulting timing is at least timerRes *
	// accuracyFactor.  Also, we mandate somewhat arbitrarily that
	// numTrials >= 4; this gives us some buffer against timer
	// variability.  Finally, don't let numTrials loop around.
	// (We're doubling it each time, so it won't take long for
	// this to happen, especially if something is wrong with the
	// benchmark and it's taking zero time, or if accuracyFactor
	// is too large, or if the timer resolution is too large.)
	const int maxNumTrials = std::numeric_limits<int>::max() / 2;
	double theTime;
	int numTrials = 2;
	do {
	  numTrials *= 2; // First value of numTrials is 4.
	  timer.start();
	  for (int trial = 0; trial < numTrials; ++trial)
	    {
	      combiner.factor_inner (numRows, numCols, R.get(), R.lda(),
				     A.get(), A.lda(), &tau[0], &work[0]);
	      combiner.apply_inner (ApplyType("N"), numRows, numCols, numCols,
				    A.get(), A.lda(), &tau[0], 
				    &Q(0, 0), Q.lda(),
				    &Q(numCols, 0), Q.lda(), 
				    &work[0]);
	    }
	  theTime = timer.stop();
	} while (theTime < minAcceptableTime && numTrials < maxNumTrials);

	return std::make_pair (numTrials, theTime);
      }


      /// \brief Benchmark TSQR::Combine on [R; A];
      ///
      /// TSQR::Combine implementations use factor_inner() to factor a
      /// vertical stack [R; A] of a square R factor and a numRows by
      /// numCols cache block A.  They then compute the explicit Q
      /// factor of [R; A] by calling apply_inner() on a ((numRows +
      /// numCols) by numCols) matrix Q consisting of the first
      /// numCols columns of the identity matrix.  This benchmark
      /// measures the time taken by repeatedly invoking
      /// factor_inner() and then apply_inner() in this way, for
      /// numTrials trials.
      /// 
      /// \param numRows [in] Number of rows in the cache block A.
      /// \param numCols [in] Number of columns in the cache block A, 
      ///   and number of rows and columns in R.
      /// \param numTrials [in] Number of timing loops.
      ///
      /// \return Cumulative time over numTrials trials.
      double
      benchmarkCacheBlock (const Ordinal numRows,
			   const Ordinal numCols,
			   const int numTrials)
      {
	if (numRows == 0 || numCols == 0)
	  throw std::invalid_argument("Benchmarking does not make sense for "
				      "a matrix with either zero rows or zero "
				      "columns.");
	TEST_FOR_EXCEPTION(numTrials < 1, std::invalid_argument,
			   "The number of trials must be positive, but "
			   "numTrials = " << numTrials << ".");

	// Random matrix generator.
	matgen_type matGen (normGenS_);

	// Generate a random R factor first.
	matrix_type R (numCols, numCols);
	std::vector<magnitude_type> sigmas (numCols);
	randomSingularValues (sigmas, numCols);
	matGen.fill_random_R (numCols, R.get(), R.lda(), &sigmas[0]);

	// Now generate a random cache block.
	matrix_type A (numRows, numCols);
	randomSingularValues (sigmas, numCols);
	matGen.fill_random_svd (numRows, numCols, A.get(), A.lda(), &sigmas[0]);

	// A place to put the Q factor.
	matrix_type Q (numRows + numCols, numCols);
	Q.fill (STS::zero());
	for (Ordinal j = 0; j < numCols; ++j)
	  Q(j,j) = STS::one();

	// TAU array (Householder reflector scaling factors).
	std::vector<Scalar> tau (numCols);
	// Work space array for factorization and applying the Q factor.
	std::vector<Scalar> work (numCols);

	// The Combine instance to benchmark.
	combine_type combiner; 

	// A few warmup runs just to avoid timing anomalies.
	const int numWarmupRuns = 3;
	for (int warmupRun = 0; warmupRun < numWarmupRuns; ++warmupRun)
	  {
	    combiner.factor_inner (numRows, numCols, R.get(), R.lda(),
				   A.get(), A.lda(), &tau[0], &work[0]);
	    combiner.apply_inner (ApplyType("N"), numRows, numCols, numCols,
				  A.get(), A.lda(), &tau[0], 
				  &Q(0, 0), Q.lda(),
				  &Q(numCols, 0), Q.lda(), 
				  &work[0]);
	  }
	//
	// The actual timing runs.
	//
	timer_type timer ("Combine cache block");
	timer.start();
	for (int trial = 0; trial < numTrials; ++trial)
	  {
	    combiner.factor_inner (numRows, numCols, R.get(), R.lda(),
				   A.get(), A.lda(), &tau[0], &work[0]);
	    combiner.apply_inner (ApplyType("N"), numRows, numCols, numCols,
				  A.get(), A.lda(), &tau[0], 
				  &Q(0, 0), Q.lda(),
				  &Q(numCols, 0), Q.lda(), 
				  &work[0]);
	  }
	return timer.stop();
      }

      /// \brief Estimate number of trials for TSQR::Combine on [R1; R2].
      ///
      /// TSQR::Combine implementations use factor_pair() to factor a
      /// stack of two square R factors [R1; R2].  They then compute
      /// the explicit Q factor of [R1; R2] by calling apply_pair() on
      /// a (2*numCols by numCols) matrix Q consisting of the first
      /// numCols columns of the 2*numCols by 2*numCols identity
      /// matrix.  This benchmark measures the time taken by
      /// repeatedly invoking factor_pair() and then apply_pair() in
      /// this way.  It returns the number of trials for such
      /// invocations which takes at least as much time as the given
      /// accuracy factor, times the timer resolution.  
      /// 
      /// \param numCols [in] Number of rows and columns in each of R1
      ///   and R2.
      /// \param accuracyFactor [in] The computed number of trials
      ///   takes at least as much time as this accuracy factor times
      ///   the timer resolution. 
      ///
      /// \return Number of trials, and cumulative time of the
      ///   benchmark over that many trials.  (The second value lets
      ///   you recycle that timing, so you don't have to run that
      ///   benchmark again.)
      std::pair<int, double>
      calibratePair (const Ordinal numCols,
		     const double accuracyFactor)
      {
	if (numCols == 0)
	  throw std::invalid_argument("Calibrating timings is impossible for "
				      "a matrix with zero columns.");
	else if (accuracyFactor < 0)
	  throw std::invalid_argument("Accuracy factor for Combine numTrials "
				      "calibration must be nonnegative.");
	// Random matrix generator.
	matgen_type matGen (normGenS_);

	// Generate R1 first.
	matrix_type R1 (numCols, numCols);
	std::vector<magnitude_type> sigmas (numCols);
	randomSingularValues (sigmas, numCols);
	matGen.fill_random_R (numCols, R1.get(), R1.lda(), &sigmas[0]);

	// Now generate R2.
	matrix_type R2 (numCols, numCols);
	randomSingularValues (sigmas, numCols);
	matGen.fill_random_R (numCols, R2.get(), R2.lda(), &sigmas[0]);

	// A place to put the Q factor of [R1; R2].
	matrix_type Q (2*numCols, numCols);
	Q.fill (STS::zero());
	for (Ordinal j = 0; j < numCols; ++j)
	  Q(j,j) = STS::one();

	// TAU array (Householder reflector scaling factors).
	std::vector<Scalar> tau (numCols);
	// Work space array for factorization and applying the Q factor.
	std::vector<Scalar> work (numCols);

	// The Combine instance to benchmark.
	combine_type combiner; 

	// A few warmup runs just to avoid timing anomalies.
	const int numWarmupRuns = 3;
	for (int warmupRun = 0; warmupRun < numWarmupRuns; ++warmupRun)
	  {
	    combiner.factor_pair (numCols, R1.get(), R1.lda(), 
				  R2.get(), R2.lda(),
				  &tau[0], &work[0]);
	    combiner.apply_pair (ApplyType("N"), numCols, numCols, 
				 R2.get(), R2.lda(), &tau[0], 
				 &Q(0, 0), Q.lda(),
				 &Q(numCols, 0), Q.lda(),
				 &work[0]);
	  }

	// How much time numTrials runs must take in order for
	// numTrials to be considered sufficiently large.
	const double minAcceptableTime = accuracyFactor * timerResolution();

	timer_type timer ("Combine pair");

	// The actual timing runs.  Repeat, doubling numTrials each
	// time, until the resulting timing is at least timerRes *
	// accuracyFactor.  Also, we mandate somewhat arbitrarily that
	// numTrials >= 4; this gives us some buffer against timer
	// variability.  Finally, don't let numTrials loop around.
	// (We're doubling it each time, so it won't take long for
	// this to happen, especially if something is wrong with the
	// benchmark and it's taking zero time, or if accuracyFactor
	// is too large, or if the timer resolution is too large.)
	const int maxNumTrials = std::numeric_limits<int>::max() / 2;
	double theTime;
	int numTrials = 2;
	do {
	  numTrials *= 2; // First value of numTrials is 4.
	  timer.start();
	  for (int trial = 0; trial < numTrials; ++trial)
	    {
	      combiner.factor_pair (numCols, R1.get(), R1.lda(), 
				    R2.get(), R2.lda(),
				    &tau[0], &work[0]);
	      combiner.apply_pair (ApplyType("N"), numCols, numCols, 
				   R2.get(), R2.lda(), &tau[0], 
				   &Q(0, 0), Q.lda(),
				   &Q(numCols, 0), Q.lda(),
				   &work[0]);
	    }
	  theTime = timer.stop();
	} while (theTime < minAcceptableTime && numTrials < maxNumTrials);

	return std::make_pair (numTrials, theTime);
      }


      /// \brief Benchmark TSQR::Combine on [R1; R2].
      ///
      /// TSQR::Combine implementations use factor_pair() to factor a
      /// stack of two square R factors [R1; R2].  They then compute
      /// the explicit Q factor of [R1; R2] by calling apply_pair() on
      /// a (2*numCols by numCols) matrix Q consisting of the first
      /// numCols columns of the 2*numCols by 2*numCols identity
      /// matrix.  This benchmark measures the time taken by
      /// repeatedly invoking factor_pair() and then apply_pair() in
      /// this way, for numTrials trials.
      /// 
      /// \param numCols [in] Number of rows and columns in each of R1
      ///   and R2.
      /// \param numTrials [in] Number of timing loops.
      ///
      /// \return Cumulative time over numTrials trials.
      double
      benchmarkPair (const Ordinal numCols,
		     const int numTrials)
      {
	if (numCols == 0)
	  throw std::invalid_argument("Benchmarking does not make sense for "
				      "a matrix with zero columns.");
	TEST_FOR_EXCEPTION(numTrials < 1, std::invalid_argument,
			   "The number of trials must be positive, but "
			   "numTrials = " << numTrials << ".");

	// Random matrix generator.
	matgen_type matGen (normGenS_);

	// Generate R1 first.
	matrix_type R1 (numCols, numCols);
	std::vector<magnitude_type> sigmas (numCols);
	randomSingularValues (sigmas, numCols);
	matGen.fill_random_R (numCols, R1.get(), R1.lda(), &sigmas[0]);

	// Now generate R2.
	matrix_type R2 (numCols, numCols);
	randomSingularValues (sigmas, numCols);
	matGen.fill_random_R (numCols, R2.get(), R2.lda(), &sigmas[0]);

	// A place to put the Q factor of [R1; R2].
	matrix_type Q (2*numCols, numCols);
	Q.fill (STS::zero());
	for (Ordinal j = 0; j < numCols; ++j)
	  Q(j,j) = STS::one();

	// TAU array (Householder reflector scaling factors).
	std::vector<Scalar> tau (numCols);
	// Work space array for factorization and applying the Q factor.
	std::vector<Scalar> work (numCols);

	// The Combine instance to benchmark.
	combine_type combiner; 

	// A few warmup runs just to avoid timing anomalies.
	const int numWarmupRuns = 3;
	for (int warmupRun = 0; warmupRun < numWarmupRuns; ++warmupRun)
	  {
	    combiner.factor_pair (numCols, R1.get(), R1.lda(), 
				  R2.get(), R2.lda(),
				  &tau[0], &work[0]);
	    combiner.apply_pair (ApplyType("N"), numCols, numCols, 
				 R2.get(), R2.lda(), &tau[0], 
				 &Q(0, 0), Q.lda(),
				 &Q(numCols, 0), Q.lda(),
				 &work[0]);
	  }
	//
	// The actual timing runs.
	//
	timer_type timer ("Combine pair");
	timer.start();
	for (int trial = 0; trial < numTrials; ++trial)
	  {
	    combiner.factor_pair (numCols, R1.get(), R1.lda(), 
				  R2.get(), R2.lda(),
				  &tau[0], &work[0]);
	    combiner.apply_pair (ApplyType("N"), numCols, numCols, 
				 R2.get(), R2.lda(), &tau[0], 
				 &Q(0, 0), Q.lda(),
				 &Q(numCols, 0), Q.lda(),
				 &work[0]);
	  }
	return timer.stop();
      }

    private:

      //! Pseudorandom normal(0,1) generator for Scalar values.
      TSQR::Random::NormalGenerator<ordinal_type, scalar_type> normGenS_;

      //! Pseudorandom normal(0,1) generator for magnitude_type values.
      TSQR::Random::NormalGenerator<ordinal_type, magnitude_type> normGenM_;

      //! Timer resolution (in seconds) for TimerType timers.
      double timerResolution_;

      /// \brief Fill sigmas with numValues random singular values.
      ///
      /// \param sigmas [out] Vector of at least numValues entries.
      ///   Resized if necessary.
      /// \param numValues [in] Number of random singular values to
      ///   generate.
      void
      randomSingularValues (std::vector<magnitude_type>& sigmas,
			    const Ordinal numValues)
      {
	// Cast to avoid compiler warnings for signed / unsigned
	// comparisons.
	typedef typename std::vector<magnitude_type>::size_type size_type;
	if (sigmas.size() < static_cast<size_type> (numValues))
	  sigmas.resize (numValues);
      
	// Relative amount by which to perturb each singular value.  The
	// perturbation will be multiplied by a normal(0,1) pseudorandom
	// number drawn from magGen.
	const magnitude_type perturbationFactor = magnitude_type(10) * STM::eps();
	const magnitude_type one = STM::one();
	for (Ordinal k = 0; k < numValues; ++k)
	  {
	    magnitude_type perturbation = perturbationFactor * normGenM_();
	    // If (1 - perturbation) is a small or nonpositive number,
	    // subtract instead.
	    if (one - perturbation <= perturbationFactor)
	      perturbation = -perturbation;
	    sigmas[k] = one - perturbation;
	  }
      }
    };

  } // namespace Test
} // namespace TSQR

#endif // __Tsqr_CombineBenchmarker_hpp
