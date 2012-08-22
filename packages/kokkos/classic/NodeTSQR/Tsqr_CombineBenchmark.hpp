//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef __TSQR_Test_CombineBenchmark_hpp
#define __TSQR_Test_CombineBenchmark_hpp

#include <Tsqr_ConfigDefs.hpp>
#include <Tsqr_CombineBenchmarker.hpp>
#include <Tsqr_CombineDefault.hpp>
#include <Tsqr_CombineNative.hpp>
#ifdef HAVE_KOKKOSCLASSIC_TSQR_FORTRAN
#  include <Tsqr_CombineFortran.hpp>
#endif // HAVE_KOKKOSCLASSIC_TSQR_FORTRAN

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>


namespace TSQR {
  namespace Test {

    /// \class CombineBenchmarkParameters
    /// \brief Parameters for the TSQR::Combine benchmarks.
    /// \author Mark Hoemmen
    ///
    /// numRows: Number of rows in the cache block A.
    ///
    /// numCols: Number of columns in the cache block A,
    ///   and number of rows and columns in the upper triangular
    ///   matrices R, R1, and R2.
    ///
    /// testReal: Whether to test real-arithmetic routines.
    ///
    /// testComplex: Whether to test complex-arithmetic routines.
    ///
    /// numTrials: If calibrate is false: the number of
    ///   trials to run each of the benchmarks.  Ignored if calibrate
    ///   is true.
    ///
    /// calibrate: Whether to calibrate the number of
    ///   trials according to the computed timer resolution.
    ///
    /// averageTimings: Whether to print average (true)
    ///   or cumulative (false) timings over all trials.
    /// 
    /// strictPerfTests: Whether to require the ratio of CombineNative
    ///   run time to CombineDefault run time to be less than
    ///   allowance.  "Require" means that we throw an exception (and
    ///   the test fails) otherwise.  CombineFortran is tested
    ///   similarly, if applicable.
    ///
    /// allowance: Allowed slowdown factor for strictPerfTests (if
    ///   applicable).
    ///
    /// seed: If useSeedValues is false, ignored; else, the
    ///   four-integer seed for the random number generator.  See the
    ///   documentation of LAPACK's _LARNV routines for requirements.
    ///
    /// useSeedValues: Whether seed (see above) is read.
    ///
    /// additionalFieldNames: Field names for any additional
    ///   data to print after each row.  May be an empty string,
    ///   in which case the number of additional fields is zero
    ///   and no additional data is printed after each row.
    ///
    /// additionalData: Any additional data to print after each row.
    ///   Same number of additional data per row as fields in
    ///   additionalFieldNames.
    ///
    /// printFieldNames: Whether to print a "%" - commented row of
    ///   comma-delimited field names before the first row of data.
    ///
    /// debug: Whether to print copious debugging output to stderr.
    ///
    struct CombineBenchmarkParameters {
      int numRows;
      int numCols;
      bool testReal;
      bool testComplex;
      int numTrials;
      bool calibrate;
      bool averageTimings;
      bool strictPerfTests;
      double allowance;
      std::vector<int> seed;
      bool useSeedValues;
      std::string additionalFieldNames;
      std::string additionalData;
      bool printFieldNames;
      bool debug;
    };

    template<class CombineType, class TimerType>
    static std::vector<double>
    benchmarkCombineType (std::ostream& out,
			  std::vector<int>& iseed,
			  const std::string& dataTypeName,
			  const std::string& combineTypeName,
			  const typename CombineType::ordinal_type numRows,
			  const typename CombineType::ordinal_type numCols,
			  const int cacheBlockNumTrials,
			  const int pairNumTrials,
			  const bool averageTimings,
			  const std::string& additionalData)
    {
      using std::endl;

      typedef typename CombineType::ordinal_type ordinal_type;
      typedef typename CombineType::scalar_type scalar_type;
      typedef typename CombineType::magnitude_type magnitude_type;
      typedef CombineBenchmarker<ordinal_type, scalar_type, CombineType, TimerType> 
	benchmarker_type;

      TEUCHOS_TEST_FOR_EXCEPTION(cacheBlockNumTrials < 1, std::invalid_argument,
			 "The number of trials for the cache block benchmark "
			 "must be positive, but you specified cacheBlockNum"
			 "Trials = " << cacheBlockNumTrials << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(pairNumTrials < 1, std::invalid_argument,
			 "The number of trials for the pair benchmark must be "
			 "positive, but you specified pairNumTrials = "
			 << pairNumTrials << ".");

      benchmarker_type b (iseed);
      std::pair<double, double> results;
      results.first = 
	b.benchmarkPair (numCols, pairNumTrials);
      results.second = 
	b.benchmarkCacheBlock (numRows, numCols, cacheBlockNumTrials);

      // Whether or not we should print the "additional data"
      // (originally supplied at command-line invocation of this
      // benchmark) after the benchmark results.  The additional data
      // option makes it easier to write parsers for benchmark
      // results, since we can include data that are known outside the
      // benchmark (when invoking the benchmark as an executable), but
      // not known (easily or at all) inside the benchmark.  A good
      // example would be environment variables, like OMP_NUM_THREADS,
      // or (for a benchmark that uses MPI, which this is not) the
      // number of MPI processes per node ("ppn").
      const bool printAdditionalData = (! additionalData.empty());

      const double pairTime = averageTimings ? 
	results.first / static_cast<double>(pairNumTrials) : 
	results.first;
      const double cacheBlockTime = averageTimings ? 
	results.second / static_cast<double>(cacheBlockNumTrials) : 
	results.second;

      out << combineTypeName 
	  << "," << "R1R2"
	  << "," << dataTypeName
	  << "," << (2*numCols)
	  << "," << numCols
	  << "," << pairNumTrials
	  << "," << pairTime;
      if (printAdditionalData)
	out << "," << additionalData;
      out << endl;
      out << combineTypeName 
	  << "," << "RA"
	  << "," << dataTypeName
	  << "," << numRows
	  << "," << numCols
	  << "," << cacheBlockNumTrials
	  << "," << cacheBlockTime;
      if (printAdditionalData)
	out << "," << additionalData;
      out << endl;

      std::vector<double> timings (2);
      timings[0] = pairTime;
      timings[1] = cacheBlockTime;
      return timings;
    }

    template<class Scalar, class TimerType>
    static void
    benchmarkAllCombineTypes (std::ostream& out,
			      const std::string& dataTypeName,
			      CombineBenchmarkParameters& params,
			      const double timerResolution)
    {
      using std::cerr;
      using std::endl;
      const bool debug = params.debug;
      const int numRows = params.numRows;
      const int numCols = params.numCols;

      TEUCHOS_TEST_FOR_EXCEPTION(timerResolution <= static_cast<double>(0), 
			 std::invalid_argument,
			 "The timer resolution must be a positive number, "
			 "but you specified timerResolution = " 
			 << timerResolution << ".");

      // If no calibration is performed, then the number of trials is
      // the same for both the cache block [R; A] benchmark and the
      // pair [R1; R2] benchmark.  Otherwise, we calibrate the number
      // of trials for each separately.  This is because we expect the
      // [R1; R2] benchmark to take much less time than the [R; A]
      // benchmark, so [R1; R2] should have more trials, in order to
      // get comparable timing accuracy without requiring too many [R;
      // A] trials.
      int pairNumTrials = params.numTrials;
      int cacheBlockNumTrials = params.numTrials;
      if (params.calibrate)
	{ // We calibrate the number of trials using the default
	  // Combine implementation.  We don't expect CombineNative or
	  // CombineFortran to be much faster than that.  
	  if (debug)
	    cerr << "Calibrating..." << endl;

	  // Calibrater gets the timer resolution.
	  typedef CombineDefault<int, Scalar> combine_type;
	  typedef CombineBenchmarker<int, Scalar, combine_type, TimerType> 
	    benchmarker_type;
	  benchmarker_type c (timerResolution, params.seed);

	  // Accuracy factor of 1000 gives us 3 digits of timer accuracy.
	  const double accuracyFactor = static_cast<double> (1000);

	  // Number of trials for factor_pair() and apply_pair().
	  std::pair<int, double> result;
	  result = c.calibratePair (numCols, accuracyFactor);
	  if (debug)
	    {
	      cerr << "- Pair number of trials: " << result.first << endl;
	      cerr << "- Pair calibration time: " << result.second << endl;
	    }
	  pairNumTrials = result.first;

	  // Number of trials for factor_inner() and apply_inner().
	  result = c.calibrateCacheBlock (numRows, numCols, accuracyFactor);
	  if (debug)
	    {
	      cerr << "- Cache block number of trials: " << result.first << endl;
	      cerr << "- Cache block calibration time: " << result.second << endl;
	    }
	  cacheBlockNumTrials = result.first;

	  // Store the updated PRNG seed in the benchmark parameters.
	  c.getSeed (params.seed);
	}

      // Always benchmark CombineDefault.  We use its timings as the
      // standard by which the other Combine implementations' timings
      // are compared.  The returned vector contains two timings: for
      // [R1; R2], and for [R; A], in that order.
      std::vector<double> defaultTimings;
      {
	typedef CombineDefault< int, Scalar > combine_type;
	std::string combineTypeName ("Default");
	defaultTimings = 
	  benchmarkCombineType<combine_type, TimerType> (out, params.seed,
							 dataTypeName, 
							 combineTypeName, 
							 numRows, 
							 numCols, 
							 cacheBlockNumTrials,
							 pairNumTrials,
							 params.averageTimings,
							 params.additionalData);
      }

      // If we're doing strict performance tests, then CombineNative
      // (and CombineFortran, if applicable) may be no slower than the
      // given allowance factor times CombineDefault's time.  For now,
      // we only look at cache block performance, since that is where
      // most of the time should be going.
      std::vector<double> nativeTimings;
      {
	typedef CombineNative<int, Scalar> combine_type;
	std::string combineTypeName ("Native");
	nativeTimings = 
	  benchmarkCombineType<combine_type, TimerType> (out, params.seed, 
							 dataTypeName, 
							 combineTypeName, 
							 numRows, 
							 numCols, 
							 cacheBlockNumTrials,
							 pairNumTrials,
							 params.averageTimings,
							 params.additionalData);
	const double slowdown = nativeTimings[1] / defaultTimings[1];
	const bool tooSlow = slowdown > params.allowance;
	// FIXME (mfh 24 May 2011) Replace std::runtime_error with a
	// more appropriately named exception.
	TEUCHOS_TEST_FOR_EXCEPTION(params.strictPerfTests && tooSlow, 
			   std::runtime_error,
			   "CombineNative is too slow!  For cache block "
			   "benchmark with numRows=" << numRows << " and numCols="
			   << numCols << ", CombineNative time (= " 
			   << nativeTimings[1] << ") / CombineDefault time (= "
			   << defaultTimings[1] << ") = " << slowdown 
			   << " > the allowed fraction " << params.allowance 
			   << ".");
      }

#ifdef HAVE_KOKKOSCLASSIC_TSQR_FORTRAN
      std::vector<double> fortranTimings;
      {
	typedef CombineFortran<Scalar> combine_type;
	std::string combineTypeName ("Fortran");
	fortranTimings = 
	  benchmarkCombineType<combine_type, TimerType> (out, params.seed, 
							 dataTypeName, 
							 combineTypeName, 
							 numRows, 
							 numCols, 
							 cacheBlockNumTrials,
							 pairNumTrials,
							 params.averageTimings,
							 params.additionalData);
	const double slowdown = fortranTimings[1] / defaultTimings[1];
	const bool tooSlow = slowdown > params.allowance;
	// FIXME (mfh 24 May 2011) Replace std::runtime_error with a
	// more appropriately named exception.
	TEUCHOS_TEST_FOR_EXCEPTION(params.strictPerfTests && tooSlow, 
			   std::runtime_error,
			   "CombineFortran is too slow!  For cache block "
			   "benchmark with numRows=" << numRows << " and numCols="
			   << numCols << ", CombineFortran time (= " 
			   << fortranTimings[1] << ") / CombineDefault time (= "
			   << defaultTimings[1] << ") = " << slowdown 
			   << " > the allowed fraction " << params.allowance 
			   << ".");
      }
#endif // HAVE_KOKKOSCLASSIC_TSQR_FORTRAN
    }


    template<class TimerType>
    static void
    benchmarkAllCombineTypesAndScalars (std::ostream& out,
					CombineBenchmarkParameters& params)
    {
      using std::cerr;
      using std::endl;
      using std::string;
      const bool debug = params.debug;

      // Compute timer resolution.
      const double timerResolution = computeTimerResolution<TimerType> ();
      if (debug)
	cerr << "Timer resolution: " << timerResolution << " seconds" << endl;

      string dataTypeName;
      if (params.testReal)
	{
	  dataTypeName = "float";
	  benchmarkAllCombineTypes<float, TimerType> (out, dataTypeName, 
						      params, timerResolution);
	  dataTypeName = "double";
	  benchmarkAllCombineTypes<double, TimerType> (out, dataTypeName,
						       params, timerResolution);
	}
      if (params.testComplex)
	{
#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
	  using std::complex;

	  dataTypeName = "complex<float>";
	  benchmarkAllCombineTypes<complex<float>, TimerType> (out, dataTypeName,
							       params, timerResolution);
	  dataTypeName = "complex<double>";
	  benchmarkAllCombineTypes<complex<double>, TimerType> (out, dataTypeName,
								params, timerResolution);

#else // Don't HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
	  throw std::logic_error("TSQR not built with complex arithmetic support");
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
	}
    }					

    /// \fn benchmarkCombine
    /// \brief Benchmark TSQR::Combine, using a timer of type TimerType.
    /// \author Mark Hoemmen
    ///
    /// Benchmarks test cache block and pair operations for all
    /// Combine implementations, over all Scalar types (modulated by
    /// testReal and testComplex).
    ///
    /// \param out [out] Output stream to which to write results.
    /// \param params [in/out] Benchmark parameters.
    template<class TimerType>
    void
    benchmarkCombine (std::ostream& out,
		      CombineBenchmarkParameters& params)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(params.numRows < 1 || params.numCols < 1, 
			 std::invalid_argument,
			 "The test matrix must have a positive number of rows "
			 "and columns, but you specified numRows = " 
			 << params.numRows << " and numCols = "
			 << params.numCols << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(! params.calibrate && params.numTrials < 1, 
			 std::invalid_argument,
			 "Since you specified no calibration is to be performed, "
			 "the number of trials must be positive, but you specified "
			 "numTrials = " << params.numTrials << ".");

      if (! params.useSeedValues)
	{ // Fill in default seed values.
	  if (params.seed.size() < 4)
	    params.seed.resize (4);
	  params.seed[0] = 0;
	  params.seed[1] = 0;
	  params.seed[2] = 0;
	  params.seed[3] = 1;
	}

      if (params.printFieldNames)
	{
	  // The row of field names begins with a '%' character, in
	  // order to help out the benchmark results parser.
	  out << "%" << "method"
	      << "," << "kernel"
	      << "," << "scalarType"
	      << "," << "numRows"
	      << "," << "numCols"
	      << "," << "numTrials"
	      << "," << "timing";
	  if (params.printFieldNames && ! params.additionalFieldNames.empty())
	    // The additionalFieldNames string should be a
	    // comma-delimited list of additional field name(s).
	    out << "," << params.additionalFieldNames;
	  out << std::endl;
	}
      benchmarkAllCombineTypesAndScalars<TimerType> (out, params);
    }

  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_CombineBenchmark_hpp
