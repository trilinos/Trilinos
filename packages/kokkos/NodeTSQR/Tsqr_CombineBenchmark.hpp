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

#ifndef __TSQR_Test_CombineBenchmark_hpp
#define __TSQR_Test_CombineBenchmark_hpp

#include <Tsqr_Config.hpp>
#include <Tsqr_CombineBenchmarker.hpp>
#include <Tsqr_CombineDefault.hpp>
#include <Tsqr_CombineNative.hpp>
#ifdef HAVE_TSQR_FORTRAN
#  include <Tsqr_CombineFortran.hpp>
#endif // HAVE_TSQR_FORTRAN

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>


namespace TSQR {
  namespace Test {

    template<class CombineType, class TimerType>
    static void
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

      TEST_FOR_EXCEPTION(numRows < 1 || numCols < 1, std::invalid_argument,
			 "The test matrix must have a positive number of rows "
			 "and columns, but you specified numRows = " << numRows 
			 << " and numCols = " << numCols << ".");
      TEST_FOR_EXCEPTION(cacheBlockNumTrials < 1, std::invalid_argument,
			 "The number of trials for the cache block benchmark "
			 "must be positive, but you specified cacheBlockNum"
			 "Trials = " << cacheBlockNumTrials << ".");
      TEST_FOR_EXCEPTION(pairNumTrials < 1, std::invalid_argument,
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

      out << combineTypeName 
	  << "," << "R1R2"
	  << "," << dataTypeName
	  << "," << (2*numCols)
	  << "," << numCols
	  << "," << pairNumTrials
	  << "," << (averageTimings ? results.first / static_cast<double>(pairNumTrials) : results.first);
      if (printAdditionalData)
	out << "," << additionalData;
      out << endl;
      out << combineTypeName 
	  << "," << "RA"
	  << "," << dataTypeName
	  << "," << numRows
	  << "," << numCols
	  << "," << cacheBlockNumTrials
	  << "," << (averageTimings ? results.second / static_cast<double>(cacheBlockNumTrials) : results.second);
      if (printAdditionalData)
	out << "," << additionalData;
      out << endl;
    }

    template<class Scalar, class TimerType>
    static void
    benchmarkAllCombineTypes (std::ostream& out,
			      std::vector<int>& iseed,
			      const std::string& dataTypeName,
			      const int numRows,
			      const int numCols,
			      const int numTrials,
			      const bool calibrate,
			      const bool averageTimings,
			      const double timerResolution,
			      const std::string& additionalData,
			      const bool debug)
    {
      using std::cerr;
      using std::endl;

      TEST_FOR_EXCEPTION(numRows < 1 || numCols < 1, std::invalid_argument,
			 "The test matrix must have a positive number of rows "
			 "and columns, but you specified numRows = " << numRows 
			 << " and numCols = " << numCols << ".");
      TEST_FOR_EXCEPTION(! calibrate && numTrials < 1, std::invalid_argument,
			 "The number of trials must be positive, but you "
			 "specified numTrials = " << numTrials << ".");
      TEST_FOR_EXCEPTION(timerResolution <= static_cast<double>(0), 
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
      int pairNumTrials = numTrials;
      int cacheBlockNumTrials = numTrials;
      if (calibrate)
	{ // We calibrate the number of trials using the default
	  // Combine implementation.  We don't expect CombineNative or
	  // CombineFortran to be much faster than that.  
	  if (debug)
	    cerr << "Calibrating..." << endl;

	  // Calibrater gets the timer resolution.
	  typedef CombineDefault<int, Scalar> combine_type;
	  typedef CombineBenchmarker<int, Scalar, combine_type, TimerType> 
	    calibrater_type;
	  calibrater_type c (timerResolution);

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
	}

      {
	typedef CombineNative<int, Scalar> combine_type;
	std::string combineTypeName ("Native");
	benchmarkCombineType<combine_type, TimerType> (out, iseed, 
						       dataTypeName, 
						       combineTypeName, 
						       numRows, 
						       numCols, 
						       cacheBlockNumTrials,
						       pairNumTrials,
						       averageTimings,
						       additionalData);
      }
#ifdef HAVE_TSQR_FORTRAN
      {
	typedef CombineFortran<Scalar> combine_type;
	std::string combineTypeName ("Fortran");
	benchmarkCombineType<combine_type, TimerType> (out, iseed, 
						       dataTypeName, 
						       combineTypeName, 
						       numRows, 
						       numCols, 
						       cacheBlockNumTrials,
						       pairNumTrials,
						       averageTimings,
						       additionalData);
      }
#endif // HAVE_TSQR_FORTRAN
      {
	typedef CombineDefault< int, Scalar > combine_type;
	std::string combineTypeName ("Default");
	benchmarkCombineType<combine_type, TimerType> (out, iseed, 
						       dataTypeName, 
						       combineTypeName, 
						       numRows, 
						       numCols, 
						       cacheBlockNumTrials,
						       pairNumTrials,
						       averageTimings,
						       additionalData);
      }
    }


    template<class TimerType>
    static void
    benchmarkAllCombineTypesAndScalars (std::ostream& out,
					std::vector<int>& iseed,
					const int numRows,
					const int numCols,
					const int numTrials,
					const bool testReal,
					const bool testComplex,
					const bool calibrate,
					const bool averageTimings,
					const std::string& additionalData,
					const bool debug)
    {
      using std::cerr;
      using std::endl;
      using std::string;
      string dataTypeName;

      TEST_FOR_EXCEPTION(numRows < 1 || numCols < 1, std::invalid_argument,
			 "The test matrix must have a positive number of rows "
			 "and columns, but you specified numRows = " << numRows 
			 << " and numCols = " << numCols << ".");
      TEST_FOR_EXCEPTION(! calibrate && numTrials < 1, std::invalid_argument,
			 "The number of trials must be positive, but you "
			 "specified numTrials = " << numTrials << ".");

      // Compute timer resolution.
      const double timerResolution = computeTimerResolution<TimerType> ();
      if (debug)
	cerr << "Timer resolution: " << timerResolution << " seconds" << endl;

      if (testReal)
	{
	  dataTypeName = "float";
	  benchmarkAllCombineTypes<float, TimerType> (out, iseed, 
						      dataTypeName, 
						      numRows, 
						      numCols, 
						      numTrials, 
						      calibrate,
						      averageTimings,
						      timerResolution,
						      additionalData, 
						      debug);
	  dataTypeName = "double";
	  benchmarkAllCombineTypes<double, TimerType> (out, iseed, 
						       dataTypeName, 
						       numRows, 
						       numCols, 
						       numTrials, 
						       calibrate,
						       averageTimings,
						       timerResolution,
						       additionalData,
						       debug);
	}
      if (testComplex)
	{
#ifdef HAVE_TSQR_COMPLEX
	  using std::complex;

	  dataTypeName = "complex<float>";
	  benchmarkAllCombineTypes<complex<float>, TimerType> (out, iseed, 
							       dataTypeName, 
							       numRows, 
							       numCols, 
							       numTrials, 
							       calibrate,
							       averageTimings,
							       timerResolution,
							       additionalData,
							       debug);
	  dataTypeName = "complex<double>";
	  benchmarkAllCombineTypes<complex<double>, TimerType> (out, iseed, 
								dataTypeName, 
								numRows, 
								numCols, 
								numTrials,
								calibrate,
								averageTimings,
								timerResolution,
								additionalData,
								debug);
#else // Don't HAVE_TSQR_COMPLEX
	  throw std::logic_error("TSQR not built with complex arithmetic support");
#endif // HAVE_TSQR_COMPLEX
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
    /// \param numRows [in] Number of rows in the cache block A.
    /// \param numCols [in] Number of columns in the cache block A,
    ///   and number of rows and columns in the upper triangular
    ///   matrices R, R1, and R2.
    /// \param testReal [in] Whether to test real-arithmetic routines.
    /// \param testComplex [in] Whether to test complex-arithmetic
    ///   routines.
    /// \param numTrials [in] If calibrate is false: the number of
    ///   trials to run each of the benchmarks.  Ignored if calibrate
    ///   is true.
    /// \param calibrate [in] Whether to calibrate the number of
    ///   trials according to the computed timer resolution.
    /// \param averageTimings [in] Whether to print average (true)
    ///   or cumulative (false) timings over all trials.
    /// \param seed [in] If useSeedValues is false, ignored; else, the
    ///   four-integer seed for the random number generator.  See the
    ///   documentation of LAPACK's _LARNV routines for requirements.
    /// \param useSeedValues [in] Whether seed (see above) is read.
    /// \param additionalFieldNames [in] Field names for additional
    ///   data to print after each row.
    /// \param additionalData [in] Additional data to print after each
    ///   row.  Same number of entries as additionalFieldNames.
    /// \param printFieldNames [in] Whether to print a "%" - commented
    ///   row of comma-delimited field names before the first row of
    ///   data.
    /// \param debug [in] Whether to print copious debugging output
    ///   to stderr.
    template<class TimerType>
    void
    benchmarkCombine (std::ostream& out,
		      const int numRows,
		      const int numCols,
		      const bool testReal,
		      const bool testComplex, 
		      const int numTrials,
		      const bool calibrate,
		      const bool averageTimings,
		      std::vector<int>& seed,
		      const bool useSeedValues,
		      const std::string& additionalFieldNames,
		      const std::string& additionalData,
		      const bool printFieldNames,
		      const bool debug)
    {
      TEST_FOR_EXCEPTION(numRows < 1 || numCols < 1, std::invalid_argument,
			 "The test matrix must have a positive number of rows "
			 "and columns, but you specified numRows = " << numRows 
			 << " and numCols = " << numCols << ".");
      TEST_FOR_EXCEPTION(! calibrate && numTrials < 1, std::invalid_argument,
			 "The number of trials must be positive, but you "
			 "specified numTrials = " << numTrials << ".");

      if (! useSeedValues)
	{
	  // Default seed values.
	  if (seed.size() < 4)
	    seed.resize (4);
	  seed[0] = 0;
	  seed[1] = 0;
	  seed[2] = 0;
	  seed[3] = 1;
	}

      if (printFieldNames)
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
	  if (printFieldNames && ! additionalFieldNames.empty())
	    // The additionalFieldNames string should be a
	    // comma-delimited list of additional field name(s).
	    out << "," << additionalFieldNames;
	  out << std::endl;
	}
      benchmarkAllCombineTypesAndScalars<TimerType> (out, seed, 
						     numRows, 
						     numCols, 
						     numTrials, 
						     testReal,
						     testComplex, 
						     calibrate,
						     averageTimings,
						     additionalData,
						     debug);
    }

  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_CombineBenchmark_hpp
