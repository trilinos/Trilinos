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

#include <Tsqr_Config.hpp>
#include <Tsqr_SeqTest.hpp>

#include <Tsqr_Random_NormalGenerator.hpp>
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
#include <fstream>
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

    /// Test the accuracy of sequential TSQR on an nrows by ncols
    /// matrix (using the given cache block size (in bytes)), and
    /// print the results to stdout.
    template< class Ordinal, class Scalar >
    static void
    verifySeqTsqrTemplate (std::ostream& out,
			   TSQR::Random::NormalGenerator< Ordinal, Scalar >& generator,
			   const std::string& datatype,
			   const std::string& shortDatatype,
			   const Ordinal nrows, 
			   const Ordinal ncols, 
			   const size_t cache_block_size,
			   const bool contiguous_cache_blocks,
			   const bool save_matrices,
			   const bool human_readable,
			   const bool b_debug,
			   const bool printFieldNames)
    {
      typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
      using std::cerr;
      using std::endl;
      using std::pair;
      using std::string;
      using std::vector;

      SequentialTsqr< Ordinal, Scalar > actor (cache_block_size);
      Ordinal numCacheBlocks;

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

      if (save_matrices)
	{
	  string filename = "A_" + shortDatatype + ".txt";
	  if (b_debug)
	    cerr << "-- Saving test problem to \"" << filename << "\"" << endl;
	  std::ofstream fileOut (filename.c_str());
	  print_local_matrix (fileOut, nrows, ncols, A.get(), A.lda());
	  fileOut.close();
	}

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

      // Count the number of cache blocks that factor() will use. 
      // This is only for diagnostic purposes.
      numCacheBlocks = 
	actor.factor_num_cache_blocks (nrows, ncols, A_copy.get(), 
				       A_copy.lda(), contiguous_cache_blocks);
      // In debug mode, report how many cache blocks factor() will use.
      if (b_debug)
	cerr << "-- Number of cache blocks factor() will use: " 
	     << numCacheBlocks << endl << endl;

      // Factor the matrix and compute the explicit Q factor
      typedef typename SequentialTsqr< Ordinal, Scalar >::FactorOutput 
	factor_output_type;
      factor_output_type factorOutput = 
	actor.factor (nrows, ncols, A_copy.get(), A_copy.lda(), 
		      R.get(), R.lda(), contiguous_cache_blocks);
      if (b_debug)
	cerr << "-- Finished SequentialTsqr::factor" << endl;

      if (save_matrices)
	{
	  string filename = "R_" + shortDatatype + ".txt";
	  if (b_debug)
	    cerr << "-- Saving R factor to \"" << filename << "\"" << endl;
	  std::ofstream fileOut (filename.c_str());
	  print_local_matrix (fileOut, ncols, ncols, R.get(), R.lda());
	  fileOut.close();
	}

      actor.explicit_Q (nrows, ncols, A_copy.get(), lda, factorOutput,
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

      if (save_matrices)
	{
	  string filename = "Q_" + shortDatatype + ".txt";
	  if (b_debug)
	    cerr << "-- Saving Q factor to \"" << filename << "\"" << endl;
	  std::ofstream fileOut (filename.c_str());
	  print_local_matrix (fileOut, nrows, ncols, Q.get(), Q.lda());
	  fileOut.close();
	}

      // Print out the R factor
      if (false && b_debug)
	{
	  cerr << endl << "-- R factor:" << endl;
	  print_local_matrix (cerr, ncols, ncols, R.get(), R.lda());
	  cerr << endl;
	}

      // Validate the factorization
      vector< magnitude_type > results =
	local_verify (nrows, ncols, A.get(), lda, Q.get(), ldq, R.get(), ldr);
      if (b_debug)
	cerr << "-- Finished local_verify" << endl;

      // Print the results
      if (human_readable)
	out << "Sequential cache-blocked TSQR:" << endl
	    << "Datatype: " << datatype << endl
	    << "Absolute residual $\\| A - QR \\|_F$: " << results[0] << endl
	    << "Absolute orthogonality $\\| I - Q^* Q \\|_F$: " << results[1] << endl
	    << "Test matrix norm $\\| A \\|_F$: " << results[2] << endl
	    << endl << endl;
      else
	{
	  if (printFieldNames)
	    {
	      const char prefix[] = "%";
	      out << prefix
		  << "method"
		  << ",scalarType"
		  << ",numRows"
		  << ",numCols"
		  << ",cacheBlockSize"
		  << ",numCacheBlocks"
		  << ",contiguousCacheBlocks"
		  << ",absFrobResid"
		  << ",absFrobOrthog"
		  << ",frobA"
		  << endl;
	    }
	  out << "SeqTSQR"
	      << "," << datatype
	      << "," << nrows
	      << "," << ncols
	      << "," << actor.cache_block_size()
	      << "," << numCacheBlocks
	      << "," << contiguous_cache_blocks 
	      << "," << results[0]
	      << "," << results[1]
	      << "," << results[2]
	      << endl;
	}
    }


    void
    verifySeqTsqr (std::ostream& out,
		   const int nrows, 
		   const int ncols, 
		   const size_t cache_block_size,
		   const bool test_complex_arithmetic,
		   const bool save_matrices,
		   const bool contiguous_cache_blocks,
		   const bool printFieldNames,
		   const bool human_readable,
		   const bool b_debug)
    {
      using TSQR::Random::NormalGenerator;
      using std::complex;
      using std::string;
      using std::vector;

      //
      // We do tests one after another, using the seed from the
      // previous test in the current test, so that the pseudorandom
      // streams used by the tests are independent.
      //

      // On output: Seed for the next pseudorandom number generator.
      vector< int > iseed(4);
      string datatype; // name of the current datatype being tested
      string shortDatatype; // one-letter version of datatype

      // First test.  The PRNG seeds itself with a default value.
      // This will be the same each time, so if you want
      // nondeterministic behavior, you should pick the seed values
      // yourself.  Only print field names (if at all) for the first
      // data type tested; field names are only printed if output is
      // not human_readable.
      NormalGenerator< int, float > normgenS;
      datatype = "float";
      shortDatatype = "S";
      verifySeqTsqrTemplate (out, normgenS, datatype, shortDatatype, nrows, ncols, 
			     cache_block_size, contiguous_cache_blocks, 
			     save_matrices, human_readable, b_debug, printFieldNames);
      // Fetch the pseudorandom seed from the previous test.
      normgenS.getSeed (iseed);
      NormalGenerator< int, double > normgenD (iseed);
      // Next test.
      datatype = "double";
      shortDatatype = "D";
      verifySeqTsqrTemplate (out, normgenD, datatype, shortDatatype, nrows, ncols, 
			     cache_block_size, contiguous_cache_blocks, 
			     save_matrices, human_readable, b_debug, false);

      if (test_complex_arithmetic)
	{
	  normgenD.getSeed (iseed);
	  NormalGenerator< int, complex<float> > normgenC (iseed);
	  datatype = "complex<float>";
	  shortDatatype = "C";
	  verifySeqTsqrTemplate (out, normgenC, datatype, shortDatatype, nrows, ncols, 
				 cache_block_size, contiguous_cache_blocks, 
				 save_matrices, human_readable, b_debug, false);
	  normgenC.getSeed (iseed);
	  NormalGenerator< int, complex<double> > normgenZ (iseed);
	  datatype = "complex<double>";
	  shortDatatype = "Z";
	  verifySeqTsqrTemplate (out, normgenZ, datatype, shortDatatype, nrows, ncols, 
				 cache_block_size, contiguous_cache_blocks, 
				 save_matrices, human_readable, b_debug, false);
	}
    }



    template< class Ordinal, class Scalar >
    static void
    verifyLapackTemplate (TSQR::Random::NormalGenerator< Ordinal, Scalar >& generator,
			  const std::string& datatype,
			  const Ordinal nrows, 
			  const Ordinal ncols, 
			  const bool human_readable,
			  const bool b_debug)
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

      if (b_debug)
	{
	  cerr << endl << "-- R factor:" << endl;
	  print_local_matrix (cerr, ncols, ncols, R.get(), R.lda());
	  cerr << endl;
	}

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
      std::vector< magnitude_type > results = 
	local_verify (nrows, ncols, A.get(), lda, Q.get(), ldq, R.get(), ldr);

      // Print the results
      if (human_readable)
	cout << "LAPACK QR (DGEQRF and DORGQR):" << endl
	     << "Datatype: " << datatype << endl
	     << "Absolute residual $\\| A - QR \\|_F$: " << results[0] << endl
	     << "Absolute orthogonality $\\| I - Q^* Q \\|_F$: " << results[1] << endl
	     << "Test matrix norm $\\| A \\|_F$: " << results[2] << endl
	     << endl << endl;
      else
	cout << "LAPACK"
	     << "," << datatype
	     << "," << nrows
	     << "," << ncols
	     << "," << size_t(0) // cache_block_size
	     << "," << false     // contiguous_cache_blocks 
	     << "," << results[0]
	     << "," << results[1]
	     << "," << results[2]
	     << endl;
    }


    void
    verifyLapack (const int nrows, 
		  const int ncols, 
		  const bool test_complex_arithmetic,
		  const bool human_readable,
		  const bool b_debug)
    {
      using TSQR::Random::NormalGenerator;
      using std::complex;
      using std::string;
      using std::vector;

      //
      // We do tests one after another, using the seed from the
      // previous test in the current test, so that the pseudorandom
      // streams used by the tests are independent.
      //

      // On output: Seed for the next pseudorandom number generator.
      vector< int > iseed(4);
      string datatype; // name of the current datatype being tested

      // First test.  The PRNG seeds itself with a default value.
      // This will be the same each time, so if you want
      // nondeterministic behavior, you should pick the seed values
      // yourself.
      NormalGenerator< int, float > normgenS;
      datatype = "float";
      verifyLapackTemplate (normgenS, datatype, nrows, ncols, 
			    human_readable, b_debug);
      // Fetch the pseudorandom seed from the previous test.
      normgenS.getSeed (iseed);
      NormalGenerator< int, double > normgenD (iseed);
      // Next test.
      datatype = "double";
      verifyLapackTemplate (normgenD, datatype, nrows, ncols, 
			    human_readable, b_debug);

      if (test_complex_arithmetic)
	{
	  normgenD.getSeed (iseed);
	  NormalGenerator< int, complex<float> > normgenC (iseed);
	  datatype = "complex<float>";
	  verifyLapackTemplate (normgenC, datatype, nrows, ncols, 
				human_readable, b_debug);
	  normgenC.getSeed (iseed);
	  NormalGenerator< int, complex<double> > normgenZ (iseed);
	  datatype = "complex<double>";
	  verifyLapackTemplate (normgenZ, datatype, nrows, ncols, 
				human_readable, b_debug);
	}
    }


    void
    benchmarkSeqTsqr (std::ostream& out,
		      const int numRows,
		      const int numCols,
		      const int numTrials,
		      const size_t cacheBlockSize,
		      const bool contiguousCacheBlocks,
		      const bool testComplex,
		      const bool printFieldNames,
		      const bool humanReadable)
    {
      typedef Teuchos::Time timer_type;
      const bool testReal = true;
      using std::string;

      // Only print field names (if at all) for the first data type tested.
      bool printedFieldNames = false;

      if (testReal)
	{
	  { // Scalar=float
	    typedef SeqTsqrBenchmarker< int, float, timer_type > benchmark_type;
	    string scalarTypeName ("float");
	    benchmark_type widget (scalarTypeName, out, humanReadable);
	    widget.benchmark (numTrials, numRows, numCols, cacheBlockSize, 
			      contiguousCacheBlocks, 
			      printFieldNames && ! printedFieldNames);
	    if (printFieldNames && ! printedFieldNames)
	      printedFieldNames = true;
	  }
	  { // Scalar=double
	    typedef SeqTsqrBenchmarker< int, double, timer_type > benchmark_type;
	    string scalarTypeName ("double");
	    benchmark_type widget (scalarTypeName, out, humanReadable);
	    widget.benchmark (numTrials, numRows, numCols, cacheBlockSize, 
			      contiguousCacheBlocks, 
			      printFieldNames && ! printedFieldNames);
	    if (printFieldNames && ! printedFieldNames)
	      printedFieldNames = true;
	  }
	}

      if (testComplex)
	{
#ifdef HAVE_TSQR_COMPLEX
	  using std::complex;
	  { // Scalar=complex<float>
	    typedef SeqTsqrBenchmarker< int, complex<float>, timer_type > benchmark_type;
	    string scalarTypeName ("complex<float>");
	    benchmark_type widget (scalarTypeName, out, humanReadable);
	    widget.benchmark (numTrials, numRows, numCols, cacheBlockSize, 
			      contiguousCacheBlocks, 
			      printFieldNames && ! printedFieldNames);
	    if (printFieldNames && ! printedFieldNames)
	      printedFieldNames = true;
	  }
	  { // Scalar=complex<double>
	    typedef SeqTsqrBenchmarker< int, complex<double>, timer_type > benchmark_type;
	    string scalarTypeName ("complex<double>");
	    benchmark_type widget (scalarTypeName, out, humanReadable);
	    widget.benchmark (numTrials, numRows, numCols, cacheBlockSize, 
			      contiguousCacheBlocks, 
			      printFieldNames && ! printedFieldNames);
	    if (printFieldNames && ! printedFieldNames)
	      printedFieldNames = true;
	  }
#else // Don't HAVE_TSQR_COMPLEX
	  throw std::logic_error("TSQR not built with complex arithmetic support");
#endif // HAVE_TSQR_COMPLEX
	}
    }



  } // namespace Test
} // namespace TSQR
