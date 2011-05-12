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
#include <Teuchos_Time.hpp>

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

    template<class Ordinal, class Scalar>
    static Ordinal
    lworkQueryLapackQr (LAPACK<Ordinal, Scalar>& lapack,
			const Ordinal nrows,
			const Ordinal ncols,
			const Ordinal lda)
    {
      typedef typename ScalarTraits<Scalar>::magnitude_type magnitude_type;
      using std::ostringstream;
      using std::endl;

      Scalar d_lwork_geqrf = Scalar (0);
      int INFO = 0;
      lapack.GEQRF (nrows, ncols, NULL, lda, NULL, &d_lwork_geqrf, -1, &INFO);
      if (INFO != 0)
	{
	  ostringstream os;
	  os << "LAPACK _GEQRF workspace size query failed: INFO = " << INFO;
	  // It's a logic error and not a runtime error, because the
	  // LWORK query should only fail if the input parameters have
	  // invalid (e.g., out of range) values.
	  throw std::logic_error (os.str());
	}

      Scalar d_lwork_orgqr = Scalar (0);
      // A workspace query appropriate for computing the explicit Q
      // factor (nrows x ncols) in place, from the QR factorization of
      // an nrows x ncols matrix with leading dimension lda.
      lapack.ORGQR (nrows, ncols, ncols, NULL, lda, NULL, &d_lwork_orgqr, -1, &INFO);
      if (INFO != 0)
	{
	  ostringstream os;
	  os << "LAPACK _ORGQR workspace size query failed: INFO = " << INFO;
	  // It's a logic error and not a runtime error, because the
	  // LWORK query should only fail if the input parameters have
	  // invalid (e.g., out of range) values.
	  throw std::logic_error (os.str());
	}

      // LAPACK workspace queries do return their results as a
      // double-precision floating-point value, but LAPACK promises
      // that that value will fit in an int.  Thus, we don't need to
      // check for valid casts to int below.  I include the checks
      // just to be "bulletproof" and also to show how to do the
      // checks for later reference.
      const magnitude_type lwork_geqrf_test = 
	static_cast< magnitude_type > (static_cast<Ordinal> (ScalarTraits<Scalar>::abs (d_lwork_geqrf)));
      if (lwork_geqrf_test != ScalarTraits<Scalar>::abs (d_lwork_geqrf))
	{
	  ostringstream os;
	  os << "LAPACK _GEQRF workspace query returned a result, " 
	     << d_lwork_geqrf << ", bigger than the max Ordinal value, " 
	     << std::numeric_limits<Ordinal>::max();
	  throw std::range_error (os.str());
	}
      const Scalar lwork_orgqr_test = 
	static_cast< magnitude_type > (static_cast<Ordinal> (ScalarTraits<Scalar>::abs ((d_lwork_orgqr))));
      if (lwork_orgqr_test != ScalarTraits<Scalar>::abs (d_lwork_orgqr))
	{
	  ostringstream os;
	  os << "LAPACK _ORGQR workspace query returned a result, " 
	     << d_lwork_orgqr << ", bigger than the max Ordinal value, " 
	     << std::numeric_limits<Ordinal>::max();
	  throw std::range_error (os.str());
	}
      return std::max (static_cast<Ordinal> (ScalarTraits<Scalar>::abs (d_lwork_geqrf)),
		       static_cast<Ordinal> (ScalarTraits<Scalar>::abs (d_lwork_orgqr)));
    }

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
			   const size_t cache_size_hint,
			   const bool contiguous_cache_blocks,
			   const bool save_matrices,
			   const std::string& additionalFieldNames,
			   const std::string& additionalData,
			   const bool printFieldNames,
			   const bool human_readable,
			   const bool b_debug)
    {
      typedef typename ScalarTraits<Scalar>::magnitude_type magnitude_type;
      using std::cerr;
      using std::endl;
      using std::pair;
      using std::string;
      using std::vector;

      SequentialTsqr< Ordinal, Scalar > actor (cache_size_hint);
      Ordinal numCacheBlocks;

      if (b_debug)
	{
	  cerr << "Sequential TSQR test problem:" << endl
	       << "* " << nrows << " x " << ncols << endl
	       << "* Cache size hint of " << actor.cache_size_hint() << " bytes" << endl;
	  if (contiguous_cache_blocks)
	    cerr << "* Contiguous cache blocks" << endl;
	}

      Matrix< Ordinal, Scalar > A (nrows, ncols);
      Matrix< Ordinal, Scalar > A_copy (nrows, ncols);
      Matrix< Ordinal, Scalar > Q (nrows, ncols);
      Matrix< Ordinal, Scalar > R (ncols, ncols);
      if (std::numeric_limits<Scalar>::has_quiet_NaN)
	{
	  A.fill (std::numeric_limits< Scalar>::quiet_NaN());
	  A_copy.fill (std::numeric_limits<Scalar>::quiet_NaN());
	  Q.fill (std::numeric_limits<Scalar>::quiet_NaN());
	  R.fill (std::numeric_limits<Scalar>::quiet_NaN());
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
	      if (std::numeric_limits<Scalar>::has_quiet_NaN)
		A2.fill (std::numeric_limits<Scalar>::quiet_NaN());

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
	    << "Scalar type: " << datatype << endl
	    << "Matrix dimensions: " << nrows << " by " << ncols << endl
	    << "Cache size hint in bytes: " << actor.cache_size_hint() << endl
	    << "Number of cache blocks: " << numCacheBlocks << endl
	    << "Contiguous cache blocks? " << contiguous_cache_blocks << endl
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
		  << ",cacheSizeHint"
		  << ",contiguousCacheBlocks"
		  << ",absFrobResid"
		  << ",absFrobOrthog"
		  << ",frobA";
	      if (! additionalFieldNames.empty())
		out << "," << additionalFieldNames;
	      out << endl;
	    }
	  out << "SeqTSQR"
	      << "," << datatype
	      << "," << nrows
	      << "," << ncols
	      << "," << actor.cache_size_hint()
	      << "," << contiguous_cache_blocks 
	      << "," << results[0]
	      << "," << results[1]
	      << "," << results[2];
	  if (! additionalData.empty())
	    out << "," << additionalData;
	  out << endl;
	}
    }


    void
    verifySeqTsqr (std::ostream& out,
		   const int nrows, 
		   const int ncols, 
		   const size_t cache_size_hint,
		   const bool test_complex_arithmetic,
		   const bool save_matrices,
		   const bool contiguous_cache_blocks,
		   const std::string& additionalFieldNames,
		   const std::string& additionalData,
		   const bool printFieldNames,
		   const bool human_readable,
		   const bool b_debug)
    {
      using TSQR::Random::NormalGenerator;
#ifdef HAVE_TSQR_COMPLEX
      using std::complex;
#endif // HAVE_TSQR_COMPLEX
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
			     cache_size_hint, contiguous_cache_blocks, 
			     save_matrices, additionalFieldNames, additionalData,
			     printFieldNames, human_readable, b_debug);
      // Fetch the pseudorandom seed from the previous test.
      normgenS.getSeed (iseed);
      NormalGenerator< int, double > normgenD (iseed);
      // Next test.
      datatype = "double";
      shortDatatype = "D";
      verifySeqTsqrTemplate (out, normgenD, datatype, shortDatatype, nrows, ncols, 
			     cache_size_hint, contiguous_cache_blocks, 
			     save_matrices, additionalFieldNames, additionalData,
			     printFieldNames, human_readable, b_debug);
#ifdef HAVE_TSQR_COMPLEX
      if (test_complex_arithmetic)
	{
	  normgenD.getSeed (iseed);
	  NormalGenerator< int, complex<float> > normgenC (iseed);
	  datatype = "complex<float>";
	  shortDatatype = "C";
	  verifySeqTsqrTemplate (out, normgenC, datatype, shortDatatype, nrows, ncols, 
				 cache_size_hint, contiguous_cache_blocks, 
				 save_matrices, additionalFieldNames, additionalData,
				 printFieldNames, human_readable, b_debug);
	  normgenC.getSeed (iseed);
	  NormalGenerator< int, complex<double> > normgenZ (iseed);
	  datatype = "complex<double>";
	  shortDatatype = "Z";
	  verifySeqTsqrTemplate (out, normgenZ, datatype, shortDatatype, nrows, ncols, 
				 cache_size_hint, contiguous_cache_blocks, 
				 save_matrices, additionalFieldNames, additionalData,
				 printFieldNames, human_readable, b_debug);
	}
#else // HAVE_TSQR_COMPLEX
      if (test_complex_arithmetic)
	throw std::logic_error ("Trilinos was not built with "
				"complex arithmetic support");
#endif // HAVE_TSQR_COMPLEX
    }



    template< class Ordinal, class Scalar >
    static void
    verifyLapackTemplate (std::ostream& out,
			  TSQR::Random::NormalGenerator< Ordinal, Scalar >& generator,
			  const std::string& datatype,
			  const Ordinal nrows, 
			  const Ordinal ncols, 
			  const std::string& additionalFieldNames,
			  const std::string& additionalData,
			  const bool printFieldNames,
			  const bool human_readable,
			  const bool b_debug)
    {
      typedef typename ScalarTraits<Scalar>::magnitude_type magnitude_type;
      using std::ostringstream;
      using std::cerr;
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
      if (std::numeric_limits<Scalar>::has_quiet_NaN)
	{
	  A.fill (std::numeric_limits< Scalar>::quiet_NaN());
	  A_copy.fill (std::numeric_limits<Scalar>::quiet_NaN());
	  Q.fill (std::numeric_limits<Scalar>::quiet_NaN());
	  R.fill (std::numeric_limits<Scalar>::quiet_NaN());
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
      std::vector<Scalar> work (lwork);
      std::vector<Scalar> tau (ncols);

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
	out << "LAPACK QR (DGEQRF and DORGQR):" << endl
	    << "Scalar type: " << datatype << endl
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
		  << ",cacheSizeHint"
		  << ",contiguousCacheBlocks"
		  << ",absFrobResid"
		  << ",absFrobOrthog"
		  << ",frobA";
	      if (! additionalFieldNames.empty())
		out << "," << additionalFieldNames;
	      out << endl;
	    }
	  out << "LAPACK"
	      << "," << datatype
	      << "," << nrows
	      << "," << ncols
	      << "," << size_t(0) // cache_size_hint
	      << "," << false     // contiguous_cache_blocks 
	      << "," << results[0]
	      << "," << results[1]
	      << "," << results[2];
	  if (! additionalData.empty())
	    out << "," << additionalData;
	  out << endl;
	}
    }


    void
    verifyLapack (std::ostream& out,
		  const int nrows, 
		  const int ncols, 
		  const bool test_complex_arithmetic,
		  const std::string& additionalFieldNames,
		  const std::string& additionalData,
		  const bool printFieldNames,
		  const bool human_readable,
		  const bool b_debug)
    {
      using TSQR::Random::NormalGenerator;
#ifdef HAVE_TSQR_COMPLEX
      using std::complex;
#endif // HAVE_TSQR_COMPLEX
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
      verifyLapackTemplate (out, normgenS, datatype, nrows, ncols, 
			    additionalFieldNames, additionalData,
			    printFieldNames, human_readable, b_debug);
      // Fetch the pseudorandom seed from the previous test.
      normgenS.getSeed (iseed);
      NormalGenerator< int, double > normgenD (iseed);
      // Next test.
      datatype = "double";
      verifyLapackTemplate (out, normgenD, datatype, nrows, ncols, 
			    additionalFieldNames, additionalData,
			    false, human_readable, b_debug);
#ifdef HAVE_TSQR_COMPLEX
      if (test_complex_arithmetic)
	{
	  normgenD.getSeed (iseed);
	  NormalGenerator< int, complex<float> > normgenC (iseed);
	  datatype = "complex<float>";
	  verifyLapackTemplate (out, normgenC, datatype, nrows, ncols,
				additionalFieldNames, additionalData, 
				false, human_readable, b_debug);
	  normgenC.getSeed (iseed);
	  NormalGenerator< int, complex<double> > normgenZ (iseed);
	  datatype = "complex<double>";
	  verifyLapackTemplate (out, normgenZ, datatype, nrows, ncols, 
				additionalFieldNames, additionalData,
				false, human_readable, b_debug);
	}
#else // HAVE_TSQR_COMPLEX
      if (test_complex_arithmetic)
	throw std::logic_error ("Trilinos was not built with "
				"complex arithmetic support");
#endif // HAVE_TSQR_COMPLEX
    }

    /// \class LapackBenchmarker
    /// \brief Template version of LAPACK QR benchmark
    ///
    /// LAPACK QR benchmark, templated on Ordinal, Scalar, and
    /// TimerType.
    template< class Ordinal, class Scalar, class TimerType >
    class LapackBenchmarker {
    public:
      typedef Ordinal ordinal_type;
      typedef Scalar scalar_type;

      /// \brief Constructor
      ///
      /// \param scalarTypeName [in] Human-readable name of the Scalar
      ///   type.
      /// \param out [out] Reference to the output stream (e.g.,
      ///   std::cout) to which to write benchmark results.
      /// \param humanReadable [in] Whether to print results to out in
      ///   a verbose human-readable way, or in a way that is easy to
      ///   parse with a script.  In either case, the results will be
      ///   printed in ASCII format.
      LapackBenchmarker (const std::string& scalarTypeName,
			 std::ostream& out = std::cout,
			 const bool humanReadable = false) :
	scalarTypeName_ (scalarTypeName),
	out_ (out), 
	humanReadable_ (humanReadable)
      {
	TSQR::Test::verifyTimerConcept< TimerType >();
      }

      void 
      benchmark (const int numTrials,
		 const Ordinal numRows,
		 const Ordinal numCols,
		 const std::string& additionalFieldNames, 
		 const std::string& additionalData,
		 const bool printFieldNames)
      {
	Matrix< Ordinal, Scalar > A (numRows, numCols);
	Matrix< Ordinal, Scalar > Q (numRows, numCols);
	Matrix< Ordinal, Scalar > R (numCols, numCols);
	const Ordinal lda = numRows;
	const Ordinal ldq = numRows;
	const Ordinal ldr = numCols;

	// Create a test problem
	nodeTestProblem (gen_, numRows, numCols, A.get(), lda, false);

	// Copy A into Q, since LAPACK QR overwrites the input.  We only
	// need Q because LAPACK's computation of the explicit Q factor
	// occurs in place.  This doesn't work with TSQR.  To give
	// LAPACK QR the fullest possible advantage over TSQR, we don't
	// allocate an A_copy here (as we would when benchmarking TSQR).
	Q.copy (A);

	// Determine the required workspace for the factorization
	const Ordinal lwork = lworkQueryLapackQr (lapack_, numRows, numCols, lda);
	std::vector<Scalar> work (lwork);
	std::vector<Scalar> tau (numCols);

	// Benchmark LAPACK's QR factorization for numTrials trials.
	//
	// Name of timer doesn't matter here; we only need the timing.
	TimerType timer("LAPACK");
	timer.start();
	for (int trialNum = 0; trialNum < numTrials; ++trialNum)
	  {
	    // Compute the QR factorization
	    int info = 0; // INFO is always an int
	    lapack_.GEQRF (numRows, numCols, Q.get(), ldq, &tau[0], &work[0], lwork, &info);
	    if (info != 0)
	      {
		std::ostringstream os;
		os << "LAPACK QR factorization (_GEQRF) failed: INFO = " << info;
		throw std::runtime_error (os.str());
	      }

	    // Extract the upper triangular factor R from Q (where it
	    // was computed in place by GEQRF), since ORGQR will
	    // overwrite all of Q with the explicit Q factor.
	    copy_upper_triangle (numRows, numCols, R.get(), ldr, Q.get(), ldq);

	    // Compute the explicit Q factor
	    lapack_.ORGQR (numRows, numCols, numCols, Q.get(), ldq,
			  &tau[0], &work[0], lwork, &info);
	    if (info != 0)
	      {
		std::ostringstream os;
		os << "LAPACK explicit Q computation (_ORGQR) failed: INFO = " << info;
		throw std::runtime_error (os.str());
	      }
	  }
	const double lapackTiming = timer.stop();
	reportResults (numTrials, numRows, numCols, lapackTiming, 
		       additionalFieldNames, additionalData, printFieldNames);
      }


    private:
      //! Wrapper around LAPACK routines.
      TSQR::LAPACK< Ordinal, Scalar > lapack_;
      
      /// \brief Pseudorandom normal(0,1) generator.  
      ///
      /// Default seed is OK, because this is a benchmark, not an
      /// accuracy test.
      TSQR::Random::NormalGenerator< ordinal_type, scalar_type > gen_;
      
      //! Human-readable string representation of the Scalar type.
      std::string scalarTypeName_;

      //! Output stream to which to print benchmark results.
      std::ostream& out_;

      /// \brief Whether results should be printed in a human-readable way,
      /// 
      /// rather than a way easily parsed by a script.
      bool humanReadable_;

      /// \brief Report benchmark results to out_
      void 
      reportResults (const int numTrials,
		     const Ordinal numRows,
		     const Ordinal numCols,
		     const double lapackTiming,
		     const std::string& additionalFieldNames, 
		     const std::string& additionalData,
		     const bool printFieldNames)
      {
	using std::endl;
	if (humanReadable_)
	  out_ << "LAPACK\'s QR factorization (_GEQRF + _ORGQR):" << endl
	       << "Scalar type = " << scalarTypeName_ << endl
	       << "# rows = " << numRows << endl
	       << "# columns = " << numCols << endl
	       << "# trials = " << numTrials << endl
	       << "Total time (s) = " << lapackTiming << endl 
	       << endl;
	else
	  {
	    if (printFieldNames)
	      {
		const char prefix[] = "%";
		out_ << prefix 
		     << "method" 
		     << ",scalarType"
		     << ",numRows"
		     << ",numCols"
		     << ",cacheSizeHint"
		     << ",contiguousCacheBlocks"
		     << ",numTrials"
		     << ",timing";
		if (! additionalFieldNames.empty())
		  out_ << "," << additionalFieldNames;
		out_ << endl;
	      }
	    // "0" refers to the cache size hint, which is not
	    // applicable in this case; we retain it for easy
	    // comparison of results with SequentialTsqr (so that the
	    // number of fields is the same in both cases).  "false"
	    // (that follows 0) refers to whether or not contiguous
	    // cache blocks were used (see TSQR::SequentialTsqr); this
	    // is also not applicable in this case.
	    out_ << "LAPACK" 
		 << "," << scalarTypeName_
		 << "," << numRows
		 << "," << numCols
		 << "," << 0
		 << "," << false
		 << "," << numTrials 
		 << "," << lapackTiming;
	    if (! additionalData.empty())
	      out_ << "," << additionalData;
	    out_ << endl;
	  }
      }
    };


    void
    benchmarkLapack (std::ostream& out,
		     const int numRows,
		     const int numCols,
		     const int numTrials,
		     const bool testComplex,
		     const std::string& additionalFieldNames,
		     const std::string& additionalData,
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
	    typedef LapackBenchmarker< int, float, timer_type > benchmark_type;
	    string scalarTypeName ("float");
	    benchmark_type widget (scalarTypeName, out, humanReadable);
	    widget.benchmark (numTrials, numRows, numCols, 
			      additionalFieldNames, additionalData,
			      printFieldNames && ! printedFieldNames);
	    if (printFieldNames && ! printedFieldNames)
	      printedFieldNames = true;
	  }
	  { // Scalar=double
	    typedef LapackBenchmarker< int, double, timer_type > benchmark_type;
	    string scalarTypeName ("double");
	    benchmark_type widget (scalarTypeName, out, humanReadable);
	    widget.benchmark (numTrials, numRows, numCols,
			      additionalFieldNames, additionalData,
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
	    typedef LapackBenchmarker< int, complex<float>, timer_type > benchmark_type;
	    string scalarTypeName ("complex<float>");
	    benchmark_type widget (scalarTypeName, out, humanReadable);
	    widget.benchmark (numTrials, numRows, numCols,
			      additionalFieldNames, additionalData,
			      printFieldNames && ! printedFieldNames);
	    if (printFieldNames && ! printedFieldNames)
	      printedFieldNames = true;
	  }
	  { // Scalar=complex<double>
	    typedef LapackBenchmarker<int, complex<double>, timer_type> benchmark_type;
	    string scalarTypeName ("complex<double>");
	    benchmark_type widget (scalarTypeName, out, humanReadable);
	    widget.benchmark (numTrials, numRows, numCols,
			      additionalFieldNames, additionalData,
			      printFieldNames && ! printedFieldNames);
	    if (printFieldNames && ! printedFieldNames)
	      printedFieldNames = true;
	  }
#else // Don't HAVE_TSQR_COMPLEX
	  throw std::logic_error ("Trilinos was not built with "
				  "complex arithmetic support");
#endif // HAVE_TSQR_COMPLEX
	}
    }



    /// \class SeqTsqrBenchmarker
    /// \brief Template version of SequentialTsqr benchmark.
    ///
    /// SequentialTsqr benchmark, templated on Ordinal, Scalar, and
    /// TimerType.
    template<class Ordinal, class Scalar, class TimerType>
    class SeqTsqrBenchmarker {
    public:
      typedef Ordinal ordinal_type;
      typedef Scalar scalar_type;

      /// \brief Constructor
      ///
      /// \param scalarTypeName [in] Human-readable name of the Scalar
      ///   type.
      /// \param out [out] Reference to the output stream (e.g.,
      ///   std::cout) to which to write benchmark results.
      /// \param humanReadable [in] Whether to print results to out in
      ///   a verbose human-readable way, or in a way that is easy to
      ///   parse with a script.  In either case, the results will be
      ///   printed in ASCII format.
      SeqTsqrBenchmarker (const std::string& scalarTypeName,
			  std::ostream& out = std::cout,
			  const bool humanReadable = false) : 
	scalarTypeName_ (scalarTypeName),
	out_ (out), 
	humanReadable_ (humanReadable)
      {
	// Make sure that TimerType satisfies the required interface.
	TSQR::Test::verifyTimerConcept<TimerType>();
      }

      void 
      benchmark (const int numTrials,
		 const Ordinal numRows,
		 const Ordinal numCols,
		 const size_t cacheSizeHint,
		 const bool contiguousCacheBlocks,
		 const std::string& additionalFieldNames,
		 const std::string& additionalData,
		 const bool printFieldNames)
      {
	SequentialTsqr<Ordinal, Scalar> actor (cacheSizeHint);

	Matrix<Ordinal, Scalar> A (numRows, numCols);
	Matrix<Ordinal, Scalar> A_copy (numRows, numCols);
	Matrix<Ordinal, Scalar> Q (numRows, numCols);
	Matrix<Ordinal, Scalar> R (numCols, numCols);
	const Ordinal lda = numRows;
	const Ordinal ldq = numRows;

	// Create a test problem
	nodeTestProblem (gen_, numRows, numCols, A.get(), lda, false);

	// Copy A into A_copy, since TSQR overwrites the input
	A_copy.copy (A);

	// Benchmark sequential TSQR for numTrials trials.
	//
	// Name of timer doesn't matter here; we only need the timing.
	TimerType timer("SeqTSQR");
	timer.start();
	for (int trialNum = 0; trialNum < numTrials; ++trialNum)
	  {
	    // Factor the matrix and extract the resulting R factor
	    typedef typename SequentialTsqr<Ordinal, Scalar>::FactorOutput 
	      factor_output_type;
	    factor_output_type factorOutput = 
	      actor.factor (numRows, numCols, A_copy.get(), lda, 
			    R.get(), R.lda(), contiguousCacheBlocks);
	    // Compute the explicit Q factor.  Unlike with LAPACK QR,
	    // this doesn't happen in place: the implicit Q factor is
	    // stored in A_copy, and the explicit Q factor is written to
	    // Q.
	    actor.explicit_Q (numRows, numCols, A_copy.get(), lda, factorOutput, 
			      numCols, Q.get(), ldq, contiguousCacheBlocks);
	  }
	const double seqTsqrTiming = timer.stop();
	reportResults (numTrials, numRows, numCols, actor.cache_size_hint(),
		       contiguousCacheBlocks, seqTsqrTiming, 
		       additionalFieldNames, additionalData, printFieldNames);
      }


    private:
      /// \brief Pseudorandom normal(0,1) generator.  
      ///
      /// Default seed is OK, because this is a benchmark, not an
      /// accuracy test.
      TSQR::Random::NormalGenerator<ordinal_type, scalar_type> gen_;
      
      //! Human-readable string representation of the Scalar type.
      std::string scalarTypeName_;

      //! Output stream to which to print benchmark results.
      std::ostream& out_;

      /// \brief Whether results should be printed in a human-readable way,
      ///
      /// as opposed to a way easily parsed by a script.
      bool humanReadable_;

      //! Report benchmark results to out_
      void 
      reportResults (const int numTrials,
		     const Ordinal numRows,
		     const Ordinal numCols,
		     const size_t actualCacheSizeHint,
		     const bool contiguousCacheBlocks,
		     const double seqTsqrTiming,
		     const std::string& additionalFieldNames,
		     const std::string& additionalData,
		     const bool printFieldNames)
      {
	using std::endl;
	if (humanReadable_)
	  out_ << "Sequential (cache-blocked) TSQR:" << endl
	       << "Scalar type = " << scalarTypeName_ << endl
	       << "# rows = " << numRows << endl
	       << "# columns = " << numCols << endl
	       << "cache size hint in bytes = " << actualCacheSizeHint << endl
	       << "contiguous cache blocks? " << contiguousCacheBlocks << endl
	       << "# trials = " << numTrials << endl
	       << "Total time (s) = " << seqTsqrTiming << endl 
	       << endl;
	else
	  {
	    if (printFieldNames)
	      {
		const char prefix[] = "%";
		out_ << prefix 
		     << "method" 
		     << ",scalarType"
		     << ",numRows"
		     << ",numCols"
		     << ",cacheSizeHint"
		     << ",contiguousCacheBlocks"
		     << ",numTrials"
		     << ",timing";
		if (! additionalFieldNames.empty())
		  out_ << "," << additionalFieldNames;
		out_ << endl;
	      }
	    out_ << "SeqTSQR" 
		 << "," << scalarTypeName_
		 << "," << numRows
		 << "," << numCols
		 << "," << actualCacheSizeHint
		 << "," << contiguousCacheBlocks
		 << "," << numTrials 
		 << "," << seqTsqrTiming;
	    if (! additionalData.empty())
	      out_ << "," << additionalData;
	    out_ << endl;
	  }
      }
    };


    void
    benchmarkSeqTsqr (std::ostream& out,
		      const int numRows,
		      const int numCols,
		      const int numTrials,
		      const size_t cacheSizeHint,
		      const bool contiguousCacheBlocks,
		      const bool testComplex,
		      const std::string& additionalFieldNames,
		      const std::string& additionalData,
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
	    typedef SeqTsqrBenchmarker<int, float, timer_type> benchmark_type;
	    string scalarTypeName ("float");
	    benchmark_type widget (scalarTypeName, out, humanReadable);
	    widget.benchmark (numTrials, numRows, numCols, cacheSizeHint, 
			      contiguousCacheBlocks, 
			      additionalFieldNames, additionalData,
			      printFieldNames && ! printedFieldNames);
	    if (printFieldNames && ! printedFieldNames)
	      printedFieldNames = true;
	  }
	  { // Scalar=double
	    typedef SeqTsqrBenchmarker< int, double, timer_type > benchmark_type;
	    string scalarTypeName ("double");
	    benchmark_type widget (scalarTypeName, out, humanReadable);
	    widget.benchmark (numTrials, numRows, numCols, cacheSizeHint, 
			      contiguousCacheBlocks, 
			      additionalFieldNames, additionalData,
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
	    widget.benchmark (numTrials, numRows, numCols, cacheSizeHint, 
			      contiguousCacheBlocks, 
			      additionalFieldNames, additionalData,
			      printFieldNames && ! printedFieldNames);
	    if (printFieldNames && ! printedFieldNames)
	      printedFieldNames = true;
	  }
	  { // Scalar=complex<double>
	    typedef SeqTsqrBenchmarker< int, complex<double>, timer_type > benchmark_type;
	    string scalarTypeName ("complex<double>");
	    benchmark_type widget (scalarTypeName, out, humanReadable);
	    widget.benchmark (numTrials, numRows, numCols, cacheSizeHint, 
			      contiguousCacheBlocks, 
			      additionalFieldNames, additionalData,
			      printFieldNames && ! printedFieldNames);
	    if (printFieldNames && ! printedFieldNames)
	      printedFieldNames = true;
	  }
#else // Don't HAVE_TSQR_COMPLEX
	  throw std::logic_error ("Trilinos was not built with "
				  "complex arithmetic support");
#endif // HAVE_TSQR_COMPLEX
	}
    }



  } // namespace Test
} // namespace TSQR
