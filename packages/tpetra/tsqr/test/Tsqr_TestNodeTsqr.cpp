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
// ************************************************************************
//@HEADER

#include "Tsqr_ConfigDefs.hpp"
#include "Teuchos_ConfigDefs.hpp" // HAVE_MPI
#ifdef HAVE_MPI
#  include "Teuchos_GlobalMPISession.hpp"
#  include "Teuchos_oblackholestream.hpp"
#endif // HAVE_MPI
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Time.hpp"

#include "Tsqr_Impl_Lapack.hpp"
#include "Tsqr_Random_NormalGenerator.hpp"
#include "Tsqr_LocalVerify.hpp"
#include "Tsqr_Matrix.hpp"
#include "Tsqr_NodeTsqrFactory.hpp"
#include "Tsqr_nodeTestProblem.hpp"
#include "Tsqr_Util.hpp"

#include <algorithm>
#include <complex>
#include <cstring> // size_t definition
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace TSQR {
  namespace Test {

    using execution_space = Kokkos::DefaultExecutionSpace;
    using memory_space = execution_space::memory_space;
    using device_type =
      Kokkos::Device<execution_space, memory_space>;

    // Command-line arguments and other test parameters.
    struct NodeTestParameters {
      NodeTestParameters () = default;

      std::string nodeTsqrType {"Default"};
      bool verify = true;
      bool benchmark = false;
      int numRows = 10000;
      int numCols = 10;
      int numTrials = 10;
      bool testReal = true;
#ifdef HAVE_KOKKOSTSQR_COMPLEX
      bool testComplex = true;
#else
      bool testComplex = false;
#endif // HAVE_KOKKOSTSQR_COMPLEX
      size_t cacheSizeHint = 0;
      bool contiguousCacheBlocks = false;
      bool printFieldNames = true;
      bool printTrilinosTestStuff = true;
      bool humanReadable = false;
      bool verbose = false;
      bool saveMatrices = false;
    };

    void
    printNodeTestParameters (std::ostream& out,
                             const NodeTestParameters& p,
                             const std::string& prefix)
    {
      using std::endl;
      out << prefix << "NodeTsqr: " << p.nodeTsqrType << endl
          << prefix << "numRows: " << p.numRows << endl
          << prefix << "numCols: " << p.numCols << endl
          << prefix << "numTrials: " << p.numTrials << endl
          << prefix << "testReal: "
          << (p.testReal ? "true" : "false") << endl
          << prefix << "testComplex: "
          << (p.testComplex ? "true" : "false") << endl
          << prefix << "cacheSizeHint: " << p.cacheSizeHint << endl
          << prefix << "contiguousCacheBlocks: "
          << (p.contiguousCacheBlocks ? "true" : "false") << endl
          << prefix << "printFieldNames: "
          << (p.printFieldNames ? "true" : "false") << endl
          << prefix << "printTrilinosTestStuff: "
          << (p.printTrilinosTestStuff ? "true" : "false") << endl
          << prefix << "humanReadable: "
          << (p.humanReadable ? "true" : "false") << endl
          << prefix << "verbose: "
          << (p.verbose ? "true" : "false") << endl
          << prefix << "saveMatrices: "
          << (p.saveMatrices ? "true" : "false") << endl;
    }

    void
    setBoolCmdLineOpt (Teuchos::CommandLineProcessor& cmdLineProc,
                       bool* variable,
                       const char trueString[],
                       const char falseString[],
                       const char docString[])
    {
      cmdLineProc.setOption (trueString, falseString, variable, docString);
    }

    // \brief Parse command-line options for this test
    //
    // \param argc [in] As usual in C(++)
    // \param argv [in] As usual in C(++)
    // \param allowedToPrint [in] Whether this (MPI) process is allowed
    //   to print to stdout/stderr.  Different per (MPI) process.
    // \param printedHelp [out] Whether this (MPI) process printed the
    //   "help" display (summary of command-line options)
    //
    // \return Encapsulation of command-line options
    static NodeTestParameters
    parseOptions (int argc,
                  char* argv[],
                  const bool allowedToPrint,
                  bool& printedHelp)
    {
      using std::cerr;
      using std::endl;

      printedHelp = false;

      // Command-line parameters, set to their default values.
      NodeTestParameters params;
      /// We really want the cache block size as a size_t, but
      /// Teuchos::CommandLineProcessor doesn't offer that option.
      /// So we read it in as an int, which means negative inputs
      /// are possible.  We check for those below in the input
      /// validation phase.
      //
      // Fetch default value of cacheSizeHint.
      int cacheSizeHintAsInt = static_cast<int> (params.cacheSizeHint);
      try {
        using Teuchos::CommandLineProcessor;
        CommandLineProcessor cmdLineProc (/* throwExceptions=*/ true,
                                          /* recognizeAllOptions=*/ false);
        const char docString[] = "This program tests TSQR::NodeTsqr, "
          "which implements the intraprocess part of TSQR.  "
          "Accuracy and performance tests are included.";
        cmdLineProc.setDocString (docString);

        setBoolCmdLineOpt (cmdLineProc, &params.verify,
                           "verify",
                           "noverify",
                           "Test accuracy");
        setBoolCmdLineOpt (cmdLineProc, &params.benchmark,
                           "benchmark",
                           "nobenchmark",
                           "Test performance");
        cmdLineProc.setOption ("numRows",
                               &params.numRows,
                               "Number of rows in the test matrix");
        cmdLineProc.setOption ("numCols",
                               &params.numCols,
                               "Number of columns in the test matrix");
        cmdLineProc.setOption ("numTrials",
                               &params.numTrials,
                               "Number of trials (only used when "
                               "\"--benchmark\"");
        setBoolCmdLineOpt (cmdLineProc, &params.testReal,
                           "testReal",
                           "noTestReal",
                           "Test real arithmetic");
        setBoolCmdLineOpt (cmdLineProc, &params.testComplex,
                           "testComplex",
                           "noTestComplex",
                           "Test complex arithmetic");
        cmdLineProc.setOption ("cacheBlockSize",
                               &cacheSizeHintAsInt,
                               "Cache size hint in bytes (0 means pick a reasonable default)");
        setBoolCmdLineOpt (cmdLineProc,
                           &params.contiguousCacheBlocks,
                           "contiguousCacheBlocks",
                           "noncontiguousCacheBlocks",
                           "Whether cache blocks should be stored contiguously");
        setBoolCmdLineOpt (cmdLineProc, &params.printFieldNames,
                           "printFieldNames",
                           "noPrintFieldNames",
                           "Print field names (for machine-readable output only)");
        setBoolCmdLineOpt (cmdLineProc, &params.printTrilinosTestStuff,
                           "printTrilinosTestStuff",
                           "noPrintTrilinosTestStuff",
                           "Print output that makes the Trilinos test framework happy, but may make benchmark results' parsing scripts unhappy.");
        setBoolCmdLineOpt (cmdLineProc, &params.humanReadable,
                           "humanReadable",
                           "machineReadable",
                           "If set, make output easy to read by humans, but harder to parse.");
        setBoolCmdLineOpt (cmdLineProc, &params.verbose,
                           "verbose",
                           "quiet",
                           "Print verbose debugging information");
        setBoolCmdLineOpt (cmdLineProc, &params.saveMatrices,
                           "saveMatrices",
                           "noSaveMatrices",
                           "If set, dump matrices to files.");
        cmdLineProc.setOption ("NodeTsqr",
                               &params.nodeTsqrType,
                               "NodeTsqr subclass type");
        cmdLineProc.parse (argc, argv);
      }
      catch (Teuchos::CommandLineProcessor::UnrecognizedOption& e) {
        if (allowedToPrint) {
          cerr << "Unrecognized command-line option: " << e.what () << endl;
        }
        throw e;
      }
      catch (Teuchos::CommandLineProcessor::HelpPrinted& e) {
        printedHelp = true;
        return params; // Don't verify parameters in this case
      }

      // Validate command-line options.  We provide default values
      // for unset options, so we don't have to validate those.
      if (params.numRows <= 0) {
        throw std::invalid_argument ("Number of rows must be positive");
      }
      else if (params.numCols <= 0) {
        throw std::invalid_argument ("Number of columns must be positive");
      }
      else if (params.numRows < params.numCols) {
        throw std::invalid_argument ("Number of rows must be >= number of columns");
      }
      else if (params.benchmark && params.numTrials < 1) {
        throw std::invalid_argument ("\"--benchmark\" option requires numTrials >= 1");
      }
      else {
        if (cacheSizeHintAsInt < 0) {
          throw std::invalid_argument ("Cache size hint must be nonnegative");
        }
        else {
          params.cacheSizeHint = size_t (cacheSizeHintAsInt);
        }
      }
      return params;
    }

    template<class Scalar>
    static int
    lworkQueryLapackQr (Impl::Lapack<Scalar>& lapack,
                        const int nrows,
                        const int ncols,
                        const int lda)
    {
      using std::ostringstream;
      using std::endl;
      using STS = Teuchos::ScalarTraits<Scalar>;
      using mag_type = typename STS::magnitudeType;

      Scalar d_lwork_geqrf {};
      lapack.compute_QR (nrows, ncols, nullptr, lda, nullptr,
                         &d_lwork_geqrf, -1);

      Scalar d_lwork_orgqr {};
      // A workspace query appropriate for computing the explicit Q
      // factor (nrows x ncols) in place, from the QR factorization of
      // an nrows x ncols matrix with leading dimension lda.
      lapack.compute_explicit_Q (nrows, ncols, ncols, nullptr, lda,
                                 nullptr, &d_lwork_orgqr, -1);

      // LAPACK workspace queries do return their results as a
      // double-precision floating-point value, but LAPACK promises
      // that that value will fit in an int.  Thus, we don't need to
      // check for valid casts to int below.  I include the checks
      // just to be "bulletproof" and also to show how to do the
      // checks for later reference.
      const mag_type lwork_geqrf_test
        (int (STS::magnitude (d_lwork_geqrf)));
      if (lwork_geqrf_test != STS::magnitude (d_lwork_geqrf)) {
        ostringstream os;
        os << "LAPACK _GEQRF workspace query returned a result, "
           << d_lwork_geqrf << ", bigger than the max int value, "
           << std::numeric_limits<int>::max ();
        throw std::range_error (os.str ());
      }
      const Scalar lwork_orgqr_test =
        mag_type (int (STS::magnitude ((d_lwork_orgqr))));
      if (lwork_orgqr_test != STS::magnitude (d_lwork_orgqr)) {
        ostringstream os;
        os << "LAPACK _UNGQR workspace query returned a result, "
           << d_lwork_orgqr << ", bigger than the max int value, "
           << std::numeric_limits<int>::max();
        throw std::range_error (os.str());
      }
      return std::max (static_cast<int> (STS::magnitude (d_lwork_geqrf)),
                       static_cast<int> (STS::magnitude (d_lwork_orgqr)));
    }

    template<class SC>
    Teuchos::RCP<
      typename ::TSQR::NodeTsqrFactory<SC, int, device_type>::node_tsqr_type
    >
    getNodeTsqr (const NodeTestParameters& p)
    {
      using fct_type = ::TSQR::NodeTsqrFactory<SC, int, device_type>;
      auto nodeTsqr = fct_type::getNodeTsqr (p.nodeTsqrType);
      TEUCHOS_ASSERT( ! nodeTsqr.is_null () );
      auto nodeTsqrParams = Teuchos::parameterList ("NodeTsqr");
      nodeTsqrParams->set ("Cache Size Hint", p.cacheSizeHint);
      nodeTsqr->setParameterList (nodeTsqrParams);
      return nodeTsqr;
    }

    static void
    printVerifyFieldNames (std::ostream& out)
    {
      const char prefix[] = "%";
      out << prefix << "method"
          << ",scalarType"
          << ",numRows"
          << ",numCols"
          << ",cacheSizeHint"
          << ",contiguousCacheBlocks"
          << ",frobA"
          << ",absFrobResid"
          << ",absFrobOrthog";
      out << std::endl;
    }

    template<class Scalar>
    static std::string
    getFileSuffix (const std::string& method)
    {
      std::string shortScalarType;
      if (std::is_same<Scalar, float>::value) {
        shortScalarType = "S";
      }
      else if (std::is_same<Scalar, double>::value) {
        shortScalarType = "D";
      }
      else if (std::is_same<Scalar, std::complex<float>>::value) {
        shortScalarType = "C";
      }
      else if (std::is_same<Scalar, std::complex<double>>::value) {
        shortScalarType = "Z";
      }
      else {
        shortScalarType = "U"; // unknown
      }
      const std::string sep ("_");
      return sep + method + sep + shortScalarType + ".txt";
    }

    // Test the accuracy of a NodeTsqr implementation on an nrows by
    // ncols matrix (using the given cache block size (in bytes)),
    // and print the results to stdout.
    template<class Scalar>
    static bool
    verifyNodeTsqrTmpl (std::ostream& out,
                        std::vector<int>& iseed,
                        const NodeTestParameters& params)
    {
      using Teuchos::TypeNameTraits;
      using std::cerr;
      using std::endl;
      using STS = Teuchos::ScalarTraits<Scalar>;
      using mag_type = typename STS::magnitudeType;
      using STM = Teuchos::ScalarTraits<mag_type>;
      const bool verbose = params.verbose;

      const std::string scalarType = TypeNameTraits<Scalar>::name ();
      const std::string fileSuffix =
        getFileSuffix<Scalar> (params.nodeTsqrType);
      if (verbose) {
        cerr << "Test NodeTsqr with Scalar=" << scalarType << endl;
      }

      const int nrows = params.numRows;
      const int ncols = params.numCols;

      Matrix<int, Scalar> A (nrows, ncols);
      Matrix<int, Scalar> A_copy (nrows, ncols);
      Matrix<int, Scalar> Q (nrows, ncols);
      Matrix<int, Scalar> R (ncols, ncols);
      if (std::numeric_limits<Scalar>::has_quiet_NaN) {
        deep_copy (A, std::numeric_limits<Scalar>::quiet_NaN ());
        deep_copy (A_copy, std::numeric_limits<Scalar>::quiet_NaN ());
        deep_copy (Q, std::numeric_limits<Scalar>::quiet_NaN ());
        deep_copy (R, std::numeric_limits<Scalar>::quiet_NaN ());
      }
      const int lda = nrows;
      const int ldq = nrows;
      const int ldr = ncols;

      if (verbose) {
        cerr << "-- Create test problem" << endl;
      }
      {
        TSQR::Random::NormalGenerator<int, Scalar> gen (iseed);
        nodeTestProblem (gen, nrows, ncols, A.data (),
                         A.stride(1), true);
        gen.getSeed (iseed); // fetch seed for the next test
      }

      if (params.saveMatrices) {
        std::string filename = std::string ("A") + fileSuffix;
        if (verbose) {
          cerr << "-- Save A to \"" << filename << "\"" << endl;
        }
        std::ofstream fileOut (filename.c_str ());
        print_local_matrix (fileOut, nrows, ncols,
                            A.data (), A.stride (1));
        fileOut.close ();
      }

      auto nodeTsqrPtr = getNodeTsqr<Scalar> (params);
      auto& actor = *nodeTsqrPtr;

      if (! params.contiguousCacheBlocks) {
        if (verbose) {
          cerr << "-- Copy A into A_copy" << endl;
        }
        deep_copy (A_copy, A);
      }
      else {
        if (verbose) {
          cerr << "-- Copy A into A_copy via cache_block" << endl;
        }
        actor.cache_block (nrows, ncols, A_copy.data (),
                           A.data (), A.stride (1));
        if (verbose) {
          cerr << "-- Verify cache_block result" << endl;
        }

        Matrix<int, Scalar> A2 (nrows, ncols);
        if (std::numeric_limits<Scalar>::has_quiet_NaN) {
          deep_copy (A2, std::numeric_limits<Scalar>::quiet_NaN ());
        }
        actor.un_cache_block (nrows, ncols, A2.data (),
                              A2.stride (1), A_copy.data ());
        const bool matrices_equal = matrix_equal (A, A2);
        TEUCHOS_TEST_FOR_EXCEPTION
          (! matrices_equal, std::logic_error, "cache_block failed!");
      }

      if (verbose) {
        cerr << "-- Fill R with zeros" << endl;
      }
      // We need to fill R with zeros, since the factorization may not
      // overwrite the strict lower triangle of R.
      deep_copy (R, Scalar {});

      if (verbose) {
        cerr << "-- Call NodeTsqr::factor" << endl;
      }
      auto factorOutput =
        actor.factor (nrows, ncols, A_copy.data(), A_copy.stride(1),
                      R.data(), R.stride(1),
                      params.contiguousCacheBlocks);
      if (params.saveMatrices) {
        std::string filename = std::string ("R") + fileSuffix;
        if (verbose) {
          cerr << "-- Save R to \"" << filename << "\"" << endl;
        }
        std::ofstream fileOut (filename.c_str ());
        print_local_matrix (fileOut, ncols, ncols,
                            R.data (), R.stride (1));
        fileOut.close ();
      }

      if (verbose) {
        cerr << "-- Call NodeTsqr::explicit_Q" << endl;
      }
      actor.explicit_Q (nrows, ncols, A_copy.data (), lda,
                        *factorOutput, ncols, Q.data (), Q.stride (1),
                        params.contiguousCacheBlocks);

      // "Un"-cache-block the output, if contiguous cache blocks were
      // used.  This is only necessary because local_verify() doesn't
      // currently support contiguous cache blocks.
      if (params.contiguousCacheBlocks) {
        // Use A_copy as temporary storage for un-cache-blocking Q.
        if (verbose) {
          cerr << "-- Call NodeTsqr::un_cache_block" << endl;
        }
        actor.un_cache_block (nrows, ncols, A_copy.data (),
                              A_copy.stride (1), Q.data ());
        deep_copy (Q, A_copy);
      }

      if (params.saveMatrices) {
        std::string filename = std::string ("Q") + fileSuffix;
        if (verbose) {
          cerr << "-- Save Q to \"" << filename << "\"" << endl;
        }
        std::ofstream fileOut (filename.c_str());
        print_local_matrix (fileOut, nrows, ncols,
                            Q.data (), Q.stride (1));
        fileOut.close ();
      }

      if (verbose) {
        cerr << "-- Call local_verify to validate the factorization"
             << endl;
      }
      auto results = local_verify (nrows, ncols, A.data (), lda,
                                   Q.data (), ldq, R.data (), ldr);

      if (verbose) {
        cerr << "-- Compute accuracy bounds and check" << endl;
      }

      // Accuracy relates to the number of floating-point operations,
      // which in turn is a function of the matrix's dimensions.
      // Avoid overflow of the local Ordinal type, by casting first to
      // a floating-point type.
      const mag_type dimsProd = mag_type(nrows) * mag_type(ncols) *
        mag_type(ncols);
      const mag_type fudgeFactor (10.0);
      // Relative residual error is ||A-Q*R|| / ||A||, or just
      // ||A-Q*R|| if ||A|| == 0.  (The result had better be zero in
      // the latter case.)  Square root of the matrix dimensions is an
      // old heuristic from Wilkinson or perhaps even an earlier
      // source.  We include a "fudge factor" so that the test won't
      // fail unless there is a really good reason.
      const mag_type relResidBound = fudgeFactor *
        STM::squareroot (dimsProd) * STS::eps ();

      // Relative residual error; avoid division by zero.
      const mag_type relResidError = results[0] /
        (results[2] == STM::zero () ? STM::one () : results[2]);

      bool success = true;
      if (relResidError > relResidBound) {
        success = false;
        if (verbose) {
          const std::string relResStr
            (results[2] == STM::zero () ? " / ||A||_F" : "");
          cerr << "*** For NodeTsqr=" << params.nodeTsqrType
               << " with Scalar=" << scalarType << ": "
               << "Residual ||A - QR||_F" << relResStr
               << " = " << relResidError << " > bound "
               << relResidBound << "." << endl;
        }
      }

      // Orthogonality of the matrix should not depend on the matrix
      // dimensions, if we measure in the 2-norm.  However, we are
      // measuring in the Frobenius norm, so it's appropriate to
      // multiply eps by the number of entries in the matrix for which
      // we compute the Frobenius norm.  We include a "fudge factor"
      // for the same reason as mentioned above.
      const mag_type orthoBound = fudgeFactor *
        mag_type (ncols) * mag_type (ncols) * STS::eps ();

      const mag_type orthoError = results[1];
      if (orthoError > orthoBound) {
        success = false;
        if (verbose) {
          cerr << "*** For NodeTsqr=" << params.nodeTsqrType
               << " with Scalar=" << scalarType << ": "
               << "Orthogonality ||I - Q^* Q||_F = " << orthoError
               << " > bound " << orthoBound << "." << endl;
        }
      }

      if (params.humanReadable) {
        out << "NodeTsqr subclass: " << params.nodeTsqrType
            << endl
            << "  - Scalar type: " << scalarType << endl
            << "  - Matrix dimensions: " << nrows << " by " << ncols
            << endl
            << "  - Cache Size Hint: " << params.cacheSizeHint
            << endl
            << "  - Contiguous cache blocks: "
            << (params.contiguousCacheBlocks ? "true" : "false")
            << endl
            << "  - Input matrix norm $\\| A \\|_F$: " << results[2]
            << endl
            << "  - Residual $\\| A - QR \\|_F$: " << results[0]
            << endl
            << "  - Orthogonality $\\| I - Q^* Q \\|_F$: "
            << results[1] << endl
            << endl;
      }
      else {
        out << params.nodeTsqrType
            << "," << scalarType
            << "," << nrows
            << "," << ncols
            << "," << params.cacheSizeHint
            << ","
            << (params.contiguousCacheBlocks ? "true" : "false")
            << "," << results[2]
            << "," << results[0]
            << "," << results[1];
        out << endl;
      }
      return success;
    }

    bool
    verifyNodeTsqr (std::ostream& out,
                    const NodeTestParameters& p)
    {
      // Seed for the next pseudorandom number generator.  We do tests
      // one after another, using the seed from the previous test in
      // the current test, so that the pseudorandom streams used by
      // the tests are independent.
      std::vector<int> iseed {{0, 0, 0, 1}};

      bool success = true;
      if (p.testReal) {
        const bool ok_S = verifyNodeTsqrTmpl<float> (out, iseed, p);
        const bool ok_D = verifyNodeTsqrTmpl<double> (out, iseed, p);
        success = success && ok_S && ok_D;
      }
      if (p.testComplex) {
#ifdef HAVE_KOKKOSTSQR_COMPLEX
        const bool ok_C =
          verifyNodeTsqrTmpl<std::complex<float>> (out, iseed, p);
        const bool ok_Z =
          verifyNodeTsqrTmpl<std::complex<double>> (out, iseed, p);
        success = success && ok_C && ok_Z;
#else // HAVE_KOKKOSTSQR_COMPLEX
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::logic_error, "TSQR was not built with complex "
           "arithmetic support.");
#endif // HAVE_KOKKOSTSQR_COMPLEX
      }
      return success;
    }

    template<class Scalar>
    static void
    verifyLapackTmpl (std::ostream& out,
                      std::vector<int>& iseed,
                      const NodeTestParameters& params)
    {
      using Teuchos::TypeNameTraits;
      using std::cerr;
      using std::endl;
      using STS = Teuchos::ScalarTraits<Scalar>;
      using mag_type = typename STS::magnitudeType;
      const bool verbose = params.verbose;

      const std::string scalarType = TypeNameTraits<Scalar>::name ();
      const std::string fileSuffix =
        getFileSuffix<Scalar> ("Lapack");
      if (verbose) {
        cerr << "Test LAPACK with Scalar=" << scalarType << endl;
      }

      const int nrows = params.numRows;
      const int ncols = params.numCols;

      Matrix<int, Scalar> A (nrows, ncols);
      Matrix<int, Scalar> A_copy (nrows, ncols);
      Matrix<int, Scalar> Q (nrows, ncols);
      Matrix<int, Scalar> R (ncols, ncols);
      if (std::numeric_limits<Scalar>::has_quiet_NaN) {
        deep_copy (A, std::numeric_limits< Scalar>::quiet_NaN ());
        deep_copy (A_copy, std::numeric_limits<Scalar>::quiet_NaN ());
        deep_copy (Q, std::numeric_limits<Scalar>::quiet_NaN ());
        deep_copy (R, std::numeric_limits<Scalar>::quiet_NaN ());
      }
      const int lda = nrows;
      const int ldq = nrows;
      const int ldr = ncols;

      if (verbose) {
        cerr << "-- Create test problem" << endl;
      }
      {
        TSQR::Random::NormalGenerator<int, Scalar> gen (iseed);
        nodeTestProblem (gen, nrows, ncols, A.data (),
                         A.stride (1), true);
        gen.getSeed (iseed); // fetch seed for the next test
      }

      if (params.saveMatrices) {
        std::string filename = std::string ("A") + fileSuffix;
        if (verbose) {
          cerr << "-- Save A to \"" << filename << "\"" << endl;
        }
        std::ofstream fileOut (filename.c_str ());
        print_local_matrix (fileOut, nrows, ncols,
                            A.data (), A.stride (1));
        fileOut.close ();
      }

      if (verbose) {
        cerr << "-- Copy A into A_copy" << endl;
      }
      deep_copy (A_copy, A);

      if (verbose) {
        cerr << "-- Do LAPACK lwork query" << endl;
      }
      Impl::Lapack<Scalar> lapack;
      const int lwork =
        lworkQueryLapackQr (lapack, nrows, ncols, A_copy.stride (1));
      if (verbose) {
        cerr << "-- lwork=" << lwork << endl;
      }
      std::vector<Scalar> work (lwork);
      std::vector<Scalar> tau (ncols);

      if (verbose) {
        cerr << "-- Fill R with zeros" << endl;
      }
      // We need to fill R with zeros, since the factorization may not
      // overwrite the strict lower triangle of R.
      deep_copy (R, Scalar {});

      if (verbose) {
        cerr << "-- Call Lapack::compute_QR" << endl;
      }
      lapack.compute_QR (nrows, ncols, A_copy.data (),
                         A_copy.stride (1), tau.data (),
                         work.data(), lwork);
      if (verbose) {
        cerr << "-- Copy R out of in-place result" << endl;
      }
      copy_upper_triangle (ncols, ncols, R.data(), ldr,
                           A_copy.data(), lda);
      if (params.saveMatrices) {
        std::string filename = std::string ("R") + fileSuffix;
        if (verbose) {
          cerr << "-- Save R to \"" << filename << "\"" << endl;
        }
        std::ofstream fileOut (filename.c_str ());
        print_local_matrix (fileOut, ncols, ncols,
                            R.data (), R.stride (1));
        fileOut.close ();
      }

      // The explicit Q factor will be computed in place, so copy the
      // result of the factorization into Q.
      deep_copy (Q, A_copy);

      if (verbose) {
        cerr << "-- Call Lapack::compute_explicit_Q" << endl;
      }
      lapack.compute_explicit_Q (nrows, ncols, ncols, Q.data (), ldq,
                                 tau.data (), work.data (), lwork);

      if (params.saveMatrices) {
        std::string filename = std::string ("Q") + fileSuffix;
        if (verbose) {
          cerr << "-- Save Q to \"" << filename << "\"" << endl;
        }
        std::ofstream fileOut (filename.c_str());
        print_local_matrix (fileOut, nrows, ncols,
                            Q.data (), Q.stride (1));
        fileOut.close ();
      }

      if (verbose) {
        cerr << "-- Call local_verify to validate the factorization"
             << endl;
      }
      auto results = local_verify (nrows, ncols, A.data (), lda,
                                   Q.data (), ldq, R.data (), ldr);

      if (params.humanReadable) {
        out << "LAPACK QR:" << endl
            << "  - Scalar type: " << scalarType << endl
            << "  - Matrix dimensions: " << nrows << " by " << ncols
            << endl
            << "  - Matrix norm $\\| A \\|_F$: "
            << results[2] << endl
            << "  - Residual $\\| A - QR \\|_F$: "
            << results[0] << endl
            << "  - Orthogonality $\\| I - Q^* Q \\|_F$: "
            << results[1] << endl
            << endl;
      }
      else {
        out << "LAPACK"
            << "," << scalarType
            << "," << nrows
            << "," << ncols
            << ",0"     // cacheSizeHint
            << ",false" // contiguousCacheBlocks
            << "," << results[2]
            << "," << results[0]
            << "," << results[1];
        out << endl;
      }
    }

    void
    verifyLapack (std::ostream& out,
                  const NodeTestParameters& p)
    {
      // We do tests one after another, using the seed from the
      // previous test in the current test, so that the pseudorandom
      // streams used by the tests are independent.

      std::vector<int> iseed {{0, 0, 0, 1}};

      if (p.testReal) {
        verifyLapackTmpl<float> (out, iseed, p);
        verifyLapackTmpl<double> (out, iseed, p);
      }
      if (p.testComplex) {
#ifdef HAVE_KOKKOSTSQR_COMPLEX
        verifyLapackTmpl<std::complex<float>> (out, iseed, p);
        verifyLapackTmpl<std::complex<double>> (out, iseed, p);
#else // HAVE_KOKKOSTSQR_COMPLEX
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::logic_error, "TSQR was not built with complex "
           "arithmetic support.");
#endif // HAVE_KOKKOSTSQR_COMPLEX
      }
    }

    static void
    printBenchmarkFieldNames (std::ostream& out)
    {
      const char prefix[] = "%";
      out << prefix << "method"
          << ",scalarType"
          << ",numRows"
          << ",numCols"
          << ",cacheSizeHint"
          << ",contiguousCacheBlocks"
          << ",numTrials"
          << ",timing" << std::endl;
    }

    template<class Scalar>
    void
    benchmarkLapackTmpl (std::ostream& out,
                         std::vector<int>& iseed,
                         const NodeTestParameters& testParams)
    {
      using std::endl;

      const int numRows = testParams.numRows;
      const int numCols = testParams.numCols;
      const int numTrials = testParams.numTrials;

      Matrix<int, Scalar> A (numRows, numCols);
      Matrix<int, Scalar> Q (numRows, numCols);
      Matrix<int, Scalar> R (numCols, numCols);
      const int lda = numRows;
      const int ldq = numRows;
      const int ldr = numCols;

      {
        using prng_type = TSQR::Random::NormalGenerator<int, Scalar>;
        prng_type gen (iseed);
        nodeTestProblem (gen, numRows, numCols,
                         A.data (), lda, false);
        gen.getSeed (iseed);
      }

      // Copy A into Q, since LAPACK QR overwrites the input.  We only
      // need Q because LAPACK's computation of the explicit Q factor
      // occurs in place.  This doesn't work with TSQR.  To give
      // LAPACK QR the fullest possible advantage over TSQR, we don't
      // allocate an A_copy here (as we would when benchmarking TSQR).
      deep_copy (Q, A);

      // Determine the required workspace for the factorization
      Impl::Lapack<Scalar> lapack;
      const int lwork =
        lworkQueryLapackQr (lapack, numRows, numCols, lda);
      std::vector<Scalar> work (lwork);
      std::vector<Scalar> tau (numCols);

      // Benchmark LAPACK's QR factorization for numTrials trials.
      Teuchos::Time timer ("LAPACK");
      timer.start ();
      for (int trialNum = 0; trialNum < numTrials; ++trialNum) {
        lapack.compute_QR (numRows, numCols, Q.data (), ldq,
                           tau.data (), work.data (), lwork);
        // Extract the upper triangular factor R from Q (where it was
        // computed in place by GEQRF), since UNGQR will overwrite all
        // of Q with the explicit Q factor.
        copy_upper_triangle (numRows, numCols, R.data (), ldr,
                             Q.data (), ldq);
        lapack.compute_explicit_Q (numRows, numCols, numCols,
                                   Q.data (), ldq, tau.data (),
                                   work.data (), lwork);
      }
      const double lapackTiming = timer.stop ();

      const std::string scalarType =
        Teuchos::TypeNameTraits<Scalar>::name ();

      if (testParams.humanReadable) {
        out << "LAPACK\'s QR factorization (_GEQRF + _UNGQR):"
            << endl << "  Scalar type = " << scalarType << endl
            << "  # rows = " << numRows << endl
            << "  # columns = " << numCols << endl
            << "  # trials = " << numTrials << endl
            << "Total time (s) = " << lapackTiming << endl
            << endl;
      }
      else {
        // "0" refers to the cache size hint, which is not applicable
        // in this case; we retain it for easy comparison of results
        // with NodeTsqr (so that the number of fields is the same in
        // both cases).  "false" (that follows 0) refers to whether or
        // not contiguous cache blocks were used (see TSQR::NodeTsqr);
        // this is also not applicable here.
        out << "LAPACK"
            << "," << scalarType
            << "," << numRows
            << "," << numCols
            << ",0"
            << ",false"
            << "," << numTrials
            << "," << lapackTiming << endl;
      }
    }

    void
    benchmarkLapack (std::ostream& out,
                     const NodeTestParameters& p)
    {
      std::vector<int> iseed {{0, 0, 0, 1}};
      if (p.testReal) {
        benchmarkLapackTmpl<float> (out, iseed, p);
        benchmarkLapackTmpl<double> (out, iseed, p);
      }
      if (p.testComplex) {
#ifdef HAVE_KOKKOSTSQR_COMPLEX
        benchmarkLapackTmpl<std::complex<float>> (out, iseed, p);
        benchmarkLapackTmpl<std::complex<double>> (out, iseed, p);
#else // Don't HAVE_KOKKOSTSQR_COMPLEX
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::logic_error,
           "TSQR was not built with complex arithmetic support.");
#endif // HAVE_KOKKOSTSQR_COMPLEX
      }
    }

    template<class Scalar>
    void
    benchmarkNodeTsqrTmpl (std::ostream& out,
                           std::vector<int>& iseed,
                           const NodeTestParameters& testParams)
    {
      using std::endl;
      auto nodeTsqrPtr = getNodeTsqr<Scalar> (testParams);
      auto& actor = *nodeTsqrPtr;

      const int numRows = testParams.numRows;
      const int numCols = testParams.numCols;
      const int numTrials = testParams.numTrials;
      const bool contiguousCacheBlocks =
        testParams.contiguousCacheBlocks;

      Matrix<int, Scalar> A (numRows, numCols);
      Matrix<int, Scalar> A_copy (numRows, numCols);
      Matrix<int, Scalar> Q (numRows, numCols);
      Matrix<int, Scalar> R (numCols, numCols);
      const int lda = numRows;
      const int ldq = numRows;

      {
        using prng_type = TSQR::Random::NormalGenerator<int, Scalar>;
        prng_type gen (iseed);
        nodeTestProblem (gen, numRows, numCols,
                         A.data (), lda, false);
        gen.getSeed (iseed);
      }
      deep_copy (A_copy, A); // need copy since TSQR overwrites

      // Benchmark sequential TSQR for numTrials trials.
      Teuchos::Time timer ("NodeTsqr");
      timer.start();
      for (int trialNum = 0; trialNum < numTrials; ++trialNum) {
        // Factor the matrix and extract the resulting R factor
        auto factorOutput =
          actor.factor (numRows, numCols, A_copy.data(), lda,
                        R.data(), R.stride(1), contiguousCacheBlocks);
        // Compute the explicit Q factor.  Unlike with LAPACK, this
        // doesn't happen in place: the implicit Q factor is stored in
        // A_copy, and the explicit Q factor is written to Q.
        actor.explicit_Q (numRows, numCols, A_copy.data (), lda,
                          *factorOutput, numCols, Q.data (), ldq,
                          contiguousCacheBlocks);
      }
      const double nodeTsqrTiming = timer.stop ();

      const std::string scalarType =
        Teuchos::TypeNameTraits<Scalar>::name ();

      if (testParams.humanReadable) {
        out << "NodeTsqr:" << endl
            << "  Scalar type = " << scalarType << endl
            << "  # rows = " << numRows << endl
            << "  # columns = " << numCols << endl
            << "  cache size hint in bytes = "
            << testParams.cacheSizeHint << endl
            << "  contiguous cache blocks? "
            << (contiguousCacheBlocks ? "true" : "false") << endl
            << "  # trials = " << numTrials << endl
            << "Total time (s) = " << nodeTsqrTiming << endl
            << endl;
      }
      else {
        out << testParams.nodeTsqrType
            << "," << scalarType
            << "," << numRows
            << "," << numCols
            << "," << testParams.cacheSizeHint
            << "," << (contiguousCacheBlocks ? "true" : "false")
            << "," << numTrials
            << "," << nodeTsqrTiming << endl;
      }
    }

    void
    benchmarkNodeTsqr (std::ostream& out,
                       const NodeTestParameters& p)
    {
      using Teuchos::TypeNameTraits;
      using LO = int;

      std::vector<int> iseed {{0, 0, 0, 1}};
      if (p.testReal) {
        benchmarkNodeTsqrTmpl<float> (out, iseed, p);
        benchmarkNodeTsqrTmpl<double> (out, iseed, p);
      }
      if (p.testComplex) {
#ifdef HAVE_KOKKOSTSQR_COMPLEX
        benchmarkNodeTsqrTmpl<std::complex<float>> (out, iseed, p);
        benchmarkNodeTsqrTmpl<std::complex<double>> (out, iseed, p);
#else // Don't HAVE_KOKKOSTSQR_COMPLEX
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::logic_error,
           "TSQR was not built with complex arithmetic support.");
#endif // HAVE_KOKKOSTSQR_COMPLEX
      }
    }
  } // namespace Test
} // namespace TSQR

int
main (int argc, char *argv[])
{
  using TSQR::Test::parseOptions;
  using std::endl;

#ifdef HAVE_MPI
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackhole);
  auto comm = Teuchos::DefaultComm<int>::getComm ();
  const int myRank = comm->getRank();
  // Only Process 0 writes to stdout.  The other processes send their
  // output to something that looks like /dev/null.
  std::ostream& out = (myRank == 0) ? std::cout : blackhole;
  // Only Process 0 performs the tests.
  const bool performingTests = (myRank == 0);
  const bool mayPrint = (myRank == 0);
#else // Don't HAVE_MPI: single-process test
  const bool performingTests = true;
  const bool mayPrint = true;
  std::ostream& out = std::cout;
#endif // HAVE_MPI

  // Fetch command-line parameters.
  bool printedHelp = false;
  auto params = parseOptions (argc, argv, mayPrint, printedHelp);
  if (printedHelp) {
    return EXIT_SUCCESS;
  }

  if (mayPrint) {
    out << "NodeTsqr verify/benchmark test options:" << endl;
    printNodeTestParameters (out, params, "  - ");
  }

  bool success = true;
  try {
    if (performingTests) {
      // We allow the same run to do both benchmark and verify.
      if (params.verify) {
        if (mayPrint && ! params.humanReadable) {
          TSQR::Test::printVerifyFieldNames (out);
        }
        TSQR::Test::verifyLapack (out, params);
        success = TSQR::Test::verifyNodeTsqr (out, params);
      }
      if (params.benchmark) {
        if (mayPrint && ! params.humanReadable) {
          TSQR::Test::printBenchmarkFieldNames (out);
        }
        TSQR::Test::benchmarkLapack (out, params);
        TSQR::Test::benchmarkNodeTsqr (out, params);
      }

      if (params.printTrilinosTestStuff) {
        // The Trilinos test framework expects a message like this.
        if (success) {
          out << "\nEnd Result: TEST PASSED" << endl;
        }
        else {
          out << "\nEnd Result: TEST FAILED" << endl;
        }
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
