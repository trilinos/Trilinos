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

#include "Teuchos_CommandLineProcessor.hpp"
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
      NodeTestParameters() = default;

      std::string nodeTsqrType {"Default"};
      bool verify = true;
      bool benchmark = false;
      int numRows = 10000;
      int numCols = 10;
      int numTrials = 10;
      bool testReal = true;
#ifdef HAVE_TPETRATSQR_COMPLEX
      bool testComplex = true;
#else
      bool testComplex = false;
#endif // HAVE_TPETRATSQR_COMPLEX
      size_t cacheSizeHint = 0;
      bool contiguousCacheBlocks = false;
      bool printFieldNames = true;
      bool printTrilinosTestStuff = true;
      bool humanReadable = false;
      bool verbose = false;
      bool saveMatrices = false;
    };

    void
    printNodeTestParameters(std::ostream& out,
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
    setBoolCmdLineOpt(Teuchos::CommandLineProcessor& cmdLineProc,
                      bool* variable,
                      const char trueString[],
                      const char falseString[],
                      const char docString[])
    {
      cmdLineProc.setOption(trueString, falseString, variable,
                            docString);
    }

    // \brief Parse command-line options for this test
    //
    // \param argc [in] As usual in C(++).
    // \param argv [in] As usual in C(++).
    // \param printedHelp [out] Whether this function printed the
    //   "help" display (summary of command-line options).
    //
    // \return Encapsulation of command-line options
    static NodeTestParameters
    parseOptions(int argc,
                 char* argv[],
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
      int cacheSizeHintAsInt = static_cast<int>(params.cacheSizeHint);
      try {
        const bool throwExceptions = true;
        const bool recognizeAllOptions = false;
        using Teuchos::CommandLineProcessor;
        CommandLineProcessor cmdLineProc(throwExceptions,
                                         recognizeAllOptions);
        const char docString[] = "This program tests TSQR::NodeTsqr, "
          "which implements the intraprocess part of TSQR.  "
          "Accuracy and performance tests are included.";
        cmdLineProc.setDocString(docString);

        setBoolCmdLineOpt(cmdLineProc, &params.verify,
                          "verify",
                          "noverify",
                          "Test accuracy");
        setBoolCmdLineOpt(cmdLineProc, &params.benchmark,
                          "benchmark",
                          "nobenchmark",
                          "Test performance");
        cmdLineProc.setOption("numRows",
                              &params.numRows,
                              "Number of rows in the test matrix");
        cmdLineProc.setOption("numCols",
                              &params.numCols,
                              "Number of columns in the test matrix");
        cmdLineProc.setOption("numTrials",
                              &params.numTrials,
                              "Number of trials (only used when "
                              "\"--benchmark\"");
        setBoolCmdLineOpt(cmdLineProc, &params.testReal,
                          "testReal",
                          "noTestReal",
                          "Test real arithmetic");
        setBoolCmdLineOpt(cmdLineProc, &params.testComplex,
                          "testComplex",
                          "noTestComplex",
                          "Test complex arithmetic");
        cmdLineProc.setOption("cacheBlockSize",
                              &cacheSizeHintAsInt,
                              "Cache size hint in bytes (0 means "
                              "pick a reasonable default)");
        setBoolCmdLineOpt(cmdLineProc,
                          &params.contiguousCacheBlocks,
                          "contiguousCacheBlocks",
                          "noncontiguousCacheBlocks",
                          "Whether cache blocks should be stored contiguously");
        setBoolCmdLineOpt(cmdLineProc, &params.printFieldNames,
                          "printFieldNames",
                          "noPrintFieldNames",
                          "Print field names (for machine-readable output only)");
        setBoolCmdLineOpt(cmdLineProc, &params.printTrilinosTestStuff,
                          "printTrilinosTestStuff",
                          "noPrintTrilinosTestStuff",
                          "Print output that makes the Trilinos test "
                          "framework happy, but may make benchmark "
                          "results' parsing scripts unhappy.");
        setBoolCmdLineOpt(cmdLineProc, &params.humanReadable,
                          "humanReadable",
                          "machineReadable",
                          "If set, make output easy to read by "
                          "humans, but harder to parse.");
        setBoolCmdLineOpt(cmdLineProc, &params.verbose,
                          "verbose",
                          "quiet",
                          "Print verbose debugging information");
        setBoolCmdLineOpt(cmdLineProc, &params.saveMatrices,
                          "saveMatrices",
                          "noSaveMatrices",
                          "If set, dump matrices to files.");
        cmdLineProc.setOption("NodeTsqr",
                              &params.nodeTsqrType,
                              "NodeTsqr subclass type");
        cmdLineProc.parse(argc, argv);
      }
      catch(Teuchos::CommandLineProcessor::UnrecognizedOption& e) {
        cerr << "Unrecognized command-line option: " << e.what()
             << endl;
        throw e;
      }
      catch(Teuchos::CommandLineProcessor::HelpPrinted& e) {
        printedHelp = true;
        return params; // Don't verify parameters in this case
      }

      // Validate command-line options.  We provide default values
      // for unset options, so we don't have to validate those.
      TEUCHOS_TEST_FOR_EXCEPTION
        (params.numRows <= 0, std::invalid_argument, "Number of "
         "rows must be positive, but you set --numRows=" <<
         params.numRows << ".");
      TEUCHOS_TEST_FOR_EXCEPTION
        (params.numCols <= 0, std::invalid_argument, "Number of "
         "columns must be positive, but you set --numCols=" <<
         params.numCols << ".");
      TEUCHOS_TEST_FOR_EXCEPTION
        (params.numRows < params.numCols, std::invalid_argument,
         "Number of rows must be >= number of columns, but you set "
         "--numRows=" << params.numRows << " and --numCols=" <<
         params.numCols << ".");
      TEUCHOS_TEST_FOR_EXCEPTION
        (params.benchmark && params.numTrials < 1,
         std::invalid_argument, "Since you set --benchmark, the "
         "number of trials must be positive, but you set --numTrials="
         << params.numTrials << ".");
      TEUCHOS_TEST_FOR_EXCEPTION
        (cacheSizeHintAsInt < 0, std::invalid_argument, "Cache size "
         "hint must be nonnegative, but you set --cacheBlockSize=" <<
         cacheSizeHintAsInt << ".");
      params.cacheSizeHint = size_t(cacheSizeHintAsInt);
      return params;
    }

    template<class Scalar>
    using kokkos_value_type = typename std::conditional<
        std::is_const<Scalar>::value,
        const typename Kokkos::ArithTraits<
          typename std::remove_const<Scalar>::type>::val_type,
        typename Kokkos::ArithTraits<Scalar>::val_type
      >::type;

    template<class LO, class Scalar>
    Kokkos::View<kokkos_value_type<Scalar>**,
                 Kokkos::LayoutLeft, Kokkos::HostSpace,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    getHostMatrixView(const MatView<LO, Scalar>& A)
    {
      using Kokkos::ALL;
      using Kokkos::subview;
      using IST = kokkos_value_type<Scalar>;
      using host_mat_view_type =
        Kokkos::View<IST**, Kokkos::LayoutLeft, Kokkos::HostSpace,
          Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

      const size_t nrows(A.extent(0));
      const size_t ncols(A.extent(1));
      const size_t lda(A.stride(1));
      IST* A_raw = reinterpret_cast<IST*>(A.data());
      host_mat_view_type A_full(A_raw, lda, ncols);
      const std::pair<size_t, size_t> rowRange(0, nrows);
      return Kokkos::subview(A_full, rowRange, Kokkos::ALL());
    }

    template<class LO, class Scalar>
    Kokkos::View<typename Kokkos::ArithTraits<Scalar>::val_type**,
                 Kokkos::LayoutLeft>
    getDeviceMatrixCopy(const MatView<LO, Scalar>& A,
                        const std::string& label)
    {
      using Kokkos::view_alloc;
      using Kokkos::WithoutInitializing;
      using IST = typename Kokkos::ArithTraits<Scalar>::val_type;
      using device_matrix_type =
        Kokkos::View<IST**, Kokkos::LayoutLeft>;

      const size_t nrows(A.extent(0));
      const size_t ncols(A.extent(1));
      device_matrix_type A_dev
        (view_alloc(label, WithoutInitializing), nrows, ncols);
      auto A_host = getHostMatrixView(A);
      Kokkos::deep_copy(A_dev, A_host);
      return A_dev;
    }


    template<template<class SC> class LapackType, class Scalar>
    static int
    lworkQueryLapackQr(LapackType<Scalar>& lapack,
                       const int nrows,
                       const int ncols,
                       const int lda)
    {
      const int lwork_geqrf =
        lapack.compute_QR_lwork(nrows, ncols, nullptr, lda);
      // A workspace query appropriate for computing the explicit Q
      // factor (nrows x ncols) in place, from the QR factorization of
      // an nrows x ncols matrix with leading dimension lda.
      const int lwork_ungqr =
        lapack.compute_explicit_Q_lwork(nrows, ncols, ncols,
                                        nullptr, lda, nullptr);
      return std::max(lwork_geqrf, lwork_ungqr);
    }

    template<class SC>
    Teuchos::RCP<
      typename ::TSQR::NodeTsqrFactory<SC, int, device_type>::node_tsqr_type
    >
    getNodeTsqr(const NodeTestParameters& p,
                const std::string& overrideNodeTsqrType = "")
    {
      const std::string nodeTsqrType = overrideNodeTsqrType == "" ?
        p.nodeTsqrType : overrideNodeTsqrType;

      TEUCHOS_ASSERT( nodeTsqrType != "" );
      using fct_type = ::TSQR::NodeTsqrFactory<SC, int, device_type>;
      auto nodeTsqr = fct_type::getNodeTsqr(nodeTsqrType);
      TEUCHOS_ASSERT( ! nodeTsqr.is_null() );
      auto nodeTsqrParams = Teuchos::parameterList("NodeTsqr");
      nodeTsqrParams->set("Cache Size Hint", p.cacheSizeHint);
      nodeTsqr->setParameterList(nodeTsqrParams);
      return nodeTsqr;
    }

    static void
    printVerifyFieldNames(std::ostream& out)
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
    getFileSuffix(const std::string& method)
    {
      std::string shortScalarType;
      if(std::is_same<Scalar, float>::value) {
        shortScalarType = "S";
      }
      else if(std::is_same<Scalar, double>::value) {
        shortScalarType = "D";
      }
      else if(std::is_same<Scalar, std::complex<float>>::value) {
        shortScalarType = "C";
      }
      else if(std::is_same<Scalar, std::complex<double>>::value) {
        shortScalarType = "Z";
      }
      else {
        shortScalarType = "U"; // unknown
      }
      const std::string sep("_");
      return sep + method + sep + shortScalarType + ".txt";
    }

    // Test the accuracy of a NodeTsqr implementation on an nrows by
    // ncols matrix (using the given cache block size (in bytes)),
    // and print the results to stdout.
    template<class Scalar>
    static bool
    verifyNodeTsqrTmpl(std::ostream& out,
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
      const std::string scalarType = TypeNameTraits<Scalar>::name();
      const std::string fileSuffix =
        getFileSuffix<Scalar>(params.nodeTsqrType);
      if(verbose) {
        cerr << "Test NodeTsqr with Scalar=" << scalarType << endl;
      }

      bool success = true;

      const int nrows = params.numRows;
      const int ncols = params.numCols;

      Matrix<int, Scalar> A(nrows, ncols);
      Matrix<int, Scalar> A_copy(nrows, ncols);
      Matrix<int, Scalar> Q(nrows, ncols);
      Matrix<int, Scalar> R(ncols, ncols);
      if(std::numeric_limits<Scalar>::has_quiet_NaN) {
        deep_copy(A, std::numeric_limits<Scalar>::quiet_NaN());
        deep_copy(A_copy, std::numeric_limits<Scalar>::quiet_NaN());
        deep_copy(Q, std::numeric_limits<Scalar>::quiet_NaN());
        deep_copy(R, std::numeric_limits<Scalar>::quiet_NaN());
      }
      const int lda = nrows;
      const int ldq = nrows;
      const int ldr = ncols;

      if(verbose) {
        cerr << "-- Create test problem" << endl;
      }
      {
        TSQR::Random::NormalGenerator<int, Scalar> gen(iseed);
        nodeTestProblem(gen, nrows, ncols, A.data(), A.stride(1),
                        true);
        gen.getSeed(iseed); // fetch seed for the next test
      }

      if(params.saveMatrices) {
        std::string filename = std::string("A") + fileSuffix;
        if(verbose) {
          cerr << "-- Save A to \"" << filename << "\"" << endl;
        }
        std::ofstream fileOut(filename.c_str());
        print_local_matrix(fileOut, nrows, ncols,
                           A.data(), A.stride(1));
        fileOut.close();
      }

      auto nodeTsqrPtr = getNodeTsqr<Scalar>(params);
      auto& actor = *nodeTsqrPtr;
      if(verbose && actor.wants_device_memory()) {
        cerr << "-- NodeTsqr claims to want device memory" << endl;
      }

      using IST = typename Kokkos::ArithTraits<Scalar>::val_type;
      using device_matrix_type =
        Kokkos::View<IST**, Kokkos::LayoutLeft>;

      auto A_h = getHostMatrixView(A.view());
      auto A_copy_h = getHostMatrixView(A_copy.view());
      auto Q_h = getHostMatrixView(Q.view());
      device_matrix_type A_d;
      device_matrix_type A_copy_d;
      device_matrix_type Q_d;
      if(actor.wants_device_memory()) {
        A_d = getDeviceMatrixCopy(A.view(), "A_d");
        // Don't copy A_copy yet; see below.
        A_copy_d = device_matrix_type("A_copy_d", nrows, ncols);
        Q_d = device_matrix_type("Q_d", nrows, ncols);
      }

      if(! params.contiguousCacheBlocks) {
        if(verbose) {
          cerr << "-- Copy A into A_copy" << endl;
        }
        deep_copy(A_copy, A);
        if(actor.wants_device_memory()) {
          deep_copy(A_copy_d, A_d);
        }
      }
      else {
        if(verbose) {
          cerr << "-- Copy A into A_copy via cache_block" << endl;
        }
        if(actor.wants_device_memory()) {
          Scalar* A_copy_d_raw =
            reinterpret_cast<Scalar*>(A_copy_d.data());
          const Scalar* A_d_raw =
            reinterpret_cast<const Scalar*>(A_d.data());
          actor.cache_block(nrows, ncols, A_copy_d_raw,
                            A_d_raw, A_d.stride(1));
          Kokkos::deep_copy(A_copy_h, A_copy_d);
        }
        else {
          actor.cache_block(nrows, ncols, A_copy.data(),
                            A.data(), A.stride(1));
        }
        if(verbose) {
          cerr << "-- Verify cache_block result" << endl;
        }

        Matrix<int, Scalar> A2(nrows, ncols);
        if(std::numeric_limits<Scalar>::has_quiet_NaN) {
          deep_copy(A2, std::numeric_limits<Scalar>::quiet_NaN());
        }
        if(actor.wants_device_memory()) {
          auto A2_h = getHostMatrixView(A2.view());
          auto A2_d = getDeviceMatrixCopy(A2.view(), "A2_d");
          Scalar* A2_d_raw = reinterpret_cast<Scalar*>(A2_d.data());
          const Scalar* A_copy_d_raw =
            reinterpret_cast<const Scalar*>(A_copy_d.data());
          actor.un_cache_block(nrows, ncols, A2_d_raw,
                               A2_d.stride(1), A_copy_d_raw);
          Kokkos::deep_copy(A2_h, A2_d);
        }
        else {
          actor.un_cache_block(nrows, ncols, A2.data(),
                               A2.stride(1), A_copy.data());
        }
        const bool matrices_equal = matrix_equal(A, A2);
        if(! matrices_equal) {
          success = false;
          if(verbose) {
            cerr << "*** cache_block failed!" << endl;
          }
        }
      }

      if(verbose) {
        cerr << "-- Fill R with zeros" << endl;
      }
      // We need to fill R with zeros, since the factorization may not
      // overwrite the strict lower triangle of R.
      deep_copy(R, Scalar {});

      if(verbose) {
        cerr << "-- Call NodeTsqr::factor" << endl;
      }
      // R is always in host memory, because that's what Belos wants.
      auto factorOutput = [&]() {
        if(actor.wants_device_memory()) {
          Scalar* A_copy_d_raw =
            reinterpret_cast<Scalar*>(A_copy_d.data());
          TEUCHOS_ASSERT( nrows == 0 || ncols == 0 ||
                          A_copy_d_raw != nullptr );
          TEUCHOS_ASSERT( size_t(A_copy_d.extent(0)) ==
                          size_t(nrows) );
          TEUCHOS_ASSERT( size_t(A_copy_d.extent(1)) ==
                          size_t(ncols) );
          auto result =
            actor.factor(nrows, ncols, A_copy_d_raw,
                         A_copy_d.stride(1),
                         R.data(), R.stride(1),
                         params.contiguousCacheBlocks);
          Kokkos::deep_copy(A_copy_h, A_copy_d);
          return result;
        }
        else {
          return actor.factor(nrows, ncols, A_copy.data(),
                              A_copy.stride(1),
                              R.data(), R.stride(1),
                              params.contiguousCacheBlocks);
        }
      }();

      if(params.saveMatrices) {
        std::string filename = std::string("R") + fileSuffix;
        if(verbose) {
          cerr << "-- Save R to \"" << filename << "\"" << endl;
        }
        std::ofstream fileOut(filename.c_str());
        print_local_matrix(fileOut, ncols, ncols,
                           R.data(), R.stride(1));
        fileOut.close();
      }

      if(verbose) {
        cerr << "-- Call NodeTsqr::explicit_Q" << endl;
      }
      if(actor.wants_device_memory()) {
        const Scalar* A_copy_d_raw =
          reinterpret_cast<const Scalar*>(A_copy_d.data());
        Scalar* Q_d_raw = reinterpret_cast<Scalar*>(Q_d.data());
        TEUCHOS_ASSERT( nrows == 0 || ncols == 0 ||
                        Q_d_raw != nullptr );
        TEUCHOS_ASSERT( size_t(Q_d.extent(0)) == size_t(nrows) );
        TEUCHOS_ASSERT( size_t(Q_d.extent(1)) == size_t(ncols) );
        actor.explicit_Q(nrows, ncols,
                         A_copy_d_raw, A_copy_d.stride(1),
                         *factorOutput, ncols,
                         Q_d_raw, Q_d.stride(1),
                         params.contiguousCacheBlocks);
        // We copy back to Q_h below, either with un_cache_block (if
        // contiguous cache blocks) or directly (if not).
      }
      else {
        actor.explicit_Q(nrows, ncols,
                         A_copy.data(), A_copy.stride(1),
                         *factorOutput, ncols,
                         Q.data(), Q.stride(1),
                         params.contiguousCacheBlocks);
      }

      // "Un"-cache-block the output, if contiguous cache blocks were
      // used.  This is only necessary because local_verify() doesn't
      // currently support contiguous cache blocks.
      if(params.contiguousCacheBlocks) {
        // Use A_copy as temporary storage for un-cache-blocking Q.
        if(verbose) {
          cerr << "-- Call NodeTsqr::un_cache_block" << endl;
        }
        if(actor.wants_device_memory()) {
          Scalar* A_copy_d_raw =
            reinterpret_cast<Scalar*>(A_copy_d.data());
          const Scalar* Q_d_raw =
            reinterpret_cast<const Scalar*>(Q_d.data());
          actor.un_cache_block(nrows, ncols, A_copy_d_raw,
                               A_copy_d.stride(1), Q_d_raw);
          Kokkos::deep_copy(Q_h, A_copy_d);
        }
        else {
          actor.un_cache_block(nrows, ncols, A_copy.data(),
                               A_copy.stride(1), Q.data());
          deep_copy(Q, A_copy);
        }
      }
      else {
        if(actor.wants_device_memory()) {
          Kokkos::deep_copy(Q_h, Q_d);
        }
      }

      if(params.saveMatrices) {
        std::string filename = std::string("Q") + fileSuffix;
        if(verbose) {
          cerr << "-- Save Q to \"" << filename << "\"" << endl;
        }
        std::ofstream fileOut(filename.c_str());
        print_local_matrix(fileOut, nrows, ncols,
                           Q.data(), Q.stride(1));
        fileOut.close();
      }

      if(verbose) {
        cerr << "-- Call local_verify to validate the factorization"
             << endl;
      }
      auto results = local_verify(nrows, ncols, A.data(), lda,
                                  Q.data(), ldq, R.data(), ldr);

      if(verbose) {
        cerr << "-- Compute accuracy bounds and check" << endl;
      }

      // Accuracy relates to the number of floating-point operations,
      // which in turn is a function of the matrix's dimensions.
      // Avoid overflow of the local Ordinal type, by casting first to
      // a floating-point type.
      const mag_type dimsProd = mag_type(nrows) * mag_type(ncols) *
        mag_type(ncols);
      const mag_type fudgeFactor(10.0);
      // Relative residual error is ||A-Q*R|| / ||A||, or just
      // ||A-Q*R|| if ||A|| == 0.  (The result had better be zero in
      // the latter case.)  Square root of the matrix dimensions is an
      // old heuristic from Wilkinson or perhaps even an earlier
      // source.  We include a "fudge factor" so that the test won't
      // fail unless there is a really good reason.
      const mag_type relResidBound = fudgeFactor *
        STM::squareroot(dimsProd) * STS::eps();

      // Relative residual error; avoid division by zero.
      const mag_type relResidError = results[0] /
        (results[2] == STM::zero() ? STM::one() : results[2]);

      if(relResidError > relResidBound) {
        success = false;
        if(verbose) {
          const std::string relResStr
            (results[2] == STM::zero() ? " / ||A||_F" : "");
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
        mag_type(ncols) * mag_type(ncols) * STS::eps();

      const mag_type orthoError = results[1];
      if(orthoError > orthoBound) {
        success = false;
        if(verbose) {
          cerr << "*** For NodeTsqr=" << params.nodeTsqrType
               << " with Scalar=" << scalarType << ": "
               << "Orthogonality ||I - Q^* Q||_F = " << orthoError
               << " > bound " << orthoBound << "." << endl;
        }
      }

      if(params.humanReadable) {
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
    verifyNodeTsqr(std::ostream& out,
                   const NodeTestParameters& p)
    {
      // Seed for the next pseudorandom number generator.  We do tests
      // one after another, using the seed from the previous test in
      // the current test, so that the pseudorandom streams used by
      // the tests are independent.
      std::vector<int> iseed{{0, 0, 0, 1}};

      bool success = true;
      if(p.testReal) {
        const bool ok_S = verifyNodeTsqrTmpl<float>(out, iseed, p);
        const bool ok_D = verifyNodeTsqrTmpl<double>(out, iseed, p);
        success = success && ok_S && ok_D;
      }
      if(p.testComplex) {
#ifdef HAVE_TPETRATSQR_COMPLEX
        const bool ok_C =
          verifyNodeTsqrTmpl<std::complex<float>>(out, iseed, p);
        const bool ok_Z =
          verifyNodeTsqrTmpl<std::complex<double>>(out, iseed, p);
        success = success && ok_C && ok_Z;
#else // HAVE_TPETRATSQR_COMPLEX
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::logic_error, "TSQR was not built with complex "
           "arithmetic support.");
#endif // HAVE_TPETRATSQR_COMPLEX
      }
      return success;
    }

    template<template<class SC> class LapackType, class Scalar>
    static void
    verifyLapackTmpl(std::ostream& out,
                     std::vector<int>& iseed,
                     LapackType<Scalar>& lapack,
                     const NodeTestParameters& params,
                     const std::string& lapackImplName)
    {
      using std::cerr;
      using std::endl;
      const bool verbose = params.verbose;

      const std::string scalarType =
        Teuchos::TypeNameTraits<Scalar>::name();
      const std::string fileSuffix = getFileSuffix<Scalar>("Lapack");

      if(verbose) {
        cerr << "Test RawQR<" << scalarType << "> implementation "
             << lapackImplName << " whose type is "
             << Teuchos::typeName(lapack) << endl;
        if(lapack.wants_device_memory()) {
          cerr << "-- RawQR subclass claims to want device memory"
               << endl;
        }
      }
      const int nrows = params.numRows;
      const int ncols = params.numCols;

      Matrix<int, Scalar> A(nrows, ncols);
      Matrix<int, Scalar> A_copy(nrows, ncols);
      Matrix<int, Scalar> Q(nrows, ncols);
      Matrix<int, Scalar> R(ncols, ncols);
      if(std::numeric_limits<Scalar>::has_quiet_NaN) {
        deep_copy(A, std::numeric_limits< Scalar>::quiet_NaN());
        deep_copy(A_copy, std::numeric_limits<Scalar>::quiet_NaN());
        deep_copy(Q, std::numeric_limits<Scalar>::quiet_NaN());
        deep_copy(R, std::numeric_limits<Scalar>::quiet_NaN());
      }
      const int lda = nrows;
      const int ldq = nrows;
      const int ldr = ncols;

      if(verbose) {
        cerr << "-- Create test problem" << endl;
      }
      {
        TSQR::Random::NormalGenerator<int, Scalar> gen(iseed);
        nodeTestProblem(gen, nrows, ncols, A.data(), A.stride(1),
                        true);
        gen.getSeed(iseed); // fetch seed for the next test
      }

      if(params.saveMatrices) {
        std::string filename = std::string("A") + fileSuffix;
        if(verbose) {
          cerr << "-- Save A to \"" << filename << "\"" << endl;
        }
        std::ofstream fileOut(filename.c_str());
        print_local_matrix(fileOut, nrows, ncols,
                           A.data(), A.stride(1));
        fileOut.close();
      }

      using IST = typename Kokkos::ArithTraits<Scalar>::val_type;
      using device_matrix_type =
        Kokkos::View<IST**, Kokkos::LayoutLeft>;

      auto A_h = getHostMatrixView(A.view());
      auto A_copy_h = getHostMatrixView(A_copy.view());
      auto Q_h = getHostMatrixView(Q.view());
      device_matrix_type A_d;
      device_matrix_type A_copy_d;
      device_matrix_type Q_d;
      if(lapack.wants_device_memory()) {
        A_d = getDeviceMatrixCopy(A.view(), "A_d");
        // Don't copy A_copy yet; see below.
        A_copy_d = device_matrix_type("A_copy_d", nrows, ncols);
        Q_d = device_matrix_type("Q_d", nrows, ncols);
      }

      if(verbose) {
        cerr << "-- Copy A into A_copy" << endl;
      }
      deep_copy(A_copy, A);
      if(lapack.wants_device_memory()) {
        deep_copy(A_copy_d, A_d);
      }

      if(verbose) {
        cerr << "-- Fill R with zeros" << endl;
      }
      // We need to do this because the factorization may not
      // overwrite the strict lower triangle of R.  R is always in
      // host memory.
      deep_copy(R, Scalar {});

      if(verbose) {
        cerr << "-- Do LAPACK lwork query" << endl;
      }
      const int lwork = [&]() {
        if(lapack.wants_device_memory()) {
          Scalar* A_copy_d_raw =
            reinterpret_cast<Scalar*>(A_copy_d.data());
          const int A_copy_d_lda(A_copy_d.stride(1));
          TEUCHOS_ASSERT( nrows == 0 || ncols == 0 ||
                          A_copy_d_raw != nullptr );
          TEUCHOS_ASSERT( size_t(A_copy_d.extent(0)) ==
                          size_t(nrows) );
          TEUCHOS_ASSERT( size_t(A_copy_d.extent(1)) ==
                          size_t(ncols) );
          return lapack.compute_QR_lwork(nrows, ncols, A_copy_d_raw,
                                         A_copy_d_lda);
        }
        else {
          Scalar* A_copy_raw = A_copy.data();
          const int A_copy_lda(A_copy.stride(1));
          return lapack.compute_QR_lwork(nrows, ncols, A_copy_raw,
                                         A_copy_lda);
        }
      }();
      if(verbose) {
        cerr << "-- lwork=" << lwork << endl;
      }
      std::vector<Scalar> work(lwork);
      std::vector<Scalar> tau(ncols);

      Kokkos::View<IST*> work_d;
      Kokkos::View<IST*> tau_d;
      if(lapack.wants_device_memory()) {
        work_d = Kokkos::View<IST*>("work_d", lwork);
        tau_d = Kokkos::View<IST*>("tau_d", ncols);
      }

      if(verbose) {
        cerr << "-- Call compute_QR" << endl;
      }

      if(lapack.wants_device_memory()) {
        Scalar* A_copy_d_raw =
          reinterpret_cast<Scalar*>(A_copy_d.data());
        Scalar* tau_d_raw = reinterpret_cast<Scalar*>(tau_d.data());
        Scalar* work_d_raw =
          reinterpret_cast<Scalar*>(work_d.data());
        TEUCHOS_ASSERT( ncols == 0 || tau_d_raw != nullptr );
        TEUCHOS_ASSERT( size_t(tau_d.extent(0)) >= size_t(ncols) );
        TEUCHOS_ASSERT( lwork == 0 || work_d_raw != nullptr );
        TEUCHOS_ASSERT( size_t(work_d.extent(0)) >= size_t(lwork) );
        TEUCHOS_ASSERT( nrows == 0 || ncols == 0 ||
                        A_copy_d_raw != nullptr );
        TEUCHOS_ASSERT( size_t(A_copy_d.extent(0)) ==
                        size_t(nrows) );
        TEUCHOS_ASSERT( size_t(A_copy_d.extent(1)) ==
                        size_t(ncols) );
        lapack.compute_QR(nrows, ncols, A_copy_d_raw,
                          A_copy_d.stride(1), tau_d_raw,
                          work_d_raw, lwork);
        Kokkos::deep_copy(A_copy_h, A_copy_d);
      }
      else {
        lapack.compute_QR(nrows, ncols, A_copy.data(),
                          A_copy.stride(1), tau.data(),
                          work.data(), lwork);
      }

      if(verbose) {
        cerr << "-- Copy R out of in-place result" << endl;
      }
      copy_upper_triangle(R, A_copy);
      if(params.saveMatrices) {
        std::string filename = std::string("R") + fileSuffix;
        if(verbose) {
          cerr << "-- Save R to \"" << filename << "\"" << endl;
        }
        std::ofstream fileOut(filename.c_str());
        print_local_matrix(fileOut, ncols, ncols,
                           R.data(), R.stride(1));
        fileOut.close();
      }

      // The explicit Q factor will be computed in place, so copy the
      // result of the factorization into Q.
      deep_copy(Q, A_copy);
      if(lapack.wants_device_memory()) {
        deep_copy(Q_d, A_copy_d);
      }

      if(verbose) {
        cerr << "-- Call Lapack::compute_explicit_Q" << endl;
      }
      if(lapack.wants_device_memory()) {
        Scalar* Q_d_raw = reinterpret_cast<Scalar*>(Q_d.data());
        const Scalar* tau_d_raw =
          reinterpret_cast<const Scalar*>(tau_d.data());
        Scalar* work_d_raw =
          reinterpret_cast<Scalar*>(work_d.data());
        lapack.compute_explicit_Q(nrows, ncols, ncols,
                                  Q_d_raw, ldq, tau_d_raw,
                                  work_d_raw, lwork);
        deep_copy(Q_h, Q_d);
      }
      else {
        lapack.compute_explicit_Q(nrows, ncols, ncols,
                                  Q.data(), ldq, tau.data(),
                                  work.data(), lwork);
      }

      if(params.saveMatrices) {
        std::string filename = std::string("Q") + fileSuffix;
        if(verbose) {
          cerr << "-- Save Q to \"" << filename << "\"" << endl;
        }
        std::ofstream fileOut(filename.c_str());
        print_local_matrix(fileOut, nrows, ncols,
                           Q.data(), Q.stride(1));
        fileOut.close();
      }

      if(verbose) {
        cerr << "-- Call local_verify to validate the factorization"
             << endl;
      }
      auto results = local_verify(nrows, ncols, A.data(), lda,
                                  Q.data(), ldq, R.data(), ldr);

      if(params.humanReadable) {
        out << lapackImplName << ":" << endl
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
        out << lapackImplName
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

    template<class Scalar>
    void
    verifyLapackImplementations(std::ostream& out,
                                std::vector<int>& iseed,
                                const NodeTestParameters& p)
    {
#if defined(HAVE_TPETRATSQR_CUBLAS) && defined(HAVE_TPETRATSQR_CUSOLVER)
      {
        // Make sure that both Lapack and CuSolver get the same
        // pseudorandom seed.
        std::vector<int> iseed_copy(iseed);
        Kokkos::View<int> info("info");
        Impl::CuSolver<Scalar> solver(info.data());
        verifyLapackTmpl(out, iseed_copy, solver, p, "CUSOLVER");
      }
#endif // HAVE_TPETRATSQR_CUBLAS && HAVE_TPETRATSQR_CUSOLVER
      {
        Impl::Lapack<Scalar> lapack;
        verifyLapackTmpl(out, iseed, lapack, p, "LAPACK");
      }
    }

    void
    verifyLapack(std::ostream& out,
                 const NodeTestParameters& p)
    {
      // We do tests one after another, using the seed from the
      // previous test in the current test, so that the pseudorandom
      // streams used by the tests are independent.
      std::vector<int> iseed {{0, 0, 0, 1}};
      if(p.testReal) {
        verifyLapackImplementations<float>(out, iseed, p);
        verifyLapackImplementations<double>(out, iseed, p);
      }
      if(p.testComplex) {
#ifdef HAVE_TPETRATSQR_COMPLEX
        verifyLapackImplementations<std::complex<float>>
          (out, iseed, p);
        verifyLapackImplementations<std::complex<double>>
          (out, iseed, p);
#else // HAVE_TPETRATSQR_COMPLEX
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::logic_error, "TSQR was not built with complex "
           "arithmetic support.");
#endif // HAVE_TPETRATSQR_COMPLEX
      }
    }

    static void
    printBenchmarkFieldNames(std::ostream& out)
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

    template<template<class SC> class LapackType, class Scalar>
    void
    benchmarkLapackTmpl(std::ostream& out,
                        std::vector<int>& iseed,
                        LapackType<Scalar>& lapack,
                        const NodeTestParameters& params,
                        const std::string& lapackImplName)
    {
      using std::endl;

      const int numRows = params.numRows;
      const int numCols = params.numCols;
      const int numTrials = params.numTrials;

      Matrix<int, Scalar> A(numRows, numCols);
      Matrix<int, Scalar> Q(numRows, numCols);
      Matrix<int, Scalar> R(numCols, numCols);
      const int lda = numRows;
      const int ldq = numRows;

      {
        using prng_type = TSQR::Random::NormalGenerator<int, Scalar>;
        prng_type gen(iseed);
        nodeTestProblem(gen, numRows, numCols, A.data(), lda, false);
        gen.getSeed(iseed);
      }

      using IST = typename Kokkos::ArithTraits<Scalar>::val_type;
      using device_matrix_type =
        Kokkos::View<IST**, Kokkos::LayoutLeft>;

      auto A_h = getHostMatrixView(A.view());
      auto Q_h = getHostMatrixView(Q.view());
      device_matrix_type A_d;
      device_matrix_type Q_d;
      if(lapack.wants_device_memory()) {
        A_d = getDeviceMatrixCopy(A.view(), "A_d");
        Q_d = device_matrix_type("Q_d", numRows, numCols);
      }

      // Copy A into Q, since LAPACK QR overwrites the input.  We only
      // need Q because LAPACK's computation of the explicit Q factor
      // occurs in place.  This doesn't work with TSQR.  To give
      // LAPACK QR the fullest possible advantage over TSQR, we don't
      // allocate an A_copy here (as we would when benchmarking TSQR).
      deep_copy(Q, A);
      if(lapack.wants_device_memory()) {
        deep_copy(Q_d, A_d);
      }

      // Determine the required workspace for the factorization
      const int lwork =
        lworkQueryLapackQr(lapack, numRows, numCols, lda);
      std::vector<Scalar> work(lwork);
      std::vector<Scalar> tau(numCols);

      Kokkos::View<IST*> work_d;
      Kokkos::View<IST*> tau_d;
      if(lapack.wants_device_memory()) {
        work_d = Kokkos::View<IST*>("work_d", lwork);
        tau_d = Kokkos::View<IST*>("tau_d", numCols);
      }

      // Benchmark LAPACK's QR factorization for numTrials trials.
      Teuchos::Time timer("LAPACK");
      timer.start();
      for(int trialNum = 0; trialNum < numTrials; ++trialNum) {
        if(lapack.wants_device_memory()) {
          Scalar* Q_raw = reinterpret_cast<Scalar*>(Q_d.data());
          Scalar* tau_raw = reinterpret_cast<Scalar*>(tau_d.data());
          Scalar* work_raw =
            reinterpret_cast<Scalar*>(work_d.data());
          lapack.compute_QR(numRows, numCols,
                            Q_raw, Q_d.stride(1),
                            tau_raw, work_raw, lwork);
        }
        else {
          lapack.compute_QR(numRows, numCols,
                            Q.data(), ldq,
                            tau.data(), work.data(), lwork);
        }

        if(lapack.wants_device_memory()) {
          // FIXME (mfh 18 Dec 2019) We should actually extract the
          // upper triangle here and copy it to host, to get a fair
          // comparison with TSQR.

          Scalar* Q_raw = reinterpret_cast<Scalar*>(Q_d.data());
          const Scalar* tau_raw =
            reinterpret_cast<const Scalar*>(tau_d.data());
          Scalar* work_raw =
            reinterpret_cast<Scalar*>(work_d.data());
          lapack.compute_explicit_Q(numRows, numCols, numCols,
                                    Q_raw, Q_d.stride(1),
                                    tau_raw, work_raw, lwork);
        }
        else {
          // Extract the upper triangular factor R from Q (where it was
          // computed in place by GEQRF), since UNGQR will overwrite all
          // of Q with the explicit Q factor.
          copy_upper_triangle(R, Q);
          lapack.compute_explicit_Q(numRows, numCols, numCols,
                                    Q.data(), ldq, tau.data(),
                                    work.data(), lwork);
        }
      }
      const double lapackTiming = timer.stop();

      const std::string scalarType =
        Teuchos::TypeNameTraits<Scalar>::name();

      if(params.humanReadable) {
        out << lapackImplName << ":" << endl
            << "  Scalar: " << scalarType << endl
            << "  numRows: " << numRows << endl
            << "  numCols: " << numCols << endl
            << "  numTrials: " << numTrials << endl
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
        out << lapackImplName
            << "," << scalarType
            << "," << numRows
            << "," << numCols
            << ",0"
            << ",false"
            << "," << numTrials
            << "," << lapackTiming << endl;
      }
    }

    template<class Scalar>
    void
    benchmarkLapackImplementations(std::ostream& out,
                                   std::vector<int>& iseed,
                                   const NodeTestParameters& p)
    {
#if defined(HAVE_TPETRATSQR_CUBLAS) && defined(HAVE_TPETRATSQR_CUSOLVER)
      {
        // Make sure that both Lapack and CuSolver get the same
        // pseudorandom seed.
        std::vector<int> iseed_copy(iseed);
        Kokkos::View<int> info("info");
        Impl::CuSolver<Scalar> solver(info.data());
        benchmarkLapackTmpl(out, iseed_copy, solver, p, "CUSOLVER");
      }
#endif // HAVE_TPETRATSQR_CUBLAS && HAVE_TPETRATSQR_CUSOLVER
      {
        Impl::Lapack<Scalar> lapack;
        benchmarkLapackTmpl(out, iseed, lapack, p, "LAPACK");
      }
    }

    void
    benchmarkLapack(std::ostream& out,
                    const NodeTestParameters& p)
    {
      std::vector<int> iseed{{0, 0, 0, 1}};
      if(p.testReal) {
        benchmarkLapackImplementations<float>(out, iseed, p);
        benchmarkLapackImplementations<double>(out, iseed, p);
      }
      if(p.testComplex) {
#ifdef HAVE_TPETRATSQR_COMPLEX
        benchmarkLapackImplementations<std::complex<float>>(out, iseed, p);
        benchmarkLapackImplementations<std::complex<double>>(out, iseed, p);
#else // Don't HAVE_TPETRATSQR_COMPLEX
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::logic_error,
           "TSQR was not built with complex arithmetic support.");
#endif // HAVE_TPETRATSQR_COMPLEX
      }
    }

    template<class Scalar>
    void
    benchmarkNodeTsqrTmpl(std::ostream& out,
                          std::vector<int>& iseed,
                          NodeTsqr<int, Scalar>& actor,
                          const NodeTestParameters& params,
                          const std::string& nodeTsqrType)
    {
      using std::endl;

      TEUCHOS_ASSERT( nodeTsqrType != "" );

      const int numRows = params.numRows;
      const int numCols = params.numCols;
      const int numTrials = params.numTrials;
      const bool contiguousCacheBlocks =
        params.contiguousCacheBlocks;

      Matrix<int, Scalar> A(numRows, numCols);
      Matrix<int, Scalar> A_copy(numRows, numCols);
      Matrix<int, Scalar> Q(numRows, numCols);
      Matrix<int, Scalar> R(numCols, numCols);

      {
        using prng_type = TSQR::Random::NormalGenerator<int, Scalar>;
        prng_type gen(iseed);
        nodeTestProblem(gen, numRows, numCols,
                        A.data(), A.stride(1), false);
        gen.getSeed(iseed);
      }
      deep_copy(A_copy, A); // need copy since TSQR overwrites

      using IST = typename Kokkos::ArithTraits<Scalar>::val_type;
      using device_matrix_type =
        Kokkos::View<IST**, Kokkos::LayoutLeft>;

      auto A_copy_h = getHostMatrixView(A_copy.view());
      auto Q_h = getHostMatrixView(Q.view());
      device_matrix_type A_copy_d;
      device_matrix_type Q_d;
      if(actor.wants_device_memory()) {
        A_copy_d = getDeviceMatrixCopy(A_copy.view(), "A_copy_d");
        Q_d = device_matrix_type("Q_d", numRows, numCols);
      }

      // Benchmark sequential TSQR for numTrials trials.
      Teuchos::Time timer("NodeTsqr");
      timer.start();
      for (int trialNum = 0; trialNum < numTrials; ++trialNum) {
        if (actor.wants_device_memory()) {
          Scalar* A_raw =
            reinterpret_cast<Scalar*>(A_copy_d.data());
          auto factorOutput =
            actor.factor(numRows, numCols,
                         A_raw, A_copy_d.stride(1),
                         R.data(), R.stride(1),
                         contiguousCacheBlocks);
          // Unlike with LAPACK, this doesn't happen in place: the
          // implicit Q factor is stored in A_copy_d, and the explicit
          // Q factor is written to Q_d.
          Scalar* Q_raw = reinterpret_cast<Scalar*>(Q_d.data());
          actor.explicit_Q(numRows, numCols,
                           A_raw, A_copy_d.stride(1),
                           *factorOutput, numCols,
                           Q_raw, Q_d.stride(1),
                           contiguousCacheBlocks);
        }
        else {
          Scalar* A_raw = A_copy.data();
          auto factorOutput =
            actor.factor(numRows, numCols,
                         A_raw, A_copy.stride(1),
                         R.data(), R.stride(1),
                         contiguousCacheBlocks);
          // Unlike with LAPACK, this doesn't happen in place: the
          // implicit Q factor is stored in A_copy, and the explicit Q
          // factor is written to Q.
          Scalar* Q_raw = Q.data();
          actor.explicit_Q(numRows, numCols,
                           A_raw, A_copy.stride(1),
                           *factorOutput, numCols,
                           Q_raw, Q.stride(1),
                           contiguousCacheBlocks);
        }
      }
      const double nodeTsqrTiming = timer.stop();

      const std::string scalarType =
        Teuchos::TypeNameTraits<Scalar>::name();

      if(params.humanReadable) {
        out << "NodeTsqr:" << endl
            << "  Implementation: " << nodeTsqrType << endl
            << "  Scalar: " << scalarType << endl
            << "  numRows: " << numRows << endl
            << "  numCols: " << numCols << endl
            << "  cache size hint (bytes): "
            << params.cacheSizeHint << endl
            << "  contiguous cache blocks? "
            << (contiguousCacheBlocks ? "true" : "false") << endl
            << "  # trials = " << numTrials << endl
            << "Total time (s) = " << nodeTsqrTiming << endl
            << endl;
      }
      else {
        out << nodeTsqrType
            << "," << scalarType
            << "," << numRows
            << "," << numCols
            << "," << params.cacheSizeHint
            << "," << (contiguousCacheBlocks ? "true" : "false")
            << "," << numTrials
            << "," << nodeTsqrTiming << endl;
      }
    }

    // If nodeTsqrType == "", use p.nodeTsqrType.
    template<class Scalar>
    void
    benchmarkNodeTsqrImplementation(std::ostream& out,
                                    const std::vector<int>& iseed,
                                    const NodeTestParameters& p,
                                    const std::string& nodeTsqrType)
    {
      TEUCHOS_ASSERT( nodeTsqrType != "" );

      // Make sure that all NodeTsqr implementations get the same
      // pseudorandom seed.  That way, if there are any data-dependent
      // performance effects (e.g., subnorms), all implementations
      // will see them.
      std::vector<int> iseed_copy(iseed);
      auto nodeTsqrPtr = getNodeTsqr<Scalar>(p, nodeTsqrType);
      benchmarkNodeTsqrTmpl(out, iseed_copy, *nodeTsqrPtr, p,
                            nodeTsqrType);
    }

    template<class Scalar>
    void
    benchmarkNodeTsqrImplementations(std::ostream& out,
                                     std::vector<int>& iseed,
                                     const NodeTestParameters& p)
    {
      TEUCHOS_ASSERT( p.nodeTsqrType != "" );

      if(p.nodeTsqrType == "all" || p.nodeTsqrType == "ALL" ||
         p.nodeTsqrType == "All") {
        const char* nodeTsqrImpls[] =
          {"CombineNodeTsqr",
#if defined(HAVE_TPETRATSQR_CUBLAS) && defined(HAVE_TPETRATSQR_CUSOLVER)
           "CuSolverNodeTsqr",
#endif
           "SequentialTsqr"};
        for(auto&& nodeTsqrType : nodeTsqrImpls) {
          benchmarkNodeTsqrImplementation<Scalar>(out, iseed, p,
                                                  nodeTsqrType);
        }
      }
      else {
        benchmarkNodeTsqrImplementation<Scalar>(out, iseed, p, p.nodeTsqrType);
      }
    }

    void
    benchmarkNodeTsqr(std::ostream& out,
                      const NodeTestParameters& p)
    {
      using Teuchos::TypeNameTraits;

      std::vector<int> iseed{{0, 0, 0, 1}};
      if(p.testReal) {
        benchmarkNodeTsqrImplementations<float>(out, iseed, p);
        benchmarkNodeTsqrImplementations<double>(out, iseed, p);
      }
      if(p.testComplex) {
#ifdef HAVE_TPETRATSQR_COMPLEX
        benchmarkNodeTsqrImplementations<std::complex<float>>
          (out, iseed, p);
        benchmarkNodeTsqrImplementations<std::complex<double>>
          (out, iseed, p);
#else // Don't HAVE_TPETRATSQR_COMPLEX
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::logic_error,
           "TSQR was not built with complex arithmetic support.");
#endif // HAVE_TPETRATSQR_COMPLEX
      }
    }
  } // namespace Test
} // namespace TSQR

int
main(int argc, char *argv[])
{
  using TSQR::Test::parseOptions;
  using std::cerr;
  using std::cout;
  using std::endl;

  // Fetch command-line parameters.
  bool printedHelp = false;
  auto params = parseOptions(argc, argv, printedHelp);
  if(printedHelp) {
    return EXIT_SUCCESS;
  }

  cout << "NodeTsqr verify/benchmark test options:" << endl;
  printNodeTestParameters(cout, params, "  - ");

  bool success = true;
  try {
    Kokkos::ScopeGuard kokkosScope(argc, argv);

    // We allow the same run to do both benchmark and verify.
    if(params.verify) {
      if(! params.humanReadable) {
        TSQR::Test::printVerifyFieldNames(cout);
      }
      TSQR::Test::verifyLapack(cout, params);
      success = TSQR::Test::verifyNodeTsqr(cout, params);
    }
    if(params.benchmark) {
      if(! params.humanReadable) {
        TSQR::Test::printBenchmarkFieldNames(cout);
      }
      TSQR::Test::benchmarkLapack(cout, params);
      TSQR::Test::benchmarkNodeTsqr(cout, params);
    }

    if(params.printTrilinosTestStuff) {
      // The Trilinos test framework expects a message like this.
      if(success) {
        cout << "\nEnd Result: TEST PASSED" << endl;
      }
      else {
        cout << "\nEnd Result: TEST FAILED" << endl;
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, cerr, success);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
