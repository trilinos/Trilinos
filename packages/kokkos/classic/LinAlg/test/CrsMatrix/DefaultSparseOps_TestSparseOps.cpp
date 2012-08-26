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

#include <Teuchos_UnitTestHarness.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_DefaultSparseOps.hpp"
#include "Kokkos_AltSparseOps.hpp"
#ifdef HAVE_KOKKOSCLASSIC_MKL
#  include "Kokkos_MklSparseOps.hpp"
#endif // HAVE_KOKKOSCLASSIC_MKL
#include "Kokkos_Version.hpp"
#include "DefaultSparseOps_TestSparseOps.hpp"


namespace {
  using Kokkos::DefaultNode;
  using Teuchos::arcp;
  using Teuchos::ArrayRCP;
  using Teuchos::null;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  //
  // mfh 28 Jun 2012: scalar_type, ordinal_type, and node_type are the
  // template parameters of interest for Kokkos local sparse kernels.
  // In practice, ordinal_type is usually int, and always should be a
  // signed integer type.  You may want to modify scalar_type if you
  // want to check correctness or performance of single-precision,
  // extended-precision, or complex arithmetic kernels.
  //

  typedef double scalar_type;
  typedef int ordinal_type;
  // mfh 28 Jun 2012: DefaultNodeType is usually TPINode, which may
  // start threads by default.  We use SerialNode to make absolutely
  // sure that this is a comparison of sequential kernels.
  //
  //typedef Kokkos::DefaultNode::DefaultNodeType node_type;
  typedef Kokkos::SerialNode node_type;

  typedef Teuchos::ScalarTraits<double> STM;

  //
  // Values of command-line arguments.
  //
  // CommandLineProcessor only accepts double for floating-point
  // values, and int for integer values, so we don't use scalar_type
  // or ordinal_type here.

  // If not "", then we assume this is the name of a Matrix Market -
  // format sparse matrix file to read.
  std::string matrixFilename;
  // The given tolerance allows some rounding error for triangular
  // solves, given that the test problems we generate should be well
  // conditioned on average.
  double tol = 1000 * Teuchos::ScalarTraits<double>::eps ();
  // Number of rows in the sparse matrices to test.
  int numRows = 5;
  // Number of columns in the sparse matrices to test.
  int numCols = 5;
  // Number of benchmark trials; timings are cumulative over all trials.
  int numTrials = 100;
  // If true, run the benchmark instead of the test.
  bool benchmark = false;
  // Number of columns in the multivectors to benchmark.
  int numVecs = 1;
  // Number of threads for the Kokkos Node instance.  -1 means use the
  // default behavior.
  int numThreads = -1;
  // Whether to force AltSparseOps to use first-touch allocation.
  bool forceFirstTouch = false;
  // Verbosity, including that of Kokkos Node initialization.  (Only
  // certain Kokkos Node types support this option.)
  bool verbose = false;
  // Whether to print copious debugging output.
  bool debug = false;
  // For AltSparseOps: Whether to use the version with loops unrolled
  // over the input and output multivectors.
  bool unroll = true;
  // For AltSparseOps: Which variant of sparse mat-vec to test; "all"
  // means test all implemented variants.
  std::string variant ("for-for");

  // This will be set after reading the command-line arguments, so
  // that we can, for the particular Kokkos Node type, set things like
  // the number of threads to run.
  RCP<node_type> node;

  // Set up the extra command-line arguments.
  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP ();
    clp.addOutputSetupOptions (true);
    clp.setOption ("matrixFilename", &matrixFilename, "If not \"\", the name "
                   "of a Matrix Market - format file from which to read the "
                   "sparse matrix A to benchmark or test.  If provided, we only"
                   "benchmark or test sparse mat-vec, not sparse triangular "
                   "solve.");
    clp.setOption ("numRows", &numRows, "Number of rows in the matrices to "
                   "test.  Ignored if a matrix file name is provided.");
    clp.setOption ("numCols", &numCols, "Number of columns in the matrices to "
                   "test.  Ignored if a matrix file name is provided.  We only "
                   "test sparse triangular solve if numRows == numCols.");
    clp.setOption ("tol", &tol, "Tolerance for relative error.");
    clp.setOption ("benchmark", "test", &benchmark, "Whether to run the "
                   "benchmark (if true) or the test (if false).");
    clp.setOption ("numTrials", &numTrials, "If running the benchmark: number "
                   "of benchmark trials.  Timings are cumulative over all "
                   "trials.");
    clp.setOption ("numVecs", &numVecs, "Number of columns in the multivectors "
                   "to benchmark.");
    clp.setOption ("numThreads", &numThreads, "Number of threads for the Kokkos "
                   "Node instance to use.  -1 means use the default behavior.");
    clp.setOption ("verbose", "quiet", &verbose, "Whether Kokkos Node "
                   "initialization should print verbose status output, if the "
                   "Node type supports this.");
    clp.setOption ("debug", "release", &debug, "Whether to print copious "
                   "debugging output.");
    clp.setOption ("unroll", "dontUnroll", &unroll, "Whether Kokkos::SeqSparse"
                   "Ops should unroll across columns of the multivectors.");
    clp.setOption ("variant", &variant, "Which algorithm variant Kokkos::Seq"
                   "SparseOps should use for sparse mat-vec.  Valid options "
                   "are \"for-for\", \"for-while\", and \"for-if\".  You may "
                   "also use \"all\", which tests all possibilities.");
    clp.setOption ("forceFirstTouch", "dontForceFirstTouch", &forceFirstTouch,
                   "Whether to force AltSparseOps to use first-touch allocation.");

    ParameterList nodeParams;
    if (numThreads != -1) {
      nodeParams.set ("Num Threads", numThreads);
    }
    nodeParams.set ("Verbose", verbose);
    node = rcp (new node_type (nodeParams));
  }

  //
  // UNIT TESTS
  //

  /// \class Tester
  /// \brief Class for testing SparseOpsType.
  ///
  /// \tparam SparseOpsType Implementation of local sparse kernels.
  ///   Examples include Kokkos::DefaultHostSparseOps.
  template<class SparseOpsType>
  class Tester {
  public:
    /// \brief Benchmark or test SparseOpsType, with default
    ///   SparseOpsType parameters.
    ///
    /// \param sparseOpsTypeName [in] Name of the SparseOpsType class;
    ///   used for generating timer labels in the benchmark.
    ///
    /// \param implicitUnitDiagTriMultCorrect [in] Whether
    ///   SparseOpsType correctly implements sparse
    ///   matrix-(multi)vector multiply with a triangular matrix with
    ///   implicitly stored unit diagonal.  "Incorrectly" means that
    ///   SparseOpsType assumes that all entries of the matrix are
    ///   stored explicitly, regardless of its Teuchos::EDiag input
    ///   value.  The SparseOpsType interface does not require this by
    ///   default, but some implementations do.
    static void
    test (const std::string& sparseOpsTypeName,
          const bool implicitUnitDiagTriMultCorrect=false)
    {
      Teuchos::ParameterList params;
      test (sparseOpsTypeName, params, implicitUnitDiagTriMultCorrect);
    }

    /// \brief Benchmark or test SparseOpsType.
    ///
    /// \param sparseOpsTypeName [in] Name of the SparseOpsType class;
    ///   used for generating timer labels in the benchmark.
    /// \param params [in/out] Parameters for SparseOpsType's constructor.
    /// \param implicitUnitDiagTriMultCorrect [in] Whether
    ///   SparseOpsType correctly implements sparse
    ///   matrix-(multi)vector multiply with a triangular matrix with
    ///   implicitly stored unit diagonal.  "Incorrectly" means that
    ///   SparseOpsType assumes that all entries of the matrix are
    ///   stored explicitly, regardless of its Teuchos::EDiag input
    ///   value.  The SparseOpsType interface does not require this by
    ///   default, but some implementations do.
    static void
    test (const std::string& sparseOpsTypeName,
          Teuchos::ParameterList& params,
          const bool implicitUnitDiagTriMultCorrect=false)
    {
      using Teuchos::FancyOStream;
      using Teuchos::getFancyOStream;
      using Teuchos::RCP;
      using Teuchos::rcpFromRef;

      RCP<FancyOStream> out = getFancyOStream (rcpFromRef (std::cout));
      TestSparseOps<SparseOpsType> tester (out, verbose, debug);
      if (benchmark) {
        std::vector<std::pair<std::string, double> > results;
        if (matrixFilename == "") {
          tester.benchmarkSparseOps (results, sparseOpsTypeName, node, params,
                                     numRows, numCols, numVecs, numTrials);
        }
        else {
          tester.benchmarkSparseMatVecFromFile (results, matrixFilename,
                                                sparseOpsTypeName, node,
                                                params, numVecs, numTrials);
        }
      }
      else {
        tester.testSparseOps (node, params, numRows, numCols, numVecs, tol,
                              implicitUnitDiagTriMultCorrect);
      }
    }
  };

  // Test sparse matrix-(multi)vector multiply and sparse triangular
  // solve, using DefaultHostSparseOps.
  TEUCHOS_UNIT_TEST( DefaultSparseOps, TestSparseOps )
  {
    using Kokkos::DefaultHostSparseOps;
    typedef DefaultHostSparseOps<scalar_type, ordinal_type, node_type> sparse_ops_type;

    const bool implicitUnitDiagTriMultCorrect = false;
    Tester<sparse_ops_type>::test ("DefaultSparseHostOps",
                                   implicitUnitDiagTriMultCorrect);
  }

#ifdef HAVE_KOKKOSCLASSIC_MKL
  // Test sparse matrix-(multi)vector multiply and sparse triangular
  // solve, using MklSparseOps.
  TEUCHOS_UNIT_TEST( MklSparseOps, TestSparseOps )
  {
    using Kokkos::MklSparseOps;
    typedef MklSparseOps<scalar_type, ordinal_type, node_type> sparse_ops_type;

    // MKL correctly implements both sparse mat-vec and sparse
    // triangular solve for the implicit unit diagonal triangular
    // matrix case.
    const bool implicitUnitDiagTriMultCorrect = true;

    // A critical test for MklSparseOps is to exercise the on-demand
    // transition from 0-based indices to 1-based indices.  We can do
    // this by starting with numVecs=1, then going to some numVecs > 1
    // and back again.  In order not to spoil tests after these, we
    // save the current numVecs value and restore it after the tests
    // (with a try-catch, so it gets restored even if the tests throw.
    const int curNumVecs = numVecs;
    try {
      numVecs = 1;
      Tester<sparse_ops_type>::test ("MklSparseOps",
                                     implicitUnitDiagTriMultCorrect);
      numVecs = 3;
      Tester<sparse_ops_type>::test ("MklSparseOps",
                                     implicitUnitDiagTriMultCorrect);
      numVecs = 1;
      Tester<sparse_ops_type>::test ("MklSparseOps",
                                     implicitUnitDiagTriMultCorrect);
      // Also test the user's numVecs, if it's neither 1 nor 3.
      if (curNumVecs != 1 && curNumVecs != 3) {
        numVecs = curNumVecs;
        Tester<sparse_ops_type>::test ("MklSparseOps",
                                       implicitUnitDiagTriMultCorrect);
      }
    }
    catch (...) {
      numVecs = curNumVecs;
      throw;
    }
    numVecs = curNumVecs; // "finally" block of the try-catch.
  }
#endif // HAVE_KOKKOSCLASSIC_MKL

  // Test sparse matrix-(multi)vector multiply and sparse triangular solve.
  // This test should work, but we omit it to save time.
#if 0
  TEUCHOS_UNIT_TEST( AltSparseOpsDefaultParameters, TestSparseOps )
  {
    using Kokkos::AltSparseOps;
    typedef AltSparseOps<scalar_type, ordinal_type, node_type> sparse_ops_type;
    const bool implicitUnitDiagTriMultCorrect = true;
    Tester<sparse_ops_type>::test ("AltSparseOps",
                                   implicitUnitDiagTriMultCorrect);
  }
#endif // 0

  // Test sparse matrix-(multi)vector multiply and sparse triangular solve.
  TEUCHOS_UNIT_TEST( AltSparseOpsNonDefaultParameters, TestSparseOps )
  {
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Kokkos::AltSparseOps;
    typedef AltSparseOps<scalar_type, ordinal_type, node_type> sparse_ops_type;
    const bool implicitUnitDiagTriMultCorrect = false;

    ParameterList params ("AltSparseOps");
    params.set ("Unroll across multivectors", unroll);
    params.set ("Force first-touch allocation", forceFirstTouch);
    if (variant == "all") {
      std::string variant = "for-for";
      params.set ("Sparse matrix-vector multiply variant", variant);
      Tester<sparse_ops_type>::test ("AltSparseOps (for-for)", params,
                                     implicitUnitDiagTriMultCorrect);
      variant = "for-while";
      params.set ("Sparse matrix-vector multiply variant", variant);
      Tester<sparse_ops_type>::test ("AltSparseOps (for-while)", params,
                                     implicitUnitDiagTriMultCorrect);
      variant = "for-if";
      params.set ("Sparse matrix-vector multiply variant", variant);
      Tester<sparse_ops_type>::test ("AltSparseOps (for-if)", params,
                                     implicitUnitDiagTriMultCorrect);
    }
    else {
      params.set ("Sparse matrix-vector multiply variant", variant);
      Tester<sparse_ops_type>::test ("AltSparseOps", params,
                                     implicitUnitDiagTriMultCorrect);
    }

    if (benchmark) {
      // Summarize timing results.  You should only call summarize()
      // for all the different SparseOps types you plan to test in
      // this executable.
      //
      // Extra endline makes the summarize() output nicer.  In
      // benchmark mode, you should run with only one MPI process (if
      // building with MPI at all).  We don't plan to run benchmark
      // mode for the check-in or continuous integration tests.
      std::cout << "\n";
      Teuchos::TimeMonitor::summarize ();
    }
  }
}

