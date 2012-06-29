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
#include "Kokkos_SeqSparseOps.hpp"
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

  // The given tolerance allows some rounding error for triangular
  // solves, given that the test problems we generate should be well
  // conditioned on average.
  double tol = 1000 * Teuchos::ScalarTraits<double>::eps ();
  // Number of rows (and columns) in the sparse matrices to test.
  int numRows = 100;
  // Number of benchmark trials; timings are cumulative over all trials.
  int numTrials = 100;
  // If true, run the benchmark instead of the test.
  bool benchmark = false;
  // Number of columns in the multivectors to benchmark.
  // Must be at least as many as the number of rows.
  int numVecs = 1;

  // Set up the extra command-line arguments.
  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP ();
    clp.addOutputSetupOptions (true);
    clp.setOption ("numRows", &numRows, "Number of rows (and columns) in the "
                   "matrices to test.");
    clp.setOption ("tol", &tol, "Tolerance for relative error.");
    clp.setOption ("benchmark", "test", &benchmark, "Whether to run the "
                   "benchmark (if true) or the test (if false).");
    clp.setOption ("numTrials", &numTrials, "If running the benchmark: number "
                   "of benchmark trials.  Timings are cumulative over all "
                   "trials.");
    clp.setOption ("numVecs", &numVecs, "Number of columns in the multivectors "
                   "to benchmark.");
  }

  //
  // UNIT TESTS
  //

  // Test sparse matrix-(multi)vector multiply and sparse triangular
  // solve, using DefaultHostSparseOps.
  TEUCHOS_UNIT_TEST( DefaultSparseOps, TestSparseOps )
  {
    using Kokkos::DefaultHostSparseOps;
    typedef DefaultHostSparseOps<scalar_type, ordinal_type, node_type> sparse_ops_type;
    RCP<node_type> node = Kokkos::DefaultNode::getDefaultNode();
    TestSparseOps<sparse_ops_type> tester;

    if (benchmark) {
      std::vector<std::pair<std::string, double> > results;
      tester.benchmarkSparseOps (results, "DefaultSparseHostOps", node,
                                 numRows, numVecs, numTrials);
    }
    else {
      tester.testSparseOps (node, numRows, numVecs, tol);
    }
  }

  // Test sparse matrix-(multi)vector multiply and sparse triangular solve.
  TEUCHOS_UNIT_TEST( SeqSparseOps, TestSparseOps )
  {
    using Kokkos::SeqSparseOps;
    typedef SeqSparseOps<scalar_type, ordinal_type, node_type> sparse_ops_type;
    RCP<node_type> node = Kokkos::DefaultNode::getDefaultNode();
    TestSparseOps<sparse_ops_type> tester;

    if (benchmark) {
      std::vector<std::pair<std::string, double> > results;
      tester.benchmarkSparseOps (results, "SeqSparseOps", node,
                                 numRows, numVecs, numTrials);
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
    else {
      // The test runs by default.
      tester.testSparseOps (node, numRows, numVecs, tol);
    }
  }
}
