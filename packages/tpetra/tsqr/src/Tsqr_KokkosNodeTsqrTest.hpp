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

#ifndef __TSQR_Test_KokkosNodeTsqrTest_hpp
#define __TSQR_Test_KokkosNodeTsqrTest_hpp

#include "Tsqr_nodeTestProblem.hpp"
#include "Tsqr_verifyTimerConcept.hpp"
#include "Tsqr_Random_NormalGenerator.hpp"
#include "Tsqr_LocalVerify.hpp"
#include "Tsqr_Matrix.hpp"
#include "Tsqr_KokkosNodeTsqr.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace TSQR {
  namespace Test {
    /// \fn verifyKokkosNodeTsqr
    /// \brief Test accuracy of KokkosNodeTsqr's QR factorization.
    ///
    /// Test the accuracy of KokkosNodeTsqr's QR factorization on a
    /// numRows by numCols matrix, and print results to stdout.
    ///
    /// \param gen [in/out] Pseudorandom number generator for the
    ///   normal(0,1) distribution.
    /// \param numRows [in] Number of rows in the test matrix.
    /// \param numCols [in] Number of columns in the test matrix.
    /// \param numPartitions [in] Number of parallel partitions (must
    ///   be a positive integer).
    /// \param cacheSizeHint [in] Cache size hint, in bytes.  Zero
    ///   means pick a reasonable default.
    /// \param contiguousCacheBlocks [in] Whether cache blocks in the
    ///   matrix to factor should be stored contiguously.
    /// \param printFieldNames [in] If humanReadable is true, this is
    ///   ignored; otherwise, whether to print a line of field names
    ///   before the line of output.
    /// \param humanReadable [in] Whether to print output that is easy
    ///   for humans to read, or instead to print output that is easy
    ///   for a script to parse.
    /// \param debug [in] Whether to print extra debugging output to
    ///   stderr.
    template<class Ordinal, class Scalar>
    void
    verifyKokkosNodeTsqr (TSQR::Random::NormalGenerator<Ordinal, Scalar>& gen,
                          const Ordinal numRows,
                          const Ordinal numCols,
                          const int numPartitions,
                          const size_t cacheSizeHint,
                          const bool contiguousCacheBlocks,
                          const bool printFieldNames,
                          const bool humanReadable,
                          const bool debug)
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;
      using Teuchos::TypeNameTraits;
      using std::cerr;
      using std::cout;
      using std::endl;
      using node_tsqr_type = TSQR::KokkosNodeTsqr<Ordinal, Scalar>;
      typedef typename node_tsqr_type::FactorOutput factor_output_type;
      typedef Teuchos::ScalarTraits<Scalar> STS;
      typedef typename STS::magnitudeType magnitude_type;
      // typedef Teuchos::Time timer_type;
      typedef Matrix<Ordinal, Scalar> matrix_type;
      typedef MatView<Ordinal, Scalar> mat_view_type;

      const std::string scalarTypeName = TypeNameTraits<Scalar>::name();

      // Set up TSQR implementation.
      RCP<ParameterList> params = parameterList ("Intranode TSQR");
      params->set ("Cache Size Hint", cacheSizeHint);
      params->set ("Num Tasks", numPartitions);
      node_tsqr_type actor (params);
      if (debug) {
        cerr << actor.description() << endl;
        if (contiguousCacheBlocks) {
          cerr << "-- Test with contiguous cache blocks" << endl;
        }
      }

      // Allocate space for test problem.
      matrix_type A (numRows, numCols);
      matrix_type A_copy (numRows, numCols);
      matrix_type Q (numRows, numCols);
      matrix_type R (numCols, numCols);
      if (std::numeric_limits<Scalar>::has_quiet_NaN) {
        deep_copy (A, std::numeric_limits<Scalar>::quiet_NaN());
        deep_copy (A_copy, std::numeric_limits<Scalar>::quiet_NaN());
        deep_copy (Q, std::numeric_limits<Scalar>::quiet_NaN());
        deep_copy (R, std::numeric_limits<Scalar>::quiet_NaN());
      }
      else {
        deep_copy (A, Scalar {});
        deep_copy (A_copy, Scalar {});
        deep_copy (Q, Scalar {});
        deep_copy (R, Scalar {});
      }
      const Ordinal lda = numRows;
      const Ordinal ldq = numRows;
      const Ordinal ldr = numCols;

      // Create a test problem
      nodeTestProblem (gen, numRows, numCols, A.data(), A.stride(1), true);

      if (debug) {
        cerr << "-- Generated test problem" << endl;
        // Don't print the matrix if it's too big.
        if (A.extent(0) <= 30) {
          cerr << "A = " << endl;
          print_local_matrix (cerr, A.extent(0), A.extent(1),
                              A.data(), A.stride(1));
          cerr << endl << endl;
        }
      }

      // Copy A into A_copy, since TSQR overwrites the input.  If
      // specified, rearrange the data in A_copy so that the data in
      // each cache block is contiguously stored.
      if (! contiguousCacheBlocks) {
        deep_copy (A_copy, A);
        if (debug) {
          cerr << "-- Copied test problem from A into A_copy" << endl;
          // Don't print the matrix if it's too big.
          if (A_copy.extent(0) <= 30) {
            cerr << "A_copy = " << endl;
            print_local_matrix (cerr, A_copy.extent(0), A_copy.extent(1),
                                A_copy.data(), A_copy.stride(1));
            cerr << endl << endl;
          }
        }
      }
      else {
        actor.cache_block (numRows, numCols, A_copy.data(), A.data(), A.stride(1));
        if (debug) {
          cerr << "-- Reorganized test matrix to have contiguous "
            "cache blocks" << endl;
          // Don't print the matrix if it's too big.
          if (A_copy.extent(0) <= 30) {
            cerr << "A_copy = " << endl;
            print_local_matrix (cerr, A_copy.extent(0), A_copy.extent(1),
                                A_copy.data(), A_copy.stride(1));
            cerr << endl << endl;
          }
        }

        // Verify cache blocking, when in debug mode.
        if (debug) {
          matrix_type A2 (numRows, numCols);
          if (std::numeric_limits<Scalar>::has_quiet_NaN) {
            deep_copy (A2, std::numeric_limits<Scalar>::quiet_NaN());
          }

          actor.un_cache_block (numRows, numCols, A2.data(), A2.stride(1), A_copy.data());
          if (matrix_equal (A, A2)) {
            if (debug)
              cerr << "-- Cache blocking test succeeded!" << endl;
          }
          else {
            if (debug) {
              cerr << "*** Cache blocking test failed! A != A2 ***"
                   << endl << endl;
              // Don't print the matrices if they are too big.
              if (A.extent(0) <= 30 && A2.extent(0) <= 30) {
                cerr << "A = " << endl;
                print_local_matrix (cerr, A.extent(0), A.extent(1),
                                    A.data(), A.stride(1));
                cerr << endl << "A2 = " << endl;
                print_local_matrix (cerr, A2.extent(0), A2.extent(1),
                                    A2.data(), A2.stride(1));
                cerr << endl;
              }
            }
            throw std::logic_error ("Cache blocking failed");
          }
        }
      }

      // Fill R with zeros, since the factorization may not
      // necessarily overwrite the strict lower triangle of R.
      if (debug) {
        cerr << "-- Filling R with zeros" << endl;
      }
      deep_copy (R, Scalar {});

      if (debug) {
        cerr << "-- Calling factor()" << endl;
      }

      // Factor the matrix and compute the explicit Q factor
      factor_output_type factor_output =
        actor.factor (numRows, numCols, A_copy.data(), A_copy.stride(1),
                      R.data(), R.stride(1), contiguousCacheBlocks);
      if (debug) {
        cerr << "-- Finished factor()" << endl;
        cerr << "-- Calling explicit_Q()" << endl;
      }

      // KokkosNodeTsqr isn't designed to be used by itself, so we
      // have to help it along by filling the top ncols x ncols
      // entries with the first ncols columns of the identity matrix.
      {
        mat_view_type Q_top =
          actor.top_block (Q.view (), contiguousCacheBlocks);
        mat_view_type Q_top_square (Q_top.extent(1), Q_top.extent(1),
                                    Q_top.data(), Q_top.stride(1));
        deep_copy (Q_top_square, Scalar {});
        for (Ordinal j = 0; j < Q_top_square.extent(1); ++j) {
          Q_top_square(j,j) = Scalar (1.0);
        }
      }
      actor.explicit_Q (numRows, numCols, A_copy.data(), A_copy.stride(1),
                        factor_output, numCols, Q.data(), Q.stride(1),
                        contiguousCacheBlocks);
      if (debug) {
        cerr << "-- Finished explicit_Q()" << endl;
      }

      // "Un"-cache-block the output Q (the explicit Q factor), if
      // contiguous cache blocks were used.  This is only necessary
      // because local_verify() doesn't currently support contiguous
      // cache blocks.
      if (contiguousCacheBlocks) {
        // Use A_copy as temporary storage for un-cache-blocking Q.
        actor.un_cache_block (numRows, numCols, A_copy.data(),
                              A_copy.stride(1), Q.data());
        deep_copy (Q, A_copy);
        if (debug) {
          cerr << "-- Un-cache-blocked output Q factor" << endl;
        }
      }

      // Print out the Q and R factors in debug mode.
      if (debug) {
        // Don't print the matrix if it's too big.
        if (Q.extent(0) <= 30) {
          cerr << endl << "-- Q factor:" << endl;
          print_local_matrix (cerr, Q.extent(0), Q.extent(1),
                              Q.data(), Q.stride(1));
          cerr << endl << endl;
        }
        cerr << endl << "-- R factor:" << endl;
        print_local_matrix (cerr, numCols, numCols, R.data(), R.stride(1));
        cerr << endl;
      }

      // Validate the factorization
      std::vector<magnitude_type> results =
        local_verify (numRows, numCols, A.data(), lda,
                      Q.data(), ldq, R.data(), ldr);
      if (debug)
        cerr << "-- Finished local_verify" << endl;

      // Print the results
      if (humanReadable) {
        cout << "KokkosNodeTsqr:" << endl
             << "Scalar type: " << scalarTypeName << endl
             << "# rows: " << numRows << endl
             << "# columns: " << numCols << endl
             << "# partitions: " << numPartitions << endl
             << "cache size hint (revised) in bytes: " << actor.cache_size_hint() << endl
             << "contiguous cache blocks? " << contiguousCacheBlocks << endl
             << "Absolute residual $\\|A - Q*R\\|_2$: "
             << results[0] << endl
             << "Absolute orthogonality $\\|I - Q^T*Q\\|_2$: "
             << results[1] << endl
             << "Test matrix norm $\\| A \\|_F$: "
             << results[2] << endl
             << endl;
      }
      else {
        if (printFieldNames) {
          const char prefix[] = "%";
          cout << prefix
               << "method"
               << ",scalarType"
               << ",numRows"
               << ",numCols"
               << ",numPartitions"
               << ",cacheSizeHint"
               << ",contiguousCacheBlocks"
               << ",absFrobResid"
               << ",absFrobOrthog"
               << ",frobA"
               << endl;
        }
        cout << "KokkosNodeTsqr"
             << "," << scalarTypeName
             << "," << numRows
             << "," << numCols
             << "," << numPartitions
             << "," << actor.cache_size_hint()
             << "," << contiguousCacheBlocks
             << "," << results[0]
             << "," << results[1]
             << "," << results[2]
             << endl;
      }
    }

    /// \fn benchmarkKokkosNodeTsqr
    /// \brief Test performance of KokkosNodeTsqr's QR factorization.
    ///
    /// Compare the performance of KokkosNodeTsqr's QR factorization
    /// to that of LAPACK's QR factorization.  Print results to
    /// stdout.
    ///
    /// \param numTrials [in] Number of times to run the benchmark;
    ///   the timing result is cumulative over all trials.  Timing
    ///   over larger numbers of trials improves certainty of the
    ///   result.
    /// \param numRows [in] Number of rows in the test matrix.
    /// \param numCols [in] Number of columns in the test matrix.
    /// \param numPartitions [in] Number of parallel partitions (must
    ///   be a positive integer).
    /// \param cacheSizeHint [in] Cache size hint, in bytes.  Zero
    ///   means pick a reasonable default.
    /// \param contiguousCacheBlocks [in] Whether cache blocks in the
    ///   matrix to factor should be stored contiguously.
    /// \param printFieldNames [in] If humanReadable is true, this is
    ///   ignored; otherwise, whether to print a line of field names
    ///   before the line of output.
    /// \param humanReadable [in] Whether to print output that is easy
    ///   for humans to read, or instead to print output that is easy
    ///   for a script to parse.
    template<class Ordinal, class Scalar>
    void
    benchmarkKokkosNodeTsqr (const int numTrials,
                             const Ordinal numRows,
                             const Ordinal numCols,
                             const int numPartitions,
                             const size_t cacheSizeHint,
                             const bool contiguousCacheBlocks,
                             const bool printFieldNames,
                             const bool humanReadable)
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;
      using Teuchos::TypeNameTraits;
      using std::cerr;
      using std::cout;
      using std::endl;
      using node_tsqr_type = TSQR::KokkosNodeTsqr<Ordinal, Scalar>;
      typedef typename node_tsqr_type::FactorOutput factor_output_type;
      typedef Teuchos::Time timer_type;
      typedef Matrix<Ordinal, Scalar> matrix_type;

      const std::string scalarTypeName = TypeNameTraits<Scalar>::name();

      // Pseudorandom normal(0,1) generator.  Default seed is OK,
      // because this is a benchmark, not an accuracy test.
      TSQR::Random::NormalGenerator<Ordinal, Scalar> gen;

      // Set up TSQR implementation.
      RCP<ParameterList> params = parameterList ("Intranode TSQR");
      params->set ("Cache Size Hint", cacheSizeHint);
      params->set ("Num Tasks", numPartitions);
      node_tsqr_type actor (params);

      // Allocate space for test problem.
      matrix_type A (numRows, numCols);
      matrix_type A_copy (numRows, numCols);
      matrix_type Q (numRows, numCols);
      matrix_type R (numCols, numCols);

      // Fill R with zeros, since the factorization may not overwrite
      // the strict lower triangle of R.
      deep_copy (R, Scalar {});

      // Create a test problem
      nodeTestProblem (gen, numRows, numCols, A.data(), A.stride(1), false);

      // Copy A into A_copy, since TSQR overwrites the input.  If
      // specified, rearrange the data in A_copy so that the data in
      // each cache block is contiguously stored.
      if (contiguousCacheBlocks) {
        actor.cache_block (numRows, numCols, A_copy.data(), A.data(), A.stride(1));
      } else {
        deep_copy (A_copy, A);
      }

      // Do a few timing runs and throw away the results, just to warm
      // up any libraries that do autotuning.
      const int numWarmupRuns = 5;
      for (int warmupRun = 0; warmupRun < numWarmupRuns; ++warmupRun) {
        // Factor the matrix in-place in A_copy, and extract the
        // resulting R factor into R.
        factor_output_type factor_output =
          actor.factor (numRows, numCols, A_copy.data(), A_copy.stride(1),
                        R.data(), R.stride(1), contiguousCacheBlocks);
        // Compute the explicit Q factor (which was stored
        // implicitly in A_copy and factor_output) and store in Q.
        // We don't need to un-cache-block the output, because we
        // aren't verifying it here.
        actor.explicit_Q (numRows, numCols, A_copy.data(), A_copy.stride(1),
                          factor_output, numCols, Q.data(), Q.stride(1),
                          contiguousCacheBlocks);
      }

      // Benchmark intranode TSQR for numTrials trials.
      //
      // Name of timer doesn't matter here; we only need the timing.
      timer_type timer("KokkosNodeTsqr");
      timer.start();
      for (int trialNum = 0; trialNum < numTrials; ++trialNum) {
        // Factor the matrix in-place in A_copy, and extract the
        // resulting R factor into R.
        factor_output_type factor_output =
          actor.factor (numRows, numCols, A_copy.data(), A_copy.stride(1),
                        R.data(), R.stride(1), contiguousCacheBlocks);
        // Compute the explicit Q factor (which was stored
        // implicitly in A_copy and factor_output) and store in Q.
        // We don't need to un-cache-block the output, because we
        // aren't verifying it here.
        actor.explicit_Q (numRows, numCols, A_copy.data(), A_copy.stride(1),
                          factor_output, numCols, Q.data(), Q.stride(1),
                          contiguousCacheBlocks);
      }
      const double timing = timer.stop();

      // Print the results
      if (humanReadable) {
        cout << "KokkosNodeTsqr cumulative timings:" << endl
             << "Scalar type: " << scalarTypeName << endl
             << "# rows = " << numRows << endl
             << "# columns = " << numCols << endl
             << "# partitions: " << numPartitions << endl
             << "Cache size hint (in bytes) = " << actor.cache_size_hint() << endl
             << "Contiguous cache blocks? " << contiguousCacheBlocks << endl
             << "# trials = " << numTrials << endl
             << "Total time (s) = " << timing << endl;
      }
      else {
        if (printFieldNames) {
          const char prefix[] = "%";
          cout << prefix
               << "method"
               << ",scalarType"
               << ",numRows"
               << ",numCols"
               << ",numPartitions"
               << ",cacheSizeHint"
               << ",contiguousCacheBlocks"
               << ",numTrials"
               << ",timing"
               << endl;
        }

        // We don't include {min,max}_seq_apply_timing() here, because
        // those times don't benefit from the accuracy of benchmarking
        // for numTrials > 1.  Thus, it's misleading to include them
        // with tbb_tsqr_timing, the total time over numTrials trials.
        cout << "KokkosNodeTsqr"
             << "," << scalarTypeName
             << "," << numRows
             << "," << numCols
             << "," << numPartitions
             << "," << actor.cache_size_hint()
             << "," << contiguousCacheBlocks
             << "," << numTrials
             << "," << timing
             << endl;
      }
    }
  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_KokkosNodeTsqrTest_hpp
