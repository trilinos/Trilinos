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

#include "Tsqr_CombineTest.hpp"

#include "Tsqr_Random_NormalGenerator.hpp"
#include "Tsqr_Random_MatrixGenerator.hpp"

#include "Tsqr_Combine.hpp"
#include "Tsqr_LocalVerify.hpp"
#include "Tsqr_Matrix.hpp"
#include "Tsqr_Util.hpp"

#include "Teuchos_Assert.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace TSQR {
  namespace Test {

    template<class Ordinal, class MagnitudeType, class NormalGenType>
    static void
    generateSingularValues (NormalGenType& magGen,
                            std::vector<MagnitudeType>& sigma,
                            const Ordinal numValues)
    {
      typedef MagnitudeType magnitude_type;
      const magnitude_type machEps =
        std::numeric_limits<magnitude_type>::epsilon();
      sigma.resize (numValues);

      // Relative amount by which to perturb each singular value.  The
      // perturbation will be multiplied by a normal(0,1) pseudorandom
      // number drawn from magGen.
      const magnitude_type perturbationFactor = magnitude_type(10) * machEps;

      sigma[0] = magnitude_type (1);
      for (Ordinal k = 1; k < numValues; ++k)
        {
          const magnitude_type perturbation = perturbationFactor * magGen();
          const magnitude_type beforePerturb = sigma[k-1] / magnitude_type(2);
          const magnitude_type candidate = beforePerturb + perturbation;

          // If adding the perturbation to beforePerturb would result
          // in a nonpositive number, subtract instead.
          if (candidate <= magnitude_type(0))
            sigma[k] = beforePerturb - perturbation;
          else
            sigma[k] = candidate;
        }
    }

    static void
    printCombineFieldNames ()
    {
      using std::cout;
      using std::endl;

      const char prefix[] = "%";
      cout << prefix
           << "method"
           << ",kernel"
           << ",scalarType"
           << ",numRows"
           << ",numCols"
           << ",absFrobResid"
           << ",absFrobOrthog"
           << ",frobA"
           << endl;
    }

    template<class MagnitudeType>
    static void
    printR1R2results (const std::string& datatype,
                      const int numCols,
                      const std::vector<MagnitudeType>& results)
    {
      using std::cout;
      using std::endl;

      cout << "Combine"
           << "," << "R1R2"
           << "," << datatype
           << "," << (2*numCols)
           << "," << numCols
           << "," << results[0]
           << "," << results[1]
           << "," << results[2]
           << endl;
    }

    template<class MagnitudeType>
    static void
    printR3Aresults (const std::string& datatype,
                     const int numRows,
                     const int numCols,
                     const std::vector<MagnitudeType>& results)
    {
      using std::cout;
      using std::endl;

      cout << "Combine"
           << "," << "R3A"
           << "," << datatype
           << "," << numRows
           << "," << numCols
           << "," << results[3]
           << "," << results[4]
           << "," << results[5]
           << endl;
    }

    template<class MagnitudeType>
    static void
    printResults (const std::string& datatype,
                  const int numRows,
                  const int numCols,
                  const std::vector<MagnitudeType>& results,
                  const bool printFieldNames)
    {
      if (printFieldNames)
        printCombineFieldNames();
      printR1R2results (datatype, numCols, results);
      printR3Aresults (datatype, numRows, numCols, results);
    }

    template<class MagnitudeType>
    static void
    printSimSeqTsqrResults (const std::string& datatype,
                            const int numRows,
                            const int numCols,
                            const std::vector<MagnitudeType>& results,
                            const bool printFieldNames)
    {
      using std::cout;
      using std::endl;

      if (printFieldNames)
        {
          const char prefix[] = "%";
          cout << prefix
               << "method"
               << ",scalarType"
               << ",numRows"
               << ",numCols"
               << ",absFrobResid"
               << ",absFrobOrthog"
               << ",frobA"
               << endl;
        }
      cout << "CombineSimSeqTsqr"
           << "," << datatype
           << "," << numRows
           << "," << numCols
           << "," << results[0]
           << "," << results[1]
           << "," << results[2]
           << endl;
    }

    template< class MatrixViewType >
    static void
    printMatrix (std::ostream& out,
                 const MatrixViewType& A)
    {
      print_local_matrix (out, A.extent(0), A.extent(1), A.data(), A.stride(1));
    }

    template<class MatrixViewType>
    static
    std::vector<
      typename Teuchos::ScalarTraits<
        typename MatrixViewType::non_const_value_type
      >::magnitudeType
    >
    localVerify (const MatrixViewType& A,
                 const MatrixViewType& Q,
                 const MatrixViewType& R)
    {
      return local_verify (A.extent(0), A.extent(1), A.data(), A.stride(1),
                           Q.data(), Q.stride(1), R.data(), R.stride(1));
    }

    /// \brief Test accuracy of TSQR::Combine
    ///
    /// 1. [R1; R2] where R1 and R2 are both ncols by ncols upper
    ///    triangular.
    ///
    /// 2. [R; A] where R is ncols by ncols upper triangular, and A is
    ///    nrows by ncols general dense.
    ///
    /// \return ($\|A - QR\|_F$, $\|I - Q^* Q\|_F$, $\|A\|_F$) for each
    ///   test problem (so, a vector of six elements).
    ///
    template<class Ordinal, class Scalar>
    static std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
    verifyCombineTemplate (TSQR::Random::NormalGenerator<Ordinal, Scalar>& gen,
                           TSQR::Random::NormalGenerator<Ordinal, typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& magGen,
                           const Ordinal numRows,
                           const Ordinal numCols,
                           const bool debug)
    {
      using TSQR::Random::MatrixGenerator;
      using TSQR::Random::NormalGenerator;
      using std::cerr;
      using std::endl;
      using std::invalid_argument;
      using std::ostringstream;
      using std::pair;
      using std::vector;

      typedef Teuchos::ScalarTraits<Scalar> STS;
      typedef typename STS::magnitudeType magnitude_type;
      typedef NormalGenerator<Ordinal, Scalar> normgen_type;
      typedef MatrixGenerator<Ordinal, Scalar, normgen_type> matgen_type;
      typedef Matrix<Ordinal, Scalar> matrix_type;
      typedef vector<magnitude_type> results_type;

      if (numRows < numCols) {
        ostringstream os;
        os << "# rows < # columns is not allowed.  You specified # rows = "
           << numRows << " and # columns = " << numCols << ".";
        throw invalid_argument (os.str());
      }
      else if (numCols == 0) {
        throw invalid_argument ("ncols == 0 is not allowed");
      }

      //
      // Generate four different sets of singular values.  Randomly
      // perturb them, but make sure all are positive.
      //
      vector< magnitude_type > sigma_R1 (numCols);
      vector< magnitude_type > sigma_R2 (numCols);
      vector< magnitude_type > sigma_R3 (numCols);
      vector< magnitude_type > sigma_A (numCols);
      generateSingularValues (magGen, sigma_R1, numCols);
      generateSingularValues (magGen, sigma_R2, numCols);
      generateSingularValues (magGen, sigma_R3, numCols);
      generateSingularValues (magGen, sigma_A, numCols);

      matrix_type R1 (numCols, numCols, Scalar(0));
      matrix_type R2 (numCols, numCols, Scalar(0));
      matrix_type R3 (numCols, numCols, Scalar(0));
      matrix_type A (numRows, numCols, Scalar(0));
      matgen_type matgen (gen);
      matgen.fill_random_R (numCols, R1.data(), R1.stride(1), &sigma_R1[0]);
      matgen.fill_random_R (numCols, R2.data(), R2.stride(1), &sigma_R2[0]);
      matgen.fill_random_R (numCols, R3.data(), R3.stride(1), &sigma_R3[0]);
      matgen.fill_random_svd (numRows, numCols, A.data(), A.stride(1), &sigma_A[0]);

      if (false && debug) {
        cerr << endl << "First test problem:" << endl;
        print_local_matrix (cerr, numCols, numCols, R1.data(), R1.stride(1));
        print_local_matrix (cerr, numCols, numCols, R2.data(), R2.stride(1));
        cerr << endl;

        cerr << endl << "Second test problem:" << endl;
        print_local_matrix (cerr, numCols, numCols, R3.data(), R3.stride(1));
        print_local_matrix (cerr, numRows, numCols, A.data(), A.stride(1));
        cerr << endl;
      }

      // Space to put the original test problem, expressed as one
      // dense matrix rather than in two blocks.  These will be deep
      // copies of the test problems, since the test problem matrices
      // will be overwritten by the factorizations.
      matrix_type A_R1R2 (Ordinal(2) * numCols, numCols, Scalar {});
      matrix_type A_R3A (numRows + numCols, numCols, Scalar {});

      // Copy [R1; R2] into A_R1R2.
      {
        auto A_R1R2_views = partition_2x1 (A_R1R2, numCols);
        deep_copy (A_R1R2_views.first, R1);
        deep_copy (A_R1R2_views.second, R2);
      }

      // Copy [R3; A] into A_R3A.
      {
        auto A_R3A_views = partition_2x1 (A_R3A, numCols);
        deep_copy (A_R3A_views.first, R3);
        deep_copy (A_R3A_views.second, A);
      }

      // Space to put the explicit Q factors.
      matrix_type Q_R1R2 (Ordinal(2) * numCols, numCols, Scalar(0));
      matrix_type Q_R3A (numRows + numCols, numCols, Scalar(0));

      // Fill the explicit Q factor matrices with the first numCols
      // columns of the identity matrix.
      for (Ordinal k = 0; k < numCols; ++k) {
        // FIXME (mfh 26 Nov 2019) Eventually we want to get away from
        // direct modification of the entries of a Matrix or MatView,
        // in favor of only doing so with a Kokkos kernel or TPL.
        Q_R1R2(k, k) = Scalar(1.0);
        Q_R3A(k, k) = Scalar(1.0);
      }

      // tau factor arrays, one for each factorization test.
      vector<Scalar> tau_R1R2 (numCols);
      vector<Scalar> tau_R3A (numCols);

      // Workspace array for factorization and applying the Q factor.
      // We recycle this workspace for all tests.
      vector<Scalar> work (numCols);

      if (debug) {
        cerr << endl << "----------------------------------------" << endl
             << "TSQR::Combine first test problem:" << endl
             << "qr( [R1; R2] ), with R1 and R2 " << numCols
             << " by " << numCols << endl << endl;
      }
      Combine<Ordinal, Scalar> combiner;
      combiner.factor_pair (R1.view(), R2.view(),
                            tau_R1R2.data(), work.data());
      combiner.apply_pair (ApplyType("N"), numCols, numCols,
                           R2.data(), R2.stride(1), tau_R1R2.data(),
                           &Q_R1R2(0, 0), Q_R1R2.stride(1),
                           &Q_R1R2(numCols, 0), Q_R1R2.stride(1),
                           work.data());
      if (debug) {
        cerr << "Results of first test problem:" << endl;
        cerr << "-- Copy of test problem:" << endl;
        print_local_matrix (cerr, A_R1R2.extent(0), A_R1R2.extent(1),
                            A_R1R2.data(), A_R1R2.stride(1));
        cerr << endl << "-- Q factor:" << endl;
        print_local_matrix (cerr, Q_R1R2.extent(0), Q_R1R2.extent(1),
                            Q_R1R2.data(), Q_R1R2.stride(1));
        cerr << endl << "-- R factor:" << endl;
        print_local_matrix (cerr, R1.extent(0), R1.extent(1),
                            R1.data(), R1.stride(1));
        cerr << endl;
      }
      const results_type firstResults =
        local_verify (A_R1R2.extent(0), A_R1R2.extent(1),
                      A_R1R2.data(), A_R1R2.stride(1),
                      Q_R1R2.data(), Q_R1R2.stride(1),
                      R1.data(), R1.stride(1));
      if (debug) {
        cerr << "\\| A - Q*R \\|_F = " << firstResults[0] << endl
             << "\\| I - Q'*Q \\|_F = " << firstResults[1] << endl
             << "\\| A \\|_A = " << firstResults[2] << endl;
        cerr << endl << "----------------------------------------" << endl
             << "TSQR::Combine second test problem:" << endl
             << "qr( [R3; A] ), with R3 " << numCols << " by " << numCols
             << " and A " << numRows << " by " << numCols << endl << endl;
      }
      combiner.factor_inner (R3.view(), A.view(),
                             tau_R3A.data(), work.data());
      combiner.apply_inner (ApplyType("N"), numRows, numCols, numCols,
                            A.data(), A.stride(1), tau_R3A.data(),
                            &Q_R3A(0, 0), Q_R3A.stride(1),
                            &Q_R3A(numCols, 0), Q_R3A.stride(1),
                            work.data());
      if (debug) {
        cerr << "Results of second test problem:" << endl;
        cerr << "-- Copy of test problem:" << endl;
        print_local_matrix (cerr, A_R3A.extent(0), A_R3A.extent(1),
                            A_R3A.data(), A_R3A.stride(1));
        cerr << endl << "-- Q factor:" << endl;
        print_local_matrix (cerr, Q_R3A.extent(0), Q_R3A.extent(1),
                            Q_R3A.data(), Q_R3A.stride(1));
        cerr << endl << "-- R factor:" << endl;
        print_local_matrix (cerr, R3.extent(0), R3.extent(1),
                            R3.data(), R3.stride(1));
        cerr << endl;
      }
      const results_type secondResults =
        local_verify (A_R3A.extent(0), A_R3A.extent(1),
                      A_R3A.data(), A_R3A.stride(1),
                      Q_R3A.data(), Q_R3A.stride(1),
                      R3.data(), R3.stride(1));
      if (debug) {
        cerr << "\\| A - Q*R \\|_F = " << secondResults[0] << endl
             << "\\| I - Q'*Q \\|_F = " << secondResults[1] << endl
             << "\\| A \\|_A = " << secondResults[2] << endl;
      }
      vector<magnitude_type> finalResults;
      finalResults.push_back (firstResults[0]);
      finalResults.push_back (firstResults[1]);
      finalResults.push_back (firstResults[2]);

      finalResults.push_back (secondResults[0]);
      finalResults.push_back (secondResults[1]);
      finalResults.push_back (secondResults[2]);
      return finalResults;
    }

    //! Simulate one combine step of Sequential TSQR
    template<class Ordinal, class Scalar>
    static std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
    verifyCombineSeqTemplate (TSQR::Random::NormalGenerator<Ordinal, Scalar>& gen,
                              TSQR::Random::NormalGenerator<Ordinal, typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& magGen,
                              const Ordinal numRows,
                              const Ordinal numCols,
                              const bool debug)
    {
      using TSQR::Random::MatrixGenerator;
      using TSQR::Random::NormalGenerator;
      using std::cerr;
      using std::endl;
      using std::invalid_argument;
      using std::ostringstream;
      using std::pair;
      using std::vector;

      typedef Teuchos::ScalarTraits<Scalar> STS;
      typedef typename STS::magnitudeType magnitude_type;
      typedef NormalGenerator< Ordinal, Scalar > normgen_type;
      typedef MatrixGenerator< Ordinal, Scalar, normgen_type > matgen_type;
      typedef Matrix<Ordinal, Scalar> matrix_type;
      typedef MatView<Ordinal, Scalar> mat_view_type;
      typedef vector<magnitude_type> results_type;

      if (numRows < numCols) {
        ostringstream os;
        os << "# rows < # columns is not allowed.  You specified # rows = "
           << numRows << " and # columns = " << numCols << ".";
        throw invalid_argument (os.str());
      }
      else if (numCols == 0) {
        throw invalid_argument ("ncols == 0 is not allowed");
      }

      // Generate two different sets of singular values.
      vector< magnitude_type > sigma_A1 (numCols);
      vector< magnitude_type > sigma_A2 (numCols);
      generateSingularValues (magGen, sigma_A1, numCols);
      generateSingularValues (magGen, sigma_A2, numCols);

      // Matrix consisting of two cache blocks.
      matrix_type A (Ordinal(2)*numRows, numCols, Scalar(0));
      // Views of the two cache blocks.
      mat_view_type A1 (numRows, numCols, &A(0,0), A.stride(1));
      mat_view_type A2 (numRows, numCols, &A(numRows,0), A.stride(1));

      // Fill the two cache blocks with random test problems.
      matgen_type matgen (gen);
      matgen.fill_random_svd (numRows, numCols, A1.data(), A1.stride(1), &sigma_A1[0]);
      matgen.fill_random_svd (numRows, numCols, A2.data(), A2.stride(1), &sigma_A2[0]);

      if (false && debug) {
        cerr << endl << "Test problem:" << endl;
        cerr << endl << "Original matrix:" << endl;
        printMatrix (cerr, A);
        cerr << endl << "First cache block:" << endl;
        printMatrix (cerr, A1);
        cerr << endl << "Second cache block:" << endl;
        printMatrix (cerr, A2);
        cerr << endl;
      }

      // Copy of the resulting test problem, stored as one dense
      // matrix rather than as two blocks.  We will use A_copy to
      // measure the residual error once we've completed the
      // factorization and computed the explicit Q factor.
      matrix_type A_copy (A);

      // Space to put the explicit Q factor.
      matrix_type Q (Ordinal(2) * numRows, numCols, Scalar(0));

      // Fill Q with the first numCols columns of the identity matrix.
      for (Ordinal k = 0; k < numCols; ++k) {
        // FIXME (mfh 26 Nov 2019) I'm assuming I can write to the
        // Matrix or MatView on host, outside of Kokkos.  TSQR always
        // assumed this, but if we want to use Kokkos, we'll need to
        // get rid of that assumption.
        Q(k, k) = Scalar(1.0);
      }

      // Two cache blocks (as views) of Q.
      mat_view_type Q1 (numRows, numCols, &Q(0,0), Q.stride(1));
      mat_view_type Q2 (numRows, numCols, &Q(numRows,0), Q.stride(1));

      // Two tau factor arrays, one for each cache block.
      vector<Scalar> tau1 (numCols);
      vector<Scalar> tau2 (numCols);

      // Workspace array for factorization and applying the Q factor.
      // We recycle this workspace for all tests.
      vector<Scalar> work (numCols);

      if (debug) {
        cerr << endl << "----------------------------------------" << endl
             << "TSQR::Combine SequentialTsqr simulation with 2 cache blocks:"
             << endl << "qr( [A1; A2] ), with A1 and A2 being each "
             << numRows << " by " << numCols << endl << endl;
      }
      Combine<Ordinal, Scalar> combiner;
      // qr( A1 )
      combiner.factor_first (A1, tau1.data(), work.data());
      // View of numCols by numCols upper triangle of A1.
      mat_view_type R1 (numCols, numCols, A1.data(), A1.stride(1));
      // qr( [R1; A2] )
      combiner.factor_inner (R1, A2, tau2.data(), work.data());
      // Extract (a deep copy of) the R factor.
      matrix_type R (R1);
      // Zero out everything below the diagonal of R.
      for (Ordinal j = 0; j < numCols; ++j) {
        for (Ordinal i = j+1; i < numCols; ++i) {
          // FIXME (mfh 26 Nov 2019) I'm assuming I can write to the
          // Matrix or MatView on host, outside of Kokkos.  TSQR
          // always assumed this in the past, but if we want to use
          // Kokkos, we'll need to get rid of that assumption.
          R(i,j) = Scalar {};
        }
      }

      // Compute the explicit Q factor, by starting with A2 and
      // (working up the matrix A,) finishing with A1.
      combiner.apply_inner (ApplyType::NoTranspose,
                            numRows, numCols, numCols,
                            A2.data(), A2.stride(1), tau2.data(),
                            Q1.data(), Q1.stride(1),
                            Q2.data(), Q2.stride(1), work.data());
      combiner.apply_first (ApplyType::NoTranspose,
                            A1, tau1.data(),
                            Q1, work.data());
      if (debug) {
        cerr << "Results of first test problem:" << endl;
        cerr << "-- Test matrix A:" << endl;
        printMatrix (cerr, A_copy);
        cerr << endl << "-- Q factor:" << endl;
        printMatrix (cerr, Q);
        cerr << endl << "-- R factor:" << endl;
        printMatrix (cerr, R);
        cerr << endl;
      }
      const results_type results = localVerify (A_copy, Q, R);

      if (debug) {
        cerr << "\\| A - Q*R \\|_F = " << results[0] << endl
             << "\\| I - Q'*Q \\|_F = " << results[1] << endl
             << "\\| A \\|_F = " << results[2] << endl;
      }
      return results;
    }

    void
    verifyCombine (const int numRows,
                   const int numCols,
                   const bool testReal,
                   const bool testComplex,
                   const bool printFieldNames,
                   const bool simulateSequentialTsqr,
                   const bool debug)
    {
      using TSQR::Random::NormalGenerator;
      using std::cerr;
#ifdef HAVE_KOKKOSTSQR_COMPLEX
      using std::complex;
#endif // HAVE_KOKKOSTSQR_COMPLEX
      using std::cout;
      using std::endl;
      using std::pair;
      using std::string;
      using std::vector;

      //
      // We do tests one after another, using the seed from the
      // previous test in the current test, so that the pseudorandom
      // streams used by the tests are independent.
      //

      // Default seed for the next pseudorandom number generator.
      // This will be the same each time, so if you want
      // nondeterministic behavior, you should modify this routine to
      // let you supply the seed values yourself.
      vector<int> iseed(4);
      iseed[0] = 0;
      iseed[1] = 0;
      iseed[2] = 0;
      iseed[3] = 1;

      // Whether to print the field (i.e., column) names for the
      // output data.
      bool doPrintFieldNames = printFieldNames;

      if (! simulateSequentialTsqr) {
        if (testReal) {
          {
            NormalGenerator<int, float> normgenS (iseed);
            const vector<float> resultsS =
              verifyCombineTemplate (normgenS, normgenS, numRows,
                                     numCols, debug);
            // Only print field names (if at all) once per run, for
            // the first data type.
            printResults (string("float"), numRows, numCols,
                          resultsS, doPrintFieldNames);
            // Print field names at most once.
            doPrintFieldNames = false;
            // Fetch the pseudorandom seed from the previous test.
            normgenS.getSeed (iseed);
          }
          {
            NormalGenerator<int, double> normgenD (iseed);
            const vector<double> resultsD =
              verifyCombineTemplate (normgenD, normgenD, numRows,
                                     numCols, debug);
            printResults (string("double"), numRows, numCols,
                          resultsD, doPrintFieldNames);
            doPrintFieldNames = false;
            normgenD.getSeed (iseed);
          }
        }

        if (testComplex)
          {
#ifdef HAVE_KOKKOSTSQR_COMPLEX
            {
              NormalGenerator<int, complex<float> > normgenC (iseed);
              NormalGenerator<int, float> normgenS (iseed);
              const vector<float> resultsC =
                verifyCombineTemplate (normgenC, normgenS, numRows,
                                       numCols, debug);
              printResults (string("complex<float>"), numRows, numCols,
                            resultsC, doPrintFieldNames);
              doPrintFieldNames = false;
              // Even though normgenC and normgenS each updated the
              // random seed independently, for now we just fetch the
              // updated seed from normgenC.  This should still
              // produce reproducible results.
              normgenC.getSeed (iseed);
            }
            {
              NormalGenerator<int, complex<double> > normgenZ (iseed);
              NormalGenerator<int, double> normgenD (iseed);
              const vector<double> resultsZ =
                verifyCombineTemplate (normgenZ, normgenD, numRows,
                                       numCols, debug);
              printResults (string("complex<double>"), numRows, numCols,
                            resultsZ, doPrintFieldNames);
              doPrintFieldNames = false;
              normgenZ.getSeed (iseed);
            }
#else // NOT HAVE_KOKKOSTSQR_COMPLEX
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                       "Trilinos was not built with "
                                       "complex arithmetic support");
#endif // HAVE_KOKKOSTSQR_COMPLEX
          }
      }
      else { // simulateSequentialTsqr
        if (testReal) {
          {
            NormalGenerator<int, float> normgenS (iseed);
            const vector<float> resultsS =
              verifyCombineSeqTemplate (normgenS, normgenS, numRows,
                                        numCols, debug);
            printSimSeqTsqrResults (string("float"), numRows, numCols,
                                    resultsS, doPrintFieldNames);
            doPrintFieldNames = false;
            normgenS.getSeed (iseed);
          }
          {
            NormalGenerator<int, double> normgenD (iseed);
            const vector<double> resultsD =
              verifyCombineSeqTemplate (normgenD, normgenD, numRows,
                                        numCols, debug);
            printSimSeqTsqrResults (string("double"), numRows, numCols,
                                    resultsD, doPrintFieldNames);
            doPrintFieldNames = false;
            normgenD.getSeed (iseed);
          }
        }

        if (testComplex) {
#ifdef HAVE_KOKKOSTSQR_COMPLEX
          {
            NormalGenerator<int, complex<float> > normgenC (iseed);
            NormalGenerator<int, float> normgenS (iseed);
            const vector<float> resultsC =
              verifyCombineSeqTemplate (normgenC, normgenS, numRows,
                                        numCols, debug);
            printSimSeqTsqrResults (string("complex<float>"), numRows, numCols,
                                    resultsC, doPrintFieldNames);
            doPrintFieldNames = false;
            normgenC.getSeed (iseed);
          }
          {
            NormalGenerator<int, complex<double> > normgenZ (iseed);
            NormalGenerator<int, double> normgenD (iseed);
            const vector<double> resultsZ =
              verifyCombineSeqTemplate (normgenZ, normgenD, numRows,
                                        numCols, debug);
            printSimSeqTsqrResults (string("complex<double>"), numRows,
                                    numCols, resultsZ, doPrintFieldNames);
            doPrintFieldNames = false;
            normgenZ.getSeed (iseed);
          }
#else // NOT HAVE_KOKKOSTSQR_COMPLEX
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                     "Trilinos was not built with "
                                     "complex arithmetic support");
#endif // HAVE_KOKKOSTSQR_COMPLEX
        }
      }
    }
  } // namespace Test
} // namespace TSQR

