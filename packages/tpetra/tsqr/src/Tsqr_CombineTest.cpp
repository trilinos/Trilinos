// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_CombineTest.hpp"

#include "Tsqr_Random_NormalGenerator.hpp"
#include "Tsqr_Random_MatrixGenerator.hpp"

#include "Tsqr_CombineFactory.hpp"
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

    template<class Ordinal, class Scalar>
    void
    fill_with_identity_columns (const MatView<Ordinal, Scalar>& A)
    {
      deep_copy (A, Scalar {});
      const Ordinal numCols = A.extent (1);
      // FIXME (mfh 08 Dec 2019) Eventually stop writing to Matrix or
      // MatView entries on host, for eventual GPU-ization.
      for (Ordinal j = 0; j < numCols; ++j) {
        A(j,j) = Scalar (1.0);
      }
    }

    template<class Ordinal, class MagnitudeType, class NormalGenType>
    void
    generateSingularValues (NormalGenType& magGen,
                            std::vector<MagnitudeType>& sigma,
                            const Ordinal numValues)
    {
      using mag_type = MagnitudeType;
      const mag_type machEps =
        std::numeric_limits<mag_type>::epsilon();
      sigma.resize (numValues);

      // Relative amount by which to perturb each singular value.  The
      // perturbation will be multiplied by a normal(0,1) pseudorandom
      // number drawn from magGen.
      const mag_type perturbationFactor = mag_type(10) * machEps;

      sigma[0] = mag_type (1);
      for (Ordinal k = 1; k < numValues; ++k) {
        const mag_type perturbation = perturbationFactor * magGen();
        const mag_type beforePerturb = sigma[k-1] / mag_type(2);
        const mag_type candidate = beforePerturb + perturbation;

        // If adding the perturbation to beforePerturb would result
        // in a nonpositive number, subtract instead.
        if (candidate <= mag_type {}) {
          sigma[k] = beforePerturb - perturbation;
        }
        else {
          sigma[k] = candidate;
        }
      }
    }

    static void
    printCombineFieldNames ()
    {
      using std::cout;
      using std::endl;

      const char prefix[] = "%";
      cout << prefix << "kernel"
           << ",combiner"
           << ",scalarType"
           << ",numRows"
           << ",numCols"
           << ",frobA"
           << ",absFrobResid"
           << ",absFrobOrthog"
           << endl;
    }

    template<class MagnitudeType>
    static void
    printR1R2results (const std::string& combinerName,
                      const std::string& scalarName,
                      const int numCols,
                      const std::vector<MagnitudeType>& results)
    {
      using std::cout;
      using std::endl;

      cout << "R1R2"
           << "," << combinerName
           << "," << scalarName
           << "," << (2*numCols)
           << "," << numCols
           << "," << results[2]
           << "," << results[0]
           << "," << results[1]
           << endl;
    }

    template<class MagnitudeType>
    static void
    printR3Aresults (const std::string& combinerName,
                     const std::string& scalarName,
                     const int numRows,
                     const int numCols,
                     const std::vector<MagnitudeType>& results)
    {
      using std::cout;
      using std::endl;

      cout << "R3A"
           << "," << combinerName
           << "," << scalarName
           << "," << numRows
           << "," << numCols
           << "," << results[5]
           << "," << results[3]
           << "," << results[4]
           << endl;
    }

    template<class MagnitudeType>
    static void
    printResults (const std::string& combinerName,
                  const std::string& scalarName,
                  const int numRows,
                  const int numCols,
                  const std::vector<MagnitudeType>& results)
    {
      printR1R2results (combinerName, scalarName, numCols, results);
      printR3Aresults (combinerName, scalarName,
                       numRows, numCols, results);
    }

    static void
    printSimSeqTsqrFieldNames ()
    {
      using std::cout;
      using std::endl;

      const char prefix[] = "%";
      cout << prefix
           << "method"
           << ",combiner"
           << ",scalarType"
           << ",numRows"
           << ",numCols"
           << ",frobA"
           << ",absFrobResid"
           << ",absFrobOrthog"
           << endl;
    }

    template<class MagnitudeType>
    static void
    printSimSeqTsqrResults (const std::string& combinerName,
                            const std::string& scalarName,
                            const int numRows,
                            const int numCols,
                            const std::vector<MagnitudeType>& results)
    {
      using std::cout;
      using std::endl;

      cout << "CombineSimSeqTsqr"
           << "," << combinerName
           << "," << scalarName
           << "," << numRows
           << "," << numCols
           << "," << results[2]
           << "," << results[0]
           << "," << results[1]
           << endl;
    }

    template< class MatrixViewType >
    static void
    printMatrix (std::ostream& out,
                 const MatrixViewType& A)
    {
      print_local_matrix (out, A.extent(0), A.extent(1),
                          A.data(), A.stride(1));
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
      return local_verify (A.extent(0), A.extent(1),
                           A.data(), A.stride(1),
                           Q.data(), Q.stride(1),
                           R.data(), R.stride(1));
    }

    /// \brief Test accuracy of TSQR::Combine
    ///
    /// 1. [R1; R2] where R1 and R2 are both ncols by ncols upper
    ///    triangular.
    ///
    /// 2. [R; A] where R is ncols by ncols upper triangular, and A is
    ///    nrows by ncols general dense.
    ///
    /// Print ($\|A - QR\|_F$, $\|I - Q^* Q\|_F$, $\|A\|_F$) for each
    /// test problem (6 numbers in total).
    ///
    template<class Ordinal,
             class Scalar,
             class CombineType>
    void
    verifyCombineTemplate (TSQR::Random::NormalGenerator<Ordinal, Scalar>& gen,
                           TSQR::Random::NormalGenerator<Ordinal, typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& magGen,
                           CombineType& combiner,
                           const std::string& combinerName,
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
      typedef typename STS::magnitudeType mag_type;
      typedef NormalGenerator<Ordinal, Scalar> normgen_type;
      typedef MatrixGenerator<Ordinal, Scalar, normgen_type> matgen_type;
      typedef Matrix<Ordinal, Scalar> matrix_type;
      typedef vector<mag_type> results_type;

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
      vector<mag_type> sigma_R1 (numCols);
      vector<mag_type> sigma_R2 (numCols);
      vector<mag_type> sigma_R3 (numCols);
      vector<mag_type> sigma_A (numCols);
      generateSingularValues (magGen, sigma_R1, numCols);
      generateSingularValues (magGen, sigma_R2, numCols);
      generateSingularValues (magGen, sigma_R3, numCols);
      generateSingularValues (magGen, sigma_A, numCols);

      matrix_type R1 (numCols, numCols, Scalar{});
      matrix_type R2 (numCols, numCols, Scalar{});
      matrix_type R3 (numCols, numCols, Scalar{});
      matrix_type A (numRows, numCols, Scalar{});
      matgen_type matgen (gen);
      matgen.fill_random_R (numCols, R1.data(), R1.stride(1), &sigma_R1[0]);
      matgen.fill_random_R (numCols, R2.data(), R2.stride(1), &sigma_R2[0]);
      matgen.fill_random_R (numCols, R3.data(), R3.stride(1), &sigma_R3[0]);
      matgen.fill_random_svd (numRows, numCols, A.data(), A.stride(1), &sigma_A[0]);

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
      matrix_type Q_R1R2 (Ordinal(2) * numCols, numCols, Scalar {});
      auto Q_R1_Q_R2 = partition_2x1 (Q_R1R2.view (), numCols);
      matrix_type Q_R3A (numCols + numRows, numCols, Scalar {});
      auto Q_R3_A = partition_2x1 (Q_R3A.view (), numCols);

      fill_with_identity_columns (Q_R1R2.view ());
      fill_with_identity_columns (Q_R3A.view ());

      // tau factor arrays, one for each factorization test.
      vector<Scalar> tau_R1R2 (numCols);
      vector<Scalar> tau_R3A (numCols);

      // Workspace array for factorization and applying the Q factor.
      // We recycle this workspace for all tests.
      const Ordinal lwork =
        combiner.work_size (numRows, numCols, numCols);
      vector<Scalar> work (lwork);

      if (debug) {
        cerr << endl << "----------------------------------------" << endl
             << "TSQR::Combine first test problem:" << endl
             << "qr( [R1; R2] ), with R1 and R2 " << numCols
             << " by " << numCols << endl << endl;
      }
      combiner.factor_pair (R1.view (), R2.view (),
                            tau_R1R2.data (), work.data (), lwork);
      combiner.apply_pair (ApplyType ("N"), R2.view (),
                           tau_R1R2.data (),
                           Q_R1_Q_R2.first, Q_R1_Q_R2.second,
                           work.data (), lwork);
      if (debug) {
        cerr << "Results of first test problem:" << endl;
        cerr << "-- Copy of test problem:" << endl;
        print_local_matrix (cerr, A_R1R2.extent (0),
                            A_R1R2.extent (1), A_R1R2.data (),
                            A_R1R2.stride (1));
        cerr << endl << "-- Q factor:" << endl;
        print_local_matrix (cerr, Q_R1R2.extent (0),
                            Q_R1R2.extent (1), Q_R1R2.data (),
                            Q_R1R2.stride (1));
        cerr << endl << "-- R factor:" << endl;
        print_local_matrix (cerr, R1.extent (0), R1.extent (1),
                            R1.data (), R1.stride (1));
        cerr << endl;
      }
      const results_type firstResults =
        local_verify (A_R1R2.extent (0), A_R1R2.extent (1),
                      A_R1R2.data (), A_R1R2.stride (1),
                      Q_R1R2.data (), Q_R1R2.stride (1),
                      R1.data (), R1.stride (1));
      if (debug) {
        cerr << "\\| A - Q*R \\|_F = " << firstResults[0] << endl
             << "\\| I - Q'*Q \\|_F = " << firstResults[1] << endl
             << "\\| A \\|_A = " << firstResults[2] << endl;
        cerr << endl << "----------------------------------------"
             << endl << "TSQR::Combine second test problem:" << endl
             << "qr( [R3; A] ), with R3 " << numCols << " by "
             << numCols << " and A " << numRows << " by " << numCols
             << endl << endl;
      }
      combiner.factor_inner (R3.view (), A.view (),
                             tau_R3A.data (), work.data (), lwork);
      combiner.apply_inner (ApplyType ("N"), A.view (),
                            tau_R3A.data (), Q_R3_A.first,
                            Q_R3_A.second, work.data (), lwork);
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
      vector<mag_type> finalResults;
      finalResults.push_back (firstResults[0]);
      finalResults.push_back (firstResults[1]);
      finalResults.push_back (firstResults[2]);

      finalResults.push_back (secondResults[0]);
      finalResults.push_back (secondResults[1]);
      finalResults.push_back (secondResults[2]);

      const std::string scalarName =
        Teuchos::TypeNameTraits<Scalar>::name ();
      printResults (combinerName, scalarName, numRows, numCols,
                    finalResults);
    }

    template<class Ordinal,
             class Scalar>
    void
    verifyCombineTemplateAllCombiners (std::vector<int>& iseed,
                                       const Ordinal numRows,
                                       const Ordinal numCols,
                                       const bool debug)
    {
      using mag_type =
        typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
      const std::string scalarName =
        Teuchos::TypeNameTraits<Scalar>::name ();

      Random::NormalGenerator<int, Scalar> normgenS (iseed);
      Random::NormalGenerator<int, mag_type> normgenM (iseed);

      using factory_type = CombineFactory<int, Scalar>;
      {
        const std::string combinerName ("Native");
        auto combiner = factory_type::create (combinerName);
        TEUCHOS_ASSERT( combiner.get () != nullptr );
        // Make sure it's the right type.
        using expected_type = CombineNative<int, Scalar>;
        expected_type* combinerPtr =
          dynamic_cast<expected_type*> (combiner.get ());
        TEUCHOS_ASSERT( combinerPtr != nullptr );
        verifyCombineTemplate (normgenS, normgenM, *combiner,
                               combinerName, numRows, numCols,
                               debug);
      }
      {
        const std::string combinerName ("Default");
        auto combiner = factory_type::create (combinerName);
        TEUCHOS_ASSERT( combiner.get () != nullptr );
        // Make sure it's the right type.
        using expected_type = CombineDefault<int, Scalar>;
        expected_type* combinerPtr =
          dynamic_cast<expected_type*> (combiner.get ());
        TEUCHOS_ASSERT( combinerPtr != nullptr );
        verifyCombineTemplate (normgenS, normgenM, *combiner,
                               combinerName, numRows, numCols,
                               debug);
      }

      // Fetch the pseudorandom seed from the previous test.
      //
      // Even though normgenS and normgenM each updated the random
      // seed independently, for now we just fetch the updated seed
      // from normgenS.  This should still produce reproducible
      // results.
      normgenS.getSeed (iseed);
    }

    //! Simulate one combine step of Sequential TSQR
    template<class Ordinal,
             class Scalar,
             class CombineType>
    std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
    verifyCombineSeqTemplate (TSQR::Random::NormalGenerator<Ordinal, Scalar>& gen,
                              TSQR::Random::NormalGenerator<Ordinal, typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& magGen,
                              CombineType& combiner,
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
      typedef typename STS::magnitudeType mag_type;
      typedef NormalGenerator< Ordinal, Scalar > normgen_type;
      typedef MatrixGenerator< Ordinal, Scalar, normgen_type > matgen_type;
      typedef Matrix<Ordinal, Scalar> matrix_type;
      typedef MatView<Ordinal, Scalar> mat_view_type;
      typedef vector<mag_type> results_type;

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
      vector<mag_type> sigma_A1 (numCols);
      vector<mag_type> sigma_A2 (numCols);
      generateSingularValues (magGen, sigma_A1, numCols);
      generateSingularValues (magGen, sigma_A2, numCols);

      // Matrix consisting of two "cache blocks."
      matrix_type A (Ordinal(2)*numRows, numCols, Scalar{});
      auto A1_A2 = partition_2x1 (A, numRows);
      // Views of the two cache blocks.
      mat_view_type A1 = A1_A2.first;
      mat_view_type A2 = A1_A2.second;

      // Fill the two cache blocks with random test problems.
      matgen_type matgen (gen);
      matgen.fill_random_svd (numRows, numCols, A1.data(),
                              A1.stride(1), sigma_A1.data ());
      matgen.fill_random_svd (numRows, numCols, A2.data(),
                              A2.stride(1), sigma_A2.data ());

      // Copy of the resulting test problem, stored as one dense
      // matrix rather than as two blocks.  We will use A_copy to
      // measure the residual error once we've completed the
      // factorization and computed the explicit Q factor.
      matrix_type A_copy (A);

      // Space to put the explicit Q factor.
      matrix_type Q (Ordinal(2) * numRows, numCols, Scalar {});
      fill_with_identity_columns (Q.view ());
      // Two "cache blocks" (as views) of Q.
      auto Q1_Q2 = partition_2x1 (Q.view (), numRows);

      // Two tau factor arrays, one for each cache block.
      vector<Scalar> tau1 (numCols);
      vector<Scalar> tau2 (numCols);

      // Workspace array for factorization and applying the Q factor.
      // We recycle this workspace for all tests.
      const Ordinal lwork =
        combiner.work_size (numRows, numCols, numCols);
      vector<Scalar> work (lwork);

      if (debug) {
        cerr << endl << "----------------------------------------"
          << endl << "TSQR::Combine SequentialTsqr simulation with 2 "
          "cache blocks:" << endl << "qr( [A1; A2] ), with A1 and A2 "
          "A2 each " << numRows << " by " << numCols << endl << endl;
      }
      // qr( A1 )
      combiner.factor_first (A1, tau1.data (), work.data (), lwork);
      // View of numCols by numCols upper triangle of A1.
      mat_view_type R1 (numCols, numCols, A1.data(), A1.stride(1));
      // qr( [R1; A2] )
      combiner.factor_inner (R1, A2, tau2.data (),
                             work.data (), lwork);
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
      combiner.apply_inner (ApplyType::NoTranspose, A2, tau2.data (),
                            Q1_Q2.first, Q1_Q2.second,
                            work.data (), lwork);
      combiner.apply_first (ApplyType::NoTranspose, A1, tau1.data (),
                            Q1_Q2.first, work.data (), lwork);
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
#ifdef HAVE_TPETRATSQR_COMPLEX
      using std::complex;
#endif // HAVE_TPETRATSQR_COMPLEX
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

      if (! simulateSequentialTsqr) {
        printCombineFieldNames ();
        if (testReal) {
          {
            using scalar_type = float;
            verifyCombineTemplateAllCombiners<int, scalar_type>
              (iseed, numRows, numCols, debug);
          }
          {
            using scalar_type = double;
            verifyCombineTemplateAllCombiners<int, scalar_type>
              (iseed, numRows, numCols, debug);
          }
        }
        if (testComplex) {
#ifdef HAVE_TPETRATSQR_COMPLEX
          {
            using scalar_type = std::complex<float>;
            verifyCombineTemplateAllCombiners<int, scalar_type>
              (iseed, numRows, numCols, debug);
          }
          {
            using scalar_type = std::complex<double>;
            verifyCombineTemplateAllCombiners<int, scalar_type>
              (iseed, numRows, numCols, debug);
          }
#else // NOT HAVE_TPETRATSQR_COMPLEX
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::logic_error, "You set testComplex=true, but "
             "Trilinos was not built with complex arithmetic support "
             "enabled.");
#endif // HAVE_TPETRATSQR_COMPLEX
        }
      }
      else { // simulateSequentialTsqr
        printSimSeqTsqrFieldNames ();
        if (testReal) {
          {
            using scalar_type = float;

            NormalGenerator<int, scalar_type> normgenS (iseed);
            auto combiner =
              CombineFactory<int, scalar_type>::create (numCols);
            const std::string combinerName ("?");
            const auto results =
              verifyCombineSeqTemplate (normgenS, normgenS, *combiner,
                                        numRows, numCols, debug);
            const std::string scalarName =
              Teuchos::TypeNameTraits<scalar_type>::name ();
            printSimSeqTsqrResults (combinerName, scalarName,
                                    numRows, numCols, results);
            normgenS.getSeed (iseed);
          }
          {
            using scalar_type = double;

            NormalGenerator<int, scalar_type> normgenS (iseed);
            auto combiner =
              CombineFactory<int, scalar_type>::create (numCols);
            const std::string combinerName ("?");
            const auto results =
              verifyCombineSeqTemplate (normgenS, normgenS, *combiner,
                                        numRows, numCols, debug);
            const std::string scalarName =
              Teuchos::TypeNameTraits<scalar_type>::name ();
            printSimSeqTsqrResults (combinerName, scalarName,
                                    numRows, numCols, results);
            normgenS.getSeed (iseed);
          }
        }

        if (testComplex) {
#ifdef HAVE_TPETRATSQR_COMPLEX
          {
            using scalar_type = complex<float>;
            using mag_type = float;

            NormalGenerator<int, scalar_type> normgenS (iseed);
            NormalGenerator<int, mag_type> normgenM (iseed);
            auto combiner =
              CombineFactory<int, scalar_type>::create (numCols);
            const std::string combinerName ("?");
            const auto results =
              verifyCombineSeqTemplate (normgenS, normgenM, *combiner,
                                        numRows, numCols, debug);
            const std::string scalarName =
              Teuchos::TypeNameTraits<scalar_type>::name ();
            printSimSeqTsqrResults (combinerName, scalarName,
                                    numRows, numCols, results);
            normgenS.getSeed (iseed);
          }
          {
            using scalar_type = complex<double>;
            using mag_type = double;

            NormalGenerator<int, scalar_type> normgenS (iseed);
            NormalGenerator<int, mag_type> normgenM (iseed);
            auto combiner =
              CombineFactory<int, scalar_type>::create (numCols);
            const std::string combinerName ("?");
            const auto results =
              verifyCombineSeqTemplate (normgenS, normgenM, *combiner,
                                        numRows, numCols, debug);
            const std::string scalarName =
              Teuchos::TypeNameTraits<scalar_type>::name ();
            printSimSeqTsqrResults (combinerName, scalarName,
                                    numRows, numCols, results);
            normgenS.getSeed (iseed);
          }
#else // NOT HAVE_TPETRATSQR_COMPLEX
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::logic_error, "Trilinos was not built with "
             "complex arithmetic support.");
#endif // HAVE_TPETRATSQR_COMPLEX
        }
      }
    }
  } // namespace Test
} // namespace TSQR
