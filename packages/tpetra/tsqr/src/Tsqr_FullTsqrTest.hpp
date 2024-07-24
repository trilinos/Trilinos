// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_TEST_FULLTSQRTEST_HPP
#define TSQR_TEST_FULLTSQRTEST_HPP

#include "Tsqr.hpp"
#include "Tsqr_NodeTsqrFactory.hpp"
#include "Tsqr_Random_NormalGenerator.hpp"
#include "Tsqr_Random_GlobalMatrix.hpp"
#include "Tsqr_SequentialTsqr.hpp"
#include "Tsqr_TestSetup.hpp"
#include "Tsqr_GlobalVerify.hpp"
#include "Tsqr_TeuchosMessenger.hpp"
#include "Tsqr_TestUtils.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_TypeNameTraits.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>

namespace TSQR {
  namespace Test {

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

    /// \class FullTsqrVerifier
    /// \brief Test (correctness and) accuracy of Tsqr for one Scalar
    ///   type.
    /// \author Mark Hoemmen
    ///
    /// \tparam Scalar Type of each matrix entry.
    ///
    /// This class is meant to be used only by FullTsqrVerifierCaller.
    /// It performs one accuracy test of Tsqr for the given Scalar
    /// type.  An accuracy test is also a correctness test.  This test
    /// computes accuracy bounds for both orthogonality and forward
    /// errors, and if those bounds are exceeded and the
    /// failIfInaccurate option is enabled, the test will throw a
    /// TsqrInaccurate exception.
    ///
    /// The test takes a Teuchos::ParameterList input.  For a
    /// ParameterList with all parameters, their default values, and
    /// documentation, see the relevant class method in
    /// FullTsqrVerifierCaller.
    template<class Scalar>
    class FullTsqrVerifier {
    public:
      using scalar_type = Scalar;
      using ordinal_type = int;
      using node_tsqr_type = NodeTsqr<ordinal_type, scalar_type>;
      using dist_tsqr_type = DistTsqr<ordinal_type, scalar_type>;
      using tsqr_type = Tsqr<ordinal_type, scalar_type>;

    private:
      static Teuchos::RCP<node_tsqr_type>
      getNodeTsqr(
        const Teuchos::RCP<Teuchos::ParameterList>& testParams,
        const bool myRank,
        const bool verbose,
        const std::string inputPrefix)
      {
        using Teuchos::RCP;
        using Teuchos::rcp;
        using Teuchos::rcp_implicit_cast;
        using std::cerr;
        using std::endl;
        using device_type =
          Kokkos::DefaultExecutionSpace::device_type;
        const char cacheSizeHintParamName[] = "Cache Size Hint";
        const std::string prefix = inputPrefix + "  ";

        auto nodeTsqrParams = Teuchos::parameterList("NodeTsqr");

        size_t cacheSizeHint = 0;
        if (testParams->isType<size_t>(cacheSizeHintParamName)) {
          cacheSizeHint =
            testParams->get<size_t>(cacheSizeHintParamName);
          nodeTsqrParams->set(cacheSizeHintParamName, cacheSizeHint);
        }
        else if (testParams->isType<int>(cacheSizeHintParamName)) {
          cacheSizeHint = static_cast<size_t>
           (testParams->get<int>(cacheSizeHintParamName));
          nodeTsqrParams->set(cacheSizeHintParamName, cacheSizeHint);
        }

        std::string nodeTsqrName("Default");
        if (testParams->isType<std::string>("NodeTsqr")) {
          nodeTsqrName = testParams->get<std::string>("NodeTsqr");
        }
        if (myRank == 0 && verbose) {
          cerr << prefix << "getNodeTsqr:" << endl
               << prefix << "  - NodeTsqr: " << nodeTsqrName << endl
               << prefix << "  - Cache Size Hint: " << cacheSizeHint
               << endl;
        }

        RCP<node_tsqr_type> nodeTsqr;
        using node_tsqr_factory_type = TSQR::NodeTsqrFactory<
          scalar_type, ordinal_type, device_type>;
        nodeTsqr = node_tsqr_factory_type::getNodeTsqr(nodeTsqrName);
        TEUCHOS_ASSERT( ! nodeTsqr.is_null() );

        if (myRank == 0 && verbose) {
          using execution_space = device_type::execution_space;
          const std::string spaceName =
            Teuchos::TypeNameTraits<execution_space>::name();
          const std::string myPrefix = prefix + "  * ";

          cerr << myPrefix << "execution_space: " << spaceName << endl
               << myPrefix << "concurrency: "
               << execution_space().concurrency() << endl
               << myPrefix << "Requested NodeTsqr subclass type: "
               << nodeTsqrName << endl
               << myPrefix << "Actual NodeTsqr subclass type: "
               << Teuchos::typeName(*nodeTsqr) << endl;
        }
        return nodeTsqr;
      }

      //! Instantiate and return a (full) Tsqr instance.
      static Teuchos::RCP<tsqr_type>
      getTsqr(const Teuchos::RCP<Teuchos::ParameterList>& testParams,
              const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
              const bool verbose)
      {
        using Teuchos::RCP;
        using Teuchos::rcp;
        using Teuchos::rcp_implicit_cast;
        using std::cerr;
        using std::endl;
        const int myRank = comm->getRank();

        const std::string prefix("  ");

        if (myRank == 0 && verbose) {
          cerr << prefix << "- Set up TSQR::Tsqr instance" << endl;
        }
        auto nodeTsqr =
          getNodeTsqr(testParams, myRank, verbose, prefix);
        auto scalarMess =
          rcp(new TeuchosMessenger<scalar_type>(comm));
        auto scalarMessBase =
          rcp_implicit_cast<MessengerBase<scalar_type>>(scalarMess);
        RCP<dist_tsqr_type> distTsqr(new dist_tsqr_type);
        distTsqr->init(scalarMessBase);

        return rcp(new tsqr_type(nodeTsqr, distTsqr));
      }

    public:
      /// \brief Verify "full" TSQR's accuracy for the Scalar type.
      ///
      /// \param comm [in] Communicator over which to run the test.
      /// \param testParams [in/out] Parameters for the test.  May
      ///   be modified by each test in turn.
      /// \param randomSeed [in/out] On input: the random seed for
      ///   LAPACK's pseudorandom number generator.  On output: the
      ///   updated random seed.
      ///
      /// \return Whether the test passed.
      static bool verify(
        const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
        const Teuchos::RCP<Teuchos::ParameterList>& testParams,
        std::vector<int>& randomSeed)
      {
        using std::cerr;
        using std::cout;
        using std::endl;
        using Teuchos::ParameterList;
        using Teuchos::parameterList;
        using Teuchos::RCP;
        using Teuchos::rcp;
        using Teuchos::rcp_implicit_cast;
        using matrix_type = Matrix<ordinal_type, scalar_type>;
        using mat_view_type = MatView<ordinal_type, scalar_type>;

        bool success = true;

        TEUCHOS_ASSERT( ! comm.is_null() );
        TEUCHOS_ASSERT( ! testParams.is_null() );

        const int myRank = comm->getRank();
        const int numProcs = comm->getSize();
        const bool verbose = testParams->get<bool>("verbose");
        const ordinal_type numRowsLocal =
          testParams->get<ordinal_type>("numRowsLocal");
        const ordinal_type numCols =
          testParams->get<ordinal_type>("numCols");
        //const int numCores = testParams->get<int>("numCores");
        const bool contiguousCacheBlocks =
          testParams->get<bool>("contiguousCacheBlocks");
        const bool testFactorExplicit =
          testParams->get<bool>("testFactorExplicit");
        const bool testRankRevealing =
          testParams->get<bool>("testRankRevealing");

        if (myRank == 0 && verbose) {
          cerr << "Full TSQR test: Scalar="
               << Teuchos::TypeNameTraits<Scalar>::name() << endl
               << "  - Command-line arguments:" << endl
               << "    * numRowsLocal: " << numRowsLocal << endl
               << "    * numCols: " << numCols << endl
               << "    * contiguousCacheBlocks: "
               << (contiguousCacheBlocks ? "true" : "false") << endl
               << "    * testFactorExplicit: "
               << (testFactorExplicit ? "true" : "false") << endl
               << "    * testRankRevealing: "
               << (testRankRevealing ? "true" : "false") << endl
               << "    * verbose: "
               << (verbose ? "true" : "false") << endl;
        }

        RCP<tsqr_type> tsqr = getTsqr(testParams, comm, verbose);
        TEUCHOS_ASSERT( ! tsqr.is_null() );

        // Space for each process's local part of the test problem.
        // A_local, A_copy, and Q_local are distributed matrices, and
        // R is replicated on all processes sharing the communicator.
        matrix_type A_local(numRowsLocal, numCols);
        matrix_type A_copy(numRowsLocal, numCols);
        matrix_type Q_local(numRowsLocal, numCols);
        matrix_type R(numCols, numCols);

        // Start by filling the test problem with zeros.
        deep_copy(A_local, Scalar {});
        deep_copy(A_copy, Scalar {});
        deep_copy(Q_local, Scalar {});
        deep_copy(R, Scalar {});

        // Create some reasonable singular values for the test problem:
        // 1, 1/2, 1/4, 1/8, ...
        using STS = Teuchos::ScalarTraits<scalar_type>;
        using magnitude_type = typename STS::magnitudeType;
        std::vector<magnitude_type> singularValues(numCols);
        using STM = Teuchos::ScalarTraits<magnitude_type>;
        {
          const magnitude_type scalingFactor = STM::one() + STM::one();
          magnitude_type curVal = STM::one();
          for (magnitude_type& singularValue : singularValues) {
            singularValue = curVal;
            curVal = curVal / scalingFactor;
          }
        }

        // Construct a normal(0,1) pseudorandom number generator with
        // the given random seed.
        using TSQR::Random::NormalGenerator;
        using generator_type = NormalGenerator<ordinal_type, scalar_type>;
        generator_type gen(randomSeed);

        // We need a Messenger for Ordinal-type data, so that we can
        // build a global random test matrix.
        auto ordinalMessenger =
          rcp_implicit_cast<MessengerBase<ordinal_type>>
            (rcp(new TeuchosMessenger<ordinal_type>(comm)));

        // We also need a Messenger for Scalar-type data.  The TSQR
        // implementation already constructed one, but it's OK to
        // construct another one; TeuchosMessenger is just a thin
        // wrapper over the Teuchos::Comm object.
        auto scalarMessenger =
          rcp_implicit_cast<MessengerBase<scalar_type>>
            (rcp(new TeuchosMessenger<scalar_type>(comm)));

        if (myRank == 0 && verbose) {
          cerr << "  - Generate test problem" << endl;
        }

        {
          // Generate a global distributed matrix (whose part local to
          // this process is in A_local) with the given singular values.
          // This part has O(P) communication for P MPI processes.
          using TSQR::Random::randomGlobalMatrix;
          mat_view_type A_local_view(A_local.extent(0),
                                     A_local.extent(1),
                                     A_local.data(),
                                     A_local.stride(1));
          const magnitude_type* const singVals = singularValues.data();
          randomGlobalMatrix(&gen, A_local_view, singVals,
                             ordinalMessenger.getRawPtr(),
                             scalarMessenger.getRawPtr());
        }
        // Save the pseudorandom number generator's seed for any later
        // tests.  The generator keeps its own copy of the seed and
        // updates it internally, so we have to ask for its copy.
        gen.getSeed(randomSeed);

        if (myRank == 0 && verbose) {
          cerr << "-- tsqr->wants_device_memory() = "
               << (tsqr->wants_device_memory() ? "true" : "false")
               << endl;
        }

        using IST =
          typename Kokkos::ArithTraits<scalar_type>::val_type;
        using device_matrix_type =
          Kokkos::View<IST**, Kokkos::LayoutLeft>;

        auto A_h = getHostMatrixView(A_local.view());
        auto A_copy_h = getHostMatrixView(A_copy.view());
        auto Q_h = getHostMatrixView(Q_local.view());
        device_matrix_type A_d;
        device_matrix_type A_copy_d;
        device_matrix_type Q_d;
        if (tsqr->wants_device_memory()) {
          A_d = getDeviceMatrixCopy(A_local.view(), "A_d");
          // Don't copy A_copy yet; see below.
          A_copy_d =
            device_matrix_type("A_copy_d", numRowsLocal, numCols);
          Q_d = device_matrix_type("Q_d", numRowsLocal, numCols);
        }

        // If specified in the test parameters, rearrange cache blocks
        // in the copy.  Otherwise, just copy the test problem into
        // A_copy.  The factorization overwrites the input matrix, so
        // we have to make a copy in order to validate the final
        // result.

        if (! contiguousCacheBlocks) {
          if (myRank == 0 && verbose) {
            cerr << "  - Copy A into A_copy" << endl;
          }
          deep_copy(A_copy, A_local);
          if (tsqr->wants_device_memory()) {
            deep_copy(A_copy_d, A_d);
          }
        }
        else {
          if (myRank == 0 && verbose) {
            cerr << "  - Copy A into A_copy via cache_block" << endl;
          }
          if (tsqr->wants_device_memory()) {
            Scalar* A_copy_d_raw =
              reinterpret_cast<Scalar*>(A_copy_d.data());
            const Scalar* A_d_raw =
              reinterpret_cast<const Scalar*>(A_d.data());
            tsqr->cache_block(numRowsLocal, numCols, A_copy_d_raw,
                              A_d_raw, A_d.stride(1));
            deep_copy(A_copy_h, A_copy_d);
          }
          else {
            tsqr->cache_block(numRowsLocal, numCols, A_copy.data(),
                              A_local.data(), A_local.stride(1));
          }
          if (myRank == 0 && verbose) {
            cerr << "  - Finished cache-blocking the test problem"
                 << endl;
          }
        }

        if (testFactorExplicit) {
          if (myRank == 0 && verbose) {
            cerr << "  - Call factorExplicitRaw" << endl;
          }
          try {
            if (tsqr->wants_device_memory()) {
              Scalar* A_raw =
                reinterpret_cast<Scalar*>(A_copy_d.data());
              Scalar* Q_raw = reinterpret_cast<Scalar*>(Q_d.data());
              tsqr->factorExplicitRaw(A_copy_d.extent(0),
                                      A_copy_d.extent(1),
                                      A_raw,
                                      A_copy_d.stride(1),
                                      Q_raw,
                                      Q_d.stride(1),
                                      R.data(), R.stride(1),
                                      contiguousCacheBlocks);
              if (myRank == 0 && verbose) {
                cerr << "  - Finished factorExplicitRaw; now "
                  "deep_copy(Q_h, Q_d)" << endl;
              }
              deep_copy(Q_h, Q_d);
            }
            else {
              Scalar* A_raw = A_copy.data();
              Scalar* Q_raw = Q_local.data();
              tsqr->factorExplicitRaw(A_copy.extent(0),
                                      A_copy.extent(1),
                                      A_raw,
                                      A_copy.stride(1),
                                      Q_raw,
                                      Q_local.stride(1),
                                      R.data(), R.stride(1),
                                      contiguousCacheBlocks);
              if (myRank == 0 && verbose) {
                cerr << "  - Finished factorExplicitRaw" << endl;
              }
            }
          }
          catch (std::exception& e) {
            std::ostringstream os;
            os << "Proc " << myRank << " threw an exception: "
               << e.what() << endl;
            cerr << os.str();
            MPI_Abort(MPI_COMM_WORLD, -1);
          }

          bool found_nonzero_in_R = false;
          for (ordinal_type j = 0; j < numCols; ++j) {
            for (ordinal_type i = 0; i < numCols; ++i) {
              if (R(i,j) != scalar_type {}) {
                found_nonzero_in_R = true;
              }
            }
          }

          if (! found_nonzero_in_R) {
            success = false;
            if (myRank == 0) {
              const std::string prefix
                (verbose ? "  - *** " : "*** ");
              const std::string scalarName =
                Teuchos::TypeNameTraits<scalar_type>::name();
              cerr << prefix << "For Scalar=" << scalarName
                   << ": R factor resulting from factorExplicitRaw "
                   << "is zero." << endl;
            }
          }
        }
        else {
          if (myRank == 0 && verbose) {
            cerr << "  - Call factor" << endl;
          }
          auto factorOutput = [&] () {
            if (tsqr->wants_device_memory()) {
              Scalar* A_raw =
                reinterpret_cast<Scalar*>(A_copy_d.data());
              auto result =
                tsqr->factor(numRowsLocal, numCols,
                             A_raw, A_copy_d.stride(1),
                             R.data(), R.stride(1),
                             contiguousCacheBlocks);
              deep_copy(A_copy_h, A_copy_d);
              return result;
            }
            else {
              Scalar* A_raw =
                reinterpret_cast<Scalar*>(A_copy_d.data());
              return tsqr->factor(numRowsLocal, numCols,
                                  A_raw, A_copy.stride(1),
                                  R.data(), R.stride(1),
                                  contiguousCacheBlocks);
            }
          } ();

          if (myRank == 0 && verbose) {
            cerr << "  - Finished factor; call explicit_Q" << endl;
          }
          if (tsqr->wants_device_memory()) {
            const Scalar* A_raw =
              reinterpret_cast<const Scalar*>(A_copy_d.data());
            Scalar* Q_raw = reinterpret_cast<Scalar*>(Q_d.data());
            tsqr->explicit_Q(numRowsLocal, numCols,
                             A_raw, A_copy_d.stride(1),
                             factorOutput, numCols,
                             Q_raw, Q_d.stride(1),
                             contiguousCacheBlocks);
            deep_copy(Q_h, Q_d);
          }
          else {
            const Scalar* A_raw = A_copy.data();
            Scalar* Q_raw = Q_local.data();
            tsqr->explicit_Q(numRowsLocal, numCols,
                             A_raw, A_copy.stride(1),
                             factorOutput, numCols,
                             Q_raw, Q_local.stride(1),
                             contiguousCacheBlocks);
          }
          if (myRank == 0 && verbose) {
            cerr << "  - Finished explicit_Q" << endl;
          }
        }

        // Optionally, test rank-revealing capability.  We do this
        // before un-cache-blocking the explicit Q factor, since
        // revealRank can work with contiguous cache blocks, and
        // modifies the Q factor if the matrix doesn't have full
        // column rank.
        if (testRankRevealing) {
          // If 2^{# columns} > machine precision, then our choice
          // of singular values will make the smallest singular
          // value < machine precision.  In that case, the SVD can't
          // promise it will distinguish between tiny and zero.  If
          // the number of columns is less than that, we can use a
          // tolerance of zero to test the purported rank with the
          // actual numerical rank.
          const magnitude_type tol = STM::zero();
          if (myRank == 0 && verbose) {
            cerr << "  - Call revealRankRaw" << endl;
          }
          const ordinal_type rank = [&] () {
            Scalar* Q_raw = tsqr->wants_device_memory() ?
              reinterpret_cast<Scalar*>(Q_d.data()) :
              Q_local.data();
            const ordinal_type ldq = tsqr->wants_device_memory() ?
              Q_d.stride(1) : Q_local.stride(1);
            return tsqr->revealRankRaw(numRowsLocal, numCols,
                                       Q_raw, ldq,
                                       R.data(), R.stride(1),
                                       tol, contiguousCacheBlocks);
          } ();
          if (myRank == 0 && verbose) {
            cerr << "  - Finished revealRankRaw" << endl;
          }
          magnitude_type two_to_the_numCols = STM::one();
          for (int k = 0; k < numCols; ++k) {
            const magnitude_type two = STM::one() + STM::one();
            two_to_the_numCols *= two;
          }
          // Throw in a factor of 10, just for more tolerance of
          // rounding error (so the test only fails if something is
          // really broken).
          if (two_to_the_numCols > magnitude_type(10) * STM::eps()) {
            TEUCHOS_TEST_FOR_EXCEPTION
              (rank != numCols, std::logic_error, "The matrix of " <<
               numCols << " columns should have full numerical rank, "
               "but Tsqr reports that it has rank " << rank << ".  "
               "Please report this bug to the Kokkos developers.");
            if (myRank == 0 && verbose) {
              cerr << "  - Tested rank-revealing capability" << endl;
            }
          }
          else {
            if (myRank == 0 && verbose) {
              cerr << "  - Not testing rank-revealing capability; "
                "too many columns" << endl;
            }
          }
        }
        // "Un"-cache-block the output, if contiguous cache blocks
        // were used.  This is only necessary because global_verify()
        // doesn't currently support contiguous cache blocks.
        if (contiguousCacheBlocks) {
          // Use A_copy(_d) as scratch for un-cache-blocking Q_local.
          if (myRank == 0 && verbose) {
            cerr << "  - Call Tsqr::un_cache_block" << endl;
          }
          if (tsqr->wants_device_memory()) {
            Scalar* A_copy_d_raw =
              reinterpret_cast<Scalar*>(A_copy_d.data());
            const Scalar* Q_d_raw =
              reinterpret_cast<const Scalar*>(Q_d.data());
            tsqr->un_cache_block(numRowsLocal, numCols,
                                 A_copy_d_raw,
                                 A_copy_d.stride(1),
                                 Q_d_raw);
            deep_copy(Q_h, A_copy_d);
          }
          else {
            tsqr->un_cache_block(numRowsLocal, numCols,
                                 A_copy.data(),
                                 A_copy.stride(1),
                                 Q_local.data());
            deep_copy(Q_local, A_copy);
          }
          if (myRank == 0 && verbose) {
            cerr << "  - Finished Tsqr::un_cache_block" << endl;
          }
        }
        else {
          if (tsqr->wants_device_memory()) {
            deep_copy(Q_h, Q_d);
          }
        }

        if (myRank == 0 && verbose) {
          cerr << "  - Call global_verify" << endl;
        }
        const auto results =
          global_verify(numRowsLocal, numCols,
                        A_local.data(), A_local.stride(1),
                        Q_local.data(), Q_local.stride(1),
                        R.data(), R.stride(1),
                        scalarMessenger.getRawPtr());
        if (myRank == 0 && verbose) {
          cerr << "  - Finished global_verify" << endl;
        }

        // Print the results on Proc 0.
        if (myRank == 0) {
          if (testParams->get<bool>("printFieldNames")) {
            cout << "%"
                 << "method"
                 << ",scalarType"
                 << ",numRowsLocal"
                 << ",numCols"
                 << ",numProcs"
                 << ",cacheSizeHint"
                 << ",contiguousCacheBlocks"
                 << ",absFrobResid"
                 << ",absFrobOrthog"
                 << ",frobA" << endl;
            // We don't need to print field names again for the other
            // tests, so set the test parameters accordingly.
            testParams->set("printFieldNames", false);
          }
          if (testParams->get<bool>("printResults")) {
            const std::string scalarName =
              Teuchos::TypeNameTraits<scalar_type>::name();
            cout << "Tsqr"
                 << "," << scalarName
                 << "," << numRowsLocal
                 << "," << numCols
                 << "," << numProcs
                 << "," << tsqr->cache_size_hint()
                 << "," << contiguousCacheBlocks
                 << "," << results[0]
                 << "," << results[1]
                 << "," << results[2]
                 << endl;
          }
        }

        // If requested, check accuracy and fail if results are not
        // sufficiently accurate.
        if (testParams->get<bool>("failIfInaccurate")) {
          // Avoid overflow of the local Ordinal type, by casting
          // first to a floating-point type.
          const magnitude_type dimsProd = magnitude_type(numRowsLocal) *
            magnitude_type(numProcs) * magnitude_type(numCols*numCols);

          // Relative residual error is ||A-Q*R|| / ||A||, or just
          // ||A-Q*R|| if ||A|| == 0.  (The result had better be zero
          // in the latter case.)  A reasonable error bound should
          // incorporate the dimensions of the matrix, since this
          // indicates the amount of rounding error.  Square root of
          // the matrix dimensions is an old heuristic from Wilkinson
          // or perhaps even an earlier source.  We include a factor
          // of 10 so that the test won't fail unless there is a
          // really good reason.
          const magnitude_type relResidBound =
            magnitude_type(10) * STM::squareroot(dimsProd) * STM::eps();

          // Orthogonality of the matrix should not depend on the
          // matrix dimensions, if we measure in the 2-norm.
          // However, we are measuring in the Frobenius norm, so
          // it's appropriate to multiply eps by the number of
          // entries in the matrix for which we compute the
          // Frobenius norm.  We include a factor of 10 for the same
          // reason as mentioned above.
          const magnitude_type orthoBound =
            magnitude_type(10*numCols*numCols) * STM::eps();

          // Avoid division by zero.
          const magnitude_type relResidError = results[0] /
            (results[2] == STM::zero() ? STM::one() : results[2]);

          if (relResidError > relResidBound) {
            success = false;
            if (myRank == 0) {
              const std::string prefix
                (verbose ? "  - *** " : "*** ");
              const std::string scalarName =
                Teuchos::TypeNameTraits<scalar_type>::name();
              const std::string relResStr
                (results[2] == STM::zero() ? " / ||A||_F" : "");
              cerr << prefix << "For Scalar=" << scalarName
                   << ": Inaccurate residual ||A - QR||_F"
                   << relResStr
                   << (results[2] == STM::zero() ? " / ||A||_F" : "")
                   << " = " << relResidError << "." << endl
                   << prefix << "It's greater than the bound "
                   << relResidBound << " by a factor of "
                   << relResidError / relResidBound << "." << endl;
            }
          }
          const magnitude_type orthoError = results[1];
          if (orthoError > orthoBound) {
            success = false;
            if (myRank == 0) {
              const std::string prefix
                (verbose ? "  - *** " : "*** ");
              const std::string scalarName =
                Teuchos::TypeNameTraits<scalar_type>::name();
              cerr << prefix << "For Scalar=" << scalarName
                   << ": Inaccurate orthogonality measure "
                   << "||I - Q^* Q||_F = " << orthoError << "."
                   << endl << prefix << "It's greater than the bound "
                   << orthoBound << " by a factor of "
                   << orthoError / orthoBound << "." << endl;
            }
          }
        } // if (the tests should fail on inaccuracy)
        return success;
      }

      /// \brief Benchmark "full" TSQR for the Scalar type.
      ///
      /// \param comm [in] Communicator over which to run the test.
      /// \param testParams [in/out] Parameters for the test.  May
      ///   be modified by each test in turn.
      /// \param randomSeed [in/out] On input: the random seed for
      ///   LAPACK's pseudorandom number generator.  On output: the
      ///   updated random seed.
      static void benchmark(
        const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
        const Teuchos::RCP<Teuchos::ParameterList>& testParams,
        std::vector<int>& randomSeed)
      {
        using std::cerr;
        using std::cout;
        using std::endl;
        using Teuchos::ParameterList;
        using Teuchos::parameterList;
        using Teuchos::RCP;
        using Teuchos::rcp;
        using Teuchos::rcp_implicit_cast;
        using matrix_type = Matrix<ordinal_type, scalar_type>;
        using mat_view_type = MatView<ordinal_type, scalar_type>;

        TEUCHOS_ASSERT( ! comm.is_null() );
        TEUCHOS_ASSERT( ! testParams.is_null() );

        const int myRank = comm->getRank();
        const int numProcs = comm->getSize();
        const bool verbose = testParams->get<bool>("verbose");
        const ordinal_type numRowsLocal =
          testParams->get<ordinal_type>("numRowsLocal");
        const ordinal_type numCols =
          testParams->get<ordinal_type>("numCols");
        const int numTrials = testParams->get<int>("numTrials");
        /* const */ bool contiguousCacheBlocks =
          testParams->get<bool>("contiguousCacheBlocks");
          // 9/2020  contiguousCacheBlocks should be const, but
          // gcc7.2+cuda9 with std=c++14 fails when const is used;
          // see https://github.com/trilinos/Trilinos/pull/8047
        const bool testFactorExplicit =
          testParams->get<bool>("testFactorExplicit");
        const bool testRankRevealing =
          testParams->get<bool>("testRankRevealing");

        if (myRank == 0 && verbose) {
          cerr << "Full TSQR test: Scalar="
               << Teuchos::TypeNameTraits<Scalar>::name() << endl
               << "  - Command-line arguments:" << endl
               << "    * numRowsLocal: " << numRowsLocal << endl
               << "    * numCols: " << numCols << endl
               << "    * numTrials: " << numTrials << endl
               << "    * contiguousCacheBlocks: "
               << (contiguousCacheBlocks ? "true" : "false") << endl
               << "    * testFactorExplicit: "
               << (testFactorExplicit ? "true" : "false") << endl
               << "    * testRankRevealing: "
               << (testRankRevealing ? "true" : "false") << endl
               << "    * verbose: "
               << (verbose ? "true" : "false") << endl;
        }

        RCP<tsqr_type> tsqr = getTsqr(testParams, comm, verbose);
        TEUCHOS_ASSERT( ! tsqr.is_null() );

        // Space for each process's local part of the test problem.
        // A_local, A_copy, and Q_local are distributed matrices, and
        // R is replicated on all processes sharing the communicator.
        matrix_type A_local(numRowsLocal, numCols);
        matrix_type A_copy(numRowsLocal, numCols);
        matrix_type Q_local(numRowsLocal, numCols);
        matrix_type R(numCols, numCols);

        // Start by filling the test problem with zeros.
        deep_copy(A_local, Scalar {});
        deep_copy(A_copy, Scalar {});
        deep_copy(Q_local, Scalar {});
        deep_copy(R, Scalar {});

        // Create some reasonable singular values for the test problem:
        // 1, 1/2, 1/4, 1/8, ...
        using STS = Teuchos::ScalarTraits<scalar_type>;
        using magnitude_type = typename STS::magnitudeType;
        std::vector<magnitude_type> singularValues(numCols);
        using STM = Teuchos::ScalarTraits<magnitude_type>;
        {
          const magnitude_type scalingFactor = STM::one() + STM::one();
          magnitude_type curVal = STM::one();
          for (magnitude_type& singularValue : singularValues) {
            singularValue = curVal;
            curVal = curVal / scalingFactor;
          }
        }

        // Construct a normal(0,1) pseudorandom number generator with
        // the given random seed.
        using TSQR::Random::NormalGenerator;
        using generator_type = NormalGenerator<ordinal_type, scalar_type>;
        generator_type gen(randomSeed);

        // We need a Messenger for Ordinal-type data, so that we can
        // build a global random test matrix.
        auto ordinalMessenger =
          rcp_implicit_cast<MessengerBase<ordinal_type>>
            (rcp(new TeuchosMessenger<ordinal_type>(comm)));

        // We also need a Messenger for Scalar-type data.  The TSQR
        // implementation already constructed one, but it's OK to
        // construct another one; TeuchosMessenger is just a thin
        // wrapper over the Teuchos::Comm object.
        auto scalarMessenger =
          rcp_implicit_cast<MessengerBase<scalar_type>>
            (rcp(new TeuchosMessenger<scalar_type>(comm)));

        if (myRank == 0 && verbose) {
          cerr << "  - Generate test problem" << endl;
        }

        {
          // Generate a global distributed matrix (whose part local to
          // this process is in A_local) with the given singular values.
          // This part has O(P) communication for P MPI processes.
          using TSQR::Random::randomGlobalMatrix;
          mat_view_type A_local_view(A_local.extent(0),
                                     A_local.extent(1),
                                     A_local.data(),
                                     A_local.stride(1));
          const magnitude_type* const singVals = singularValues.data();
          randomGlobalMatrix(&gen, A_local_view, singVals,
                             ordinalMessenger.getRawPtr(),
                             scalarMessenger.getRawPtr());
        }
        // Save the pseudorandom number generator's seed for any later
        // tests.  The generator keeps its own copy of the seed and
        // updates it internally, so we have to ask for its copy.
        gen.getSeed(randomSeed);

        if (myRank == 0 && verbose) {
          cerr << "-- tsqr->wants_device_memory() = "
               << (tsqr->wants_device_memory() ? "true" : "false")
               << endl;
        }

        using IST =
          typename Kokkos::ArithTraits<scalar_type>::val_type;
        using device_matrix_type =
          Kokkos::View<IST**, Kokkos::LayoutLeft>;

        auto A_h = getHostMatrixView(A_local.view());
        auto A_copy_h = getHostMatrixView(A_copy.view());
        auto Q_h = getHostMatrixView(Q_local.view());
        device_matrix_type A_d;
        device_matrix_type A_copy_d;
        device_matrix_type Q_d;
        if (tsqr->wants_device_memory()) {
          A_d = getDeviceMatrixCopy(A_local.view(), "A_d");
          A_copy_d = getDeviceMatrixCopy(A_local.view(), "A_copy_d");
          Q_d = device_matrix_type("Q_d", numRowsLocal, numCols);
        }

        //
        // Time (cache_block, un_cache_block) repeatedly.
        //
        Teuchos::Time cacheBlockTimer("cache_block");
        if (contiguousCacheBlocks) {
          {
            Teuchos::TimeMonitor timeMon(cacheBlockTimer);

            for (int trialNum = 0; trialNum < numTrials; ++trialNum) {
              if (tsqr->wants_device_memory()) {
                // cache_block result goes into A_copy_d.
                Scalar* A_copy_d_raw =
                  reinterpret_cast<Scalar*>(A_copy_d.data());
                Scalar* A_d_raw =
                  reinterpret_cast<Scalar*>(A_d.data());
                tsqr->cache_block(numRowsLocal, numCols, A_copy_d_raw,
                                  A_d_raw, A_d.stride(1));
                // un_cache_block result goes into A_d.
                tsqr->un_cache_block(numRowsLocal, numCols,
                                     A_d_raw,
                                     A_d.stride(1),
                                     A_copy_d_raw);
              }
              else {
                // cache_block result goes into A_copy.
                tsqr->cache_block(numRowsLocal, numCols, A_copy.data(),
                                  A_local.data(), A_local.stride(1));
                // un_cache_block result goes into A_local.
                tsqr->un_cache_block(numRowsLocal, numCols,
                                     A_local.data(),
                                     A_local.stride(1),
                                     A_copy.data());
              }
            } // for each trial
          } // end timing region for cache_block + un_cache_block

          // Finish with an untimed cache_block, so that the code that
          // follows gets cache-blocked input.

          if (tsqr->wants_device_memory()) {
            // cache_block result goes into A_copy_d.
            Scalar* A_copy_d_raw =
              reinterpret_cast<Scalar*>(A_copy_d.data());
            const Scalar* A_d_raw =
              reinterpret_cast<const Scalar*>(A_d.data());
            tsqr->cache_block(numRowsLocal, numCols, A_copy_d_raw,
                              A_d_raw, A_d.stride(1));
          }
          else {
            // cache_block result goes into A_copy.
            tsqr->cache_block(numRowsLocal, numCols, A_copy.data(),
                              A_local.data(), A_local.stride(1));
          }
        } // if contiguousCacheBlocks

        Teuchos::Time fullTsqrTimer("FullTsqr");
        {
          Teuchos::TimeMonitor timeMon(fullTsqrTimer);

          for (int trialNum = 0; trialNum < numTrials; ++trialNum) {
            if (testFactorExplicit) {
              if (tsqr->wants_device_memory()) {
                // factorExplicitRaw result goes into Q_d.
                Scalar* A_raw = contiguousCacheBlocks ?
                  reinterpret_cast<Scalar*>(A_copy_d.data()) :
                  reinterpret_cast<Scalar*>(A_d.data());
                Scalar* Q_raw = reinterpret_cast<Scalar*>(Q_d.data());
                tsqr->factorExplicitRaw(A_copy_d.extent(0),
                                        A_copy_d.extent(1),
                                        A_raw,
                                        A_copy_d.stride(1),
                                        Q_raw,
                                        Q_d.stride(1),
                                        R.data(), R.stride(1),
                                        contiguousCacheBlocks);
              }
              else {
                // factorExplicitRaw result goes into Q_local.
                Scalar* A_raw = contiguousCacheBlocks ?
                  A_copy.data() :
                  A_local.data();
                Scalar* Q_raw = Q_local.data();
                tsqr->factorExplicitRaw(A_copy.extent(0),
                                        A_copy.extent(1),
                                        A_raw,
                                        A_copy.stride(1),
                                        Q_raw,
                                        Q_local.stride(1),
                                        R.data(), R.stride(1),
                                        contiguousCacheBlocks);
              }
            }
            else { // call factor, then explicit_Q
              //
              // factor overwrites its input with part of the result.
              //
              auto factorOutput = [&] () {
                if (tsqr->wants_device_memory()) {
                  Scalar* A_raw = contiguousCacheBlocks ?
                    reinterpret_cast<Scalar*>(A_copy_d.data()) :
                    reinterpret_cast<Scalar*>(A_d.data());
                  auto result =
                    tsqr->factor(numRowsLocal, numCols,
                                 A_raw, A_copy_d.stride(1),
                                 R.data(), R.stride(1),
                                 contiguousCacheBlocks);
                  deep_copy(A_copy_h, A_copy_d);
                  return result;
                }
                else {
                  Scalar* A_raw = contiguousCacheBlocks ?
                    A_copy.data() :
                    A_local.data();
                  return tsqr->factor(numRowsLocal, numCols,
                                      A_raw, A_copy.stride(1),
                                      R.data(), R.stride(1),
                                      contiguousCacheBlocks);
                }
              } ();

              if (tsqr->wants_device_memory()) {
                // explicit_Q result goes into Q_d.
                const Scalar* A_raw = contiguousCacheBlocks ?
                  reinterpret_cast<Scalar*>(A_copy_d.data()) :
                  reinterpret_cast<Scalar*>(A_d.data());
                Scalar* Q_raw = reinterpret_cast<Scalar*>(Q_d.data());
                tsqr->explicit_Q(numRowsLocal, numCols,
                                 A_raw, A_copy_d.stride(1),
                                 factorOutput, numCols,
                                 Q_raw, Q_d.stride(1),
                                 contiguousCacheBlocks);
              }
              else {
                // explicit_Q result goes into Q_local.
                const Scalar* A_raw = contiguousCacheBlocks ?
                  A_copy.data() :
                  A_local.data();
                Scalar* Q_raw = Q_local.data();
                tsqr->explicit_Q(numRowsLocal, numCols,
                                 A_raw, A_copy.stride(1),
                                 factorOutput, numCols,
                                 Q_raw, Q_local.stride(1),
                                 contiguousCacheBlocks);
              }
            } // call factor, then explicit_Q

            // Optionally, test rank-revealing capability.
            // revealRank can work with contiguous cache blocks.
            if (testRankRevealing) {
              const magnitude_type tol = STM::zero();
              Scalar* Q_raw = tsqr->wants_device_memory() ?
                reinterpret_cast<Scalar*>(Q_d.data()) :
                Q_local.data();
              const ordinal_type ldq = tsqr->wants_device_memory() ?
                Q_d.stride(1) : Q_local.stride(1);
              const auto rank =
                tsqr->revealRankRaw(numRowsLocal, numCols,
                                    Q_raw, ldq,
                                    R.data(), R.stride(1),
                                    tol, contiguousCacheBlocks);
              if (myRank == 0 && verbose) {
                cerr << "  - Rank: " << rank << endl;
              }
            }
          } // for each trial
        } // end timing region for "full" TSQR

        // We don't need to un_cache_block the output.  If you have
        // cache-blocked input, then Tsqr assumes that you're mainly
        // concerned about performance for cache-blocked output.
        // We've already timed cache_block and un_cache_block
        // separately above.

        if (myRank == 0) {
          if (testParams->get<bool>("printFieldNames")) {
            cout << "%"
                 << "method"
                 << ",nodeTsqr"
                 << ",scalarType"
                 << ",numRowsLocal"
                 << ",numCols"
                 << ",numTrials"
                 << ",numProcs"
                 << ",cacheSizeHint"
                 << ",contiguousCacheBlocks"
                 << ",testRankRevealing";

            if (contiguousCacheBlocks) {
              cout << ",cacheBlockTime";
            }
            cout << ",fullTsqrTime" << endl;

            // We don't need to print field names again for the other
            // tests, so set the test parameters accordingly.
            testParams->set("printFieldNames", false);
          }
          if (testParams->get<bool>("printResults")) {
            const std::string scalarName =
              Teuchos::TypeNameTraits<scalar_type>::name();
            cout << "FullTsqr"
                 << "," << testParams->get<std::string>("NodeTsqr")
                 << "," << scalarName
                 << "," << numRowsLocal
                 << "," << numCols
                 << "," << numTrials
                 << "," << numProcs
                 << "," << tsqr->cache_size_hint()
                 << "," << (contiguousCacheBlocks ? "true" : "false")
                 << "," << (testRankRevealing ? "true" : "false");

            if (contiguousCacheBlocks) {
              cout << "," << cacheBlockTimer.totalElapsedTime();
            }
            cout << "," << fullTsqrTimer.totalElapsedTime() << endl;
          }
        }
      }
    };

    /// \class FullTsqrVerifierCallerImpl
    /// \brief Implementation of verify and benchmark in
    ///   FullTsqrVerifierCaller.
    /// \author Mark Hoemmen
    template<class ... ScalarTypes>
    class FullTsqrVerifierCallerImpl {
    public:
      static bool verify(
        const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
        const Teuchos::RCP<Teuchos::ParameterList>& testParams,
        std::vector<int>& randomSeed);

      static void benchmark(
        const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
        const Teuchos::RCP<Teuchos::ParameterList>& testParams,
        std::vector<int>& randomSeed);
    };

    // Partial specialization for FirstScalarType, RestScalarTypes...
    template<class FirstScalarType, class ... RestScalarTypes>
    class FullTsqrVerifierCallerImpl<
      FirstScalarType, RestScalarTypes...> {
    public:
      static bool verify(
        const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
        const Teuchos::RCP<Teuchos::ParameterList>& testParams,
        std::vector<int>& randomSeed)
      {
        using first_type = FullTsqrVerifier<FirstScalarType>;
        using rest_type = FullTsqrVerifierCallerImpl<RestScalarTypes...>;
        const bool success1 =
          first_type::verify(comm, testParams, randomSeed);
        const bool success2 =
          rest_type::verify(comm, testParams, randomSeed);
        return success1 && success2;
      }

      static void benchmark(
        const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
        const Teuchos::RCP<Teuchos::ParameterList>& testParams,
        std::vector<int>& randomSeed)
      {
        using first_type = FullTsqrVerifier<FirstScalarType>;
        using rest_type = FullTsqrVerifierCallerImpl<RestScalarTypes...>;
        first_type::benchmark(comm, testParams, randomSeed);
        rest_type::benchmark(comm, testParams, randomSeed);
      }
    };

    // Full specialization for the empty list of Scalar types.
    template<>
    class FullTsqrVerifierCallerImpl<> {
    public:
      static bool verify(
        const Teuchos::RCP<const Teuchos::Comm<int>>& /* comm */,
        const Teuchos::RCP<Teuchos::ParameterList>& /* testParams */,
        std::vector<int>& /* randomSeed */)
      {
        return true;
      }

      static void benchmark(
        const Teuchos::RCP<const Teuchos::Comm<int>>& /* comm */,
        const Teuchos::RCP<Teuchos::ParameterList>& /* testParams */,
        std::vector<int>& /* randomSeed */)
      {}
    };

    /// \class FullTsqrVerifierCaller
    /// \brief Invokes FullTsqrVerifier::verify() over all Scalar types
    ///   in a type list.
    /// \author Mark Hoemmen
    ///
    /// Use this class to test the full TSQR implementation in Tsqr.
    /// It will test Tsqr over a list of Scalar types that you define.
    /// Represent the list as std::tuple.
    class FullTsqrVerifierCaller {
    public:
      /// \typedef ordinal_type
      /// \brief The (local) Ordinal type to use for TSQR.
      using ordinal_type = int;

      /// \brief Return a valid parameter list for verifying Tsqr.
      ///
      /// Call this once to get a valid parameter list with all the
      /// defaults filled in.  This list is valid for all the Scalar
      /// types which TsqrVerifierCaller::run tests.
      Teuchos::RCP<const Teuchos::ParameterList>
      getValidParameterList() const
      {
        using Teuchos::ParameterList;
        using Teuchos::parameterList;
        using Teuchos::RCP;

        RCP<ParameterList> plist = parameterList("FullTsqrVerifier");

        const size_t cacheSizeHint = 0;
        // const int numCores = 1;
        const ordinal_type numRowsLocal = 100;
        const ordinal_type numCols = 10;
        const int numTrials = 100;
        const bool contiguousCacheBlocks = false;
        const bool testFactorExplicit = true;
        const bool testRankRevealing = true;
        const bool printFieldNames = true;
        const bool printResults = true;
        const bool failIfInaccurate = true;
        const std::string nodeTsqr("Default");
        const bool verbose = false;

        // Parameters for configuring Tsqr itself.
        plist->set("Cache Size Hint", cacheSizeHint,
                   "Cache size hint in bytes.  "
                   "Zero means TSQR picks a reasonable default.");

        // Parameters for testing Tsqr.
        plist->set("numRowsLocal", numRowsLocal,
                   "Number of rows per (MPI) process in the test "
                   "matrix.  Must be >= the number of columns.");
        plist->set("numCols", numCols,
                   "Number of columns in the test matrix.");
        plist->set("numTrials", numTrials,
                   "Number of trials; only used when the "
                   "command-line option \"--benchmark\" is set).");
        plist->set("contiguousCacheBlocks", contiguousCacheBlocks,
                   "Whether to test the factorization with "
                   "contiguously stored cache blocks.");
        plist->set("testFactorExplicit", testFactorExplicit,
                   "Whether to test TSQR's factorExplicit() (a "
                   "hopefully faster path than calling factor() and "
                   "explicit_Q() in sequence).");
        plist->set("testRankRevealing", testRankRevealing,
                   "Whether to test TSQR's rank-revealing capability.");
        plist->set("printFieldNames", printFieldNames,
                   "Whether to print field names (this is only done "
                   "once, for all Scalar types tested).");
        plist->set("printResults", printResults,
                   "Whether to print test results.");
        plist->set("failIfInaccurate", failIfInaccurate,
                   "Whether to fail the test if the factorization "
                   "is not sufficiently accurate.");
        plist->set("NodeTsqr", nodeTsqr, "NodeTsqr subclass to use; "
                   "\"Default\" means let TSQR pick it");
        plist->set("verbose", verbose,
                    "Whether to print verbose debugging output.");
        return plist;
      }

      /// \brief Run TsqrVerifier<ScalarType>::verify() for every
      ///   ScalarType in ScalarTypes...
      ///
      /// \tparam ScalarTypes Zero or more Scalar types to test.
      ///
      /// \param testParams [in/out] List of parameters for all tests
      ///   to run.  Call getValidParameterList() to get a valid list
      ///   of parameters with default values and documentation.
      template<class ... ScalarTypes>
      bool verify(
        const Teuchos::RCP<Teuchos::ParameterList>& testParams)
      {
        using impl_type = FullTsqrVerifierCallerImpl<ScalarTypes...>;
        return impl_type::verify(comm_, testParams, randomSeed_);
      }

      /// \brief Run TsqrVerifier<ScalarType>::benchmark() for every
      ///   ScalarType in ScalarTypes...
      ///
      /// \tparam ScalarTypes Zero or more Scalar types to test.
      ///
      /// \param testParams [in/out] List of parameters for all tests
      ///   to run.  Call getValidParameterList() to get a valid list
      ///   of parameters with default values and documentation.
      template<class ... ScalarTypes>
      void benchmark(
        const Teuchos::RCP<Teuchos::ParameterList>& testParams)
      {
        using impl_type = FullTsqrVerifierCallerImpl<ScalarTypes...>;
        impl_type::benchmark(comm_, testParams, randomSeed_);
      }

      /// \brief Full constructor.
      ///
      /// \param comm [in] Communicator (with one or more processes)
      ///   over which to perform tests.
      ///
      /// \param randomSeed [in] The seed for LAPACK's pseudorandom
      ///   number generator.  An array of four integers, satisfying
      ///   the requirements of LAPACK's _LARNV routines.  The array
      ///   elements must be in [0,4095], and the last element
      ///   (iseed[3]) must be odd.  Call \c defaultRandomSeed() for a
      ///   constant default value (if you want the same results each
      ///   time; not "random" but reproducible).
      FullTsqrVerifierCaller(
        const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
        const std::vector<int>& randomSeed)
        : comm_(comm),
          randomSeed_(validateRandomSeed(randomSeed))
      {}

      /// \brief One-argument constructor.
      ///
      /// Fills in defaults for the other arguments that the full
      /// constructor would take.
      ///
      /// \param comm [in] Communicator (with one or more processes)
      ///   over which to perform tests.
      FullTsqrVerifierCaller(
        const Teuchos::RCP<const Teuchos::Comm<int>>& comm)
        : comm_(comm),
          randomSeed_(defaultRandomSeed())
      {}

      //! Validate the given random seed.
      static std::vector<int>
      validateRandomSeed(const std::vector<int>& seed)
      {
        TEUCHOS_TEST_FOR_EXCEPTION
          (seed.size() < 4, std::invalid_argument, "Invalid random "
           "seed: Need an array of four integers, but you gave us "
           << seed.size() << " of them.");
        for (size_t k = 0; k < seed.size(); ++k) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (seed[k] < 0 || seed[k] > 4095, std::invalid_argument,
             "seed[" << k << "]=" << seed[k] << " is invalid.  "
             "Each of the four seeds must be in [0, 4095].");
        }
        TEUCHOS_TEST_FOR_EXCEPTION
          (seed[3] % 2 != 1, std::invalid_argument, "seed[3]="
           << seed[3] << " is invalid: it must be odd.");
        return seed;
      }

      //! Default random seed.
      static std::vector<int> defaultRandomSeed() {
        return {0, 0, 0, 1};
      }

    private:
      /// \brief Communicator over which to perform tests.
      ///
      /// This communicator may include one or more processes.
      /// MPI is not required (it may be a "serial communicator").
      Teuchos::RCP<const Teuchos::Comm<int>> comm_;

      /// \brief The seed for LAPACK's pseudorandom number generator.
      ///
      /// Array of four integers, satisfying the requirements of
      /// LAPACK's _LARNV routines.  The array elements must be in
      /// [0,4095], and the last element (iseed[3]) must be odd.
      std::vector<int> randomSeed_;
    };

  } // namespace Test
} // namespace TSQR

#endif // TSQR_TEST_FULLTSQRTEST_HPP
