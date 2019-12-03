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

#ifndef TSQR_TEST_FULLTSQRTEST_HPP
#define TSQR_TEST_FULLTSQRTEST_HPP

#include "Tsqr.hpp"
#include "Tsqr_NodeTsqrFactory.hpp"
#include "Tsqr_Random_NormalGenerator.hpp"
#include "Tsqr_Random_GlobalMatrix.hpp"
#include "Tsqr_TestSetup.hpp"
#include "Tsqr_GlobalVerify.hpp"
#include "Tsqr_TeuchosMessenger.hpp"
#include "Tsqr_TestUtils.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace TSQR {
  namespace Test {

    /// \class TsqrInaccurate
    /// \brief Signals that a TSQR test failed due to insufficient accuracy.
    class TsqrInaccurate : public std::exception {
    public:
      TsqrInaccurate (const std::string& msg) : msg_ (msg) {}
      const char* what() const throw() override { return msg_.c_str(); }
      ~TsqrInaccurate() throw() override = default;

    private:
      std::string msg_;
    };

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

      //! Instantiate and return a (full) Tsqr instance.
      static Teuchos::RCP<tsqr_type>
      getTsqr (const Teuchos::RCP<Teuchos::ParameterList>& testParams,
               const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
               const bool verbose)
      {
        using Teuchos::ParameterList;
        using Teuchos::rcp_implicit_cast;
        using Teuchos::RCP;
        using Teuchos::rcp;
        using std::endl;
        const char cacheSizeHintParamName[] = "Cache Size Hint";
        const int myRank = comm->getRank ();

        if (myRank == 0 && verbose) {
          std::cerr << "Setting up TSQR::Tsqr instance" << std::endl;
        }
        auto nodeTsqrParams = Teuchos::parameterList ("NodeTsqr");

        if (testParams->isType<size_t> (cacheSizeHintParamName)) {
          const size_t cacheSizeHint =
            testParams->get<size_t> (cacheSizeHintParamName);
          nodeTsqrParams->set (cacheSizeHintParamName, cacheSizeHint);
        }
        else if (testParams->isType<int> (cacheSizeHintParamName)) {
          const size_t cacheSizeHint
            (testParams->get<int> (cacheSizeHintParamName));
          nodeTsqrParams->set (cacheSizeHintParamName, cacheSizeHint);
        }

        //const int numTasks = testParams->get<int> ("numTasks");
        //tsqrParams->set ("Num Tasks", numCores);

        using device_type =
          Kokkos::DefaultExecutionSpace::device_type;
        using node_tsqr_factory_type = TSQR::NodeTsqrFactory<
          scalar_type, ordinal_type, device_type>;
        auto nodeTsqr = node_tsqr_factory_type::getNodeTsqr ();
        TEUCHOS_ASSERT( ! nodeTsqr.is_null () );
        if (myRank == 0 && verbose) {
          using execution_space = device_type::execution_space;
          const std::string spaceName =
            Teuchos::TypeNameTraits<execution_space>::name ();
          std::cerr << "execution_space: " << spaceName << endl
                    << "concurrency: "
                    << execution_space ().concurrency () << endl
                    << "NodeTsqr subclass type: "
                    << Teuchos::typeName (*nodeTsqr) << endl;
        }

        RCP<TeuchosMessenger<scalar_type>> scalarMess =
          rcp (new TeuchosMessenger<scalar_type> (comm));
        RCP<MessengerBase<scalar_type>> scalarMessBase =
          rcp_implicit_cast<MessengerBase<scalar_type>> (scalarMess);
        RCP<dist_tsqr_type> distTsqr = rcp (new dist_tsqr_type);
        distTsqr->init (scalarMessBase);

        return rcp (new tsqr_type (nodeTsqr, distTsqr));
      }

    public:
      /// \brief Run the test for the Scalar type.
      ///
      /// \param comm [in] Communicator over which to run the test.
      /// \param testParams [in/out] Parameters for the test.  May
      ///   be modified by each test in turn.
      /// \param randomSeed [in/out] On input: the random seed for
      ///   LAPACK's pseudorandom number generator.  On output: the
      ///   updated random seed.
      static void
      run (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
           const Teuchos::RCP<Teuchos::ParameterList>& testParams,
           std::vector<int>& randomSeed,
           const bool verbose)
      {
        using std::cerr;
        using std::cout;
        using std::endl;
        using Teuchos::arcp;
        using Teuchos::ParameterList;
        using Teuchos::parameterList;
        using Teuchos::RCP;
        using Teuchos::rcp;
        using Teuchos::rcp_const_cast;
        using Teuchos::rcp_implicit_cast;
        typedef Matrix<ordinal_type, scalar_type> matrix_type;
        typedef MatView<ordinal_type, scalar_type> mat_view_type;
        typedef typename tsqr_type::FactorOutput factor_output_type;

        const int myRank = Teuchos::rank (*comm);
        const int numProcs = Teuchos::size (*comm);

        // Construct TSQR implementation instance.
        RCP<tsqr_type> tsqr = getTsqr (testParams, comm, verbose);
        TEUCHOS_ASSERT( ! tsqr.is_null () );

        // Fetch test parameters from the input parameter list.
        const ordinal_type numRowsLocal =
          testParams->get<ordinal_type> ("numRowsLocal");
        const ordinal_type numCols =
          testParams->get<ordinal_type> ("numCols");
        //const int numCores = testParams->get<int> ("numCores");
        const bool contiguousCacheBlocks =
          testParams->get<bool> ("contiguousCacheBlocks");
        const bool testFactorExplicit =
          testParams->get<bool> ("testFactorExplicit");
        const bool testRankRevealing =
          testParams->get<bool> ("testRankRevealing");
        const bool debug = testParams->get<bool> ("debug");

        if (debug) {
          comm->barrier ();
          if (myRank == 0) {
            cerr << "Full TSQR test command-line arguments:" << endl
                 << "  numRowsLocal: " << numRowsLocal << endl
                 << "  numCols: " << numCols << endl
              // << "  numCores: " << numCores << endl
                 << "  contiguousCacheBlocks: "
                 << (contiguousCacheBlocks ? "true" : "false") << endl
                 << "  testFactorExplicit: "
                 << (testFactorExplicit ? "true" : "false") << endl
                 << "  testRankRevealing: "
                 << (testRankRevealing ? "true" : "false") << endl
                 << "  debug: "
                 << (debug ? "true" : "false") << endl;
          }
        }

        // Space for each process's local part of the test problem.
        // A_local, A_copy, and Q_local are distributed matrices, and
        // R is replicated on all processes sharing the communicator.
        matrix_type A_local (numRowsLocal, numCols);
        matrix_type A_copy (numRowsLocal, numCols);
        matrix_type Q_local (numRowsLocal, numCols);
        matrix_type R (numCols, numCols);

        // Start out by filling the test problem with zeros.
        deep_copy (A_local, Scalar {});
        deep_copy (A_copy, Scalar {});
        deep_copy (Q_local, Scalar {});
        deep_copy (R, Scalar {});

        // Create some reasonable singular values for the test problem:
        // 1, 1/2, 1/4, 1/8, ...
        using STS = Teuchos::ScalarTraits<scalar_type>;
        using magnitude_type = typename STS::magnitudeType;
        std::vector<magnitude_type> singularValues (numCols);
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
        generator_type gen (randomSeed);

        // We need a Messenger for Ordinal-type data, so that we can
        // build a global random test matrix.
        RCP<MessengerBase<ordinal_type>> ordinalMessenger =
          rcp_implicit_cast<MessengerBase<ordinal_type>> (rcp (new TeuchosMessenger<ordinal_type> (comm)));

        // We also need a Messenger for Scalar-type data.  The TSQR
        // implementation already constructed one, but it's OK to
        // construct another one; TeuchosMessenger is just a thin
        // wrapper over the Teuchos::Comm object.
        RCP<MessengerBase<scalar_type>> scalarMessenger =
          rcp_implicit_cast<MessengerBase<scalar_type>> (rcp (new TeuchosMessenger<scalar_type> (comm)));

        if (debug) {
          comm->barrier ();
          if (myRank == 0) {
            cerr << "Generate test problem" << endl;
          }
        }

        {
          // Generate a global distributed matrix (whose part local to
          // this process is in A_local) with the given singular values.
          // This part has O(P) communication for P MPI processes.
          using TSQR::Random::randomGlobalMatrix;
          mat_view_type A_local_view (A_local.extent(0),
                                      A_local.extent(1),
                                      A_local.data(), A_local.stride(1));
          const magnitude_type* const singVals = singularValues.data();
          randomGlobalMatrix<mat_view_type, generator_type> (&gen, A_local_view, singVals,
                                                             ordinalMessenger.getRawPtr(),
                                                             scalarMessenger.getRawPtr());
        }
        // Save the pseudorandom number generator's seed for any later
        // tests.  The generator keeps its own copy of the seed and
        // updates it internally, so we have to ask for its copy.
        gen.getSeed (randomSeed);

        // If specified in the test parameters, rearrange cache blocks
        // in the copy.  Otherwise, just copy the test problem into
        // A_copy.  The factorization overwrites the input matrix, so
        // we have to make a copy in order to validate the final
        // result.
        if (contiguousCacheBlocks) {
          if (debug) {
            comm->barrier ();
            if (myRank == 0) {
              cerr << "Cache-block the test problem" << endl;
            }
          }
          tsqr->cache_block (numRowsLocal, numCols, A_copy.data(),
                             A_local.data(), A_local.stride(1));
          if (debug) {
            comm->barrier ();
            if (myRank == 0) {
              cerr << "Finished cache-blocking the test problem"
                   << endl;
            }
          }
        }
        else {
          if (debug) {
            comm->barrier ();
            if (myRank == 0) {
              cerr << "Copy the test problem (no cache blocking)"
                   << endl;
            }
          }
          deep_copy (A_copy, A_local);
        }

        if (testFactorExplicit) {
          if (debug) {
            comm->barrier ();
            if (myRank == 0) {
              cerr << "Call factorExplicitRaw" << endl;
            }
          }
          tsqr->factorExplicitRaw (A_copy.extent (0), A_copy.extent (1),
                                   A_copy.data (), A_copy.stride (1),
                                   Q_local.data (), Q_local.stride (1),
                                   R.data (), R.stride (1),
                                   contiguousCacheBlocks);
          if (debug) {
            comm->barrier ();
            if (myRank == 0) {
              cerr << "Finished factorExplicitRaw" << endl;
            }
          }
        }
        else {
          if (debug) {
            comm->barrier ();
            if (myRank == 0) {
              cerr << "Call factor" << endl;
            }
          }
          factor_output_type factorOutput =
            tsqr->factor (numRowsLocal, numCols, A_copy.data(),
                          A_copy.stride(1), R.data(), R.stride(1),
                          contiguousCacheBlocks);
          if (debug) {
            comm->barrier ();
            if (myRank == 0) {
              cerr << "Finished factor; call explicit_Q" << endl;
            }
          }
          // Compute the explicit Q factor in Q_local.
          tsqr->explicit_Q (numRowsLocal, numCols, A_copy.data(),
                            A_copy.stride(1), factorOutput, numCols,
                            Q_local.data(), Q_local.stride(1),
                            contiguousCacheBlocks);
          if (debug) {
            comm->barrier ();
            if (myRank == 0) {
              cerr << "Finished explicit_Q" << endl;
            }
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
          const ordinal_type rank =
            tsqr->revealRankRaw (Q_local.extent (0), Q_local.extent (1),
                                 Q_local.data (), Q_local.stride (1),
                                 R.data (), R.stride (1), tol,
                                 contiguousCacheBlocks);

          magnitude_type two_to_the_numCols = STM::one();
          for (int k = 0; k < numCols; ++k) {
            const magnitude_type two = STM::one() + STM::one();
            two_to_the_numCols *= two;
          }
          // Throw in a factor of 10, just for more tolerance of
          // rounding error (so the test only fails if something is
          // really broken).
          if (two_to_the_numCols > magnitude_type(10) * STM::eps ()) {
            TEUCHOS_TEST_FOR_EXCEPTION(
              rank != numCols, std::logic_error, "The matrix of " << numCols
              << " columns should have full numerical rank, but Tsqr reports "
              "that it has rank " << rank << ".  Please report this bug to "
              "the Kokkos developers.");
            if (debug) {
              Teuchos::barrier (*comm);
              if (myRank == 0)
                cerr << "-- Tested rank-revealing capability" << endl;
            }
          }
          else {
            if (debug) {
              Teuchos::barrier (*comm);
              if (myRank == 0)
                cerr << "-- Not testing rank-revealing capability; too many columns" << endl;
            }
          }
        }
        // "Un"-cache-block the output, if contiguous cache blocks
        // were used.  This is only necessary because global_verify()
        // doesn't currently support contiguous cache blocks.
        if (contiguousCacheBlocks) {
          // We can use A_copy as scratch space for
          // un-cache-blocking Q_local, since we're done using
          // A_copy for other things.
          tsqr->un_cache_block (numRowsLocal, numCols, A_copy.data(),
                                A_copy.stride(1), Q_local.data());
          // Overwrite Q_local with the un-cache-blocked Q factor.
          deep_copy (Q_local, A_copy);
          if (debug) {
            Teuchos::barrier (*comm);
            if (myRank == 0)
              cerr << "-- Finished Tsqr::un_cache_block" << endl;
          }
        }

        // Test accuracy of the factorization.
        const std::vector<magnitude_type> results =
          global_verify (numRowsLocal, numCols, A_local.data(), A_local.stride(1),
                         Q_local.data(), Q_local.stride(1), R.data(), R.stride(1),
                         scalarMessenger.getRawPtr());
        if (debug) {
          Teuchos::barrier (*comm);
          if (myRank == 0)
            cerr << "-- Finished global_verify" << endl;
        }

        // Print the results on Proc 0.
        if (myRank == 0) {
          if (testParams->get<bool> ("printFieldNames")) {
            cout << "%"
                 << "method"
                 << ",scalarType"
                 << ",numRowsLocal"
                 << ",numCols"
                 << ",numProcs"
              // << ",numCores"
                 << ",cacheSizeHint"
                 << ",contiguousCacheBlocks"
                 << ",absFrobResid"
                 << ",absFrobOrthog"
                 << ",frobA" << endl;
            // We don't need to print field names again for the other
            // tests, so set the test parameters accordingly.
            testParams->set ("printFieldNames", false);
          }
          if (testParams->get<bool> ("printResults")) {
            cout << "Tsqr"
                 << "," << Teuchos::TypeNameTraits<scalar_type>::name()
                 << "," << numRowsLocal
                 << "," << numCols
                 << "," << numProcs
              // << "," << numCores
                 << "," << tsqr->cache_size_hint()
                 << "," << contiguousCacheBlocks
                 << "," << results[0]
                 << "," << results[1]
                 << "," << results[2]
                 << endl;
          }
        } // if (myRank == 0)

        // If requested, check accuracy and fail if results are not
        // sufficiently accurate.
        if (testParams->get<bool> ("failIfInaccurate")) {
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
          const magnitude_type relResidError =
            results[0] / (results[2] == STM::zero() ? STM::one() : results[2]);
          TEUCHOS_TEST_FOR_EXCEPTION(
            relResidError > relResidBound, TsqrInaccurate, "Full Tsqr "
            "has an inaccurate relative residual ||A - QR||_F"
            << (results[2] == STM::zero() ? " / ||A||_F" : "")
            << " = " << relResidError << ", which is greater than the bound "
            << relResidBound << " by a factor of "
            << relResidError / relResidBound << ".");
          const magnitude_type orthoError = results[1];
          TEUCHOS_TEST_FOR_EXCEPTION(
            orthoError > orthoBound, TsqrInaccurate,
            "Full Tsqr has an inaccurate orthogonality measure ||I - Q^* Q||_F"
            << results[1] << " = " << orthoError << ", which is greater than "
            "the bound " << orthoBound << " by a factor of "
            << orthoError / orthoBound << ".");
        } // if (the tests should fail on inaccuracy)
      }
    };

    /// \class FullTsqrVerifierCallerImpl
    /// \brief This class implements a "function template specialization."
    /// \author Mark Hoemmen
    ///
    /// We want to make FullTsqrVerifierCaller::run() a template
    /// function, with a partial specialization for Cons<CarType,
    /// CdrType> and a full specialization for NullType.  However,
    /// function templates can't have partial specializations, at
    /// least not in the version of the C++ standard currently
    /// supported by Trilinos.  Thus, I've taken the advice of Herb
    /// Sutter (C/C++ Users Journal, 19(7), July 2001), which can be
    /// read online here:
    ///
    /// http://www.gotw.ca/publications/mill17.htm
    ///
    /// Namely, I've implemented the function template via a class
    /// template.  This class is an implementation detail and not
    /// meant to be used anywhere else other than in
    /// FullTsqrVerifierCaller::run().
    template<class TypeListType>
    class FullTsqrVerifierCallerImpl {
    public:
      static void
      run (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
           const Teuchos::RCP<Teuchos::ParameterList>& testParams,
           std::vector<int>& randomSeed,
           const bool verbose);
    };

    //
    // Partial specialization for Cons<CarType, CdrType>.
    //
    template<class CarType, class CdrType>
    class FullTsqrVerifierCallerImpl<TSQR::Test::Cons<CarType, CdrType> > {
    public:
      static void
      run (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
           const Teuchos::RCP<Teuchos::ParameterList>& testParams,
           std::vector<int>& randomSeed,
           const bool verbose)
      {
        using car_type = CarType;
        using cdr_type = CdrType;
        FullTsqrVerifier<car_type>::run (comm, testParams,
                                         randomSeed, verbose);
        FullTsqrVerifierCallerImpl<cdr_type>::run (comm, testParams,
                                                   randomSeed, verbose);
      }
    };

    //
    // Full specialization for NullCons.
    //
    template<>
    class FullTsqrVerifierCallerImpl<TSQR::Test::NullCons> {
    public:
      static void
      run (const Teuchos::RCP<const Teuchos::Comm<int> >&,
           const Teuchos::RCP<Teuchos::ParameterList>&,
           std::vector<int>&,
           const bool /* verbose */)
      {
        // We're at the end of the type list, so do nothing.
      }
    };

    /// \class FullTsqrVerifierCaller
    /// \brief Invokes FullTsqrVerifier::run() over all Scalar types in a type list.
    /// \author Mark Hoemmen
    ///
    /// Use this class to test the full TSQR implementation in Tsqr.
    /// It will test Tsqr over a list of Scalar types that you define,
    /// using \c Cons and \c NullCons.
    class FullTsqrVerifierCaller {
    public:
      /// \typedef ordinal_type
      /// \brief The (local) Ordinal type to use for TSQR.
      typedef int ordinal_type;

      /// \brief Return a valid parameter list for verifying Tsqr.
      ///
      /// Call this once to get a valid parameter list with all the
      /// defaults filled in.  This list is valid for all the Scalar
      /// types which TsqrVerifierCaller::run tests.
      Teuchos::RCP<const Teuchos::ParameterList>
      getValidParameterList () const
      {
        using Teuchos::ParameterList;
        using Teuchos::parameterList;
        using Teuchos::RCP;

        RCP<ParameterList> plist = parameterList ("FullTsqrVerifier");

        const size_t cacheSizeHint = 0;
        // const int numCores = 1;
        const ordinal_type numRowsLocal = 100;
        const ordinal_type numCols = 10;
        const bool contiguousCacheBlocks = false;
        const bool testFactorExplicit = true;
        const bool testRankRevealing = true;
        const bool printFieldNames = true;
        const bool printResults = true;
        const bool failIfInaccurate = true;
        const bool debug = false;

        // Parameters for configuring Tsqr itself.
        plist->set ("Cache Size Hint", cacheSizeHint,
                    "Cache size hint in bytes.  "
                    "Zero means TSQR picks a reasonable default.");
        // plist->set ("Num Tasks", numCores,
        //             "Number of partition(s) to use for TbbTsqr (if "
        //             "applicable).  Must be a positive integer.");

        // Parameters for testing Tsqr.
        plist->set ("numRowsLocal", numRowsLocal,
                    "Number of rows per (MPI) process in the test "
                    "matrix.  Must be >= the number of columns.");
        plist->set ("numCols", numCols,
                    "Number of columns in the test matrix.");
        plist->set ("contiguousCacheBlocks", contiguousCacheBlocks,
                    "Whether to test the factorization with "
                    "contiguously stored cache blocks.");
        plist->set ("testFactorExplicit", testFactorExplicit,
                    "Whether to test TSQR's factorExplicit() (a "
                    "hopefully faster path than calling factor() and "
                    "explicit_Q() in sequence).");
        plist->set ("testRankRevealing", testRankRevealing,
                    "Whether to test TSQR's rank-revealing capability.");
        plist->set ("printFieldNames", printFieldNames,
                    "Whether to print field names (this is only done "
                    "once, for all Scalar types tested).");
        plist->set ("printResults", printResults,
                    "Whether to print test results.");
        plist->set ("failIfInaccurate", failIfInaccurate,
                    "Whether to fail the test if the factorization "
                    "is not sufficiently accurate.");
        plist->set ("debug", debug,
                    "Whether to print debugging output.");
        return plist;
      }

      /// \brief Run TsqrVerifier<T>::run() for every type in the type
      ///   list.
      ///
      /// TypeListType should be either a NullCons (representing an
      /// empty type list, in which case this function does nothing),
      /// or a Cons (whose CarType is a Scalar type to test, and whose
      /// CdrType is either a NullCons or a Cons).
      ///
      /// \param testParams [in/out] List of parameters for all tests
      ///   to run.  Call getValidParameterList() to get a valid list
      ///   of parameters with default values and documentation.
      ///
      template<class TypeListType>
      void
      run (const Teuchos::RCP<Teuchos::ParameterList>& testParams,
           const bool verbose)
      {
        // Using a class with a static method is a way to implement
        // "partial specialization of function templates" (which by
        // itself is not allowed in C++).
        using impl_type = FullTsqrVerifierCallerImpl<TypeListType>;
        impl_type::run (comm_, testParams, randomSeed_, verbose);
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
      FullTsqrVerifierCaller (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                              const std::vector<int>& randomSeed) :
        comm_ (comm),
        randomSeed_ (validateRandomSeed (randomSeed))
      {}

      /// \brief One-argument constructor.
      ///
      /// Fills in defaults for the other arguments that the full
      /// constructor would take.
      ///
      /// \param comm [in] Communicator (with one or more processes)
      ///   over which to perform tests.
      FullTsqrVerifierCaller (const Teuchos::RCP<const Teuchos::Comm<int> >& comm) :
        comm_ (comm),
        randomSeed_ (defaultRandomSeed ())
      {}

      //! Validate the given random seed.
      static std::vector<int>
      validateRandomSeed (const std::vector<int>& seed)
      {
        TEUCHOS_TEST_FOR_EXCEPTION
          (seed.size () < 4, std::invalid_argument, "Invalid random "
           "seed: Need an array of four integers, but you gave us "
           << seed.size () << " of them.");
        for (size_t k = 0; k < seed.size (); ++k) {
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
      static std::vector<int>
      defaultRandomSeed ()
      {
        std::vector<int> seed (4);
        seed[0] = 0;
        seed[1] = 0;
        seed[2] = 0;
        seed[3] = 1;
        return seed;
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

