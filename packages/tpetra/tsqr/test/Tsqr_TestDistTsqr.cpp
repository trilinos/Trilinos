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

#include "Tsqr_Random_NormalGenerator.hpp"
#include "Tsqr_generateStack.hpp"
#include "Tsqr_DistTsqr.hpp"
#include "Tsqr_GlobalTimeStats.hpp"
#include "Tsqr_GlobalVerify.hpp"
#include "Tsqr_printGlobalMatrix.hpp"

#include "Tsqr_Test_MpiAndKokkosScope.hpp"
#include "Tsqr_TeuchosMessenger.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include <algorithm>
#ifdef HAVE_TPETRATSQR_COMPLEX
#  include <complex>
#endif // HAVE_TPETRATSQR_COMPLEX
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace TSQR {
  namespace Test {
    /// \class DistTsqrVerifier
    /// \brief Generic version of \c DistTsqr accuracy test.
    template<class Ordinal, class Scalar>
    class DistTsqrVerifier {
      TSQR::Random::NormalGenerator<Ordinal, Scalar> gen_;
      Teuchos::RCP<MessengerBase<Ordinal> > const ordinalComm_;
      Teuchos::RCP<MessengerBase<Scalar> > const scalarComm_;
      std::string scalarTypeName_;
      std::ostream& out_;
      std::ostream& err_;
      const bool testFactorExplicit_, testFactorImplicit_;
      const bool humanReadable_, printMatrices_, debug_;

    public:
      using ordinal_type = Ordinal;
      using scalar_type = Scalar;
      using mag_type =
        typename Teuchos::ScalarTraits<scalar_type>::magnitudeType;
      using result_type = std::vector<mag_type>;

      /// \brief Constructor, with custom seed value
      ///
      /// \param scalarComm [in/out] Communicator object over which to
      ///   test.
      /// \param seed [in] 4-element vector; the random seed input of
      ///   TSQR::Random::NormalGenerator (which see, since there are
      ///   restrictions on the set of valid seeds)
      /// \param scalarTypeName [in] Human-readable name of the Scalar
      ///   template type parameter
      /// \param out [out] Output stream to which to write results
      /// \param err [out] Output stream to which to write any
      ///   debugging outputs (if applicable) or errors
      /// \param testFactorExplicit [in] Whether to test
      ///   DistTsqr::factorExplicit()
      /// \param testFactorImplicit [in] Whether to test
      ///   DistTsqr::factor() and DistTsqr::explicit_Q()
      /// \param humanReadable [in] Whether printed results should be
      ///   easy for humans to read (vs. easy for parsers to parse)
      /// \param debug [in] Whether to write verbose debug output to
      ///   err
      DistTsqrVerifier(const Teuchos::RCP<MessengerBase<Ordinal> >& ordinalComm,
                       const Teuchos::RCP<MessengerBase<Scalar> >& scalarComm,
                       const std::vector<int>& seed,
                       const std::string& scalarTypeName,
                       std::ostream& out,
                       std::ostream& err,
                       const bool testFactorExplicit,
                       const bool testFactorImplicit,
                       const bool humanReadable,
                       const bool printMatrices,
                       const bool debug) :
        gen_(seed),
        ordinalComm_(ordinalComm),
        scalarComm_(scalarComm),
        scalarTypeName_(scalarTypeName),
        out_(out),
        err_(err),
        testFactorExplicit_(testFactorExplicit),
        testFactorImplicit_(testFactorImplicit),
        humanReadable_(humanReadable),
        printMatrices_(printMatrices),
        debug_(debug)
      {}

      /// \brief Constructor, with default seed value
      ///
      /// This constructor sets a default seed (for the pseudorandom
      /// number generator), which is the same seed (0,0,0,1) each
      /// time.
      ///
      /// \param scalarComm [in/out] Communicator object over which to
      ///   test.
      /// \param scalarTypeName [in] Human-readable name of the Scalar
      ///   template type parameter
      /// \param out [out] Output stream to which to write results
      /// \param err [out] Output stream to which to write any
      ///   debugging outputs (if applicable) or errors
      /// \param testFactorExplicit [in] Whether to test
      ///   DistTsqr::factorExplicit()
      /// \param testFactorImplicit [in] Whether to test
      ///   DistTsqr::factor() and DistTsqr::explicit_Q()
      /// \param humanReadable [in] Whether printed results should be
      ///   easy for humans to read (vs. easy for parsers to parse)
      /// \param debug [in] Whether to write verbose debug output to
      ///   err
      DistTsqrVerifier(const Teuchos::RCP<MessengerBase<Ordinal> >& ordinalComm,
                       const Teuchos::RCP<MessengerBase<Scalar> >& scalarComm,
                       const std::string& scalarTypeName,
                       std::ostream& out,
                       std::ostream& err,
                       const bool testFactorExplicit,
                       const bool testFactorImplicit,
                       const bool humanReadable,
                       const bool printMatrices,
                       const bool debug) :
        ordinalComm_(ordinalComm),
        scalarComm_(scalarComm),
        scalarTypeName_(scalarTypeName),
        out_(out),
        err_(err),
        testFactorExplicit_(testFactorExplicit),
        testFactorImplicit_(testFactorImplicit),
        humanReadable_(humanReadable),
        printMatrices_(printMatrices),
        debug_(debug)
      {}

      /// \brief Get seed vector for pseudorandom number generator
      ///
      /// Fill seed (changing size of vector as necessary) with the
      /// seed vector used by the pseudorandom number generator.  You
      /// can use this to resume the pseudorandom number stream from
      /// where you last were.
      void
      getSeed(std::vector<int>& seed) const
      {
        gen_.getSeed(seed);
      }

      /// \brief Run the DistTsqr accuracy test
      ///
      /// \param numCols [in] Number of columns in the matrix to test.
      ///   Number of rows := (# MPI processors) * ncols.
      void
      verify(const Ordinal numCols,
             const std::string& additionalFieldNames,
             const std::string& additionalData,
             const bool printFieldNames)
      {
        using std::endl;

        const int myRank = scalarComm_->rank();
        if(debug_) {
          scalarComm_->barrier();
          if(myRank == 0) {
            err_ << "Verifying DistTsqr:" << endl;
          }
          scalarComm_->barrier();
        }

        // Generate test problem.
        Matrix<Ordinal, Scalar> A_local, Q_local, R;
        testProblem(A_local, Q_local, R, numCols);
        if(debug_) {
          scalarComm_->barrier();
          if(myRank == 0) {
            err_ << "-- Generated test problem." << endl;
          }
          scalarComm_->barrier();
        }

        // Set up TSQR implementation.
        DistTsqr<Ordinal, Scalar> par;
        par.init (scalarComm_);
        if(debug_) {
          scalarComm_->barrier();
          if(myRank == 0) {
            err_ << "-- DistTsqr object initialized" << endl << endl;
          }
        }

        // Whether we've printed field names (i.e., column headers)
        // yet.  Only matters for non-humanReadable output.
        bool printedFieldNames = false;

        // Test DistTsqr::factor() and DistTsqr::explicit_Q().
        if(testFactorImplicit_) {
          // Factor the matrix A (copied into R, which will be
          // overwritten on output)
          typedef typename DistTsqr<Ordinal, Scalar>::FactorOutput
            factor_output_type;
          factor_output_type factorOutput = par.factor (R.view());
          if(debug_) {
            scalarComm_->barrier();
            if(myRank == 0) {
              err_ << "-- Finished DistTsqr::factor" << endl;
            }
          }
          // Compute the explicit Q factor
          par.explicit_Q(numCols, Q_local.data(), Q_local.stride(1),
                         factorOutput);
          if(debug_) {
            scalarComm_->barrier();
            if(myRank == 0) {
              err_ << "-- Finished DistTsqr::explicit_Q" << endl;
            }
          }
          // Verify the factorization
          auto result =
            global_verify(numCols, numCols, A_local.data(),
                          A_local.stride(1), Q_local.data(),
                          Q_local.stride(1), R.data(), R.stride(1),
                          scalarComm_.get());
          if(debug_) {
            scalarComm_->barrier();
            if(myRank == 0) {
              err_ << "-- Finished global_verify" << endl;
            }
          }
          reportResults("DistTsqr", numCols, result,
                        additionalFieldNames, additionalData,
                        printFieldNames && (! printedFieldNames));
          if(printFieldNames && (! printedFieldNames)) {
            printedFieldNames = true;
          }
        }

        // Test DistTsqr::factorExplicit()
        if(testFactorExplicit_) {
          // Factor the matrix and compute the explicit Q factor, both
          // in a single operation.
          par.factorExplicit(R.view(), Q_local.view());
          if(debug_) {
            scalarComm_->barrier();
            if(myRank == 0) {
              err_ << "-- Finished DistTsqr::factorExplicit" << endl;
            }
          }

          if(printMatrices_) {
            if(myRank == 0) {
              err_ << std::endl << "Computed Q factor:" << std::endl;
            }
            printGlobalMatrix(err_, Q_local, scalarComm_.get(),
                              ordinalComm_.get());
            if(myRank == 0) {
              err_ << std::endl << "Computed R factor:" << std::endl;
              print_local_matrix (err_, R.extent(0), R.extent(1),
                                  R.data(), R.stride(1));
              err_ << std::endl;
            }
          }

          // Verify the factorization
          result_type result =
            global_verify(numCols, numCols, A_local.data(),
                          A_local.stride(1), Q_local.data(),
                          Q_local.stride(1), R.data(), R.stride(1),
                          scalarComm_.get());
          if(debug_) {
            scalarComm_->barrier();
            if(myRank == 0) {
              err_ << "-- Finished global_verify" << endl;
            }
          }
          reportResults("DistTsqrRB", numCols, result,
                        additionalFieldNames, additionalData,
                        printFieldNames && (! printedFieldNames));
          if(printFieldNames && (! printedFieldNames)) {
            printedFieldNames = true;
          }
        }
      }

    private:
      /// Report verification results.  Call on ALL MPI processes, not
      /// just Process 0.
      ///
      /// \param method [in] String to print before reporting results
      /// \param numCols [in] Number of columns in the matrix tested.
      /// \param result [in] (relative residual, orthogonality)
      void
      reportResults (const std::string& method,
                     const Ordinal numCols,
                     const result_type& result,
                     const std::string& additionalFieldNames,
                     const std::string& additionalData,
                     const bool printFieldNames)
      {
        using std::endl;

        const int numProcs = scalarComm_->size();
        const int myRank = scalarComm_->rank();

        if(myRank == 0) {
          if(humanReadable_) {
            out_ << method << " accuracy results:" << endl
                 << "Scalar: " << scalarTypeName_ << endl
                 << "numCols: " << numCols << endl
                 << "Number of (MPI) processes: " << numProcs << endl
                 << "Absolute residual $\\| A - Q R \\|_2: "
                 << result[0] << endl
                 << "Absolute orthogonality $\\| I - Q^* Q \\|_2$: "
                 << result[1] << endl
                 << "Test matrix norm $\\| A \\|_F$: "
                 << result[2] << endl;
          }
          else {
            // Use scientific notation for floating-point numbers
            out_ << std::scientific;

            if(printFieldNames) {
              out_ << "%method,scalarType,numCols,numProcs"
                ",absFrobResid,absFrobOrthog,frobA";
              if(! additionalFieldNames.empty())
                out_ << "," << additionalFieldNames;
              out_ << endl;
            }

            out_ << method
                 << "," << scalarTypeName_
                 << "," << numCols
                 << "," << numProcs
                 << "," << result[0]
                 << "," << result[1]
                 << "," << result[2];
            if(! additionalData.empty()) {
              out_ << "," << additionalData;
            }
            out_ << endl;
          }
        }
      }

      void
      testProblem(Matrix<Ordinal, Scalar>& A_local,
                  Matrix<Ordinal, Scalar>& Q_local,
                  Matrix<Ordinal, Scalar>& R,
                  const Ordinal numCols)
      {
        const Ordinal numRowsLocal = numCols;

        // A_local: Space for the matrix A to factor -- local to each
        //   processor.
        //
        // A_global: Global matrix (only nonempty on Proc 0); only
        //   used temporarily.
        Matrix<Ordinal, Scalar> A_global;

        // This modifies A_local on all procs, and A_global on Proc 0.
        par_tsqr_test_problem(gen_, A_local, A_global, numCols, scalarComm_);

        if(printMatrices_) {
          const int myRank = scalarComm_->rank();
          if(myRank == 0) {
            err_ << "Input matrix A:" << std::endl;
          }
          printGlobalMatrix(err_, A_local, scalarComm_.get(),
                            ordinalComm_.get());
          if(myRank == 0) {
            err_ << std::endl;
          }
        }

        // Copy the test problem input into R, since the factorization
        // will overwrite it in place with the final R factor.
        R.reshape(numCols, numCols);
        deep_copy(R, Scalar{});
        deep_copy(R, A_local);

        // Prepare space in which to construct the explicit Q factor
        // (local component on this processor)
        Q_local.reshape(numRowsLocal, numCols);
        deep_copy(Q_local, Scalar {});
      }
    };

    /// \class DistTsqrBenchmarker
    /// \brief Generic version of DistTsqr performance test.
    template< class Ordinal, class Scalar>
    class DistTsqrBenchmarker {
      TSQR::Random::NormalGenerator<Ordinal, Scalar> gen_;
      Teuchos::RCP<MessengerBase<Scalar>> scalarComm_;
      Teuchos::RCP<MessengerBase<double>> doubleComm_;
      std::string scalarTypeName_;

      std::ostream& out_;
      std::ostream& err_;
      const bool testFactorExplicit_;
      const bool testFactorImplicit_;
      const bool humanReadable_;
      const bool debug_;

    public:
      using ordinal_type = Ordinal;
      using scalar_type = Scalar;
      using timer_type = Teuchos::Time;

      /// \brief Constructor, with custom seed value
      ///
      /// \param scalarComm [in/out] Communicator object over which
      ///   to test.
      /// \param doubleComm [in/out] Communicator object for doubles,
      ///   used for finding the min and max of timing results over
      ///   all the MPI processes.
      /// \param seed [in] 4-element vector; the random seed input of
      ///   TSQR::Random::NormalGenerator (which see, since there are
      ///   restrictions on the set of valid seeds)
      /// \param scalarTypeName [in] Human-readable name of the Scalar
      ///   template type parameter
      /// \param out [out] Output stream to which to write results
      /// \param err [out] Output stream to which to write any
      ///   debugging outputs (if applicable) or errors
      /// \param testFactorExplicit [in] Whether to test
      ///   DistTsqr::factorExplicit()
      /// \param testFactorImplicit [in] Whether to test
      ///   DistTsqr::factor() and DistTsqr::explicit_Q()
      /// \param humanReadable [in] Whether printed results should be
      ///   easy for humans to read (vs. easy for parsers to parse)
      /// \param debug [in] Whether to write verbose debug output to
      ///   err
      DistTsqrBenchmarker(const Teuchos::RCP<MessengerBase<Scalar>>& scalarComm,
                          const Teuchos::RCP<MessengerBase<double>>& doubleComm,
                          const std::vector<int>& seed,
                          const std::string& scalarTypeName,
                          std::ostream& out,
                          std::ostream& err,
                          const bool testFactorExplicit,
                          const bool testFactorImplicit,
                          const bool humanReadable,
                          const bool debug) :
        gen_(seed),
        scalarComm_(scalarComm),
        doubleComm_(doubleComm),
        scalarTypeName_(scalarTypeName),
        out_(out),
        err_(err),
        testFactorExplicit_(testFactorExplicit),
        testFactorImplicit_(testFactorImplicit),
        humanReadable_(humanReadable),
        debug_(debug)
      {}

      /// \brief Constructor, with default seed value
      ///
      /// This constructor sets a default seed (for the pseudorandom
      /// number generator), which is the same seed (0,0,0,1) each
      /// time.
      ///
      /// \param scalarComm [in/out] Communicator object over which
      ///   to test.
      /// \param doubleComm [in/out] Communicator object for doubles,
      ///   used for finding the min and max of timing results over
      ///   all the MPI processes.
      /// \param scalarTypeName [in] Human-readable name of the Scalar
      ///   template type parameter
      /// \param out [out] Output stream to which to write results
      /// \param err [out] Output stream to which to write any
      ///   debugging outputs (if applicable) or errors
      /// \param testFactorExplicit [in] Whether to test
      ///   DistTsqr::factorExplicit()
      /// \param testFactorImplicit [in] Whether to test
      ///   DistTsqr::factor() and DistTsqr::explicit_Q()
      /// \param humanReadable [in] Whether printed results should be
      ///   easy for humans to read (vs. easy for parsers to parse)
      /// \param debug [in] Whether to write verbose debug output to
      ///   err
      DistTsqrBenchmarker(const Teuchos::RCP<MessengerBase<Scalar>>& scalarComm,
                          const Teuchos::RCP<MessengerBase<double>>& doubleComm,
                          const std::string& scalarTypeName,
                          std::ostream& out,
                          std::ostream& err,
                          const bool testFactorExplicit,
                          const bool testFactorImplicit,
                          const bool humanReadable,
                          const bool debug) :
        scalarComm_(scalarComm),
        doubleComm_(doubleComm),
        scalarTypeName_(scalarTypeName),
        out_(out),
        err_(err),
        testFactorExplicit_(testFactorExplicit),
        testFactorImplicit_(testFactorImplicit),
        humanReadable_(humanReadable),
        debug_(debug)
      {}

      /// \brief Get seed vector for pseudorandom number generator
      ///
      /// Fill seed (changing size of vector as necessary) with the
      /// seed vector used by the pseudorandom number generator.  You
      /// can use this to resume the pseudorandom number stream from
      /// where you last were.
      void
      getSeed(std::vector<int>& seed) const
      {
        gen_.getSeed(seed);
      }

      /// \brief Run the DistTsqr benchmark
      ///
      /// \param numTrials [in] Number of times to repeat the computation
      ///   in a single timing run
      /// \param numCols [in] Number of columns in the matrix to test.
      ///   Number of rows := (# MPI processors) * ncols
      void
      benchmark(const int numTrials,
                const Ordinal numCols,
                const std::string& additionalFieldNames,
                const std::string& additionalData,
                const bool printFieldNames)
      {
        using std::endl;

        // Set up test problem.
        Matrix<Ordinal, Scalar> A_local, Q_local, R;
        testProblem(A_local, Q_local, R, numCols);

        // Set up TSQR implementation.
        DistTsqr<Ordinal, Scalar> par;
        par.init(scalarComm_);

        // Whether we've printed field names (i.e., column headers)
        // yet.  Only matters for non-humanReadable output.
        bool printedFieldNames = false;

        if(testFactorImplicit_) {
          std::string timerName("DistTsqr");

          // Throw away some number of runs, because some MPI libraries
          // (recent versions of OpenMPI at least) do autotuning for the
          // first few collectives calls.
          const int numThrowAwayRuns = 5;
          for(int runNum = 0; runNum < numThrowAwayRuns; ++runNum) {
            auto factorOutput = par.factor(R.view());
            par.explicit_Q(numCols, Q_local.data(),
                           Q_local.stride(1), factorOutput);
          }

          // Now do the actual timing runs.  Benchmark DistTsqr
          // (factor() and explicit_Q()) for numTrials trials.
          timer_type timer (timerName);
          timer.start();
          for(int trialNum = 0; trialNum < numTrials; ++trialNum) {
            auto factorOutput = par.factor(R.view());
            par.explicit_Q(numCols, Q_local.data(),
                           Q_local.stride(1), factorOutput);
          }
          // Cumulative timing on this MPI process.  "Cumulative"
          // means the elapsed time of numTrials executions.
          const double localCumulativeTiming = timer.stop();

          // reportResults() must be called on all processes, since this
          // figures out the min and max timings over all processes.
          reportResults(timerName, numTrials, numCols,
                        localCumulativeTiming, additionalFieldNames,
                        additionalData,
                        printFieldNames && (! printedFieldNames));
          if(printFieldNames && (! printedFieldNames)) {
            printedFieldNames = true;
          }
        }

        if(testFactorExplicit_) {
          std::string timerName ("DistTsqrRB");

          // Throw away some number of runs, because some MPI libraries
          // (recent versions of OpenMPI at least) do autotuning for the
          // first few collectives calls.
          const int numThrowAwayRuns = 5;
          for(int runNum = 0; runNum < numThrowAwayRuns; ++runNum) {
            par.factorExplicit(R.view(), Q_local.view());
          }

          // Benchmark DistTsqr::factorExplicit() for numTrials trials.
          timer_type timer(timerName);
          timer.start();
          for(int trialNum = 0; trialNum < numTrials; ++trialNum) {
            par.factorExplicit(R.view(), Q_local.view());
          }
          // Cumulative timing on this MPI process.
          // "Cumulative" means the elapsed time of numTrials executions.
          const double localCumulativeTiming = timer.stop();

          // Report cumulative (not per-invocation) timing results
          reportResults(timerName, numTrials, numCols, localCumulativeTiming,
                        additionalFieldNames, additionalData,
                        printFieldNames && (! printedFieldNames));
          if(printFieldNames && (! printedFieldNames)) {
            printedFieldNames = true;
          }

          // Per-invocation timings (for factorExplicit() benchmark
          // only).  localTimings were computed on this MPI process;
          // globalTimings are statistical summaries of those over
          // all MPI processes.  We only collect that data for
          // factorExplicit().
          std::vector<TimeStats> localTimings;
          std::vector<TimeStats> globalTimings;
          par.getFactorExplicitTimings(localTimings);
          for(size_t k = 0; k < localTimings.size(); ++k) {
            globalTimings.push_back
              (globalTimeStats(*doubleComm_, localTimings[k]));
          }
          std::vector<std::string> timingLabels;
          par.getFactorExplicitTimingLabels(timingLabels);

          if(humanReadable_) {
            out_ << timerName << " per-invocation benchmark results:" << endl;
          }
          const std::string labelLabel("label,scalarType");
          for (size_t k = 0; k < timingLabels.size(); ++k) {
            // Only print column headers (i.e., field names) once, if at all.
            const bool printHeaders = (k == 0) && printFieldNames;
            globalTimings[k].print (out_, humanReadable_,
                                    timingLabels[k] + "," + scalarTypeName_,
                                    labelLabel, printHeaders);
          }
        }
      }

    private:
      /// Report timing results to the given output stream
      ///
      /// \param method [in] String to print before reporting results
      /// \param numTrials [in] Number of times to repeat the computation
      ///   in a single timing run
      /// \param numCols [in] Number of columns in the matrix to test.
      ///   Number of rows := (# MPI processors) * ncols
      /// \param timing [in] Total benchmark time, as measured on this
      ///   MPI process.  This may differ on each process; we report
      ///   the min and the max.
      ///
      /// \warning Call on ALL MPI processes, not just Rank 0!
      void
      reportResults(const std::string& method,
                    const int numTrials,
                    const ordinal_type numCols,
                    const double localTiming,
                    const std::string& additionalFieldNames,
                    const std::string& additionalData,
                    const bool printFieldNames)
      {
        using std::endl;

        // Find min and max timing over all MPI processes
        TimeStats localStats;
        localStats.update (localTiming);
        TimeStats globalStats = globalTimeStats (*doubleComm_, localStats);

        // Only Rank 0 prints the final results.
        const bool printResults = (doubleComm_->rank() == 0);
        if(printResults) {
          const int numProcs = doubleComm_->size();
          if(humanReadable_) {
            out_ << method << " cumulative benchmark results "
                 << "(total time over all trials):" << endl
                 << "Scalar: " << scalarTypeName_ << endl
                 << "numCols: " << numCols << endl
                 << "MPI comm size: " << numProcs << endl
                 << "numTrials: " << numTrials << endl
                 << "Min timing (s): " << globalStats.min() << endl
                 << "Mean timing (s): " << globalStats.mean() << endl
                 << "Max timing (s): " << globalStats.max() << endl
                 << endl;
          }
          else {
            // Use scientific notation for floating-point numbers
            out_ << std::scientific;

            if(printFieldNames) {
              out_ << "%method,scalarType,numCols,numProcs,numTrials"
                   << ",minTiming,meanTiming,maxTiming";
              if(! additionalFieldNames.empty()) {
                out_ << "," << additionalFieldNames;
              }
              out_ << endl;
            }

            out_ << method
                 << "," << scalarTypeName_
                 << "," << numCols
                 << "," << numProcs
                 << "," << numTrials
                 << "," << globalStats.min()
                 << "," << globalStats.mean()
                 << "," << globalStats.max();
            if(! additionalData.empty()) {
              out_ << "," << additionalData;
            }
            out_ << endl;
          }
        }
      }

      void
      testProblem(Matrix<Ordinal, Scalar>& A_local,
                  Matrix<Ordinal, Scalar>& Q_local,
                  Matrix<Ordinal, Scalar>& R,
                  const Ordinal numCols)
      {
        const Ordinal numRowsLocal = numCols;

        // A_local: Space for the matrix A to factor -- local to each
        //   (MPI) process.
        //
        // A_global: Global matrix (only nonempty on Proc 0); only
        //   used temporarily.
        Matrix<Ordinal, Scalar> A_global;

        // This modifies A_local on all procs, and A_global on Proc 0.
        par_tsqr_test_problem(gen_, A_local, A_global, numCols,
                              scalarComm_);

        // Copy the test problem input into R, since the factorization
        // will overwrite it in place with the final R factor.
        R.reshape(numCols, numCols);
        deep_copy(R, A_local);

        // Prepare space in which to construct the explicit Q factor
        // (local component on this processor)
        Q_local.reshape(numRowsLocal, numCols);
        deep_copy(Q_local, Scalar {});
      }
    };
  } // namespace Test
} // namespace TSQR

template<class Ordinal, class Scalar>
class MessengerPairMaker {
public:
  using ordinal_type = Ordinal;
  using scalar_type = Scalar;

  using pair_type = std::pair<
    Teuchos::RCP<TSQR::MessengerBase<ordinal_type>>,
    Teuchos::RCP<TSQR::MessengerBase<scalar_type>>
    >;

  static pair_type
  makePair(const Teuchos::RCP<const Teuchos::Comm<int>>& comm)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_implicit_cast;
    using TSQR::MessengerBase;
    using TSQR::TeuchosMessenger;

    auto derivedOrdinalComm =
      rcp(new TeuchosMessenger<ordinal_type>(comm));
    auto ordinalComm =
      rcp_implicit_cast<MessengerBase<ordinal_type>>(derivedOrdinalComm);
    auto derivedScalarComm =
      rcp (new TeuchosMessenger<scalar_type>(comm));
    auto scalarComm =
      rcp_implicit_cast<MessengerBase<scalar_type>>(derivedScalarComm);

    return {ordinalComm, scalarComm};
  }
};

#define TSQR_TEST_DIST_TSQR( ScalarType, typeString )                   \
  do {                                                                  \
    using TSQR::Test::DistTsqrVerifier;                                 \
    using LO = int;                                                     \
    using SC = ScalarType;                                              \
    using verifier_type = DistTsqrVerifier<LO, SC>;                     \
                                                                        \
    std::string scalarTypeName (typeString);                            \
    auto messPair = MessengerPairMaker<LO, SC>::makePair (comm);        \
    verifier_type verifier (messPair.first, messPair.second, seed,      \
                            scalarTypeName, out, err,                   \
                            testFactorExplicit, testFactorImplicit,     \
                            humanReadable, printMatrices, debug);       \
    verifier.verify (numCols, params.additionalFieldNames,              \
                     params.additionalData, params.printFieldNames);    \
    verifier.getSeed (seed);                                            \
  } while (false)


#define TSQR_BENCHMARK_DIST_TSQR( theType, typeString )                 \
  do {                                                                  \
    using TSQR::Test::DistTsqrBenchmarker;                              \
    using Teuchos::RCP;                                                 \
    using SC = theType;                                                 \
    using base_messenger_type = TSQR::MessengerBase<SC>;                \
    using derived_messenger_type = TSQR::TeuchosMessenger<SC>;          \
    using derived_messenger_ptr = RCP<derived_messenger_type>;          \
    using benchmarker_type = DistTsqrBenchmarker<int, SC>;              \
                                                                        \
    std::string scalarTypeName (typeString);                            \
    derived_messenger_ptr scalarCommDerived                             \
      (new derived_messenger_type (comm));                              \
    auto scalarComm =                                                   \
      rcp_implicit_cast<base_messenger_type> (scalarCommDerived);       \
    benchmarker_type benchmarker (scalarComm, doubleComm, seed,         \
                                  scalarTypeName, out, err,             \
                                  testFactorExplicit,                   \
                                  testFactorImplicit,                   \
                                  humanReadable, debug);                \
    benchmarker.benchmark (numTrials, numCols,                          \
                           params.additionalFieldNames,                 \
                           params.additionalData,                       \
                           params.printFieldNames);                     \
    benchmarker.getSeed (seed);                                         \
  } while (false)

/// \class DistTsqrTestParameters
/// \brief Encapsulates values of command-line parameters
struct DistTsqrTestParameters {
  int numCols = 10;
  int numTrials = 10;
  bool verify = true;
  bool benchmark = false;
  bool testReal = true;
#ifdef HAVE_TPETRATSQR_COMPLEX
  bool testComplex = true;
#else
  bool testComplex = false;
#endif // HAVE_TPETRATSQR_COMPLEX
  bool testFactorExplicit = true;
  bool testFactorImplicit = true;
  bool printFieldNames = true;
  bool printTrilinosTestStuff = true;
  bool humanReadable = false;
  bool printMatrices = false;
  bool debug = false;

  std::string additionalFieldNames;
  std::string additionalData;
};

static void
verify(Teuchos::RCP<const Teuchos::Comm<int>> comm,
       const DistTsqrTestParameters& params,
       std::ostream& out,
       std::ostream& err,
       std::vector<int>& seed,
       const bool useSeed)
{
  const bool testReal = params.testReal;
  const bool testComplex = params.testComplex;
  const int numCols = params.numCols;
  const bool testFactorExplicit = params.testFactorExplicit;
  const bool testFactorImplicit = params.testFactorImplicit;
  const bool humanReadable = params.humanReadable;
  const bool printMatrices = params.printMatrices;
  const bool debug = params.debug;

  if(! useSeed) {
    seed.resize(4);
    seed[0] = 0;
    seed[1] = 0;
    seed[2] = 0;
    seed[3] = 1;
  }
  if(testReal) {
    TSQR_TEST_DIST_TSQR( float, "float" );
    TSQR_TEST_DIST_TSQR( double, "double" );
  }
  if(testComplex) {
#ifdef HAVE_TPETRATSQR_COMPLEX
    using std::complex;

    TSQR_TEST_DIST_TSQR( complex<float>, "complex<float>" );
    TSQR_TEST_DIST_TSQR( complex<double>, "complex<double>" );

#else // Don't HAVE_TPETRATSQR_COMPLEX
    throw std::logic_error("TSQR was not built with complex "
                           "arithmetic support");
#endif // HAVE_TPETRATSQR_COMPLEX
  }
}


static void
benchmark(Teuchos::RCP<const Teuchos::Comm<int>> comm,
          const DistTsqrTestParameters& params,
          std::ostream& out,
          std::ostream& err,
          std::vector<int>& seed,
          const bool useSeed)
{
  const bool testReal = params.testReal;
  const bool testComplex = params.testComplex;
  const int numCols = params.numCols;
  const int numTrials = params.numTrials;
  const bool testFactorExplicit = params.testFactorExplicit;
  const bool testFactorImplicit = params.testFactorImplicit;
  const bool humanReadable = params.humanReadable;
  const bool debug = params.debug;

  if(! useSeed) {
    seed.resize(4);
    seed[0] = 0;
    seed[1] = 0;
    seed[2] = 0;
    seed[3] = 1;
  }
  using Teuchos::rcp;
  auto doubleCommSub =
    rcp(new TSQR::TeuchosMessenger<double>(comm));
  using TSQR::MessengerBase;
  using Teuchos::rcp_implicit_cast;
  auto doubleComm =
    rcp_implicit_cast<MessengerBase<double>>(doubleCommSub);

  if(testReal) {
    TSQR_BENCHMARK_DIST_TSQR( float, "float" );
    TSQR_BENCHMARK_DIST_TSQR( double, "double" );
  }
  if(testComplex) {
#ifdef HAVE_TPETRATSQR_COMPLEX
    using std::complex;

    TSQR_BENCHMARK_DIST_TSQR( complex<float>, "complex<float>" );
    TSQR_BENCHMARK_DIST_TSQR( complex<double>, "complex<double>" );

#else // Don't HAVE_TPETRATSQR_COMPLEX
    throw std::logic_error("TSQR was not built with complex "
                           "arithmetic support");
#endif // HAVE_TPETRATSQR_COMPLEX
  }
}


/// \brief Parse command-line options for this test
///
/// \param argc [in] As usual in C(++)
/// \param argv [in] As usual in C(++)
/// \param allowedToPrint [in] Whether this (MPI) process is allowed
///   to print to stdout/stderr.  Different per (MPI) process.
/// \param printedHelp [out] Whether this (MPI) process printed the
///   "help" display (summary of command-line options)
///
/// \return Encapsulation of command-line options
static DistTsqrTestParameters
parseOptions(int argc,
             char* argv[],
             std::ostream& err,
             bool& printedHelp)
{
  using std::endl;
  printedHelp = false;

  // Command-line parameters, set to their default values.
  DistTsqrTestParameters params {};
  try {
    constexpr bool throwExceptions = true;
    constexpr bool recognizeAllOptions = true;
    using CLP = Teuchos::CommandLineProcessor;
    CLP cmdLineProc(throwExceptions, recognizeAllOptions);

    const char docString[] = "This program tests TSQR::DistTsqr, which "
      "implements the internode-parallel part of TSQR (TSQR::Tsqr).  "
      "Accuracy and performance tests are included.";
    cmdLineProc.setDocString(docString);
    cmdLineProc.setOption("verify",
        "noverify",
        &params.verify,
        "Test accuracy");
    cmdLineProc.setOption("benchmark",
        "nobenchmark",
        &params.benchmark,
        "Test performance");
    cmdLineProc.setOption("implicit",
        "noimplicit",
        &params.testFactorImplicit,
        "Test DistTsqr\'s factor() and explicit_Q()");
    cmdLineProc.setOption("explicit",
        "noexplicit",
        &params.testFactorExplicit,
        "Test DistTsqr\'s factorExplicit()");
    cmdLineProc.setOption("field-names",
        &params.additionalFieldNames,
        "Any additional field name(s) (comma-delimited "
        "string) to add to the benchmark output.  Empty "
        "by default.  Good for things known when invoking "
        "the benchmark executable, but not (easily) known "
        "inside the benchmark -- e.g., environment "
        "variables.");
    cmdLineProc.setOption("output-data",
        &params.additionalData,
        "Any additional data to add to the output, "
        "corresponding to the above field name(s). "
        "Empty by default.");
    cmdLineProc.setOption("print-field-names",
        "no-print-field-names",
        &params.printFieldNames,
        "Print field names (for machine-readable output only)");
    cmdLineProc.setOption("print-trilinos-test-stuff",
        "no-print-trilinos-test-stuff",
        &params.printTrilinosTestStuff,
        "Print output that makes the Trilinos test "
        "framework happy (but makes benchmark results "
        "parsing scripts unhappy)");
    cmdLineProc.setOption("print-matrices",
        "no-print-matrices",
        &params.printMatrices,
        "Print global test matrices and computed results to stderr");
    cmdLineProc.setOption("debug",
        "nodebug",
        &params.debug,
        "Print debugging information");
    cmdLineProc.setOption("human-readable",
        "machine-readable",
        &params.humanReadable,
        "If set, make output easy to read by humans "
        "(but hard to parse)");
    cmdLineProc.setOption("ncols",
        &params.numCols,
        "Number of columns in the test matrix");
    cmdLineProc.setOption("ntrials",
        &params.numTrials,
        "Number of trials (only used when \"--benchmark\"");
    cmdLineProc.setOption("real",
        "noreal",
        &params.testReal,
        "Test real arithmetic routines");
    cmdLineProc.setOption("complex",
        "nocomplex",
        &params.testComplex,
        "Test complex arithmetic routines (only set to true if "
        "complex arithmetic support was enabled at configure "
        "time)");
    cmdLineProc.parse (argc, argv);
  }
  catch (Teuchos::CommandLineProcessor::UnrecognizedOption& e) {
    err << "Unrecognized command-line option: " << e.what() << endl;
    throw e;
  }
  catch (Teuchos::CommandLineProcessor::HelpPrinted& e) {
    printedHelp = true;
  }

  // Validate command-line options.  We provide default values
  // for unset options, so we don't have to validate those.
  TEUCHOS_TEST_FOR_EXCEPTION
    (params.numCols <= 0, std::invalid_argument,
     "You set --numCols=" << params.numCols << ".  The number of "
     "columns in the matrix to test must be positive.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (params.benchmark && params.numTrials < 1, std::invalid_argument,
     "\"--benchmark\" option requires positive --numTrials, but you "
     "set --numTrials=" << params.numTrials << ".");
#ifndef HAVE_TPETRATSQR_COMPLEX
    TEUCHOS_TEST_FOR_EXCEPTION
      (params.testComplex, std::invalid_argument, "Complex "
       "arithmetic support was not enabled at configure time, "
       "but you set --testComplex.");
#endif // HAVE_TPETRATSQR_COMPLEX
  return params;
}

int
main(int argc, char *argv[])
{
  TSQR::Test::MpiAndKokkosScope testScope(&argc, &argv);
  auto comm = testScope.getComm();
  std::ostream& out = testScope.outStream();
  std::ostream& err = testScope.errStream();

  // Fetch command-line parameters.
  bool printedHelp = false;
  auto params = parseOptions(argc, argv, err, printedHelp);
  if(printedHelp) {
    return EXIT_SUCCESS;
  }
  bool success = false;
  constexpr bool actually_print_caught_exceptions = true;
  try {
    if(params.verify) {
      std::vector<int> seed(4);
      const bool useSeed = false;
      verify(comm, params, out, err, seed, useSeed);
    }

    if(params.benchmark) {
      std::vector<int> seed(4);
      const bool useSeed = false;
      benchmark(comm, params, out, err, seed, useSeed);
    }

    success = true;

    if(params.printTrilinosTestStuff) {
      // The Trilinos test framework expects a message like this.
      out << "\nEnd Result: TEST PASSED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS
    (actually_print_caught_exceptions, err, success);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
