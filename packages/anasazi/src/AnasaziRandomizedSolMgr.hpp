// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_RANDOMIZED_SOLMGR_HPP
#define ANASAZI_RANDOMIZED_SOLMGR_HPP

/*! \file AnasaziRandomizedSolMgr.hpp
  \brief The Anasazi::RandomizedSolMgr approximates largest eigenvalues/eigenvectors
  by performing a simple Rayleigh-Ritz projection over a random block of vectors. 
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"

#include "AnasaziBasicSort.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziICGSOrthoManager.hpp"

#include "AnasaziStatusTestMaxIters.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziStatusTestOutput.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziOutputStreamTraits.hpp"
#include "AnasaziSolverUtils.hpp"

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_LAPACK.hpp"

/*! \class Anasazi::RandomizedSolMgr
  \brief The Anasazi::RandomizedSolMgr approximates largest eigenvalues/eigenvectors
  by performing a simple Rayleigh-Ritz projection over a random block of vectors. 

  This algorithm is well-known. The implementation is adapted from
  that found in N. Halko et al, "Finding Structure with Randomness:
  Probabilistic Algorithms for Constructing Approximate Matrix 
  Decompositions", SIAM Review, 53 (2011), pp. 217-288, 
  https://doi.org/10.1137/090771806

  \ingroup anasazi_solver_framework

  \author Jennifer A. Loe, Erik G. Boman, Heather M. Switzer
  */

namespace Anasazi {

  namespace Experimental {

    template<class ScalarType, class MV, class OP>
      class RandomizedSolMgr : public SolverManager<ScalarType,MV,OP> {

        private:
          typedef int OT; 
          typedef MultiVecTraits<ScalarType,MV> MVT;
          typedef OperatorTraits<ScalarType,MV,OP>  OPT;
          typedef Teuchos::ScalarTraits<ScalarType> SCT;
          typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MT;
          const ScalarType ONE  = SCT::one();

        public:

          //!@name Constructors/Destructor
          //@{

          /*! \brief Basic constructor for RandomizedSolMgr.
           *
           * This constructor accepts the Eigenproblem to be solved in addition
           * to a parameter list of options for the solver manager. These options include the following:
           *   - "Which" - a \c string specifying the desired eigenvalues: LM or LR. This method does not support finding the smallest eigenvalues. Default: LM
           *   - "Block Size" - an \c int specifying the block size to be used by the underlying LOBPCG solver. Default: problem->getNEV()
           *   - "Maximum Iterations" - an \c int specifying the maximum number of iterations the underlying solver is allowed to perform. For the randomized solver, this corresponds to the number of times we multiply A by the initial vectors. Default: 5
           *   - "Orthogonalization" - a \c string specifying the desired orthogonalization: DGKS, ICGS, and SVQB. Default: "SVQB"
           *   - "Orthogonalization Frequency" - a \c int specifying how many iterations should pass before the system. Default: 0
           *   - "Residual Frequency" - a \c int specifying how many iterations should pass before checking the residuals of the system. Default: 0
           *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Anasazi::Errors
           *   - "Output Stream" - a reference-counted pointer to the formatted output stream where all
           *                      solver output is sent.  Default: Teuchos::getFancyOStream ( Teuchos::rcpFromRef (std::cout) )
           *   - "Output Processor" - an \c int specifying the MPI processor that will print solver/timer details.  Default: 0
           *   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: machine precision
           */
          RandomizedSolMgr( const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
              Teuchos::ParameterList &pl );

          //! Destructor.
          virtual ~RandomizedSolMgr() {};
          //@}

          //! @name Accessor methods
          //@{

          const Eigenproblem<ScalarType,MV,OP>& getProblem() const {
            return *problem_;
          }

          int getNumIters() const {
            return numIters_;
          }

          int getNumFailed() const {
            return numFailed_;
          }

          //@}

          //! @name Solver application methods
          //@{

          /*! \brief This method performs possibly repeated calls to the underlying eigensolver's iterate() routine
           * until the problem has been solved (as decided by the solver manager) or the solver manager decides to
           * quit.
           *
           * \returns ::ReturnType specifying:
           *    - ::Converged: the eigenproblem was solved to the specification required by the solver manager.
           *    - ::Unconverged: the eigenproblem was not solved to the specification desired by the solver manager
           */
          ReturnType solve();
          //@}

        private:
          Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > problem_;
          Teuchos::RCP<Teuchos::FancyOStream> osp_;
          std::string whch_;
          MT tol_;
          int osProc_;
          int verb_;
          Teuchos::RCP<Teuchos::Time> timerOrtho_;
          Teuchos::RCP<Teuchos::Time> timerSolve_;
          Teuchos::RCP<Teuchos::Time> timerOp_;
          std::string ortho_;
          int orthoFreq_;
          int resFreq_;
          int blockSize_;
          int maxIters_;
          int numIters_;
          bool trackResNorms_;
          int numFailed_;
      };


    ////////////////////////////////////////////////////////////////////////////////////////
    template<class ScalarType, class MV, class OP>
      RandomizedSolMgr<ScalarType,MV,OP>::RandomizedSolMgr(
          const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
          Teuchos::ParameterList &pl ) :
        problem_(problem),
        whch_("LM"),
        tol_(1e-6),
        osProc_(0),
        verb_(Anasazi::Errors),
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        timerOrtho_(Teuchos::TimeMonitor::getNewTimer("Anasazi: Randomized::Orthogonalization")),
        timerSolve_(Teuchos::TimeMonitor::getNewTimer("Anasazi: Randomized::solve()")),
        timerOp_(Teuchos::TimeMonitor::getNewTimer("Anasazi: Randomized::Operation Op*x")),
#endif
        ortho_("SVQB"),
        orthoFreq_(0),
        resFreq_(0),
        blockSize_(0),
        maxIters_(5),
        numIters_(0),
        trackResNorms_(true)
        {
          TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null,              std::invalid_argument, "Problem not given to solver manager.");
          TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isProblemSet(),              std::invalid_argument, "Problem not set.");
          TEUCHOS_TEST_FOR_EXCEPTION(problem_->getInitVec() == Teuchos::null,std::invalid_argument, "Problem does not contain initial vectors to clone from.");

          whch_ = pl.get("Which",whch_);
          TEUCHOS_TEST_FOR_EXCEPTION(whch_ != "LM" && whch_ != "LR",
              AnasaziError,
              "RandomizedSolMgr: \"Which\" parameter must be LM or LR. Other options not available."); 

          tol_ = pl.get("Convergence Tolerance",tol_);
          TEUCHOS_TEST_FOR_EXCEPTION(tol_ <= 0,
              AnasaziError,
              "RandomizedSolMgr: \"Tolerance\" parameter must be strictly positive.");

          // Create a formatted output stream to print to.
          // See if user requests output processor.
          osProc_ = pl.get("Output Processor", osProc_);

          // If not passed in by user, it will be chosen based upon operator type.
          if (pl.isParameter("Output Stream")) {
            osp_ = Teuchos::getParameter<Teuchos::RCP<Teuchos::FancyOStream> >(pl,"Output Stream");
          }
          else {
            osp_ = OutputStreamTraits<OP>::getOutputStream (*problem_->getOperator(), osProc_);
          }

          // verbosity level
          if (pl.isParameter("Verbosity")) {
            if (Teuchos::isParameterType<int>(pl,"Verbosity")) {
              verb_ = pl.get("Verbosity", verb_);
            } else {
              verb_ = (int)Teuchos::getParameter<Anasazi::MsgType>(pl,"Verbosity");
            }
          }

          // Orthogonalization type
          ortho_ = pl.get("Orthogonalization","SVQB");

          blockSize_= pl.get("Block Size",problem_->getNEV());
          TEUCHOS_TEST_FOR_EXCEPTION(blockSize_ <= 0,
              AnasaziError,
              "RandomizedSolMgr: \"Block Size\" parameter must be strictly positive.");

          maxIters_ = pl.get("Maximum Iterations",maxIters_);
          trackResNorms_ = pl.get("Track Residuals",true);

          // How many iterations between orthogonalizations
          if (pl.isParameter("Orthogonalization Frequency")) {
            if (Teuchos::isParameterType<int>(pl,"Orthogonalization Frequency")) {
              orthoFreq_ = pl.get("Orthogonalization Frequency", orthoFreq_);
            } else {
              orthoFreq_ = (int)Teuchos::getParameter<Anasazi::MsgType>(pl,"Orthogonalization Frequency");
            }
          }

          // How many iterations between checking the residuals
          if (pl.isParameter("Residual Frequency")) {
            if (Teuchos::isParameterType<int>(pl,"Residual Frequency")) {
              resFreq_ = pl.get("Residual Frequency", resFreq_);
            } else {
              resFreq_ = (int)Teuchos::getParameter<Anasazi::MsgType>(pl,"Residual Frequency");
            }
          }
        }

    ////////////////////////////////////////////////////////////////////////////////////////
    template<class ScalarType, class MV, class OP>
      ReturnType
      RandomizedSolMgr<ScalarType,MV,OP>::solve() {

        // Sort manager
        Teuchos::RCP<BasicSort<MT> > sorter = Teuchos::rcp( new BasicSort<MT> );
        sorter->setSortType(whch_);
        std::vector<int> order(blockSize_);     /* Permutation array for sorting the eigenvectors */ 
        SolverUtils<ScalarType,MV,OP> msutils;

        // Output manager
        Teuchos::RCP<OutputManager<ScalarType> > printer = Teuchos::rcp( new OutputManager<ScalarType>(verb_,osp_) );

        // Eigensolution manager
        Eigensolution<ScalarType,MV> sol;
        sol.numVecs = 0;
        problem_->setSolution(sol);	/* In case there is an exception thrown */

        // ortho manager 
        Teuchos::RCP<Anasazi::OrthoManager<ScalarType,MV> > orthoMgr;
        int rank;
        if (ortho_=="SVQB") {
          orthoMgr = Teuchos::rcp( new Anasazi::SVQBOrthoManager<ScalarType,MV,OP>());
        } else if (ortho_=="DGKS") {
          orthoMgr = Teuchos::rcp( new Anasazi::BasicOrthoManager<ScalarType,MV,OP>());
        } else if (ortho_=="ICGS") {
          orthoMgr = Teuchos::rcp( new Anasazi::ICGSOrthoManager<ScalarType,MV,OP>());
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(ortho_!="SVQB"&&ortho_!="DGKS"&&ortho_!="ICGS",std::logic_error,"Anasazi::RandomSolver Invalid orthogonalization type.");
        }

        if(blockSize_ < problem_->getNEV()){ 
          printer->stream(Warnings) << "Warning! Block size smaller than number evals. Increasing Block Size to num evals." << std::endl;
          blockSize_ = problem_->getNEV();
        }

        /* Grab some Multivector to Clone
         * in practice, getInitVec() should always provide this, but it is possible to use a 
         * Eigenproblem with nothing in getInitVec() by manually initializing with initialize(); 
         * in case of that strange scenario, we will try to Clone from V_; first resort to getInitVec(), 
         * because we would like to clear the storage associated with V_ so we have room for the new V_ */
        Teuchos::RCP<MV> randVecs;
        TEUCHOS_TEST_FOR_EXCEPTION(problem_->getInitVec() == Teuchos::null,std::invalid_argument,
            "Anasazi::Randomized: eigenproblem did not specify initial vectors to clone from.");
        if(MVT::GetNumberVecs(*(problem_->getInitVec()))==blockSize_){
          randVecs = MVT::CloneCopy(*(problem_->getInitVec()));
        } else {
          randVecs = MVT::Clone(*(problem_->getInitVec()),blockSize_);
          MVT::MvRandom(*randVecs);
        }

        { 	/* Ortho Timer */
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor lcltimer( *timerOrtho_ );
#endif
          rank = orthoMgr->normalize(*randVecs);
          if( rank < blockSize_ ) printer->stream(Warnings) << "Warning! Anasazi::RandomSolver Random vectors did not have full rank!" << std::endl;
        } 	/* End Ortho Timer */

        /* Set up variables for residual computation ------------------------------ */
        int i, j; 	// Loop variables 

        /* For computing H = Q^TAQ. (RR projection) */ 
        Teuchos::RCP<MV> TmpVecs = MVT::Clone(*randVecs,blockSize_);
        Teuchos::SerialDenseMatrix<OT,ScalarType> H (blockSize_, blockSize_);

        /* For solving the projected eigenvalue problem. */
        Teuchos::LAPACK<OT,ScalarType> lapack;
        Teuchos::SerialDenseMatrix<OT,ScalarType> evects (blockSize_, blockSize_);
        std::vector<MT> evals_real(blockSize_);
        std::vector<MT> evals_imag(blockSize_);

        /* Size of workspace and workspace for DGEEV */
        int info = -1000; 
        ScalarType* vlr = 0; 
        const int ldv = 1; 
        int lwork = -1;
        std::vector<ScalarType> work(1);
        std::vector<MT> rwork(2*blockSize_);
        numIters_ = 0;

        /* For computing the residuals of the eigenproblem */
        int numev;
        std::vector<Value<ScalarType>> EigenVals(blockSize_); 
        Teuchos::RCP<MV> Avecs, evecs;
        Teuchos::SerialDenseMatrix<OT,ScalarType> T (blockSize_, blockSize_ );
        std::vector<MT> normV( blockSize_ );
        bool converged = false;

        { // Solve Timer 
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor slvtimer(*timerSolve_);
#endif
          sol.numVecs = blockSize_;
          sol.Evals.resize(sol.numVecs);

          // Perform multiplies by A and Rayleigh-Ritz
          for( i = 0; i < maxIters_; i++ ){
            if (converged == true) {
              numFailed_ = 0;
              numIters_ = i-1;
              break;
            }

            { /* Begin Operator Timer */
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
              Teuchos::TimeMonitor lcltimer( *timerOp_ );
#endif
              OPT::Apply( *(problem_->getOperator()), *randVecs, *randVecs );
            } /* End Operator Timer */

            if ((orthoFreq_ > 0 && i % orthoFreq_ == 0) || (resFreq_ > 0 && i % resFreq_ == 0)) {
              { /* Start Ortho Timer */
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
                Teuchos::TimeMonitor lcltimer( *timerOrtho_ );
#endif
                rank = orthoMgr->normalize(*randVecs);
                if( rank < blockSize_ ) printer->stream(Warnings) << "Warning! Anasazi::RandomSolver Random vectors did not have full rank!" << std::endl;
              } /* End Ortho Timer */
            } /* End  Ortho */

            if (resFreq_ > 0 && i % resFreq_ == 0) {

              // Build the H matrix to run Rayleigh-Ritz on
              OPT::Apply( *(problem_->getOperator()), *randVecs, *TmpVecs ); 
              MVT::MvTransMv(ONE, *randVecs, *TmpVecs, H);

              // Run GEEV once to find the correct size for rwork
              lapack.GEEV('N','V',blockSize_,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);
              lwork = std::abs (static_cast<int> (Teuchos::ScalarTraits<ScalarType>::real (work[0])));
              work.resize( lwork );

              // Run GEEV a second time to solve for Harmonic Ritz Values:
              lapack.GEEV('N','V',blockSize_,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);
              if(info != 0) printer->stream(Warnings) << "Warning! Anasazi::RandomSolver GEEV solve possible failure: info = " << info << std::endl;
              lwork = -1;

              // sort the eigenvalues and permute the eigenvectors appropriately
              sorter->sort(evals_real,evals_imag,Teuchos::rcpFromRef(order),blockSize_);
              msutils.permuteVectors(order, evects);

              for( j = 0; j < blockSize_; j++){
                EigenVals[j].realpart = evals_real[j];
                EigenVals[j].imagpart = evals_imag[j];
              }

              // Project Evects back up to large problem. 
              MVT::MvTimesMatAddMv(ONE, *randVecs, evects, 0.0, *TmpVecs);

              // Copy only the eigenvectors we asked for to soln. 
              sol.Evecs = MVT::CloneCopy(*TmpVecs, Teuchos::Range1D(0,sol.numVecs-1));
              sol.numVecs = blockSize_;
              sol.Evals = EigenVals;

              // Extract evects/evals from solution
              evecs = sol.Evecs;
              numev = sol.numVecs;       

              // Check residuals for convergence
              if (numev > 0 ) {
                for ( j = 0; j < numev; ++j ) T(j, j) = sol.Evals[j].realpart;

                Avecs = MVT::Clone(*evecs, numev);
                OPT::Apply(*(problem_->getOperator()), *evecs, *Avecs);

                MVT::MvTimesMatAddMv(-ONE, *evecs, T, ONE, *Avecs);   /* Residuals = A*evecs - evecs*lambda */
                MVT::MvNorm(*Avecs, normV);

                numFailed_ = 0;
                converged = true;
                for ( j = 0; j < numev; ++j )
                {
                  if ( SCT::magnitude(sol.Evals[j].realpart) != SCT::zero() ) normV[j] = SCT::magnitude(normV[j]/sol.Evals[j].realpart);
                  if (normV[j] > tol_) {
                    numFailed_++;
                    converged = false;
                    break;
                  }
                }
              } // End residual computation              
            } // End Rayleigh-Ritz solve
          } // End subspace iterations

          if(converged == false)
          {
            { /* Begin Ortho Timer */
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
              Teuchos::TimeMonitor lcltimer( *timerOrtho_ );
#endif
              rank = orthoMgr->normalize(*randVecs);
              if( rank < blockSize_ ) printer->stream(Warnings) << "Warning! Anasazi::RandomSolver Random vectors did not have full rank!" << std::endl;
            } /* End Ortho Timer */

            OPT::Apply( *(problem_->getOperator()), *randVecs, *TmpVecs ); 
            MVT::MvTransMv(ONE, *randVecs, *TmpVecs, H);

            /* --------------------------------------------------------------------------
             * Parameters for DGEEV:
             *   'N'                 = Don't compute left eigenvectors. 
             *   'V'                 = Compute right eigenvectors. 
             *   blockSize           = Dimension of H (numEvals)
             *   H.values            = H matrix (Q'AQ)
             *   H.stride            = Leading dimension of H (numEvals)
             *   evals_real.data()   = Array to store evals, real parts
             *   evals_imag.data()   = Array to store evals, imag parts
             *   vlr                 = Stores left evects, so don't need this
             *   ldv                 = Leading dimension of vlr
             *   evects              = Array to store right evects
             *   evects.stride       = Leading dimension of evects
             *   work                = Work array
             *   lwork               = -1 means to query for array size
             *   rwork               = Not referenced because ST is not complex
             * -------------------------------------------------------------------------- */
            //Find workspace size for DGEEV:
            lapack.GEEV('N','V',blockSize_,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);
            if(info != 0) printer->stream(IterationDetails) << "Warning!! Anasazi::RandomSolver GEEV solve possible failure: info = " << info << std::endl;

            lwork = std::abs (static_cast<int> (Teuchos::ScalarTraits<ScalarType>::real (work[0])));
            work.resize( lwork );

            // Solve for Harmonic Ritz Values:
            lapack.GEEV('N','V',blockSize_,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);
            if(info != 0) printer->stream(IterationDetails) << "Warning!! Anasazi::RandomSolver GEEV solve possible failure: info = " << info << std::endl;

            // Sort the eigenvalues and permute the eigenvectors appropriately
            sorter->sort(evals_real,evals_imag,Teuchos::rcpFromRef(order),blockSize_);
            msutils.permuteVectors(order, evects);

            for( j = 0; j < blockSize_; j++){
              EigenVals[j].realpart = evals_real[j];
              EigenVals[j].imagpart = evals_imag[j];
            }
            sol.Evals = EigenVals;

            // Project Evects back up to large problem and permute
            MVT::MvTimesMatAddMv(ONE,*randVecs,evects, 0.0,*TmpVecs);

            //------Post-Solve Processing----------------------------
            //Copy only the eigenvectors we asked for to soln. 
            sol.numVecs = blockSize_;
            sol.Evecs = MVT::CloneCopy(*TmpVecs, Teuchos::Range1D(0,sol.numVecs-1));
            numIters_ = maxIters_;

          } // End converged == false
        } // End solve timer

        // Send the solution to the eigenproblem
        problem_->setSolution(sol);
        printer->stream(Debug) << "Returning " << sol.numVecs << " eigenpairs to eigenproblem." << std::endl;

        // print timing information
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        if ( printer->isVerbosity( TimingDetails ) ) Teuchos::TimeMonitor::summarize( printer->stream( TimingDetails ) );
#endif

        if (converged) return Converged;  // Return from RandomizedSolMgr::solve()
        return Unconverged; // Return from RandomizedSolMgr::solve()

      } // End Solve function
  } // end Experimental namespace
} // end Anasazi namespace

#endif /* ANASAZI_RANDOMIZED_SOLMGR_HPP */
