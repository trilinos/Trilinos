// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_TraceMinBase_SOLMGR_HPP
#define ANASAZI_TraceMinBase_SOLMGR_HPP

/*! \file AnasaziTraceMinBaseSolMgr.hpp
 *  \brief The Anasazi::TraceMinBaseSolMgr provides an abstract base class for the TraceMin series of solver managers.
*/

#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziOutputStreamTraits.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziICGSOrthoManager.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziSolverUtils.hpp"
#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziStatusTestOutput.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziStatusTestSpecTrans.hpp"
#include "AnasaziStatusTestWithOrdering.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziTraceMinBase.hpp"
#include "AnasaziTraceMinTypes.hpp"
#include "AnasaziTypes.hpp"

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_FancyOStream.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace Anasazi {
namespace Experimental {


/*! \class TraceMinBaseSolMgr
 *
 *  \brief The Anasazi::TraceMinBaseSolMgr provides an abstract base class for the TraceMin series of solver managers.
 *
 * This solver manager implements a hard-locking mechanism, whereby eigenpairs designated to be locked are moved from the eigensolver and placed in
 * auxilliary storage. The eigensolver is then restarted and continues to iterate, orthogonal to the locked eigenvectors.
 *
 * The solver manager provides to the solver a StatusTestCombo object constructed as follows:<br>
 *    &nbsp;&nbsp;&nbsp;<tt>combo = globaltest OR lockingtest OR debugtest</tt><br>
 * where
 *    - \c globaltest terminates computation when global convergence has been detected.<br>
 *      It is encapsulated in a StatusTestWithOrdering object, to ensure that computation is terminated
 *      only after the most significant eigenvalues/eigenvectors have met the convergence criteria.<br>
 *      If not specified via setGlobalStatusTest(), \c globaltest is a StatusTestResNorm object which tests the
 *      2-norms of the direct residuals relative to the Ritz values.
 *    - \c lockingtest halts TraceMinBase::iterate() in order to deflate converged eigenpairs for locking.<br>
 *      It will query the underlying TraceMinBase eigensolver to determine when eigenvectors should be locked.<br>
 *      If not specified via setLockingStatusTest(), \c lockingtest is a StatusTestResNorm object.
 *    - \c debugtest allows a user to specify additional monitoring of the iteration, encapsulated in a StatusTest object<br>
 *      If not specified via setDebugStatusTest(), \c debugtest is ignored.<br> 
 *      In most cases, it should return ::Failed; if it returns ::Passed, solve() will throw an AnasaziError exception.
 *
 * Additionally, the solver manager will terminate solve() after a specified number of restarts or iterations.
 * 
 * Much of this behavior is controlled via parameters and options passed to the
 * solver manager. For more information, see TraceMinBaseSolMgr().

 \ingroup anasazi_solver_framework

 \author Alicia Klinvex
 */

template<class ScalarType, class MV, class OP>
class TraceMinBaseSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;
    
  public:

  //! @name Constructors/Destructor
  //@{ 

  /*! \brief Basic constructor for TraceMinBaseSolMgr. 
   *
   * This constructor accepts the Eigenproblem to be solved in addition
   * to a parameter list of options for the solver manager. These options include the following:
   *   - Solver parameters
   *      - \c "Which" - a \c string that specifies whether we want the largest eigenvalues "LM" or the smallest "SM". Default: "SM"
   *      - \c "Verbosity" - a sum of MsgType specifying the verbosity. Default: ::Errors
   *      - \c "Output Stream" - a reference-counted pointer to the formatted output stream where all
   *                             solver output is sent.  Default: Teuchos::getFancyOStream ( Teuchos::rcpFromRef (std::cout) )
   *      - \c "Output Processor" - an \c int specifying the MPI processor that will print solver/timer details.  Default: 0
   *      - \c "Maximum Restarts" - a \c int specifying the maximum number of restarts the underlying solver is allowed to perform. Default: 20
   *      - \c "Saddle Solver Type" - a \c string specifying how to solve the saddle point problem arising at each iteration.
   *           Options are "Projected Krylov", "Schur Complement", and "Block Diagonal Preconditioned Minres". Default: "Projected Krylov"
   *            - \c "Projected Krylov": Uses projected-minres to solve the problem.
   *            - \c "Schur Complement": Explicitly forms the (inexact) Schur complement using minres. 
   *            - \c "Block Diagonal Preconditioned Minres": Uses a block preconditioner on the entire saddle point problem.  For more information, please see "Overview of Anasazi and its newest eigensolver, TraceMin" on the main Anasazi page.
   *           We recommend using "Projected Krylov" in the absence of preconditioning.  If you want to use a preconditioner, "Block Diagonal Preconditioned Minres" is recommended.
   *           "Schur Complement" mainly exists for special use cases.
   *      - Ritz shift parameters
   *         - \c "When To Shift" - a \c string specifying when Ritz shifts should be performed. Options are "Never", "After Trace Levels", and "Always". Default: "Always"
   *            - \c "Never": Do not perform Ritz shifts.  This option produces guaranteed convergence but converges linearly.  Not recommended.
   *            - \c "After Trace Levels": Do not perform Ritz shifts until the trace of \f$X^TKX\f$ has stagnated (i.e. the relative change in trace has become small).  
   *                 The \c MagnitudeType specifying how small the relative change in trace must become may be provided via the parameter \c "Trace Threshold", whose default value is 0.02.
   *            - \c "Always": Always attempt to use Ritz shifts.
   *         - \c "How To Choose Shift" - a \c string specifying how to choose the Ritz shifts (assuming Ritz shifts are being used).  
   *              Options are "Largest Converged", "Adjusted Ritz Values", and "Ritz Values". Default: "Adjusted Ritz Values"
   *            - \c "Largest Converged": Ritz shifts are chosen to be the largest converged eigenvalue.  Until an eigenvalue converges, the Ritz shifts are all 0.
   *            - \c "Adjusted Ritz Values": Ritz shifts are chosen based on the Ritz values and their associated residuals in such a way as to guarantee global convergence.  
   *                 This method is described in "The trace minimization method for the symmetric generalized eigenvalue problem."
   *            - \c "Ritz Values": Ritz shifts are chosen to equal the Ritz values.  This does NOT guarantee global convergence.
   *         - \c "Use Multiple Shifts" - a \c bool specifying whether to use one or many Ritz shifts (assuming shifting is enabled). Default: true
   *   - Convergence parameters (if using default convergence test; see setGlobalStatusTest())
   *      - \c "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: machine precision.
   *      - \c "Relative Convergence Tolerance" - a \c bool specifying whether residuals norms should be scaled by their eigenvalues for the purposing of deciding convergence. Default: true
   *      - \c "Convergence Norm" - a \c string specifying the norm for convergence testing: "2" or "M".  Default: "2"
   *   - Locking parameters (if using default locking test; see setLockingStatusTest())
   *      - \c "Use Locking" - a \c bool specifying whether the algorithm should employ locking of converged eigenpairs. Default: true
   *      - \c "Max Locked" - a \c int specifying the maximum number of eigenpairs to be locked. Default: problem->getNEV()
   *      - \c "Locking Quorum" - a \c int specifying the number of eigenpairs that must meet the locking criteria before locking actually occurs. Default: 1
   *      - \c "Locking Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide locking. Default: 0.1*convergence tolerance
   *      - \c "Relative Locking Tolerance" - a \c bool specifying whether residuals norms should be scaled by their eigenvalues for the purposing of deciding locking. Default: true
   *      - \c "Locking Norm" - a \c string specifying the norm for locking testing: "2" or "M". Default: "2" 
   *
   * Anasazi's trace minimization solvers are still in development, and we plan to add additional features in the future, including the ability to perform spectral transformations.
   */
  TraceMinBaseSolMgr( const RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
                             Teuchos::ParameterList &pl );

  //! Destructor.
  virtual ~TraceMinBaseSolMgr() {};
  //@}
  
  //! @name Accessor methods
  //@{ 

  //! Return the eigenvalue problem.
  const Eigenproblem<ScalarType,MV,OP>& getProblem() const {
    return *problem_;
  }

  //! Get the iteration count for the most recent call to \c solve().
  int getNumIters() const { 
    return numIters_; 
  }

  /*! \brief Return the timers for this object. 
   *
   * The timers are ordered as follows:
   *   - time spent in solve() routine
   *   - time spent restarting
   *   - time spent locking converged eigenvectors
   */
   Teuchos::Array<RCP<Teuchos::Time> > getTimers() const {
     return Teuchos::tuple(_timerSolve, _timerRestarting, _timerLocking);
   }

  //@}

  //! @name Solver application methods
  //@{ 
    
  /*! \brief This method performs possibly repeated calls to the underlying eigensolver's iterate() routine
   * until the problem has been solved (as decided by the solver manager) or the solver manager decides to 
   * quit.
   *
   * This method calls TraceMinBase::iterate(), which will return either because a specially constructed status test evaluates to ::Passed
   * or an exception is thrown.
   *
   * A return from TraceMinBase::iterate() signifies one of the following scenarios:
   *    - the maximum number of restarts/iterations has been exceeded. In this scenario, the solver manager will place\n
   *      all converged eigenpairs into the eigenproblem and return ::Unconverged.
   *    - the locking conditions have been met. In this scenario, some of the current eigenpairs will be removed\n
   *      from the eigensolver and placed into auxiliary storage. The eigensolver will be restarted with the remaining part of the Krylov subspace\n
   *      and some random information to replace the removed subspace.
   *    - the subspace is full, and we need to remove some vectors.  The eigensolver will be restarted with the most 
   *      significant part of the Krylov subspace.
   *    - global convergence has been met. In this case, the most significant NEV eigenpairs in the solver and locked storage  \n
   *      have met the convergence criterion. (Here, NEV refers to the number of eigenpairs requested by the Eigenproblem.)    \n
   *      In this scenario, the solver manager will return ::Converged.
   *
   * \returns ::ReturnType specifying:
   *     - ::Converged: the eigenproblem was solved to the specification required by the solver manager.
   *     - ::Unconverged: the eigenproblem was not solved to the specification desired by the solver manager.
  */
  ReturnType solve();

  //! Set the status test defining global convergence.
  void setGlobalStatusTest(const RCP< StatusTest<ScalarType,MV,OP> > &global);

  //! Get the status test defining global convergence.
  const RCP< StatusTest<ScalarType,MV,OP> > & getGlobalStatusTest() const;

  //! Set the status test defining locking.
  void setLockingStatusTest(const RCP< StatusTest<ScalarType,MV,OP> > &locking);

  //! Get the status test defining locking.
  const RCP< StatusTest<ScalarType,MV,OP> > & getLockingStatusTest() const;

  //! Set the status test for debugging.
  void setDebugStatusTest(const RCP< StatusTest<ScalarType,MV,OP> > &debug);

  //! Get the status test for debugging.
  const RCP< StatusTest<ScalarType,MV,OP> > & getDebugStatusTest() const;

  //@}

  protected:
  RCP<Eigenproblem<ScalarType,MV,OP> > problem_;

  int numIters_;

  // Block variables
  int blockSize_, numBlocks_, numRestartBlocks_;

  // Output variables
  RCP<OutputManager<ScalarType> > printer_;

  // Convergence variables
  MagnitudeType convTol_;
  bool relConvTol_;
  enum ResType convNorm_;

  // Locking variables
  MagnitudeType lockTol_;
  int maxLocked_, lockQuorum_;
  bool useLocking_, relLockTol_;
  enum ResType lockNorm_;

  // Shifting variables
  enum WhenToShiftType whenToShift_;
  MagnitudeType traceThresh_, shiftTol_;
  enum HowToShiftType howToShift_;
  bool useMultipleShifts_, relShiftTol_, considerClusters_;
  std::string shiftNorm_;
  
  // Other variables
  int maxKrylovIter_;
  std::string ortho_, which_;
  enum SaddleSolType saddleSolType_;
  bool projectAllVecs_, projectLockedVecs_, computeAllRes_, useRHSR_, useHarmonic_, noSort_;
  MagnitudeType alpha_;

  // Timers
  RCP<Teuchos::Time> _timerSolve, _timerRestarting, _timerLocking;

  // Status tests
  RCP<StatusTest<ScalarType,MV,OP> > globalTest_;
  RCP<StatusTest<ScalarType,MV,OP> > lockingTest_; 
  RCP<StatusTest<ScalarType,MV,OP> > debugTest_;

  // TraceMin specific functions
  void copyPartOfState(const TraceMinBaseState<ScalarType,MV>& oldState, TraceMinBaseState<ScalarType,MV>& newState, const std::vector<int> indToCopy) const;

  void setParameters(Teuchos::ParameterList &pl) const;

  void printParameters(std::ostream &os) const;

  virtual RCP< TraceMinBase<ScalarType,MV,OP> > createSolver( 
            const RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > &sorter,
            const RCP<StatusTest<ScalarType,MV,OP> >      &outputtest,
            const RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
            Teuchos::ParameterList &plist
          ) =0;

  virtual bool needToRestart(const RCP< TraceMinBase<ScalarType,MV,OP> > solver) =0;

  virtual bool performRestart(int &numRestarts, RCP< TraceMinBase<ScalarType,MV,OP> > solver) =0;
};


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor
template<class ScalarType, class MV, class OP>
TraceMinBaseSolMgr<ScalarType,MV,OP>::TraceMinBaseSolMgr( 
        const RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
        Teuchos::ParameterList &pl ) : 
  problem_(problem)
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  , _timerSolve(Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinBaseSolMgr::solve()")),
  _timerRestarting(Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinBaseSolMgr restarting")),
  _timerLocking(Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinBaseSolMgr locking"))
#endif
{
  TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null,              std::invalid_argument, "Problem not given to solver manager.");
  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isProblemSet(),              std::invalid_argument, "Problem not set.");
  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isHermitian(),               std::invalid_argument, "Problem not symmetric.");
  TEUCHOS_TEST_FOR_EXCEPTION(problem_->getInitVec() == Teuchos::null,std::invalid_argument, "Problem does not contain initial vectors to clone from.");

  std::string strtmp;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Block parameters

  // TODO: the default is different for TraceMin and TraceMin-Davidson
  // block size: default is nev()
//  blockSize_ = pl.get("Block Size",problem_->getNEV());
//  TEUCHOS_TEST_FOR_EXCEPTION(blockSize_ <= 0, std::invalid_argument,
//                     "Anasazi::TraceMinBaseSolMgr: \"Block Size\" must be strictly positive.");

  // TODO: add Num Blocks as a parameter to both child classes, since they have different default values
//  numBlocks_ = pl.get("Num Blocks",5);
//  TEUCHOS_TEST_FOR_EXCEPTION(numBlocks_ < 2, std::invalid_argument,
//                     "Anasazi::TraceMinBaseSolMgr: \"Num Blocks\" must be >= 2.");

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Output parameters

  // Create a formatted output stream to print to.
  // See if user requests output processor.
  int osProc = pl.get("Output Processor", 0);
     
  // If not passed in by user, it will be chosen based upon operator type.
  Teuchos::RCP<Teuchos::FancyOStream> osp;

  if (pl.isParameter("Output Stream")) {
    osp = Teuchos::getParameter<Teuchos::RCP<Teuchos::FancyOStream> >(pl,"Output Stream");
  }
  else {
    osp = OutputStreamTraits<OP>::getOutputStream (*problem_->getOperator(), osProc);
  }

  int verbosity = Anasazi::Errors;
  if (pl.isParameter("Verbosity")) {
    if (Teuchos::isParameterType<int>(pl,"Verbosity")) {
      verbosity = pl.get("Verbosity", verbosity);
    } else {
      verbosity = (int)Teuchos::getParameter<Anasazi::MsgType>(pl,"Verbosity");
    }
  }
  printer_ = rcp( new OutputManager<ScalarType>(verbosity,osp) );

  // TODO: Add restart parameters to TraceMin-Davidson

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Convergence parameters
  convTol_ = pl.get("Convergence Tolerance",MT::prec());
  TEUCHOS_TEST_FOR_EXCEPTION(convTol_ < 0, std::invalid_argument,
                     "Anasazi::TraceMinBaseSolMgr: \"Convergence Tolerance\" must be nonnegative.");

  relConvTol_ = pl.get("Relative Convergence Tolerance",true);
  strtmp = pl.get("Convergence Norm",std::string("2"));
  if (strtmp == "2") {
    convNorm_ = RES_2NORM;
  }
  else if (strtmp == "M") {
    convNorm_ = RES_ORTH;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
        "Anasazi::TraceMinBaseSolMgr: Invalid Convergence Norm.");
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Locking parameters
  useLocking_ = pl.get("Use Locking",true);
  relLockTol_ = pl.get("Relative Locking Tolerance",true);
  lockTol_ = pl.get("Locking Tolerance",convTol_/10);

  TEUCHOS_TEST_FOR_EXCEPTION(relConvTol_ != relLockTol_, std::invalid_argument,
         "Anasazi::TraceMinBaseSolMgr: \"Relative Convergence Tolerance\" and \"Relative Locking Tolerance\" have different values.  If you set one, you should always set the other.");

  TEUCHOS_TEST_FOR_EXCEPTION(lockTol_ < 0, std::invalid_argument,
                     "Anasazi::TraceMinBaseSolMgr: \"Locking Tolerance\" must be nonnegative.");

  strtmp = pl.get("Locking Norm",std::string("2"));
  if (strtmp == "2") {
    lockNorm_ = RES_2NORM;
  }
  else if (strtmp == "M") {
    lockNorm_ = RES_ORTH;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
        "Anasazi::TraceMinBaseSolMgr: Invalid Locking Norm.");
  }

  // max locked: default is nev(), must satisfy maxLocked_ + blockSize_ >= nev
  if (useLocking_) {
    maxLocked_ = pl.get("Max Locked",problem_->getNEV());
    TEUCHOS_TEST_FOR_EXCEPTION(maxLocked_ <= 0, std::invalid_argument,
           "Anasazi::TraceMinBaseSolMgr: \"Max Locked\" must be strictly positive.");
  }
  else {
    maxLocked_ = 0;
  }

  if (useLocking_) {
    lockQuorum_ = pl.get("Locking Quorum",1);
    TEUCHOS_TEST_FOR_EXCEPTION(lockQuorum_ <= 0, std::invalid_argument,
                       "Anasazi::TraceMinBaseSolMgr: \"Locking Quorum\" must be strictly positive.");
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Ritz shift parameters

  // When to shift - what triggers a shift?
  strtmp = pl.get("When To Shift", "Always");

  if(strtmp == "Never")
    whenToShift_ = NEVER_SHIFT;
  else if(strtmp == "After Trace Levels")
    whenToShift_ = SHIFT_WHEN_TRACE_LEVELS;
  else if(strtmp == "Residual Becomes Small")
    whenToShift_ = SHIFT_WHEN_RESID_SMALL;
  else if(strtmp == "Always")
    whenToShift_ = ALWAYS_SHIFT;
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
           "Anasazi::TraceMinBaseSolMgr: Invalid value for \"When To Shift\"; valid options are \"Never\", \"After Trace Levels\", \"Residual Becomes Small\", \"Always\".");  


  // How small does the % change in trace have to get before shifting?
  traceThresh_ = pl.get("Trace Threshold", 0.02);

  TEUCHOS_TEST_FOR_EXCEPTION(traceThresh_ < 0, std::invalid_argument,
                       "Anasazi::TraceMinBaseSolMgr: \"Trace Threshold\" must be nonnegative.");

  // Shift threshold - if the residual of an eigenpair is less than this, then shift
  shiftTol_ = pl.get("Shift Tolerance", 0.1);

  TEUCHOS_TEST_FOR_EXCEPTION(shiftTol_ < 0, std::invalid_argument,
                       "Anasazi::TraceMinBaseSolMgr: \"Shift Tolerance\" must be nonnegative.");

  // Use relative convergence tolerance - scale by eigenvalue?
  relShiftTol_ = pl.get("Relative Shift Tolerance", true);

  // Which norm to use in determining whether to shift
  shiftNorm_ = pl.get("Shift Norm", "2");

  TEUCHOS_TEST_FOR_EXCEPTION(shiftNorm_ != "2" && shiftNorm_ != "M", std::invalid_argument,
         "Anasazi::TraceMinBaseSolMgr: Invalid value for \"Shift Norm\"; valid options are \"2\" and \"M\".");  

  noSort_ = pl.get("No Sorting", false);

  // How to choose shift
  strtmp = pl.get("How To Choose Shift", "Adjusted Ritz Values");
  
  if(strtmp == "Largest Converged")
    howToShift_ = LARGEST_CONVERGED_SHIFT;
  else if(strtmp == "Adjusted Ritz Values")
    howToShift_ = ADJUSTED_RITZ_SHIFT;
  else if(strtmp == "Ritz Values")
    howToShift_ = RITZ_VALUES_SHIFT;
  else if(strtmp == "Experimental Shift")
    howToShift_ = EXPERIMENTAL_SHIFT;
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
           "Anasazi::TraceMinBaseSolMgr: Invalid value for \"How To Choose Shift\"; valid options are \"Largest Converged\", \"Adjusted Ritz Values\", \"Ritz Values\".");

  // Consider clusters - if all eigenvalues are in one cluster, it's not expecially safe to shift
  considerClusters_ = pl.get("Consider Clusters", true);

  // Use multiple shifts
  useMultipleShifts_ = pl.get("Use Multiple Shifts", true);

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Other parameters

  // which orthogonalization to use
  ortho_ = pl.get("Orthogonalization", "SVQB");
  TEUCHOS_TEST_FOR_EXCEPTION(ortho_ != "DGKS" && ortho_ != "SVQB" && ortho_ != "ICGS", std::invalid_argument,
         "Anasazi::TraceMinBaseSolMgr: Invalid value for \"Orthogonalization\"; valid options are \"DGKS\", \"SVQB\", \"ICGS\".");  

  strtmp = pl.get("Saddle Solver Type", "Projected Krylov");
  
  if(strtmp == "Projected Krylov")
    saddleSolType_ = PROJECTED_KRYLOV_SOLVER;
  else if(strtmp == "Schur Complement")
    saddleSolType_ = SCHUR_COMPLEMENT_SOLVER;
  else if(strtmp == "Block Diagonal Preconditioned Minres")
    saddleSolType_ = BD_PREC_MINRES;
  else if(strtmp == "HSS Preconditioned Gmres")
    saddleSolType_ = HSS_PREC_GMRES;
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
           "Anasazi::TraceMinBaseSolMgr: Invalid value for \"Saddle Solver Type\"; valid options are \"Projected Krylov\", \"Schur Complement\", and \"Block Diagonal Preconditioned Minres\".");

  projectAllVecs_ = pl.get("Project All Vectors", true);
  projectLockedVecs_ = pl.get("Project Locked Vectors", true);
  computeAllRes_ = pl.get("Compute All Residuals", true);
  useRHSR_ = pl.get("Use Residual as RHS", false);
  alpha_ = pl.get("HSS: alpha", 1.0);

  TEUCHOS_TEST_FOR_EXCEPTION(projectLockedVecs_ && ! projectAllVecs_, std::invalid_argument,
         "Anasazi::TraceMinBaseSolMgr: If you want to project out the locked vectors, you should really project out ALL the vectors of X.");

  // Maximum number of inner iterations
  maxKrylovIter_ = pl.get("Maximum Krylov Iterations", 200);
  TEUCHOS_TEST_FOR_EXCEPTION(maxKrylovIter_ < 1, std::invalid_argument,
         "Anasazi::TraceMinBaseSolMgr: \"Maximum Krylov Iterations\" must be greater than 0.");
		 
  // Which eigenvalues we want to get
  which_ = pl.get("Which", "SM");
  TEUCHOS_TEST_FOR_EXCEPTION(which_ != "SM" && which_ != "LM", std::invalid_argument,
         "Anasazi::TraceMinBaseSolMgr: Invalid value for \"Which\"; valid options are \"SM\" and \"LM\".");

  // Test whether we are shifting without an operator K
  // This is a really bad idea
  TEUCHOS_TEST_FOR_EXCEPTION(problem_->getOperator() == Teuchos::null && whenToShift_ != NEVER_SHIFT, std::invalid_argument,
         "Anasazi::TraceMinBaseSolMgr: It is an exceptionally bad idea to use Ritz shifts when finding the largest eigenpairs of a standard eigenvalue problem.  If we don't use Ritz shifts, it may take extra iterations to converge, but we NEVER have to solve a single linear system.  Using Ritz shifts forces us to solve systems of the form (I + sigma A)x=f, and it probably doesn't benefit us enough to outweigh the extra cost.  We may add support for this feature in the future, but for now, please set \"When To Shift\" to \"Never\".");

#ifdef BELOS_PSEUDO_BLOCK_GMRES_SOLMGR_HPP
  // Test whether we are using a projected preconditioner with multiple Ritz shifts
  // We can't currently do this for reasons that are complicated and are explained in the user manual
  TEUCHOS_TEST_FOR_EXCEPTION(problem_->getPrec() != Teuchos::null && saddleSolType_ == PROJECTED_KRYLOV_SOLVER && useMultipleShifts_, std::invalid_argument,
         "Anasazi::TraceMinBaseSolMgr: When you use the projected Krylov solver with preconditioning, the preconditioner must be projected as well.  In theory, if the preconditioner is SPD, the projected preconditioner will be SPSD, but in practice, it can have small negative eigenvalues, presumably due to machine arithmetic.  This means we can't use TraceMin's built-in MINRES, and we are forced to use Belos for now.  When you use multiple Ritz shifts, you are essentially using a different operator to solve each linear system.  Belos can't handle this right now, but we're working on a solution.  For now, please set \"Use Multiple Shifts\" to false.");
#else
  // Test whether we are using a projected preconditioner without Belos.
  // P Prec P should be positive definite if Prec is positive-definite, 
  // but it tends not to be in practice, presumably due to machine arithmetic
  // As a result, we have to use pseudo-block gmres for now.
  // Make sure it's available.
  TEUCHOS_TEST_FOR_EXCEPTION(problem_->getPrec() != Teuchos::null && saddleSolType_ == PROJECTED_KRYLOV_SOLVER, std::invalid_argument,
         "Anasazi::TraceMinBaseSolMgr: When you use the projected Krylov solver with preconditioning, the preconditioner must be projected as well.  In theory, if the preconditioner is SPD, the projected preconditioner will be SPSD, but in practice, it can have small negative eigenvalues, presumably due to machine arithmetic.  This means we can't use TraceMin's built-in MINRES, and we are forced to use Belos for now.  You didn't install Belos.  You have three options to correct this problem:\n1. Reinstall Trilinos with Belos enabled\n2. Don't use a preconditioner\n3. Choose a different method for solving the saddle-point problem (Recommended)");


#endif

  
}


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// solve()
template<class ScalarType, class MV, class OP>
ReturnType 
TraceMinBaseSolMgr<ScalarType,MV,OP>::solve() 
{
  typedef SolverUtils<ScalarType,MV,OP> msutils;

  const int nev = problem_->getNEV();

#ifdef TEUCHOS_DEBUG
    RCP<Teuchos::FancyOStream>
      out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(printer_->stream(Debug)));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
    *out << "Entering Anasazi::TraceMinBaseSolMgr::solve()\n";
#endif

  //////////////////////////////////////////////////////////////////////////////////////
  // Sort manager
  RCP<BasicSort<MagnitudeType> > sorter = rcp( new BasicSort<MagnitudeType>("SM") );

  //////////////////////////////////////////////////////////////////////////////////////
  // Handle the spectral transformation if necessary
  // TODO: Make sure we undo this before returning...
  if(which_ == "LM")
  {
    RCP<const OP> swapHelper = problem_->getOperator();
    problem_->setOperator(problem_->getM());
    problem_->setM(swapHelper);
    problem_->setProblem();
  } 

  //////////////////////////////////////////////////////////////////////////////////////
  // Status tests
  //
  // convergence
  RCP<StatusTest<ScalarType,MV,OP> > convtest;
  if (globalTest_ == Teuchos::null) {
    if(which_ == "SM")
      convtest = rcp( new StatusTestResNorm<ScalarType,MV,OP>(convTol_,nev,convNorm_,relConvTol_) );
    else
      convtest = rcp( new StatusTestSpecTrans<ScalarType,MV,OP>(convTol_,nev,convNorm_,relConvTol_,true,problem_->getOperator()) );
  }
  else {
    convtest = globalTest_;
  }
  RCP<StatusTestWithOrdering<ScalarType,MV,OP> > ordertest 
    = rcp( new StatusTestWithOrdering<ScalarType,MV,OP>(convtest,sorter,nev) );
  // locking
  RCP<StatusTest<ScalarType,MV,OP> > locktest;
  if (useLocking_) {
    if (lockingTest_ == Teuchos::null) {
      if(which_ == "SM")
        locktest = rcp( new StatusTestResNorm<ScalarType,MV,OP>(lockTol_,lockQuorum_,lockNorm_,relLockTol_) );
      else
        locktest = rcp( new StatusTestSpecTrans<ScalarType,MV,OP>(lockTol_,lockQuorum_,lockNorm_,relLockTol_,true,problem_->getOperator()) );
    }
    else {
      locktest = lockingTest_;
    }
  }
  // for a non-short-circuited OR test, the order doesn't matter
  Teuchos::Array<RCP<StatusTest<ScalarType,MV,OP> > > alltests;
  alltests.push_back(ordertest);
  if (locktest != Teuchos::null) alltests.push_back(locktest);
  if (debugTest_ != Teuchos::null) alltests.push_back(debugTest_);

  RCP<StatusTestCombo<ScalarType,MV,OP> > combotest
    = rcp( new StatusTestCombo<ScalarType,MV,OP>( StatusTestCombo<ScalarType,MV,OP>::OR, alltests) );
  // printing StatusTest
  RCP<StatusTestOutput<ScalarType,MV,OP> > outputtest;
  if ( printer_->isVerbosity(Debug) ) {
    outputtest = rcp( new StatusTestOutput<ScalarType,MV,OP>( printer_,combotest,1,Passed+Failed+Undefined ) );
  }
  else {
    outputtest = rcp( new StatusTestOutput<ScalarType,MV,OP>( printer_,combotest,1,Passed ) );
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Orthomanager
  RCP<MatOrthoManager<ScalarType,MV,OP> > ortho; 
  if (ortho_=="SVQB") {
    ortho = rcp( new SVQBOrthoManager<ScalarType,MV,OP>(problem_->getM()) );
  } else if (ortho_=="DGKS") {
    ortho = rcp( new BasicOrthoManager<ScalarType,MV,OP>(problem_->getM()) );
  } else if (ortho_=="ICGS") {
    ortho = rcp( new ICGSOrthoManager<ScalarType,MV,OP>(problem_->getM()) );
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Anasazi::TraceMinBaseSolMgr::solve(): Invalid orthogonalization type.");
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Parameter list
  Teuchos::ParameterList plist;
  setParameters(plist);

  //////////////////////////////////////////////////////////////////////////////////////
  // TraceMinBase solver
  RCP<TraceMinBase<ScalarType,MV,OP> > tm_solver 
    = createSolver(sorter,outputtest,ortho,plist);
  // set any auxiliary vectors defined in the problem
  RCP< const MV > probauxvecs = problem_->getAuxVecs();
  if (probauxvecs != Teuchos::null) {
    tm_solver->setAuxVecs( Teuchos::tuple< RCP<const MV> >(probauxvecs) );
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Storage
  // 
  // lockvecs will contain eigenvectors that have been determined "locked" by the status test
  int curNumLocked = 0;
  RCP<MV> lockvecs;
  // lockvecs is used to hold the locked eigenvectors, as well as for temporary storage when locking.
  // when locking, we will lock some number of vectors numnew, where numnew <= maxlocked - curlocked
  // we will produce numnew random vectors, which will go into the space with the new basis.
  // we will also need numnew storage for the image of these random vectors under A and M; 
  // columns [curlocked+1,curlocked+numnew] will be used for this storage
  if (maxLocked_ > 0) {
    lockvecs = MVT::Clone(*problem_->getInitVec(),maxLocked_);
  }
  std::vector<MagnitudeType> lockvals;
  //
  // Restarting occurs under two scenarios: when the basis is full and after locking.
  //
  // For the former, a new basis of size blockSize*numRestartBlocks is generated using the current basis
  // and the most significant primitive Ritz vectors (projected eigenvectors).
  //     [S,L] = eig(KK)
  //     S = [Sr St]   // some for "r"estarting, some are "t"runcated
  //     newV = V*Sr
  //     KK_new = newV'*K*newV = Sr'*V'*K*V*Sr = Sr'*KK*Sr
  //  Therefore, the only multivector operation needed is for the generation of newV.
  //
  //  * If the multiplication is explicit, it requires a workspace of blockSize*numRestartBlocks vectors. 
  //    This space must be specifically allocated for that task, as we don't have any space of that size.
  //    It (workMV) will be allocated at the beginning of solve()
  //  * Optionally, the multiplication can be performed implicitly, via a Householder QR factorization of 
  //    Sr. This can be done in situ, using the basis multivector contained in the solver. This requires
  //    that we cast away the const on the multivector returned from getState(). Workspace for this approach
  //    is a single vector. the solver's internal storage must be preserved (X,MX,KX,R), requiring us to 
  //    allocate this vector.
  //
  // For the latter (restarting after locking), the new basis is the same size as existing basis. If numnew
  // vectors are locked, they are deflated from the current basis and replaced with randomly generated 
  // vectors.
  //     [S,L] = eig(KK)
  //     S = [Sl Su]  // partitioned: "l"ocked and "u"nlocked
  //     newL = V*Sl = X(locked)
  //     defV = V*Su
  //     augV = rand(numnew)  // orthogonal to oldL,newL,defV,auxvecs
  //     newV = [defV augV]
  //     Kknew = newV'*K*newV = [Su'*KK*Su    defV'*K*augV]
  //                            [augV'*K*defV augV'*K*augV]
  //     locked = [oldL newL]
  // Clearly, this operation is more complicated than the previous.
  // Here is a list of the significant computations that need to be performed:
  // - newL will be put into space in lockvecs, but will be copied from getState().X at the end
  // - defV,augV will be stored in workspace the size of the current basis.
  // - M*augV and K*augV are needed; they will be stored in lockvecs. As a result, newL will
  //   not be put into lockvecs until the end.
  //
  // Therefore, we must allocate workMV when ((maxRestarts_ > 0) || (useLocking_ == true)) && inSituRestart == false
  // It will be allocated to size (numBlocks-1)*blockSize
  //

  // some consts and utils
  const ScalarType ONE = SCT::one();
  const ScalarType ZERO = SCT::zero();

  // go ahead and initialize the solution to nothing in case we throw an exception
  Eigensolution<ScalarType,MV> sol;
  sol.numVecs = 0;
  problem_->setSolution(sol);

  int numRestarts = 0;

  // enter solve() iterations
  {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*_timerSolve);
#endif

    // tell tm_solver to iterate
    while (1) {
      try {
        tm_solver->iterate();

        ////////////////////////////////////////////////////////////////////////////////////
        //
        //
        ////////////////////////////////////////////////////////////////////////////////////
        if (debugTest_ != Teuchos::null && debugTest_->getStatus() == Passed) {
          throw AnasaziError("Anasazi::TraceMinBaseSolMgr::solve(): User-specified debug status test returned Passed.");
        }
        ////////////////////////////////////////////////////////////////////////////////////
        //
        // check convergence next
        //
        ////////////////////////////////////////////////////////////////////////////////////
        else if (ordertest->getStatus() == Passed ) {
          // we have convergence
          // ordertest->whichVecs() tells us which vectors from lockvecs and solver state are the ones we want
          // ordertest->howMany() will tell us how many
          break;
        }
        ////////////////////////////////////////////////////////////////////////////////////
        //
        // check locking if we didn't converge or restart
        //
        ////////////////////////////////////////////////////////////////////////////////////
        else if (locktest != Teuchos::null && locktest->getStatus() == Passed) {

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor lcktimer(*_timerLocking);
#endif

          // 
          // get current state
          TraceMinBaseState<ScalarType,MV> state = tm_solver->getState();
          const int curdim = state.curDim;

          //
          // get number,indices of vectors to be locked
          TEUCHOS_TEST_FOR_EXCEPTION(locktest->howMany() <= 0,std::logic_error,
              "Anasazi::TraceMinBaseSolMgr::solve(): status test mistake: howMany() non-positive.");
          TEUCHOS_TEST_FOR_EXCEPTION(locktest->howMany() != (int)locktest->whichVecs().size(),std::logic_error,
              "Anasazi::TraceMinBaseSolMgr::solve(): status test mistake: howMany() not consistent with whichVecs().");
          // we should have room to lock vectors; otherwise, locking should have been deactivated
          TEUCHOS_TEST_FOR_EXCEPTION(curNumLocked == maxLocked_,std::logic_error,
              "Anasazi::TraceMinBaseSolMgr::solve(): status test mistake: locking not deactivated.");
          //
          // don't lock more than maxLocked_; we didn't allocate enough space.
          std::vector<int> tmp_vector_int;
          if (curNumLocked + locktest->howMany() > maxLocked_) {
            // just use the first of them
			for(int i=0; i<maxLocked_-curNumLocked; i++)
			  tmp_vector_int.push_back(locktest->whichVecs()[i]);
//            tmp_vector_int.insert(tmp_vector_int.begin(),locktest->whichVecs().begin(),locktest->whichVecs().begin()+maxLocked_-curNumLocked);
          }
          else {
            tmp_vector_int = locktest->whichVecs();
          }

          const std::vector<int> lockind(tmp_vector_int);
          const int numNewLocked = lockind.size();
          //
          // generate indices of vectors left unlocked
          // curind = [0,...,curdim-1] = UNION( lockind, unlockind )
          const int numUnlocked = curdim-numNewLocked;
          tmp_vector_int.resize(curdim);
          for (int i=0; i<curdim; i++) tmp_vector_int[i] = i;
          const std::vector<int> curind(tmp_vector_int);       // curind = [0 ... curdim-1]
          tmp_vector_int.resize(numUnlocked); 
          std::set_difference(curind.begin(),curind.end(),lockind.begin(),lockind.end(),tmp_vector_int.begin());
          const std::vector<int> unlockind(tmp_vector_int);    // unlockind = [0 ... curdim-1] - lockind
          tmp_vector_int.clear();

          //
          // debug printing
          if (printer_->isVerbosity(Debug)) {
            printer_->print(Debug,"Locking vectors: ");
            for (unsigned int i=0; i<lockind.size(); i++) {printer_->stream(Debug) << " " << lockind[i];}
            printer_->print(Debug,"\n");
            printer_->print(Debug,"Not locking vectors: ");
            for (unsigned int i=0; i<unlockind.size(); i++) {printer_->stream(Debug) << " " << unlockind[i];}
            printer_->print(Debug,"\n");
          }

          // Copy eigenvalues we want to lock into lockvals
          std::vector<Value<ScalarType> > allvals = tm_solver->getRitzValues();
          for(unsigned int i=0; i<allvals.size(); i++)
            printer_->stream(Debug) << "Ritz value[" << i << "] = " << allvals[i].realpart << std::endl;
          for (int i=0; i<numNewLocked; i++) {
            lockvals.push_back(allvals[lockind[i]].realpart);
          }

          // Copy vectors we want to lock into lockvecs
          RCP<const MV> newLocked = MVT::CloneView(*tm_solver->getRitzVectors(),lockind);
          std::vector<int> indlock(numNewLocked);
          for (int i=0; i<numNewLocked; i++) indlock[i] = curNumLocked+i;
          if(useHarmonic_)
          {
            RCP<MV> tempMV = MVT::CloneCopy(*newLocked);
            ortho->normalizeMat(*tempMV);
            MVT::SetBlock(*tempMV,indlock,*lockvecs);
          }
          else
          {
            MVT::SetBlock(*newLocked,indlock,*lockvecs);
          }

          // Tell the StatusTestWithOrdering that things have been locked
          // This is VERY important
          // If this set of lines is removed, the code does not terminate correctly
          if(noSort_)
          {
            for(unsigned int aliciaInd=0; aliciaInd<lockvals.size(); aliciaInd++)
            {
              lockvals[aliciaInd] = 0.0;
            }
          }
          ordertest->setAuxVals(lockvals);

          // Set the auxiliary vectors so that we remain orthogonal to the ones we locked
          // Remember to include any aux vecs provided by the user
          curNumLocked += numNewLocked;

          if(ordertest->getStatus() == Passed)  break;

          std::vector<int> curlockind(curNumLocked);
          for (int i=0; i<curNumLocked; i++) curlockind[i] = i;
          RCP<const MV> curlocked = MVT::CloneView(*lockvecs,curlockind);

          Teuchos::Array< RCP<const MV> > aux;
          if (probauxvecs != Teuchos::null) aux.push_back(probauxvecs);
          aux.push_back(curlocked);
          tm_solver->setAuxVecs(aux);

          // Disable locking if we have locked the maximum number of things
          printer_->stream(Debug) << "curNumLocked: " << curNumLocked << std::endl;
          printer_->stream(Debug) << "maxLocked: " << maxLocked_ << std::endl;
          if (curNumLocked == maxLocked_) {
            // disabled locking now
            combotest->removeTest(locktest);
            locktest = Teuchos::null;
            printer_->stream(Debug) << "Removed locking test\n";
          }

          int newdim = numRestartBlocks_*blockSize_;
          TraceMinBaseState<ScalarType,MV> newstate;
          if(newdim <= numUnlocked)
          {
            if(useHarmonic_)
            {
              std::vector<int> desiredSubscripts(newdim);
              for(int i=0; i<newdim; i++)
              {
                desiredSubscripts[i] = unlockind[i];
                printer_->stream(Debug) << "H desiredSubscripts[" << i << "] = " << desiredSubscripts[i] << std::endl;
              }
              newstate.V = MVT::CloneView(*tm_solver->getRitzVectors(),desiredSubscripts);
              newstate.curDim = newdim;
            }
            else
            {
              std::vector<int> desiredSubscripts(newdim);
              for(int i=0; i<newdim; i++)
              {
                desiredSubscripts[i] = unlockind[i];
                printer_->stream(Debug) << "desiredSubscripts[" << i << "] = " << desiredSubscripts[i] << std::endl;
              }

              copyPartOfState(state, newstate, desiredSubscripts);
            }
          }
          else
          {
            // TODO: Come back to this and make it more efficient

            // Replace the converged eigenvectors with random ones
            int nrandom = newdim-numUnlocked;
  
            RCP<const MV> helperMV;
            RCP<MV> totalV = MVT::Clone(*tm_solver->getRitzVectors(),newdim);

            // Holds old things that we're keeping
            tmp_vector_int.resize(numUnlocked);
            for(int i=0; i<numUnlocked; i++) tmp_vector_int[i] = i;
            RCP<MV> oldV = MVT::CloneViewNonConst(*totalV,tmp_vector_int);

            // Copy over the old things
            if(useHarmonic_)
              helperMV = MVT::CloneView(*tm_solver->getRitzVectors(),unlockind);
            else
              helperMV = MVT::CloneView(*state.V,unlockind);
            MVT::Assign(*helperMV,*oldV);

            // Holds random vectors we're generating
            tmp_vector_int.resize(nrandom);
            for(int i=0; i<nrandom; i++) tmp_vector_int[i] = i+numUnlocked;
            RCP<MV> newV = MVT::CloneViewNonConst(*totalV,tmp_vector_int);

            // Create random things
            MVT::MvRandom(*newV);

            newstate.V = totalV;
            newstate.curDim = newdim;

            // Copy Ritz shifts
            RCP< std::vector<ScalarType> > helperRS = rcp( new std::vector<ScalarType>(blockSize_) );
            for(unsigned int i=0; i<(unsigned int)blockSize_; i++) 
            {
              if(i < unlockind.size() && unlockind[i] < blockSize_)
                (*helperRS)[i] = (*state.ritzShifts)[unlockind[i]];
              else
                (*helperRS)[i] = ZERO;
            }
            newstate.ritzShifts  = helperRS;
          }

          // Determine the largest safe shift
          newstate.largestSafeShift = std::abs(lockvals[0]);
          for(size_t i=0; i<lockvals.size(); i++)
            newstate.largestSafeShift = std::max(std::abs(lockvals[i]), newstate.largestSafeShift);

          // Prepare new state, removing converged vectors
          // TODO: Init will perform some unnecessary calculations; look into it
          // TODO: The residual norms should be part of the state
          newstate.NEV = state.NEV - numNewLocked;
          tm_solver->initialize(newstate);
        } // end of locking
        ////////////////////////////////////////////////////////////////////////////////////
        //
        // check for restarting before locking: if we need to lock, it will happen after the restart
        //
        ////////////////////////////////////////////////////////////////////////////////////
        else if ( needToRestart(tm_solver) ) {
          // if performRestart returns false, we exceeded the maximum number of restarts
          if(performRestart(numRestarts, tm_solver) == false)
            break;
        } // end of restarting
        ////////////////////////////////////////////////////////////////////////////////////
        //
        // we returned from iterate(), but none of our status tests Passed.
        // something is wrong, and it is probably our fault.
        //
        ////////////////////////////////////////////////////////////////////////////////////
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Anasazi::TraceMinBaseSolMgr::solve(): Invalid return from tm_solver::iterate().");
        }
      }
      catch (const AnasaziError &err) {
        printer_->stream(Errors) 
          << "Anasazi::TraceMinBaseSolMgr::solve() caught unexpected exception from Anasazi::TraceMinBase::iterate() at iteration " << tm_solver->getNumIters() << std::endl
          << err.what() << std::endl
          << "Anasazi::TraceMinBaseSolMgr::solve() returning Unconverged with no solutions." << std::endl;
        return Unconverged;
      }
    }

    sol.numVecs = ordertest->howMany();
    if (sol.numVecs > 0) {
      sol.Evecs = MVT::Clone(*problem_->getInitVec(),sol.numVecs);
      sol.Espace = sol.Evecs;
      sol.Evals.resize(sol.numVecs);
      std::vector<MagnitudeType> vals(sol.numVecs);

      // copy them into the solution
      std::vector<int> which = ordertest->whichVecs();
      // indices between [0,blockSize) refer to vectors/values in the solver
      // indices between [-curNumLocked,-1] refer to locked vectors/values [0,curNumLocked)
      // everything has already been ordered by the solver; we just have to partition the two references
      std::vector<int> inlocked(0), insolver(0);
      for (unsigned int i=0; i<which.size(); i++) {
        if (which[i] >= 0) {
          TEUCHOS_TEST_FOR_EXCEPTION(which[i] >= blockSize_,std::logic_error,"Anasazi::TraceMinBaseSolMgr::solve(): positive indexing mistake from ordertest.");
          insolver.push_back(which[i]);
        }
        else {
          // sanity check
          TEUCHOS_TEST_FOR_EXCEPTION(which[i] < -curNumLocked,std::logic_error,"Anasazi::TraceMinBaseSolMgr::solve(): negative indexing mistake from ordertest.");
          inlocked.push_back(which[i] + curNumLocked);
        }
      }

      TEUCHOS_TEST_FOR_EXCEPTION(insolver.size() + inlocked.size() != (unsigned int)sol.numVecs,std::logic_error,"Anasazi::TraceMinBaseSolMgr::solve(): indexing mistake.");

      // set the vecs,vals in the solution
      if (insolver.size() > 0) {
        // set vecs
        int lclnum = insolver.size();
        std::vector<int> tosol(lclnum);
        for (int i=0; i<lclnum; i++) tosol[i] = i;
        RCP<const MV> v = MVT::CloneView(*tm_solver->getRitzVectors(),insolver);
        MVT::SetBlock(*v,tosol,*sol.Evecs);
        // set vals
        std::vector<Value<ScalarType> > fromsolver = tm_solver->getRitzValues();
        for (unsigned int i=0; i<insolver.size(); i++) {
          vals[i] = fromsolver[insolver[i]].realpart;
        }
      }

      // get the vecs,vals from locked storage
      if (inlocked.size() > 0) {
        int solnum = insolver.size();
        // set vecs
        int lclnum = inlocked.size();
        std::vector<int> tosol(lclnum);
        for (int i=0; i<lclnum; i++) tosol[i] = solnum + i;
        RCP<const MV> v = MVT::CloneView(*lockvecs,inlocked);
        MVT::SetBlock(*v,tosol,*sol.Evecs);
        // set vals
        for (unsigned int i=0; i<inlocked.size(); i++) {
          vals[i+solnum] = lockvals[inlocked[i]];
        }
      }

      // undo the spectral transformation if necessary
      // if we really passed the solver Bx = \lambda A x, invert the eigenvalues
      if(which_ == "LM")
      {
        for(size_t i=0; i<vals.size(); i++)
          vals[i] = ONE/vals[i];
      }

      // sort the eigenvalues and permute the eigenvectors appropriately
      {
        std::vector<int> order(sol.numVecs);
        sorter->sort(vals,Teuchos::rcpFromRef(order),sol.numVecs);
        // store the values in the Eigensolution
        for (int i=0; i<sol.numVecs; i++) {
          sol.Evals[i].realpart = vals[i];
          sol.Evals[i].imagpart = MT::zero();
        }
        // now permute the eigenvectors according to order
        msutils::permuteVectors(sol.numVecs,order,*sol.Evecs);
      }

      // setup sol.index, remembering that all eigenvalues are real so that index = {0,...,0}
      sol.index.resize(sol.numVecs,0);
    }
  }

  // print final summary
  tm_solver->currentStatus(printer_->stream(FinalSummary));

  printParameters(printer_->stream(FinalSummary));

  // print timing information
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  if ( printer_->isVerbosity( TimingDetails ) ) {
    Teuchos::TimeMonitor::summarize( printer_->stream( TimingDetails ) );
  }
#endif

  problem_->setSolution(sol);
  printer_->stream(Debug) << "Returning " << sol.numVecs << " eigenpairs to eigenproblem." << std::endl;

  // get the number of iterations taken for this call to solve().
  numIters_ = tm_solver->getNumIters();

  if (sol.numVecs < nev) {
    return Unconverged; // return from TraceMinBaseSolMgr::solve() 
  }
  return Converged; // return from TraceMinBaseSolMgr::solve() 
}


template <class ScalarType, class MV, class OP>
void 
TraceMinBaseSolMgr<ScalarType,MV,OP>::setGlobalStatusTest(
    const RCP< StatusTest<ScalarType,MV,OP> > &global) 
{
  globalTest_ = global;
}

template <class ScalarType, class MV, class OP>
const RCP< StatusTest<ScalarType,MV,OP> > & 
TraceMinBaseSolMgr<ScalarType,MV,OP>::getGlobalStatusTest() const 
{
  return globalTest_;
}

template <class ScalarType, class MV, class OP>
void 
TraceMinBaseSolMgr<ScalarType,MV,OP>::setDebugStatusTest(
    const RCP< StatusTest<ScalarType,MV,OP> > &debug)
{
  debugTest_ = debug;
}

template <class ScalarType, class MV, class OP>
const RCP< StatusTest<ScalarType,MV,OP> > & 
TraceMinBaseSolMgr<ScalarType,MV,OP>::getDebugStatusTest() const
{
  return debugTest_;
}

template <class ScalarType, class MV, class OP>
void 
TraceMinBaseSolMgr<ScalarType,MV,OP>::setLockingStatusTest(
    const RCP< StatusTest<ScalarType,MV,OP> > &locking) 
{
  lockingTest_ = locking;
}

template <class ScalarType, class MV, class OP>
const RCP< StatusTest<ScalarType,MV,OP> > & 
TraceMinBaseSolMgr<ScalarType,MV,OP>::getLockingStatusTest() const 
{
  return lockingTest_;
}

template <class ScalarType, class MV, class OP>
void TraceMinBaseSolMgr<ScalarType,MV,OP>::copyPartOfState(const TraceMinBaseState<ScalarType,MV>& oldState, TraceMinBaseState<ScalarType,MV>& newState, const std::vector<int> indToCopy) const
{
  const ScalarType ONE = Teuchos::ScalarTraits<MagnitudeType>::one();
  const ScalarType ZERO = Teuchos::ScalarTraits<MagnitudeType>::zero();

  newState.curDim = indToCopy.size();
  std::vector<int> fullIndices(oldState.curDim);
  for(int i=0; i<oldState.curDim; i++) fullIndices[i] = i;

  // Initialize with X.  
  // Note that we didn't compute enough vectors of X, but we can very easily using the Ritz vectors.
  // That's why they're part of the state.
  // Note that there will ALWAYS be enough vectors

  // Helpful vectors for computing views and whatnot
  std::vector<int> oldIndices;
  std::vector<int> newIndices;
  for(int i=0; i<newState.curDim; i++)
  {
    if(indToCopy[i] < blockSize_)
      oldIndices.push_back(indToCopy[i]);
    else
      newIndices.push_back(indToCopy[i]);
  }

  int olddim = oldIndices.size();
  int newdim = newIndices.size();

  // If there are no new vectors being copied
  if(computeAllRes_)
  {
    newState.V  = MVT::CloneView(*oldState.X, indToCopy);
    newState.R  = MVT::CloneView(*oldState.R, indToCopy);
    newState.X = newState.V;

    if(problem_->getOperator() != Teuchos::null)
    {
      newState.KV = MVT::CloneView(*oldState.KX, indToCopy);
      newState.KX = newState.KV;
    }
    else
    {
      newState.KV = Teuchos::null;
      newState.KX = Teuchos::null;
    }

    if(problem_->getM() != Teuchos::null)
    {
      newState.MopV = MVT::CloneView(*oldState.MX, indToCopy);
      newState.MX = newState.MopV;
    }
    else
    {
      newState.MopV = Teuchos::null;
      newState.MX = Teuchos::null;
    }
  }
  else if(newdim == 0)
  {
    std::vector<int> blockind(blockSize_);
    for(int i=0; i<blockSize_; i++)
      blockind[i] = i;

    // Initialize with X
    newState.V  = MVT::CloneView(*oldState.X, blockind);
    newState.KV = MVT::CloneView(*oldState.KX, blockind);
    newState.R  = MVT::CloneView(*oldState.R, blockind);
    newState.X = MVT::CloneView(*newState.V, blockind);
    newState.KX = MVT::CloneView(*newState.KV, blockind);

    if(problem_->getM() != Teuchos::null)
    {
      newState.MopV = MVT::CloneView(*oldState.MX, blockind);
      newState.MX = MVT::CloneView(*newState.MopV, blockind);
    }
    else
    {
      newState.MopV = Teuchos::null;
      newState.MX = Teuchos::null;
    }
  }
  else
  {
    // More helpful vectors
    std::vector<int> oldPart(olddim);
    for(int i=0; i<olddim; i++) oldPart[i] = i;
    std::vector<int> newPart(newdim);
    for(int i=0; i<newdim; i++) newPart[i] = olddim+i;

    // Helpful multivectors for views and whatnot
    RCP<MV> helper = MVT::Clone(*oldState.V,newState.curDim);
    RCP<MV> oldHelper = MVT::CloneViewNonConst(*helper,oldPart);
    RCP<MV> newHelper = MVT::CloneViewNonConst(*helper,newPart);
    RCP<const MV> viewHelper;

    // Get the parts of the Ritz vectors we are interested in.
    Teuchos::SerialDenseMatrix<int,ScalarType> newRV(oldState.curDim,newdim);
    for(int r=0; r<oldState.curDim; r++)
    {
      for(int c=0; c<newdim; c++)
        newRV(r,c) = (*oldState.RV)(r,newIndices[c]);
    }

    // We're going to compute X as V*RitzVecs
    viewHelper = MVT::CloneView(*oldState.V,fullIndices);
    MVT::MvTimesMatAddMv(ONE,*viewHelper,newRV,ZERO,*newHelper);
    viewHelper = MVT::CloneView(*oldState.X,oldIndices);
    MVT::Assign(*viewHelper,*oldHelper);
    newState.V = MVT::CloneCopy(*helper);

    // Also compute KX as KV*RitzVecs
    viewHelper = MVT::CloneView(*oldState.KV,fullIndices);
    MVT::MvTimesMatAddMv(ONE,*viewHelper,newRV,ZERO,*newHelper);
    viewHelper = MVT::CloneView(*oldState.KX,oldIndices);
    MVT::Assign(*viewHelper,*oldHelper);
    newState.KV = MVT::CloneCopy(*helper);

    // Do the same with MX if necessary
    if(problem_->getM() != Teuchos::null)
    {
      viewHelper = MVT::CloneView(*oldState.MopV,fullIndices);
      MVT::MvTimesMatAddMv(ONE,*viewHelper,newRV,ZERO,*newHelper);
      viewHelper = MVT::CloneView(*oldState.MX,oldIndices);
      MVT::Assign(*viewHelper,*oldHelper);
      newState.MopV = MVT::CloneCopy(*helper);
    }
    else
      newState.MopV = newState.V;

    // Get X, MX, KX
    std::vector<int> blockVec(blockSize_);
    for(int i=0; i<blockSize_; i++) blockVec[i] = i;
    newState.X = MVT::CloneView(*newState.V,blockVec);
    newState.KX = MVT::CloneView(*newState.KV,blockVec);
    newState.MX = MVT::CloneView(*newState.MopV,blockVec);

    // Update the residuals
    if(blockSize_-oldIndices.size() > 0)
    {
      // There are vectors we have not computed the residual for yet
      newPart.resize(blockSize_-oldIndices.size());
      helper = MVT::Clone(*oldState.V,blockSize_);
      oldHelper = MVT::CloneViewNonConst(*helper,oldPart);
      newHelper = MVT::CloneViewNonConst(*helper,newPart);

      RCP<MV> scaledMV = MVT::CloneCopy(*newState.MX,newPart);
      RCP<const MV> localKV = MVT::CloneView(*newState.KX,newPart);
      std::vector<ScalarType> scalarVec(blockSize_-oldIndices.size());
      for(unsigned int i=0; i<(unsigned int)blockSize_-oldIndices.size(); i++) scalarVec[i] = (*oldState.T)[newPart[i]];
      MVT::MvScale(*scaledMV,scalarVec);
            
      helper = MVT::Clone(*oldState.V,blockSize_);
      oldHelper = MVT::CloneViewNonConst(*helper,oldPart);
      newHelper = MVT::CloneViewNonConst(*helper,newPart);
      MVT::MvAddMv(ONE,*localKV,-ONE,*scaledMV,*newHelper);
      viewHelper = MVT::CloneView(*oldState.R,oldIndices);
      MVT::Assign(*viewHelper,*oldHelper);
      newState.R = MVT::CloneCopy(*helper);
    }
    else
      newState.R = oldState.R;
  }

  // Since we are setting V:=X, V is orthonormal
  newState.isOrtho = true;

  // Get the first eigenvalues
  RCP< std::vector<ScalarType> > helperT = rcp( new std::vector<ScalarType>(newState.curDim) );
  for(int i=0; i<newState.curDim; i++) (*helperT)[i] = (*oldState.T)[indToCopy[i]];
  newState.T  = helperT;

  // X'KX is diag(T)
  RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > newKK = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newState.curDim,newState.curDim) );
  for(int i=0; i<newState.curDim; i++)
    (*newKK)(i,i) = (*newState.T)[i];
  newState.KK = newKK;

  // The associated Ritz vectors are I
  RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > newRV = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newState.curDim,newState.curDim) );
  for(int i=0; i<newState.curDim; i++)
    (*newRV)(i,i) = ONE;
  newState.RV = newRV;

  // Get the Ritz shifts
  RCP< std::vector<ScalarType> > helperRS = rcp( new std::vector<ScalarType>(blockSize_) );
  for(int i=0; i<blockSize_; i++) 
  {
    if(indToCopy[i] < blockSize_)
      (*helperRS)[i] = (*oldState.ritzShifts)[indToCopy[i]];
    else
      (*helperRS)[i] = ZERO;
  }
  newState.ritzShifts  = helperRS;
}


template <class ScalarType, class MV, class OP>
void TraceMinBaseSolMgr<ScalarType,MV,OP>::setParameters(Teuchos::ParameterList &pl) const
{
  pl.set("Block Size", blockSize_);
  pl.set("Num Blocks", numBlocks_);
  pl.set("Num Restart Blocks", numRestartBlocks_);
  pl.set("When To Shift", whenToShift_);
  pl.set("Trace Threshold", traceThresh_);
  pl.set("Shift Tolerance", shiftTol_);
  pl.set("Relative Shift Tolerance", relShiftTol_);
  pl.set("Shift Norm", shiftNorm_);
  pl.set("How To Choose Shift", howToShift_);
  pl.set("Consider Clusters", considerClusters_);
  pl.set("Use Multiple Shifts", useMultipleShifts_);
  pl.set("Saddle Solver Type", saddleSolType_);
  pl.set("Project All Vectors", projectAllVecs_);
  pl.set("Project Locked Vectors", projectLockedVecs_);
  pl.set("Compute All Residuals", computeAllRes_);
  pl.set("Use Residual as RHS", useRHSR_);
  pl.set("Use Harmonic Ritz Values", useHarmonic_);
  pl.set("Maximum Krylov Iterations", maxKrylovIter_);
  pl.set("HSS: alpha", alpha_);
}


template <class ScalarType, class MV, class OP>
void TraceMinBaseSolMgr<ScalarType,MV,OP>::printParameters(std::ostream &os) const
{
  os << "\n\n\n";
  os << "========================================\n";
  os << "========= TraceMin parameters ==========\n";
  os << "========================================\n";
  os << "=========== Block parameters ===========\n";
  os << "Block Size: " << blockSize_ << std::endl;
  os << "Num Blocks: " << numBlocks_ << std::endl;
  os << "Num Restart Blocks: " << numRestartBlocks_ << std::endl;
  os << "======== Convergence parameters ========\n";
  os << "Convergence Tolerance: " << convTol_ << std::endl;
  os << "Relative Convergence Tolerance: " << relConvTol_ << std::endl;
  os << "========== Locking parameters ==========\n";
  os << "Use Locking: " << useLocking_ << std::endl;
  os << "Locking Tolerance: " << lockTol_ << std::endl;
  os << "Relative Locking Tolerance: " << relLockTol_ << std::endl;
  os << "Max Locked: " << maxLocked_ << std::endl;
  os << "Locking Quorum: " << lockQuorum_ << std::endl;
  os << "========== Shifting parameters =========\n";
  os << "When To Shift: ";
  if(whenToShift_ == NEVER_SHIFT) os << "Never\n";
  else if(whenToShift_ == SHIFT_WHEN_TRACE_LEVELS) os << "After Trace Levels\n";
  else if(whenToShift_ == SHIFT_WHEN_RESID_SMALL) os << "Residual Becomes Small\n";
  else if(whenToShift_ == ALWAYS_SHIFT) os << "Always\n";
  os << "Consider Clusters: " << considerClusters_ << std::endl;
  os << "Trace Threshohld: " << traceThresh_ << std::endl;
  os << "Shift Tolerance: " << shiftTol_ << std::endl;
  os << "Relative Shift Tolerance: " << relShiftTol_ << std::endl;
  os << "How To Choose Shift: ";
  if(howToShift_ == LARGEST_CONVERGED_SHIFT) os << "Largest Converged\n";
  else if(howToShift_ == ADJUSTED_RITZ_SHIFT) os << "Adjusted Ritz Values\n";
  else if(howToShift_ == RITZ_VALUES_SHIFT) os << "Ritz Values\n";
  os << "Use Multiple Shifts: " << useMultipleShifts_ << std::endl;
  os << "=========== Other parameters ===========\n";
  os << "Orthogonalization: " << ortho_ << std::endl;
  os << "Saddle Solver Type: ";
  if(saddleSolType_ == PROJECTED_KRYLOV_SOLVER) os << "Projected Krylov\n";
  else if(saddleSolType_ == SCHUR_COMPLEMENT_SOLVER) os << "Schur Complement\n";
  os << "Project All Vectors: " << projectAllVecs_ << std::endl;
  os << "Project Locked Vectors: " << projectLockedVecs_ << std::endl;
  os << "Compute All Residuals: " << computeAllRes_ << std::endl;
  os << "========================================\n\n\n";
}


}} // end Anasazi namespace

#endif /* ANASAZI_TraceMinBase_SOLMGR_HPP */
