// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_RTR_SOLMGR_HPP
#define ANASAZI_RTR_SOLMGR_HPP

/*! \file AnasaziRTRSolMgr.hpp
  \brief The Anasazi::RTRSolMgr provides a simple solver manager over the IRTR
  eigensolvers.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziSolverUtils.hpp"

#include "AnasaziIRTR.hpp"
#include "AnasaziSIRTR.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziICGSOrthoManager.hpp"
#include "AnasaziStatusTestMaxIters.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziStatusTestWithOrdering.hpp"
#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziStatusTestOutput.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziOutputStreamTraits.hpp"

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_FancyOStream.hpp"

/*! \class Anasazi::RTRSolMgr
  \brief The Anasazi::RTRSolMgr provides a simple solver
  manager over the RTR eigensolver. For more information, see the discussion for RTRBase.

  \ingroup anasazi_solver_framework

  \author Chris Baker
*/

namespace Anasazi {

template<class ScalarType, class MV, class OP>
class RTRSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;

  public:

  //! @name Constructors/Destructor
  //@{

  /*! \brief Basic constructor for RTRSolMgr.
   *
   * This constructor accepts the Eigenproblem to be solved in addition
   * to a parameter list of options for the solver manager. These options include the following:
   *   - Solver parameters
   *      - \c "Skinny Solver" - a \c bool specifying whether a non-caching ("skinny") solver implementation is used. Determines whether the underlying solver is IRTR or SIRTR.
   *      - \c "Which" - a \c string specifying the desired eigenvalues: SR or LR, i.e., smallest or largest algebraic eigenvalues.
   *      - \c "Block Size" - an \c int specifying the block size to be used by the underlying RTR solver. Default: problem->getNEV()
   *      - \c "Verbosity" - a sum of MsgType specifying the verbosity. Default: ::Errors
   *      - \c "Output Stream" - a reference-counted pointer to the formatted output stream where all
   *                             solver output is sent.  Default: Teuchos::getFancyOStream ( Teuchos::rcpFromRef (std::cout) )
   *      - \c "Output Processor" - an \c int specifying the MPI processor that will print solver/timer details.  Default: 0
   *   - Convergence parameters
   *      - \c "Maximum Iterations" - an \c int specifying the maximum number of iterations the underlying solver is allowed to perform. Default: 100
   *      - \c "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: machine precision.
   *      - \c "Relative Convergence Tolerance" - a \c bool specifying whether residuals norms should be scaled by their eigenvalues for the purposing of deciding convergence. Default: true
   *      - \c "Convergence Norm" - a \c string specifying the norm for convergence testing: "2" or "M"
   */
  RTRSolMgr( const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
             Teuchos::ParameterList &pl );

  //! Destructor.
  virtual ~RTRSolMgr() {};
  //@}

  //! @name Accessor methods
  //@{

  //! Return the eigenvalue problem.
  const Eigenproblem<ScalarType,MV,OP>& getProblem() const {
    return *problem_;
  }

  /*! \brief Return the timers for this object.
   *
   * The timers are ordered as follows:
   *   - time spent in solve() routine
   */
   Teuchos::Array<Teuchos::RCP<Teuchos::Time> > getTimers() const {
     return Teuchos::tuple(_timerSolve);
   }

  //! Get the iteration count for the most recent call to solve.
  int getNumIters() const {
    return numIters_;
  }


  //@}

  //! @name Solver application methods
  //@{

  /*! \brief This method performs possibly repeated calls to the underlying eigensolver's iterate() routine
   * until the problem has been solved (as decided by the solver manager) or the solver manager decides to
   * quit.
   *
   * \returns ::ReturnType specifying:
   *     - ::Converged: the eigenproblem was solved to the specification required by the solver manager.
   *     - ::Unconverged: the eigenproblem was not solved to the specification desired by the solver manager.
  */
  ReturnType solve();
  //@}

  private:
  Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > problem_;
  std::string whch_;
  std::string ortho_;

  bool skinny_;
  MagnitudeType convtol_;
  int maxIters_;
  bool relconvtol_;
  enum ResType convNorm_;
  int numIters_;
  int numICGS_;
  int blkSize_;

  Teuchos::RCP<Teuchos::Time> _timerSolve;
  Teuchos::RCP<OutputManager<ScalarType> > printer_;
  Teuchos::ParameterList& pl_;
};


////////////////////////////////////////////////////////////////////////////////////////
template<class ScalarType, class MV, class OP>
RTRSolMgr<ScalarType,MV,OP>::RTRSolMgr(
        const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
        Teuchos::ParameterList &pl ) :
  problem_(problem),
  whch_("SR"),
  skinny_(true),
  convtol_(MT::prec()),
  maxIters_(100),
  relconvtol_(true),
  numIters_(-1),
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  _timerSolve(Teuchos::TimeMonitor::getNewTimer("Anasazi: RTRSolMgr::solve()")),
#endif
  pl_(pl)
{
  TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null,              std::invalid_argument, "Problem not given to solver manager.");
  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isProblemSet(),              std::invalid_argument, "Problem not set.");
  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isHermitian(),               std::invalid_argument, "Problem not symmetric.");
  TEUCHOS_TEST_FOR_EXCEPTION(problem_->getInitVec() == Teuchos::null,std::invalid_argument, "Problem does not contain initial vectors to clone from.");

  std::string strtmp;

  whch_ = pl_.get("Which","SR");
  TEUCHOS_TEST_FOR_EXCEPTION(whch_ != "SR" && whch_ != "LR",
      std::invalid_argument, "Anasazi::RTRSolMgr: Invalid sorting string. RTR solvers compute only LR or SR.");

  // convergence tolerance
  convtol_ = pl_.get("Convergence Tolerance",convtol_);
  relconvtol_ = pl_.get("Relative Convergence Tolerance",relconvtol_);
  strtmp = pl_.get("Convergence Norm",std::string("2"));
  if (strtmp == "2") {
    convNorm_ = RES_2NORM;
  }
  else if (strtmp == "M") {
    convNorm_ = RES_ORTH;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "Anasazi::RTRSolMgr: Invalid Convergence Norm.");
  }


  // maximum number of (outer) iterations
  maxIters_ = pl_.get("Maximum Iterations",maxIters_);

  // skinny solver or not
  skinny_ = pl_.get("Skinny Solver",skinny_);

  // number if ICGS iterations
  numICGS_ = pl_.get("Num ICGS",2);

  // block size
  blkSize_ = pl_.get("Block Size", problem_->getNEV());

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
  if (pl_.isParameter("Verbosity")) {
    if (Teuchos::isParameterType<int>(pl_,"Verbosity")) {
      verbosity = pl_.get("Verbosity", verbosity);
    } else {
      verbosity = (int)Teuchos::getParameter<Anasazi::MsgType>(pl_,"Verbosity");
    }
  }
  printer_ = Teuchos::rcp( new OutputManager<ScalarType>(verbosity,osp) );

}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType
RTRSolMgr<ScalarType,MV,OP>::solve() {

  using std::endl;

  // typedef SolverUtils<ScalarType,MV,OP> msutils; // unused

  const int nev = problem_->getNEV();

  // clear num iters
  numIters_ = -1;

#ifdef TEUCHOS_DEBUG
    Teuchos::RCP<Teuchos::FancyOStream>
      out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(printer_->stream(Debug)));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
    *out << "Entering Anasazi::RTRSolMgr::solve()\n";
#endif

  //////////////////////////////////////////////////////////////////////////////////////
  // Sort manager
  Teuchos::RCP<BasicSort<MagnitudeType> > sorter = Teuchos::rcp( new BasicSort<MagnitudeType>(whch_) );

  //////////////////////////////////////////////////////////////////////////////////////
  // Status tests
  //
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > maxtest;
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > ordertest;
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > combotest;
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > convtest;
  // maximum iters
  if (maxIters_ > 0) {
    maxtest = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>(maxIters_) );
  }
  else {
    maxtest = Teuchos::null;
  }
  //
  // convergence
  convtest = Teuchos::rcp( new StatusTestResNorm<ScalarType,MV,OP>(convtol_,nev,convNorm_,relconvtol_) );
  ordertest = Teuchos::rcp( new StatusTestWithOrdering<ScalarType,MV,OP>(convtest,sorter,nev) );
  //
  // combo
  Teuchos::Array<Teuchos::RCP<StatusTest<ScalarType,MV,OP> > > alltests;
  alltests.push_back(ordertest);
  if (maxtest != Teuchos::null) alltests.push_back(maxtest);
  combotest = Teuchos::rcp( new StatusTestCombo<ScalarType,MV,OP>( StatusTestCombo<ScalarType,MV,OP>::OR, alltests) );
  //
  // printing StatusTest
  Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputtest;
  if ( printer_->isVerbosity(Debug) ) {
    outputtest = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer_,combotest,1,Passed+Failed+Undefined ) );
  }
  else {
    outputtest = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer_,combotest,1,Passed ) );
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Orthomanager
  Teuchos::RCP<ICGSOrthoManager<ScalarType,MV,OP> > ortho
    = Teuchos::rcp( new ICGSOrthoManager<ScalarType,MV,OP>(problem_->getM(),numICGS_) );

  //////////////////////////////////////////////////////////////////////////////////////
  // create an RTR solver
  // leftmost or rightmost?
  bool leftMost = true;
  if (whch_ == "LR" || whch_ == "LM") {
    leftMost = false;
  }
  pl_.set<bool>("Leftmost",leftMost);
  Teuchos::RCP<RTRBase<ScalarType,MV,OP> > rtr_solver;
  if (skinny_ == false) {
    // "hefty" IRTR
    rtr_solver = Teuchos::rcp( new  IRTR<ScalarType,MV,OP>(problem_,sorter,printer_,outputtest,ortho,pl_) );
  }
  else {
    // "skinny" IRTR
    rtr_solver = Teuchos::rcp( new SIRTR<ScalarType,MV,OP>(problem_,sorter,printer_,outputtest,ortho,pl_) );
  }
  // set any auxiliary vectors defined in the problem
  Teuchos::RCP< const MV > probauxvecs = problem_->getAuxVecs();
  if (probauxvecs != Teuchos::null) {
    rtr_solver->setAuxVecs( Teuchos::tuple< Teuchos::RCP<const MV> >(probauxvecs) );
  }

  TEUCHOS_TEST_FOR_EXCEPTION(rtr_solver->getBlockSize() < problem_->getNEV(),std::logic_error,
            "Anasazi::RTRSolMgr requires block size >= requested number of eigenvalues.");

  int numfound = 0;
  Teuchos::RCP<MV> foundvecs;
  std::vector<MagnitudeType> foundvals;

  // tell the solver to iterate
  try {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*_timerSolve);
#endif
    rtr_solver->iterate();
    numIters_ = rtr_solver->getNumIters();
  }
  catch (const std::exception &e) {
    // we are a simple solver manager. we don't catch exceptions. set solution empty, then rethrow.
    printer_->stream(Anasazi::Errors) << "Exception: " << e.what() << endl;
    Eigensolution<ScalarType,MV> sol;
    sol.numVecs = 0;
    problem_->setSolution(sol);
    throw;
  }

  // check the status tests
  if (convtest->getStatus() == Passed || (maxtest != Teuchos::null && maxtest->getStatus() == Passed))
  {
    int num = convtest->howMany();
    if (num > 0) {
      std::vector<int> ind = convtest->whichVecs();
      // copy the converged eigenvectors
      foundvecs = MVT::CloneCopy(*rtr_solver->getRitzVectors(),ind);
      // copy the converged eigenvalues
      foundvals.resize(num);
      std::vector<Value<ScalarType> > all = rtr_solver->getRitzValues();
      for (int i=0; i<num; i++) {
        foundvals[i] = all[ind[i]].realpart;
      }
      numfound = num;
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Anasazi::RTRSolMgr::solve(): solver returned without satisfy status test.");
  }

  // create contiguous storage for all eigenvectors, eigenvalues
  Eigensolution<ScalarType,MV> sol;
  sol.numVecs = numfound;
  sol.Evecs = foundvecs;
  sol.Espace = sol.Evecs;
  sol.Evals.resize(sol.numVecs);
  for (int i=0; i<sol.numVecs; i++) {
    sol.Evals[i].realpart = foundvals[i];
  }
  // all real eigenvalues: set index vectors [0,...,numfound-1]
  sol.index.resize(numfound,0);

  // print final summary
  rtr_solver->currentStatus(printer_->stream(FinalSummary));

  // print timing information
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  if ( printer_->isVerbosity( TimingDetails ) ) {
    Teuchos::TimeMonitor::summarize( printer_->stream( TimingDetails ) );
  }
#endif

  // send the solution to the eigenproblem
  problem_->setSolution(sol);
  printer_->stream(Debug) << "Returning " << sol.numVecs << " eigenpairs to eigenproblem." << endl;

  // return from SolMgr::solve()
  if (sol.numVecs < nev) return Unconverged;
  return Converged;
}




} // end Anasazi namespace

#endif /* ANASAZI_RTR_SOLMGR_HPP */
