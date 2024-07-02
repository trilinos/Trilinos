// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_BLOCK_KRYLOV_SCHUR_SOLMGR_HPP
#define ANASAZI_BLOCK_KRYLOV_SCHUR_SOLMGR_HPP

/// \file AnasaziBlockKrylovSchurSolMgr.hpp
/// \brief The Anasazi::BlockKrylovSchurSolMgr class provides a user
///   interface for the block Krylov-Schur eigensolver.

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziSolverUtils.hpp"

#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziICGSOrthoManager.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziStatusTestWithOrdering.hpp"
#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziStatusTestOutput.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziOutputStreamTraits.hpp"

#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_FancyOStream.hpp"

/** \example BlockKrylovSchur/BlockKrylovSchurEpetraEx.cpp
    \brief Use Anasazi::BlockKrylovSchurSolMgr to solve a standard
      (not generalized) eigenvalue problem, using Epetra data
      structures.
*/

/// \example BlockKrylovSchur/BlockKrylovSchurEpetraExGenAmesos.cpp
/// \brief Compute smallest eigenvalues of a generalized eigenvalue
///   problem, using block Krylov-Schur with Epetra and an Amesos direct
///   solver.
///
/// This example computes the eigenvalues of smallest magnitude of a
/// generalized eigenvalue problem \f$K x = \lambda M x\f$, using
/// Anasazi's implementation of the block Krylov-Schur method, with
/// Epetra linear algebra and a direct solver from the Amesos package.
///
/// Anasazi computes the smallest-magnitude eigenvalues using a
/// shift-and-invert strategy.  For simplicity, this example uses a
/// shift of zero.  It illustrates the general pattern for using
/// Anasazi for this problem:
///
///   1. Construct an "operator" A such that \f$Az = K^{-1} M z\f$.
///   2. Use Anasazi to solve \f$Az = \sigma z\f$, which is a spectral
///      transformation of the original problem \f$K x = \lambda M x\f$.
///   3. The eigenvalues \f$\lambda\f$ of the original problem are the
///      inverses of the eigenvalues \f$\sigma\f$ of the transformed
///      problem.
///
/// In the example, the "operator A such that \f$A z = K^{-1} M z\f$"
/// is a subclass of Epetra_Operator.  The Apply method of that
/// operator takes the vector b, and computes \f$x = K^{-1} M b\f$.
/// It does so first by applying the matrix M, and then by solving the
/// linear system \f$K x = M b\f$ for x.  Trilinos implements many
/// different ways to solve linear systems.  The example uses the
/// sparse direct solver KLU to do so.  Trilinos' Amesos package has
/// an interface to KLU.

/** \example BlockKrylovSchur/BlockKrylovSchurEpetraExGenAztecOO.cpp
    \brief Use Anasazi::BlockKrylovSchurSolMgr to solve a generalized
      eigenvalue problem, using Epetra data stuctures and the AztecOO
      package of iterative linear solvers and preconditioners.
*/

/** \example BlockKrylovSchur/BlockKrylovSchurEpetraExGenBelos.cpp
    \brief Use Anasazi::BlockKrylovSchurSolMgr to solve a generalized
      eigenvalue problem, using Epetra data stuctures and the Belos
      iterative linear solver package.
*/

/** \example BlockKrylovSchur/BlockKrylovSchurEpetraExSVD.cpp
    \brief Use Anasazi::BlockKrylovSchurSolMgr to compute a singular
      value decomposition (SVD), using Epetra data structures.
*/

namespace Anasazi {


/*! \class BlockKrylovSchurSolMgr
 *
 *  \brief The Anasazi::BlockKrylovSchurSolMgr provides a flexible solver manager over the BlockKrylovSchur eigensolver.
 *
 * The solver manager provides to the solver a StatusTestCombo object constructed as follows:<br>
 *    &nbsp;&nbsp;&nbsp;<tt>combo = globaltest OR debugtest</tt><br>
 * where
 *    - \c globaltest terminates computation when global convergence has been detected.<br>
 *      It is encapsulated in a StatusTestWithOrdering object, to ensure that computation is terminated
 *      only after the most significant eigenvalues/eigenvectors have met the convergence criteria.<br>
 *      If not specified via setGlobalStatusTest(), this test is a StatusTestResNorm object which tests the
 *      2-norms of the Ritz residuals relative to the Ritz values.
 *    - \c debugtest allows a user to specify additional monitoring of the iteration, encapsulated in a StatusTest object<br>
 *      If not specified via setDebugStatusTest(), \c debugtest is ignored.<br>
 *      In most cases, it should return ::Failed; if it returns ::Passed, solve() will throw an AnasaziError exception.
 *
 * Additionally, the solver manager will terminate solve() after a specified number of restarts.
 *
 * Much of this behavior is controlled via parameters and options passed to the
 * solver manager. For more information, see BlockKrylovSchurSolMgr().

 \ingroup anasazi_solver_framework

 \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
 */

template<class ScalarType, class MV, class OP>
class BlockKrylovSchurSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;

  public:

  //! @name Constructors/Destructor
  //@{

  /*! \brief Basic constructor for BlockKrylovSchurSolMgr.
   *
   * This constructor accepts the Eigenproblem to be solved in addition
   * to a parameter list of options for the solver manager. These options include the following:
   *   - Solver parameters
   *      - \c "Which" - a \c string specifying the desired eigenvalues: SM, LM, SR or LR. Default: "LM"
   *      - \c "Block Size" - an \c int specifying the block size to be used by the underlying block Krylov-Schur solver. Default: 1
   *      - \c "Num Blocks" - an \c int specifying the number of blocks allocated for the Krylov basis. Default: 3*nev
   *      - \c "Extra NEV Blocks" - an \c int specifying the number of extra blocks the solver should keep in addition to those
             required to compute the number of eigenvalues requested.  Default: 0
   *      - \c "Maximum Restarts" - an \c int specifying the maximum number of restarts the underlying solver is allowed to perform. Default: 20
   *      - \c "Orthogonalization" - a \c string specifying the desired orthogonalization:  DGKS, ICGS, and SVQB. Default: "SVQB"
   *      - \c "Verbosity" - a sum of MsgType specifying the verbosity. Default: Anasazi::Errors
   *      - \c "Output Stream" - a reference-counted pointer to the formatted output stream where all
   *                             solver output is sent.  Default: Teuchos::getFancyOStream ( Teuchos::rcpFromRef (std::cout) )
   *      - \c "Output Processor" - an \c int specifying the MPI processor that will print solver/timer details.  Default: 0 
   *   - Convergence parameters
   *      - \c "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: machine precision.
   *      - \c "Relative Convergence Tolerance" - a \c bool specifying whether residuals norms should be scaled by their eigenvalues for the purposing of deciding convergence. Default: true
   */
  BlockKrylovSchurSolMgr( const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
                             Teuchos::ParameterList &pl );

  //! Destructor.
  virtual ~BlockKrylovSchurSolMgr() {};
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

  /*! \brief Return the Ritz values from the most recent solve.
   */
  std::vector<Value<ScalarType> > getRitzValues() const {
    std::vector<Value<ScalarType> > ret( ritzValues_ );
    return ret;
  }

  /*! \brief Return the timers for this object.
   *
   * The timers are ordered as follows:
   *   - time spent in solve() routine
   *   - time spent restarting
   */
   Teuchos::Array<Teuchos::RCP<Teuchos::Time> > getTimers() const {
     return Teuchos::tuple(timerSolve_, timerRestarting_);
   }

  //@}

  //! @name Solver application methods
  //@{

  /*! \brief This method performs possibly repeated calls to the underlying eigensolver's iterate() routine
   * until the problem has been solved (as decided by the solver manager) or the solver manager decides to
   * quit.
   *
   * This method calls BlockKrylovSchur::iterate(), which will return either because a specially constructed status test evaluates to ::Passed
   * or an exception is thrown.
   *
   * A return from BlockKrylovSchur::iterate() signifies one of the following scenarios:
   *    - the maximum number of restarts has been exceeded. In this scenario, the solver manager will place\n
   *      all converged eigenpairs into the eigenproblem and return ::Unconverged.
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
  void setGlobalStatusTest(const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > &global);

  //! Get the status test defining global convergence.
  const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > & getGlobalStatusTest() const;

  //! Set the status test for debugging.
  void setDebugStatusTest(const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > &debug);

  //! Get the status test for debugging.
  const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > & getDebugStatusTest() const;

  //@}

  private:
  Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > problem_;
  Teuchos::RCP<SortManager<MagnitudeType> > sort_;

  std::string whch_, ortho_;
  MagnitudeType ortho_kappa_;

  MagnitudeType convtol_;
  int maxRestarts_;
  bool relconvtol_,conjSplit_;
  int blockSize_, numBlocks_, stepSize_, nevBlocks_, xtra_nevBlocks_;
  int numIters_;
  int verbosity_;
  bool inSituRestart_, dynXtraNev_;
  int osProc_;

  std::vector<Value<ScalarType> > ritzValues_;

  int printNum_;
  Teuchos::RCP<Teuchos::Time> timerSolve_, timerRestarting_;

  Teuchos::RCP<Teuchos::FancyOStream> osp_;

  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > globalTest_;
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > debugTest_;

};


// Constructor
template<class ScalarType, class MV, class OP>
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::BlockKrylovSchurSolMgr(
        const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
        Teuchos::ParameterList &pl ) :
  problem_(problem),
  whch_("LM"),
  ortho_("SVQB"),
  ortho_kappa_(-1.0),
  convtol_(0),
  maxRestarts_(20),
  relconvtol_(true),
  conjSplit_(false),
  blockSize_(0),
  numBlocks_(0),
  stepSize_(0),
  nevBlocks_(0),
  xtra_nevBlocks_(0),
  numIters_(0),
  verbosity_(Anasazi::Errors),
  inSituRestart_(false),
  dynXtraNev_(false),
  osProc_(0),
  printNum_(-1)
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  ,timerSolve_(Teuchos::TimeMonitor::getNewTimer("Anasazi: BlockKrylovSchurSolMgr::solve()")),
  timerRestarting_(Teuchos::TimeMonitor::getNewTimer("Anasazi: BlockKrylovSchurSolMgr restarting"))
#endif
{
  TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null,               std::invalid_argument, "Problem not given to solver manager.");
  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isProblemSet(),               std::invalid_argument, "Problem not set.");
  TEUCHOS_TEST_FOR_EXCEPTION(problem_->getInitVec() == Teuchos::null, std::invalid_argument, "Problem does not contain initial vectors to clone from.");

  const int nev = problem_->getNEV();

  // convergence tolerance
  convtol_ = pl.get("Convergence Tolerance",MT::prec());
  relconvtol_ = pl.get("Relative Convergence Tolerance",relconvtol_);

  // maximum number of restarts
  maxRestarts_ = pl.get("Maximum Restarts",maxRestarts_);

  // block size: default is 1
  blockSize_ = pl.get("Block Size",1);
  TEUCHOS_TEST_FOR_EXCEPTION(blockSize_ <= 0, std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: \"Block Size\" must be strictly positive.");

  // set the number of blocks we need to save to compute the nev eigenvalues of interest.
  xtra_nevBlocks_ = pl.get("Extra NEV Blocks",0);
  if (nev%blockSize_) {
    nevBlocks_ = nev/blockSize_ + 1;
  } else {
    nevBlocks_ = nev/blockSize_;
  }

  // determine if we should use the dynamic scheme for selecting the current number of retained eigenvalues.
  // NOTE:  This employs a technique similar to ARPACK in that it increases the number of retained eigenvalues
  //        by one for every converged eigenpair to accelerate convergence.
  if (pl.isParameter("Dynamic Extra NEV")) {
    if (Teuchos::isParameterType<bool>(pl,"Dynamic Extra NEV")) {
      dynXtraNev_ = pl.get("Dynamic Extra NEV",dynXtraNev_);
    } else {
      dynXtraNev_ = ( Teuchos::getParameter<int>(pl,"Dynamic Extra NEV") != 0 );
    }
  }

  numBlocks_ = pl.get("Num Blocks",3*nevBlocks_);
  TEUCHOS_TEST_FOR_EXCEPTION(numBlocks_ <= nevBlocks_, std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: \"Num Blocks\" must be strictly positive and large enough to compute the requested eigenvalues.");

  TEUCHOS_TEST_FOR_EXCEPTION(static_cast<ptrdiff_t>(numBlocks_)*blockSize_ > MVT::GetGlobalLength(*problem_->getInitVec()),
                     std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: Potentially impossible orthogonality requests. Reduce basis size.");

  // step size: the default is maxRestarts_*numBlocks_, so that Ritz values are only computed every restart.
  if (maxRestarts_) {
    stepSize_ = pl.get("Step Size", (maxRestarts_+1)*(numBlocks_+1));
  } else {
    stepSize_ = pl.get("Step Size", numBlocks_+1);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(stepSize_ < 1, std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: \"Step Size\" must be strictly positive.");

  // get the sort manager
  if (pl.isParameter("Sort Manager")) {
    sort_ = Teuchos::getParameter<Teuchos::RCP<Anasazi::SortManager<MagnitudeType> > >(pl,"Sort Manager");
  } else {
    // which values to solve for
    whch_ = pl.get("Which",whch_);
    TEUCHOS_TEST_FOR_EXCEPTION(whch_ != "SM" && whch_ != "LM" && whch_ != "SR" && whch_ != "LR" && whch_ != "SI" && whch_ != "LI",
                       std::invalid_argument, "Invalid sorting string.");
    sort_ = Teuchos::rcp( new BasicSort<MagnitudeType>(whch_) );
  }

  // which orthogonalization to use
  ortho_ = pl.get("Orthogonalization",ortho_);
  if (ortho_ != "DGKS" && ortho_ != "SVQB" && ortho_ != "ICGS") {
    ortho_ = "SVQB";
  }

  // which orthogonalization constant to use
  ortho_kappa_ = pl.get("Orthogonalization Constant",ortho_kappa_);

  // Create a formatted output stream to print to.
  // See if user requests output processor.
  osProc_ = pl.get("Output Processor", osProc_);
  
  // If not passed in by user, it will be chosen based upon operator type.
  if (pl.isParameter("Output Stream")) {
    osp_ = Teuchos::getParameter<Teuchos::RCP<Teuchos::FancyOStream> >(pl,"Output Stream");
  }
  else {
    osp_ = OutputStreamTraits<OP>::getOutputStream (*problem_->getOperator(),osProc_);
  }

  // verbosity level
  if (pl.isParameter("Verbosity")) {
    if (Teuchos::isParameterType<int>(pl,"Verbosity")) {
      verbosity_ = pl.get("Verbosity", verbosity_);
    } else {
      verbosity_ = (int)Teuchos::getParameter<Anasazi::MsgType>(pl,"Verbosity");
    }
  }

  // restarting technique: V*Q or applyHouse(V,H,tau)
  if (pl.isParameter("In Situ Restarting")) {
    if (Teuchos::isParameterType<bool>(pl,"In Situ Restarting")) {
      inSituRestart_ = pl.get("In Situ Restarting",inSituRestart_);
    } else {
      inSituRestart_ = ( Teuchos::getParameter<int>(pl,"In Situ Restarting") != 0 );
    }
  }

  printNum_ = pl.get<int>("Print Number of Ritz Values",printNum_);
}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::solve() {

  const int nev = problem_->getNEV();
  ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

  Teuchos::BLAS<int,ScalarType> blas;
  Teuchos::LAPACK<int,ScalarType> lapack;

  //////////////////////////////////////////////////////////////////////////////////////
  // Output manager
  Teuchos::RCP<OutputManager<ScalarType> > printer = Teuchos::rcp( new OutputManager<ScalarType>(verbosity_,osp_) );

  //////////////////////////////////////////////////////////////////////////////////////
  // Status tests
  //
  // convergence
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > convtest;
  if (globalTest_ == Teuchos::null) {
    convtest = Teuchos::rcp( new StatusTestResNorm<ScalarType,MV,OP>(convtol_,nev,RITZRES_2NORM,relconvtol_) );
  }
  else {
    convtest = globalTest_;
  }
  Teuchos::RCP<StatusTestWithOrdering<ScalarType,MV,OP> > ordertest
    = Teuchos::rcp( new StatusTestWithOrdering<ScalarType,MV,OP>(convtest,sort_,nev) );
  // for a non-short-circuited OR test, the order doesn't matter
  Teuchos::Array<Teuchos::RCP<StatusTest<ScalarType,MV,OP> > > alltests;
  alltests.push_back(ordertest);

  if (debugTest_ != Teuchos::null) alltests.push_back(debugTest_);

  Teuchos::RCP<StatusTestCombo<ScalarType,MV,OP> > combotest
    = Teuchos::rcp( new StatusTestCombo<ScalarType,MV,OP>( StatusTestCombo<ScalarType,MV,OP>::OR, alltests) );
  // printing StatusTest
  Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputtest;
  if ( printer->isVerbosity(Debug) ) {
    outputtest = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer,combotest,1,Passed+Failed+Undefined ) );
  }
  else {
    outputtest = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer,combotest,1,Passed ) );
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Orthomanager
  Teuchos::RCP<OrthoManager<ScalarType,MV> > ortho;
  if (ortho_=="SVQB") {
    ortho = Teuchos::rcp( new SVQBOrthoManager<ScalarType,MV,OP>(problem_->getM()) );
  } else if (ortho_=="DGKS") {
    if (ortho_kappa_ <= 0) {
      ortho = Teuchos::rcp( new BasicOrthoManager<ScalarType,MV,OP>(problem_->getM()) );
    } else {
      ortho = Teuchos::rcp( new BasicOrthoManager<ScalarType,MV,OP>(problem_->getM(),ortho_kappa_) );
    }
  } else if (ortho_=="ICGS") {
    ortho = Teuchos::rcp( new ICGSOrthoManager<ScalarType,MV,OP>(problem_->getM()) );
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(ortho_!="SVQB"&&ortho_!="DGKS"&&ortho_!="ICGS",std::logic_error,"Anasazi::BlockKrylovSchurSolMgr::solve(): Invalid orthogonalization type.");
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Parameter list
  Teuchos::ParameterList plist;
  plist.set("Block Size",blockSize_);
  plist.set("Num Blocks",numBlocks_);
  plist.set("Step Size",stepSize_);
  if (printNum_ == -1) {
    plist.set("Print Number of Ritz Values",nevBlocks_*blockSize_);
  }
  else {
    plist.set("Print Number of Ritz Values",printNum_);
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // BlockKrylovSchur solver
  Teuchos::RCP<BlockKrylovSchur<ScalarType,MV,OP> > bks_solver
    = Teuchos::rcp( new BlockKrylovSchur<ScalarType,MV,OP>(problem_,sort_,printer,outputtest,ortho,plist) );
  // set any auxiliary vectors defined in the problem
  Teuchos::RCP< const MV > probauxvecs = problem_->getAuxVecs();
  if (probauxvecs != Teuchos::null) {
    bks_solver->setAuxVecs( Teuchos::tuple< Teuchos::RCP<const MV> >(probauxvecs) );
  }

  // Create workspace for the Krylov basis generated during a restart
  // Need at most (nevBlocks_*blockSize_+1) for the updated factorization and another block for the current factorization residual block (F).
  //  ---> (nevBlocks_*blockSize_+1) + blockSize_
  // If Hermitian, this becomes nevBlocks_*blockSize_ + blockSize_
  // we only need this if there is the possibility of restarting, ex situ

  // Maximum allowable extra vectors that BKS can keep to accelerate convergence
  int maxXtraBlocks = 0;
  if ( dynXtraNev_ ) maxXtraBlocks = ( ( bks_solver->getMaxSubspaceDim() - nev ) / blockSize_ ) / 2;

  Teuchos::RCP<MV> workMV;
  if (maxRestarts_ > 0) {
    if (inSituRestart_==true) {
      // still need one work vector for applyHouse()
      workMV = MVT::Clone( *problem_->getInitVec(), 1 );
    }
    else { // inSituRestart == false
      if (problem_->isHermitian()) {
        workMV = MVT::Clone( *problem_->getInitVec(), (nevBlocks_+maxXtraBlocks)*blockSize_ + blockSize_ );
      } else {
        // Increase workspace by 1 because of potential complex conjugate pairs on the Ritz vector boundary
        workMV = MVT::Clone( *problem_->getInitVec(), (nevBlocks_+maxXtraBlocks)*blockSize_ + 1 + blockSize_ );
      }
    }
  } else {
    workMV = Teuchos::null;
  }

  // go ahead and initialize the solution to nothing in case we throw an exception
  Eigensolution<ScalarType,MV> sol;
  sol.numVecs = 0;
  problem_->setSolution(sol);

  int numRestarts = 0;
  int cur_nevBlocks = 0;

  // enter solve() iterations
  {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*timerSolve_);
#endif

    // tell bks_solver to iterate
    while (1) {
      try {
        bks_solver->iterate();

        ////////////////////////////////////////////////////////////////////////////////////
        //
        // check convergence first
        //
        ////////////////////////////////////////////////////////////////////////////////////
        if ( ordertest->getStatus() == Passed ) {
          // we have convergence
          // ordertest->whichVecs() tells us which vectors from solver state are the ones we want
          // ordertest->howMany() will tell us how many
          break;
        }
        ////////////////////////////////////////////////////////////////////////////////////
        //
        // check for restarting, i.e. the subspace is full
        //
        ////////////////////////////////////////////////////////////////////////////////////
        // this is for the Hermitian case, or non-Hermitian conjugate split situation.
        // --> for the Hermitian case the current subspace dimension needs to match the maximum subspace dimension
        // --> for the non-Hermitian case:
        //     --> if a conjugate pair was detected in the previous restart then the current subspace dimension needs to match the
        //         maximum subspace dimension (the BKS solver keeps one extra vector if the problem is non-Hermitian).
        //     --> if a conjugate pair was not detected in the previous restart then the current subspace dimension will be one less
        //         than the maximum subspace dimension.
        else if ( (bks_solver->getCurSubspaceDim() == bks_solver->getMaxSubspaceDim()) ||
                  (!problem_->isHermitian() && !conjSplit_ && (bks_solver->getCurSubspaceDim()+1 == bks_solver->getMaxSubspaceDim())) ) {

          // Update the Schur form of the projected eigenproblem, then sort it.
          if (!bks_solver->isSchurCurrent()) {
            bks_solver->computeSchurForm( true );

            // Check for convergence, just in case we wait for every restart to check
            outputtest->checkStatus( &*bks_solver );
          }

          // Don't bother to restart if we've converged or reached the maximum number of restarts
          if ( numRestarts >= maxRestarts_ || ordertest->getStatus() == Passed) {
            break; // break from while(1){bks_solver->iterate()}
          }

          // Start restarting timer and increment counter
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor restimer(*timerRestarting_);
#endif
          numRestarts++;

          int numConv = ordertest->howMany();
          cur_nevBlocks = nevBlocks_*blockSize_;

          // Add in extra blocks for restarting if either static or dynamic boundaries are being used.
          int moreNevBlocks = std::min( maxXtraBlocks, std::max( numConv/blockSize_, xtra_nevBlocks_) );
          if ( dynXtraNev_ )
            cur_nevBlocks += moreNevBlocks * blockSize_;
          else if ( xtra_nevBlocks_ )
            cur_nevBlocks += xtra_nevBlocks_ * blockSize_;
/*
            int cur_numConv = numConv;
            while ( (cur_nevBlocks < (_nevBlocks + maxXtraVecs)) && cur_numConv > 0 ) {
              cur_nevBlocks++;
              cur_numConv--;
*/

          printer->stream(Debug) << " Performing restart number " << numRestarts << " of " << maxRestarts_ << std::endl << std::endl;
          printer->stream(Debug) << "   - Current NEV blocks is " << cur_nevBlocks << ", the minimum is " << nevBlocks_*blockSize_ << std::endl;

          // Get the most current Ritz values before we continue.
          ritzValues_ = bks_solver->getRitzValues();

          // Get the state.
          BlockKrylovSchurState<ScalarType,MV> oldState = bks_solver->getState();

          // Get the current dimension of the factorization
          int curDim = oldState.curDim;

          // Determine if the storage for the nev eigenvalues of interest splits a complex conjugate pair.
          std::vector<int> ritzIndex = bks_solver->getRitzIndex();
          if (ritzIndex[cur_nevBlocks-1]==1) {
            conjSplit_ = true;
            cur_nevBlocks++;
          } else {
            conjSplit_ = false;
          }

          // Print out a warning to the user if complex eigenvalues were found on the boundary of the restart subspace
          // and the eigenproblem is Hermitian.  This solver is not prepared to handle this situation.
          if (problem_->isHermitian() && conjSplit_)
          {
            printer->stream(Warnings)
              << " Eigenproblem is Hermitian, complex eigenvalues have been detected, and eigenvalues of interest split a conjugate pair!!!"
              << std::endl
              << " Block Krylov-Schur eigensolver cannot guarantee correct behavior in this situation, please turn Hermitian flag off!!!"
              << std::endl;
          }

          // Update the Krylov-Schur decomposition

          // Get a view of the Schur vectors of interest.
          Teuchos::SerialDenseMatrix<int,ScalarType> Qnev(Teuchos::View, *(oldState.Q), curDim, cur_nevBlocks);

          // Get a view of the current Krylov basis.
          std::vector<int> curind( curDim );
          for (int i=0; i<curDim; i++) { curind[i] = i; }
          Teuchos::RCP<const MV> basistemp = MVT::CloneView( *(oldState.V), curind );

          // Compute the new Krylov basis: Vnew = V*Qnev
          //
          // this will occur ex situ in workspace allocated for this purpose (tmpMV)
          // or in situ in the solver's memory space.
          //
          // we will also set a pointer for the location that the current factorization residual block (F),
          // currently located after the current basis in oldstate.V, will be moved to
          //
          Teuchos::RCP<MV> newF;
          if (inSituRestart_) {
            //
            // get non-const pointer to solver's basis so we can work in situ
            Teuchos::RCP<MV> solverbasis = Teuchos::rcp_const_cast<MV>(oldState.V);
            Teuchos::SerialDenseMatrix<int,ScalarType> copyQnev(Teuchos::Copy, Qnev);
            //
            // perform Householder QR of copyQnev = Q [D;0], where D is unit diag. We will want D below.
            std::vector<ScalarType> tau(cur_nevBlocks), work(cur_nevBlocks);
            int info;
            lapack.GEQRF(curDim,cur_nevBlocks,copyQnev.values(),copyQnev.stride(),&tau[0],&work[0],work.size(),&info);
            TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,
                               "Anasazi::BlockKrylovSchurSolMgr::solve(): error calling GEQRF during restarting.");
            // we need to get the diagonal of D
            std::vector<ScalarType> d(cur_nevBlocks);
            for (int j=0; j<copyQnev.numCols(); j++) {
              d[j] = copyQnev(j,j);
            }
            if (printer->isVerbosity(Debug)) {
              Teuchos::SerialDenseMatrix<int,ScalarType> R(Teuchos::Copy,copyQnev,cur_nevBlocks,cur_nevBlocks);
              for (int j=0; j<R.numCols(); j++) {
                R(j,j) = SCT::magnitude(R(j,j)) - 1.0;
                for (int i=j+1; i<R.numRows(); i++) {
                  R(i,j) = zero;
                }
              }
              printer->stream(Debug) << "||Triangular factor of Su - I||: " << R.normFrobenius() << std::endl;
            }
            //
            // perform implicit V*Qnev
            // this actually performs V*[Qnev Qtrunc*M] = [newV truncV], for some unitary M
            // we are interested in only the first cur_nevBlocks vectors of the result
            curind.resize(curDim);
            for (int i=0; i<curDim; i++) curind[i] = i;
            {
              Teuchos::RCP<MV> oldV = MVT::CloneViewNonConst(*solverbasis,curind);
              SolverUtils<ScalarType,MV,OP>::applyHouse(cur_nevBlocks,*oldV,copyQnev,tau,workMV);
            }
            // multiply newV*D
            // get pointer to new basis
            curind.resize(cur_nevBlocks);
            for (int i=0; i<cur_nevBlocks; i++) { curind[i] = i; }
            {
              Teuchos::RCP<MV> newV = MVT::CloneViewNonConst( *solverbasis, curind );
              MVT::MvScale(*newV,d);
            }
            // get pointer to new location for F
            curind.resize(blockSize_);
            for (int i=0; i<blockSize_; i++) { curind[i] = cur_nevBlocks + i; }
            newF = MVT::CloneViewNonConst( *solverbasis, curind );
          }
          else {
            // get pointer to first part of work space
            curind.resize(cur_nevBlocks);
            for (int i=0; i<cur_nevBlocks; i++) { curind[i] = i; }
            Teuchos::RCP<MV> tmp_newV = MVT::CloneViewNonConst(*workMV, curind );
            // perform V*Qnev
            MVT::MvTimesMatAddMv( one, *basistemp, Qnev, zero, *tmp_newV );
            tmp_newV = Teuchos::null;
            // get pointer to new location for F
            curind.resize(blockSize_);
            for (int i=0; i<blockSize_; i++) { curind[i] = cur_nevBlocks + i; }
            newF = MVT::CloneViewNonConst( *workMV, curind );
          }

          // Move the current factorization residual block (F) to the last block of newV.
          curind.resize(blockSize_);
          for (int i=0; i<blockSize_; i++) { curind[i] = curDim + i; }
          Teuchos::RCP<const MV> oldF = MVT::CloneView( *(oldState.V), curind );
          for (int i=0; i<blockSize_; i++) { curind[i] = i; }
          MVT::SetBlock( *oldF, curind, *newF );
          newF = Teuchos::null;

          // Update the Krylov-Schur quasi-triangular matrix.
          //
          // Create storage for the new Schur matrix of the Krylov-Schur factorization
          // Copy over the current quasi-triangular factorization of oldState.H which is stored in oldState.S.
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > newH =
            Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::Copy, *(oldState.S), cur_nevBlocks+blockSize_, cur_nevBlocks) );
          //
          // Get a view of the B block of the current factorization
          Teuchos::SerialDenseMatrix<int,ScalarType> oldB(Teuchos::View, *(oldState.H), blockSize_, blockSize_, curDim, curDim-blockSize_);
          //
          // Get a view of the a block row of the Schur vectors.
          Teuchos::SerialDenseMatrix<int,ScalarType> subQ(Teuchos::View, *(oldState.Q), blockSize_, cur_nevBlocks, curDim-blockSize_);
          //
          // Get a view of the new B block of the updated Krylov-Schur factorization
          Teuchos::SerialDenseMatrix<int,ScalarType> newB(Teuchos::View, *newH,  blockSize_, cur_nevBlocks, cur_nevBlocks);
          //
          // Compute the new B block.
          blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, blockSize_, cur_nevBlocks, blockSize_, one,
                     oldB.values(), oldB.stride(), subQ.values(), subQ.stride(), zero, newB.values(), newB.stride() );


          //
          // Set the new state and initialize the solver.
          BlockKrylovSchurState<ScalarType,MV> newstate;
          if (inSituRestart_) {
            newstate.V = oldState.V;
          } else {
            newstate.V = workMV;
          }
          newstate.H = newH;
          newstate.curDim = cur_nevBlocks;
          bks_solver->initialize(newstate);

        } // end of restarting
        ////////////////////////////////////////////////////////////////////////////////////
        //
        // we returned from iterate(), but none of our status tests Passed.
        // something is wrong, and it is probably our fault.
        //
        ////////////////////////////////////////////////////////////////////////////////////
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Anasazi::BlockKrylovSchurSolMgr::solve(): Invalid return from bks_solver::iterate().");
        }
      }
      catch (const AnasaziError &err) {
        printer->stream(Errors)
          << "Anasazi::BlockKrylovSchurSolMgr::solve() caught unexpected exception from Anasazi::BlockKrylovSchur::iterate() at iteration " << bks_solver->getNumIters() << std::endl
          << err.what() << std::endl
          << "Anasazi::BlockKrylovSchurSolMgr::solve() returning Unconverged with no solutions." << std::endl;
        return Unconverged;
      }
    }

    //
    // free temporary space
    workMV = Teuchos::null;

    // Get the most current Ritz values before we return
    ritzValues_ = bks_solver->getRitzValues();

    sol.numVecs = ordertest->howMany();
    printer->stream(Debug) << "ordertest->howMany() : " << sol.numVecs << std::endl;
    std::vector<int> whichVecs = ordertest->whichVecs();

    // Place any converged eigenpairs in the solution container.
    if (sol.numVecs > 0) {

      // Next determine if there is a conjugate pair on the boundary and resize.
      std::vector<int> tmpIndex = bks_solver->getRitzIndex();
      for (int i=0; i<(int)ritzValues_.size(); ++i) {
        printer->stream(Debug) << ritzValues_[i].realpart << " + i " << ritzValues_[i].imagpart << ", Index = " << tmpIndex[i] << std::endl;
      }
      printer->stream(Debug) << "Number of converged eigenpairs (before) = " << sol.numVecs << std::endl;
      for (int i=0; i<sol.numVecs; ++i) {
        printer->stream(Debug) << "whichVecs[" << i << "] = " << whichVecs[i] << ", tmpIndex[" << whichVecs[i] << "] = " << tmpIndex[whichVecs[i]] << std::endl;
      }
      if (tmpIndex[whichVecs[sol.numVecs-1]]==1) {
        printer->stream(Debug) << "There is a conjugate pair on the boundary, resizing sol.numVecs" << std::endl;
        whichVecs.push_back(whichVecs[sol.numVecs-1]+1);
        sol.numVecs++;
        for (int i=0; i<sol.numVecs; ++i) {
          printer->stream(Debug) << "whichVecs[" << i << "] = " << whichVecs[i] << ", tmpIndex[" << whichVecs[i] << "] = " << tmpIndex[whichVecs[i]] << std::endl;
        }
      }

      bool keepMore = false;
      int numEvecs = sol.numVecs;
      printer->stream(Debug) << "Number of converged eigenpairs (after) = " << sol.numVecs << std::endl;
      printer->stream(Debug) << "whichVecs[sol.numVecs-1] > sol.numVecs-1 : " << whichVecs[sol.numVecs-1] << " > " << sol.numVecs-1 << std::endl;
      if (whichVecs[sol.numVecs-1] > (sol.numVecs-1)) {
        keepMore = true;
        numEvecs = whichVecs[sol.numVecs-1]+1;  // Add 1 to fix zero-based indexing
        printer->stream(Debug) << "keepMore = true; numEvecs = " << numEvecs << std::endl;
      }

      // Next set the number of Ritz vectors that the iteration must compute and compute them.
      bks_solver->setNumRitzVectors(numEvecs);
      bks_solver->computeRitzVectors();

      // If the leading Ritz pairs are the converged ones, get the information
      // from the iteration to the solution container. Otherwise copy the necessary
      // information using 'whichVecs'.
      if (!keepMore) {
        sol.index = bks_solver->getRitzIndex();
        sol.Evals = bks_solver->getRitzValues();
        sol.Evecs = MVT::CloneCopy( *(bks_solver->getRitzVectors()) );
      }

      // Resize based on the number of solutions being returned and set the number of Ritz
      // vectors for the iteration to compute.
      sol.Evals.resize(sol.numVecs);
      sol.index.resize(sol.numVecs);

      // If the converged Ritz pairs are not the leading ones, copy over the information directly.
      if (keepMore) {
        std::vector<Anasazi::Value<ScalarType> > tmpEvals = bks_solver->getRitzValues();
        for (int vec_i=0; vec_i<sol.numVecs; ++vec_i) {
          sol.index[vec_i] = tmpIndex[whichVecs[vec_i]];
          sol.Evals[vec_i] = tmpEvals[whichVecs[vec_i]];
        }
        sol.Evecs = MVT::CloneCopy( *(bks_solver->getRitzVectors()), whichVecs );
      }

      // Set the solution space to be the Ritz vectors at this time.
      sol.Espace = sol.Evecs;
    }
  }

  // print final summary
  bks_solver->currentStatus(printer->stream(FinalSummary));

  // print timing information
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  if ( printer->isVerbosity( TimingDetails ) ) {
    Teuchos::TimeMonitor::summarize( printer->stream( TimingDetails ) );
  }
#endif

  problem_->setSolution(sol);
  printer->stream(Debug) << "Returning " << sol.numVecs << " eigenpairs to eigenproblem." << std::endl;

  // get the number of iterations performed during this solve.
  numIters_ = bks_solver->getNumIters();

  if (sol.numVecs < nev) {
    return Unconverged; // return from BlockKrylovSchurSolMgr::solve()
  }
  return Converged; // return from BlockKrylovSchurSolMgr::solve()
}


template <class ScalarType, class MV, class OP>
void
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::setGlobalStatusTest(
    const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > &global)
{
  globalTest_ = global;
}

template <class ScalarType, class MV, class OP>
const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > &
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::getGlobalStatusTest() const
{
  return globalTest_;
}

template <class ScalarType, class MV, class OP>
void
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::setDebugStatusTest(
    const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > &debug)
{
  debugTest_ = debug;
}

template <class ScalarType, class MV, class OP>
const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > &
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::getDebugStatusTest() const
{
  return debugTest_;
}

} // end Anasazi namespace

#endif /* ANASAZI_BLOCK_KRYLOV_SCHUR_SOLMGR_HPP */
