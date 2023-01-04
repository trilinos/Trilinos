
// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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

//TODO Check these includes!!
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

  The algorithms is well-known. See \\TODO citations. 

  \ingroup anasazi_solver_framework

  \author Jennifer A. Loe, Erik G. Boman
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
   *   - "Which" - a \c string specifying the desired eigenvalues: SM, LM, SR or LR. Default: SR
   *  NOTE: LM, SM, LR options are not yet implemented. 
   *   - "Block Size" - an \c int specifying the block size to be used by the underlying LOBPCG solver. Default: problem->getNEV()
   *   - "Maximum Iterations" - an \c int specifying the maximum number of iterations the underlying solver is allowed to perform. For the randomized solver, this corresponds to the number of times we multiply A by the initial vectors. Default: 5
   *   - "Orthogonalization" - a \c string specifying the desired orthogonalization: DGKS, ICGS, and SVQB. Default: "SVQB"
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
  std::string ortho_;
  int blockSize_;
  int maxIters_;
  int numIters_;
  bool trackResNorms_;
  Teuchos::RCP<Teuchos::Time> timerOp_, timerOrtho_, timerSolve_;
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
    timerOp_(Teuchos::TimeMonitor::getNewTimer("Anasazi: Randomized::Operation Op*x")),
    timerOrtho_(Teuchos::TimeMonitor::getNewTimer("Anasazi: Randomized::Orthogonalization")),
    timerSolve_(Teuchos::TimeMonitor::getNewTimer("Anasazi: Randomized::solve()")),
#endif
  ortho_("SVQB"),
  blockSize_(0),
  maxIters_(5),
  numIters_(0),
  trackResNorms_(true)
{
  TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null,              std::invalid_argument, "Problem not given to solver manager.");
  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isProblemSet(),              std::invalid_argument, "Problem not set.");
  TEUCHOS_TEST_FOR_EXCEPTION(problem_->getInitVec() == Teuchos::null,std::invalid_argument, "Problem does not contain initial vectors to clone from.");

  whch_ = pl.get("Which","LM");
  TEUCHOS_TEST_FOR_EXCEPTION(whch_ != "LM", //&& whch_ != "SM" && whch_ != "SR" && whch_ != "LR",
                     AnasaziError,
                     "RandomizedSolMgr: \"Which\" parameter must be LM. Other options not currently implemented."); //TODO Other options??

  tol_ = pl.get("Convergence Tolerance",tol_);
  TEUCHOS_TEST_FOR_EXCEPTION(tol_ <= 0,
                     AnasaziError,
                     "RandomizedSolMgr: \"Tolerance\" parameter must be strictly postiive.");

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
}



////////////////////////////////////////////////////////////////////////////////////////
template<class ScalarType, class MV, class OP>
ReturnType
RandomizedSolMgr<ScalarType,MV,OP>::solve() {
//std::cout << "DEBUG: In the solve function." << std::endl;

  // sort manager
  Teuchos::RCP<BasicSort<MT> > sorter = Teuchos::rcp( new BasicSort<MT> );
  // output manager
  Teuchos::RCP<OutputManager<ScalarType> > printer = Teuchos::rcp( new OutputManager<ScalarType>(verb_,osp_) );
  { //Timer 
  if(blockSize_ < problem_->getNEV()){ //TODO: Fix couts to proper ostream.
    std::cout << "Block size smaller than number evals. Increasing Block Size to num evals." << std::endl;
    blockSize_ = problem_->getNEV();
  }
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*timerSolve_);
#endif
  Teuchos::RCP<MV> randVecs;
    // grab some Multivector to Clone
    // in practice, getInitVec() should always provide this, but it is possible to use a 
    // Eigenproblem with nothing in getInitVec() by manually initializing with initialize(); 
    // in case of that strange scenario, we will try to Clone from V_; first resort to getInitVec(), 
    // because we would like to clear the storage associated with V_ so we have room for the new V_
  TEUCHOS_TEST_FOR_EXCEPTION(problem_->getInitVec() == Teuchos::null,std::invalid_argument,
    "Anasazi::Randomized: eigenproblem did not specify initial vectors to clone from.");
  //std::cout << "DEBUG: past teuchos check for getInitVec." << std::endl;
  if(MVT::GetNumberVecs(*(problem_->getInitVec()))==blockSize_){
    randVecs = MVT::CloneCopy(*(problem_->getInitVec()));
  }
  else{
    randVecs = MVT::Clone(*(problem_->getInitVec()),blockSize_);
    MVT::MvRandom(*randVecs);
  }
  //TEUCHOS_TEST_FOR_EXCEPTION(problem_->getA() == Teuchos::null,std::invalid_argument,
    //"Anasazi::Randomized: There is no A to get.");
  //std::cout << "DEBUG: got past check of getA." << std::endl;

  // Perform multiplies by A and Rayleigh-Ritz
  for( int i = 0; i < maxIters_; i++ ){
    {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor lcltimer( *timerOp_ );
#endif
    OPT::Apply( *(problem_->getOperator()), *randVecs, *randVecs );
    }
  }
  //std::cout  << "DEBUG: Past A multiply." << std::endl;
  // Set up Orthomanager and orthonormalize random vecs. 
  Teuchos::RCP<Anasazi::OrthoManager<ScalarType,MV> > orthoMgr;
  if (ortho_=="SVQB") {
    orthoMgr = Teuchos::rcp( new Anasazi::SVQBOrthoManager<ScalarType,MV,OP>());
  } else if (ortho_=="DGKS") {
    //if (ortho_kappa_ <= 0) {
    //  orthoMgr= Teuchos::rcp( new Anasazi::BasicOrthoManager<ScalarType,MV,OP>();
    //} else {
      orthoMgr = Teuchos::rcp( new Anasazi::BasicOrthoManager<ScalarType,MV,OP>());
    //}   
  } else if (ortho_=="ICGS") {
    orthoMgr = Teuchos::rcp( new Anasazi::ICGSOrthoManager<ScalarType,MV,OP>());
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(ortho_!="SVQB"&&ortho_!="DGKS"&&ortho_!="ICGS",std::logic_error,"Anasazi::RandomSolver Invalid orthogonalization type.");
  }
  //Orthogonalize the vectors
  int rank;
  {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor lcltimer( *timerOrtho_ );
#endif
  rank = orthoMgr->normalize(*randVecs);
  }
  if( rank < blockSize_ ){
    std::cout << "Warning! Anasazi::RandomSolver Random vectors did not have full rank!" << std::endl;
  }

  //std::cout  << "DEBUG: Past orthog." << std::endl;
  //Compute H = Q^TAQ. (RR projection) 
  Teuchos::RCP<MV> TmpVecs = MVT::Clone(*randVecs,blockSize_);
  Teuchos::SerialDenseMatrix<OT,ScalarType> H (blockSize_, blockSize_);

  OPT::Apply( *(problem_->getOperator()), *randVecs, *TmpVecs ); 
  MVT::MvTransMv(ONE, *randVecs, *TmpVecs, H);

  // Solve projected eigenvalue problem.
  Teuchos::LAPACK<OT,ScalarType> lapack;
  Teuchos::SerialDenseMatrix<OT,ScalarType> evects (blockSize_, blockSize_);
  std::vector<MT> evals_real(blockSize_);
  std::vector<MT> evals_imag(blockSize_);

  // Size of workspace and workspace for DGEEV
  int info = -1000; 
  ScalarType* vlr = 0; 
  const int ldv = 1; 

  int lwork = -1;
  std::vector<ScalarType> work(1);
  std::vector<MT> rwork(2*blockSize_);

  //Find workspace size for DGEEV:

  // 'N' no left eigenvectors. 
  // 'V' to compute right eigenvectors. 
  // blockSize = dimension of H
  // H matrix
  // H.stride = leading dimension of H
  // Array to store evals, real parts
  // Array to store evals, imag parts
  // vlr -> stores left evects, so don't need this
  // lead dim of vlr
  // evects =  array to store right evects
  // evects.stride = lead dim ovf evects
  // work = work array
  // lwork -1 means to query for array size
  // rwork - not referenced because ScalarType is not complex
  lapack.GEEV('N','V',blockSize_,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);
  lwork = std::abs (static_cast<int> (Teuchos::ScalarTraits<ScalarType>::real (work[0])));
  work.resize( lwork );
  // Solve for Harmonic Ritz Values:
  lapack.GEEV('N','V',blockSize_,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);

  if(info != 0){
  printer->stream(IterationDetails) << "Warning!! Anasazi::RandomSolver GEEV solve possible failure: info = " << info << std::endl;
  }
  //std::cout << "DEBUG: Past small eval probl." << std::endl;
  // Compute the eigenvalues and eigenvectors from the original eigenproblem
  Eigensolution<ScalarType,MV> sol;
  sol.numVecs = problem_->getNEV();
  sol.Evals.resize(sol.numVecs);
  // sort the eigenvalues and permute the eigenvectors appropriately
  std::vector<int> order(blockSize_);
  sorter->sort(evals_real,evals_imag,Teuchos::rcpFromRef(order),blockSize_);
  //std::cout << "DEBUG: past sort statement." << std::endl;
  for( int i = 0; i < sol.numVecs; i++){
    sol.Evals[i].realpart = evals_real[i];
    sol.Evals[i].imagpart = evals_imag[i];
  }
  //std::cout << "DEBUG: finished sorting evals." << std::endl;
  // Project Evects back up to large problem. 
  MVT::MvTimesMatAddMv(ONE,*randVecs,evects,0.0,*TmpVecs);
  //std::cout << "DEBUG: Past computing full evects. " << std::endl;

//------Post-Solve Processing----------------------------
  // now permute the eigenvectors according to order
  SolverUtils<ScalarType,MV,OP> msutils;
  msutils.permuteVectors(blockSize_,order,*TmpVecs);
  //Copy only the eigenvectors we asked for to soln. 
  sol.Evecs = MVT::CloneCopy(*TmpVecs, Teuchos::Range1D(0,sol.numVecs-1));

  //std::cout << "DEBUG: Past permuting eigenvectors. " << std::endl;
  // print final summary
  //lobpcg_solver->currentStatus(printer->stream(FinalSummary));

  // print timing information
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  if ( printer->isVerbosity( TimingDetails ) ) {
    Teuchos::TimeMonitor::summarize( printer->stream( TimingDetails ) );
  }
#endif
  //std::cout << "DEBUG: In eigensolver; setting solution." << std::endl;
  // send the solution to the eigenproblem
  problem_->setSolution(sol);
  printer->stream(Debug) << "Returning " << sol.numVecs << " eigenpairs to eigenproblem." << std::endl;
  //std::cout << "DEBUG: Set solution successfully. Returning from eigensolver." << std::endl;

  // return from SolMgr::solve()
  //if (sol.numVecs < nev) return Unconverged;
  } //End solve timer
  return Converged;
}

} // end Experimental  namespace
} // end Anasazi namespace

#endif /* ANASAZI_RANDOMIZED_SOLMGR_HPP */
