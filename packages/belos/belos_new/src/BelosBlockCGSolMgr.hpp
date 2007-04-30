// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef BELOS_BLOCK_CG_SOLMGR_HPP
#define BELOS_BLOCK_CG_SOLMGR_HPP

/*! \file BelosBlockCGSolMgr.hpp
 *  \brief The Belos::BlockCGSolMgr provides a solver manager for the BlockCG linear solver.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosCGIter.hpp"
#include "BelosDGKSOrthoManager.hpp"
#include "BelosICGSOrthoManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutput.hpp"
#include "BelosOutputManager.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_TimeMonitor.hpp"

/** \example BlockCG/BlockCGEpetraEx.cpp
    This is an example of how to use the Belos::BlockCGSolMgr solver manager.
*/

/*! \class Belos::BlockCGSolMgr
 *
 *  \brief The Belos::BlockCGSolMgr provides a powerful and fully-featured solver manager over the BlockCG linear solver.

 \ingroup belos_solver_framework

 \author Heidi Thornquist, Chris Baker, and Teri Barth
 */

namespace Belos {
  
  //! @name BlockCGSolMgr Exceptions
  //@{
  
  /** \brief BlockCGSolMgrLinearProblemFailure is thrown when the linear problem is
   * not setup (i.e. setProblem() was not called) when solve() is called.
   *
   * This exception is thrown from the BlockCGSolMgr::solve() method.
   *
   */
  class BlockCGSolMgrLinearProblemFailure : public BelosError {public:
    BlockCGSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  /** \brief BlockCGSolMgrOrthoFailure is thrown when the orthogonalization manager is
   * unable to generate orthonormal columns from the initial basis vectors.
   *
   * This exception is thrown from the BlockCGSolMgr::solve() method.
   *
   */
  class BlockCGSolMgrOrthoFailure : public BelosError {public:
    BlockCGSolMgrOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  template<class ScalarType, class MV, class OP>
  class BlockCGSolMgr : public SolverManager<ScalarType,MV,OP> {
    
  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;
    
  public:
    
    //! @name Constructors/Destructor
    //@{ 
    
    /*! \brief Basic constructor for BlockCGSolMgr.
     *
     * This constructor accepts the LinearProblem to be solved in addition
     * to a parameter list of options for the solver manager. These options include the following:
     *   - "Block Size" - a \c int specifying the block size to be used by the underlying block Krylov-Schur solver. Default: 1
     *   - "Adaptive Block Size" - a \c bool specifying whether the block size can be modified throughout the solve. Default: true
     *   - "Num Blocks" - a \c int specifying the number of blocks allocated for the Krylov basis. Default: 3*nev
     *   - "Maximum Iterations" - a \c int specifying the maximum number of iterations the underlying solver is allowed to perform. Default: 300
     *   - "Maximum Restarts" - a \c int specifying the maximum number of restarts the underlying solver is allowed to perform. Default: 20
     *   - "Orthogonalization" - a \c string specifying the desired orthogonalization:  DGKS and ICGS. Default: "DGKS"
     *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Belos::Errors
     *   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: machine precision.
     *   - "Relative Convergence Tolerance" - a \c bool specifying whether residuals norms should be scaled for the purposing of deciding convergence. Default: true
     */
    BlockCGSolMgr( const Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> > &problem,
		      Teuchos::ParameterList &pl );
    
    //! Destructor.
    virtual ~BlockCGSolMgr() {};
    //@}
    
    //! @name Accessor methods
    //@{ 
    
    const LinearProblem<ScalarType,MV,OP>& getProblem() const {
      return *problem_;
    }
    
    /*! \brief Return the timers for this object. 
     *
     * The timers are ordered as follows:
     *   - time spent in solve() routine
     */
    Teuchos::Array<Teuchos::RefCountPtr<Teuchos::Time> > getTimers() const {
      return tuple(timerSolve_);
    }
    
    //@}
    
    //! @name Set methods
    //@{
    
    void setProblem( const Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> > &problem ) { problem_ = problem; }
    
    void setParameters( Teuchos::ParameterList &list ) {}
    
    //@}
    
    //! @name Solver application methods
    //@{ 
    
    /*! \brief This method performs possibly repeated calls to the underlying linear solver's iterate() routine
     * until the problem has been solved (as decided by the solver manager) or the solver manager decides to 
     * quit.
     *
     * This method calls BlockCGIter::iterate(), which will return either because a specially constructed status test evaluates to 
     * ::Passed or an exception is thrown.
     *
     * A return from BlockCGIter::iterate() signifies one of the following scenarios:
     *    - the maximum number of restarts has been exceeded. In this scenario, the current solutions to the linear system
     *      will be placed in the linear problem and return ::Unconverged.
     *    - global convergence has been met. In this case, the current solutions to the linear system will be placed in the linear
     *      problem and the solver manager will return ::Converged
     *
     * \returns ::ReturnType specifying:
     *     - ::Converged: the linear problem was solved to the specification required by the solver manager.
     *     - ::Unconverged: the linear problem was not solved to the specification desired by the solver manager.
     */
    ReturnType solve();
    
    //@}
    
    /** \name Overridden from Teuchos::Describable */
    //@{
    
    /** \brief Method to return description of the block CG solver manager */
    std::string description() const;
    
    //@}
    
  private:
    // Linear problem.
    Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> > problem_;
    
    // Output manager.
    Teuchos::RefCountPtr<OutputManager<ScalarType> > printer_;

    // Status test.
    Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > maxIterTest_;
    Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > convTest_;
    Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > sTest_;
    Teuchos::RefCountPtr<StatusTestOutput<ScalarType,MV,OP> > outputTest_;

    // Orthogonalization manager.
    Teuchos::RefCountPtr<MatOrthoManager<ScalarType,MV,OP> > ortho_; 
    
    string orthoType_; 
    MagnitudeType ortho_kappa_;
    
    MagnitudeType convtol_;
    int maxRestarts_, maxIters_;
    bool adaptiveBlockSize_;
    int blockSize_, numBlocks_;
    int verbosity_, output_freq_;
    
    // Timers.
    Teuchos::RefCountPtr<Teuchos::Time> timerSolve_;
  };


// Constructor
template<class ScalarType, class MV, class OP>
BlockCGSolMgr<ScalarType,MV,OP>::BlockCGSolMgr( 
						     const Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> > &problem,
						     Teuchos::ParameterList &pl ) : 
  problem_(problem),
  orthoType_("DGKS"),
  ortho_kappa_(-1.0),
  convtol_(0),
  maxRestarts_(20),
  adaptiveBlockSize_(true),
  blockSize_(0),
  numBlocks_(0),
  verbosity_(Belos::Errors),
  output_freq_(-1),
  timerSolve_(Teuchos::TimeMonitor::getNewTimer("BlockCGSolMgr::solve()"))
{
  TEST_FOR_EXCEPTION(problem_ == Teuchos::null, std::invalid_argument, "Problem not given to solver manager.");
  
  // convergence tolerance
  convtol_ = pl.get("Convergence Tolerance",MT::prec());
  
  // maximum number of restarts
  maxRestarts_ = pl.get("Maximum Restarts",maxRestarts_);
  
  // maximum number of iterations
  maxIters_ = pl.get("Maximum Iterations",maxIters_);

  // block size: default is 1
  blockSize_ = pl.get("Block Size",1);
  TEST_FOR_EXCEPTION(blockSize_ <= 0, std::invalid_argument,
                     "Belos::BlockCGSolMgr: \"Block Size\" must be strictly positive.");
  adaptiveBlockSize_ = pl.get("Adaptive Block Size",adaptiveBlockSize_);

  numBlocks_ = pl.get("Num Blocks",25);
  TEST_FOR_EXCEPTION(numBlocks_ <= 0, std::invalid_argument,
                     "Belos::BlockCGSolMgr: \"Num Blocks\" must be strictly positive.");

  // which orthogonalization to use
  orthoType_ = pl.get("Orthogonalization",orthoType_);
  if (orthoType_ != "DGKS" && orthoType_ != "ICGS") {
    orthoType_ = "DGKS";
  }

  // which orthogonalization constant to use
  ortho_kappa_ = pl.get("Orthogonalization Constant",ortho_kappa_);
  
  // verbosity level
  if (pl.isParameter("Verbosity")) {
    if (Teuchos::isParameterType<int>(pl,"Verbosity")) {
      verbosity_ = pl.get("Verbosity", verbosity_);
    } else {
      verbosity_ = (int)Teuchos::getParameter<Belos::MsgType>(pl,"Verbosity");
    }
  }
  
  // frequency level
  if (verbosity_ & Belos::StatusTestDetails) {
    if (pl.isParameter("Output Frequency")) {
      output_freq_ = pl.get("Output Frequency", output_freq_);
    }
  }

  // Create output manager
  printer_ = Teuchos::rcp( new OutputManager<ScalarType>(verbosity_) );
  
  // Create status tests

  // Convergence
  typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
  typedef Belos::StatusTestResNorm<ScalarType,MV,OP>  StatusTestResNorm_t;
  
  // Basic test checks maximum iterations and native residual.
  maxIterTest_ = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>( maxIters_ ) );

  // Implicit residual test, using the native residual to determine if convergence was achieved.
  Teuchos::RefCountPtr<StatusTestResNorm_t> impConvTest = Teuchos::rcp( new StatusTestResNorm_t( convtol_ ) );
  impConvTest->defineScaleForm( StatusTestResNorm_t::NormOfPrecInitRes, Belos::TwoNorm );
  
  // Explicit residual test once the native residual is below the tolerance
  Teuchos::RefCountPtr<StatusTestResNorm_t> expConvTest  = Teuchos::rcp( new StatusTestResNorm_t( convtol_ ) );
  expConvTest->defineResForm( StatusTestResNorm_t::Explicit, Belos::TwoNorm );
  
  convTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::SEQ, impConvTest, expConvTest ) );

  sTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::OR, maxIterTest_, convTest_ ) );
  
  if (output_freq_ > 0) {
    outputTest_ = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer_, 
                                                                        sTest_, 
									output_freq_, 
									Passed+Failed+Undefined ) ); 
  }
  else {
    outputTest_ = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer_, 
									sTest_, 1 ) );
  }

  // Create orthogonalization manager
  if (orthoType_=="DGKS") {
    if (ortho_kappa_ <= 0) {
      ortho_ = Teuchos::rcp( new DGKSOrthoManager<ScalarType,MV,OP>() );
    }
    else {
      ortho_ = Teuchos::rcp( new DGKSOrthoManager<ScalarType,MV,OP>() );
      Teuchos::rcp_dynamic_cast<DGKSOrthoManager<ScalarType,MV,OP> >(ortho_)->setDepTol( ortho_kappa_ );
    }
  }
  else if (orthoType_=="ICGS") {
    ortho_ = Teuchos::rcp( new ICGSOrthoManager<ScalarType,MV,OP>() );
  } 
  else {
    TEST_FOR_EXCEPTION(orthoType_!="ICGS"&&orthoType_!="DGKS",std::logic_error,
		       "Belos::BlockCGSolMgr(): Invalid orthogonalization type.");
  }  

}

  
// solve()
template<class ScalarType, class MV, class OP>
ReturnType BlockCGSolMgr<ScalarType,MV,OP>::solve() {

  Teuchos::BLAS<int,ScalarType> blas;
  Teuchos::LAPACK<int,ScalarType> lapack;
  
  TEST_FOR_EXCEPTION(!problem_->isProblemSet(),BlockCGSolMgrLinearProblemFailure,
                     "Belos::BlockCGSolMgr::solve(): Linear problem is not ready, setProblem() has not been called.");

  // Create indices for the linear systems to be solved.
  int startPtr = 0;
  int numRHS2Solve = MVT::GetNumberVecs( *(problem_->getRHS()) );
  int numCurrRHS = ( numRHS2Solve < blockSize_) ? numRHS2Solve : blockSize_;

  std::vector<int> currIdx;
  //  If an adaptive block size is allowed then only the linear systems that need to be solved are solved.
  //  Otherwise, the index set is generated that informs the linear problem that some linear systems are augmented.
  if ( adaptiveBlockSize_ ) {
    blockSize_ = numCurrRHS;
    currIdx.resize( numCurrRHS  );
    for (int i=0; i<numCurrRHS; ++i) 
      { currIdx[i] = startPtr+i; }
    
  }
  else {
    currIdx.resize( blockSize_ );
    for (int i=0; i<numCurrRHS; ++i) 
      { currIdx[i] = startPtr+i; }
    for (int i=numCurrRHS; i<blockSize_; ++i) 
      { currIdx[i] = -1; }
  }

  // Inform the linear problem of the current linear system to solve.
  problem_->setLSIndex( currIdx );

  //////////////////////////////////////////////////////////////////////////////////////
  // Parameter list
  Teuchos::ParameterList plist;
  plist.set("Block Size",blockSize_);
  plist.set("Num Blocks",numBlocks_);
  
  // Reset the status test.  
  outputTest_->reset();

  // Assume convergence is achieved, then let any failed convergence set this to false.
  bool isConverged = true;	

  //////////////////////////////////////////////////////////////////////////////////////
  // BlockCG solver

  Teuchos::RefCountPtr<CGIter<ScalarType,MV,OP> > block_cg_iter
    = Teuchos::rcp( new CGIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,plist) );
  

  // Enter solve() iterations
  {
    Teuchos::TimeMonitor slvtimer(*timerSolve_);

    while ( numRHS2Solve > 0 ) {

      // Reset the number of iterations.
      block_cg_iter->resetNumIters();

      // Reset the number of calls that the status test output knows about.
      outputTest_->resetNumCalls();

      // Create the first block in the current Krylov basis.
      Teuchos::RefCountPtr<MV> r_0 = problem_->getCurrResVec();

      // Set the new state and initialize the solver.
      CGIterState<ScalarType,MV> newstate;
      newstate.r = r_0;
      block_cg_iter->initialize(newstate);

      while(1) {
	
	// tell block_cg_iter to iterate
	try {
	  block_cg_iter->iterate();
	  
	  ////////////////////////////////////////////////////////////////////////////////////
	  //
	  // check convergence first
	  //
	  ////////////////////////////////////////////////////////////////////////////////////
	  if ( convTest_->getStatus() == Passed ) {
	    // we have convergence
	    break;  // break from while(1){block_cg_iter->iterate()}
	  }
	  ////////////////////////////////////////////////////////////////////////////////////
	  //
	  // check for maximum iterations
	  //
	  ////////////////////////////////////////////////////////////////////////////////////
	  else if ( maxIterTest_->getStatus() == Passed ) {
	    // we don't have convergence
	    isConverged = false;
	    break;  // break from while(1){block_cg_iter->iterate()}
	  }

	  ////////////////////////////////////////////////////////////////////////////////////
	  //
	  // we returned from iterate(), but none of our status tests Passed.
	  // something is wrong, and it is probably our fault.
	  //
	  ////////////////////////////////////////////////////////////////////////////////////

	  else {
	    TEST_FOR_EXCEPTION(true,std::logic_error,
			       "Belos::BlockCGSolMgr::solve(): Invalid return from BlockCGIter::iterate().");
	  }
	}
	catch (std::exception e) {
	  printer_->stream(Errors) << "Error! Caught exception in BlockCGIter::iterate() at iteration " 
				  << block_cg_iter->getNumIters() << endl 
				  << e.what() << endl;
	  throw;
	}
      }
      
      // Compute the current solution.
      // Update the linear problem.
      Teuchos::RefCountPtr<MV> update = block_cg_iter->getCurrentUpdate();
      problem_->updateSolution( update, true );

      // Inform the linear problem that we are finished with this block linear system.
      problem_->setCurrLS();
      
      // Update indices for the linear systems to be solved.
      startPtr += numCurrRHS;
      numRHS2Solve -= numCurrRHS;
      if ( numRHS2Solve > 0 ) {
	numCurrRHS = ( numRHS2Solve < blockSize_) ? numRHS2Solve : blockSize_;

	if ( adaptiveBlockSize_ ) {
          blockSize_ = numCurrRHS;
	  currIdx.resize( numCurrRHS  );
	  for (int i=0; i<numCurrRHS; ++i) 
	    { currIdx[i] = startPtr+i; }	  
	}
	else {
	  currIdx.resize( blockSize_ );
	  for (int i=0; i<numCurrRHS; ++i) 
	    { currIdx[i] = startPtr+i; }
	  for (int i=numCurrRHS; i<blockSize_; ++i) 
	    { currIdx[i] = -1; }
	}
	// Set the next indices.
	problem_->setLSIndex( currIdx );
      }
      else {
        currIdx.resize( numRHS2Solve );
      }
      
    }// while ( numRHS2Solve > 0 )
    
  }
 
  // print timing information
  Teuchos::TimeMonitor::summarize( printer_->stream(TimingDetails) );
  
  if (!isConverged) {
    return Unconverged; // return from BlockCGSolMgr::solve() 
  }
  return Converged; // return from BlockCGSolMgr::solve() 
}

//  This method requires the solver manager to return a string that describes itself.
template<class ScalarType, class MV, class OP>
std::string BlockCGSolMgr<ScalarType,MV,OP>::description() const
{
  std::ostringstream oss;
  oss << "Belos::BlockCGSolMgr<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
  oss << "{";
  oss << "Ortho Type='"<<orthoType_<<"\'";
  oss << "}";
  return oss.str();
}
  
} // end Belos namespace

#endif /* BELOS_BLOCK_CG_SOLMGR_HPP */
