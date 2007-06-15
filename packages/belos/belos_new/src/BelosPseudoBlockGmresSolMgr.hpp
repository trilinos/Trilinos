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

#ifndef BELOS_PSEUDO_BLOCK_GMRES_SOLMGR_HPP
#define BELOS_PSEUDO_BLOCK_GMRES_SOLMGR_HPP

/*! \file BelosPseudoBlockGmresSolMgr.hpp
 *  \brief The Belos::PseudoBlockGmresSolMgr provides a solver manager for the BlockGmres linear solver.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosPseudoBlockGmresIter.hpp"
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

/** \example BlockGmres/BlockGmresEpetraEx.cpp
    This is an example of how to use the Belos::PseudoBlockGmresSolMgr solver manager.
*/

/*! \class Belos::PseudoBlockGmresSolMgr
 *
 *  \brief The Belos::PseudoBlockGmresSolMgr provides a powerful and fully-featured solver manager over the BlockGmres linear solver.

 \ingroup belos_solver_framework

 \author Heidi Thornquist, Chris Baker, and Teri Barth
 */

namespace Belos {
  
  //! @name PseudoBlockGmresSolMgr Exceptions
  //@{
  
  /** \brief PseudoBlockGmresSolMgrLinearProblemFailure is thrown when the linear problem is
   * not setup (i.e. setProblem() was not called) when solve() is called.
   *
   * This exception is thrown from the PseudoBlockGmresSolMgr::solve() method.
   *
   */
  class PseudoBlockGmresSolMgrLinearProblemFailure : public BelosError {public:
    PseudoBlockGmresSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  /** \brief PseudoBlockGmresSolMgrOrthoFailure is thrown when the orthogonalization manager is
   * unable to generate orthonormal columns from the initial basis vectors.
   *
   * This exception is thrown from the PseudoBlockGmresSolMgr::solve() method.
   *
   */
  class PseudoBlockGmresSolMgrOrthoFailure : public BelosError {public:
    PseudoBlockGmresSolMgrOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  template<class ScalarType, class MV, class OP>
  class PseudoBlockGmresSolMgr : public SolverManager<ScalarType,MV,OP> {
    
  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;
    
  public:
    
    //! @name Constructors/Destructor
    //@{ 
    
    /*! \brief Basic constructor for PseudoBlockGmresSolMgr.
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
    PseudoBlockGmresSolMgr( const Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> > &problem,
		      Teuchos::ParameterList &pl );
    
    //! Destructor.
    virtual ~PseudoBlockGmresSolMgr() {};
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
     * This method calls PseudoBlockGmresIter::iterate(), which will return either because a specially constructed status test evaluates to 
     * ::Passed or an exception is thrown.
     *
     * A return from PseudoBlockGmresIter::iterate() signifies one of the following scenarios:
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
    
    /** \brief Method to return description of the block GMRES solver manager */
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
    Teuchos::RefCountPtr<StatusTestResNorm<ScalarType,MV,OP> > impConvTest_, expConvTest_;
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
    int defQuorum_;
    typename StatusTestResNorm<ScalarType,MV,OP>::ScaleType impResScale_, expResScale_;       
 
    // Timers.
    string label_;
    Teuchos::RefCountPtr<Teuchos::Time> timerSolve_;
  };


// Constructor
template<class ScalarType, class MV, class OP>
PseudoBlockGmresSolMgr<ScalarType,MV,OP>::PseudoBlockGmresSolMgr( 
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
  defQuorum_(1),
  impResScale_(Belos::StatusTestResNorm<ScalarType,MV,OP>::NormOfPrecInitRes),
  expResScale_(Belos::StatusTestResNorm<ScalarType,MV,OP>::NormOfInitRes),
  label_("Belos")
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
                     "Belos::PseudoBlockGmresSolMgr: \"Block Size\" must be strictly positive.");
  adaptiveBlockSize_ = pl.get("Adaptive Block Size",adaptiveBlockSize_);

  numBlocks_ = pl.get("Num Blocks",25);
  TEST_FOR_EXCEPTION(numBlocks_ <= 0, std::invalid_argument,
                     "Belos::PseudoBlockGmresSolMgr: \"Num Blocks\" must be strictly positive.");

  label_ = pl.get("Timer Label", label_);
  string solveLabel = label_ + ": PseudoBlockGmresSolMgr total solve time";
  timerSolve_ = Teuchos::TimeMonitor::getNewTimer(solveLabel);
  
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

  // Get the deflation quorum, or number of converged systems before deflation is allowed
  if (pl.isParameter("Deflation Quorum")) {
    defQuorum_ = pl.get("Deflation Quorum", defQuorum_);
  }
  TEST_FOR_EXCEPTION(defQuorum_ > blockSize_, std::invalid_argument,
                     "Belos::PseudoBlockGmresSolMgr: \"Deflation Quorum\" cannot be larger than \"Block Size\".");

  // Convergence
  typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
  typedef Belos::StatusTestResNorm<ScalarType,MV,OP>  StatusTestResNorm_t;

  if (pl.isParameter("Implicit Residual Scaling")) {
    string impResScalingType = Teuchos::getParameter<string>( pl, "Implicit Residual Scaling" );
    if (impResScalingType == "Norm of Initial Residual") {
      impResScale_ = StatusTestResNorm_t::NormOfInitRes;
    } else if (impResScalingType == "Norm of Preconditioned Initial Residual") {
      impResScale_ = StatusTestResNorm_t::NormOfPrecInitRes;
    } else if (impResScalingType == "Norm of RHS") {
      impResScale_ = StatusTestResNorm_t::NormOfRHS;
    } else if (impResScalingType == "None") {
      impResScale_ = StatusTestResNorm_t::None;
    } else 
      TEST_FOR_EXCEPTION( true ,std::logic_error,
			  "Belos::PseudoBlockGmresSolMgr(): Invalid implicit residual scaling type.");
  }
  
  if (pl.isParameter("Explicit Residual Scaling")) {
    string expResScalingType = Teuchos::getParameter<string>( pl, "Explicit Residual Scaling" );
    if (expResScalingType == "Norm of Initial Residual") {
      expResScale_ = StatusTestResNorm_t::NormOfInitRes;
    } else if (expResScalingType == "Norm of Preconditioned Initial Residual") {
      expResScale_ = StatusTestResNorm_t::NormOfPrecInitRes;
    } else if (expResScalingType == "Norm of RHS") {
      expResScale_ = StatusTestResNorm_t::NormOfRHS;
    } else if (expResScalingType == "None") {
      expResScale_ = StatusTestResNorm_t::None;
    } else 
      TEST_FOR_EXCEPTION( true ,std::logic_error,
			  "Belos::PseudoBlockGmresSolMgr(): Invalid explicit residual scaling type.");
  }
  
  // Basic test checks maximum iterations and native residual.
  maxIterTest_ = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>( maxIters_ ) );

  // Implicit residual test, using the native residual to determine if convergence was achieved.
  impConvTest_ = Teuchos::rcp( new StatusTestResNorm_t( convtol_, defQuorum_ ) );
  impConvTest_->defineScaleForm( impResScale_, Belos::TwoNorm );
  
  // Explicit residual test once the native residual is below the tolerance
  expConvTest_  = Teuchos::rcp( new StatusTestResNorm_t( convtol_, defQuorum_ ) );
  expConvTest_->defineResForm( StatusTestResNorm_t::Explicit, Belos::TwoNorm );
  expConvTest_->defineScaleForm( expResScale_, Belos::TwoNorm );
  
  convTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::SEQ, impConvTest_, expConvTest_ ) );

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
      ortho_ = Teuchos::rcp( new DGKSOrthoManager<ScalarType,MV,OP>( label_ ) );
    }
    else {
      ortho_ = Teuchos::rcp( new DGKSOrthoManager<ScalarType,MV,OP>( label_ ) );
      Teuchos::rcp_dynamic_cast<DGKSOrthoManager<ScalarType,MV,OP> >(ortho_)->setDepTol( ortho_kappa_ );
    }
  }
  else if (orthoType_=="ICGS") {
    ortho_ = Teuchos::rcp( new ICGSOrthoManager<ScalarType,MV,OP>( label_ ) );
  } 
  else {
    TEST_FOR_EXCEPTION(orthoType_!="ICGS"&&orthoType_!="DGKS",std::logic_error,
		       "Belos::PseudoBlockGmresSolMgr(): Invalid orthogonalization type.");
  }  

}

  
// solve()
template<class ScalarType, class MV, class OP>
ReturnType PseudoBlockGmresSolMgr<ScalarType,MV,OP>::solve() {

  Teuchos::BLAS<int,ScalarType> blas;
  Teuchos::LAPACK<int,ScalarType> lapack;
  
  TEST_FOR_EXCEPTION(!problem_->isProblemSet(),PseudoBlockGmresSolMgrLinearProblemFailure,
                     "Belos::PseudoBlockGmresSolMgr::solve(): Linear problem is not ready, setProblem() has not been called.");

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
  plist.set("Num Blocks",numBlocks_);
  
  // Reset the status test.  
  outputTest_->reset();

  // Assume convergence is achieved, then let any failed convergence set this to false.
  bool isConverged = true;	

  //////////////////////////////////////////////////////////////////////////////////////
  // BlockGmres solver

  Teuchos::RefCountPtr<PseudoBlockGmresIter<ScalarType,MV,OP> > block_gmres_iter
    = Teuchos::rcp( new PseudoBlockGmresIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,ortho_,plist) );  

  // Enter solve() iterations
  {
    Teuchos::TimeMonitor slvtimer(*timerSolve_);

    while ( numRHS2Solve > 0 ) {

      // Reset the active / converged vectors from this block
      std::vector<int> convRHSIdx;
      std::vector<int> currRHSIdx( currIdx );
      currRHSIdx.resize(numCurrRHS);

      // Set the current number of blocks with the pseudo Gmres iteration.
      block_gmres_iter->setNumBlocks( numBlocks_ );

      // Reset the number of iterations.
      block_gmres_iter->resetNumIters();

      // Reset the number of calls that the status test output knows about.
      outputTest_->resetNumCalls();

      // Get a new state struct and initialize the solver.
      PseudoBlockGmresIterState<ScalarType,MV> newState;

      // Create the first block in the current Krylov basis for each right-hand side.
      std::vector<int> index(1);
      Teuchos::RefCountPtr<MV> tmpV, R_0 = MVT::Clone( *(problem_->getCurrResVec()), blockSize_ );
      problem_->computeCurrResVec( &*R_0 );
      newState.V.resize( blockSize_ );
      newState.Z.resize( blockSize_ );
      for (int i=0; i<blockSize_; ++i) {
	index[0]=i;
	tmpV = MVT::CloneCopy( *R_0, index );
	
	// Get a matrix to hold the orthonormalization coefficients.
	Teuchos::RefCountPtr<Teuchos::SerialDenseVector<int,ScalarType> > tmpZ
	  = Teuchos::rcp( new Teuchos::SerialDenseVector<int,ScalarType>( 1 ));
      
	// Orthonormalize the new V_0
	int rank = ortho_->normalize( *tmpV, tmpZ );
	TEST_FOR_EXCEPTION(rank != 1, PseudoBlockGmresSolMgrOrthoFailure,
			   "Belos::PseudoBlockGmresSolMgr::solve(): Failed to compute initial block of orthonormal vectors.");

	newState.V[i] = tmpV;
	newState.Z[i] = tmpZ;
      }

      newState.curDim = 0;
      block_gmres_iter->initialize(newState);
      int numRestarts = 0;

      while(1) {
	
	// tell block_gmres_iter to iterate
	try {
	  block_gmres_iter->iterate();
	  
	  ////////////////////////////////////////////////////////////////////////////////////
	  //
	  // check convergence first
	  //
	  ////////////////////////////////////////////////////////////////////////////////////
	  if ( convTest_->getStatus() == Passed ) {

	    // Figure out which linear systems converged.
	    std::vector<int> convIdx = expConvTest_->convIndices();

	    // If the number of converged linear systems is equal to the
            // number of current linear systems, then we are done with this block.
	    if (convIdx.size() == currRHSIdx.size())
	      break;  // break from while(1){block_gmres_iter->iterate()}

	    // Get a new state struct and initialize the solver.
	    PseudoBlockGmresIterState<ScalarType,MV> newState;

	    // Inform the linear problem that we are finished with this current linear system.
	    problem_->setCurrLS();

	    // Get the state.
	    PseudoBlockGmresIterState<ScalarType,MV> oldState = block_gmres_iter->getState();
	    int curDim = oldState.curDim;
	    std::vector<int> index( curDim );
	    for (int i=0; i<curDim; ++i) { index[i] = i; }	      

	    // Get a new state struct and reset currRHSIdx to have the right-hand sides that 
	    // are left to converge for this block.
	    int have = 0;
	    std::vector<int> oldRHSIdx( currRHSIdx );
	    std::vector<int> defRHSIdx;
	    for (unsigned int i=0; i<currRHSIdx.size(); ++i) {
	      bool found = false;
	      for (unsigned int j=0; j<convIdx.size(); ++j) {
		if (currRHSIdx[i] == convIdx[j]) {
		  found = true;
		  break;
		}
	      }
	      if (found) {
		defRHSIdx.push_back( i );
	      }
	      else {
		newState.V.push_back( Teuchos::rcp_const_cast<MV>( oldState.V[i] ) );
		newState.Z.push_back( Teuchos::rcp_const_cast<Teuchos::SerialDenseVector<int,ScalarType> >( oldState.Z[i] ) );
		newState.H.push_back( Teuchos::rcp_const_cast<Teuchos::SerialDenseMatrix<int,ScalarType> >( oldState.H[i] ) );
		newState.sn.push_back( Teuchos::rcp_const_cast<Teuchos::SerialDenseVector<int,ScalarType> >( oldState.sn[i] ) );
		newState.cs.push_back( Teuchos::rcp_const_cast<Teuchos::SerialDenseVector<int,MagnitudeType> >(oldState.cs[i] ) );
		currRHSIdx[have] = currRHSIdx[i];
		have++;
	      }
	    }
	    defRHSIdx.resize(currRHSIdx.size()-have);
	    currRHSIdx.resize(have);

	    // Compute the current solution that needs to be deflated if this solver has taken any steps.
	    if (curDim) {
	      Teuchos::RefCountPtr<MV> update = block_gmres_iter->getCurrentUpdate();
	      Teuchos::RefCountPtr<MV> defUpdate = MVT::CloneView( *update, defRHSIdx );
	      
	      // Set the deflated indices so we can update the solution.
	      problem_->setLSIndex( convIdx );
	      
	      // Update the linear problem.
	      problem_->updateSolution( defUpdate, true );
	    }
	    
	    // Set the remaining indices after deflation.
	    problem_->setLSIndex( currRHSIdx );
	    
	    // Set the dimension of the subspace, which is the same as the old subspace size.
	    newState.curDim = curDim;
	    
	    // Initialize the solver with the deflated system.
	    block_gmres_iter->initialize(newState);
	  }
	  ////////////////////////////////////////////////////////////////////////////////////
	  //
	  // check for maximum iterations
	  //
	  ////////////////////////////////////////////////////////////////////////////////////
	  else if ( maxIterTest_->getStatus() == Passed ) {
	    // we don't have convergence
	    isConverged = false;
	    break;  // break from while(1){block_gmres_iter->iterate()}
	  }
	  ////////////////////////////////////////////////////////////////////////////////////
	  //
	  // check for restarting, i.e. the subspace is full
	  //
	  ////////////////////////////////////////////////////////////////////////////////////
	  else if ( block_gmres_iter->getCurSubspaceDim() == block_gmres_iter->getMaxSubspaceDim() ) {
	
	    if ( numRestarts >= maxRestarts_ ) {
	      isConverged = false;
	      break; // break from while(1){block_gmres_iter->iterate()}
	    }
	    numRestarts++;

	    printer_->stream(Debug) << " Performing restart number " << numRestarts << " of " << maxRestarts_ << endl << endl;
	    
	    // Update the linear problem.
	    Teuchos::RefCountPtr<MV> update = block_gmres_iter->getCurrentUpdate();
	    problem_->updateSolution( update, true );
	    
	    // Get the state.
	    PseudoBlockGmresIterState<ScalarType,MV> oldState = block_gmres_iter->getState();
	    
	    // Set the new state.
	    PseudoBlockGmresIterState<ScalarType,MV> newstate;
	    newstate.V.resize(currRHSIdx.size());
	    newstate.Z.resize(currRHSIdx.size());

	    // Compute the restart vectors
	    // NOTE: Force the linear problem to update the current residual since the solution was updated.
	    Teuchos::RefCountPtr<MV> R_0 = MVT::Clone( *(problem_->getCurrResVec()), currRHSIdx.size() );
	    problem_->computeCurrResVec( &*R_0 );
	    std::vector<int> index(1);
	    for (unsigned int i=0; i<currRHSIdx.size(); ++i) {
	      index[0] = i;
	
	      tmpV = MVT::CloneCopy( *R_0, index );
	
	      // Get a matrix to hold the orthonormalization coefficients.
	      Teuchos::RefCountPtr<Teuchos::SerialDenseVector<int,ScalarType> > tmpZ
		= Teuchos::rcp( new Teuchos::SerialDenseVector<int,ScalarType>( 1 ));
	      
	      // Orthonormalize the new V_0
	      int rank = ortho_->normalize( *tmpV, tmpZ );
	      TEST_FOR_EXCEPTION(rank != 1 ,PseudoBlockGmresSolMgrOrthoFailure,
				 "Belos::PseudoBlockGmresSolMgr::solve(): Failed to compute initial block of orthonormal vectors after the restart.");
	      
	      newstate.V[i] = tmpV;
	      newstate.Z[i] = tmpZ;
	    }

	    // Initialize the solver.
	    newstate.curDim = 0;
	    block_gmres_iter->initialize(newstate);

	  } // end of restarting

	  ////////////////////////////////////////////////////////////////////////////////////
	  //
	  // we returned from iterate(), but none of our status tests Passed.
	  // something is wrong, and it is probably our fault.
	  //
	  ////////////////////////////////////////////////////////////////////////////////////

	  else {
	    TEST_FOR_EXCEPTION(true,std::logic_error,
			       "Belos::PseudoBlockGmresSolMgr::solve(): Invalid return from PseudoBlockGmresIter::iterate().");
	  }
	}
        catch (PseudoBlockGmresIterOrthoFailure e) {
     
	  // Try to recover the most recent least-squares solution
	  block_gmres_iter->updateLSQR( block_gmres_iter->getCurSubspaceDim() );

	  // Check to see if the most recent least-squares solution yielded convergence.
	  sTest_->checkStatus( &*block_gmres_iter );
	  if (convTest_->getStatus() != Passed)
	    isConverged = false;
	  break;
        }
	catch (std::exception e) {
	  printer_->stream(Errors) << "Error! Caught exception in PseudoBlockGmresIter::iterate() at iteration " 
				  << block_gmres_iter->getNumIters() << endl 
				  << e.what() << endl;
	  throw;
	}
      }
      
      // Compute the current solution.
      // Update the linear problem.
      Teuchos::RefCountPtr<MV> update = block_gmres_iter->getCurrentUpdate();
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
    return Unconverged; // return from PseudoBlockGmresSolMgr::solve() 
  }
  return Converged; // return from PseudoBlockGmresSolMgr::solve() 
}

//  This method requires the solver manager to return a string that describes itself.
template<class ScalarType, class MV, class OP>
std::string PseudoBlockGmresSolMgr<ScalarType,MV,OP>::description() const
{
  std::ostringstream oss;
  oss << "Belos::PseudoBlockGmresSolMgr<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
  oss << "{";
  oss << "Ortho Type='"<<orthoType_<<"\'";
  oss << "}";
  return oss.str();
}
  
} // end Belos namespace

#endif /* BELOS_PSEUDO_BLOCK_GMRES_SOLMGR_HPP */
