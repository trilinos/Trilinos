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

#ifndef BELOS_BLOCK_GMRES_SOLMGR_HPP
#define BELOS_BLOCK_GMRES_SOLMGR_HPP

/*! \file BelosBlockGmresSolMgr.hpp
 *  \brief The Belos::BlockGmresSolMgr provides a solver manager for the BlockGmres linear solver.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosBlockGmresIter.hpp"
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
    This is an example of how to use the Belos::BlockGmresSolMgr solver manager.
*/

/*! \class Belos::BlockGmresSolMgr
 *
 *  \brief The Belos::BlockGmresSolMgr provides a powerful and fully-featured solver manager over the BlockGmres linear solver.

 \ingroup belos_solver_framework

 \author Heidi Thornquist, Chris Baker, and Teri Barth
 */

namespace Belos {
  
  //! @name BlockGmresSolMgr Exceptions
  //@{
  
  /** \brief BlockGmresSolMgrLinearProblemFailure is thrown when the linear problem is
   * not setup (i.e. setProblem() was not called) when solve() is called.
   *
   * This exception is thrown from the BlockGmresSolMgr::solve() method.
   *
   */
  class BlockGmresSolMgrLinearProblemFailure : public BelosError {public:
    BlockGmresSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  /** \brief BlockGmresSolMgrOrthoFailure is thrown when the orthogonalization manager is
   * unable to generate orthonormal columns from the initial basis vectors.
   *
   * This exception is thrown from the BlockGmresSolMgr::solve() method.
   *
   */
  class BlockGmresSolMgrOrthoFailure : public BelosError {public:
    BlockGmresSolMgrOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  template<class ScalarType, class MV, class OP>
  class BlockGmresSolMgr : public SolverManager<ScalarType,MV,OP> {
    
  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;
    
  public:
    
    //! @name Constructors/Destructor
    //@{ 
   
    /*! \brief Empty constructor for BlockGmresSolMgr.
     * This constructor takes no arguments and sets the default values for the solver.
     * The linear problem must be passed in using setProblem() before solve() is called on this object.
     * The solver values can be changed using setParameters().
     */
     BlockGmresSolMgr();
 
    /*! \brief Basic constructor for BlockGmresSolMgr.
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
    BlockGmresSolMgr( const Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> > &problem,
		      const Teuchos::RefCountPtr<Teuchos::ParameterList> &pl );
    
    //! Destructor.
    virtual ~BlockGmresSolMgr() {};
    //@}
    
    //! @name Accessor methods
    //@{ 
    
    /*! \brief Get current linear problem being solved for in this object.
     */
    const LinearProblem<ScalarType,MV,OP>& getProblem() const {
      return *problem_;
    }

    /*! \brief Get a parameter list containing the valid parameters for this object.
     */
    Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const { return defaultParams_; }

    /*! \brief Get a parameter list containing the current parameters for this object.
     */
    Teuchos::RefCountPtr<const Teuchos::ParameterList> getCurrentParameters() const { return params_; }
 
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
    
    void setParameters( const Teuchos::RefCountPtr<Teuchos::ParameterList> &params );
    
    //@}
    
    //! @name Solver application methods
    //@{ 
    
    /*! \brief This method performs possibly repeated calls to the underlying linear solver's iterate() routine
     * until the problem has been solved (as decided by the solver manager) or the solver manager decides to 
     * quit.
     *
     * This method calls BlockGmresIter::iterate(), which will return either because a specially constructed status test evaluates to 
     * ::Passed or an exception is thrown.
     *
     * A return from BlockGmresIter::iterate() signifies one of the following scenarios:
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

    // Method to set the default parameters.
    void setDefaultParams();

    // Method to convert string to enumerated type for residual.
    typename StatusTestResNorm<ScalarType,MV,OP>::ScaleType convertStringToScaleType( string& scaleType ) {
      typedef Belos::StatusTestResNorm<ScalarType,MV,OP>  StatusTestResNorm_t;
      if (scaleType == "Norm of Initial Residual") {
	return StatusTestResNorm_t::NormOfInitRes;
      } else if (scaleType == "Norm of Preconditioned Initial Residual") {
	return StatusTestResNorm_t::NormOfPrecInitRes;
      } else if (scaleType == "Norm of RHS") {
	return StatusTestResNorm_t::NormOfRHS;
      } else if (scaleType == "None") {
	return StatusTestResNorm_t::None;
      } else 
	TEST_FOR_EXCEPTION( true ,std::logic_error,
			    "Belos::BlockGmresSolMgr(): Invalid residual scaling type.");
    }

    // Linear problem.
    Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> > problem_;
    
    // Output manager.
    Teuchos::RefCountPtr<OutputManager<ScalarType> > printer_;
    Teuchos::RefCountPtr<ostream> outputStream_;

    // Status test.
    Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > sTest_;
    Teuchos::RefCountPtr<StatusTestMaxIters<ScalarType,MV,OP> > maxIterTest_;
    Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > convTest_;
    Teuchos::RefCountPtr<StatusTestResNorm<ScalarType,MV,OP> > expConvTest_, impConvTest_;
    Teuchos::RefCountPtr<StatusTestOutput<ScalarType,MV,OP> > outputTest_;

    // Orthogonalization manager.
    Teuchos::RefCountPtr<MatOrthoManager<ScalarType,MV,OP> > ortho_; 
    
    // Current parameter list.
    Teuchos::RefCountPtr<ParameterList> params_, defaultParams_;

    // Default solver values.
    static const MagnitudeType convtol_default_;
    static const MagnitudeType orthoKappa_default_;
    static const int maxRestarts_default_;
    static const int maxIters_default_;
    static const bool adaptiveBlockSize_default_;
    static const bool showMaxResNormOnly_default_;
    static const int blockSize_default_;
    static const int numBlocks_default_;
    static const int verbosity_default_;
    static const int outputFreq_default_;
    static const std::string impResScale_default_; 
    static const std::string expResScale_default_; 
    static const std::string label_default_;
    static const std::string orthoType_default_;
    static const Teuchos::RefCountPtr<ostream> outputStream_default_;

    // Current solver values.
    MagnitudeType convtol_, orthoKappa_;
    int maxRestarts_, maxIters_;
    int blockSize_, numBlocks_, verbosity_, outputFreq_;
    bool adaptiveBlockSize_, showMaxResNormOnly_;
    std::string orthoType_; 
    std::string impResScale_, expResScale_;
    
    // Timers.
    std::string label_;
    Teuchos::RefCountPtr<Teuchos::Time> timerSolve_;

    // Internal state variables.
    bool isSet_;
  };


// Default solver values.
template<class ScalarType, class MV, class OP>
const typename Teuchos::ScalarTraits<ScalarType>::magnitudeType BlockGmresSolMgr<ScalarType,MV,OP>::convtol_default_ = 1e-8;

template<class ScalarType, class MV, class OP>
const typename Teuchos::ScalarTraits<ScalarType>::magnitudeType BlockGmresSolMgr<ScalarType,MV,OP>::orthoKappa_default_ = -1.0;

template<class ScalarType, class MV, class OP>
const int BlockGmresSolMgr<ScalarType,MV,OP>::maxRestarts_default_ = 20;

template<class ScalarType, class MV, class OP>
const int BlockGmresSolMgr<ScalarType,MV,OP>::maxIters_default_ = 1000;

template<class ScalarType, class MV, class OP>
const bool BlockGmresSolMgr<ScalarType,MV,OP>::adaptiveBlockSize_default_ = true;

template<class ScalarType, class MV, class OP>
const bool BlockGmresSolMgr<ScalarType,MV,OP>::showMaxResNormOnly_default_ = false;

template<class ScalarType, class MV, class OP>
const int BlockGmresSolMgr<ScalarType,MV,OP>::blockSize_default_ = 1;

template<class ScalarType, class MV, class OP>
const int BlockGmresSolMgr<ScalarType,MV,OP>::numBlocks_default_ = 25;

template<class ScalarType, class MV, class OP>
const int BlockGmresSolMgr<ScalarType,MV,OP>::verbosity_default_ = Belos::Errors;

template<class ScalarType, class MV, class OP>
const int BlockGmresSolMgr<ScalarType,MV,OP>::outputFreq_default_ = -1;

template<class ScalarType, class MV, class OP>
const std::string BlockGmresSolMgr<ScalarType,MV,OP>::impResScale_default_ = "Norm of Preconditioned Initial Residual";

template<class ScalarType, class MV, class OP>
const std::string BlockGmresSolMgr<ScalarType,MV,OP>::expResScale_default_ = "Norm of Initial Residual";

template<class ScalarType, class MV, class OP>
const std::string BlockGmresSolMgr<ScalarType,MV,OP>::label_default_ = "Belos";

template<class ScalarType, class MV, class OP>
const std::string BlockGmresSolMgr<ScalarType,MV,OP>::orthoType_default_ = "DGKS";

template<class ScalarType, class MV, class OP>
const Teuchos::RefCountPtr<ostream> BlockGmresSolMgr<ScalarType,MV,OP>::outputStream_default_ = Teuchos::rcp(&std::cout,false);


// Empty Constructor
template<class ScalarType, class MV, class OP>
BlockGmresSolMgr<ScalarType,MV,OP>::BlockGmresSolMgr() :
  outputStream_(outputStream_default_),
  convtol_(convtol_default_),
  orthoKappa_(orthoKappa_default_),
  maxRestarts_(maxRestarts_default_),
  maxIters_(maxIters_default_),
  blockSize_(blockSize_default_),
  numBlocks_(numBlocks_default_),
  verbosity_(verbosity_default_),
  outputFreq_(outputFreq_default_),
  adaptiveBlockSize_(adaptiveBlockSize_default_),
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  orthoType_(orthoType_default_),
  impResScale_(impResScale_default_),
  expResScale_(expResScale_default_),
  label_(label_default_),
  isSet_(false)
{
  // Set the default parameter list.
  setDefaultParams();

  // Set the current parameter list to the default parameter list, don't set parameters 
  // just in case the user decides to set them later with a call to setParameters().
  params_ = defaultParams_;
}


// Basic Constructor
template<class ScalarType, class MV, class OP>
BlockGmresSolMgr<ScalarType,MV,OP>::BlockGmresSolMgr( 
						     const Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> > &problem,
						     const Teuchos::RefCountPtr<Teuchos::ParameterList> &pl ) : 
  problem_(problem),
  outputStream_(outputStream_default_),
  convtol_(convtol_default_),
  orthoKappa_(orthoKappa_default_),
  maxRestarts_(maxRestarts_default_),
  maxIters_(maxIters_default_),
  blockSize_(blockSize_default_),
  numBlocks_(numBlocks_default_),
  verbosity_(verbosity_default_),
  outputFreq_(outputFreq_default_),
  adaptiveBlockSize_(adaptiveBlockSize_default_),
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  orthoType_(orthoType_default_),
  impResScale_(impResScale_default_),
  expResScale_(expResScale_default_),
  label_(label_default_),
  isSet_(false)
{
  TEST_FOR_EXCEPTION(problem_ == Teuchos::null, std::invalid_argument, "Problem not given to solver manager.");

  // Set the default parameter list.
  setDefaultParams();

  // If the parameter list pointer is null, then set the current parameters to the default parameter list.
  if (pl == Teuchos::null) {
    params_ = defaultParams_;
  }
  else {
    // Set the parameters using the list that was passed in.
    setParameters( pl );  
  }
}


template<class ScalarType, class MV, class OP>
void BlockGmresSolMgr<ScalarType,MV,OP>::setParameters( const Teuchos::RefCountPtr<Teuchos::ParameterList> &params )
{
  // Create the internal parameter list if ones doesn't already exist.
  if (params_ == Teuchos::null) {
    params_ = Teuchos::rcp( new Teuchos::ParameterList() );
  }

  // Check for maximum number of restarts
  if (params->isParameter("Maximum Restarts")) {
    maxRestarts_ = params->get("Maximum Restarts",maxRestarts_default_);

    // Update parameter in our list.
    params_->set("Maximum Restarts", maxRestarts_);
  }

  // Check for maximum number of iterations
  if (params->isParameter("Maximum Iterations")) {
    maxIters_ = params->get("Maximum Iterations",maxIters_default_);

    // Update parameter in our list and in status test.
    params_->set("Maximum Iterations", maxIters_);
    if (maxIterTest_!=Teuchos::null)
      maxIterTest_->setMaxIters( maxIters_ );
  }

  // Check for blocksize
  if (params->isParameter("Block Size")) {
    blockSize_ = params->get("Block Size",blockSize_default_);    
    TEST_FOR_EXCEPTION(blockSize_ <= 0, std::invalid_argument,
		       "Belos::BlockGmresSolMgr: \"Block Size\" must be strictly positive.");

    // Update parameter in our list.
    params_->set("Block Size", blockSize_);
  }

  // Check if the blocksize should be adaptive
  if (params->isParameter("Adapative Block Size")) {
    adaptiveBlockSize_ = params->get("Adaptive Block Size",adaptiveBlockSize_default_);
    
    // Update parameter in our list.
    params_->set("Adaptive Block Size", adaptiveBlockSize_);
  }

  // Check for the maximum number of blocks.
  if (params->isParameter("Num Blocks")) {
    numBlocks_ = params->get("Num Blocks",numBlocks_default_);
    TEST_FOR_EXCEPTION(numBlocks_ <= 0, std::invalid_argument,
		       "Belos::BlockGmresSolMgr: \"Num Blocks\" must be strictly positive.");

    // Update parameter in our list.
    params_->set("Num Blocks", numBlocks_);
  }

  // Check to see if the timer label changed.
  if (params->isParameter("Timer Label")) {
    string tempLabel = params->get("Timer Label", label_default_);

    // Update parameter in our list and solver timer
    if (tempLabel != label_) {
      label_ = tempLabel;
      params_->set("Timer Label", label_);
      string solveLabel = label_ + ": BlockGmresSolMgr total solve time";
      timerSolve_ = Teuchos::TimeMonitor::getNewTimer(solveLabel);
    }
  }

  // Check if the orthogonalization changed.
  if (params->isParameter("Orthogonalization")) {
    string tempOrthoType = params->get("Orthogonalization",orthoType_default_);
    TEST_FOR_EXCEPTION( tempOrthoType != "DGKS" && tempOrthoType != "ICGS", std::invalid_argument,
			"Belos::BlockGmresSolMgr: \"Orthogonalization\" must be either \"DGKS\" or \"ICGS\".");
    if (tempOrthoType != orthoType_) {
      orthoType_ = tempOrthoType;
      // Create orthogonalization manager
      if (orthoType_=="DGKS") {
	if (orthoKappa_ <= 0) {
	  ortho_ = Teuchos::rcp( new DGKSOrthoManager<ScalarType,MV,OP>( label_ ) );
	}
	else {
	  ortho_ = Teuchos::rcp( new DGKSOrthoManager<ScalarType,MV,OP>( label_ ) );
	  Teuchos::rcp_dynamic_cast<DGKSOrthoManager<ScalarType,MV,OP> >(ortho_)->setDepTol( orthoKappa_ );
	}
      }
      else if (orthoType_=="ICGS") {
	ortho_ = Teuchos::rcp( new ICGSOrthoManager<ScalarType,MV,OP>( label_ ) );
      } 
    }  
  }

  // Check which orthogonalization constant to use.
  if (params->isParameter("Orthogonalization Constant")) {
    orthoKappa_ = params->get("Orthogonalization Constant",orthoKappa_default_);

    // Update parameter in our list.
    params_->set("Orthogonalization Constant",orthoKappa_);
    if (orthoType_=="DGKS") {
      if (orthoKappa_ > 0 && ortho_ != Teuchos::null) {
	Teuchos::rcp_dynamic_cast<DGKSOrthoManager<ScalarType,MV,OP> >(ortho_)->setDepTol( orthoKappa_ );
      }
    } 
  }

  // Check for a change in verbosity level
  if (params->isParameter("Verbosity")) {
    if (Teuchos::isParameterType<int>(*params,"Verbosity")) {
      verbosity_ = params->get("Verbosity", verbosity_default_);
    } else {
      verbosity_ = (int)Teuchos::getParameter<Belos::MsgType>(*params,"Verbosity");
    }

    // Update parameter in our list.
    params_->set("Verbosity", verbosity_);
    if (printer_ != Teuchos::null)
      printer_->setVerbosity(verbosity_);
  }

  // output stream
  if (params->isParameter("Output Stream")) {
    outputStream_ = Teuchos::getParameter<Teuchos::RefCountPtr<ostream> >(*params,"Output Stream");

    // Update parameter in our list.
    params_->set("Output Stream", outputStream_);
    if (printer_ != Teuchos::null)
      printer_->setOStream( outputStream_ );
  }

  // frequency level
  if (verbosity_ & Belos::StatusTestDetails) {
    if (params->isParameter("Output Frequency")) {
      outputFreq_ = params->get("Output Frequency", outputFreq_default_);
    }

    // Update parameter in out list and output status test.
    params_->set("Output Frequency", outputFreq_);
    if (outputTest_ != Teuchos::null)
      outputTest_->setOutputFrequency( outputFreq_ );
  }

  // Create output manager if we need to.
  if (printer_ == Teuchos::null) {
    printer_ = Teuchos::rcp( new OutputManager<ScalarType>(verbosity_, outputStream_) );
  }  
  
  // Convergence
  typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
  typedef Belos::StatusTestResNorm<ScalarType,MV,OP>  StatusTestResNorm_t;

  // Check for convergence tolerance
  if (params->isParameter("Convergence Tolerance")) {
    convtol_ = params->get("Convergence Tolerance",convtol_default_);

    // Update parameter in our list and residual tests.
    params_->set("Convergence Tolerance", convtol_);
    if (impConvTest_ != Teuchos::null)
      impConvTest_->setTolerance( convtol_ );
    if (expConvTest_ != Teuchos::null)
      expConvTest_->setTolerance( convtol_ );
  }
 
  // Check for a change in scaling, if so we need to build new residual tests.
  bool newImpResTest = false, newExpResTest = false; 
  if (params->isParameter("Implicit Residual Scaling")) {
    string tempImpResScale = Teuchos::getParameter<string>( *params, "Implicit Residual Scaling" );

    // Only update the scaling if it's different.
    if (impResScale_ != tempImpResScale) {
      typename StatusTestResNorm_t::ScaleType impResScaleType = convertStringToScaleType( tempImpResScale );
      impResScale_ = tempImpResScale;

      // Update parameter in our list and residual tests
      params_->set("Implicit Residual Scaling", impResScale_);
      if (impConvTest_ != Teuchos::null) {
        try { 
          impConvTest_->defineScaleForm( impResScaleType, Belos::TwoNorm );
        }
        catch (exception& e) { 
	  newImpResTest = true;
        }
      }
    }      
  }
  
  if (params->isParameter("Explicit Residual Scaling")) {
    string tempExpResScale = Teuchos::getParameter<string>( *params, "Explicit Residual Scaling" );

    // Only update the scaling if it's different.
    if (expResScale_ != tempExpResScale) {
      typename StatusTestResNorm_t::ScaleType expResScaleType = convertStringToScaleType( tempExpResScale );
      expResScale_ = tempExpResScale;

      // Update parameter in our list and residual tests
      params_->set("Explicit Residual Scaling", expResScale_);
      if (expConvTest_ != Teuchos::null) {
        try { 
          expConvTest_->defineScaleForm( expResScaleType, Belos::TwoNorm );
        }
        catch (exception& e) {
	  newExpResTest = true;
        }
      }
    }      
  }


  if (params->isParameter("Show Maximum Residual Norm Only")) {
    showMaxResNormOnly_ = Teuchos::getParameter<bool>(*params,"Show Maximum Residual Norm Only");

    // Update parameter in our list and residual tests
    params_->set("Show Maximum Residual Norm Only", showMaxResNormOnly_);
    if (impConvTest_ != Teuchos::null)
      impConvTest_->setShowMaxResNormOnly( showMaxResNormOnly_ );
    if (expConvTest_ != Teuchos::null)
      expConvTest_->setShowMaxResNormOnly( showMaxResNormOnly_ );
  }

  // Create status tests if we need to.

  // Basic test checks maximum iterations and native residual.
  if (maxIterTest_ == Teuchos::null)
    maxIterTest_ = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>( maxIters_ ) );

  // Implicit residual test, using the native residual to determine if convergence was achieved.
  if (impConvTest_ == Teuchos::null || newImpResTest) {
    impConvTest_ = Teuchos::rcp( new StatusTestResNorm_t( convtol_ ) );
    impConvTest_->defineScaleForm( convertStringToScaleType(impResScale_), Belos::TwoNorm );
    impConvTest_->setShowMaxResNormOnly( showMaxResNormOnly_ );
  }

  // Explicit residual test once the native residual is below the tolerance
  if (expConvTest_ == Teuchos::null || newExpResTest) {
    expConvTest_ = Teuchos::rcp( new StatusTestResNorm_t( convtol_ ) );
    expConvTest_->defineResForm( StatusTestResNorm_t::Explicit, Belos::TwoNorm );
    expConvTest_->defineScaleForm( convertStringToScaleType(expResScale_), Belos::TwoNorm );
    expConvTest_->setShowMaxResNormOnly( showMaxResNormOnly_ );
  }

  if (convTest_ == Teuchos::null || newImpResTest || newExpResTest)
    convTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::SEQ, impConvTest_, expConvTest_ ) );

  if (sTest_ == Teuchos::null || newImpResTest || newExpResTest)
    sTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::OR, maxIterTest_, convTest_ ) );
  
  if (outputTest_ == Teuchos::null || newImpResTest || newExpResTest) {
    if (outputFreq_ > 0) {
      outputTest_ = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer_, 
									  sTest_, 
									  outputFreq_, 
									  Passed+Failed+Undefined ) ); 
    }
    else {
      outputTest_ = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer_, 
									  sTest_, 1 ) );
    }
  }

  // Create orthogonalization manager if we need to.
  if (ortho_ == Teuchos::null) {
    if (orthoType_=="DGKS") {
      if (orthoKappa_ <= 0) {
	ortho_ = Teuchos::rcp( new DGKSOrthoManager<ScalarType,MV,OP>( label_ ) );
      }
      else {
	ortho_ = Teuchos::rcp( new DGKSOrthoManager<ScalarType,MV,OP>( label_ ) );
	Teuchos::rcp_dynamic_cast<DGKSOrthoManager<ScalarType,MV,OP> >(ortho_)->setDepTol( orthoKappa_ );
      }
    }
    else if (orthoType_=="ICGS") {
      ortho_ = Teuchos::rcp( new ICGSOrthoManager<ScalarType,MV,OP>( label_ ) );
    } 
    else {
      TEST_FOR_EXCEPTION(orthoType_!="ICGS"&&orthoType_!="DGKS",std::logic_error,
			 "Belos::BlockGmresSolMgr(): Invalid orthogonalization type.");
    }  
  }

  // Create the timer if we need to.
  if (timerSolve_ == Teuchos::null) {
    string solveLabel = label_ + ": BlockGmresSolMgr total solve time";
    timerSolve_ = Teuchos::TimeMonitor::getNewTimer(solveLabel);
  }

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}

    
template<class ScalarType, class MV, class OP>
void BlockGmresSolMgr<ScalarType,MV,OP>::setDefaultParams()
{
  defaultParams_ = Teuchos::rcp( new Teuchos::ParameterList() );
  
  // Set all the valid parameters and their default values.
  defaultParams_->set("Convergence Tolerance", convtol_default_);
  defaultParams_->set("Maximum Restarts", maxRestarts_default_);
  defaultParams_->set("Maximum Iterations", maxIters_default_);
  defaultParams_->set("Num Blocks", numBlocks_default_);
  defaultParams_->set("Block Size", blockSize_default_);
  defaultParams_->set("Adaptive Block Size", adaptiveBlockSize_default_);
  defaultParams_->set("Verbosity", verbosity_default_);
  defaultParams_->set("Output Frequency", outputFreq_default_);  
  defaultParams_->set("Output Stream", outputStream_default_);
  defaultParams_->set("Show Maximum Residual Norm Only", showMaxResNormOnly_default_);
  defaultParams_->set("Implicit Residual Scaling", impResScale_default_);
  defaultParams_->set("Explicit Residual Scaling", expResScale_default_);
  defaultParams_->set("Timer Label", label_default_);
  //  defaultParams_->set("Restart Timers", restartTimers_);
  defaultParams_->set("Orthogonalization", orthoType_default_);
  defaultParams_->set("Orthogonalization Constant",orthoKappa_default_);
}

  
// solve()
template<class ScalarType, class MV, class OP>
ReturnType BlockGmresSolMgr<ScalarType,MV,OP>::solve() {

  // Set the current parameters if they were not set before.
  // NOTE:  This may occur if the user generated the solver manager with the default constructor and 
  // then didn't set any parameters using setParameters().
  if (!isSet_) { setParameters( params_ ); }

  Teuchos::BLAS<int,ScalarType> blas;
  Teuchos::LAPACK<int,ScalarType> lapack;
  
  TEST_FOR_EXCEPTION(problem_ == Teuchos::null,BlockGmresSolMgrLinearProblemFailure,
                     "Belos::BlockGmresSolMgr::solve(): Linear problem is not a valid object.");

  TEST_FOR_EXCEPTION(!problem_->isProblemSet(),BlockGmresSolMgrLinearProblemFailure,
                     "Belos::BlockGmresSolMgr::solve(): Linear problem is not ready, setProblem() has not been called.");

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

  int dim = MVT::GetVecLength( *(problem_->getRHS()) );  
  if (blockSize_*numBlocks_ > dim) {
    int tmpNumBlocks = 0;
    if (blockSize_ == 1)
      tmpNumBlocks = dim / blockSize_;  // Allow for a good breakdown.
    else
      tmpNumBlocks = ( dim - blockSize_) / blockSize_;  // Allow for restarting.
    printer_->stream(Warnings) << 
      "Warning! Requested Krylov subspace dimension is larger that operator dimension!" << endl <<
      " The maximum number of blocks allowed for the Krylov subspace will be adjusted to " << tmpNumBlocks << endl;
    plist.set("Num Blocks",tmpNumBlocks);
  } 
  else 
    plist.set("Num Blocks",numBlocks_);
  
  // Reset the status test.  
  outputTest_->reset();

  // Assume convergence is achieved, then let any failed convergence set this to false.
  bool isConverged = true;	

  //////////////////////////////////////////////////////////////////////////////////////
  // BlockGmres solver

  Teuchos::RefCountPtr<BlockGmresIter<ScalarType,MV,OP> > block_gmres_iter
    = Teuchos::rcp( new BlockGmresIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,ortho_,plist) );
  
  // Enter solve() iterations
  {
    Teuchos::TimeMonitor slvtimer(*timerSolve_);

    while ( numRHS2Solve > 0 ) {

      // Set the current number of blocks and blocksize with the Gmres iteration.
      if (blockSize_*numBlocks_ > dim) {
        int tmpNumBlocks = 0;
        if (blockSize_ == 1)
          tmpNumBlocks = dim / blockSize_;  // Allow for a good breakdown.
        else
          tmpNumBlocks = ( dim - blockSize_) / blockSize_;  // Allow for restarting.
        block_gmres_iter->setSize( blockSize_, tmpNumBlocks );
      }
      else
        block_gmres_iter->setSize( blockSize_, numBlocks_ );

      // Reset the number of iterations.
      block_gmres_iter->resetNumIters();

      // Reset the number of calls that the status test output knows about.
      outputTest_->resetNumCalls();

      // Create the first block in the current Krylov basis.
      Teuchos::RefCountPtr<MV> V_0 = MVT::Clone( *(problem_->getRHS()), blockSize_ );
      problem_->computeCurrResVec( &*V_0 );

      // Get a matrix to hold the orthonormalization coefficients.
      Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > Z_0 = 
        rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(blockSize_, blockSize_) );
      
      // Orthonormalize the new V_0
      int rank = ortho_->normalize( *V_0, Z_0 );
      TEST_FOR_EXCEPTION(rank != blockSize_,BlockGmresSolMgrOrthoFailure,
			 "Belos::BlockGmresSolMgr::solve(): Failed to compute initial block of orthonormal vectors.");
     
      // Set the new state and initialize the solver.
      BlockGmresIterState<ScalarType,MV> newstate;
      newstate.V = V_0;
      newstate.Z = Z_0;
      newstate.curDim = 0;
      block_gmres_iter->initialize(newstate);
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
	    // we have convergence
	    break;  // break from while(1){block_gmres_iter->iterate()}
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
	    BlockGmresIterState<ScalarType,MV> oldState = block_gmres_iter->getState();
	    
	    // Compute the restart vector.
	    // Get a view of the current Krylov basis.
	    Teuchos::RefCountPtr<MV> V_0  = MVT::Clone( *(oldState.V), blockSize_ );
	    problem_->computeCurrResVec( &*V_0 );

	    // Get a view of the first block of the Krylov basis.
            Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > Z_0 = 
              rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(blockSize_, blockSize_) );
	    
	    // Orthonormalize the new V_0
	    int rank = ortho_->normalize( *V_0, Z_0 );
	    TEST_FOR_EXCEPTION(rank != blockSize_,BlockGmresSolMgrOrthoFailure,
			       "Belos::BlockGmresSolMgr::solve(): Failed to compute initial block of orthonormal vectors after restart.");

	    // Set the new state and initialize the solver.
	    BlockGmresIterState<ScalarType,MV> newstate;
	    newstate.V = V_0;
	    newstate.Z = Z_0;
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
			       "Belos::BlockGmresSolMgr::solve(): Invalid return from BlockGmresIter::iterate().");
	  }
	}
        catch (BlockGmresIterOrthoFailure e) {
	  // If the block size is not one, it's not considered a lucky breakdown.
	  if (blockSize_ != 1) {
	    printer_->stream(Errors) << "Error! Caught exception in BlockGmresIter::iterate() at iteration " 
				    << block_gmres_iter->getNumIters() << endl 
				    << e.what() << endl;
	    if (convTest_->getStatus() != Passed)
	      isConverged = false;
            break;
          } 
          else {
            // If the block size is one, try to recover the most recent least-squares solution
  	    block_gmres_iter->updateLSQR( block_gmres_iter->getCurSubspaceDim() );

	    // Check to see if the most recent least-squares solution yielded convergence.
	    sTest_->checkStatus( &*block_gmres_iter );
	    if (convTest_->getStatus() != Passed)
	      isConverged = false;
	    break;
          }
        }
        catch (std::exception e) {
	  printer_->stream(Errors) << "Error! Caught exception in BlockGmresIter::iterate() at iteration " 
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

  // print final summary
  sTest_->print( printer_->stream(FinalSummary) );
 
  // print timing information
  Teuchos::TimeMonitor::summarize( printer_->stream(TimingDetails) );
  
  if (!isConverged) {
    return Unconverged; // return from BlockGmresSolMgr::solve() 
  }
  return Converged; // return from BlockGmresSolMgr::solve() 
}

//  This method requires the solver manager to return a string that describes itself.
template<class ScalarType, class MV, class OP>
std::string BlockGmresSolMgr<ScalarType,MV,OP>::description() const
{
  std::ostringstream oss;
  oss << "Belos::BlockGmresSolMgr<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
  oss << "{";
  oss << "Ortho Type='"<<orthoType_<<"\', Block Size=" <<blockSize_;
  oss << ", Num Blocks=" <<numBlocks_<< ", Max Restarts=" << maxRestarts_;
  oss << "}";
  return oss.str();
}
  
} // end Belos namespace

#endif /* BELOS_BLOCK_GMRES_SOLMGR_HPP */
