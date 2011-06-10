//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

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
#include "BelosBlockCGIter.hpp"
#include "BelosDGKSOrthoManager.hpp"
#include "BelosICGSOrthoManager.hpp"
#include "BelosIMGSOrthoManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestGenResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutputFactory.hpp"
#include "BelosOutputManager.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

/** \example BlockCG/BlockCGEpetraExFile.cpp
    This is an example of how to use the Belos::BlockCGSolMgr solver manager.
*/
/** \example BlockCG/BlockPrecCGEpetraExFile.cpp
    This is an example of how to use the Belos::BlockCGSolMgr solver manager with an Ifpack preconditioner.
*/

/*! \class Belos::BlockCGSolMgr
 *
 *  \brief The Belos::BlockCGSolMgr provides a powerful and fully-featured solver manager over the CG and BlockCG linear solver.

 \ingroup belos_solver_framework

 \author Heidi Thornquist, Chris Baker, and Teri Barth
 */

namespace Belos {
  
  //! @name BlockCGSolMgr Exceptions
  //@{
  
  /** \brief BlockCGSolMgrLinearProblemFailure is thrown when the linear problem is
   * not setup (i.e. setProblem() was not called) when solve() is called.
   *
   * This std::exception is thrown from the BlockCGSolMgr::solve() method.
   *
   */
  class BlockCGSolMgrLinearProblemFailure : public BelosError {public:
    BlockCGSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  /** \brief BlockCGSolMgrOrthoFailure is thrown when the orthogonalization manager is
   * unable to generate orthonormal columns from the initial basis vectors.
   *
   * This std::exception is thrown from the BlockCGSolMgr::solve() method.
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

    /*! \brief Empty constructor for BlockCGSolMgr.
     * This constructor takes no arguments and sets the default values for the solver.
     * The linear problem must be passed in using setProblem() before solve() is called on this object.
     * The solver values can be changed using setParameters().
     */
     BlockCGSolMgr();

    /*! \brief Basic constructor for BlockCGSolMgr.
     *
     * This constructor accepts the LinearProblem to be solved in addition
     * to a parameter list of options for the solver manager. These options include the following:
     *   - "Block Size" - an \c int specifying the block size to be used by the underlying block 
     *                    conjugate-gradient solver. Default: 1
     *   - "Adaptive Block Size" - a \c bool specifying whether the block size can be modified 
     *                             throughout the solve. Default: true
     *   - "Maximum Iterations" - an \c int specifying the maximum number of iterations the 
     *                            underlying solver is allowed to perform. Default: 1000
     *   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms 
     *                               must reach to decide convergence. Default: 1e-8.
     *   - "Orthogonalization" - a \c std::string specifying the desired orthogonalization:  
     *                           DGKS ,ICGS, and IMGS. Default: "DGKS"
     *   - "Orthogonalization Constant" - a \c MagnitudeType used by DGKS orthogonalization to 
     *                                    determine whether another step of classical Gram-Schmidt 
     *                                    is necessary.  Default: -1 (use DGKS default)
     *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Belos::Errors
     *   - "Output Style" - a OutputType specifying the style of output. Default: Belos::General
     *   - "Output Stream" - a reference-counted pointer to the output stream where all
     *                       solver output is sent.  Default: Teuchos::rcp(&std::cout,false)
     *   - "Output Frequency" - an \c int specifying how often convergence information should be 
     *                          outputted.  Default: -1 (never)
     *   - "Show Maximum Residual Norm Only" - a \c bool specifying whether that only the maximum
     *                                         relative residual norm is printed if convergence 
     *                                         information is printed. Default: false
     *   - "Timer Label" - a \c std::string to use as a prefix for the timer labels.  Default: "Belos"
     */
    BlockCGSolMgr( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
		   const Teuchos::RCP<Teuchos::ParameterList> &pl );
    
    //! Destructor.
    virtual ~BlockCGSolMgr() {};
    //@}
    
    //! @name Accessor methods
    //@{ 
    
    const LinearProblem<ScalarType,MV,OP>& getProblem() const {
      return *problem_;
    }

    /*! \brief Get a parameter list containing the valid parameters for this object.
     */
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    
    /*! \brief Get a parameter list containing the current parameters for this object.
     */
    Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const { return params_; }
    
    /*! \brief Return the timers for this object. 
     *
     * The timers are ordered as follows:
     *   - time spent in solve() routine
     */
    Teuchos::Array<Teuchos::RCP<Teuchos::Time> > getTimers() const {
      return tuple(timerSolve_);
    }

    //! Get the iteration count for the most recent call to \c solve().
    int getNumIters() const {
      return numIters_;
    }

    /*! \brief Return whether a loss of accuracy was detected by this solver during the most current solve.
     */
    bool isLOADetected() const { return false; }
 
    //@}
    
    //! @name Set methods
    //@{
   
    //! Set the linear problem that needs to be solved. 
    void setProblem( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem ) { problem_ = problem; }
   
    //! Set the parameters the solver manager should use to solve the linear problem. 
    void setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params );
    
    //@}
   
    //! @name Reset methods
    //@{
    /*! \brief Performs a reset of the solver manager specified by the \c ResetType.  This informs the
     *  solver manager that the solver should prepare for the next call to solve by resetting certain elements
     *  of the iterative solver strategy.
     */
    void reset( const ResetType type ) { if ((type & Belos::Problem) && !Teuchos::is_null(problem_)) problem_->setProblem(); }
    //@}
 
    //! @name Solver application methods
    //@{ 
    
    /*! \brief This method performs possibly repeated calls to the underlying linear solver's 
     *         iterate() routine until the problem has been solved (as decided by the solver manager) 
     *         or the solver manager decides to quit.
     *
     * This method calls BlockCGIter::iterate() or CGIter::iterate(), which will return either because a 
     * specially constructed status test evaluates to ::Passed or an std::exception is thrown.
     *
     * A return from BlockCGIter::iterate() signifies one of the following scenarios:
     *    - the maximum number of iterations has been exceeded. In this scenario, the current solutions 
     *      to the linear system will be placed in the linear problem and return ::Unconverged.
     *    - global convergence has been met. In this case, the current solutions to the linear system 
     *      will be placed in the linear problem and the solver manager will return ::Converged
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
    Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;
    
    // Output manager.
    Teuchos::RCP<OutputManager<ScalarType> > printer_;
    Teuchos::RCP<std::ostream> outputStream_;

    // Status test.
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > sTest_;
    Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > maxIterTest_;
    Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> > convTest_;
    Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputTest_;

    // Orthogonalization manager.
    Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > ortho_; 
    
    // Current parameter list.
    Teuchos::RCP<Teuchos::ParameterList> params_;
    
    // Default solver values.
    static const MagnitudeType convtol_default_;
    static const MagnitudeType orthoKappa_default_;
    static const int maxIters_default_;
    static const bool adaptiveBlockSize_default_;
    static const bool showMaxResNormOnly_default_;
    static const int blockSize_default_;
    static const int verbosity_default_;
    static const int outputStyle_default_;
    static const int outputFreq_default_;
    static const std::string label_default_;
    static const std::string orthoType_default_;
    static const Teuchos::RCP<std::ostream> outputStream_default_;

    // Current solver values.
    MagnitudeType convtol_, orthoKappa_;
    int maxIters_, numIters_;
    int blockSize_, verbosity_, outputStyle_, outputFreq_;
    bool adaptiveBlockSize_, showMaxResNormOnly_;
    std::string orthoType_; 
    
    // Timers.
    std::string label_;
    Teuchos::RCP<Teuchos::Time> timerSolve_;

    // Internal state variables.
    bool isSet_;
  };


// Default solver values.
template<class ScalarType, class MV, class OP>
const typename BlockCGSolMgr<ScalarType,MV,OP>::MagnitudeType BlockCGSolMgr<ScalarType,MV,OP>::convtol_default_ = 1e-8;

template<class ScalarType, class MV, class OP>
const typename BlockCGSolMgr<ScalarType,MV,OP>::MagnitudeType BlockCGSolMgr<ScalarType,MV,OP>::orthoKappa_default_ = -1.0;

template<class ScalarType, class MV, class OP>
const int BlockCGSolMgr<ScalarType,MV,OP>::maxIters_default_ = 1000;

template<class ScalarType, class MV, class OP>
const bool BlockCGSolMgr<ScalarType,MV,OP>::adaptiveBlockSize_default_ = true;

template<class ScalarType, class MV, class OP>
const bool BlockCGSolMgr<ScalarType,MV,OP>::showMaxResNormOnly_default_ = false;

template<class ScalarType, class MV, class OP>
const int BlockCGSolMgr<ScalarType,MV,OP>::blockSize_default_ = 1;

template<class ScalarType, class MV, class OP>
const int BlockCGSolMgr<ScalarType,MV,OP>::verbosity_default_ = Belos::Errors;

template<class ScalarType, class MV, class OP>
const int BlockCGSolMgr<ScalarType,MV,OP>::outputStyle_default_ = Belos::General;

template<class ScalarType, class MV, class OP>
const int BlockCGSolMgr<ScalarType,MV,OP>::outputFreq_default_ = -1;

template<class ScalarType, class MV, class OP>
const std::string BlockCGSolMgr<ScalarType,MV,OP>::label_default_ = "Belos";

template<class ScalarType, class MV, class OP>
const std::string BlockCGSolMgr<ScalarType,MV,OP>::orthoType_default_ = "DGKS";

template<class ScalarType, class MV, class OP>
const Teuchos::RCP<std::ostream> BlockCGSolMgr<ScalarType,MV,OP>::outputStream_default_ = Teuchos::rcp(&std::cout,false);


// Empty Constructor
template<class ScalarType, class MV, class OP>
BlockCGSolMgr<ScalarType,MV,OP>::BlockCGSolMgr() :
  outputStream_(outputStream_default_),
  convtol_(convtol_default_),
  orthoKappa_(orthoKappa_default_),
  maxIters_(maxIters_default_),
  blockSize_(blockSize_default_),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  adaptiveBlockSize_(adaptiveBlockSize_default_),
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  orthoType_(orthoType_default_),
  label_(label_default_),
  isSet_(false)
{}


// Basic Constructor
template<class ScalarType, class MV, class OP>
BlockCGSolMgr<ScalarType,MV,OP>::BlockCGSolMgr( 
						     const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
						     const Teuchos::RCP<Teuchos::ParameterList> &pl ) : 
  problem_(problem),
  outputStream_(outputStream_default_),
  convtol_(convtol_default_),
  orthoKappa_(orthoKappa_default_),
  maxIters_(maxIters_default_),
  blockSize_(blockSize_default_),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  adaptiveBlockSize_(adaptiveBlockSize_default_),
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  orthoType_(orthoType_default_),
  label_(label_default_),
  isSet_(false)
{
  TEST_FOR_EXCEPTION(problem_ == Teuchos::null, std::invalid_argument, "Problem not given to solver manager.");

  // If the parameter list pointer is null, then set the current parameters to the default parameter list.
  if ( !is_null(pl) ) {
    setParameters( pl );  
  }
}

template<class ScalarType, class MV, class OP>
void BlockCGSolMgr<ScalarType,MV,OP>::setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params )
{
  // Create the internal parameter list if ones doesn't already exist.
  if (params_ == Teuchos::null) {
    params_ = Teuchos::rcp( new Teuchos::ParameterList(*getValidParameters()) );
  }
  else {
    params->validateParameters(*getValidParameters());
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
		       "Belos::BlockCGSolMgr: \"Block Size\" must be strictly positive.");

    // Update parameter in our list.
    params_->set("Block Size", blockSize_);
  }

  // Check if the blocksize should be adaptive
  if (params->isParameter("Adaptive Block Size")) {
    adaptiveBlockSize_ = params->get("Adaptive Block Size",adaptiveBlockSize_default_);
    
    // Update parameter in our list.
    params_->set("Adaptive Block Size", adaptiveBlockSize_);
  }

  // Check to see if the timer label changed.
  if (params->isParameter("Timer Label")) {
    std::string tempLabel = params->get("Timer Label", label_default_);

    // Update parameter in our list and solver timer
    if (tempLabel != label_) {
      label_ = tempLabel;
      params_->set("Timer Label", label_);
      std::string solveLabel = label_ + ": BlockCGSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerSolve_ = Teuchos::TimeMonitor::getNewTimer(solveLabel);
#endif
    }
  }

  // Check if the orthogonalization changed.
  if (params->isParameter("Orthogonalization")) {
    std::string tempOrthoType = params->get("Orthogonalization",orthoType_default_);
    TEST_FOR_EXCEPTION( tempOrthoType != "DGKS" && tempOrthoType != "ICGS" && tempOrthoType != "IMGS", 
                        std::invalid_argument,
			"Belos::BlockCGSolMgr: \"Orthogonalization\" must be either \"DGKS\", \"ICGS\", or \"IMGS\".");
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
      else if (orthoType_=="IMGS") {
	ortho_ = Teuchos::rcp( new IMGSOrthoManager<ScalarType,MV,OP>( label_ ) );
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

  // Check for a change in output style
  if (params->isParameter("Output Style")) {
    if (Teuchos::isParameterType<int>(*params,"Output Style")) {
      outputStyle_ = params->get("Output Style", outputStyle_default_);
    } else {
      outputStyle_ = (int)Teuchos::getParameter<Belos::OutputType>(*params,"Output Style");
    }

    // Update parameter in our list.
    params_->set("Output Style", outputStyle_);
    outputTest_ = Teuchos::null;
  }

  // output stream
  if (params->isParameter("Output Stream")) {
    outputStream_ = Teuchos::getParameter<Teuchos::RCP<std::ostream> >(*params,"Output Stream");

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
  typedef Belos::StatusTestGenResNorm<ScalarType,MV,OP>  StatusTestResNorm_t;

  // Check for convergence tolerance
  if (params->isParameter("Convergence Tolerance")) {
    convtol_ = params->get("Convergence Tolerance",convtol_default_);

    // Update parameter in our list and residual tests.
    params_->set("Convergence Tolerance", convtol_);
    if (convTest_ != Teuchos::null)
      convTest_->setTolerance( convtol_ );
  }
  
  if (params->isParameter("Show Maximum Residual Norm Only")) {
    showMaxResNormOnly_ = Teuchos::getParameter<bool>(*params,"Show Maximum Residual Norm Only");

    // Update parameter in our list and residual tests
    params_->set("Show Maximum Residual Norm Only", showMaxResNormOnly_);
    if (convTest_ != Teuchos::null)
      convTest_->setShowMaxResNormOnly( showMaxResNormOnly_ );
  }

  // Create status tests if we need to.

  // Basic test checks maximum iterations and native residual.
  if (maxIterTest_ == Teuchos::null)
    maxIterTest_ = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>( maxIters_ ) );
  
  // Implicit residual test, using the native residual to determine if convergence was achieved.
  if (convTest_ == Teuchos::null)
    convTest_ = Teuchos::rcp( new StatusTestResNorm_t( convtol_, 1 ) );
  
  if (sTest_ == Teuchos::null)
    sTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::OR, maxIterTest_, convTest_ ) );
  
  if (outputTest_ == Teuchos::null) {

    // Create the status test output class.
    // This class manages and formats the output from the status test.
    StatusTestOutputFactory<ScalarType,MV,OP> stoFactory( outputStyle_ );
    outputTest_ = stoFactory.create( printer_, sTest_, outputFreq_, Passed+Failed+Undefined );

    // Set the solver string for the output test
    std::string solverDesc = " Block CG ";
    outputTest_->setSolverDesc( solverDesc );

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
    else if (orthoType_=="IMGS") {
      ortho_ = Teuchos::rcp( new IMGSOrthoManager<ScalarType,MV,OP>( label_ ) );
    } 
    else {
      TEST_FOR_EXCEPTION(orthoType_!="ICGS"&&orthoType_!="DGKS"&&orthoType_!="IMGS",std::logic_error,
			 "Belos::BlockCGSolMgr(): Invalid orthogonalization type.");
    }  
  }

  // Create the timer if we need to.
  if (timerSolve_ == Teuchos::null) {
    std::string solveLabel = label_ + ": BlockCGSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    timerSolve_ = Teuchos::TimeMonitor::getNewTimer(solveLabel);
#endif
  }

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}

    
template<class ScalarType, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList> 
BlockCGSolMgr<ScalarType,MV,OP>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  
  // Set all the valid parameters and their default values.
  if(is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Convergence Tolerance", convtol_default_,
      "The relative residual tolerance that needs to be achieved by the\n"
      "iterative solver in order for the linear system to be declared converged.");
    pl->set("Maximum Iterations", maxIters_default_,
      "The maximum number of block iterations allowed for each\n"
      "set of RHS solved.");
    pl->set("Block Size", blockSize_default_,
      "The number of vectors in each block.");
    pl->set("Adaptive Block Size", adaptiveBlockSize_default_,
      "Whether the solver manager should adapt to the block size\n"
      "based on the number of RHS to solve.");
    pl->set("Verbosity", verbosity_default_,
      "What type(s) of solver information should be outputted\n"
      "to the output stream.");
    pl->set("Output Style", outputStyle_default_,
      "What style is used for the solver information outputted\n"
      "to the output stream.");
    pl->set("Output Frequency", outputFreq_default_,
      "How often convergence information should be outputted\n"
      "to the output stream.");  
    pl->set("Output Stream", outputStream_default_,
      "A reference-counted pointer to the output stream where all\n"
      "solver output is sent.");
    pl->set("Show Maximum Residual Norm Only", showMaxResNormOnly_default_,
      "When convergence information is printed, only show the maximum\n"
      "relative residual norm when the block size is greater than one.");
    pl->set("Timer Label", label_default_,
      "The string to use as a prefix for the timer labels.");
    //  pl->set("Restart Timers", restartTimers_);
    pl->set("Orthogonalization", orthoType_default_,
      "The type of orthogonalization to use: DGKS, ICGS, or IMGS.");
    pl->set("Orthogonalization Constant",orthoKappa_default_,
      "The constant used by DGKS orthogonalization to determine\n"
      "whether another step of classical Gram-Schmidt is necessary.");
    validPL = pl;
  }
  return validPL;
}

  
// solve()
template<class ScalarType, class MV, class OP>
ReturnType BlockCGSolMgr<ScalarType,MV,OP>::solve() {

  // Set the current parameters if they were not set before.
  // NOTE:  This may occur if the user generated the solver manager with the default constructor and 
  // then didn't set any parameters using setParameters().
  if (!isSet_) {
    setParameters(Teuchos::parameterList(*getValidParameters()));
  }

  Teuchos::BLAS<int,ScalarType> blas;
  Teuchos::LAPACK<int,ScalarType> lapack;
  
  TEST_FOR_EXCEPTION(!problem_->isProblemSet(),BlockCGSolMgrLinearProblemFailure,
                     "Belos::BlockCGSolMgr::solve(): Linear problem is not ready, setProblem() has not been called.");

  // Create indices for the linear systems to be solved.
  int startPtr = 0;
  int numRHS2Solve = MVT::GetNumberVecs( *(problem_->getRHS()) );
  int numCurrRHS = ( numRHS2Solve < blockSize_) ? numRHS2Solve : blockSize_;

  std::vector<int> currIdx, currIdx2;
  //  If an adaptive block size is allowed then only the linear systems that need to be solved are solved.
  //  Otherwise, the index set is generated that informs the linear problem that some linear systems are augmented.
  if ( adaptiveBlockSize_ ) {
    blockSize_ = numCurrRHS;
    currIdx.resize( numCurrRHS  );
    currIdx2.resize( numCurrRHS  );
    for (int i=0; i<numCurrRHS; ++i) 
      { currIdx[i] = startPtr+i; currIdx2[i]=i; }
    
  }
  else {
    currIdx.resize( blockSize_ );
    currIdx2.resize( blockSize_ );
    for (int i=0; i<numCurrRHS; ++i) 
      { currIdx[i] = startPtr+i; currIdx2[i]=i; }
    for (int i=numCurrRHS; i<blockSize_; ++i) 
      { currIdx[i] = -1; currIdx2[i] = i; }
  }

  // Inform the linear problem of the current linear system to solve.
  problem_->setLSIndex( currIdx );

  //////////////////////////////////////////////////////////////////////////////////////
  // Parameter list
  Teuchos::ParameterList plist;
  plist.set("Block Size",blockSize_);
  
  // Reset the status test.  
  outputTest_->reset();

  // Assume convergence is achieved, then let any failed convergence set this to false.
  bool isConverged = true;	

  //////////////////////////////////////////////////////////////////////////////////////
  // BlockCG solver

  Teuchos::RCP<CGIteration<ScalarType,MV,OP> > block_cg_iter;
  if (blockSize_ == 1)
    block_cg_iter = Teuchos::rcp( new CGIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,plist) );
  else
    block_cg_iter = Teuchos::rcp( new BlockCGIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,ortho_,plist) );
  

  // Enter solve() iterations
  {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*timerSolve_);
#endif

    while ( numRHS2Solve > 0 ) {
      //
      // Reset the active / converged vectors from this block
      std::vector<int> convRHSIdx;
      std::vector<int> currRHSIdx( currIdx );
      currRHSIdx.resize(numCurrRHS);

      // Reset the number of iterations.
      block_cg_iter->resetNumIters();

      // Reset the number of calls that the status test output knows about.
      outputTest_->resetNumCalls();

      // Get the current residual for this block of linear systems.
      Teuchos::RCP<MV> R_0 = MVT::CloneViewNonConst( *(Teuchos::rcp_const_cast<MV>(problem_->getInitResVec())), currIdx );

      // Set the new state and initialize the solver.
      CGIterationState<ScalarType,MV> newstate;
      newstate.R = R_0;
      block_cg_iter->initializeCG(newstate);

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
	    // We have convergence of at least one linear system.

	    // Figure out which linear systems converged.
	    std::vector<int> convIdx = Teuchos::rcp_dynamic_cast<StatusTestGenResNorm<ScalarType,MV,OP> >(convTest_)->convIndices();

	    // If the number of converged linear systems is equal to the
            // number of current linear systems, then we are done with this block.
	    if (convIdx.size() == currRHSIdx.size())
	      break;  // break from while(1){block_cg_iter->iterate()}

	    // Inform the linear problem that we are finished with this current linear system.
	    problem_->setCurrLS();

	    // Reset currRHSIdx to have the right-hand sides that are left to converge for this block.
	    int have = 0;
            std::vector<int> unconvIdx(currRHSIdx.size());
	    for (unsigned int i=0; i<currRHSIdx.size(); ++i) {
	      bool found = false;
	      for (unsigned int j=0; j<convIdx.size(); ++j) {
		if (currRHSIdx[i] == convIdx[j]) {
		  found = true;
		  break;
		}
	      }
	      if (!found) {
                currIdx2[have] = currIdx2[i];
		currRHSIdx[have++] = currRHSIdx[i];
	      }
              else {
              } 
	    }
	    currRHSIdx.resize(have);
	    currIdx2.resize(have);

	    // Set the remaining indices after deflation.
	    problem_->setLSIndex( currRHSIdx );

	    // Get the current residual vector.
	    std::vector<MagnitudeType> norms;
            R_0 = MVT::CloneCopy( *(block_cg_iter->getNativeResiduals(&norms)),currIdx2 );
	    for (int i=0; i<have; ++i) { currIdx2[i] = i; }

	    // Set the new blocksize for the solver.
	    block_cg_iter->setBlockSize( have );

	    // Set the new state and initialize the solver.
	    CGIterationState<ScalarType,MV> defstate;
	    defstate.R = R_0;
	    block_cg_iter->initializeCG(defstate);
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
			       "Belos::BlockCGSolMgr::solve(): Invalid return from CGIteration::iterate().");
	  }
	}
	catch (const std::exception &e) {
	  printer_->stream(Errors) << "Error! Caught std::exception in CGIteration::iterate() at iteration " 
				   << block_cg_iter->getNumIters() << std::endl 
				   << e.what() << std::endl;
	  throw;
	}
      }
      
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
	  currIdx2.resize( numCurrRHS  );
	  for (int i=0; i<numCurrRHS; ++i) 
	    { currIdx[i] = startPtr+i; currIdx2[i] = i; }	  
	}
	else {
	  currIdx.resize( blockSize_ );
	  currIdx2.resize( blockSize_ );
	  for (int i=0; i<numCurrRHS; ++i) 
	    { currIdx[i] = startPtr+i; currIdx2[i] = i; }
	  for (int i=numCurrRHS; i<blockSize_; ++i) 
	    { currIdx[i] = -1; currIdx2[i] = i; }
	}
	// Set the next indices.
	problem_->setLSIndex( currIdx );

	// Set the new blocksize for the solver.
	block_cg_iter->setBlockSize( blockSize_ );	
      }
      else {
        currIdx.resize( numRHS2Solve );
      }
      
    }// while ( numRHS2Solve > 0 )
    
  }

  // print final summary
  sTest_->print( printer_->stream(FinalSummary) );
 
  // print timing information
#ifdef BELOS_TEUCHOS_TIME_MONITOR
  // Calling summarize() can be expensive, so don't call unless the
  // user wants to print out timing details.  summarize() will do all
  // the work even if it's passed a "black hole" output stream.
  if (verbosity_ & TimingDetails)
    Teuchos::TimeMonitor::summarize( printer_->stream(TimingDetails) );
#endif
 
  // get iteration information for this solve
  numIters_ = maxIterTest_->getNumIters();
 
  if (!isConverged) {
    return Unconverged; // return from BlockCGSolMgr::solve() 
  }
  return Converged; // return from BlockCGSolMgr::solve() 
}

//  This method requires the solver manager to return a std::string that describes itself.
template<class ScalarType, class MV, class OP>
std::string BlockCGSolMgr<ScalarType,MV,OP>::description() const
{
  std::ostringstream oss;
  oss << "Belos::BlockCGSolMgr<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
  oss << "{";
  oss << "Ortho Type='"<<orthoType_<<"\', Block Size=" << blockSize_;
  oss << "}";
  return oss.str();
}
  
} // end Belos namespace

#endif /* BELOS_BLOCK_CG_SOLMGR_HPP */
