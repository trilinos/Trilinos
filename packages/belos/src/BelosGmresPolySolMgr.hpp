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

#ifndef BELOS_GMRES_POLY_SOLMGR_HPP
#define BELOS_GMRES_POLY_SOLMGR_HPP

/*! \file BelosGmresPolySolMgr.hpp
 *  \brief The Belos::GmresPolySolMgr provides a solver manager for the hybrid block GMRES linear solver.
 */

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosGmresPolyOp.hpp"
#include "BelosGmresIteration.hpp"
#include "BelosBlockGmresIter.hpp"
#include "BelosDGKSOrthoManager.hpp"
#include "BelosICGSOrthoManager.hpp"
#include "BelosIMGSOrthoManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestGenResNorm.hpp"
#include "BelosStatusTestImpResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutputFactory.hpp"
#include "BelosOutputManager.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_as.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

/** \example BlockGmres/BlockGmresPolyEpetraExFile.cpp
    This is an example of how to use the Belos::GmresPolySolMgr solver manager.
*/

/*! \class Belos::GmresPolySolMgr
 *
 *  \brief The Belos::GmresPolySolMgr provides a powerful and fully-featured solver manager over the hybrid block GMRES linear solver.

 \ingroup belos_solver_framework

 \author Heidi Thornquist
*/

namespace Belos {
  
//! @name GmresPolySolMgr Exceptions
//@{
  
/** \brief GmresPolySolMgrLinearProblemFailure is thrown when the linear problem is
 * not setup (i.e. setProblem() was not called) when solve() is called.
 *
 * This std::exception is thrown from the GmresPolySolMgr::solve() method.
 *
 */
class GmresPolySolMgrLinearProblemFailure : public BelosError {public:
  GmresPolySolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
/** \brief GmresPolySolMgrPolynomialFailure is thrown when their is a problem generating
 * the GMRES polynomial for this linear problem.
 *
 * This std::exception is thrown from the GmresPolySolMgr::solve() method.
 *
 */
class GmresPolySolMgrPolynomialFailure : public BelosError {public:
  GmresPolySolMgrPolynomialFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
/** \brief GmresPolySolMgrOrthoFailure is thrown when the orthogonalization manager is
 * unable to generate orthonormal columns from the initial basis vectors.
 *
 * This std::exception is thrown from the GmresPolySolMgr::solve() method.
 *
 */
class GmresPolySolMgrOrthoFailure : public BelosError {public:
  GmresPolySolMgrOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
template<class ScalarType, class MV, class OP>
class GmresPolySolMgr : public SolverManager<ScalarType,MV,OP> {
    
private:
  typedef MultiVecTraits<ScalarType,MV> MVT;
  typedef MultiVecTraitsExt<ScalarType,MV> MVText;
  typedef OperatorTraits<ScalarType,MV,OP> OPT;
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
  typedef Teuchos::ScalarTraits<MagnitudeType> MT;
    
public:
    
  //! @name Constructors/Destructor
  //@{ 
   
  /*! \brief Empty constructor for GmresPolySolMgr.
   * This constructor takes no arguments and sets the default values for the solver.
   * The linear problem must be passed in using setProblem() before solve() is called on this object.
   * The solver values can be changed using setParameters().
   */
  GmresPolySolMgr();
 
  /*! \brief Basic constructor for GmresPolySolMgr.
   *
   * This constructor accepts the LinearProblem to be solved in addition
   * to a parameter list of options for the solver manager. These options include the following:
   *   - "Block Size" - a \c int specifying the block size to be used by the underlying block GMRES solver. Default: 1
   *   - "Num Blocks" - a \c int specifying the number of blocks allocated for the Krylov basis. Default: 300
   *   - "Maximum Iterations" - a \c int specifying the maximum number of iterations the underlying solver is allowed to perform. Default: 1000
   *   - "Maximum Restarts" - a \c int specifying the maximum number of restarts the underlying solver is allowed to perform. Default: 20
   *   - "Orthogonalization" - a \c std::string specifying the desired orthogonalization:  DGKS, ICGS, and IMGS. Default: "DGKS"
   *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Belos::Errors
   *   - "Output Style" - a OutputType specifying the style of output. Default: Belos::General
   *   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: 1e-8
   */
  GmresPolySolMgr( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
    const Teuchos::RCP<Teuchos::ParameterList> &pl );
    
  //! Destructor.
  virtual ~GmresPolySolMgr() {};
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
    return Teuchos::tuple(timerSolve_, timerPoly_);
  }
  
  //! Get the iteration count for the most recent call to \c solve().
  int getNumIters() const {
    return numIters_;
  }
 
  /*! \brief Return whether a loss of accuracy was detected by this solver during the most current solve.
      \note This flag will be reset the next time solve() is called.
   */
  bool isLOADetected() const { return loaDetected_; }
 
  //@}
    
  //! @name Set methods
  //@{
  
  //! Set the linear problem that needs to be solved.  
  void setProblem( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem ) { problem_ = problem; isSTSet_ = false; }

  //! Set the parameters the solver manager should use to solve the linear problem.  
  void setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params );
    
  //@}
  
  //! @name Reset methods
  //@{
  /*! \brief Performs a reset of the solver manager specified by the \c ResetType.  This informs the
   *  solver manager that the solver should prepare for the next call to solve by resetting certain elements
   *  of the iterative solver strategy.
  */
  void reset( const ResetType type ) { 
    if ((type & Belos::Problem) && !Teuchos::is_null(problem_)) {
      problem_->setProblem(); 
      isPolyBuilt_ = false;  // Rebuild the GMRES polynomial
    }
  }
  //@}
  
  //! @name Solver application methods
  //@{ 
    
  /*! \brief This method performs possibly repeated calls to the underlying linear solver's iterate() routine
   * until the problem has been solved (as decided by the solver manager) or the solver manager decides to 
   * quit.
   *
   * This method calls BlockGmresIter::iterate(), which will return either because a specially constructed status test evaluates to 
   * ::Passed or an std::exception is thrown.
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
    
  /** \brief Method to return description of the hybrid block GMRES solver manager */
  std::string description() const;
    
  //@}
    
private:

  // Method for checking current status test against defined linear problem.
  bool checkStatusTest();

  // Method for generating GMRES polynomial.
  bool generatePoly();

  // Linear problem.
  Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;
    
  // Output manager.
  Teuchos::RCP<OutputManager<ScalarType> > printer_;
  Teuchos::RCP<std::ostream> outputStream_;

  // Status test.
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > sTest_;
  Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > maxIterTest_;
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > convTest_;
  Teuchos::RCP<StatusTestResNorm<ScalarType,MV,OP> > expConvTest_, impConvTest_;
  Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputTest_;

  // Orthogonalization manager.
  Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > ortho_; 
    
  // Current parameter list.
  Teuchos::RCP<Teuchos::ParameterList> params_;

  // Default solver values.
  static const MagnitudeType polytol_default_;
  static const MagnitudeType convtol_default_;
  static const MagnitudeType orthoKappa_default_;
  static const int maxDegree_default_;
  static const int maxRestarts_default_;
  static const int maxIters_default_;
  static const bool strictConvTol_default_;
  static const bool showMaxResNormOnly_default_;
  static const int blockSize_default_;
  static const int numBlocks_default_;
  static const int verbosity_default_;
  static const int outputStyle_default_;
  static const int outputFreq_default_;
  static const std::string impResScale_default_; 
  static const std::string expResScale_default_; 
  static const std::string label_default_;
  static const std::string orthoType_default_;
  static const Teuchos::RCP<std::ostream> outputStream_default_;

  // Current solver values.
  MagnitudeType polytol_, convtol_, orthoKappa_;
  int maxDegree_, maxRestarts_, maxIters_, numIters_;
  int blockSize_, numBlocks_, verbosity_, outputStyle_, outputFreq_;
  bool strictConvTol_, showMaxResNormOnly_;
  std::string orthoType_; 
  std::string impResScale_, expResScale_;

  // Polynomial storage
  int poly_dim_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int, ScalarType> > poly_H_, poly_y_;
  Teuchos::RCP<Teuchos::SerialDenseVector<int, ScalarType> > poly_r0_;
  Teuchos::RCP<Belos::GmresPolyOp<ScalarType, MV, OP> > poly_Op_;  
  
  // Timers.
  std::string label_;
  Teuchos::RCP<Teuchos::Time> timerSolve_, timerPoly_;

  // Internal state variables.
  bool isPolyBuilt_;	
  bool isSet_, isSTSet_, expResTest_;
  bool loaDetected_;
};


// Default solver values.
template<class ScalarType, class MV, class OP>
const typename GmresPolySolMgr<ScalarType,MV,OP>::MagnitudeType GmresPolySolMgr<ScalarType,MV,OP>::polytol_default_ = 1e-12;

template<class ScalarType, class MV, class OP>
const typename GmresPolySolMgr<ScalarType,MV,OP>::MagnitudeType GmresPolySolMgr<ScalarType,MV,OP>::convtol_default_ = 1e-8;

template<class ScalarType, class MV, class OP>
const typename GmresPolySolMgr<ScalarType,MV,OP>::MagnitudeType GmresPolySolMgr<ScalarType,MV,OP>::orthoKappa_default_ = -1.0;

template<class ScalarType, class MV, class OP>
const int GmresPolySolMgr<ScalarType,MV,OP>::maxDegree_default_ = 25;

template<class ScalarType, class MV, class OP>
const int GmresPolySolMgr<ScalarType,MV,OP>::maxRestarts_default_ = 20;

template<class ScalarType, class MV, class OP>
const int GmresPolySolMgr<ScalarType,MV,OP>::maxIters_default_ = 1000;

template<class ScalarType, class MV, class OP>
const bool GmresPolySolMgr<ScalarType,MV,OP>::strictConvTol_default_ = false;

template<class ScalarType, class MV, class OP>
const bool GmresPolySolMgr<ScalarType,MV,OP>::showMaxResNormOnly_default_ = false;

template<class ScalarType, class MV, class OP>
const int GmresPolySolMgr<ScalarType,MV,OP>::blockSize_default_ = 1;

template<class ScalarType, class MV, class OP>
const int GmresPolySolMgr<ScalarType,MV,OP>::numBlocks_default_ = 300;

template<class ScalarType, class MV, class OP>
const int GmresPolySolMgr<ScalarType,MV,OP>::verbosity_default_ = Belos::Errors;

template<class ScalarType, class MV, class OP>
const int GmresPolySolMgr<ScalarType,MV,OP>::outputStyle_default_ = Belos::General;

template<class ScalarType, class MV, class OP>
const int GmresPolySolMgr<ScalarType,MV,OP>::outputFreq_default_ = -1;

template<class ScalarType, class MV, class OP>
const std::string GmresPolySolMgr<ScalarType,MV,OP>::impResScale_default_ = "Norm of RHS";

template<class ScalarType, class MV, class OP>
const std::string GmresPolySolMgr<ScalarType,MV,OP>::expResScale_default_ = "Norm of RHS";

template<class ScalarType, class MV, class OP>
const std::string GmresPolySolMgr<ScalarType,MV,OP>::label_default_ = "Belos";

template<class ScalarType, class MV, class OP>
const std::string GmresPolySolMgr<ScalarType,MV,OP>::orthoType_default_ = "DGKS";

template<class ScalarType, class MV, class OP>
const Teuchos::RCP<std::ostream> GmresPolySolMgr<ScalarType,MV,OP>::outputStream_default_ = Teuchos::rcp(&std::cout,false);


// Empty Constructor
template<class ScalarType, class MV, class OP>
GmresPolySolMgr<ScalarType,MV,OP>::GmresPolySolMgr() :
  outputStream_(outputStream_default_),
  polytol_(polytol_default_),
  convtol_(convtol_default_),
  orthoKappa_(orthoKappa_default_),
  maxDegree_(maxDegree_default_),
  maxRestarts_(maxRestarts_default_),
  maxIters_(maxIters_default_),
  blockSize_(blockSize_default_),
  numBlocks_(numBlocks_default_),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  strictConvTol_(strictConvTol_default_),
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  orthoType_(orthoType_default_),
  impResScale_(impResScale_default_),
  expResScale_(expResScale_default_),
  label_(label_default_),
  isPolyBuilt_(false),
  isSet_(false),
  isSTSet_(false),
  expResTest_(false),
  loaDetected_(false)
{}


// Basic Constructor
template<class ScalarType, class MV, class OP>
GmresPolySolMgr<ScalarType,MV,OP>::GmresPolySolMgr( 
  const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
  const Teuchos::RCP<Teuchos::ParameterList> &pl ) : 
  problem_(problem),
  outputStream_(outputStream_default_),
  polytol_(polytol_default_),
  convtol_(convtol_default_),
  orthoKappa_(orthoKappa_default_),
  maxDegree_(maxDegree_default_),
  maxRestarts_(maxRestarts_default_),
  maxIters_(maxIters_default_),
  blockSize_(blockSize_default_),
  numBlocks_(numBlocks_default_),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  strictConvTol_(strictConvTol_default_), 
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  orthoType_(orthoType_default_),
  impResScale_(impResScale_default_),
  expResScale_(expResScale_default_),
  label_(label_default_),
  isPolyBuilt_(false),
  isSet_(false),
  isSTSet_(false),
  expResTest_(false),
  loaDetected_(false)
{

  TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null, std::invalid_argument, "Problem not given to solver manager.");

  // If the parameter list pointer is null, then set the current parameters to the default parameter list.
  if ( !is_null(pl) ) {
    setParameters( pl );  
  }

}


template<class ScalarType, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList>
GmresPolySolMgr<ScalarType,MV,OP>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Polynomial Tolerance", polytol_default_,
      "The relative residual tolerance that used to construct the GMRES polynomial.");
    pl->set("Maximum Degree", maxDegree_default_,
      "The maximum degree allowed for any GMRES polynomial.");
    pl->set("Convergence Tolerance", convtol_default_,
      "The relative residual tolerance that needs to be achieved by the\n"
      "iterative solver in order for the linear system to be declared converged." );
    pl->set("Maximum Restarts", maxRestarts_default_,
      "The maximum number of restarts allowed for each\n"
      "set of RHS solved.");
    pl->set("Maximum Iterations", maxIters_default_,
      "The maximum number of block iterations allowed for each\n"
      "set of RHS solved.");
    pl->set("Num Blocks", numBlocks_default_,
      "The maximum number of blocks allowed in the Krylov subspace\n"
      "for each set of RHS solved.");
    pl->set("Block Size", blockSize_default_,
      "The number of vectors in each block.  This number times the\n"
      "number of blocks is the total Krylov subspace dimension.");
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
    pl->set("Strict Convergence", strictConvTol_default_,
      "After polynomial is applied, whether solver should try to achieve\n"
      "the relative residual tolerance.");
    pl->set("Show Maximum Residual Norm Only", showMaxResNormOnly_default_,
      "When convergence information is printed, only show the maximum\n"
      "relative residual norm when the block size is greater than one.");
    pl->set("Implicit Residual Scaling", impResScale_default_,
      "The type of scaling used in the implicit residual convergence test.");
    pl->set("Explicit Residual Scaling", expResScale_default_,
      "The type of scaling used in the explicit residual convergence test.");
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


template<class ScalarType, class MV, class OP>
void GmresPolySolMgr<ScalarType,MV,OP>::setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params )
{

  // Create the internal parameter list if ones doesn't already exist.
  if (params_ == Teuchos::null) {
    params_ = Teuchos::rcp( new Teuchos::ParameterList(*getValidParameters()) );
  }
  else {
    params->validateParameters(*getValidParameters());
  }

  // Check for maximum polynomial degree
  if (params->isParameter("Maximum Degree")) {
    maxDegree_ = params->get("Maximum Degree",maxDegree_default_);

    // Update parameter in our list.
    params_->set("Maximum Degree", maxDegree_);
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
    TEUCHOS_TEST_FOR_EXCEPTION(blockSize_ <= 0, std::invalid_argument,
      "Belos::GmresPolySolMgr: \"Block Size\" must be strictly positive.");

    // Update parameter in our list.
    params_->set("Block Size", blockSize_);
  }

  // Check for the maximum number of blocks.
  if (params->isParameter("Num Blocks")) {
    numBlocks_ = params->get("Num Blocks",numBlocks_default_);
    TEUCHOS_TEST_FOR_EXCEPTION(numBlocks_ <= 0, std::invalid_argument,
      "Belos::GmresPolySolMgr: \"Num Blocks\" must be strictly positive.");

    // Update parameter in our list.
    params_->set("Num Blocks", numBlocks_);
  }

  // Check to see if the timer label changed.
  if (params->isParameter("Timer Label")) {
    std::string tempLabel = params->get("Timer Label", label_default_);

    // Update parameter in our list and solver timer
    if (tempLabel != label_) {
      label_ = tempLabel;
      params_->set("Timer Label", label_);
      std::string solveLabel = label_ + ": GmresPolySolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerSolve_ = Teuchos::TimeMonitor::getNewCounter(solveLabel);
#endif
      std::string polyLabel = label_ + ": GmresPolySolMgr polynomial creation time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerPoly_ = Teuchos::TimeMonitor::getNewCounter(polyLabel);
#endif
      if (ortho_ != Teuchos::null) {
        ortho_->setLabel( label_ );
      }
    }
  }

  // Check if the orthogonalization changed.
  if (params->isParameter("Orthogonalization")) {
    std::string tempOrthoType = params->get("Orthogonalization",orthoType_default_);
    TEUCHOS_TEST_FOR_EXCEPTION( tempOrthoType != "DGKS" && tempOrthoType != "ICGS" && tempOrthoType != "IMGS", 
                        std::invalid_argument,
			"Belos::GmresPolySolMgr: \"Orthogonalization\" must be either \"DGKS\", \"ICGS\", or \"IMGS\".");
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

    // Reconstruct the convergence test if the explicit residual test is not being used.
    params_->set("Output Style", outputStyle_);
    if (outputTest_ != Teuchos::null) {
      isSTSet_ = false;
    }
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

  // Check for polynomial convergence tolerance
  if (params->isParameter("Polynomial Tolerance")) {
    polytol_ = params->get("Polynomial Tolerance",polytol_default_);

    // Update parameter in our list and residual tests.
    params_->set("Polynomial Tolerance", polytol_);
  }

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
 
  // Check if user requires solver to reach convergence tolerance
  if (params->isParameter("Strict Convergence")) {
    strictConvTol_ = params->get("Strict Convergence",strictConvTol_default_);

    // Update parameter in our list and residual tests
    params_->set("Strict Convergence", strictConvTol_);
  }
 
  // Check for a change in scaling, if so we need to build new residual tests.
  if (params->isParameter("Implicit Residual Scaling")) {
    std::string tempImpResScale = Teuchos::getParameter<std::string>( *params, "Implicit Residual Scaling" );

    // Only update the scaling if it's different.
    if (impResScale_ != tempImpResScale) {
      Belos::ScaleType impResScaleType = convertStringToScaleType( tempImpResScale );
      impResScale_ = tempImpResScale;

      // Update parameter in our list and residual tests
      params_->set("Implicit Residual Scaling", impResScale_);
      if (impConvTest_ != Teuchos::null) {
        try { 
          impConvTest_->defineScaleForm( impResScaleType, Belos::TwoNorm );
        }
        catch (std::exception& e) { 
          // Make sure the convergence test gets constructed again.
          isSTSet_ = false;
        }
      }
    }      
  }
  
  if (params->isParameter("Explicit Residual Scaling")) {
    std::string tempExpResScale = Teuchos::getParameter<std::string>( *params, "Explicit Residual Scaling" );

    // Only update the scaling if it's different.
    if (expResScale_ != tempExpResScale) {
      Belos::ScaleType expResScaleType = convertStringToScaleType( tempExpResScale );
      expResScale_ = tempExpResScale;

      // Update parameter in our list and residual tests
      params_->set("Explicit Residual Scaling", expResScale_);
      if (expConvTest_ != Teuchos::null) {
        try { 
          expConvTest_->defineScaleForm( expResScaleType, Belos::TwoNorm );
        }
        catch (std::exception& e) {
          // Make sure the convergence test gets constructed again.
          isSTSet_ = false;
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
      TEUCHOS_TEST_FOR_EXCEPTION(orthoType_!="ICGS"&&orthoType_!="DGKS"&&orthoType_!="IMGS",std::logic_error,
        "Belos::GmresPolySolMgr(): Invalid orthogonalization type.");
    }  
  }

  // Create the timers if we need to.
  if (timerSolve_ == Teuchos::null) {
    std::string solveLabel = label_ + ": GmresPolySolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    timerSolve_ = Teuchos::TimeMonitor::getNewCounter(solveLabel);
#endif
  }
  
  if (timerPoly_ == Teuchos::null) {
    std::string polyLabel = label_ + ": GmresPolySolMgr polynomial creation time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    timerPoly_ = Teuchos::TimeMonitor::getNewCounter(polyLabel);
#endif
  }

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}

// Check the status test versus the defined linear problem
template<class ScalarType, class MV, class OP>
bool GmresPolySolMgr<ScalarType,MV,OP>::checkStatusTest() {

  typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
  typedef Belos::StatusTestGenResNorm<ScalarType,MV,OP>  StatusTestGenResNorm_t;
  typedef Belos::StatusTestImpResNorm<ScalarType,MV,OP>  StatusTestImpResNorm_t;

  // Basic test checks maximum iterations and native residual.
  maxIterTest_ = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>( maxIters_ ) );

  // If there is a left preconditioner, we create a combined status test that checks the implicit
  // and then explicit residual norm to see if we have convergence.
  if (!Teuchos::is_null(problem_->getLeftPrec())) {
    expResTest_ = true;
  }

  if (expResTest_) {
   
    // Implicit residual test, using the native residual to determine if convergence was achieved.
    Teuchos::RCP<StatusTestGenResNorm_t> tmpImpConvTest =
      Teuchos::rcp( new StatusTestGenResNorm_t( convtol_ ) );
    tmpImpConvTest->defineScaleForm( convertStringToScaleType(impResScale_), Belos::TwoNorm );
    tmpImpConvTest->setShowMaxResNormOnly( showMaxResNormOnly_ );
    impConvTest_ = tmpImpConvTest;

    // Explicit residual test once the native residual is below the tolerance
    Teuchos::RCP<StatusTestGenResNorm_t> tmpExpConvTest =
      Teuchos::rcp( new StatusTestGenResNorm_t( convtol_ ) );
    tmpExpConvTest->defineResForm( StatusTestGenResNorm_t::Explicit, Belos::TwoNorm );
    tmpExpConvTest->defineScaleForm( convertStringToScaleType(expResScale_), Belos::TwoNorm );
    tmpExpConvTest->setShowMaxResNormOnly( showMaxResNormOnly_ );
    expConvTest_ = tmpExpConvTest;

    // The convergence test is a combination of the "cheap" implicit test and explicit test.
    convTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::SEQ, impConvTest_, expConvTest_ ) );
  }
  else {

    // Implicit residual test, using the native residual to determine if convergence was achieved.
    // Use test that checks for loss of accuracy.
    Teuchos::RCP<StatusTestImpResNorm_t> tmpImpConvTest =
      Teuchos::rcp( new StatusTestImpResNorm_t( convtol_ ) );
    tmpImpConvTest->defineScaleForm( convertStringToScaleType(impResScale_), Belos::TwoNorm );
    tmpImpConvTest->setShowMaxResNormOnly( showMaxResNormOnly_ );
    impConvTest_ = tmpImpConvTest;

    // Set the explicit and total convergence test to this implicit test that checks for accuracy loss.
    expConvTest_ = impConvTest_;
    convTest_ = impConvTest_;
  }

  sTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::OR, maxIterTest_, convTest_ ) );

  // Create the status test output class.
  // This class manages and formats the output from the status test.
  StatusTestOutputFactory<ScalarType,MV,OP> stoFactory( outputStyle_ );
  outputTest_ = stoFactory.create( printer_, sTest_, outputFreq_, Passed+Failed+Undefined );

  // Set the solver string for the output test
  std::string solverDesc = " Gmres Polynomial ";
  outputTest_->setSolverDesc( solverDesc );


  // The status test is now set.
  isSTSet_ = true;

  return false;
}

template<class ScalarType, class MV, class OP>
bool GmresPolySolMgr<ScalarType,MV,OP>::generatePoly()
{
  // Create a copy of the linear problem that has a zero initial guess and random RHS.
  Teuchos::RCP<MV> newX  = MVT::Clone( *(problem_->getLHS()), 1 );
  Teuchos::RCP<MV> newB  = MVT::Clone( *(problem_->getRHS()), 1 );
  MVT::MvInit( *newX, SCT::zero() );
  MVT::MvRandom( *newB );
  Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > newProblem =
    Teuchos::rcp( new LinearProblem<ScalarType,MV,OP>( problem_->getOperator(), newX, newB ) );
  newProblem->setLeftPrec( problem_->getLeftPrec() );
  newProblem->setRightPrec( problem_->getRightPrec() );
  newProblem->setLabel("Belos GMRES Poly Generation");
  newProblem->setProblem();
  std::vector<int> idx(1,0);       // Must set the index to be the first vector (0)!
  newProblem->setLSIndex( idx );

  // Create a parameter list for the GMRES iteration.
  Teuchos::ParameterList polyList;

  // Tell the block solver that the block size is one.
  polyList.set("Num Blocks",maxDegree_);
  polyList.set("Block Size",1);
  polyList.set("Keep Hessenberg", true);

  // Create a simple status test that either reaches the relative residual tolerance or maximum polynomial size.
  Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > maxItrTst =
    Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>( maxDegree_ ) );

  // Implicit residual test, using the native residual to determine if convergence was achieved.
  Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> > convTst = 
    Teuchos::rcp( new StatusTestGenResNorm<ScalarType,MV,OP>( polytol_ ) );
  convTst->defineScaleForm( convertStringToScaleType(impResScale_), Belos::TwoNorm );

  // Convergence test that stops the iteration when either are satisfied.
  Teuchos::RCP<StatusTestCombo<ScalarType,MV,OP> > polyTest = 
    Teuchos::rcp( new StatusTestCombo<ScalarType,MV,OP>( StatusTestCombo<ScalarType,MV,OP>::OR, maxItrTst, convTst ) );

  // Create Gmres iteration object to perform one cycle of Gmres.
  Teuchos::RCP<BlockGmresIter<ScalarType,MV,OP> > gmres_iter;
  gmres_iter = Teuchos::rcp( new BlockGmresIter<ScalarType,MV,OP>(newProblem,printer_,polyTest,ortho_,polyList) );

  // Create the first block in the current Krylov basis (residual).
  Teuchos::RCP<MV> V_0 = MVT::Clone( *(newProblem->getRHS()), 1 );
  newProblem->computeCurrPrecResVec( &*V_0 );

  // Get a matrix to hold the orthonormalization coefficients.
  poly_r0_ = Teuchos::rcp( new Teuchos::SerialDenseVector<int,ScalarType>(1) );

  // Orthonormalize the new V_0
  int rank = ortho_->normalize( *V_0, poly_r0_ );
  TEUCHOS_TEST_FOR_EXCEPTION(rank != 1,GmresPolySolMgrOrthoFailure,
    "Belos::GmresPolySolMgr::generatePoly(): Failed to compute initial block of orthonormal vectors for polynomial generation.");

  // Set the new state and initialize the solver.
  GmresIterationState<ScalarType,MV> newstate;
  newstate.V = V_0;
  newstate.z = poly_r0_;
  newstate.curDim = 0;
  gmres_iter->initializeGmres(newstate);

  // Perform Gmres iteration
  bool polyConverged = false;
  try {
    gmres_iter->iterate();

    // Check convergence first
    if ( convTst->getStatus() == Passed ) {
      // we have convergence
      polyConverged = true;
    }
  }
  catch (GmresIterationOrthoFailure e) {
    // Try to recover the most recent least-squares solution
    gmres_iter->updateLSQR( gmres_iter->getCurSubspaceDim() );

    // Check to see if the most recent least-squares solution yielded convergence.
    polyTest->checkStatus( &*gmres_iter );
    if (convTst->getStatus() == Passed)
      polyConverged = true;
    }
    catch (std::exception e) {
    printer_->stream(Errors) << "Error! Caught exception in BlockGmresIter::iterate() at iteration "
                             << gmres_iter->getNumIters() << std::endl
                             << e.what() << std::endl;
    throw;
  }

  // Get the solution for this polynomial, use in comparison below
  Teuchos::RCP<MV> currX = gmres_iter->getCurrentUpdate();

  // Record polynomial info, get current GMRES state
  GmresIterationState<ScalarType,MV> gmresState = gmres_iter->getState();

  // If the polynomial has no dimension, the tolerance is too low, return false
  poly_dim_ = gmresState.curDim;
  if (poly_dim_ == 0) {
    return false;
  }
  //
  //  Make a view and then copy the RHS of the least squares problem.
  //
  poly_y_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::Copy, *gmresState.z, poly_dim_, 1 ) );
  poly_H_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( *gmresState.H ) );
  //
  // Solve the least squares problem.
  //
  const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  Teuchos::BLAS<int,ScalarType> blas;
  blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
             Teuchos::NON_UNIT_DIAG, poly_dim_, 1, one,
             gmresState.R->values(), gmresState.R->stride(), 
             poly_y_->values(), poly_y_->stride() );
  //
  // Generate the polynomial operator
  //
  poly_Op_ = Teuchos::rcp( 
               new Belos::GmresPolyOp<ScalarType,MV,OP>( problem_, poly_H_, poly_y_, poly_r0_ ) );
 
  return true;
}
  
// solve()
template<class ScalarType, class MV, class OP>
ReturnType GmresPolySolMgr<ScalarType,MV,OP>::solve() {

  // Set the current parameters if they were not set before.
  // NOTE:  This may occur if the user generated the solver manager with the default constructor and 
  // then didn't set any parameters using setParameters().
  if (!isSet_) {
    setParameters(Teuchos::parameterList(*getValidParameters()));
  }

  Teuchos::BLAS<int,ScalarType> blas;
  Teuchos::LAPACK<int,ScalarType> lapack;
  
  TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null,GmresPolySolMgrLinearProblemFailure,
    "Belos::GmresPolySolMgr::solve(): Linear problem is not a valid object.");

  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isProblemSet(),GmresPolySolMgrLinearProblemFailure,
    "Belos::GmresPolySolMgr::solve(): Linear problem is not ready, setProblem() has not been called.");

  if (!isSTSet_ || (!expResTest_ && !Teuchos::is_null(problem_->getLeftPrec())) ) {
    TEUCHOS_TEST_FOR_EXCEPTION( checkStatusTest(),GmresPolySolMgrLinearProblemFailure,
      "Belos::GmresPolySolMgr::solve(): Linear problem and requested status tests are incompatible.");
  }

  // If the GMRES polynomial has not been constructed for this matrix, preconditioner pair, generate it
  if (!isPolyBuilt_) {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*timerPoly_);
#endif
    isPolyBuilt_ = generatePoly();
    TEUCHOS_TEST_FOR_EXCEPTION( !isPolyBuilt_ && poly_dim_==0, GmresPolySolMgrPolynomialFailure,
      "Belos::GmresPolySolMgr::generatePoly(): Failed to generate a non-trivial polynomial, reduce polynomial tolerance.");
    TEUCHOS_TEST_FOR_EXCEPTION( !isPolyBuilt_, GmresPolySolMgrPolynomialFailure,
      "Belos::GmresPolySolMgr::generatePoly(): Failed to generate polynomial that satisfied requirements.");
  } 

  // Assume convergence is achieved if user does not require strict convergence.
  bool isConverged = true;	  

  // Solve the linear system using the polynomial
  {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*timerSolve_);
#endif
    
    // Apply the polynomial to the current linear system
    poly_Op_->Apply( *problem_->getRHS(), *problem_->getLHS() );
    
    // Reset the problem to acknowledge the updated solution
    problem_->setProblem(); 
    
    // If we have to strictly adhere to the requested convergence tolerance, set up a standard GMRES solver.
    if (strictConvTol_) {
      
      // Create indices for the linear systems to be solved.
      int startPtr = 0;
      int numRHS2Solve = MVT::GetNumberVecs( *(problem_->getRHS()) );
      int numCurrRHS = ( numRHS2Solve < blockSize_) ? numRHS2Solve : blockSize_;
      
      std::vector<int> currIdx;
      //  If an adaptive block size is allowed then only the linear systems that need to be solved are solved.
      //  Otherwise, the index set is generated that informs the linear problem that some linear systems are augmented.
      currIdx.resize( blockSize_ );
      for (int i=0; i<numCurrRHS; ++i) 
	{ currIdx[i] = startPtr+i; }
      for (int i=numCurrRHS; i<blockSize_; ++i) 
	{ currIdx[i] = -1; }
      
      // Inform the linear problem of the current linear system to solve.
      problem_->setLSIndex( currIdx );
      
      //////////////////////////////////////////////////////////////////////////////////////
      // Parameter list
      Teuchos::ParameterList plist;
      plist.set("Block Size",blockSize_);
      
      ptrdiff_t dim = MVText::GetGlobalLength( *(problem_->getRHS()) );  
      if (blockSize_*static_cast<ptrdiff_t>(numBlocks_) > dim) {
	ptrdiff_t tmpNumBlocks = 0;
	if (blockSize_ == 1)
	  tmpNumBlocks = dim / blockSize_;  // Allow for a good breakdown.
	else
	  tmpNumBlocks = ( dim - blockSize_) / blockSize_;  // Allow for restarting.
	printer_->stream(Warnings) << 
	  "Warning! Requested Krylov subspace dimension is larger than operator dimension!" << std::endl <<
	  " The maximum number of blocks allowed for the Krylov subspace will be adjusted to " << tmpNumBlocks << std::endl;
	plist.set("Num Blocks",Teuchos::asSafe<int>(tmpNumBlocks));
      } 
      else 
	plist.set("Num Blocks",numBlocks_);
      
      // Reset the status test.  
      outputTest_->reset();
      loaDetected_ = false;
      
      // Assume convergence is achieved, then let any failed convergence set this to false.
      isConverged = true;	
      
      //////////////////////////////////////////////////////////////////////////////////////
      // BlockGmres solver
      
      Teuchos::RCP<GmresIteration<ScalarType,MV,OP> > block_gmres_iter;
      
      block_gmres_iter = Teuchos::rcp( new BlockGmresIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,ortho_,plist) );
      
      // Enter solve() iterations
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
	Teuchos::RCP<MV> V_0;
	V_0 = MVT::CloneCopy( *(problem_->getInitPrecResVec()), currIdx );
	
	
	// Get a matrix to hold the orthonormalization coefficients.
	Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > z_0 = 
	  Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( blockSize_, blockSize_ ) );
	
	// Orthonormalize the new V_0
	int rank = ortho_->normalize( *V_0, z_0 );
	TEUCHOS_TEST_FOR_EXCEPTION(rank != blockSize_,GmresPolySolMgrOrthoFailure,
			   "Belos::GmresPolySolMgr::solve(): Failed to compute initial block of orthonormal vectors.");
	
	// Set the new state and initialize the solver.
	GmresIterationState<ScalarType,MV> newstate;
	newstate.V = V_0;
	newstate.z = z_0;
	newstate.curDim = 0;
	block_gmres_iter->initializeGmres(newstate);
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
	      if ( expConvTest_->getLOADetected() ) {
              // we don't have convergence
		loaDetected_ = true;
		isConverged = false;
	      }
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
	      
	      printer_->stream(Debug) << " Performing restart number " << numRestarts << " of " << maxRestarts_ << std::endl << std::endl;
	      
	      // Update the linear problem.
	      Teuchos::RCP<MV> update = block_gmres_iter->getCurrentUpdate();
	      problem_->updateSolution( update, true );
	      
	      // Get the state.
	      GmresIterationState<ScalarType,MV> oldState = block_gmres_iter->getState();
	      
	      // Compute the restart std::vector.
	      // Get a view of the current Krylov basis.
	      Teuchos::RCP<MV> V_0  = MVT::Clone( *(oldState.V), blockSize_ );
	      problem_->computeCurrPrecResVec( &*V_0 );
	      
	      // Get a view of the first block of the Krylov basis.
	      Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > z_0 = 
	          Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( blockSize_, blockSize_ ) );
	      
	      // Orthonormalize the new V_0
	      int rank = ortho_->normalize( *V_0, z_0 );
	      TEUCHOS_TEST_FOR_EXCEPTION(rank != blockSize_,GmresPolySolMgrOrthoFailure,
				 "Belos::GmresPolySolMgr::solve(): Failed to compute initial block of orthonormal vectors after restart.");
	      
	      // Set the new state and initialize the solver.
	      GmresIterationState<ScalarType,MV> newstate;
	      newstate.V = V_0;
	      newstate.z = z_0;
	      newstate.curDim = 0;
	      block_gmres_iter->initializeGmres(newstate);
	      
	    } // end of restarting
	    
	    ////////////////////////////////////////////////////////////////////////////////////
	    //
	    // we returned from iterate(), but none of our status tests Passed.
	    // something is wrong, and it is probably our fault.
	    //
	    ////////////////////////////////////////////////////////////////////////////////////
	    
	    else {
	      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
				 "Belos::GmresPolySolMgr::solve(): Invalid return from BlockGmresIter::iterate().");
	    }
	  }
	  catch (const GmresIterationOrthoFailure &e) {
	    // If the block size is not one, it's not considered a lucky breakdown.
	    if (blockSize_ != 1) {
	      printer_->stream(Errors) << "Error! Caught std::exception in BlockGmresIter::iterate() at iteration " 
				       << block_gmres_iter->getNumIters() << std::endl 
				       << e.what() << std::endl;
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
	  catch (const std::exception &e) {
	    printer_->stream(Errors) << "Error! Caught std::exception in BlockGmresIter::iterate() at iteration " 
				     << block_gmres_iter->getNumIters() << std::endl 
				     << e.what() << std::endl;
	    throw;
	  }
	}
	
	// Compute the current solution.
	// Update the linear problem.
	// Attempt to get the current solution from the residual status test, if it has one.
	if ( !Teuchos::is_null(expConvTest_->getSolution()) ) {
	  Teuchos::RCP<MV> newX = expConvTest_->getSolution();
	  Teuchos::RCP<MV> curX = problem_->getCurrLHSVec();
	  MVT::MvAddMv( 0.0, *newX, 1.0, *newX, *curX );
	}
	else {
	  Teuchos::RCP<MV> update = block_gmres_iter->getCurrentUpdate();
	  problem_->updateSolution( update, true );
	}  
	
	// Inform the linear problem that we are finished with this block linear system.
	problem_->setCurrLS();
	
	// Update indices for the linear systems to be solved.
	startPtr += numCurrRHS;
	numRHS2Solve -= numCurrRHS;
	if ( numRHS2Solve > 0 ) {
	  numCurrRHS = ( numRHS2Solve < blockSize_) ? numRHS2Solve : blockSize_;
	  
	  currIdx.resize( blockSize_ );
	  for (int i=0; i<numCurrRHS; ++i) 
	    { currIdx[i] = startPtr+i; }
	  for (int i=numCurrRHS; i<blockSize_; ++i) 
	    { currIdx[i] = -1; }
	  
	  // Set the next indices.
	  problem_->setLSIndex( currIdx );
	}
	else {
	  currIdx.resize( numRHS2Solve );
	}
	
      }// while ( numRHS2Solve > 0 )
      
      // print final summary
      sTest_->print( printer_->stream(FinalSummary) );
      
    } // if (strictConvTol_)
  } // timing block

  // print timing information
#ifdef BELOS_TEUCHOS_TIME_MONITOR
  // Calling summarize() can be expensive, so don't call unless the
  // user wants to print out timing details.  summarize() will do all
  // the work even if it's passed a "black hole" output stream.
  if (verbosity_ & TimingDetails)
    Teuchos::TimeMonitor::summarize( printer_->stream(TimingDetails) );
#endif
  
  if (!isConverged || loaDetected_) {
    return Unconverged; // return from GmresPolySolMgr::solve() 
  }
  return Converged; // return from GmresPolySolMgr::solve() 
}
  
//  This method requires the solver manager to return a std::string that describes itself.
template<class ScalarType, class MV, class OP>
std::string GmresPolySolMgr<ScalarType,MV,OP>::description() const
  {
  std::ostringstream oss;
  oss << "Belos::GmresPolySolMgr<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
  oss << "{";
  oss << "Ortho Type='"<<orthoType_<<"\', Block Size=" <<blockSize_;
  oss << ", Num Blocks=" <<numBlocks_<< ", Max Restarts=" << maxRestarts_;
  oss << "}";
  return oss.str();
}
  
} // end Belos namespace

#endif /* BELOS_GMRES_POLY_SOLMGR_HPP */
