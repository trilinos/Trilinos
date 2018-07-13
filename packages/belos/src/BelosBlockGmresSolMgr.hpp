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

#ifndef BELOS_BLOCK_GMRES_SOLMGR_HPP
#define BELOS_BLOCK_GMRES_SOLMGR_HPP

/*! \file BelosBlockGmresSolMgr.hpp
 *  \brief The Belos::BlockGmresSolMgr provides a solver manager for the BlockGmres linear solver.
 */

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosGmresIteration.hpp"
#include "BelosBlockGmresIter.hpp"
#include "BelosBlockFGmresIter.hpp"
#include "BelosDGKSOrthoManager.hpp"
#include "BelosICGSOrthoManager.hpp"
#include "BelosIMGSOrthoManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestGenResNorm.hpp"
#include "BelosStatusTestImpResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutput.hpp"
#include "BelosStatusTestOutputFactory.hpp"
#include "BelosOutputManager.hpp"
#include "Teuchos_BLAS.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

/** \example BlockGmres/BlockGmresEpetraExFile.cpp
    This is an example of how to use the Belos::BlockGmresSolMgr solver manager.
*/
/** \example BlockGmres/BlockPrecGmresEpetraExFile.cpp
    This is an example of how to use the Belos::BlockGmresSolMgr solver manager with an Ifpack preconditioner.
*/
/** \example BlockGmres/BlockFlexGmresEpetraExFile.cpp
    This is an example of how to use the Belos::BlockGmresSolMgr solver manager with flexible Gmres.
*/

namespace Belos {

//! @name BlockGmresSolMgr Exceptions
//@{

/** \brief BlockGmresSolMgrLinearProblemFailure is thrown when the linear problem is
 * not setup (i.e. setProblem() was not called) when solve() is called.
 *
 * This std::exception is thrown from the BlockGmresSolMgr::solve() method.
 *
 */
class BlockGmresSolMgrLinearProblemFailure : public BelosError {public:
  BlockGmresSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

/** \brief BlockGmresSolMgrOrthoFailure is thrown when the orthogonalization manager is
 * unable to generate orthonormal columns from the initial basis vectors.
 *
 * This std::exception is thrown from the BlockGmresSolMgr::solve() method.
 *
 */
class BlockGmresSolMgrOrthoFailure : public BelosError {public:
  BlockGmresSolMgrOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

/*! \class BlockGmresSolMgr
 * \brief Interface to Block GMRES and Flexible GMRES.
 * \ingroup belos_solver_framework
 * \author Heidi Thornquist, Chris Baker, and Teri Barth
 *
 * This class provides an interface to Block GMRES, for solving linear
 * systems with one or more right-hand sides.  Our Block GMRES
 * implementation also has an option (the Boolean "Flexible Gmres"
 * parameter) to use the Flexible variant of GMRES.  Flexible GMRES
 * allows the preconditioner (which must be a right preconditioner) to
 * change from iteration to iteration.
 *
 * If you are a new Belos user and just want standard GMRES, use \c
 * PseudoBlockGmresSolMgr.  If you want Flexible GMRES, use this class
 * with the "Flexible Gmres" parameter set to true.
 */
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
   *   - "Block Size" - a \c int specifying the block size to be used by the underlying block GMRES solver. Default: 1
   *   - "Adaptive Block Size" - a \c bool specifying whether the block size can be modified throughout the solve. Default: true
   *   - "Num Blocks" - a \c int specifying the number of blocks allocated for the Krylov basis. Default: 300
   *   - "Maximum Iterations" - a \c int specifying the maximum number of iterations the underlying solver is allowed to perform. Default: 1000
   *   - "Maximum Restarts" - a \c int specifying the maximum number of restarts the underlying solver is allowed to perform. Default: 20
   *   - "Orthogonalization" - a \c std::string specifying the desired orthogonalization:  DGKS, ICGS, and IMGS. Default: "DGKS"
   *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Belos::Errors
   *   - "Output Style" - a OutputType specifying the style of output. Default: Belos::General
   *   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: 1e-8
   *            \param pl [in] ParameterList with construction information
   *                    \htmlonly
   *                    <iframe src="belos_BlockGmres.xml" width=100% scrolling="no" frameborder="0">
   *                    </iframe>
   *                    <hr />
   *                    \endhtmlonly
   */
  BlockGmresSolMgr( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
    const Teuchos::RCP<Teuchos::ParameterList> &pl );

  //! Destructor.
  virtual ~BlockGmresSolMgr() {};

  //! clone for Inverted Injection (DII)
  Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const override {
    return Teuchos::rcp(new BlockGmresSolMgr<ScalarType,MV,OP>);
  }
  //@}

  //! @name Accessor methods
  //@{

  /*! \brief Get current linear problem being solved for in this object.
   */
  const LinearProblem<ScalarType,MV,OP>& getProblem() const override {
    return *problem_;
  }

  /*! \brief Get a parameter list containing the valid parameters for this object.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

  /*! \brief Get a parameter list containing the current parameters for this object.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const override { return params_; }

  /*! \brief Return the timers for this object.
   *
   * The timers are ordered as follows:
   *   - time spent in solve() routine
   */
  Teuchos::Array<Teuchos::RCP<Teuchos::Time> > getTimers() const {
    return Teuchos::tuple(timerSolve_);
  }

  /// \brief Tolerance achieved by the last \c solve() invocation.
  ///
  /// This is the maximum over all right-hand sides' achieved
  /// convergence tolerances, and is set whether or not the solve
  /// actually managed to achieve the desired convergence tolerance.
  ///
  /// \warning This result may not be meaningful if there was a loss
  ///   of accuracy during the solve.  You should first call \c
  ///   isLOADetected() to check for a loss of accuracy during the
  ///   last solve.
  MagnitudeType achievedTol() const override {
    return achievedTol_;
  }

  //! Get the iteration count for the most recent call to \c solve().
  int getNumIters() const override {
    return numIters_;
  }

  /*! \brief Return whether a loss of accuracy was detected by this solver during the most current solve.
      \note This flag will be reset the next time solve() is called.
   */
  bool isLOADetected() const override { return loaDetected_; }

  //@}

  //! @name Set methods
  //@{

  //! Set the linear problem that needs to be solved.
  void setProblem( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem ) override { problem_ = problem; isSTSet_ = false; }

  //! Set the parameters the solver manager should use to solve the linear problem.
  void setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params ) override;

  //! Set a debug status test that will be checked at the same time as the top-level status test.
  void setDebugStatusTest( const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &debugStatusTest ) override;

  //@}

  //! @name Reset methods
  //@{
  /*! \brief Performs a reset of the solver manager specified by the \c ResetType.  This informs the
   *  solver manager that the solver should prepare for the next call to solve by resetting certain elements
   *  of the iterative solver strategy.
  */
  void reset( const ResetType type ) override { if ((type & Belos::Problem) && !Teuchos::is_null(problem_)) problem_->setProblem(); }
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
  ReturnType solve() override;

  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{

  /// \brief Print the object with the given verbosity level to a FancyOStream.
  ///
  /// \param out [out] Output stream to which to print.
  ///
  /// \param verbLevel [in] Verbosity level.  The default verbosity
  ///   (verbLevel=Teuchos::VERB_DEFAULT) is Teuchos::VERB_LOW.
  void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const override;

  //! Return a one-line description of this object.
  std::string description () const override;

  //@}

private:

  // Method for checking current status test against defined linear problem.
  bool checkStatusTest();

  // Linear problem.
  Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;

  // Output manager.
  Teuchos::RCP<OutputManager<ScalarType> > printer_;
  Teuchos::RCP<std::ostream> outputStream_;

  // Status test.
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > debugStatusTest_;
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
  static constexpr int maxRestarts_default_ = 20;
  static constexpr int maxIters_default_ = 1000;
  static constexpr bool adaptiveBlockSize_default_ = true;
  static constexpr bool showMaxResNormOnly_default_ = false;
  static constexpr bool flexibleGmres_default_ = false;
  static constexpr bool expResTest_default_ = false;
  static constexpr int blockSize_default_ = 1;
  static constexpr int numBlocks_default_ = 300;
  static constexpr int verbosity_default_ = Belos::Errors;
  static constexpr int outputStyle_default_ = Belos::General;
  static constexpr int outputFreq_default_ = -1;
  static constexpr const char * impResScale_default_ = "Norm of Preconditioned Initial Residual";
  static constexpr const char * expResScale_default_ = "Norm of Initial Residual";
  static constexpr const char * label_default_ = "Belos";
  static constexpr const char * orthoType_default_ = "DGKS";
  static constexpr std::ostream * outputStream_default_ = &std::cout;

  // Current solver values.
  MagnitudeType convtol_, orthoKappa_, achievedTol_;
  int maxRestarts_, maxIters_, numIters_;
  int blockSize_, numBlocks_, verbosity_, outputStyle_, outputFreq_;
  bool adaptiveBlockSize_, showMaxResNormOnly_, isFlexible_, expResTest_;
  std::string orthoType_;
  std::string impResScale_, expResScale_;

  // Timers.
  std::string label_;
  Teuchos::RCP<Teuchos::Time> timerSolve_;

  // Internal state variables.
  bool isSet_, isSTSet_;
  bool loaDetected_;
};


// Empty Constructor
template<class ScalarType, class MV, class OP>
BlockGmresSolMgr<ScalarType,MV,OP>::BlockGmresSolMgr() :
  outputStream_(Teuchos::rcp(outputStream_default_,false)),
  convtol_(DefaultSolverParameters::convTol),
  orthoKappa_(DefaultSolverParameters::orthoKappa),
  achievedTol_(Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>::zero()),
  maxRestarts_(maxRestarts_default_),
  maxIters_(maxIters_default_),
  numIters_(0),
  blockSize_(blockSize_default_),
  numBlocks_(numBlocks_default_),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  adaptiveBlockSize_(adaptiveBlockSize_default_),
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  isFlexible_(flexibleGmres_default_),
  expResTest_(expResTest_default_),
  orthoType_(orthoType_default_),
  impResScale_(impResScale_default_),
  expResScale_(expResScale_default_),
  label_(label_default_),
  isSet_(false),
  isSTSet_(false),
  loaDetected_(false)
{}


// Basic Constructor
template<class ScalarType, class MV, class OP>
BlockGmresSolMgr<ScalarType,MV,OP>::
BlockGmresSolMgr (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                  const Teuchos::RCP<Teuchos::ParameterList> &pl) :
  problem_(problem),
  outputStream_(Teuchos::rcp(outputStream_default_,false)),
  convtol_(DefaultSolverParameters::convTol),
  orthoKappa_(DefaultSolverParameters::orthoKappa),
  achievedTol_(Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>::zero()),
  maxRestarts_(maxRestarts_default_),
  maxIters_(maxIters_default_),
  numIters_(0),
  blockSize_(blockSize_default_),
  numBlocks_(numBlocks_default_),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  adaptiveBlockSize_(adaptiveBlockSize_default_),
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  isFlexible_(flexibleGmres_default_),
  expResTest_(expResTest_default_),
  orthoType_(orthoType_default_),
  impResScale_(impResScale_default_),
  expResScale_(expResScale_default_),
  label_(label_default_),
  isSet_(false),
  isSTSet_(false),
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
BlockGmresSolMgr<ScalarType,MV,OP>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

    // The static_cast is to resolve an issue with older clang versions which
    // would cause the constexpr to link fail. With c++17 the problem is resolved.
    pl->set("Convergence Tolerance", static_cast<MagnitudeType>(DefaultSolverParameters::convTol),
      "The relative residual tolerance that needs to be achieved by the\n"
      "iterative solver in order for the linear system to be declared converged." );
    pl->set("Maximum Restarts", static_cast<int>(maxRestarts_default_),
      "The maximum number of restarts allowed for each\n"
      "set of RHS solved.");
    pl->set("Maximum Iterations", static_cast<int>(maxIters_default_),
      "The maximum number of block iterations allowed for each\n"
      "set of RHS solved.");
    pl->set("Num Blocks", static_cast<int>(numBlocks_default_),
      "The maximum number of blocks allowed in the Krylov subspace\n"
      "for each set of RHS solved.");
    pl->set("Block Size", static_cast<int>(blockSize_default_),
      "The number of vectors in each block.  This number times the\n"
      "number of blocks is the total Krylov subspace dimension.");
    pl->set("Adaptive Block Size", static_cast<bool>(adaptiveBlockSize_default_),
      "Whether the solver manager should adapt the block size\n"
      "based on the number of RHS to solve.");
    pl->set("Verbosity", static_cast<int>(verbosity_default_),
      "What type(s) of solver information should be outputted\n"
      "to the output stream.");
    pl->set("Output Style", static_cast<int>(outputStyle_default_),
      "What style is used for the solver information outputted\n"
      "to the output stream.");
    pl->set("Output Frequency", static_cast<int>(outputFreq_default_),
      "How often convergence information should be outputted\n"
      "to the output stream.");
    pl->set("Output Stream", Teuchos::rcp(outputStream_default_,false),
      "A reference-counted pointer to the output stream where all\n"
      "solver output is sent.");
    pl->set("Show Maximum Residual Norm Only", static_cast<bool>(showMaxResNormOnly_default_),
      "When convergence information is printed, only show the maximum\n"
      "relative residual norm when the block size is greater than one.");
    pl->set("Flexible Gmres", static_cast<bool>(flexibleGmres_default_),
      "Whether the solver manager should use the flexible variant\n"
      "of GMRES.");
    pl->set("Explicit Residual Test", static_cast<bool>(expResTest_default_),
      "Whether the explicitly computed residual should be used in the convergence test.");
    pl->set("Implicit Residual Scaling", static_cast<const char *>(impResScale_default_),
      "The type of scaling used in the implicit residual convergence test.");
    pl->set("Explicit Residual Scaling", static_cast<const char *>(expResScale_default_),
      "The type of scaling used in the explicit residual convergence test.");
    pl->set("Timer Label", static_cast<const char *>(label_default_),
      "The string to use as a prefix for the timer labels.");
    pl->set("Orthogonalization", static_cast<const char *>(orthoType_default_),
      "The type of orthogonalization to use: DGKS, ICGS, or IMGS.");
    pl->set("Orthogonalization Constant",static_cast<MagnitudeType>(DefaultSolverParameters::orthoKappa),
      "The constant used by DGKS orthogonalization to determine\n"
      "whether another step of classical Gram-Schmidt is necessary.");
    validPL = pl;
  }
  return validPL;
}


template<class ScalarType, class MV, class OP>
void BlockGmresSolMgr<ScalarType,MV,OP>::setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params )
{

  // Create the internal parameter list if ones doesn't already exist.
  if (params_ == Teuchos::null) {
    params_ = Teuchos::rcp( new Teuchos::ParameterList(*getValidParameters()) );
  }
  else {
    params->validateParameters(*getValidParameters());
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
      "Belos::BlockGmresSolMgr: \"Block Size\" must be strictly positive.");

    // Update parameter in our list.
    params_->set("Block Size", blockSize_);
  }

  // Check if the blocksize should be adaptive
  if (params->isParameter("Adaptive Block Size")) {
    adaptiveBlockSize_ = params->get("Adaptive Block Size",adaptiveBlockSize_default_);

    // Update parameter in our list.
    params_->set("Adaptive Block Size", adaptiveBlockSize_);
  }

  // Check for the maximum number of blocks.
  if (params->isParameter("Num Blocks")) {
    numBlocks_ = params->get("Num Blocks",numBlocks_default_);
    TEUCHOS_TEST_FOR_EXCEPTION(numBlocks_ <= 0, std::invalid_argument,
      "Belos::BlockGmresSolMgr: \"Num Blocks\" must be strictly positive.");

    // Update parameter in our list.
    params_->set("Num Blocks", numBlocks_);
  }

  // Check to see if the timer label changed.
  if (params->isParameter("Timer Label")) {
    std::string tempLabel = params->get("Timer Label", label_default_);

    // Update parameter in our list, solver timer, and orthogonalization label
    if (tempLabel != label_) {
      label_ = tempLabel;
      params_->set("Timer Label", label_);
      std::string solveLabel = label_ + ": BlockGmresSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerSolve_ = Teuchos::TimeMonitor::getNewCounter(solveLabel);
#endif
      if (ortho_ != Teuchos::null) {
        ortho_->setLabel( label_ );
      }
    }
  }

  // Determine whether this solver should be "flexible".
  if (params->isParameter("Flexible Gmres")) {
    isFlexible_ = Teuchos::getParameter<bool>(*params,"Flexible Gmres");
    params_->set("Flexible Gmres", isFlexible_);
    if (isFlexible_ && expResTest_) {
      // Use an implicit convergence test if the Gmres solver is flexible
      isSTSet_ = false;
    }
  }


  // Check if the orthogonalization changed.
  if (params->isParameter("Orthogonalization")) {
    std::string tempOrthoType = params->get("Orthogonalization",orthoType_default_);
    TEUCHOS_TEST_FOR_EXCEPTION( tempOrthoType != "DGKS" && tempOrthoType != "ICGS" && tempOrthoType != "IMGS",
                        std::invalid_argument,
                        "Belos::BlockGmresSolMgr: \"Orthogonalization\" must be either \"DGKS\", \"ICGS\", or \"IMGS\".");
    if (tempOrthoType != orthoType_) {
      orthoType_ = tempOrthoType;
      params_->set("Orthogonalization", orthoType_);
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
    if (params->isType<MagnitudeType> ("Orthogonalization Constant")) {
      orthoKappa_ = params->get ("Orthogonalization Constant",
                                 static_cast<MagnitudeType> (DefaultSolverParameters::orthoKappa));
    }
    else {
      orthoKappa_ = params->get ("Orthogonalization Constant",
                                 DefaultSolverParameters::orthoKappa);
    }

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

  // Check for convergence tolerance
  if (params->isParameter("Convergence Tolerance")) {
    if (params->isType<MagnitudeType> ("Convergence Tolerance")) {
      convtol_ = params->get ("Convergence Tolerance",
                              static_cast<MagnitudeType> (DefaultSolverParameters::convTol));
    }
    else {
      convtol_ = params->get ("Convergence Tolerance", DefaultSolverParameters::convTol);
    }

    // Update parameter in our list and residual tests.
    params_->set("Convergence Tolerance", convtol_);
    if (impConvTest_ != Teuchos::null)
      impConvTest_->setTolerance( convtol_ );
    if (expConvTest_ != Teuchos::null)
      expConvTest_->setTolerance( convtol_ );
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

  if (params->isParameter("Explicit Residual Test")) {
    expResTest_ = Teuchos::getParameter<bool>( *params,"Explicit Residual Test" );

    // Reconstruct the convergence test if the explicit residual test is not being used.
    params_->set("Explicit Residual Test", expResTest_);
    if (expConvTest_ == Teuchos::null) {
      isSTSet_ = false;
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
    params_->set("Orthogonalization", orthoType_);
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
        "Belos::BlockGmresSolMgr(): Invalid orthogonalization type.");
    }
  }

  // Create the timer if we need to.
  if (timerSolve_ == Teuchos::null) {
    std::string solveLabel = label_ + ": BlockGmresSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    timerSolve_ = Teuchos::TimeMonitor::getNewCounter(solveLabel);
#endif
  }

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}

// Check the status test versus the defined linear problem
template<class ScalarType, class MV, class OP>
bool BlockGmresSolMgr<ScalarType,MV,OP>::checkStatusTest() {

  typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
  typedef Belos::StatusTestGenResNorm<ScalarType,MV,OP>  StatusTestGenResNorm_t;
  typedef Belos::StatusTestImpResNorm<ScalarType,MV,OP>  StatusTestImpResNorm_t;

  // Basic test checks maximum iterations and native residual.
  maxIterTest_ = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>( maxIters_ ) );

  // If there is a left preconditioner, we create a combined status test that checks the implicit
  // and then explicit residual norm to see if we have convergence.
  if (!Teuchos::is_null(problem_->getLeftPrec()) && !isFlexible_) {
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

    if (isFlexible_) {
      // Implicit residual test, using the native residual to determine if convergence was achieved.
      Teuchos::RCP<StatusTestGenResNorm_t> tmpImpConvTest =
        Teuchos::rcp( new StatusTestGenResNorm_t( convtol_ ) );
      tmpImpConvTest->defineScaleForm( convertStringToScaleType(impResScale_), Belos::TwoNorm );
      tmpImpConvTest->setShowMaxResNormOnly( showMaxResNormOnly_ );
      impConvTest_ = tmpImpConvTest;
    }
    else {
      // Implicit residual test, using the native residual to determine if convergence was achieved.
      // Use test that checks for loss of accuracy.
      Teuchos::RCP<StatusTestImpResNorm_t> tmpImpConvTest =
        Teuchos::rcp( new StatusTestImpResNorm_t( convtol_ ) );
      tmpImpConvTest->defineScaleForm( convertStringToScaleType(impResScale_), Belos::TwoNorm );
      tmpImpConvTest->setShowMaxResNormOnly( showMaxResNormOnly_ );
      impConvTest_ = tmpImpConvTest;
    }

    // Set the explicit and total convergence test to this implicit test that checks for accuracy loss.
    expConvTest_ = impConvTest_;
    convTest_ = impConvTest_;
  }

  // Create the status test.
  sTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::OR, maxIterTest_, convTest_ ) );

  // Add debug status test, if one is provided by the user
  if (nonnull(debugStatusTest_) ) {
    // Add debug convergence test
    Teuchos::rcp_dynamic_cast<StatusTestCombo_t>(sTest_)->addStatusTest( debugStatusTest_ );
  }

  // Create the status test output class.
  // This class manages and formats the output from the status test.
  StatusTestOutputFactory<ScalarType,MV,OP> stoFactory( outputStyle_ );
  outputTest_ = stoFactory.create( printer_, sTest_, outputFreq_, Passed+Failed+Undefined );

  // Set the solver string for the output test
  std::string solverDesc = " Block Gmres ";
  if (isFlexible_)
    solverDesc = "Flexible" + solverDesc;
  outputTest_->setSolverDesc( solverDesc );

  // The status test is now set.
  isSTSet_ = true;

  return false;
}

template<class ScalarType, class MV, class OP>
void BlockGmresSolMgr<ScalarType,MV,OP>::setDebugStatusTest(
  const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &debugStatusTest
  )
{
  debugStatusTest_ = debugStatusTest;
}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType BlockGmresSolMgr<ScalarType,MV,OP>::solve() {

  // Set the current parameters if they were not set before.
  // NOTE:  This may occur if the user generated the solver manager with the default constructor and
  // then didn't set any parameters using setParameters().
  if (!isSet_) {
    setParameters(Teuchos::parameterList(*getValidParameters()));
  }

  Teuchos::BLAS<int,ScalarType> blas;

  TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null,BlockGmresSolMgrLinearProblemFailure,
    "Belos::BlockGmresSolMgr::solve(): Linear problem is not a valid object.");

  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isProblemSet(),BlockGmresSolMgrLinearProblemFailure,
    "Belos::BlockGmresSolMgr::solve(): Linear problem is not ready, setProblem() has not been called.");

  if (isFlexible_) {
    TEUCHOS_TEST_FOR_EXCEPTION(problem_->getRightPrec()==Teuchos::null,BlockGmresSolMgrLinearProblemFailure,
      "Belos::BlockGmresSolMgr::solve(): Linear problem does not have a preconditioner required for flexible GMRES, call setRightPrec().");
  }

  if (!isSTSet_ || (!expResTest_ && !Teuchos::is_null(problem_->getLeftPrec())) ) {
    TEUCHOS_TEST_FOR_EXCEPTION( checkStatusTest(),BlockGmresSolMgrLinearProblemFailure,
      "Belos::BlockGmresSolMgr::solve(): Linear problem and requested status tests are incompatible.");
  }

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

  ptrdiff_t dim = MVT::GetGlobalLength( *(problem_->getRHS()) );
  if (blockSize_*static_cast<ptrdiff_t>(numBlocks_) > dim) {
    int tmpNumBlocks = 0;
    if (blockSize_ == 1)
      tmpNumBlocks = dim / blockSize_;  // Allow for a good breakdown.
    else
      tmpNumBlocks = ( dim - blockSize_) / blockSize_;  // Allow for restarting.
    printer_->stream(Warnings) <<
      "Belos::BlockGmresSolMgr::solve():  Warning! Requested Krylov subspace dimension is larger than operator dimension!"
      << std::endl << " The maximum number of blocks allowed for the Krylov subspace will be adjusted to " << tmpNumBlocks << std::endl;
    plist.set("Num Blocks",tmpNumBlocks);
  }
  else
    plist.set("Num Blocks",numBlocks_);

  // Reset the status test.
  outputTest_->reset();
  loaDetected_ = false;

  // Assume convergence is achieved, then let any failed convergence set this to false.
  bool isConverged = true;

  //////////////////////////////////////////////////////////////////////////////////////
  // BlockGmres solver

  Teuchos::RCP<GmresIteration<ScalarType,MV,OP> > block_gmres_iter;

  if (isFlexible_)
    block_gmres_iter = Teuchos::rcp( new BlockFGmresIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,ortho_,plist) );
  else
    block_gmres_iter = Teuchos::rcp( new BlockGmresIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,ortho_,plist) );

  // Enter solve() iterations
  {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*timerSolve_);
#endif

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
      if (isFlexible_) {
        // Load the correct residual if the system is augmented
        if (currIdx[blockSize_-1] == -1) {
          V_0 = MVT::Clone( *(problem_->getInitResVec()), blockSize_ );
          problem_->computeCurrResVec( &*V_0 );
        }
        else {
          V_0 = MVT::CloneCopy( *(problem_->getInitResVec()), currIdx );
        }
      }
      else {
        // Load the correct residual if the system is augmented
        if (currIdx[blockSize_-1] == -1) {
          V_0 = MVT::Clone( *(problem_->getInitPrecResVec()), blockSize_ );
          problem_->computeCurrPrecResVec( &*V_0 );
        }
        else {
          V_0 = MVT::CloneCopy( *(problem_->getInitPrecResVec()), currIdx );
        }
      }

      // Get a matrix to hold the orthonormalization coefficients.
      Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > z_0 =
        Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( blockSize_, blockSize_ ) );

      // Orthonormalize the new V_0
      int rank = ortho_->normalize( *V_0, z_0 );
      TEUCHOS_TEST_FOR_EXCEPTION(rank != blockSize_,BlockGmresSolMgrOrthoFailure,
        "Belos::BlockGmresSolMgr::solve(): Failed to compute initial block of orthonormal vectors.");

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
              printer_->stream(Warnings) <<
                "Belos::BlockGmresSolMgr::solve(): Warning! Solver has experienced a loss of accuracy!" << std::endl;
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
            if (isFlexible_) {
              // Update the solution manually, since the preconditioning doesn't need to be undone.
              Teuchos::RCP<MV> curX = problem_->getCurrLHSVec();
              MVT::MvAddMv( 1.0, *curX, 1.0, *update, *curX );
            }
            else
              problem_->updateSolution( update, true );

            // Get the state.
            GmresIterationState<ScalarType,MV> oldState = block_gmres_iter->getState();

            // Compute the restart std::vector.
            // Get a view of the current Krylov basis.
            V_0  = MVT::Clone( *(oldState.V), blockSize_ );
            if (isFlexible_)
              problem_->computeCurrResVec( &*V_0 );
            else
              problem_->computeCurrPrecResVec( &*V_0 );

            // Get a view of the first block of the Krylov basis.
            z_0 = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( blockSize_, blockSize_ ) );

            // Orthonormalize the new V_0
            rank = ortho_->normalize( *V_0, z_0 );
            TEUCHOS_TEST_FOR_EXCEPTION(rank != blockSize_,BlockGmresSolMgrOrthoFailure,
              "Belos::BlockGmresSolMgr::solve(): Failed to compute initial block of orthonormal vectors after restart.");

            // Set the new state and initialize the solver.
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
              "Belos::BlockGmresSolMgr::solve(): Invalid return from BlockGmresIter::iterate().");
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
      if (isFlexible_) {
        // Update the solution manually, since the preconditioning doesn't need to be undone.
        Teuchos::RCP<MV> update = block_gmres_iter->getCurrentUpdate();
        Teuchos::RCP<MV> curX = problem_->getCurrLHSVec();
        // Update the solution only if there is a valid update from the iteration
        if (update != Teuchos::null)
          MVT::MvAddMv( 1.0, *curX, 1.0, *update, *curX );
      }
      else {
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
#ifdef BELOS_TEUCHOS_TIME_MONITOR
  // Calling summarize() can be expensive, so don't call unless the
  // user wants to print out timing details.  summarize() will do all
  // the work even if it's passed a "black hole" output stream.
  if (verbosity_ & TimingDetails)
    Teuchos::TimeMonitor::summarize( printer_->stream(TimingDetails) );
#endif

  // get iteration information for this solve
  numIters_ = maxIterTest_->getNumIters();

  // Save the convergence test value ("achieved tolerance") for this
  // solve.  This requires a bit more work than for BlockCGSolMgr,
  // since for this solver, convTest_ may either be a single residual
  // norm test, or a combination of two residual norm tests.  In the
  // latter case, the master convergence test convTest_ is a SEQ combo
  // of the implicit resp. explicit tests.  If the implicit test never
  // passes, then the explicit test won't ever be executed.  This
  // manifests as expConvTest_->getTestValue()->size() < 1.  We deal
  // with this case by using the values returned by
  // impConvTest_->getTestValue().
  {
    // We'll fetch the vector of residual norms one way or the other.
    const std::vector<MagnitudeType>* pTestValues = NULL;
    if (expResTest_) {
      pTestValues = expConvTest_->getTestValue();
      if (pTestValues == NULL || pTestValues->size() < 1) {
        pTestValues = impConvTest_->getTestValue();
      }
    }
    else {
      // Only the implicit residual norm test is being used.
      pTestValues = impConvTest_->getTestValue();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(pTestValues == NULL, std::logic_error,
      "Belos::BlockGmresSolMgr::solve(): The implicit convergence test's "
      "getTestValue() method returned NULL.  Please report this bug to the "
      "Belos developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(pTestValues->size() < 1, std::logic_error,
      "Belos::BlockGmresSolMgr::solve(): The implicit convergence test's "
      "getTestValue() method returned a vector of length zero.  Please report "
      "this bug to the Belos developers.");

    // FIXME (mfh 12 Dec 2011) Does pTestValues really contain the
    // achieved tolerances for all vectors in the current solve(), or
    // just for the vectors from the last deflation?
    achievedTol_ = *std::max_element (pTestValues->begin(), pTestValues->end());
  }

  if (!isConverged || loaDetected_) {
    return Unconverged; // return from BlockGmresSolMgr::solve()
  }
  return Converged; // return from BlockGmresSolMgr::solve()
}


template<class ScalarType, class MV, class OP>
std::string BlockGmresSolMgr<ScalarType,MV,OP>::description() const
{
  std::ostringstream out;
  out << "\"Belos::BlockGmresSolMgr\": {";
  if (this->getObjectLabel () != "") {
    out << "Label: " << this->getObjectLabel () << ", ";
  }
  out << "Flexible: " << (isFlexible_ ? "true" : "false")
      << ", Num Blocks: " << numBlocks_
      << ", Maximum Iterations: " << maxIters_
      << ", Maximum Restarts: " << maxRestarts_
      << ", Convergence Tolerance: " << convtol_
      << "}";
  return out.str ();
}


template<class ScalarType, class MV, class OP>
void
BlockGmresSolMgr<ScalarType, MV, OP>::
describe (Teuchos::FancyOStream &out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::TypeNameTraits;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  // using Teuchos::VERB_MEDIUM;
  // using Teuchos::VERB_HIGH;
  // using Teuchos::VERB_EXTREME;
  using std::endl;

  // Set default verbosity if applicable.
  const Teuchos::EVerbosityLevel vl =
    (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

  if (vl != VERB_NONE) {
    Teuchos::OSTab tab0 (out);

    out << "\"Belos::BlockGmresSolMgr\":" << endl;
    Teuchos::OSTab tab1 (out);
    out << "Template parameters:" << endl;
    {
      Teuchos::OSTab tab2 (out);
      out << "ScalarType: " << TypeNameTraits<ScalarType>::name () << endl
          << "MV: " << TypeNameTraits<MV>::name () << endl
          << "OP: " << TypeNameTraits<OP>::name () << endl;
    }
    if (this->getObjectLabel () != "") {
      out << "Label: " << this->getObjectLabel () << endl;
    }
    out << "Flexible: " << (isFlexible_ ? "true" : "false") << endl
        << "Num Blocks: " << numBlocks_ << endl
        << "Maximum Iterations: " << maxIters_ << endl
        << "Maximum Restarts: " << maxRestarts_ << endl
        << "Convergence Tolerance: " << convtol_ << endl;
  }
}


} // end Belos namespace

#endif /* BELOS_BLOCK_GMRES_SOLMGR_HPP */
