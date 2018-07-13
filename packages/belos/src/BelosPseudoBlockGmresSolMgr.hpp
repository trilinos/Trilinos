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
#include "BelosIMGSOrthoManager.hpp"
#ifdef HAVE_BELOS_TSQR
#  include "BelosTsqrOrthoManager.hpp"
#endif // HAVE_BELOS_TSQR
#include "BelosStatusTestFactory.hpp"
#include "BelosStatusTestOutputFactory.hpp"
#include "BelosOutputManager.hpp"
#include "Teuchos_BLAS.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

/** \example BlockGmres/PseudoBlockGmresEpetraExFile.cpp
    This is an example of how to use the Belos::PseudoBlockGmresSolMgr solver manager.
*/
/** \example BlockGmres/PseudoBlockPrecGmresEpetraExFile.cpp
    This is an example of how to use the Belos::PseudoBlockGmresSolMgr solver manager with an Ifpack preconditioner.
*/

namespace Belos {

  //! @name PseudoBlockGmresSolMgr Exceptions
  //@{

  /** \brief PseudoBlockGmresSolMgrLinearProblemFailure is thrown when the linear problem is
   * not setup (i.e. setProblem() was not called) when solve() is called.
   *
   * This std::exception is thrown from the PseudoBlockGmresSolMgr::solve() method.
   *
   */
  class PseudoBlockGmresSolMgrLinearProblemFailure : public BelosError {public:
    PseudoBlockGmresSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief PseudoBlockGmresSolMgrOrthoFailure is thrown when the orthogonalization manager is
   * unable to generate orthonormal columns from the initial basis vectors.
   *
   * This std::exception is thrown from the PseudoBlockGmresSolMgr::solve() method.
   *
   */
  class PseudoBlockGmresSolMgrOrthoFailure : public BelosError {public:
    PseudoBlockGmresSolMgrOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /*! \class PseudoBlockGmresSolMgr
   * \brief Interface to standard and "pseudoblock" GMRES.
   * \author Heidi Thornquist, Chris Baker, and Teri Barth
   * \ingroup belos_solver_framework
   *
   * This class provides an interface to the following iterative solvers:
   * - GMRES, for linear systems with one right-hand side
   * - The "pseudoblock" variant of GMRES, for linear systems
   *   with multiple right-hand sides
   *
   * If you are a new Belos user and just want standard GMRES, use
   * this class.  If you want Flexible GMRES, use \c BlockGmresSolMgr
   * with the appropriate option set.
   *
   * "Pseudoblock" GMRES is a way to improve performance when solving
   * systems with multiple right-hand sides, without changing the
   * convergence characteristics.  It is equivalent in terms of
   * convergence to running a separate instance of (standard) GMRES
   * for each right-hand side, but should often be faster.  When
   * solving for multiple right-hand sides, "Block GMRES" (as
   * implemented by \c BlockGmresSolMgr) is a different algorithm with
   * different convergence characteristics than Pseudoblock GMRES.
   */
  template<class ScalarType, class MV, class OP>
  class PseudoBlockGmresSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;

  public:

    //! @name Constructors and destructor
    //@{

    /*! \brief Empty constructor.
     *
     * This constructor takes no arguments.  It sets default solver
     * parameters, which you may change by calling setParameters().
     * Before you may call solve(), you must first give the solver a
     * linear problem to solve, by calling setProblem().
     */
    PseudoBlockGmresSolMgr();

    /*! \brief Constructor that takes the problem to solve, and a list
     *    of solver options.
     *
     * \param problem [in/out] The linear problem to be solved.
     * \param pl [in/out] A list of solver options.
     *
     * Belos' solvers accept many different options.  You may accept
     * their default values, or set any of them yourself.  We will
     * explain the options by category.
     *
     * The following options govern the number of iterations and
     * restarts:
     * - "Num Blocks" (\c int): The restart length.  The number of
     *   vectors (or blocks, in the case of multiple right-hand sides)
     *   allocated for the Krylov basis.  Its default value is 300.
     * - "Maximum Iterations" (\c int): The maximum number of
     *   iterations the solver is allowed to perform.  This does
     *   <i>not</i> include computing the initial residual, but it
     *   <i>does</i> include iterations before and after any restarts.
     *   Its default value is 1000.
     * - "Maximum Restarts" (\c int): The maximum number of restarts.
     *   This does <i>not</i> include the first "Num Blocks"
     *   iterations (before the first restart).  Its default value is
     *   20.
     *
     * We do not currently perform any sanity checks for these
     * options.  This may affect you if you set some of them but let
     * others keep their default values.  For example, if you set "Num
     * Blocks" to 2 and "Maximum Iterations" to 100, but don't set
     * "Maximum Restarts", you will only get 40 = 20*2 total
     * iterations, rather than 100.  Thus, if you set one of these
     * parameters, you should always set them all.
     *
     * When solving with multiple right-hand sides, the "Block Size"
     * (\c int) parameter controls the number of right-hand sides for
     * which the solver solves at once.  This setting controls both
     * performance and total memory use.  Doubling it (approximately)
     * doubles the total amount of memory used by the solver, but
     * might make the solves faster by reducing synchronization
     * overhead and improving memory bandwidth utilization.  The gain
     * from increasing this tends to level off quickly.  Making this
     * setting too large may actually hurt performance.
     *
     * These options govern convergence and the numerical algorithm:
     * - "Convergence Tolerance" (\c MagnitudeType): The level that
     *   residual norms must reach in order for the solver to stop
     *   iterating.
     * - "Implicit Residual Scaling" (\c std::string): How to scale
     *   the implicit residual norm.  The default is the norm of the
     *   preconditioned initial residual.
     * - "Explicit Residual Scaling" (\c std::string): How to scale
     *   the explicit residual norm.  The default is the norm of the
     *   (unpreconditioned) initial residual.
     * - "Deflation Quorum" (\c int): When solving with multiple
     *   right-hand sides: the number of right-hand sides that must
     *   have converged to the given tolerance, before the solver will
     *   consider all the systems converged.  If -1, then the solver
     *   will require that all the right-hand sides have converged
     *   before declaring all the systems converged.  This must be no
     *   bigger than the "Block Size" parameter.
     * - "Orthogonalization" (\c std::string): The desired
     *   orthogonalization method.  Currently accepted values are
     *   "DGKS", "ICGS", and "IMGS".  Please refer to Belos'
     *   documentation for more details.
     *
     * For an explanation of "implicit" vs. "explicit" residuals,
     * please see the documentation of isLOADetected().  The
     * difference matters if using left preconditioning.  Otherwise,
     * it is not so important to most users.
     *
     * The residual scaling parameters ("Implicit Residual Scaling"
     * and "Explicit Residual Scaling") accept the following values:
     * - "Norm of Initial Residual"
     * - "Norm of Preconditioned Initial Residual"
     * - "Norm of RHS" (RHS stands for "right-hand side")
     * - "None" (no scaling factor)
     *
     * GMRES always uses the 2 norm (square root of sum of squares of
     * magnitudes of entries) to measure convergence.
     *
     * Belos' solvers let users control intermediate "status" output.
     * This output tells you the current iteration and the values of
     * current convergence criteria.  The following parameters control
     * output.  The default values are fine for users who only care
     * about the final result and don't want to see status output.
     * - "Verbosity": a sum of \c MsgType enum values specifying the
     *   verbosity. Default: Belos::Errors.
     * - "Output Frequency" (\c int): How often (in terms of number of
     *   iterations) to print intermediate status output.  The default
     *   (-1) means not to print intermediate status output at all.
     * - "Output Style" (\c OutputType): The style of output.
     *   Accepted values are General and Brief.  Default: General.
     * - "Output Stream" (<tt>Teuchos::RCP<std::ostream></tt>): A
     *   pointer to an output stream to which the solver will write
     *   status output.  The default is a pointer to
     *   <tt>std::cout</tt>.  Currently, if Trilinos was built with
     *   MPI support, only the MPI process with rank 0 in
     *   MPI_COMM_WORLD will print to this output stream.
     * - "Show Maximum Residual Norm Only": When solving for multiple
     *   right-hand sides, this controls whether output shows residual
     *   norms for all the right-hand sides, or just the current
     *   maximum residual norm over all right-hand sides.
     *          \param pl [in] ParameterList with construction information
     *                  \htmlonly
     *                  <iframe src="belos_PseudoBlockGmres.xml" width=100% scrolling="no" frameborder="0">
     *                  </iframe>
     *                  <hr />
     *                  \endhtmlonly
     */
    PseudoBlockGmresSolMgr( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                            const Teuchos::RCP<Teuchos::ParameterList> &pl );

    //! Destructor.
    virtual ~PseudoBlockGmresSolMgr() {};

    //! clone for Inverted Injection (DII)
    Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const override {
      return Teuchos::rcp(new PseudoBlockGmresSolMgr<ScalarType,MV,OP>);
    }
    //@}

    //! @name Accessor methods
    //@{

    const LinearProblem<ScalarType,MV,OP>& getProblem() const override {
      return *problem_;
    }

    //! A list of valid default parameters for this solver.
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

    //! The current parameters for this solver.
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

    //! Iteration count for the most recent call to \c solve().
    int getNumIters() const override {
      return numIters_;
    }

    /// \brief Whether a "loss of accuracy" was detected during the last solve().
    ///
    /// This solver uses two different residual norms to predict
    /// convergence: "implicit" (also called "native") and "explicit"
    /// (also called "exact," not to be confused with "exact
    /// arithmetic").  The "implicit" residuals are computed by the
    /// solver via a recurrence relation (the Arnoldi relation, in the
    /// case of GMRES).  The "explicit" residuals are computed
    /// directly as $B - A X_k$.  Implicit residuals are much cheaper
    /// to compute, since they are available almost "for free" from
    /// the recurrence relation.  In contrast, computing exact
    /// residuals requires computing the current approximate solution
    /// \f$X_k\f$, applying the global operator \f$A\f$ to \f$X_k\f$,
    /// and then computing the norm of the resulting vector(s) via a
    /// global reduction.  Thus, GMRES favors using the cheaper
    /// implicit residuals to predict convergence.  Users typically
    /// want convergence with respect to explicit residuals, though.
    ///
    /// Implicit and explicit residuals may differ due to rounding
    /// error.  However, the difference between implicit and explicit
    /// residuals matters most when using a left (or split)
    /// preconditioner.  In that case, the implicit residuals are
    /// those of the left-preconditioned problem \f$M_L^{-1} A X =
    /// M_L^{-1} B\f$ instead of the original problem \f$A X = B\f$.
    /// The implicit residual norms may thus differ significantly from
    /// the explicit residual norms, even if one could compute without
    /// rounding error.
    ///
    /// When using a left preconditioner, this solver tries to detect
    /// if the implicit residuals have converged but the explicit
    /// residuals have not.  In that case, it will reduce the
    /// convergence tolerance and iterate a little while longer to
    /// attempt to reduce the explicit residual norm.  However, if
    /// that doesn't work, it declares a "loss of accuracy" for the
    /// affected right-hand side(s), and stops iterating on them.
    /// (Not all right-hand sides may have experienced a loss of
    /// accuracy.)  Thus, the affected right-hand sides may or may not
    /// have converged to the desired residual norm tolerance.
    /// Calling this method tells you whether a "loss of accuracy"
    /// (LOA) occurred during the last \c solve() invocation.
    ///
    /// When <i>not</i> using a left preconditioner, this solver will
    /// iterate until both the implicit and explicit residuals
    /// converge.  (It does not start testing the explicit residuals
    /// until the implicit residuals have converged.  This avoids
    /// whenever possible the cost of computing explicit residuals.)
    /// Implicit and explicit residuals may differ due to rounding
    /// error, even though they are identical when no rounding error
    /// occurs.  In this case, the algorithm does <i>not</i> report a
    /// "loss of accuracy," since it continues iterating until the
    /// explicit residuals converge.
    ///
    /// \note Calling \c solve() again resets the flag that reports
    ///   whether a loss of accuracy was detected.  Thus, you should
    ///   call this method immediately after calling \c solve().
    bool isLOADetected() const override { return loaDetected_; }

    //@}

    //! @name Set methods
    //@{

    //! Set the linear problem to solve.
    void setProblem (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem) override {
      problem_ = problem;
    }

    //! Set the parameters the solver manager should use to solve the linear problem.
    void setParameters (const Teuchos::RCP<Teuchos::ParameterList> &params) override;

    /// \brief Set a custom status test.
    ///
    /// A custom status test is not required.  If you decide to set
    /// one, the current implementation will apply it sequentially
    /// (short-circuiting OR, like the || operator in C++) after
    /// Pseudoblock GMRES' standard convergence test.
    virtual void setUserConvStatusTest(
      const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &userConvStatusTest,
      const typename StatusTestCombo<ScalarType,MV,OP>::ComboType &comboType =
          StatusTestCombo<ScalarType,MV,OP>::SEQ
      ) override;

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
     * This method calls PseudoBlockGmresIter::iterate(), which will return either because a specially constructed status test evaluates to
     * ::Passed or an std::exception is thrown.
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

    /// \brief Check current status tests against current linear problem.
    ///
    /// (Re)create all the status tests, based on the current solve
    /// parameters and the current linear problem to solve.  This is
    /// necessary whenever the linear problem is set or changed via \c
    /// setProblem(), because the residual norm test to use depends on
    /// whether or not the (new) linear problem defines a left
    /// preconditioner.  Furthermore, include the user's custom
    /// convergence test if they set one via \c
    /// setUserConvStatusTest().
    ///
    /// \return False if we were able to (re)create all the status
    ///   tests correctly, else true.  The \c solve() routine may call
    ///   this method.  If it does, it checks the return value.
    bool checkStatusTest();

    //! The current linear problem to solve.
    Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;

    // Output manager.
    Teuchos::RCP<OutputManager<ScalarType> > printer_;
    Teuchos::RCP<std::ostream> outputStream_;

    // Status tests.
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > userConvStatusTest_;
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > debugStatusTest_;
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > sTest_;
    Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > maxIterTest_;
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > convTest_;
    Teuchos::RCP<StatusTestResNorm<ScalarType,MV,OP> > impConvTest_, expConvTest_;
    Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputTest_;
    typename StatusTestCombo<ScalarType,MV,OP>::ComboType comboType_;
    Teuchos::RCP<std::map<std::string, Teuchos::RCP<StatusTest<ScalarType, MV, OP> > > > taggedTests_;

    // Orthogonalization manager.
    Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > ortho_;

     // Current parameter list.
    Teuchos::RCP<Teuchos::ParameterList> params_;

    // Default solver values.
    static constexpr int maxRestarts_default_ = 20;
    static constexpr int maxIters_default_ = 1000;
    static constexpr bool showMaxResNormOnly_default_ = false;
    static constexpr int blockSize_default_ = 1;
    static constexpr int numBlocks_default_ = 300;
    static constexpr int verbosity_default_ = Belos::Errors;
    static constexpr int outputStyle_default_ = Belos::General;
    static constexpr int outputFreq_default_ = -1;
    static constexpr int defQuorum_default_ = 1;
    static constexpr const char * impResScale_default_ = "Norm of Preconditioned Initial Residual";
    static constexpr const char * expResScale_default_ = "Norm of Initial Residual";
    static constexpr const char * label_default_ = "Belos";
    static constexpr const char * orthoType_default_ = "DGKS";
    static constexpr std::ostream * outputStream_default_ = &std::cout;

    // Current solver values.
    MagnitudeType convtol_, orthoKappa_, achievedTol_;
    int maxRestarts_, maxIters_, numIters_;
    int blockSize_, numBlocks_, verbosity_, outputStyle_, outputFreq_, defQuorum_;
    bool showMaxResNormOnly_;
    std::string orthoType_;
    std::string impResScale_, expResScale_;
    MagnitudeType resScaleFactor_;

    // Timers.
    std::string label_;
    Teuchos::RCP<Teuchos::Time> timerSolve_;

    // Internal state variables.
    bool isSet_, isSTSet_, expResTest_;
    bool loaDetected_;
  };


// Empty Constructor
template<class ScalarType, class MV, class OP>
PseudoBlockGmresSolMgr<ScalarType,MV,OP>::PseudoBlockGmresSolMgr() :
  outputStream_(Teuchos::rcp(outputStream_default_,false)),
  taggedTests_(Teuchos::null),
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
  defQuorum_(defQuorum_default_),
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  orthoType_(orthoType_default_),
  impResScale_(impResScale_default_),
  expResScale_(expResScale_default_),
  resScaleFactor_(DefaultSolverParameters::resScaleFactor),
  label_(label_default_),
  isSet_(false),
  isSTSet_(false),
  expResTest_(false),
  loaDetected_(false)
{}

// Basic Constructor
template<class ScalarType, class MV, class OP>
PseudoBlockGmresSolMgr<ScalarType,MV,OP>::
PseudoBlockGmresSolMgr (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                        const Teuchos::RCP<Teuchos::ParameterList> &pl) :
  problem_(problem),
  outputStream_(Teuchos::rcp(outputStream_default_,false)),
  taggedTests_(Teuchos::null),
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
  defQuorum_(defQuorum_default_),
  showMaxResNormOnly_(showMaxResNormOnly_default_),
  orthoType_(orthoType_default_),
  impResScale_(impResScale_default_),
  expResScale_(expResScale_default_),
  resScaleFactor_(DefaultSolverParameters::resScaleFactor),
  label_(label_default_),
  isSet_(false),
  isSTSet_(false),
  expResTest_(false),
  loaDetected_(false)
{
  TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null, std::invalid_argument, "Problem not given to solver manager.");

  // If the parameter list pointer is null, then set the current parameters to the default parameter list.
  if (!is_null(pl)) {
    // Set the parameters using the list that was passed in.
    setParameters( pl );
  }
}

template<class ScalarType, class MV, class OP>
void
PseudoBlockGmresSolMgr<ScalarType,MV,OP>::
setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // Create the internal parameter list if one doesn't already exist.
  if (params_ == Teuchos::null) {
    params_ = parameterList (*getValidParameters ());
  } else {
    // TAW: 3/8/2016: do not validate sub parameter lists as they
    //                might not have a pre-defined structure
    //                e.g. user-specified status tests
    // The Belos Pseudo Block GMRES parameters on the first level are
    // not affected and verified.
    params->validateParameters (*getValidParameters (), 0);
  }

  // Check for maximum number of restarts
  if (params->isParameter ("Maximum Restarts")) {
    maxRestarts_ = params->get ("Maximum Restarts", maxRestarts_default_);

    // Update parameter in our list.
    params_->set ("Maximum Restarts", maxRestarts_);
  }

  // Check for maximum number of iterations
  if (params->isParameter ("Maximum Iterations")) {
    maxIters_ = params->get ("Maximum Iterations", maxIters_default_);

    // Update parameter in our list and in status test.
    params_->set ("Maximum Iterations", maxIters_);
    if (! maxIterTest_.is_null ()) {
      maxIterTest_->setMaxIters (maxIters_);
    }
  }

  // Check for blocksize
  if (params->isParameter ("Block Size")) {
    blockSize_ = params->get ("Block Size", blockSize_default_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      blockSize_ <= 0, std::invalid_argument,
      "Belos::PseudoBlockGmresSolMgr::setParameters: "
      "The \"Block Size\" parameter must be strictly positive, "
      "but you specified a value of " << blockSize_ << ".");

    // Update parameter in our list.
    params_->set ("Block Size", blockSize_);
  }

  // Check for the maximum number of blocks.
  if (params->isParameter ("Num Blocks")) {
    numBlocks_ = params->get ("Num Blocks", numBlocks_default_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      numBlocks_ <= 0, std::invalid_argument,
      "Belos::PseudoBlockGmresSolMgr::setParameters: "
      "The \"Num Blocks\" parameter must be strictly positive, "
      "but you specified a value of " << numBlocks_ << ".");

    // Update parameter in our list.
    params_->set ("Num Blocks", numBlocks_);
  }

  // Check to see if the timer label changed.
  if (params->isParameter ("Timer Label")) {
    const std::string tempLabel = params->get ("Timer Label", label_default_);

    // Update parameter in our list and solver timer
    if (tempLabel != label_) {
      label_ = tempLabel;
      params_->set ("Timer Label", label_);
      const std::string solveLabel =
        label_ + ": PseudoBlockGmresSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerSolve_ = Teuchos::TimeMonitor::getNewCounter (solveLabel);
#endif // BELOS_TEUCHOS_TIME_MONITOR
      if (ortho_ != Teuchos::null) {
        ortho_->setLabel( label_ );
      }
    }
  }

  // Check if the orthogonalization changed.
  if (params->isParameter ("Orthogonalization")) {
    std::string tempOrthoType = params->get ("Orthogonalization", orthoType_default_);
#ifdef HAVE_BELOS_TSQR
    TEUCHOS_TEST_FOR_EXCEPTION(
      tempOrthoType != "DGKS" && tempOrthoType != "ICGS" &&
      tempOrthoType != "IMGS" && tempOrthoType != "TSQR",
      std::invalid_argument,
      "Belos::PseudoBlockGmresSolMgr::setParameters: "
      "The \"Orthogonalization\" parameter must be one of \"DGKS\", \"ICGS\", "
      "\"IMGS\", or \"TSQR\".");
#else
    TEUCHOS_TEST_FOR_EXCEPTION(
      tempOrthoType != "DGKS" && tempOrthoType != "ICGS" &&
      tempOrthoType != "IMGS",
      std::invalid_argument,
      "Belos::PseudoBlockGmresSolMgr::setParameters: "
      "The \"Orthogonalization\" parameter must be one of \"DGKS\", \"ICGS\", "
      "or \"IMGS\".");
#endif // HAVE_BELOS_TSQR

    if (tempOrthoType != orthoType_) {
      orthoType_ = tempOrthoType;
      params_->set("Orthogonalization", orthoType_);
      // Create orthogonalization manager
      if (orthoType_ == "DGKS") {
        typedef DGKSOrthoManager<ScalarType, MV, OP> ortho_type;
        if (orthoKappa_ <= 0) {
          ortho_ = rcp (new ortho_type (label_));
        }
        else {
          ortho_ = rcp (new ortho_type (label_));
          rcp_dynamic_cast<ortho_type> (ortho_)->setDepTol (orthoKappa_);
        }
      }
      else if (orthoType_ == "ICGS") {
        typedef ICGSOrthoManager<ScalarType, MV, OP> ortho_type;
        ortho_ = rcp (new ortho_type (label_));
      }
      else if (orthoType_ == "IMGS") {
        typedef IMGSOrthoManager<ScalarType, MV, OP> ortho_type;
        ortho_ = rcp (new ortho_type (label_));
      }
#ifdef HAVE_BELOS_TSQR
      else if (orthoType_ == "TSQR") {
        typedef TsqrMatOrthoManager<ScalarType, MV, OP> ortho_type;
        ortho_ = rcp (new ortho_type (label_));
      }
#endif // HAVE_BELOS_TSQR
    }
  }

  // Check which orthogonalization constant to use.
  if (params->isParameter ("Orthogonalization Constant")) {
    if (params->isType<MagnitudeType> ("Orthogonalization Constant")) {
      orthoKappa_ = params->get ("Orthogonalization Constant",
                                 static_cast<MagnitudeType> (DefaultSolverParameters::orthoKappa));
    }
    else {
      orthoKappa_ = params->get ("Orthogonalization Constant",
                                 DefaultSolverParameters::orthoKappa);
    }

    // Update parameter in our list.
    params_->set ("Orthogonalization Constant", orthoKappa_);
    if (orthoType_ == "DGKS") {
      if (orthoKappa_ > 0 && ! ortho_.is_null ()) {
        typedef DGKSOrthoManager<ScalarType, MV, OP> ortho_type;
        rcp_dynamic_cast<ortho_type> (ortho_)->setDepTol (orthoKappa_);
      }
    }
  }

  // Check for a change in verbosity level
  if (params->isParameter ("Verbosity")) {
    if (Teuchos::isParameterType<int> (*params, "Verbosity")) {
      verbosity_ = params->get ("Verbosity", verbosity_default_);
    } else {
      verbosity_ = (int) Teuchos::getParameter<Belos::MsgType> (*params, "Verbosity");
    }

    // Update parameter in our list.
    params_->set ("Verbosity", verbosity_);
    if (! printer_.is_null ()) {
      printer_->setVerbosity (verbosity_);
    }
  }

  // Check for a change in output style.
  if (params->isParameter ("Output Style")) {
    if (Teuchos::isParameterType<int> (*params, "Output Style")) {
      outputStyle_ = params->get ("Output Style", outputStyle_default_);
    } else {
      outputStyle_ = (int) Teuchos::getParameter<Belos::OutputType> (*params, "Output Style");
    }

    // Update parameter in our list.
    params_->set ("Output Style", verbosity_);
    if (! outputTest_.is_null ()) {
      isSTSet_ = false;
    }

  }

  // Check if user has specified his own status tests
  if (params->isSublist ("User Status Tests")) {
    Teuchos::ParameterList userStatusTestsList = params->sublist("User Status Tests",true);
    if ( userStatusTestsList.numParams() > 0 ) {
      std::string userCombo_string = params->get<std::string>("User Status Tests Combo Type", "SEQ");
      Teuchos::RCP<StatusTestFactory<ScalarType,MV,OP> > testFactory = Teuchos::rcp(new StatusTestFactory<ScalarType,MV,OP>());
      setUserConvStatusTest( testFactory->buildStatusTests(userStatusTestsList), testFactory->stringToComboType(userCombo_string) );
      taggedTests_ = testFactory->getTaggedTests();
      isSTSet_ = false;
    }
  }

  // output stream
  if (params->isParameter ("Output Stream")) {
    outputStream_ = Teuchos::getParameter<Teuchos::RCP<std::ostream> > (*params, "Output Stream");

    // Update parameter in our list.
    params_->set("Output Stream", outputStream_);
    if (! printer_.is_null ()) {
      printer_->setOStream (outputStream_);
    }
  }

  // frequency level
  if (verbosity_ & Belos::StatusTestDetails) {
    if (params->isParameter ("Output Frequency")) {
      outputFreq_ = params->get ("Output Frequency", outputFreq_default_);
    }

    // Update parameter in out list and output status test.
    params_->set ("Output Frequency", outputFreq_);
    if (! outputTest_.is_null ()) {
      outputTest_->setOutputFrequency (outputFreq_);
    }
  }

  // Create output manager if we need to.
  if (printer_.is_null ()) {
    printer_ = rcp (new OutputManager<ScalarType> (verbosity_, outputStream_));
  }

  // Convergence

  // Check for convergence tolerance
  if (params->isParameter ("Convergence Tolerance")) {
    if (params->isType<MagnitudeType> ("Convergence Tolerance")) {
      convtol_ = params->get ("Convergence Tolerance",
                              static_cast<MagnitudeType> (DefaultSolverParameters::convTol));
    }
    else {
      convtol_ = params->get ("Convergence Tolerance", DefaultSolverParameters::convTol);
    }

    // Update parameter in our list and residual tests.
    params_->set ("Convergence Tolerance", convtol_);
    if (! impConvTest_.is_null ()) {
      impConvTest_->setTolerance (convtol_);
    }
    if (! expConvTest_.is_null ()) {
      expConvTest_->setTolerance (convtol_);
    }
  }

  // Grab the user defined residual scaling
  bool userDefinedResidualScalingUpdated = false;
  if (params->isParameter ("User Defined Residual Scaling")) {
    MagnitudeType tempResScaleFactor = DefaultSolverParameters::resScaleFactor;
    if (params->isType<MagnitudeType> ("User Defined Residual Scaling")) {
      tempResScaleFactor = params->get ("User Defined Residual Scaling",
                                        static_cast<MagnitudeType> (DefaultSolverParameters::resScaleFactor));
    }
    else {
      tempResScaleFactor = params->get ("User Defined Residual Scaling",
                                        DefaultSolverParameters::resScaleFactor);
    }

    // Only update the scaling if it's different.
    if (resScaleFactor_ != tempResScaleFactor) {
      resScaleFactor_ = tempResScaleFactor;
      userDefinedResidualScalingUpdated = true;
    }

    if(userDefinedResidualScalingUpdated)
    {
      if (! params->isParameter ("Implicit Residual Scaling") && ! impConvTest_.is_null ()) {
        try {
          if(impResScale_ == "User Provided")
            impConvTest_->defineScaleForm (Belos::UserProvided, Belos::TwoNorm, resScaleFactor_);
        }
        catch (std::exception& e) {
          // Make sure the convergence test gets constructed again.
          isSTSet_ = false;
        }
      }
      if (! params->isParameter ("Explicit Residual Scaling") && ! expConvTest_.is_null ()) {
        try {
          if(expResScale_ == "User Provided")
            expConvTest_->defineScaleForm (Belos::UserProvided, Belos::TwoNorm, resScaleFactor_);
        }
        catch (std::exception& e) {
          // Make sure the convergence test gets constructed again.
          isSTSet_ = false;
        }
      }
    }
  }

  // Check for a change in scaling, if so we need to build new residual tests.
  if (params->isParameter ("Implicit Residual Scaling")) {
    const std::string tempImpResScale =
      Teuchos::getParameter<std::string> (*params, "Implicit Residual Scaling");

    // Only update the scaling if it's different.
    if (impResScale_ != tempImpResScale) {
      Belos::ScaleType impResScaleType = convertStringToScaleType (tempImpResScale);
      impResScale_ = tempImpResScale;

      // Update parameter in our list and residual tests
      params_->set ("Implicit Residual Scaling", impResScale_);
      if (! impConvTest_.is_null ()) {
        try {
          if(impResScale_ == "User Provided")
            impConvTest_->defineScaleForm (impResScaleType, Belos::TwoNorm, resScaleFactor_);
          else
            impConvTest_->defineScaleForm (impResScaleType, Belos::TwoNorm);
        }
        catch (std::exception& e) {
          // Make sure the convergence test gets constructed again.
          isSTSet_ = false;
        }
      }
    }
    else if (userDefinedResidualScalingUpdated) {
      Belos::ScaleType impResScaleType = convertStringToScaleType (impResScale_);

      if (! impConvTest_.is_null ()) {
        try {
          if(impResScale_ == "User Provided")
            impConvTest_->defineScaleForm (impResScaleType, Belos::TwoNorm, resScaleFactor_);
        }
        catch (std::exception& e) {
          // Make sure the convergence test gets constructed again.
          isSTSet_ = false;
        }
      }
    }
  }

  if (params->isParameter ("Explicit Residual Scaling")) {
    const std::string tempExpResScale =
      Teuchos::getParameter<std::string> (*params, "Explicit Residual Scaling");

    // Only update the scaling if it's different.
    if (expResScale_ != tempExpResScale) {
      Belos::ScaleType expResScaleType = convertStringToScaleType (tempExpResScale);
      expResScale_ = tempExpResScale;

      // Update parameter in our list and residual tests
      params_->set ("Explicit Residual Scaling", expResScale_);
      if (! expConvTest_.is_null ()) {
        try {
          if(expResScale_ == "User Provided")
            expConvTest_->defineScaleForm (expResScaleType, Belos::TwoNorm, resScaleFactor_);
          else
            expConvTest_->defineScaleForm (expResScaleType, Belos::TwoNorm);
        }
        catch (std::exception& e) {
          // Make sure the convergence test gets constructed again.
          isSTSet_ = false;
        }
      }
    }
    else if (userDefinedResidualScalingUpdated) {
      Belos::ScaleType expResScaleType = convertStringToScaleType (expResScale_);

      if (! expConvTest_.is_null ()) {
        try {
          if(expResScale_ == "User Provided")
            expConvTest_->defineScaleForm (expResScaleType, Belos::TwoNorm, resScaleFactor_);
        }
        catch (std::exception& e) {
          // Make sure the convergence test gets constructed again.
          isSTSet_ = false;
        }
      }
    }
  }


  if (params->isParameter ("Show Maximum Residual Norm Only")) {
    showMaxResNormOnly_ =
      Teuchos::getParameter<bool> (*params, "Show Maximum Residual Norm Only");

    // Update parameter in our list and residual tests
    params_->set ("Show Maximum Residual Norm Only", showMaxResNormOnly_);
    if (! impConvTest_.is_null ()) {
      impConvTest_->setShowMaxResNormOnly (showMaxResNormOnly_);
    }
    if (! expConvTest_.is_null ()) {
      expConvTest_->setShowMaxResNormOnly (showMaxResNormOnly_);
    }
  }

  // Create status tests if we need to.

  // Get the deflation quorum, or number of converged systems before deflation is allowed
  if (params->isParameter("Deflation Quorum")) {
    defQuorum_ = params->get("Deflation Quorum", defQuorum_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      defQuorum_ > blockSize_, std::invalid_argument,
      "Belos::PseudoBlockGmresSolMgr::setParameters: "
      "The \"Deflation Quorum\" parameter (= " << defQuorum_ << ") must not be "
      "larger than \"Block Size\" (= " << blockSize_ << ").");
    params_->set ("Deflation Quorum", defQuorum_);
    if (! impConvTest_.is_null ()) {
      impConvTest_->setQuorum (defQuorum_);
    }
    if (! expConvTest_.is_null ()) {
      expConvTest_->setQuorum (defQuorum_);
    }
  }

  // Create orthogonalization manager if we need to.
  if (ortho_.is_null ()) {
    params_->set("Orthogonalization", orthoType_);
    if (orthoType_ == "DGKS") {
      typedef DGKSOrthoManager<ScalarType, MV, OP> ortho_type;
      if (orthoKappa_ <= 0) {
        ortho_ = rcp (new ortho_type (label_));
      }
      else {
        ortho_ = rcp (new ortho_type (label_));
        rcp_dynamic_cast<ortho_type> (ortho_)->setDepTol (orthoKappa_);
      }
    }
    else if (orthoType_ == "ICGS") {
      typedef ICGSOrthoManager<ScalarType, MV, OP> ortho_type;
      ortho_ = rcp (new ortho_type (label_));
    }
    else if (orthoType_ == "IMGS") {
      typedef IMGSOrthoManager<ScalarType, MV, OP> ortho_type;
      ortho_ = rcp (new ortho_type (label_));
    }
#ifdef HAVE_BELOS_TSQR
    else if (orthoType_ == "TSQR") {
      typedef TsqrMatOrthoManager<ScalarType, MV, OP> ortho_type;
      ortho_ = rcp (new ortho_type (label_));
    }
#endif // HAVE_BELOS_TSQR
    else {
#ifdef HAVE_BELOS_TSQR
      TEUCHOS_TEST_FOR_EXCEPTION(
        orthoType_ != "ICGS" && orthoType_ != "DGKS" &&
        orthoType_ != "IMGS" && orthoType_ != "TSQR",
        std::logic_error,
        "Belos::PseudoBlockGmresSolMgr::setParameters(): "
        "Invalid orthogonalization type \"" << orthoType_ << "\".");
#else
      TEUCHOS_TEST_FOR_EXCEPTION(
        orthoType_ != "ICGS" && orthoType_ != "DGKS" &&
        orthoType_ != "IMGS",
        std::logic_error,
        "Belos::PseudoBlockGmresSolMgr::setParameters(): "
        "Invalid orthogonalization type \"" << orthoType_ << "\".");
#endif // HAVE_BELOS_TSQR
    }
  }

  // Create the timer if we need to.
  if (timerSolve_ == Teuchos::null) {
    std::string solveLabel = label_ + ": PseudoBlockGmresSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    timerSolve_ = Teuchos::TimeMonitor::getNewCounter (solveLabel);
#endif
  }

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}


template<class ScalarType, class MV, class OP>
void
PseudoBlockGmresSolMgr<ScalarType,MV,OP>::setUserConvStatusTest(
  const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &userConvStatusTest,
  const typename StatusTestCombo<ScalarType,MV,OP>::ComboType &comboType
  )
{
  userConvStatusTest_ = userConvStatusTest;
  comboType_ = comboType;
}

template<class ScalarType, class MV, class OP>
void
PseudoBlockGmresSolMgr<ScalarType,MV,OP>::setDebugStatusTest(
  const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &debugStatusTest
  )
{
  debugStatusTest_ = debugStatusTest;
}



template<class ScalarType, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList>
PseudoBlockGmresSolMgr<ScalarType,MV,OP>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    // Set all the valid parameters and their default values.

    // The static_cast is to resolve an issue with older clang versions which
    // would cause the constexpr to link fail. With c++17 the problem is resolved.
    pl= Teuchos::rcp( new Teuchos::ParameterList() );
    pl->set("Convergence Tolerance", static_cast<MagnitudeType>(DefaultSolverParameters::convTol),
      "The relative residual tolerance that needs to be achieved by the\n"
      "iterative solver in order for the linear system to be declared converged.");
    pl->set("Maximum Restarts", static_cast<int>(maxRestarts_default_),
      "The maximum number of restarts allowed for each\n"
      "set of RHS solved.");
    pl->set("Maximum Iterations", static_cast<int>(maxIters_default_),
      "The maximum number of block iterations allowed for each\n"
      "set of RHS solved.");
    pl->set("Num Blocks", static_cast<int>(numBlocks_default_),
      "The maximum number of vectors allowed in the Krylov subspace\n"
      "for each set of RHS solved.");
    pl->set("Block Size", static_cast<int>(blockSize_default_),
      "The number of RHS solved simultaneously.");
    pl->set("Verbosity", static_cast<int>(verbosity_default_),
      "What type(s) of solver information should be outputted\n"
      "to the output stream.");
    pl->set("Output Style", static_cast<int>(outputStyle_default_),
      "What style is used for the solver information outputted\n"
      "to the output stream.");
    pl->set("Output Frequency", static_cast<int>(outputFreq_default_),
      "How often convergence information should be outputted\n"
      "to the output stream.");
    pl->set("Deflation Quorum", static_cast<int>(defQuorum_default_),
      "The number of linear systems that need to converge before\n"
      "they are deflated.  This number should be <= block size.");
    pl->set("Output Stream", Teuchos::rcp(outputStream_default_,false),
      "A reference-counted pointer to the output stream where all\n"
      "solver output is sent.");
    pl->set("Show Maximum Residual Norm Only", static_cast<bool>(showMaxResNormOnly_default_),
      "When convergence information is printed, only show the maximum\n"
      "relative residual norm when the block size is greater than one.");
    pl->set("Implicit Residual Scaling", static_cast<const char *>(impResScale_default_),
      "The type of scaling used in the implicit residual convergence test.");
    pl->set("Explicit Residual Scaling", static_cast<const char *>(expResScale_default_),
      "The type of scaling used in the explicit residual convergence test.");
    pl->set("Timer Label", static_cast<const char *>(label_default_),
      "The string to use as a prefix for the timer labels.");
    pl->set("Orthogonalization", static_cast<const char *>(orthoType_default_),
      "The type of orthogonalization to use: DGKS, ICGS, IMGS.");
    pl->set("Orthogonalization Constant",static_cast<MagnitudeType>(DefaultSolverParameters::orthoKappa),
      "The constant used by DGKS orthogonalization to determine\n"
      "whether another step of classical Gram-Schmidt is necessary.");
    pl->sublist("User Status Tests");
    pl->set("User Status Tests Combo Type", "SEQ",
        "Type of logical combination operation of user-defined\n"
        "and/or solver-specific status tests.");
    validPL = pl;
  }
  return validPL;
}

// Check the status test versus the defined linear problem
template<class ScalarType, class MV, class OP>
bool PseudoBlockGmresSolMgr<ScalarType,MV,OP>::checkStatusTest() {

  typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
  typedef Belos::StatusTestGenResNorm<ScalarType,MV,OP>  StatusTestGenResNorm_t;
  typedef Belos::StatusTestImpResNorm<ScalarType,MV,OP>  StatusTestImpResNorm_t;

  // Basic test checks maximum iterations and native residual.
  maxIterTest_ = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>( maxIters_ ) );

  // If there is a left preconditioner, we create a combined status test that checks the implicit
  // and then explicit residual norm to see if we have convergence.
  if ( !Teuchos::is_null(problem_->getLeftPrec()) ) {
    expResTest_ = true;
  }

  if (expResTest_) {

    // Implicit residual test, using the native residual to determine if convergence was achieved.
    Teuchos::RCP<StatusTestGenResNorm_t> tmpImpConvTest =
      Teuchos::rcp( new StatusTestGenResNorm_t( convtol_, defQuorum_ ) );
    if(impResScale_ == "User Provided")
      tmpImpConvTest->defineScaleForm( convertStringToScaleType(impResScale_), Belos::TwoNorm, resScaleFactor_ );
    else
      tmpImpConvTest->defineScaleForm( convertStringToScaleType(impResScale_), Belos::TwoNorm );

    tmpImpConvTest->setShowMaxResNormOnly( showMaxResNormOnly_ );
    impConvTest_ = tmpImpConvTest;

    // Explicit residual test once the native residual is below the tolerance
    Teuchos::RCP<StatusTestGenResNorm_t> tmpExpConvTest =
      Teuchos::rcp( new StatusTestGenResNorm_t( convtol_, defQuorum_ ) );
    tmpExpConvTest->defineResForm( StatusTestGenResNorm_t::Explicit, Belos::TwoNorm );
    if(expResScale_ == "User Provided")
      tmpExpConvTest->defineScaleForm( convertStringToScaleType(expResScale_), Belos::TwoNorm, resScaleFactor_ );
    else
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
      Teuchos::rcp( new StatusTestImpResNorm_t( convtol_, defQuorum_ ) );
    if(impResScale_ == "User Provided")
      tmpImpConvTest->defineScaleForm( convertStringToScaleType(impResScale_), Belos::TwoNorm, resScaleFactor_ );
    else
      tmpImpConvTest->defineScaleForm( convertStringToScaleType(impResScale_), Belos::TwoNorm );
    tmpImpConvTest->setShowMaxResNormOnly( showMaxResNormOnly_ );
    impConvTest_ = tmpImpConvTest;

    // Set the explicit and total convergence test to this implicit test that checks for accuracy loss.
    expConvTest_ = impConvTest_;
    convTest_ = impConvTest_;
  }

  if (nonnull(userConvStatusTest_) ) {
    // Override the overall convergence test with the users convergence test
    convTest_ = Teuchos::rcp(
      new StatusTestCombo_t( comboType_, convTest_, userConvStatusTest_ ) );
    // brief output style not compatible with more general combinations
    //outputStyle_ = Belos::General;
    // NOTE: Above, you have to run the other convergence tests also because
    // the logic in this class depends on it.  This is very unfortunate.
  }

  sTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::OR, maxIterTest_, convTest_ ) );

  // Add debug status test, if one is provided by the user
  if (nonnull(debugStatusTest_) ) {
    // Add debug convergence test
    Teuchos::rcp_dynamic_cast<StatusTestCombo_t>(sTest_)->addStatusTest( debugStatusTest_ );
  }

  // Create the status test output class.
  // This class manages and formats the output from the status test.
  StatusTestOutputFactory<ScalarType,MV,OP> stoFactory( outputStyle_, taggedTests_ );
  outputTest_ = stoFactory.create( printer_, sTest_, outputFreq_, Passed+Failed+Undefined );

  // Set the solver string for the output test
  std::string solverDesc = " Pseudo Block Gmres ";
  outputTest_->setSolverDesc( solverDesc );


  // The status test is now set.
  isSTSet_ = true;

  return false;
}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType PseudoBlockGmresSolMgr<ScalarType,MV,OP>::solve() {

  // Set the current parameters if they were not set before.
  // NOTE:  This may occur if the user generated the solver manager with the default constructor and
  // then didn't set any parameters using setParameters().
  if (!isSet_) { setParameters( params_ ); }

  Teuchos::BLAS<int,ScalarType> blas;

  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isProblemSet(),PseudoBlockGmresSolMgrLinearProblemFailure,
                     "Belos::PseudoBlockGmresSolMgr::solve(): Linear problem is not ready, setProblem() has not been called.");

  // Check if we have to create the status tests.
  if (!isSTSet_ || (!expResTest_ && !Teuchos::is_null(problem_->getLeftPrec())) ) {
    TEUCHOS_TEST_FOR_EXCEPTION( checkStatusTest(),PseudoBlockGmresSolMgrLinearProblemFailure,
      "Belos::BlockGmresSolMgr::solve(): Linear problem and requested status tests are incompatible.");
  }

  // Create indices for the linear systems to be solved.
  int startPtr = 0;
  int numRHS2Solve = MVT::GetNumberVecs( *(problem_->getRHS()) );
  int numCurrRHS = ( numRHS2Solve < blockSize_) ? numRHS2Solve : blockSize_;

  std::vector<int> currIdx( numCurrRHS );
  blockSize_ = numCurrRHS;
  for (int i=0; i<numCurrRHS; ++i)
    { currIdx[i] = startPtr+i; }

  // Inform the linear problem of the current linear system to solve.
  problem_->setLSIndex( currIdx );

  //////////////////////////////////////////////////////////////////////////////////////
  // Parameter list
  Teuchos::ParameterList plist;
  plist.set("Num Blocks",numBlocks_);

  // Reset the status test.
  outputTest_->reset();
  loaDetected_ = false;

  // Assume convergence is achieved, then let any failed convergence set this to false.
  bool isConverged = true;

  //////////////////////////////////////////////////////////////////////////////////////
  // BlockGmres solver

  Teuchos::RCP<PseudoBlockGmresIter<ScalarType,MV,OP> > block_gmres_iter
    = Teuchos::rcp( new PseudoBlockGmresIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,ortho_,plist) );

  // Enter solve() iterations
  {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*timerSolve_);
#endif

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
      Teuchos::RCP<MV> tmpV, R_0 = MVT::CloneCopy( *(problem_->getInitPrecResVec()), currIdx );
      newState.V.resize( blockSize_ );
      newState.Z.resize( blockSize_ );
      for (int i=0; i<blockSize_; ++i) {
        index[0]=i;
        tmpV = MVT::CloneViewNonConst( *R_0, index );

        // Get a matrix to hold the orthonormalization coefficients.
        Teuchos::RCP<Teuchos::SerialDenseVector<int,ScalarType> > tmpZ
          = Teuchos::rcp( new Teuchos::SerialDenseVector<int,ScalarType>( 1 ));

        // Orthonormalize the new V_0
        int rank = ortho_->normalize( *tmpV, tmpZ );
        TEUCHOS_TEST_FOR_EXCEPTION(rank != 1, PseudoBlockGmresSolMgrOrthoFailure,
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

            if ( expConvTest_->getLOADetected() ) {
              // Oops!  There was a loss of accuracy (LOA) for one or
              // more right-hand sides.  That means the implicit
              // (a.k.a. "native") residuals claim convergence,
              // whereas the explicit (hence expConvTest_, i.e.,
              // "explicit convergence test") residuals have not
              // converged.
              //
              // We choose to deal with this situation by deflating
              // out the affected right-hand sides and moving on.
              loaDetected_ = true;
              printer_->stream(Warnings) <<
                "Belos::PseudoBlockGmresSolMgr::solve(): Warning! Solver has experienced a loss of accuracy!" << std::endl;
              isConverged = false;
            }

            // Figure out which linear systems converged.
            std::vector<int> convIdx = expConvTest_->convIndices();

            // If the number of converged linear systems is equal to the
            // number of current linear systems, then we are done with this block.
            if (convIdx.size() == currRHSIdx.size())
              break;  // break from while(1){block_gmres_iter->iterate()}

            // Get a new state struct and initialize the solver.
            PseudoBlockGmresIterState<ScalarType,MV> defState;

            // Inform the linear problem that we are finished with this current linear system.
            problem_->setCurrLS();

            // Get the state.
            PseudoBlockGmresIterState<ScalarType,MV> oldState = block_gmres_iter->getState();
            int curDim = oldState.curDim;

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
                defState.V.push_back( Teuchos::rcp_const_cast<MV>( oldState.V[i] ) );
                defState.Z.push_back( Teuchos::rcp_const_cast<Teuchos::SerialDenseVector<int,ScalarType> >( oldState.Z[i] ) );
                defState.H.push_back( Teuchos::rcp_const_cast<Teuchos::SerialDenseMatrix<int,ScalarType> >( oldState.H[i] ) );
                defState.sn.push_back( Teuchos::rcp_const_cast<Teuchos::SerialDenseVector<int,ScalarType> >( oldState.sn[i] ) );
                defState.cs.push_back( Teuchos::rcp_const_cast<Teuchos::SerialDenseVector<int,MagnitudeType> >(oldState.cs[i] ) );
                currRHSIdx[have] = currRHSIdx[i];
                have++;
              }
            }
            defRHSIdx.resize(currRHSIdx.size()-have);
            currRHSIdx.resize(have);

            // Compute the current solution that needs to be deflated if this solver has taken any steps.
            if (curDim) {
              Teuchos::RCP<MV> update = block_gmres_iter->getCurrentUpdate();
              Teuchos::RCP<MV> defUpdate = MVT::CloneViewNonConst( *update, defRHSIdx );

              // Set the deflated indices so we can update the solution.
              problem_->setLSIndex( convIdx );

              // Update the linear problem.
              problem_->updateSolution( defUpdate, true );
            }

            // Set the remaining indices after deflation.
            problem_->setLSIndex( currRHSIdx );

            // Set the dimension of the subspace, which is the same as the old subspace size.
            defState.curDim = curDim;

            // Initialize the solver with the deflated system.
            block_gmres_iter->initialize(defState);
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
            PseudoBlockGmresIterState<ScalarType,MV> oldState = block_gmres_iter->getState();

            // Set the new state.
            PseudoBlockGmresIterState<ScalarType,MV> newstate;
            newstate.V.resize(currRHSIdx.size());
            newstate.Z.resize(currRHSIdx.size());

            // Compute the restart vectors
            // NOTE: Force the linear problem to update the current residual since the solution was updated.
            R_0 = MVT::Clone( *(problem_->getInitPrecResVec()), currRHSIdx.size() );
            problem_->computeCurrPrecResVec( &*R_0 );
            for (unsigned int i=0; i<currRHSIdx.size(); ++i) {
              index[0] = i;  // index(1) vector declared on line 891

              tmpV = MVT::CloneViewNonConst( *R_0, index );

              // Get a matrix to hold the orthonormalization coefficients.
              Teuchos::RCP<Teuchos::SerialDenseVector<int,ScalarType> > tmpZ
                = Teuchos::rcp( new Teuchos::SerialDenseVector<int,ScalarType>( 1 ));

              // Orthonormalize the new V_0
              int rank = ortho_->normalize( *tmpV, tmpZ );
              TEUCHOS_TEST_FOR_EXCEPTION(rank != 1 ,PseudoBlockGmresSolMgrOrthoFailure,
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
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                "Belos::PseudoBlockGmresSolMgr::solve(): Invalid return from PseudoBlockGmresIter::iterate().");
          }
        }
        catch (const PseudoBlockGmresIterOrthoFailure &e) {

          // Try to recover the most recent least-squares solution
          block_gmres_iter->updateLSQR( block_gmres_iter->getCurSubspaceDim() );

          // Check to see if the most recent least-squares solution yielded convergence.
          sTest_->checkStatus( &*block_gmres_iter );
          if (convTest_->getStatus() != Passed)
            isConverged = false;
          break;
        }
        catch (const std::exception &e) {
          printer_->stream(Errors) << "Error! Caught std::exception in PseudoBlockGmresIter::iterate() at iteration "
                                   << block_gmres_iter->getNumIters() << std::endl
                                   << e.what() << std::endl;
          throw;
        }
      }

      // Compute the current solution.
      // Update the linear problem.
      if (nonnull(userConvStatusTest_)) {
        //std::cout << "\nnonnull(userConvStatusTest_)\n";
        Teuchos::RCP<MV> update = block_gmres_iter->getCurrentUpdate();
        problem_->updateSolution( update, true );
      }
      else if (nonnull(expConvTest_->getSolution())) {
        //std::cout << "\nexpConvTest_->getSolution()\n";
        Teuchos::RCP<MV> newX = expConvTest_->getSolution();
        Teuchos::RCP<MV> curX = problem_->getCurrLHSVec();
        MVT::MvAddMv( 0.0, *newX, 1.0, *newX, *curX );
      }
      else {
        //std::cout << "\nblock_gmres_iter->getCurrentUpdate()\n";
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

        blockSize_ = numCurrRHS;
        currIdx.resize( numCurrRHS  );
        for (int i=0; i<numCurrRHS; ++i)
        { currIdx[i] = startPtr+i; }

        // Adapt the status test quorum if we need to.
        if (defQuorum_ > blockSize_) {
          if (impConvTest_ != Teuchos::null)
            impConvTest_->setQuorum( blockSize_ );
          if (expConvTest_ != Teuchos::null)
            expConvTest_->setQuorum( blockSize_ );
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
  // solve.  For this solver, convTest_ may either be a single
  // residual norm test, or a combination of two residual norm tests.
  // In the latter case, the master convergence test convTest_ is a
  // SEQ combo of the implicit resp. explicit tests.  If the implicit
  // test never passes, then the explicit test won't ever be executed.
  // This manifests as expConvTest_->getTestValue()->size() < 1.  We
  // deal with this case by using the values returned by
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
      "Belos::PseudoBlockGmresSolMgr::solve(): The implicit convergence test's "
      "getTestValue() method returned NULL.  Please report this bug to the "
      "Belos developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(pTestValues->size() < 1, std::logic_error,
      "Belos::PseudoBlockGmresSolMgr::solve(): The implicit convergence test's "
      "getTestValue() method returned a vector of length zero.  Please report "
      "this bug to the Belos developers.");

    // FIXME (mfh 12 Dec 2011) Does pTestValues really contain the
    // achieved tolerances for all vectors in the current solve(), or
    // just for the vectors from the last deflation?
    achievedTol_ = *std::max_element (pTestValues->begin(), pTestValues->end());
  }

  if (!isConverged || loaDetected_) {
    return Unconverged; // return from PseudoBlockGmresSolMgr::solve()
  }
  return Converged; // return from PseudoBlockGmresSolMgr::solve()
}


template<class ScalarType, class MV, class OP>
std::string PseudoBlockGmresSolMgr<ScalarType,MV,OP>::description () const
{
  std::ostringstream out;

  out << "\"Belos::PseudoBlockGmresSolMgr\": {";
  if (this->getObjectLabel () != "") {
    out << "Label: " << this->getObjectLabel () << ", ";
  }
  out << "Num Blocks: " << numBlocks_
      << ", Maximum Iterations: " << maxIters_
      << ", Maximum Restarts: " << maxRestarts_
      << ", Convergence Tolerance: " << convtol_
      << "}";
  return out.str ();
}


template<class ScalarType, class MV, class OP>
void
PseudoBlockGmresSolMgr<ScalarType, MV, OP>::
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

    out << "\"Belos::PseudoBlockGmresSolMgr\":" << endl;
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
    out << "Num Blocks: " << numBlocks_ << endl
        << "Maximum Iterations: " << maxIters_ << endl
        << "Maximum Restarts: " << maxRestarts_ << endl
        << "Convergence Tolerance: " << convtol_ << endl;
  }
}

} // end Belos namespace

#endif /* BELOS_PSEUDO_BLOCK_GMRES_SOLMGR_HPP */
