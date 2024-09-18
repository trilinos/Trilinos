// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_RCG_SOLMGR_HPP
#define BELOS_RCG_SOLMGR_HPP

/*! \file BelosRCGSolMgr.hpp
 *  \brief The Belos::RCGSolMgr provides a solver manager for the RCG (Recycling Conjugate Gradient) linear solver.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosRCGIter.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestGenResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutputFactory.hpp"
#include "BelosOutputManager.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_as.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

/** \example epetra/example/RCG/RCGEpetraExFile.cpp
    This is an example of how to use the Belos::RCGSolMgr solver manager in Epetra.
*/
/** \example tpetra/example/RCG/RCGTpetraExFile.cpp
    This is an example of how to use the Belos::RCGSolMgr solver manager in Tpetra.
*/

/*! \class Belos::RCGSolMgr
\brief Implementation of the RCG (Recycling Conjugate Gradient) iterative linear solver.
\ingroup belos_solver_framework
\author Michael Parks and Heidi Thornquist

\section Belos_RCGSol_summary Summary

This class implements the RCG (Recycling CG) iterative linear
solver.  This solver is suited for solving sequences of related linear
systems \f$A_i x_i = b_i\f$, where each matrix \f$A_i\f$ is symmetric
positive definite.

\section Belos_RCG_real ScalarType must be real

This RCG implementation currently only supports real-valued (not
complex-valued) ScalarType types.  You may check whether ScalarType is
complex using the following code:
\code
if (Teuchos::ScalarTraits<ScalarType>::isComplex) {
  // ScalarType is complex valued.
} else {
  // ScalarType is real valued.
}
\endcode

This is not a limitation of the RCG method itself, just of the current
implementation.  If there is sufficient interest, we can remedy this
deficiency.  For now, if you attempt to invoke the constructor when
<tt>ScalarType</tt> is complex, the constructor will throw an
exception.  This is why this class inherits from
Details::RealSolverManager.  RCGSolMgr can still compile if
<tt>ScalarType</tt> is complex, but you will not be able to construct
a RCGSolMgr instance in that case, due to the aforementioned run-time
error that the constructor raises.  We do this so that the class will
still compile, whether ScalarType is real or complex.  This helps make
SolverFactory valid to compile, whether ScalarType is real or complex.

*/

namespace Belos {

  //! @name RCGSolMgr Exceptions
  //@{

  /** \brief RCGSolMgrLinearProblemFailure is thrown when the linear problem is
   * not setup (i.e. setProblem() was not called) when solve() is called.
   *
   * This exception is thrown from the RCGSolMgr::solve() method.
   *
   */
  class RCGSolMgrLinearProblemFailure : public BelosError {public:
    RCGSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief RCGSolMgrLAPACKFailure is thrown when a nonzero value is retuned
   * from an LAPACK call.
   *
   * This exception is thrown from the RCGSolMgr::solve() method.
   *
   */
  class RCGSolMgrLAPACKFailure : public BelosError {public:
    RCGSolMgrLAPACKFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  //@}


  // Partial specialization for unsupported ScalarType types.
  // This contains a stub implementation.
  template<class ScalarType, class MV, class OP, class DM = Teuchos::SerialDenseMatrix<int,ScalarType>,
           const bool supportsScalarType =
             Belos::Details::LapackSupportsScalar<ScalarType>::value &&
             ! Teuchos::ScalarTraits<ScalarType>::isComplex>
  class RCGSolMgr :
    public Details::SolverManagerRequiresRealLapack<ScalarType, MV, OP, DM,
                                                    Belos::Details::LapackSupportsScalar<ScalarType>::value &&
                                                    ! Teuchos::ScalarTraits<ScalarType>::isComplex>
  {
    static const bool scalarTypeIsSupported =
      Belos::Details::LapackSupportsScalar<ScalarType>::value &&
      ! Teuchos::ScalarTraits<ScalarType>::isComplex;
    typedef Details::SolverManagerRequiresRealLapack<ScalarType, MV, OP, DM,
                                                     scalarTypeIsSupported> base_type;

  public:
    RCGSolMgr () :
      base_type ()
    {}
    RCGSolMgr (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
               const Teuchos::RCP<Teuchos::ParameterList> &pl) :
      base_type ()
    {}
    virtual ~RCGSolMgr () {}

    //! clone for Inverted Injection (DII)
    Teuchos::RCP<SolverManager<ScalarType, MV, OP, DM> > clone () const override {
      return Teuchos::rcp(new RCGSolMgr<ScalarType,MV,OP,DM,supportsScalarType>);
    }
  };

  // Partial specialization for real ScalarType.
  // This contains the actual working implementation of RCG.
  // See discussion in the class documentation above.
  template<class ScalarType, class MV, class OP, class DM>
  class RCGSolMgr<ScalarType, MV, OP, DM, true> :
    public Details::SolverManagerRequiresRealLapack<ScalarType, MV, OP, DM, true> {
  private:
    typedef MultiVecTraits<ScalarType,MV,DM> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef DenseMatTraits<ScalarType,DM>    DMT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;

  public:

    //! @name Constructors/Destructor
    //@{

    /*! \brief Empty constructor for RCGSolMgr.
     * This constructor takes no arguments and sets the default values for the solver.
     * The linear problem must be passed in using setProblem() before solve() is called on this object.
     * The solver values can be changed using setParameters().
     */
     RCGSolMgr();

    /*! \brief Basic constructor for RCGSolMgr.
     *
     * This constructor accepts the LinearProblem to be solved in addition
     * to a parameter list of options for the solver manager. These options include the following:
     *   - "Num Blocks" - a \c int specifying length of a cycle (and thus number of max number of search vectors kept). Default: 25
     *   - "Num Recycled Blocks" - a \c int specifying the number of vectors selected for recycling. Default: 3
     *   - "Maximum Iterations" - an \c int specifying the maximum number of iterations the
     *                            underlying solver is allowed to perform. Default: 1000
     *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Belos::Errors
     *   - "Output Style" - a OutputType specifying the style of output. Default: Belos::General
     *   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: 1e-8.
     *   - "Output Stream" - a reference-counted pointer to the output stream where all
     *                       solver output is sent.  Default: Teuchos::rcp(&std::cout,false)
     *   - "Output Frequency" - an \c int specifying how often convergence information should be
     *                          outputted.  Default: -1 (never)
     *   - "Show Maximum Residual Norm Only" - a \c bool specifying whether that only the maximum
     *                                         relative residual norm is printed if convergence
     *                                         information is printed. Default: false
     *   - "Timer Label" - a \c std::string to use as a prefix for the timer labels.  Default: "Belos"
     */

    RCGSolMgr( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                   const Teuchos::RCP<Teuchos::ParameterList> &pl );

    //! Destructor.
    virtual ~RCGSolMgr() {};

    //! clone for Inverted Injection (DII)
    Teuchos::RCP<SolverManager<ScalarType, MV, OP, DM> > clone () const override {
      return Teuchos::rcp(new RCGSolMgr<ScalarType,MV,OP,DM>);
    }
    //@}

    //! @name Accessor methods
    //@{

    const LinearProblem<ScalarType,MV,OP>& getProblem() const override {
      return *problem_;
    }

    /*! \brief Get a parameter list containing the valid parameters for this object. */
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

    /*! \brief Get a parameter list containing the current parameters for this object. */
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
    /// This is set whether or not the solve actually managed to
    /// achieve the desired convergence tolerance.
    MagnitudeType achievedTol() const override {
      return achievedTol_;
    }

    //! Get the iteration count for the most recent call to \c solve().
    int getNumIters() const override {
      return numIters_;
    }

    /*! \brief Return whether a loss of accuracy was detected by this solver during the most current solve. */
    bool isLOADetected() const override { return false; }

    //@}

    //! @name Set methods
    //@{

    //! Set the linear problem that needs to be solved.
    void setProblem( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem ) override { problem_ = problem; }

    //! Set the parameters the solver manager should use to solve the linear problem.
    void setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params ) override;

    //@}

    //! @name Reset method
    //@{
    /*! \brief Performs a reset of the solver manager specified by the \c ResetType.  This informs the
     *  solver manager that the solver should prepare for the next call to solve by resetting certain elements
     *  of the iterative solver strategy. Belos::Problem forces a call to setProblem on the linear problem, and
     *  Belos::RecycleSubspace causes the solver manager to "forget" the recycle space generated by previous calls to the solver.
     *  In the latter case, the next call to solve() will act as if the solver has never been called before.
     */
    void reset( const ResetType type ) override {
      if ((type & Belos::Problem) && !Teuchos::is_null(problem_)) problem_->setProblem();
      else if (type & Belos::RecycleSubspace) existU_ = false;
    }
    //@}

    //! @name Solver application methods
    //@{

    /*! \brief This method performs possibly repeated calls to the underlying linear solver's
     *         iterate() routine until the problem has been solved (as decided by the solver manager)
     *         or the solver manager decides to quit.
     *
     * This method calls RCGIter::iterate(), which will return either because a
     * specially constructed status test evaluates to ::Passed or an std::exception is thrown.
     *
     * A return from RCGIter::iterate() signifies one of the following scenarios:
     *    - the maximum number of iterations has been exceeded. In this scenario, the current solutions
     *      to the linear system will be placed in the linear problem and return ::Unconverged.
     *    - global convergence has been met. In this case, the current solutions to the linear system
     *      will be placed in the linear problem and the solver manager will return ::Converged
     *
     * \returns ::ReturnType specifying:
     *     - ::Converged: the linear problem was solved to the specification required by the solver manager.
     *     - ::Unconverged: the linear problem was not solved to the specification desired by the solver manager.
     */
    ReturnType solve() override;

    //@}

    /** \name Overridden from Teuchos::Describable */
    //@{

    /** \brief Method to return description of the RCG solver manager */
    std::string description() const override;

    //@}

  private:

    // Called by all constructors; Contains init instructions common to all constructors
    void init();

    //  Computes harmonic eigenpairs of projected matrix created during one cycle.
    //  Y contains the harmonic Ritz vectors corresponding to the recycleBlocks eigenvalues of smallest magnitude.
    void getHarmonicVecs(const DM& F,
                         const DM& G,
                         DM& Y);

    // Sort list of n floating-point numbers and return permutation vector
    void sort(std::vector<ScalarType>& dlist, int n, std::vector<int>& iperm);

    // Initialize solver state storage
    void initializeStateStorage();

    // Linear problem.
    Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;

    // Output manager.
    Teuchos::RCP<OutputManager<ScalarType> > printer_;
    Teuchos::RCP<std::ostream> outputStream_;

    // Status test.
    Teuchos::RCP<StatusTest<ScalarType,MV,OP,DM> > sTest_;
    Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP,DM> > maxIterTest_;
    Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP,DM> > convTest_;
    Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP,DM> > outputTest_;

    // Current parameter list.
    Teuchos::RCP<Teuchos::ParameterList> params_;

    // Default solver values.
    static constexpr int maxIters_default_ = 1000;
    static constexpr int blockSize_default_ = 1;
    static constexpr int numBlocks_default_ = 25;
    static constexpr int recycleBlocks_default_ = 3;
    static constexpr bool showMaxResNormOnly_default_ = false;
    static constexpr int verbosity_default_ = Belos::Errors;
    static constexpr int outputStyle_default_ = Belos::General;
    static constexpr int outputFreq_default_ = -1;
    static constexpr const char * label_default_ = "Belos";

    //
    // Current solver values.
    //

    //! Convergence tolerance (read from parameter list).
    MagnitudeType convtol_;

    /// \brief Tolerance achieved by the last \c solve() invocation.
    ///
    /// This is set whether or not the solve actually managed to
    /// achieve the desired convergence tolerance.
    MagnitudeType achievedTol_;

    //! Maximum iteration count (read from parameter list).
    int maxIters_;

    //! Number of iterations taken by the last \c solve() invocation.
    int numIters_;

    int numBlocks_, recycleBlocks_;
    bool showMaxResNormOnly_;
    int verbosity_, outputStyle_, outputFreq_;

    /////////////////////////////////////////////////////////////////////////
    // Solver State Storage
    /////////////////////////////////////////////////////////////////////////
    // Search vectors
    Teuchos::RCP<MV> P_;
    //
    // A times current search direction
    Teuchos::RCP<MV> Ap_;
    //
    // Residual vector
    Teuchos::RCP<MV> r_;
    //
    // Preconditioned residual
    Teuchos::RCP<MV> z_;
    //
    // Flag indicating that the recycle space should be used
    bool existU_;
    //
    // Flag indicating that the updated recycle space has been created
    bool existU1_;
    //
    // Recycled subspace and its image
    Teuchos::RCP<MV> U_, AU_;
    //
    // Recycled subspace for next system and its image
    Teuchos::RCP<MV> U1_;
    //
    // Coefficients arising in RCG iteration
    Teuchos::RCP<std::vector<ScalarType> > Beta_;
    Teuchos::RCP<std::vector<ScalarType> > Alpha_;
    Teuchos::RCP<std::vector<ScalarType> > D_;
    //
    // Solutions to local least-squares problems
    Teuchos::RCP<DM> Delta_;
    //
    // The matrix U^T A U
    Teuchos::RCP<DM> UTAU_;
    //
    // LU factorization of U^T A U
    Teuchos::RCP<DM> LUUTAU_;
    //
    // Data from LU factorization of UTAU
    Teuchos::RCP<std::vector<int> > ipiv_;
    //
    // The matrix (AU)^T AU
    Teuchos::RCP<DM> AUTAU_;
    //
    // The scalar r'*z
    Teuchos::RCP<std::vector<ScalarType> > rTz_old_;
    //
    // Matrices needed for calculation of harmonic Ritz eigenproblem
    Teuchos::RCP<DM> F_,G_,Y_;
    //
    // Matrices needed for updating recycle space
    Teuchos::RCP<DM> L2_,DeltaL2_,AU1TUDeltaL2_;
    Teuchos::RCP<DM> AU1TAU1_, AU1TU1_, AU1TAP_;
    Teuchos::RCP<DM> FY_,GY_;
    Teuchos::RCP<DM> APTAP_;
    Teuchos::RCP<MV> U1Y1_, PY2_;
    Teuchos::RCP<DM> AUTAP_, AU1TU_;
    ScalarType dold;
    /////////////////////////////////////////////////////////////////////////

    // Timers.
    std::string label_;
    Teuchos::RCP<Teuchos::Time> timerSolve_;

    // Internal state variables.
    bool params_Set_;
  };


// Empty Constructor
template<class ScalarType, class MV, class OP, class DM>
RCGSolMgr<ScalarType,MV,OP,DM,true>::RCGSolMgr():
  achievedTol_(0.0),
  numIters_(0)
{
  init();
}

// Basic Constructor
template<class ScalarType, class MV, class OP, class DM>
RCGSolMgr<ScalarType,MV,OP,DM,true>::RCGSolMgr(
                                                     const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                                                     const Teuchos::RCP<Teuchos::ParameterList> &pl ) :
  problem_(problem),
  achievedTol_(0.0),
  numIters_(0)
{
  init();
  TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null, std::invalid_argument, "Problem not given to solver manager.");

  // If the parameter list pointer is null, then set the current parameters to the default parameter list.
  if ( !is_null(pl) ) {
    setParameters( pl );
  }
}

// Common instructions executed in all constructors
template<class ScalarType, class MV, class OP, class DM>
void RCGSolMgr<ScalarType,MV,OP,DM,true>::init()
{
  outputStream_ = Teuchos::rcpFromRef(std::cout);
  convtol_ = DefaultSolverParameters::convTol;
  maxIters_ = maxIters_default_;
  numBlocks_ = numBlocks_default_;
  recycleBlocks_ = recycleBlocks_default_;
  verbosity_ = verbosity_default_;
  outputStyle_= outputStyle_default_;
  outputFreq_= outputFreq_default_;
  showMaxResNormOnly_ = showMaxResNormOnly_default_;
  label_ = label_default_;
  params_Set_ = false;
  P_ = Teuchos::null;
  Ap_ = Teuchos::null;
  r_ = Teuchos::null;
  z_ = Teuchos::null;
  existU_ = false;
  existU1_ = false;
  U_ = Teuchos::null;
  AU_ = Teuchos::null;
  U1_ = Teuchos::null;
  Alpha_ = Teuchos::null;
  Beta_ = Teuchos::null;
  D_ = Teuchos::null;
  Delta_ = Teuchos::null;
  UTAU_ = Teuchos::null;
  LUUTAU_ = Teuchos::null;
  ipiv_ = Teuchos::null;
  AUTAU_ = Teuchos::null;
  rTz_old_ = Teuchos::null;
  F_ = Teuchos::null;
  G_ = Teuchos::null;
  Y_ = Teuchos::null;
  L2_ = Teuchos::null;
  DeltaL2_ = Teuchos::null;
  AU1TUDeltaL2_ = Teuchos::null;
  AU1TAU1_ = Teuchos::null;
  AU1TU1_ = Teuchos::null;
  AU1TAP_ = Teuchos::null;
  FY_ = Teuchos::null;
  GY_ = Teuchos::null;
  APTAP_ = Teuchos::null;
  U1Y1_ = Teuchos::null;
  PY2_ = Teuchos::null;
  AUTAP_ = Teuchos::null;
  AU1TU_ = Teuchos::null;
  dold = 0.;
}

template<class ScalarType, class MV, class OP, class DM>
void RCGSolMgr<ScalarType,MV,OP,DM,true>::setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params )
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

  // Check for the maximum number of blocks.
  if (params->isParameter("Num Blocks")) {
    numBlocks_ = params->get("Num Blocks",numBlocks_default_);
    TEUCHOS_TEST_FOR_EXCEPTION(numBlocks_ <= 0, std::invalid_argument,
                       "Belos::RCGSolMgr: \"Num Blocks\" must be strictly positive.");

    // Update parameter in our list.
    params_->set("Num Blocks", numBlocks_);
  }

  // Check for the maximum number of blocks.
  if (params->isParameter("Num Recycled Blocks")) {
    recycleBlocks_ = params->get("Num Recycled Blocks",recycleBlocks_default_);
    TEUCHOS_TEST_FOR_EXCEPTION(recycleBlocks_ <= 0, std::invalid_argument,
                       "Belos::RCGSolMgr: \"Num Recycled Blocks\" must be strictly positive.");

    TEUCHOS_TEST_FOR_EXCEPTION(recycleBlocks_ >= numBlocks_, std::invalid_argument,
                       "Belos::RCGSolMgr: \"Num Recycled Blocks\" must be less than \"Num Blocks\".");

    // Update parameter in our list.
    params_->set("Num Recycled Blocks", recycleBlocks_);
  }

  // Check to see if the timer label changed.
  if (params->isParameter("Timer Label")) {
    std::string tempLabel = params->get("Timer Label", label_default_);

    // Update parameter in our list and solver timer
    if (tempLabel != label_) {
      label_ = tempLabel;
      params_->set("Timer Label", label_);
      std::string solveLabel = label_ + ": RCGSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerSolve_ = Teuchos::TimeMonitor::getNewCounter(solveLabel);
#endif
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
  typedef Belos::StatusTestCombo<ScalarType,MV,OP,DM>  StatusTestCombo_t;
  typedef Belos::StatusTestGenResNorm<ScalarType,MV,OP,DM>  StatusTestResNorm_t;

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
    maxIterTest_ = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP,DM>( maxIters_ ) );

  // Implicit residual test, using the native residual to determine if convergence was achieved.
  if (convTest_ == Teuchos::null)
    convTest_ = Teuchos::rcp( new StatusTestResNorm_t( convtol_, 1 ) );

  if (sTest_ == Teuchos::null)
    sTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::OR, maxIterTest_, convTest_ ) );

  if (outputTest_ == Teuchos::null) {

    // Create the status test output class.
    // This class manages and formats the output from the status test.
    StatusTestOutputFactory<ScalarType,MV,OP,DM> stoFactory( outputStyle_ );
    outputTest_ = stoFactory.create( printer_, sTest_, outputFreq_, Passed+Failed+Undefined );

    // Set the solver string for the output test
    std::string solverDesc = " Recycling CG ";
    outputTest_->setSolverDesc( solverDesc );
  }

  // Create the timer if we need to.
  if (timerSolve_ == Teuchos::null) {
    std::string solveLabel = label_ + ": RCGSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    timerSolve_ = Teuchos::TimeMonitor::getNewCounter(solveLabel);
#endif
  }

  // Inform the solver manager that the current parameters were set.
  params_Set_ = true;
}


template<class ScalarType, class MV, class OP, class DM>
Teuchos::RCP<const Teuchos::ParameterList>
RCGSolMgr<ScalarType,MV,OP,DM,true>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;

  // Set all the valid parameters and their default values.
  if(is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Convergence Tolerance", static_cast<MagnitudeType>(DefaultSolverParameters::convTol),
      "The relative residual tolerance that needs to be achieved by the\n"
      "iterative solver in order for the linear system to be declared converged.");
    pl->set("Maximum Iterations", static_cast<int>(maxIters_default_),
      "The maximum number of block iterations allowed for each\n"
      "set of RHS solved.");
    pl->set("Block Size", static_cast<int>(blockSize_default_),
      "Block Size Parameter -- currently must be 1 for RCG");
    pl->set("Num Blocks", static_cast<int>(numBlocks_default_),
      "The length of a cycle (and this max number of search vectors kept)\n");
    pl->set("Num Recycled Blocks", static_cast<int>(recycleBlocks_default_),
      "The number of vectors in the recycle subspace.");
    pl->set("Verbosity", static_cast<int>(verbosity_default_),
      "What type(s) of solver information should be outputted\n"
      "to the output stream.");
    pl->set("Output Style", static_cast<int>(outputStyle_default_),
      "What style is used for the solver information outputted\n"
      "to the output stream.");
    pl->set("Output Frequency", static_cast<int>(outputFreq_default_),
      "How often convergence information should be outputted\n"
      "to the output stream.");
    pl->set("Output Stream", Teuchos::rcpFromRef(std::cout),
      "A reference-counted pointer to the output stream where all\n"
      "solver output is sent.");
    pl->set("Show Maximum Residual Norm Only", static_cast<bool>(showMaxResNormOnly_default_),
      "When convergence information is printed, only show the maximum\n"
      "relative residual norm when the block size is greater than one.");
    pl->set("Timer Label", static_cast<const char *>(label_default_),
      "The string to use as a prefix for the timer labels.");
    validPL = pl;
  }
  return validPL;
}

// initializeStateStorage
template<class ScalarType, class MV, class OP, class DM>
void RCGSolMgr<ScalarType,MV,OP,DM,true>::initializeStateStorage() {

    // Check if there is any multivector to clone from.
    Teuchos::RCP<const MV> rhsMV = problem_->getRHS();
    if (rhsMV == Teuchos::null) {
      // Nothing to do
      return;
    }
    else {

      // Initialize the state storage
      TEUCHOS_TEST_FOR_EXCEPTION(static_cast<ptrdiff_t>(numBlocks_) > MVT::GetGlobalLength(*rhsMV),std::invalid_argument,
                         "Belos::RCGSolMgr::initializeStateStorage(): Cannot generate a Krylov basis with dimension larger the operator!");

      // If the subspace has not been initialized before, generate it using the RHS from lp_.
      if (P_ == Teuchos::null) {
        P_ = MVT::Clone( *rhsMV, numBlocks_+2 );
      }
      else {
        // Generate P_ by cloning itself ONLY if more space is needed.
        if (MVT::GetNumberVecs(*P_) < numBlocks_+2) {
          Teuchos::RCP<const MV> tmp = P_;
          P_ = MVT::Clone( *tmp, numBlocks_+2 );
        }
      }

      // Generate Ap_ only if it doesn't exist
      if (Ap_ == Teuchos::null)
        Ap_ = MVT::Clone( *rhsMV, 1 );

      // Generate r_ only if it doesn't exist
      if (r_ == Teuchos::null)
        r_ = MVT::Clone( *rhsMV, 1 );

      // Generate z_ only if it doesn't exist
      if (z_ == Teuchos::null)
        z_ = MVT::Clone( *rhsMV, 1 );

      // If the recycle space has not been initialized before, generate it using the RHS from lp_.
      if (U_ == Teuchos::null) {
        U_ = MVT::Clone( *rhsMV, recycleBlocks_ );
      }
      else {
        // Generate U_ by cloning itself ONLY if more space is needed.
        if (MVT::GetNumberVecs(*U_) < recycleBlocks_) {
          Teuchos::RCP<const MV> tmp = U_;
          U_ = MVT::Clone( *tmp, recycleBlocks_ );
        }
      }

      // If the recycle space has not be initialized before, generate it using the RHS from lp_.
      if (AU_ == Teuchos::null) {
        AU_ = MVT::Clone( *rhsMV, recycleBlocks_ );
      }
      else {
        // Generate AU_ by cloning itself ONLY if more space is needed.
        if (MVT::GetNumberVecs(*AU_) < recycleBlocks_) {
          Teuchos::RCP<const MV> tmp = AU_;
          AU_ = MVT::Clone( *tmp, recycleBlocks_ );
        }
      }

      // If the recycle space has not been initialized before, generate it using the RHS from lp_.
      if (U1_ == Teuchos::null) {
        U1_ = MVT::Clone( *rhsMV, recycleBlocks_ );
      }
      else {
        // Generate U1_ by cloning itself ONLY if more space is needed.
        if (MVT::GetNumberVecs(*U1_) < recycleBlocks_) {
          Teuchos::RCP<const MV> tmp = U1_;
          U1_ = MVT::Clone( *tmp, recycleBlocks_ );
        }
      }

      // Generate Alpha_ only if it doesn't exist, otherwise resize it.
      if (Alpha_ == Teuchos::null)
        Alpha_ = Teuchos::rcp( new std::vector<ScalarType>( numBlocks_, 1 ) );
      else {
        if ( (int)Alpha_->size() != numBlocks_ )
          Alpha_->resize( numBlocks_, 1 );
      }

      // Generate Beta_ only if it doesn't exist, otherwise resize it.
      if (Beta_ == Teuchos::null)
        Beta_ = Teuchos::rcp( new std::vector<ScalarType>( numBlocks_ + 1 ) );
      else {
        if ( ((int)Beta_->size() != (numBlocks_+1)) )
          Beta_->resize( numBlocks_ + 1 );
      }

      // Generate D_ only if it doesn't exist, otherwise resize it.
      if (D_ == Teuchos::null)
        D_ = Teuchos::rcp( new std::vector<ScalarType>( numBlocks_ ) );
      else {
        if ( (int)D_->size() != numBlocks_ ) 
          D_->resize( numBlocks_ );
      }

      // Generate Delta_ only if it doesn't exist, otherwise resize it.
      if (Delta_ == Teuchos::null)
        Delta_ = DMT::Create( recycleBlocks_, numBlocks_ + 1 );
      else {
        if ( (DMT::GetNumRows(*Delta_) != recycleBlocks_) || (DMT::GetNumCols(*Delta_)!= (numBlocks_ + 1)) )
          DMT::Reshape( *Delta_, recycleBlocks_, numBlocks_ + 1 );
      }

      // Generate UTAU_ only if it doesn't exist, otherwise resize it.
      if (UTAU_ == Teuchos::null)
        UTAU_ = DMT::Create( recycleBlocks_, recycleBlocks_ );
      else {
        if ( (DMT::GetNumRows(*UTAU_) != recycleBlocks_) || (DMT::GetNumCols(*UTAU_) != recycleBlocks_) )
          DMT::Reshape( *UTAU_, recycleBlocks_, recycleBlocks_ );
      }

      // Generate LUUTAU_ only if it doesn't exist, otherwise resize it.
      if (LUUTAU_ == Teuchos::null)
        LUUTAU_ = DMT::Create( recycleBlocks_, recycleBlocks_ );
      else {
        if ( (DMT::GetNumRows(*LUUTAU_) != recycleBlocks_) || (DMT::GetNumCols(*LUUTAU_) != recycleBlocks_) )
          DMT::Reshape( *LUUTAU_, recycleBlocks_, recycleBlocks_ );
      }

      // Generate ipiv_ only if it doesn't exist, otherwise resize it.
      if (ipiv_ == Teuchos::null)
        ipiv_ = Teuchos::rcp( new std::vector<int>(recycleBlocks_) );
      else {
        if ( (int)ipiv_->size() != recycleBlocks_ ) // if ipiv not correct size, always resize it
          ipiv_->resize(recycleBlocks_);
      }

      // Generate AUTAU_ only if it doesn't exist, otherwise resize it.
      if (AUTAU_ == Teuchos::null)
        AUTAU_ = DMT::Create( recycleBlocks_, recycleBlocks_ );
      else {
        if ( (DMT::GetNumRows(*AUTAU_) != recycleBlocks_) || (DMT::GetNumCols(*AUTAU_) != recycleBlocks_) )
          DMT::Reshape( *AUTAU_, recycleBlocks_, recycleBlocks_ );
      }

      // Generate rTz_old_ only if it doesn't exist
      if (rTz_old_ == Teuchos::null)
        rTz_old_ = Teuchos::rcp( new std::vector<ScalarType>(1) );
      else {
        if ( (rTz_old_->size() != 1) )
          rTz_old_->resize( 1 );
      }

      // Generate F_ only if it doesn't exist
      if (F_ == Teuchos::null)
        F_ = DMT::Create( numBlocks_+recycleBlocks_, numBlocks_+recycleBlocks_ );
      else {
        if ( (DMT::GetNumRows(*F_) != (numBlocks_+recycleBlocks_)) || (DMT::GetNumCols(*F_) != numBlocks_+recycleBlocks_) )
          DMT::Reshape( *F_, numBlocks_+recycleBlocks_, numBlocks_+recycleBlocks_ );
      }

      // Generate G_ only if it doesn't exist
      if (G_ == Teuchos::null)
        G_ = DMT::Create( numBlocks_+recycleBlocks_, numBlocks_+recycleBlocks_ );
      else {
        if ( (DMT::GetNumRows(*G_) != (numBlocks_+recycleBlocks_)) || (DMT::GetNumCols(*G_) != numBlocks_+recycleBlocks_) )
          DMT::Reshape( *G_, numBlocks_+recycleBlocks_, numBlocks_+recycleBlocks_ );
      }

      // Generate Y_ only if it doesn't exist
      if (Y_ == Teuchos::null)
        Y_ = DMT::Create( numBlocks_+recycleBlocks_, recycleBlocks_ );
      else {
        if ( (DMT::GetNumRows(*Y_) != (numBlocks_+recycleBlocks_)) || (DMT::GetNumCols(*Y_) != recycleBlocks_) )
          DMT::Reshape( *Y_, numBlocks_+recycleBlocks_, recycleBlocks_ );
      }

      // Generate L2_ only if it doesn't exist
      if (L2_ == Teuchos::null)
        L2_ = DMT::Create( numBlocks_+1, numBlocks_ );
      else {
        if ( (DMT::GetNumRows(*L2_) != (numBlocks_+1)) || (DMT::GetNumCols(*L2_) != numBlocks_) )
          DMT::Reshape( *L2_, numBlocks_+1, numBlocks_ );
      }

      // Generate DeltaL2_ only if it doesn't exist
      if (DeltaL2_ == Teuchos::null)
        DeltaL2_ = DMT::Create( recycleBlocks_, numBlocks_ );
      else {
        if ( (DMT::GetNumRows(*DeltaL2_) != recycleBlocks_) || (DMT::GetNumCols(*DeltaL2_) != numBlocks_ ) )
          DMT::Reshape( *DeltaL2_, recycleBlocks_, numBlocks_ );
      }

      // Generate AU1TUDeltaL2_ only if it doesn't exist
      if (AU1TUDeltaL2_ == Teuchos::null)
        AU1TUDeltaL2_ = DMT::Create( recycleBlocks_, numBlocks_ );
      else {
        if ( (DMT::GetNumRows(*AU1TUDeltaL2_) != recycleBlocks_) || (DMT::GetNumCols(*AU1TUDeltaL2_) != numBlocks_ ) )
          DMT::Reshape( *AU1TUDeltaL2_, recycleBlocks_, numBlocks_ );
      }

      // Generate AU1TAU1_ only if it doesn't exist
      if (AU1TAU1_ == Teuchos::null)
        AU1TAU1_ = DMT::Create( recycleBlocks_, recycleBlocks_ );
      else {
        if ( (DMT::GetNumRows(*AU1TAU1_) != recycleBlocks_) || (DMT::GetNumCols(*AU1TAU1_) != recycleBlocks_ ) )
          DMT::Reshape( *AU1TAU1_, recycleBlocks_, recycleBlocks_ );
      }

      // Generate GY_ only if it doesn't exist
      if (GY_ == Teuchos::null)
        GY_ = DMT::Create( numBlocks_ + recycleBlocks_, recycleBlocks_ );
      else {
        if ( (DMT::GetNumRows(*GY_) != (numBlocks_ + recycleBlocks_)) || (DMT::GetNumCols(*GY_) != recycleBlocks_ ) )
          DMT::Reshape( *GY_, numBlocks_+recycleBlocks_, recycleBlocks_ );
      }

      // Generate AU1TU1_ only if it doesn't exist
      if (AU1TU1_ == Teuchos::null)
        AU1TU1_ = DMT::Create( recycleBlocks_, recycleBlocks_ );
      else {
        if ( (DMT::GetNumRows(*AU1TU1_) != recycleBlocks_) || (DMT::GetNumCols(*AU1TU1_) != recycleBlocks_ ) )
          DMT::Reshape( *AU1TU1_, recycleBlocks_, recycleBlocks_ );
      }

      // Generate FY_ only if it doesn't exist
      if (FY_ == Teuchos::null)
        FY_ = DMT::Create( numBlocks_ + recycleBlocks_, recycleBlocks_ );
      else {
        if ( (DMT::GetNumRows(*FY_) != (numBlocks_ + recycleBlocks_)) || (DMT::GetNumCols(*FY_) != recycleBlocks_ ) )
          DMT::Reshape( *FY_, numBlocks_+recycleBlocks_, recycleBlocks_ );
      }

      // Generate AU1TAP_ only if it doesn't exist
      if (AU1TAP_ == Teuchos::null)
        AU1TAP_ = DMT::Create( recycleBlocks_, numBlocks_ );
      else {
        if ( (DMT::GetNumRows(*AU1TAP_) != recycleBlocks_) || (DMT::GetNumCols(*AU1TAP_) != numBlocks_ ) )
          DMT::Reshape( *AU1TAP_, recycleBlocks_, numBlocks_ );
      }

      // Generate APTAP_ only if it doesn't exist
      if (APTAP_ == Teuchos::null)
        APTAP_ = DMT::Create( numBlocks_, numBlocks_ );
      else {
        if ( (DMT::GetNumRows(*APTAP_) != numBlocks_) || (DMT::GetNumCols(*APTAP_) != (numBlocks_) ) )
          DMT::Reshape( *APTAP_, numBlocks_, numBlocks_ );
      }

      // If the subspace has not been initialized before, generate it using the RHS from lp_.
      if (U1Y1_ == Teuchos::null) {
        U1Y1_ = MVT::Clone( *rhsMV, recycleBlocks_ );
      }
      else {
        // Generate U1Y1_ by cloning itself ONLY if more space is needed.
        if (MVT::GetNumberVecs(*U1Y1_) < recycleBlocks_) {
          Teuchos::RCP<const MV> tmp = U1Y1_;
          U1Y1_ = MVT::Clone( *tmp, recycleBlocks_ );
        }
      }

      // If the subspace has not been initialized before, generate it using the RHS from lp_.
      if (PY2_ == Teuchos::null) {
        PY2_ = MVT::Clone( *rhsMV, recycleBlocks_ );
      }
      else {
        // Generate PY2_ by cloning itself ONLY if more space is needed.
        if (MVT::GetNumberVecs(*PY2_) < recycleBlocks_) {
          Teuchos::RCP<const MV> tmp = PY2_;
          PY2_ = MVT::Clone( *tmp, recycleBlocks_ );
        }
      }

      // Generate AUTAP_ only if it doesn't exist
      if (AUTAP_ == Teuchos::null)
        AUTAP_ = DMT::Create( recycleBlocks_, numBlocks_ );
      else {
        if ( (DMT::GetNumRows(*AUTAP_) != recycleBlocks_) || (DMT::GetNumCols(*AUTAP_) != numBlocks_ ) )
          DMT::Reshape( *AUTAP_, recycleBlocks_, numBlocks_ );
      }

      // Generate AU1TU_ only if it doesn't exist
      if (AU1TU_ == Teuchos::null)
        AU1TU_ = DMT::Create( recycleBlocks_, recycleBlocks_ );
      else {
        if ( (DMT::GetNumRows(*AU1TU_) != recycleBlocks_) || (DMT::GetNumCols(*AU1TU_) != recycleBlocks_ ) )
          DMT::Reshape( *AU1TU_, recycleBlocks_, recycleBlocks_ );
      }


    }
}

template<class ScalarType, class MV, class OP, class DM>
ReturnType RCGSolMgr<ScalarType,MV,OP,DM,true>::solve() {

  Teuchos::BLAS<int,ScalarType> blas;
  Teuchos::LAPACK<int,ScalarType> lapack;
  ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

  std::vector<int> index(1), rindex(recycleBlocks_), nindex(numBlocks_);
  for (int i=0; i<recycleBlocks_; ++i) { rindex[i] = i; }
  for (int i=0; i<numBlocks_; ++i) { nindex[i] = i; }

  // Count of number of cycles performed on current rhs
  int cycle = 0;

  // Set the current parameters if they were not set before.
  // NOTE:  This may occur if the user generated the solver manager with the default constructor and
  // then didn't set any parameters using setParameters().
  if (!params_Set_) {
    setParameters(Teuchos::parameterList(*getValidParameters()));
  }

  TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null,RCGSolMgrLinearProblemFailure,
                     "Belos::RCGSolMgr::solve(): Linear problem is not a valid object.");
  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isProblemSet(),RCGSolMgrLinearProblemFailure,
                     "Belos::RCGSolMgr::solve(): Linear problem is not ready, setProblem() has not been called.");
  TEUCHOS_TEST_FOR_EXCEPTION((problem_->getLeftPrec() != Teuchos::null)&&(problem_->getRightPrec() != Teuchos::null),
                     RCGSolMgrLinearProblemFailure,
                     "Belos::RCGSolMgr::solve(): RCG does not support split preconditioning, only set left or right preconditioner.");

  // Grab the preconditioning object
  Teuchos::RCP<OP> precObj;
  if (problem_->getLeftPrec() != Teuchos::null) {
    precObj = Teuchos::rcp_const_cast<OP>(problem_->getLeftPrec());
  }
  else if (problem_->getRightPrec() != Teuchos::null) {
    precObj = Teuchos::rcp_const_cast<OP>(problem_->getRightPrec());
  }

  // Create indices for the linear systems to be solved.
  int numRHS2Solve = MVT::GetNumberVecs( *(problem_->getRHS()) );
  std::vector<int> currIdx(1);
  currIdx[0] = 0;

  // Inform the linear problem of the current linear system to solve.
  problem_->setLSIndex( currIdx );

  // Check the number of blocks and change them if necessary.
  ptrdiff_t dim = MVT::GetGlobalLength( *(problem_->getRHS()) );
  if (numBlocks_ > dim) {
    numBlocks_ = Teuchos::asSafe<int>(dim);
    params_->set("Num Blocks", numBlocks_);
    printer_->stream(Warnings) <<
      "Warning! Requested Krylov subspace dimension is larger than operator dimension!" << std::endl <<
      " The maximum number of blocks allowed for the Krylov subspace will be adjusted to " << numBlocks_ << std::endl;
  }

  // Initialize storage for all state variables
  initializeStateStorage();

  // Parameter list
  Teuchos::ParameterList plist;
  plist.set("Num Blocks",numBlocks_);
  plist.set("Recycled Blocks",recycleBlocks_);

  // Reset the status test.
  outputTest_->reset();

  // Assume convergence is achieved, then let any failed convergence set this to false.
  bool isConverged = true;

  // Compute AU = A*U, UTAU = U'*AU, AUTAU = (AU)'*(AU)
  if (existU_) {
    Teuchos::RCP<const MV> Utmp  = MVT::CloneView( *U_,  rindex );
    Teuchos::RCP<MV> AUtmp = MVT::CloneViewNonConst( *AU_, rindex );
    // Initialize AU
    problem_->applyOp( *Utmp, *AUtmp );
    // Initialize UTAU
    MVT::MvTransMv( one, *Utmp, *AUtmp, *UTAU_ );
    // Initialize AUTAU  ( AUTAU = AU'*(M\AU) )
    if ( precObj != Teuchos::null ) {
      Teuchos::RCP<MV> PCAU = MVT::CloneViewNonConst( *U1_, rindex ); // use U1 as temp storage
      OPT::Apply( *precObj, *AUtmp, *PCAU );
      MVT::MvTransMv( one, *AUtmp, *PCAU, *AUTAU_ );
    } else {
      MVT::MvTransMv( one, *AUtmp, *AUtmp, *AUTAU_ );
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // RCG solver

  Teuchos::RCP<RCGIter<ScalarType,MV,OP,DM> > rcg_iter;
  rcg_iter = Teuchos::rcp( new RCGIter<ScalarType,MV,OP,DM>(problem_,printer_,outputTest_,plist) );

  // Enter solve() iterations
  {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*timerSolve_);
#endif

    while ( numRHS2Solve > 0 ) {

      // Debugging output to tell use if recycle space exists and will be used
      if (printer_->isVerbosity( Debug ) ) {
        if (existU_) printer_->print( Debug, "Using recycle space generated from previous call to solve()." );
        else printer_->print( Debug, "No recycle space exists." );
      }

      // Reset the number of iterations.
      rcg_iter->resetNumIters();

      // Set the current number of recycle blocks and subspace dimension with the RCG iteration.
      rcg_iter->setSize( recycleBlocks_, numBlocks_ );

      // Reset the number of calls that the status test output knows about.
      outputTest_->resetNumCalls();

      // indicate that updated recycle space has not yet been generated for this linear system
      existU1_ = false;

      // reset cycle count
      cycle = 0;

      // Get the current residual
      problem_->computeCurrResVec( &*r_ );

      // If U exists, find best soln over this space first
      if (existU_) {
        // Solve linear system UTAU * y = (U'*r)
        Teuchos::RCP<DM> Utr = DMT::Create(recycleBlocks_,1);
        Teuchos::RCP<const MV> Utmp  = MVT::CloneView( *U_,  rindex );
        MVT::MvTransMv( one, *Utmp, *r_, *Utr );
         
        DMT::SyncHostToDevice(*LUUTAU_);
        DMT::Assign(*LUUTAU_,*UTAU_);
        DMT::SyncDeviceToHost( *LUUTAU_ );
        DMT::SyncDeviceToHost( *Utr );
        int info = 0;
        lapack.GESV(recycleBlocks_, 1, DMT::GetRawHostPtr(*LUUTAU_), DMT::GetStride(*LUUTAU_), 
                    &(*ipiv_)[0], DMT::GetRawHostPtr(*Utr), DMT::GetStride(*Utr), &info);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, RCGSolMgrLAPACKFailure,
                           "Belos::RCGSolMgr::solve(): LAPACK GESV failed to compute a solution.");
        DMT::SyncHostToDevice( *Utr );

        // Update solution (x = x + U*y)
        MVT::MvTimesMatAddMv( one, *Utmp, *Utr, one, *problem_->getCurrLHSVec() );

        // Update residual ( r = r - AU*y )
        Teuchos::RCP<const MV> AUtmp  = MVT::CloneView( *AU_, rindex );
        MVT::MvTimesMatAddMv( -one, *AUtmp, *Utr, one, *r_ );
      }

      if ( precObj != Teuchos::null ) {
        OPT::Apply( *precObj, *r_, *z_ );
      } else {
        z_ = r_;
      }

      // rTz_old = r'*z
      MVT::MvDot( *r_, *z_, *rTz_old_ );

      if ( existU_ ) {
        // mu = UTAU\(AU'*z);
        Teuchos::RCP<DM> mu = DMT::Subview(*Delta_, recycleBlocks_, 1);
        Teuchos::RCP<const MV> AUtmp  = MVT::CloneView( *AU_, rindex );
        MVT::MvTransMv( one, *AUtmp, *z_, *mu );

        DMT::SyncDeviceToHost( *Delta_ );
        char TRANS = 'N';
        int info;
        lapack.GETRS( TRANS, recycleBlocks_, 1, DMT::GetConstRawHostPtr(*LUUTAU_), DMT::GetStride(*LUUTAU_), 
                      &(*ipiv_)[0], DMT::GetRawHostPtr(*mu), DMT::GetStride(*mu), &info );
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, RCGSolMgrLAPACKFailure,
                           "Belos::RCGSolMgr::solve(): LAPACK GETRS failed to compute a solution.");
        DMT::SyncHostToDevice( *mu );

        // p  = z - U*mu;
        index.resize( 1 );
        index[0] = 0;
        Teuchos::RCP<MV> Ptmp  = MVT::CloneViewNonConst( *P_, index );
        MVT::Assign(*z_,*Ptmp);
        MVT::MvTimesMatAddMv( -one, *U_, *mu, one, *Ptmp );
      } else {
        // p = z;
        index.resize( 1 );
        index[0] = 0;
        Teuchos::RCP<MV> Ptmp  = MVT::CloneViewNonConst( *P_, index );
        MVT::Assign(*z_,*Ptmp);
      }

      // Set the new state and initialize the solver.
      RCGIterState<ScalarType,MV,DM> newstate;

      // Create RCP views here
      index.resize( numBlocks_+1 );
      for (int ii=0; ii<(numBlocks_+1); ++ii) { index[ii] = ii; }
      newstate.P  = MVT::CloneViewNonConst( *P_,  index );
      newstate.U  = MVT::CloneViewNonConst( *U_,  rindex );
      newstate.AU  = MVT::CloneViewNonConst( *AU_,  rindex );
      newstate.Alpha = Alpha_;
      newstate.Beta = Beta_;
      newstate.Beta_i = 0;
      newstate.D = D_;
      newstate.Delta = Delta_;
      newstate.LUUTAU = LUUTAU_;
      // assign the rest of the values in the struct
      newstate.curDim = 1; // We have initialized the first search vector
      newstate.Ap = Ap_;
      newstate.r = r_;
      newstate.z = z_;
      newstate.existU = existU_;
      newstate.ipiv = ipiv_;
      newstate.rTz_old = rTz_old_;

      rcg_iter->initialize(newstate);

      while(1) {

        // tell rcg_iter to iterate
        try {
          rcg_iter->iterate();

          ////////////////////////////////////////////////////////////////////////////////////
          //
          // check convergence first
          //
          ////////////////////////////////////////////////////////////////////////////////////
          if ( convTest_->getStatus() == Passed ) {
	    // We have convergence
            break; // break from while(1){rcg_iter->iterate()}
          }
          ////////////////////////////////////////////////////////////////////////////////////
          //
          // check for maximum iterations
          //
          ////////////////////////////////////////////////////////////////////////////////////
          else if ( maxIterTest_->getStatus() == Passed ) {
            // we don't have convergence
            isConverged = false;
            break; // break from while(1){rcg_iter->iterate()}
          }
          ////////////////////////////////////////////////////////////////////////////////////
          //
          // check if cycle complete; update for next cycle
          //
          ////////////////////////////////////////////////////////////////////////////////////
          else if ( rcg_iter->getCurSubspaceDim() == rcg_iter->getMaxSubspaceDim() ) {
            // index into P_ of last search vector generated this cycle
            int lastp = -1;
            // index into Beta_ of last entry generated this cycle
            int lastBeta = -1;
            if (recycleBlocks_ > 0) {
              if (!existU_) {
                if (cycle == 0) { // No U, no U1

                   Teuchos::RCP<DM> Ftmp = DMT::Subview( *F_, numBlocks_, numBlocks_ );
                   Teuchos::RCP<DM> Gtmp = DMT::Subview( *G_, numBlocks_, numBlocks_ );
                   DMT::PutScalar( *Ftmp, zero );
                   DMT::PutScalar( *Gtmp, zero );
                   DMT::SyncDeviceToHost( *F_ );
                   DMT::SyncDeviceToHost( *G_ );
		   for (int ii=0;ii<numBlocks_;ii++) {
                     DMT::Value(*Gtmp,ii,ii) = ((*D_)[ii] / (*Alpha_)[ii])*(1 + (*Beta_)[ii]);
                     if (ii > 0) {
                       DMT::Value(*Gtmp,ii-1,ii) = -(*D_)[ii]/(*Alpha_)[ii-1];
                       DMT::Value(*Gtmp,ii,ii-1) = -(*D_)[ii]/(*Alpha_)[ii-1];
                     }
                     DMT::Value(*Ftmp,ii,ii) = (*D_)[ii];
                   }
                   DMT::SyncHostToDevice( *F_ );
                   DMT::SyncHostToDevice( *G_ );

                   // compute harmonic Ritz vectors
                   DMT::SyncDeviceToHost( *Y_ );
                   Teuchos::RCP<DM> Ytmp = DMT::Subview( *Y_, numBlocks_, recycleBlocks_ );
                   getHarmonicVecs(*Ftmp,*Gtmp,*Ytmp);
                   DMT::SyncHostToDevice( *Y_ );

                   // U1 = [P(:,1:end-1)*Y];
                   Teuchos::RCP<const MV> Ptmp = MVT::CloneView( *P_,  nindex );
                   Teuchos::RCP<MV> U1tmp  = MVT::CloneViewNonConst( *U1_,  rindex );
                   MVT::MvTimesMatAddMv( one, *Ptmp, *Ytmp, zero, *U1tmp );

                   // Precompute some variables for next cycle
                   DMT::SyncDeviceToHost(*GY_);
                   DMT::SyncDeviceToHost(*AU1TAU1_);
                   DMT::SyncDeviceToHost(*FY_);
                   DMT::SyncDeviceToHost(*AU1TU1_);

                   // AU1TAU1     = Y'*G*Y;
                   Teuchos::RCP<DM> GYtmp = DMT::Subview( *GY_, numBlocks_, recycleBlocks_ );
                   //GYtmp->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,*Gtmp,*Ytmp,zero);
                   blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, numBlocks_, recycleBlocks_, numBlocks_, 
                              one, DMT::GetConstRawHostPtr(*Gtmp), DMT::GetStride(*Gtmp),
                              DMT::GetConstRawHostPtr(*Ytmp), DMT::GetStride(*Ytmp),
                              zero, DMT::GetRawHostPtr(*GYtmp), DMT::GetStride(*GYtmp));
                   //AU1TAU1_->multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,*Ytmp,*GYtmp,zero);
                   blas.GEMM( Teuchos::TRANS, Teuchos::NO_TRANS, recycleBlocks_, recycleBlocks_, numBlocks_,
                              one, DMT::GetConstRawHostPtr(*Ytmp), DMT::GetStride(*Ytmp),
                              DMT::GetConstRawHostPtr(*GYtmp), DMT::GetStride(*GYtmp),
                              zero, DMT::GetRawHostPtr(*AU1TAU1_), DMT::GetStride(*AU1TAU1_));


                   // AU1TU1      = Y'*F*Y;
                   Teuchos::RCP<DM> FYtmp = DMT::Subview( *FY_, numBlocks_, recycleBlocks_ );
                   //FYtmp->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,*Ftmp,*Ytmp,zero);
                   blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, numBlocks_, recycleBlocks_, numBlocks_, 
                              one, DMT::GetConstRawHostPtr(*Ftmp), DMT::GetStride(*Ftmp),
                              DMT::GetConstRawHostPtr(*Ytmp), DMT::GetStride(*Ytmp),
                              zero, DMT::GetRawHostPtr(*FYtmp), DMT::GetStride(*FYtmp));
                   //AU1TU1_->multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,*Ytmp,*FYtmp,zero);
                   blas.GEMM( Teuchos::TRANS, Teuchos::NO_TRANS, recycleBlocks_, recycleBlocks_, numBlocks_,
                              one, DMT::GetConstRawHostPtr(*Ytmp), DMT::GetStride(*Ytmp),
                              DMT::GetConstRawHostPtr(*FYtmp), DMT::GetStride(*FYtmp),
                              zero, DMT::GetRawHostPtr(*AU1TU1_), DMT::GetStride(*AU1TU1_));

                   DMT::SyncHostToDevice(*AU1TAU1_);
                   DMT::SyncHostToDevice(*AU1TU1_);
                   DMT::SyncHostToDevice(*AU1TAP_);
                   
                   Teuchos::RCP<DM> AU1TAPtmp = DMT::Subview( *AU1TAP_, recycleBlocks_, 1 );
                   // Must reinitialize AU1TAP; can become dense later
                   DMT::PutScalar( *AU1TAPtmp, zero );
                   // AU1TAP(:,1) = Y(end,:)' * (-1/Alpha(end));
		   DMT::SyncDeviceToHost( *AU1TAP_ );
		   ScalarType alphatmp = -1.0 / (*Alpha_)[numBlocks_-1];
                   for (int ii=0; ii<recycleBlocks_; ++ii) {
                      DMT::Value(*AU1TAPtmp,ii,0) = DMT::ValueConst(*Ytmp,numBlocks_-1,ii) * alphatmp;
                   }
                   DMT::SyncHostToDevice(*AU1TAP_);

                   // indicate that updated recycle space now defined
                   existU1_ = true;

                   // Indicate the size of the P, Beta structures generated this cycle
                   lastp = numBlocks_;
                   lastBeta = numBlocks_-1;

                } // if (cycle == 0)
                else { // No U, but U1 guaranteed to exist now

                   // Finish computation of subblocks
                   // AU1TAP = AU1TAP * D(1);
                   DMT::Scale(*AU1TAP_,(*D_)[0]);

                   DMT::PutScalar(*APTAP_,zero);
                   DMT::SyncDeviceToHost(*APTAP_);
                   for (int ii=0; ii<numBlocks_; ii++) {
                     DMT::Value(*APTAP_,ii,ii) = ((*D_)[ii] / (*Alpha_)[ii])*(1 + (*Beta_)[ii+1]);
                     if (ii > 0) {
                       DMT::Value(*APTAP_,ii-1,ii) = -(*D_)[ii]/(*Alpha_)[ii-1];
                       DMT::Value(*APTAP_,ii,ii-1) = -(*D_)[ii]/(*Alpha_)[ii-1];
                     }
                   }
                   DMT::SyncHostToDevice(*APTAP_);

                   // F = [AU1TU1 zeros(k,m); zeros(m,k) diag(D)];
                   DMT::PutScalar(*F_,zero);
                   Teuchos::RCP<DM> F11 = DMT::Subview( *F_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::RCP<DM> F22 = DMT::Subview( *F_, numBlocks_, numBlocks_, recycleBlocks_, recycleBlocks_ );
                   DMT::Assign(*F11,*AU1TU1_);
                   DMT::SyncDeviceToHost(*F_);
                   for(int ii=0;ii<numBlocks_;ii++) {
                     DMT::Value(*F22,ii,ii) = (*D_)[ii];
                   }
                   DMT::SyncHostToDevice(*F_);

                   // G = [AU1TAU1 AU1TAP; AU1TAP' APTAP];
                   Teuchos::RCP<DM> G11 = DMT::Subview( *G_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::RCP<DM> G12 = DMT::Subview( *G_, recycleBlocks_, numBlocks_, 0, recycleBlocks_ );
                   Teuchos::RCP<DM> G21 = DMT::Subview( *G_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   Teuchos::RCP<DM> G22 = DMT::Subview( *G_, numBlocks_, numBlocks_, recycleBlocks_, recycleBlocks_ );
                   DMT::Assign(*G11, *AU1TAU1_);
                   DMT::Assign(*G12, *AU1TAP_);
                   DMT::Assign(*G22, *APTAP_);
                   DMT::SyncDeviceToHost( *G_ );
                   // G21 = G12'; (no transpose operator exists for DM; Do copy manually)
                   for (int ii=0;ii<recycleBlocks_;++ii)
                     for (int jj=0;jj<numBlocks_;++jj)
                       DMT::Value(*G21,jj,ii) = DMT::ValueConst(*G12,ii,jj);
                   DMT::SyncHostToDevice( *G_ );

                   // compute harmonic Ritz vectors
                   getHarmonicVecs(*F_,*G_,*Y_);
                   DMT::SyncHostToDevice( *Y_ );

                   // U1 = [U1 P(:,2:end-1)]*Y;
                   index.resize( numBlocks_ );
                   for (int ii=0; ii<numBlocks_; ++ii) { index[ii] = ii+1; }
                   Teuchos::RCP<const MV> Ptmp = MVT::CloneView( *P_,  index );
                   Teuchos::RCP<MV> PY2tmp = MVT::CloneViewNonConst( *PY2_,  rindex );
                   Teuchos::RCP<MV> U1tmp  = MVT::CloneViewNonConst( *U1_,  rindex );
                   Teuchos::RCP<MV> U1Y1tmp  = MVT::CloneViewNonConst( *U1Y1_,  rindex );
                   Teuchos::RCP<const DM> Y1 = DMT::SubviewConst( *Y_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::RCP<const DM> Y2 = DMT::SubviewConst( *Y_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   MVT::MvTimesMatAddMv( one, *Ptmp, *Y2, zero, *PY2tmp );
                   MVT::MvTimesMatAddMv( one, *U1tmp, *Y1, zero, *U1Y1tmp );
                   MVT::MvAddMv(one,*U1Y1tmp, one, *PY2tmp, *U1tmp);

                   // Precompute some variables for next cycle
                   DMT::SyncDeviceToHost(*GY_);
                   DMT::SyncDeviceToHost(*AU1TAU1_);

                   // AU1TAU1     = Y'*G*Y;
                   //GY_->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,*G_,*Y_,zero);
                   blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, numBlocks_+recycleBlocks_, 
                              recycleBlocks_, numBlocks_+recycleBlocks_, 
                              one, DMT::GetConstRawHostPtr(*G_), DMT::GetStride(*G_),
                              DMT::GetConstRawHostPtr(*Y_), DMT::GetStride(*Y_),
                              zero, DMT::GetRawHostPtr(*GY_), DMT::GetStride(*GY_));
                   //AU1TAU1_->multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,*Y_,*GY_,zero);
                   blas.GEMM( Teuchos::TRANS, Teuchos::NO_TRANS, recycleBlocks_, recycleBlocks_, 
                              numBlocks_+recycleBlocks_,
                              one, DMT::GetConstRawHostPtr(*Y_), DMT::GetStride(*Y_),
                              DMT::GetConstRawHostPtr(*GY_), DMT::GetStride(*GY_),
                              zero, DMT::GetRawHostPtr(*AU1TAU1_), DMT::GetStride(*AU1TAU1_));

                   DMT::SyncHostToDevice(*GY_);
                   DMT::SyncHostToDevice(*AU1TAU1_);

                   // AU1TAP      = zeros(k,m);
                   // AU1TAP(:,1) = Y(end,:)' * (-1/Alpha(end));
                   DMT::PutScalar(*AU1TAP_,zero);
                   ScalarType alphatmp = -1.0 / (*Alpha_)[numBlocks_-1];
                   DMT::SyncDeviceToHost(*AU1TAP_);
                   for (int ii=0; ii<recycleBlocks_; ++ii) {
                      DMT::Value(*AU1TAP_,ii,0) = DMT::ValueConst(*Y_,numBlocks_+recycleBlocks_-1,ii) * alphatmp;
                   }
                   DMT::SyncHostToDevice(*AU1TAP_);

                   DMT::SyncDeviceToHost(*FY_);
                   DMT::SyncDeviceToHost(*AU1TU1_);

                   // AU1TU1      = Y'*F*Y;
                   //FY_->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,*F_,*Y_,zero);
                   blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, numBlocks_+recycleBlocks_, 
                              recycleBlocks_, numBlocks_+recycleBlocks_, 
                              one, DMT::GetConstRawHostPtr(*F_), DMT::GetStride(*F_),
                              DMT::GetConstRawHostPtr(*Y_), DMT::GetStride(*Y_),
                              zero, DMT::GetRawHostPtr(*FY_), DMT::GetStride(*FY_));
                   //AU1TU1_->multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,*Y_,*FY_,zero);
                   blas.GEMM( Teuchos::TRANS, Teuchos::NO_TRANS, recycleBlocks_, recycleBlocks_, 
                              numBlocks_+recycleBlocks_,
                              one, DMT::GetConstRawHostPtr(*Y_), DMT::GetStride(*Y_),
                              DMT::GetConstRawHostPtr(*FY_), DMT::GetStride(*FY_),
                              zero, DMT::GetRawHostPtr(*AU1TU1_), DMT::GetStride(*AU1TU1_));

                   DMT::SyncHostToDevice(*FY_);
                   DMT::SyncHostToDevice(*AU1TU1_);

                   // Indicate the size of the P, Beta structures generated this cycle
                   lastp = numBlocks_+1;
                   lastBeta = numBlocks_;

                } // if (cycle != 1)
              } // if (!existU_)
              else { // U exists
                if (cycle == 0) { // No U1, but U exists
                   DMT::PutScalar(*APTAP_,zero);
                   DMT::SyncDeviceToHost(*APTAP_);
                   for (int ii=0; ii<numBlocks_; ii++) {
                     DMT::Value(*APTAP_,ii,ii) = ((*D_)[ii] / (*Alpha_)[ii])*(1 + (*Beta_)[ii]);
                     if (ii > 0) {
                       DMT::Value(*APTAP_,ii-1,ii) = -(*D_)[ii]/(*Alpha_)[ii-1];
                       DMT::Value(*APTAP_,ii,ii-1) = -(*D_)[ii]/(*Alpha_)[ii-1];
                     }
                   }
                   DMT::SyncHostToDevice(*APTAP_);

                   DMT::PutScalar(*L2_,zero);
                   DMT::SyncDeviceToHost(*L2_);
                   for(int ii=0;ii<numBlocks_;ii++) {
                     DMT::Value(*L2_,ii,ii)   = 1./(*Alpha_)[ii];
                     DMT::Value(*L2_,ii+1,ii) = -1./(*Alpha_)[ii];
                   }
                   DMT::SyncHostToDevice(*L2_);

                   // AUTAP = UTAU*Delta*L2;
                   DMT::SyncDeviceToHost(*Delta_);
		   DMT::SyncDeviceToHost(*DeltaL2_);
                   DMT::SyncDeviceToHost(*AUTAP_);
                   DMT::SyncDeviceToHost(*UTAU_);

                   //DeltaL2_->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,*Delta_,*L2_,zero);
                   blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, recycleBlocks_, numBlocks_, numBlocks_+1,
                              one, DMT::GetConstRawHostPtr(*Delta_), DMT::GetStride(*Delta_),
                              DMT::GetConstRawHostPtr(*L2_), DMT::GetStride(*L2_),
                              zero, DMT::GetRawHostPtr(*DeltaL2_), DMT::GetStride(*DeltaL2_));
                   //AUTAP_->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,*UTAU_,*DeltaL2_,zero);
                   blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, recycleBlocks_, numBlocks_, recycleBlocks_, 
                              one, DMT::GetConstRawHostPtr(*UTAU_), DMT::GetStride(*UTAU_),
                              DMT::GetConstRawHostPtr(*DeltaL2_), DMT::GetStride(*DeltaL2_),
                              zero, DMT::GetRawHostPtr(*AUTAP_), DMT::GetStride(*AUTAP_));

                   DMT::SyncHostToDevice(*DeltaL2_);
                   DMT::SyncHostToDevice(*AUTAP_);

                   // F = [UTAU zeros(k,m); zeros(m,k) diag(D)];
                   DMT::PutScalar(*F_,zero);
                   Teuchos::RCP<DM> F11 = DMT::Subview( *F_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::RCP<DM> F22 = DMT::Subview( *F_, numBlocks_, numBlocks_, recycleBlocks_, recycleBlocks_ );
                   DMT::Assign(*F11,*UTAU_);
                   DMT::SyncDeviceToHost(*F_);
                   for(int ii=0;ii<numBlocks_;ii++) {
                     DMT::Value(*F22,ii,ii) = (*D_)[ii];
                   }
                   DMT::SyncHostToDevice(*F_);

                   // G = [AUTAU AUTAP; AUTAP' APTAP];
                   Teuchos::RCP<DM> G11 = DMT::Subview( *G_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::RCP<DM> G12 = DMT::Subview( *G_, recycleBlocks_, numBlocks_, 0, recycleBlocks_ );
                   Teuchos::RCP<DM> G21 = DMT::Subview( *G_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   Teuchos::RCP<DM> G22 = DMT::Subview( *G_, numBlocks_, numBlocks_, recycleBlocks_, recycleBlocks_ );
                   DMT::Assign(*G11,*AUTAU_);
                   DMT::Assign(*G12,*AUTAP_);
                   DMT::Assign(*G22,*APTAP_);
                   DMT::SyncDeviceToHost(*G_);
                   // G21 = G12'; (no transpose operator exists for DM; Do copy manually)
                   for (int ii=0;ii<recycleBlocks_;++ii)
                     for (int jj=0;jj<numBlocks_;++jj)
                       DMT::Value(*G21,jj,ii) = DMT::ValueConst(*G12,ii,jj);
                   DMT::SyncHostToDevice(*G_);

                   // compute harmonic Ritz vectors
                   getHarmonicVecs(*F_,*G_,*Y_);
                   DMT::SyncHostToDevice(*Y_);

                   // U1 = [U P(:,1:end-1)]*Y;
                   Teuchos::RCP<const MV> Utmp = MVT::CloneView( *U_,  rindex );
                   Teuchos::RCP<const MV> Ptmp = MVT::CloneView( *P_,  nindex );
                   Teuchos::RCP<MV> PY2tmp = MVT::CloneViewNonConst( *PY2_,  rindex );
                   Teuchos::RCP<MV> UY1tmp  = MVT::CloneViewNonConst( *U1Y1_,  rindex );
                   Teuchos::RCP<MV> U1tmp  = MVT::CloneViewNonConst( *U1_,  rindex );
                   Teuchos::RCP<const DM> Y1 = DMT::SubviewConst( *Y_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::RCP<const DM> Y2 = DMT::SubviewConst( *Y_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   MVT::MvTimesMatAddMv( one, *Ptmp, *Y2, zero, *PY2tmp );
                   MVT::MvTimesMatAddMv( one, *Utmp, *Y1, zero, *UY1tmp );
                   MVT::MvAddMv(one,*UY1tmp, one, *PY2tmp, *U1tmp);

                   // Precompute some variables for next cycle
                   DMT::SyncDeviceToHost(*GY_);
                   DMT::SyncDeviceToHost(*AU1TAU1_);
                   DMT::SyncDeviceToHost(*FY_);
                   DMT::SyncDeviceToHost(*AU1TU1_);

                   // AU1TAU1     = Y'*G*Y;
                   //GY_->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,*G_,*Y_,zero);
                   blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, numBlocks_+recycleBlocks_, 
                              recycleBlocks_, numBlocks_+recycleBlocks_, 
                              one, DMT::GetConstRawHostPtr(*G_), DMT::GetStride(*G_),
                              DMT::GetConstRawHostPtr(*Y_), DMT::GetStride(*Y_),
                              zero, DMT::GetRawHostPtr(*GY_), DMT::GetStride(*GY_));
                   //AU1TAU1_->multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,*Y_,*GY_,zero);
                   blas.GEMM( Teuchos::TRANS, Teuchos::NO_TRANS, recycleBlocks_, recycleBlocks_, 
                              numBlocks_+recycleBlocks_,
                              one, DMT::GetConstRawHostPtr(*Y_), DMT::GetStride(*Y_),
                              DMT::GetConstRawHostPtr(*GY_), DMT::GetStride(*GY_),
                              zero, DMT::GetRawHostPtr(*AU1TAU1_), DMT::GetStride(*AU1TAU1_));

                   // AU1TU1      = Y'*F*Y;
                   //FY_->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,*F_,*Y_,zero);
                   blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, numBlocks_+recycleBlocks_, 
                              recycleBlocks_, numBlocks_+recycleBlocks_, 
                              one, DMT::GetConstRawHostPtr(*F_), DMT::GetStride(*F_),
                              DMT::GetConstRawHostPtr(*Y_), DMT::GetStride(*Y_),
                              zero, DMT::GetRawHostPtr(*FY_), DMT::GetStride(*FY_));
                   //AU1TU1_->multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,*Y_,*FY_,zero);
                   blas.GEMM( Teuchos::TRANS, Teuchos::NO_TRANS, recycleBlocks_, recycleBlocks_, 
                              numBlocks_+recycleBlocks_,
                              one, DMT::GetConstRawHostPtr(*Y_), DMT::GetStride(*Y_),
                              DMT::GetConstRawHostPtr(*FY_), DMT::GetStride(*FY_),
                              zero, DMT::GetRawHostPtr(*AU1TU1_), DMT::GetStride(*AU1TU1_));

                   DMT::SyncHostToDevice(*GY_);
                   DMT::SyncHostToDevice(*AU1TAU1_);
                   DMT::SyncHostToDevice(*FY_);
                   DMT::SyncHostToDevice(*AU1TU1_);

                   // AU1TU   = UTAU;
                   DMT::Assign(*AU1TU_,*UTAU_);

                   // dold    = D(end);
                   dold = (*D_)[numBlocks_-1];

                   // indicate that updated recycle space now defined
                   existU1_ = true;

                   // Indicate the size of the P, Beta structures generated this cycle
                   lastp = numBlocks_;
                   lastBeta = numBlocks_-1;
                }
                else { // Have U and U1

                   DMT::SyncDeviceToHost(*APTAP_);
                   for (int ii=0; ii<numBlocks_; ii++) {
                     DMT::Value(*APTAP_,ii,ii) = ((*D_)[ii] / (*Alpha_)[ii])*(1 + (*Beta_)[ii+1]);
                     if (ii > 0) {
                       DMT::Value(*APTAP_,ii-1,ii) = -(*D_)[ii]/(*Alpha_)[ii-1];
                       DMT::Value(*APTAP_,ii,ii-1) = -(*D_)[ii]/(*Alpha_)[ii-1];
                     }
                   }
                   DMT::SyncHostToDevice(*APTAP_);

                   DMT::SyncDeviceToHost(*L2_);
                   for(int ii=0;ii<numBlocks_;ii++) {
                     DMT::Value(*L2_,ii,ii)   = 1./(*Alpha_)[ii];
                     DMT::Value(*L2_,ii+1,ii) = -1./(*Alpha_)[ii];
                   }
                   DMT::SyncHostToDevice(*L2_);

                   DMT::SyncDeviceToHost(*Delta_);
                   DMT::SyncDeviceToHost(*DeltaL2_);
                   DMT::SyncDeviceToHost(*AU1TUDeltaL2_);
                   DMT::SyncDeviceToHost(*AU1TAP_);

                   // M(end,1) = dold*(-Beta(1)/Alpha(1));
                   // AU1TAP = Y'*[AU1TU*Delta*L2; M];
                   //DeltaL2_->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,*Delta_,*L2_,zero);
                   blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, recycleBlocks_, numBlocks_, numBlocks_+1,
                              one, DMT::GetConstRawHostPtr(*Delta_), DMT::GetStride(*Delta_),
                              DMT::GetConstRawHostPtr(*L2_), DMT::GetStride(*L2_),
                              zero, DMT::GetRawHostPtr(*DeltaL2_), DMT::GetStride(*DeltaL2_));
                   //AU1TUDeltaL2_->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,*AU1TU_,*DeltaL2_,zero);
                   blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, recycleBlocks_, numBlocks_, recycleBlocks_, 
                              one, DMT::GetConstRawHostPtr(*AU1TU_), DMT::GetStride(*AU1TU_),
                              DMT::GetConstRawHostPtr(*DeltaL2_), DMT::GetStride(*DeltaL2_),
                              zero, DMT::GetRawHostPtr(*AU1TUDeltaL2_), DMT::GetStride(*AU1TUDeltaL2_));

		   DMT::SyncDeviceToHost( *Y_);
                   Teuchos::RCP<const DM> Y1 = DMT::SubviewConst( *Y_, recycleBlocks_, recycleBlocks_ );
		   Teuchos::RCP<const DM> Y2 = DMT::SubviewConst( *Y_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );

		   //AU1TAP_->multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,*Y1,*AU1TUDeltaL2_,zero);
                   blas.GEMM( Teuchos::TRANS, Teuchos::NO_TRANS, recycleBlocks_, numBlocks_, recycleBlocks_, 
                              one, DMT::GetConstRawHostPtr(*Y1), DMT::GetStride(*Y1),
                              DMT::GetConstRawHostPtr(*AU1TUDeltaL2_), DMT::GetStride(*AU1TUDeltaL2_),
                              zero, DMT::GetRawHostPtr(*AU1TAP_), DMT::GetStride(*AU1TAP_));
                   ScalarType val = dold * (-(*Beta_)[0]/(*Alpha_)[0]);
                   for(int ii=0;ii<recycleBlocks_;ii++) {
                     DMT::Value(*AU1TAP_,ii,0) += DMT::ValueConst(*Y2,numBlocks_-1,ii)*val;
                   }
                   DMT::SyncHostToDevice(*AU1TAP_);

                   // AU1TU = Y1'*AU1TU
                   Teuchos::RCP<DM> Y1TAU1TU = DMT::Subview( *GY_, recycleBlocks_, recycleBlocks_ );
                   //Y1TAU1TU->multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,*Y1,*AU1TU_,zero);
                   blas.GEMM( Teuchos::TRANS, Teuchos::NO_TRANS, recycleBlocks_, recycleBlocks_, recycleBlocks_, 
                              one, DMT::GetConstRawHostPtr(*Y1), DMT::GetStride(*Y1),
                              DMT::GetConstRawHostPtr(*AU1TU_), DMT::GetStride(*AU1TU_),
                              zero, DMT::GetRawHostPtr(*Y1TAU1TU), DMT::GetStride(*Y1TAU1TU));
                   DMT::SyncHostToDevice(*GY_);
                   DMT::Assign(*AU1TU_,*Y1TAU1TU);

                   // F = [AU1TU1 zeros(k,m); zeros(m,k) diag(D)];
                   DMT::PutScalar(*F_,zero);
                   Teuchos::RCP<DM> F11 = DMT::Subview( *F_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::RCP<DM> F22 = DMT::Subview( *F_, numBlocks_, numBlocks_, recycleBlocks_, recycleBlocks_ );
                   DMT::Assign(*F11,*AU1TU1_);
                   DMT::SyncDeviceToHost(*F_);
                   for(int ii=0;ii<numBlocks_;ii++) {
                     DMT::Value(*F22,ii,ii) = (*D_)[ii];
                   }
                   DMT::SyncHostToDevice(*F_);

                   // G = [AU1TAU1 AU1TAP; AU1TAP' APTAP];
                   Teuchos::RCP<DM> G11 = DMT::Subview( *G_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::RCP<DM> G12 = DMT::Subview( *G_, recycleBlocks_, numBlocks_, 0, recycleBlocks_ );
                   Teuchos::RCP<DM> G21 = DMT::Subview( *G_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   Teuchos::RCP<DM> G22 = DMT::Subview( *G_, numBlocks_, numBlocks_, recycleBlocks_, recycleBlocks_ );
                   DMT::Assign(*G11,*AU1TAU1_);
                   DMT::Assign(*G12,*AU1TAP_);
                   DMT::Assign(*G22,*APTAP_);
                   DMT::SyncDeviceToHost(*G_);
                   // G21 = G12'; (no transpose operator exists for DM; Do copy manually)
                   for (int ii=0;ii<recycleBlocks_;++ii)
                     for (int jj=0;jj<numBlocks_;++jj)
                       DMT::Value(*G21,jj,ii) = DMT::ValueConst(*G12,ii,jj);
                   DMT::SyncHostToDevice(*G_);

                   // compute harmonic Ritz vectors
                   getHarmonicVecs(*F_,*G_,*Y_);
                   DMT::SyncHostToDevice(*Y_);

                   // U1 = [U1 P(:,2:end-1)]*Y;
                   index.resize( numBlocks_ );
                   for (int ii=0; ii<numBlocks_; ++ii) { index[ii] = ii+1; }
                   Teuchos::RCP<const MV> Ptmp = MVT::CloneView( *P_,  index );
                   Teuchos::RCP<MV> PY2tmp = MVT::CloneViewNonConst( *PY2_,  rindex );
                   Teuchos::RCP<MV> U1tmp  = MVT::CloneViewNonConst( *U1_,  rindex );
                   Teuchos::RCP<MV> U1Y1tmp  = MVT::CloneViewNonConst( *U1Y1_,  rindex );
                   MVT::MvTimesMatAddMv( one, *Ptmp, *Y2, zero, *PY2tmp );
                   MVT::MvTimesMatAddMv( one, *U1tmp, *Y1, zero, *U1Y1tmp );
                   MVT::MvAddMv(one,*U1Y1tmp, one, *PY2tmp, *U1tmp);

                   // Precompute some variables for next cycle
                   DMT::SyncDeviceToHost(*GY_);
                   DMT::SyncDeviceToHost(*AU1TAU1_);
                   DMT::SyncDeviceToHost(*FY_);
                   DMT::SyncDeviceToHost(*AU1TU1_);

                   // AU1TAU1     = Y'*G*Y;
                   //GY_->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,*G_,*Y_,zero);
                   blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, numBlocks_+recycleBlocks_, 
                              recycleBlocks_, numBlocks_+recycleBlocks_, 
                              one, DMT::GetConstRawHostPtr(*G_), DMT::GetStride(*G_),
                              DMT::GetConstRawHostPtr(*Y_), DMT::GetStride(*Y_),
                              zero, DMT::GetRawHostPtr(*GY_), DMT::GetStride(*GY_));
                   //AU1TAU1_->multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,*Y_,*GY_,zero);
                   blas.GEMM( Teuchos::TRANS, Teuchos::NO_TRANS, recycleBlocks_, recycleBlocks_, 
                              numBlocks_+recycleBlocks_,
                              one, DMT::GetConstRawHostPtr(*Y_), DMT::GetStride(*Y_),
                              DMT::GetConstRawHostPtr(*GY_), DMT::GetStride(*GY_),
                              zero, DMT::GetRawHostPtr(*AU1TAU1_), DMT::GetStride(*AU1TAU1_));

                   // AU1TU1      = Y'*F*Y;
                   //FY_->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,*F_,*Y_,zero);
                   blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, numBlocks_+recycleBlocks_, 
                              recycleBlocks_, numBlocks_+recycleBlocks_, 
                              one, DMT::GetConstRawHostPtr(*F_), DMT::GetStride(*F_),
                              DMT::GetConstRawHostPtr(*Y_), DMT::GetStride(*Y_),
                              zero, DMT::GetRawHostPtr(*FY_), DMT::GetStride(*FY_));
                   //AU1TU1_->multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,*Y_,*FY_,zero);
                   blas.GEMM( Teuchos::TRANS, Teuchos::NO_TRANS, recycleBlocks_, recycleBlocks_, 
                              numBlocks_+recycleBlocks_,
                              one, DMT::GetConstRawHostPtr(*Y_), DMT::GetStride(*Y_),
                              DMT::GetConstRawHostPtr(*FY_), DMT::GetStride(*FY_),
                              zero, DMT::GetRawHostPtr(*AU1TU1_), DMT::GetStride(*AU1TU1_));

                   DMT::SyncHostToDevice(*GY_);
                   DMT::SyncHostToDevice(*AU1TAU1_);
                   DMT::SyncHostToDevice(*FY_);
                   DMT::SyncHostToDevice(*AU1TU1_);

                   // dold    = D(end);
                   dold = (*D_)[numBlocks_-1];

                   // Indicate the size of the P, Beta structures generated this cycle
                   lastp = numBlocks_+1;
                   lastBeta = numBlocks_;

                }
              }
            } // if (recycleBlocks_ > 0)

            // Cleanup after end of cycle

            // P = P(:,end-1:end);
            index.resize( 2 );
            index[0] = lastp-1; index[1] = lastp;
            Teuchos::RCP<const MV> Ptmp2 = MVT::CloneView( *P_,  index );
            index[0] = 0;       index[1] = 1;
            MVT::SetBlock(*Ptmp2,index,*P_);

            // Beta = Beta(end);
            (*Beta_)[0] = (*Beta_)[lastBeta];

            // Delta = Delta(:,end);
            if (existU_) { // Delta only initialized if U exists
              Teuchos::RCP<DM> mu1 = DMT::Subview( *Delta_, recycleBlocks_, 1, 0, 0 );
              Teuchos::RCP<DM> mu2 = DMT::Subview( *Delta_, recycleBlocks_, 1, 0, numBlocks_ );
              DMT::Assign(*mu1,*mu2);
            }

            // Now reinitialize state variables for next cycle
            newstate.P = Teuchos::null;
            index.resize( numBlocks_+1 );
            for (int ii=0; ii<(numBlocks_+1); ++ii) { index[ii] = ii+1; }
            newstate.P  = MVT::CloneViewNonConst( *P_,  index );

            newstate.Beta = Beta_;
            newstate.Beta_i = 1;

            newstate.Delta = Delta_;

            newstate.curDim = 1; // We have initialized the first search vector

            // Pass to iteration object
            rcg_iter->initialize(newstate);

            // increment cycle count
            cycle = cycle + 1;

          }
          ////////////////////////////////////////////////////////////////////////////////////
          //
          // we returned from iterate(), but none of our status tests Passed.
          // something is wrong, and it is probably our fault.
          //
          ////////////////////////////////////////////////////////////////////////////////////
          else {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                               "Belos::RCGSolMgr::solve(): Invalid return from RCGIter::iterate().");
          }
        }
        catch (const StatusTestNaNError& e) {
          // A NaN was detected in the solver.  Set the solution to zero and return unconverged.
          achievedTol_ = MT::one();
          Teuchos::RCP<MV> X = problem_->getLHS();
          MVT::MvInit( *X, SCT::zero() );
          printer_->stream(Warnings) << "Belos::RCGSolMgr::solve(): Warning! NaN has been detected!" 
                                     << std::endl;
          return Unconverged;
        }
        catch (const std::exception &e) {
          printer_->stream(Errors) << "Error! Caught std::exception in RCGIter::iterate() at iteration "
                                   << rcg_iter->getNumIters() << std::endl
                                   << e.what() << std::endl;
          throw;
        }
      }

      // Inform the linear problem that we are finished with this block linear system.
      problem_->setCurrLS();

      // Update indices for the linear systems to be solved.
      numRHS2Solve -= 1;
      if ( numRHS2Solve > 0 ) {
        currIdx[0]++;
        // Set the next indices.
        problem_->setLSIndex( currIdx );
      }
      else {
        currIdx.resize( numRHS2Solve );
      }

      // Update the recycle space for the next linear system
      if (existU1_) { // be sure updated recycle space was created
        // U = U1
        MVT::SetBlock(*U1_,rindex,*U_);
        // Set flag indicating recycle space is now defined
        existU_ = true;
        if (numRHS2Solve > 0) { // also update AU, UTAU, and AUTAU
          // Free pointers in newstate
          newstate.P = Teuchos::null;
          newstate.Ap = Teuchos::null;
          newstate.r = Teuchos::null;
          newstate.z = Teuchos::null;
          newstate.U = Teuchos::null;
          newstate.AU = Teuchos::null;
          newstate.Alpha = Teuchos::null;
          newstate.Beta = Teuchos::null;
          newstate.Beta_i = 0;
          newstate.D = Teuchos::null;
          newstate.Delta = Teuchos::null;
          newstate.LUUTAU = Teuchos::null;
          newstate.ipiv = Teuchos::null;
          newstate.rTz_old = Teuchos::null;

          // Reinitialize AU, UTAU, AUTAU
          Teuchos::RCP<const MV> Utmp  = MVT::CloneView( *U_,  rindex );
          Teuchos::RCP<MV> AUtmp = MVT::CloneViewNonConst( *AU_, rindex );
          // Initialize AU
          problem_->applyOp( *Utmp, *AUtmp );
          // Initialize UTAU
          MVT::MvTransMv( one, *Utmp, *AUtmp, *UTAU_ );
          // Initialize AUTAU  ( AUTAU = AU'*(M\AU) )
          if ( precObj != Teuchos::null ) {
            Teuchos::RCP<MV> LeftPCAU = MVT::CloneViewNonConst( *U1_, rindex ); // use U1 as temp storage
            OPT::Apply( *precObj, *AUtmp, *LeftPCAU );
            MVT::MvTransMv( one, *AUtmp, *LeftPCAU, *AUTAU_ );
          } else {
            MVT::MvTransMv( one, *AUtmp, *AUtmp, *AUTAU_ );
          }
        } // if (numRHS2Solve > 0)

      } // if (existU1)
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

  // Save the convergence test value ("achieved tolerance") for this solve.
  {
    using Teuchos::rcp_dynamic_cast;
    typedef StatusTestGenResNorm<ScalarType,MV,OP,DM> conv_test_type;
    // testValues is nonnull and not persistent.
    const std::vector<MagnitudeType>* pTestValues =
      rcp_dynamic_cast<conv_test_type>(convTest_)->getTestValue();

    TEUCHOS_TEST_FOR_EXCEPTION(pTestValues == NULL, std::logic_error,
      "Belos::RCGSolMgr::solve(): The convergence test's getTestValue() "
      "method returned NULL.  Please report this bug to the Belos developers.");

    TEUCHOS_TEST_FOR_EXCEPTION(pTestValues->size() < 1, std::logic_error,
      "Belos::RCGSolMgr::solve(): The convergence test's getTestValue() "
      "method returned a vector of length zero.  Please report this bug to the "
      "Belos developers.");

    // FIXME (mfh 12 Dec 2011) Does pTestValues really contain the
    // achieved tolerances for all vectors in the current solve(), or
    // just for the vectors from the last deflation?
    achievedTol_ = *std::max_element (pTestValues->begin(), pTestValues->end());
  }

  if (!isConverged) {
    return Unconverged; // return from RCGSolMgr::solve()
  }
  return Converged; // return from RCGSolMgr::solve()
}

//  Compute the harmonic eigenpairs of the projected, dense system.
template<class ScalarType, class MV, class OP, class DM>
void RCGSolMgr<ScalarType,MV,OP,DM,true>::getHarmonicVecs(const DM& F,
                                                          const DM& G,
                                                          DM& Y ) {

  // order of F,G
  int n = DMT::GetNumCols(F);

  // The LAPACK interface
  Teuchos::LAPACK<int,ScalarType> lapack;

  // Magnitude of harmonic Ritz values
  std::vector<MagnitudeType> w(n);

  // Sorted order of harmonic Ritz values
  std::vector<int> iperm(n);

  // Compute k smallest harmonic Ritz pairs
  // SUBROUTINE DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, INFO )
  int itype = 1; // solve A*x = (lambda)*B*x
  char jobz='V'; // compute eigenvalues and eigenvectors
  char uplo='U'; // since F,G symmetric, reference only their upper triangular data
  std::vector<ScalarType> work(1);
  int lwork = -1;
  int info = 0;
  // since SYGV destroys workspace, create copies of F,G
  Teuchos::RCP<DM> F2 = DMT::CreateCopy( F );
  Teuchos::RCP<DM> G2 = DMT::CreateCopy( G );

  DMT::SyncDeviceToHost(*F2);
  DMT::SyncDeviceToHost(*G2);

  // query for optimal workspace size
  lapack.SYGV(itype, jobz, uplo, n, DMT::GetRawHostPtr(*G2), DMT::GetStride(*G2), 
              DMT::GetRawHostPtr(*F2), DMT::GetStride(*F2), &w[0], &work[0], lwork, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, RCGSolMgrLAPACKFailure,
                     "Belos::RCGSolMgr::solve(): LAPACK SYGV failed to query optimal work size.");
  lwork = (int)work[0];
  work.resize(lwork);
  lapack.SYGV(itype, jobz, uplo, n, DMT::GetRawHostPtr(*G2), DMT::GetStride(*G2), 
              DMT::GetRawHostPtr(*F2), DMT::GetStride(*F2), &w[0], &work[0], lwork, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, RCGSolMgrLAPACKFailure,
                     "Belos::RCGSolMgr::solve(): LAPACK SYGV failed to compute eigensolutions.");


  // Construct magnitude of each harmonic Ritz value
  this->sort(w,n,iperm);

  // Select recycledBlocks_ smallest eigenvectors
  for( int i=0; i<recycleBlocks_; i++ ) {
    for( int j=0; j<n; j++ ) {
      DMT::Value(Y,j,i) = DMT::ValueConst(*G2,j,iperm[i]);
    }
  }

}

// This method sorts list of n floating-point numbers and return permutation vector
template<class ScalarType, class MV, class OP, class DM>
void RCGSolMgr<ScalarType,MV,OP,DM,true>::sort(std::vector<ScalarType>& dlist, int n, std::vector<int>& iperm)
{
  int l, r, j, i, flag;
  int    RR2;
  double dRR, dK;

  // Initialize the permutation vector.
  for(j=0;j<n;j++)
    iperm[j] = j;

  if (n <= 1) return;

  l    = n / 2 + 1;
  r    = n - 1;
  l    = l - 1;
  dRR  = dlist[l - 1];
  dK   = dlist[l - 1];

  RR2 = iperm[l - 1];
  while (r != 0) {
    j = l;
    flag = 1;

    while (flag == 1) {
      i = j;
      j = j + j;

      if (j > r + 1)
        flag = 0;
      else {
        if (j < r + 1)
          if (dlist[j] > dlist[j - 1]) j = j + 1;

        if (dlist[j - 1] > dK) {
          dlist[i - 1] = dlist[j - 1];
          iperm[i - 1] = iperm[j - 1];
        }
        else {
          flag = 0;
        }
      }
    }
    dlist[i - 1] = dRR;
    iperm[i - 1] = RR2;
    if (l == 1) {
      dRR  = dlist [r];
      RR2 = iperm[r];
      dK = dlist[r];
      dlist[r] = dlist[0];
      iperm[r] = iperm[0];
      r = r - 1;
    }
    else {
      l   = l - 1;
      dRR  = dlist[l - 1];
      RR2  = iperm[l - 1];
      dK   = dlist[l - 1];
    }
  }
  dlist[0] = dRR;
  iperm[0] = RR2;
}

//  This method requires the solver manager to return a std::string that describes itself.
template<class ScalarType, class MV, class OP, class DM>
std::string RCGSolMgr<ScalarType,MV,OP,DM,true>::description() const
{
  std::ostringstream oss;
  oss << "Belos::RCGSolMgr<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
  return oss.str();
}

} // end Belos namespace

#endif /* BELOS_RCG_SOLMGR_HPP */
