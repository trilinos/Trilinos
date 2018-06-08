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
#include "Teuchos_as.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

/*! \class Belos::RCGSolMgr
\brief Implementation of the RCG (Recycling Conjugate Gradient) iterative linear solver.
\ingroup belos_solver_framework
\author Michael Parks and Heidi Thornquist

\section Belos_GCRODR_summary Summary

This class implements the GCRODR (Recycling GMRES) iterative linear
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

  /** \brief RCGSolMgrRecyclingFailure is thrown when any problem occurs in using/creating
   * the recycling subspace.
   *
   * This exception is thrown from the RCGSolMgr::solve() method.
   *
   */
  class RCGSolMgrRecyclingFailure : public BelosError {public:
    RCGSolMgrRecyclingFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  //@}


  // Partial specialization for unsupported ScalarType types.
  // This contains a stub implementation.
  template<class ScalarType, class MV, class OP,
           const bool supportsScalarType =
             Belos::Details::LapackSupportsScalar<ScalarType>::value &&
             ! Teuchos::ScalarTraits<ScalarType>::isComplex>
  class RCGSolMgr :
    public Details::SolverManagerRequiresRealLapack<ScalarType, MV, OP,
                                                    Belos::Details::LapackSupportsScalar<ScalarType>::value &&
                                                    ! Teuchos::ScalarTraits<ScalarType>::isComplex>
  {
    static const bool scalarTypeIsSupported =
      Belos::Details::LapackSupportsScalar<ScalarType>::value &&
      ! Teuchos::ScalarTraits<ScalarType>::isComplex;
    typedef Details::SolverManagerRequiresRealLapack<ScalarType, MV, OP,
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
    Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const override {
      return Teuchos::rcp(new RCGSolMgr<ScalarType,MV,OP,supportsScalarType>);
    }
  };

  // Partial specialization for real ScalarType.
  // This contains the actual working implementation of RCG.
  // See discussion in the class documentation above.
  template<class ScalarType, class MV, class OP>
  class RCGSolMgr<ScalarType, MV, OP, true> :
    public Details::SolverManagerRequiresRealLapack<ScalarType, MV, OP, true> {
  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
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
    Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const override {
      return Teuchos::rcp(new RCGSolMgr<ScalarType,MV,OP>);
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
    void getHarmonicVecs(const Teuchos::SerialDenseMatrix<int,ScalarType> &F,
                         const Teuchos::SerialDenseMatrix<int,ScalarType> &G,
                         Teuchos::SerialDenseMatrix<int,ScalarType>& Y);

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
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > sTest_;
    Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > maxIterTest_;
    Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> > convTest_;
    Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputTest_;

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
    static constexpr std::ostream * outputStream_default_ = &std::cout;

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
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > Alpha_;
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > Beta_;
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > D_;
    //
    // Solutions to local least-squares problems
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > Delta_;
    //
    // The matrix U^T A U
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > UTAU_;
    //
    // LU factorization of U^T A U
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > LUUTAU_;
    //
    // Data from LU factorization of UTAU
    Teuchos::RCP<std::vector<int> > ipiv_;
    //
    // The matrix (AU)^T AU
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > AUTAU_;
    //
    // The scalar r'*z
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > rTz_old_;
    //
    // Matrices needed for calculation of harmonic Ritz eigenproblem
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > F_,G_,Y_;
    //
    // Matrices needed for updating recycle space
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > L2_,DeltaL2_,AU1TUDeltaL2_;
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > AU1TAU1_, AU1TU1_, AU1TAP_;
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > FY_,GY_;
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > APTAP_;
    Teuchos::RCP<MV> U1Y1_, PY2_;
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > AUTAP_, AU1TU_;
    ScalarType dold;
    /////////////////////////////////////////////////////////////////////////

    // Timers.
    std::string label_;
    Teuchos::RCP<Teuchos::Time> timerSolve_;

    // Internal state variables.
    bool params_Set_;
  };


// Empty Constructor
template<class ScalarType, class MV, class OP>
RCGSolMgr<ScalarType,MV,OP,true>::RCGSolMgr():
  achievedTol_(0.0),
  numIters_(0)
{
  init();
}

// Basic Constructor
template<class ScalarType, class MV, class OP>
RCGSolMgr<ScalarType,MV,OP,true>::RCGSolMgr(
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
template<class ScalarType, class MV, class OP>
void RCGSolMgr<ScalarType,MV,OP,true>::init()
{
  outputStream_ = Teuchos::rcp(outputStream_default_,false);
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

template<class ScalarType, class MV, class OP>
void RCGSolMgr<ScalarType,MV,OP,true>::setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params )
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
  typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
  typedef Belos::StatusTestGenResNorm<ScalarType,MV,OP>  StatusTestResNorm_t;

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


template<class ScalarType, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList>
RCGSolMgr<ScalarType,MV,OP,true>::getValidParameters() const
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
    pl->set("Output Stream", Teuchos::rcp(outputStream_default_,false),
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
template<class ScalarType, class MV, class OP>
void RCGSolMgr<ScalarType,MV,OP,true>::initializeStateStorage() {

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
        Alpha_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( numBlocks_, 1 ) );
      else {
        if ( (Alpha_->numRows() != numBlocks_) || (Alpha_->numCols() != 1) )
          Alpha_->reshape( numBlocks_, 1 );
      }

      // Generate Beta_ only if it doesn't exist, otherwise resize it.
      if (Beta_ == Teuchos::null)
        Beta_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( numBlocks_ + 1, 1 ) );
      else {
        if ( (Beta_->numRows() != (numBlocks_+1)) || (Beta_->numCols() != 1) )
          Beta_->reshape( numBlocks_ + 1, 1 );
      }

      // Generate D_ only if it doesn't exist, otherwise resize it.
      if (D_ == Teuchos::null)
        D_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( numBlocks_ , 1 ) );
      else {
        if ( (D_->numRows() != numBlocks_) || (D_->numCols() != 1) )
          D_->reshape( numBlocks_, 1 );
      }

      // Generate Delta_ only if it doesn't exist, otherwise resize it.
      if (Delta_ == Teuchos::null)
        Delta_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( recycleBlocks_, numBlocks_ + 1 ) );
      else {
        if ( (Delta_->numRows() != recycleBlocks_) || (Delta_->numCols() != (numBlocks_ + 1)) )
          Delta_->reshape( recycleBlocks_, numBlocks_ + 1 );
      }

      // Generate UTAU_ only if it doesn't exist, otherwise resize it.
      if (UTAU_ == Teuchos::null)
        UTAU_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( recycleBlocks_, recycleBlocks_ ) );
      else {
        if ( (UTAU_->numRows() != recycleBlocks_) || (UTAU_->numCols() != recycleBlocks_) )
          UTAU_->reshape( recycleBlocks_, recycleBlocks_ );
      }

      // Generate LUUTAU_ only if it doesn't exist, otherwise resize it.
      if (LUUTAU_ == Teuchos::null)
        LUUTAU_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( recycleBlocks_, recycleBlocks_ ) );
      else {
        if ( (LUUTAU_->numRows() != recycleBlocks_) || (LUUTAU_->numCols() != recycleBlocks_) )
          LUUTAU_->reshape( recycleBlocks_, recycleBlocks_ );
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
        AUTAU_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( recycleBlocks_, recycleBlocks_ ) );
      else {
        if ( (AUTAU_->numRows() != recycleBlocks_) || (AUTAU_->numCols() != recycleBlocks_) )
          AUTAU_->reshape( recycleBlocks_, recycleBlocks_ );
      }

      // Generate rTz_old_ only if it doesn't exist
      if (rTz_old_ == Teuchos::null)
        rTz_old_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( 1, 1 ) );
      else {
        if ( (rTz_old_->numRows() != 1) || (rTz_old_->numCols() != 1) )
          rTz_old_->reshape( 1, 1 );
      }

      // Generate F_ only if it doesn't exist
      if (F_ == Teuchos::null)
        F_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( numBlocks_+recycleBlocks_, numBlocks_+recycleBlocks_ ) );
      else {
        if ( (F_->numRows() != (numBlocks_+recycleBlocks_)) || (F_->numCols() != numBlocks_+recycleBlocks_) )
          F_->reshape( numBlocks_+recycleBlocks_, numBlocks_+recycleBlocks_ );
      }

      // Generate G_ only if it doesn't exist
      if (G_ == Teuchos::null)
        G_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( numBlocks_+recycleBlocks_, numBlocks_+recycleBlocks_ ) );
      else {
        if ( (G_->numRows() != (numBlocks_+recycleBlocks_)) || (G_->numCols() != numBlocks_+recycleBlocks_) )
          G_->reshape( numBlocks_+recycleBlocks_, numBlocks_+recycleBlocks_ );
      }

      // Generate Y_ only if it doesn't exist
      if (Y_ == Teuchos::null)
        Y_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( numBlocks_+recycleBlocks_, recycleBlocks_ ) );
      else {
        if ( (Y_->numRows() != (numBlocks_+recycleBlocks_)) || (Y_->numCols() != numBlocks_+recycleBlocks_) )
          Y_->reshape( numBlocks_+recycleBlocks_, numBlocks_+recycleBlocks_ );
      }

      // Generate L2_ only if it doesn't exist
      if (L2_ == Teuchos::null)
        L2_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( numBlocks_+1, numBlocks_ ) );
      else {
        if ( (L2_->numRows() != (numBlocks_+1)) || (L2_->numCols() != numBlocks_) )
          L2_->reshape( numBlocks_+1, numBlocks_ );
      }

      // Generate DeltaL2_ only if it doesn't exist
      if (DeltaL2_ == Teuchos::null)
        DeltaL2_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( recycleBlocks_, numBlocks_ ) );
      else {
        if ( (DeltaL2_->numRows() != (recycleBlocks_)) || (DeltaL2_->numCols() != (numBlocks_) ) )
          DeltaL2_->reshape( recycleBlocks_, numBlocks_ );
      }

      // Generate AU1TUDeltaL2_ only if it doesn't exist
      if (AU1TUDeltaL2_ == Teuchos::null)
        AU1TUDeltaL2_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( recycleBlocks_, numBlocks_ ) );
      else {
        if ( (AU1TUDeltaL2_->numRows() != (recycleBlocks_)) || (AU1TUDeltaL2_->numCols() != (numBlocks_) ) )
          AU1TUDeltaL2_->reshape( recycleBlocks_, numBlocks_ );
      }

      // Generate AU1TAU1_ only if it doesn't exist
      if (AU1TAU1_ == Teuchos::null)
        AU1TAU1_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( recycleBlocks_, recycleBlocks_ ) );
      else {
        if ( (AU1TAU1_->numRows() != (recycleBlocks_)) || (AU1TAU1_->numCols() != (recycleBlocks_) ) )
          AU1TAU1_->reshape( recycleBlocks_, recycleBlocks_ );
      }

      // Generate GY_ only if it doesn't exist
      if (GY_ == Teuchos::null)
        GY_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( numBlocks_ + recycleBlocks_, recycleBlocks_ ) );
      else {
        if ( (GY_->numRows() != (numBlocks_ + recycleBlocks_)) || (GY_->numCols() != (recycleBlocks_) ) )
          GY_->reshape( numBlocks_+recycleBlocks_, recycleBlocks_ );
      }

      // Generate AU1TU1_ only if it doesn't exist
      if (AU1TU1_ == Teuchos::null)
        AU1TU1_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( recycleBlocks_, recycleBlocks_ ) );
      else {
        if ( (AU1TU1_->numRows() != (recycleBlocks_)) || (AU1TU1_->numCols() != (recycleBlocks_) ) )
          AU1TU1_->reshape( recycleBlocks_, recycleBlocks_ );
      }

      // Generate FY_ only if it doesn't exist
      if (FY_ == Teuchos::null)
        FY_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( numBlocks_ + recycleBlocks_, recycleBlocks_ ) );
      else {
        if ( (FY_->numRows() != (numBlocks_ + recycleBlocks_)) || (FY_->numCols() != (recycleBlocks_) ) )
          FY_->reshape( numBlocks_+recycleBlocks_, recycleBlocks_ );
      }

      // Generate AU1TAP_ only if it doesn't exist
      if (AU1TAP_ == Teuchos::null)
        AU1TAP_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( recycleBlocks_, numBlocks_ ) );
      else {
        if ( (AU1TAP_->numRows() != (recycleBlocks_)) || (AU1TAP_->numCols() != (numBlocks_) ) )
          AU1TAP_->reshape( recycleBlocks_, numBlocks_ );
      }

      // Generate APTAP_ only if it doesn't exist
      if (APTAP_ == Teuchos::null)
        APTAP_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( numBlocks_, numBlocks_ ) );
      else {
        if ( (APTAP_->numRows() != (numBlocks_)) || (APTAP_->numCols() != (numBlocks_) ) )
          APTAP_->reshape( numBlocks_, numBlocks_ );
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
        AUTAP_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( recycleBlocks_, numBlocks_ ) );
      else {
        if ( (AUTAP_->numRows() != (recycleBlocks_)) || (AUTAP_->numCols() != (numBlocks_) ) )
          AUTAP_->reshape( recycleBlocks_, numBlocks_ );
      }

      // Generate AU1TU_ only if it doesn't exist
      if (AU1TU_ == Teuchos::null)
        AU1TU_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( recycleBlocks_, recycleBlocks_ ) );
      else {
        if ( (AU1TU_->numRows() != (recycleBlocks_)) || (AU1TU_->numCols() != (recycleBlocks_) ) )
          AU1TU_->reshape( recycleBlocks_, recycleBlocks_ );
      }


    }
}

template<class ScalarType, class MV, class OP>
ReturnType RCGSolMgr<ScalarType,MV,OP,true>::solve() {

  Teuchos::LAPACK<int,ScalarType> lapack;
  std::vector<int> index(recycleBlocks_);
  ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

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
    index.resize(recycleBlocks_);
    for (int i=0; i<recycleBlocks_; ++i) { index[i] = i; }
    Teuchos::RCP<const MV> Utmp  = MVT::CloneView( *U_,  index );
    index.resize(recycleBlocks_);
    for (int i=0; i<recycleBlocks_; ++i) { index[i] = i; }
    Teuchos::RCP<MV> AUtmp = MVT::CloneViewNonConst( *AU_, index );
    // Initialize AU
    problem_->applyOp( *Utmp, *AUtmp );
    // Initialize UTAU
    Teuchos::SerialDenseMatrix<int,ScalarType> UTAUtmp( Teuchos::View, *UTAU_, recycleBlocks_, recycleBlocks_ );
    MVT::MvTransMv( one, *Utmp, *AUtmp, UTAUtmp );
    // Initialize AUTAU  ( AUTAU = AU'*(M\AU) )
    Teuchos::SerialDenseMatrix<int,ScalarType> AUTAUtmp( Teuchos::View, *AUTAU_, recycleBlocks_, recycleBlocks_ );
    if ( precObj != Teuchos::null ) {
      index.resize(recycleBlocks_);
      for (int i=0; i<recycleBlocks_; ++i) { index[i] = i; }
      index.resize(recycleBlocks_);
      for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
      Teuchos::RCP<MV> PCAU = MVT::CloneViewNonConst( *U1_, index ); // use U1 as temp storage
      OPT::Apply( *precObj, *AUtmp, *PCAU );
      MVT::MvTransMv( one, *AUtmp, *PCAU, AUTAUtmp );
    } else {
      MVT::MvTransMv( one, *AUtmp, *AUtmp, AUTAUtmp );
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // RCG solver

  Teuchos::RCP<RCGIter<ScalarType,MV,OP> > rcg_iter;
  rcg_iter = Teuchos::rcp( new RCGIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,plist) );

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
        Teuchos::SerialDenseMatrix<int,ScalarType> Utr(recycleBlocks_,1);
        index.resize(recycleBlocks_);
        for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
        Teuchos::RCP<const MV> Utmp  = MVT::CloneView( *U_,  index );
        MVT::MvTransMv( one, *Utmp, *r_, Utr );
        Teuchos::SerialDenseMatrix<int,ScalarType> UTAUtmp( Teuchos::View, *UTAU_, recycleBlocks_, recycleBlocks_ );
        Teuchos::SerialDenseMatrix<int,ScalarType> LUUTAUtmp( Teuchos::View, *LUUTAU_, recycleBlocks_, recycleBlocks_ );
        LUUTAUtmp.assign(UTAUtmp);
        int info = 0;
        lapack.GESV(recycleBlocks_, 1, LUUTAUtmp.values(), LUUTAUtmp.stride(), &(*ipiv_)[0], Utr.values(), Utr.stride(), &info);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, RCGSolMgrLAPACKFailure,
                           "Belos::RCGSolMgr::solve(): LAPACK GESV failed to compute a solution.");

        // Update solution (x = x + U*y)
        MVT::MvTimesMatAddMv( one, *Utmp, Utr, one, *problem_->getCurrLHSVec() );

        // Update residual ( r = r - AU*y )
        index.resize(recycleBlocks_);
        for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
        Teuchos::RCP<const MV> AUtmp  = MVT::CloneView( *AU_, index );
        MVT::MvTimesMatAddMv( -one, *AUtmp, Utr, one, *r_ );
      }

      if ( precObj != Teuchos::null ) {
        OPT::Apply( *precObj, *r_, *z_ );
      } else {
        z_ = r_;
      }

      // rTz_old = r'*z
      MVT::MvTransMv( one, *r_, *z_, *rTz_old_ );

      if ( existU_ ) {
        // mu = UTAU\(AU'*z);
        Teuchos::SerialDenseMatrix<int,ScalarType> mu( Teuchos::View, *Delta_, recycleBlocks_, 1);
        index.resize(recycleBlocks_);
        for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
        Teuchos::RCP<const MV> AUtmp  = MVT::CloneView( *AU_, index );
        MVT::MvTransMv( one, *AUtmp, *z_, mu );
        Teuchos::SerialDenseMatrix<int,ScalarType> LUUTAUtmp( Teuchos::View, *LUUTAU_, recycleBlocks_, recycleBlocks_ );
        char TRANS = 'N';
        int info;
        lapack.GETRS( TRANS, recycleBlocks_, 1, LUUTAUtmp.values(), LUUTAUtmp.stride(), &(*ipiv_)[0], mu.values(), mu.stride(), &info );
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, RCGSolMgrLAPACKFailure,
                           "Belos::RCGSolMgr::solve(): LAPACK GETRS failed to compute a solution.");
        // p  = z - U*mu;
        index.resize( 1 );
        index[0] = 0;
        Teuchos::RCP<MV> Ptmp  = MVT::CloneViewNonConst( *P_, index );
        MVT::MvAddMv(one,*z_,zero,*z_,*Ptmp);
        MVT::MvTimesMatAddMv( -one, *U_, mu, one, *Ptmp );
      } else {
        // p = z;
        index.resize( 1 );
        index[0] = 0;
        Teuchos::RCP<MV> Ptmp  = MVT::CloneViewNonConst( *P_, index );
        MVT::MvAddMv(one,*z_,zero,*z_,*Ptmp);
      }

      // Set the new state and initialize the solver.
      RCGIterState<ScalarType,MV> newstate;

      // Create RCP views here
      index.resize( numBlocks_+1 );
      for (int ii=0; ii<(numBlocks_+1); ++ii) { index[ii] = ii; }
      newstate.P  = MVT::CloneViewNonConst( *P_,  index );
      index.resize( recycleBlocks_ );
      for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
      newstate.U  = MVT::CloneViewNonConst( *U_,  index );
      index.resize( recycleBlocks_ );
      for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
      newstate.AU  = MVT::CloneViewNonConst( *AU_,  index );
      newstate.Alpha = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::View, *Alpha_, numBlocks_, 1 ) );
      newstate.Beta = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::View, *Beta_, numBlocks_, 1 ) );
      newstate.D = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::View, *D_, numBlocks_, 1 ) );
      newstate.Delta = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::View, *Delta_, recycleBlocks_, numBlocks_, 0, 1 ) );
      newstate.LUUTAU = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::View, *LUUTAU_, recycleBlocks_, recycleBlocks_ ) );
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

                   Teuchos::SerialDenseMatrix<int,ScalarType> Ftmp( Teuchos::View, *F_, numBlocks_, numBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Gtmp( Teuchos::View, *G_, numBlocks_, numBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Dtmp( Teuchos::View, *D_, numBlocks_, 1 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Alphatmp( Teuchos::View, *Alpha_, numBlocks_, 1 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Betatmp( Teuchos::View, *Beta_, numBlocks_, 1 );
                   Ftmp.putScalar(zero);
                   Gtmp.putScalar(zero);
                   for (int ii=0;ii<numBlocks_;ii++) {
                     Gtmp(ii,ii) = (Dtmp(ii,0) / Alphatmp(ii,0))*(1 + Betatmp(ii,0));
                     if (ii > 0) {
                       Gtmp(ii-1,ii) = -Dtmp(ii,0)/Alphatmp(ii-1,0);
                       Gtmp(ii,ii-1) = -Dtmp(ii,0)/Alphatmp(ii-1,0);
                     }
                     Ftmp(ii,ii) = Dtmp(ii,0);
                   }

                   // compute harmonic Ritz vectors
                   Teuchos::SerialDenseMatrix<int,ScalarType> Ytmp( Teuchos::View, *Y_, numBlocks_, recycleBlocks_ );
                   getHarmonicVecs(Ftmp,Gtmp,Ytmp);

                   // U1 = [P(:,1:end-1)*Y];
                   index.resize( numBlocks_ );
                   for (int ii=0; ii<numBlocks_; ++ii) { index[ii] = ii; }
                   Teuchos::RCP<const MV> Ptmp = MVT::CloneView( *P_,  index );
                   index.resize( recycleBlocks_ );
                   for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
                   Teuchos::RCP<MV> U1tmp  = MVT::CloneViewNonConst( *U1_,  index );
                   MVT::MvTimesMatAddMv( one, *Ptmp, Ytmp, zero, *U1tmp );

                   // Precompute some variables for next cycle

                   // AU1TAU1     = Y'*G*Y;
                   Teuchos::SerialDenseMatrix<int,ScalarType> GYtmp( Teuchos::View, *GY_, numBlocks_, recycleBlocks_ );
                   GYtmp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,Gtmp,Ytmp,zero);
                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TAU1tmp( Teuchos::View, *AU1TAU1_, recycleBlocks_, recycleBlocks_ );
                   AU1TAU1tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,Ytmp,GYtmp,zero);

                   // AU1TU1      = Y'*F*Y;
                   Teuchos::SerialDenseMatrix<int,ScalarType> FYtmp( Teuchos::View, *FY_, numBlocks_, recycleBlocks_ );
                   FYtmp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,Ftmp,Ytmp,zero);
                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TU1tmp( Teuchos::View, *AU1TU1_, recycleBlocks_, recycleBlocks_ );
                   AU1TU1tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,Ytmp,FYtmp,zero);

                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TAPtmp( Teuchos::View, *AU1TAP_, recycleBlocks_, 1 );
                   // Must reinitialize AU1TAP; can become dense later
                   AU1TAPtmp.putScalar(zero);
                   // AU1TAP(:,1) = Y(end,:)' * (-1/Alpha(end));
                   ScalarType alphatmp = -1.0 / Alphatmp(numBlocks_-1,0);
                   for (int ii=0; ii<recycleBlocks_; ++ii) {
                      AU1TAPtmp(ii,0) = Ytmp(numBlocks_-1,ii) * alphatmp;
                   }

                   // indicate that updated recycle space now defined
                   existU1_ = true;

                   // Indicate the size of the P, Beta structures generated this cycle
                   lastp = numBlocks_;
                   lastBeta = numBlocks_-1;

                } // if (cycle == 0)
                else { // No U, but U1 guaranteed to exist now

                   // Finish computation of subblocks
                   // AU1TAP = AU1TAP * D(1);
                   Teuchos::SerialDenseMatrix<int,ScalarType> Dtmp( Teuchos::View, *D_, numBlocks_, 1 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TAPtmp( Teuchos::View, *AU1TAP_, recycleBlocks_, numBlocks_ );
                   AU1TAPtmp.scale(Dtmp(0,0));

                   Teuchos::SerialDenseMatrix<int,ScalarType> Alphatmp( Teuchos::View, *Alpha_, numBlocks_, 1 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Betatmp( Teuchos::View, *Beta_, numBlocks_+1, 1 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> APTAPtmp( Teuchos::View, *APTAP_, numBlocks_, numBlocks_ );
                   APTAPtmp.putScalar(zero);
                   for (int ii=0; ii<numBlocks_; ii++) {
                     APTAPtmp(ii,ii) = (Dtmp(ii,0) / Alphatmp(ii,0))*(1 + Betatmp(ii+1,0));
                     if (ii > 0) {
                       APTAPtmp(ii-1,ii) = -Dtmp(ii,0)/Alphatmp(ii-1,0);
                       APTAPtmp(ii,ii-1) = -Dtmp(ii,0)/Alphatmp(ii-1,0);
                     }
                   }

                   // F = [AU1TU1 zeros(k,m); zeros(m,k) diag(D)];
                   Teuchos::SerialDenseMatrix<int,ScalarType> Ftmp( Teuchos::View, *F_, (numBlocks_+recycleBlocks_), (numBlocks_+recycleBlocks_) );
                   Teuchos::SerialDenseMatrix<int,ScalarType> F11( Teuchos::View, *F_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> F12( Teuchos::View, *F_, recycleBlocks_, numBlocks_, 0, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> F21( Teuchos::View, *F_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> F22( Teuchos::View, *F_, numBlocks_, numBlocks_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TU1tmp( Teuchos::View, *AU1TU1_, recycleBlocks_, recycleBlocks_ );
                   F11.assign(AU1TU1tmp);
                   F12.putScalar(zero);
                   F21.putScalar(zero);
                   F22.putScalar(zero);
                   for(int ii=0;ii<numBlocks_;ii++) {
                     F22(ii,ii) = Dtmp(ii,0);
                   }

                   // G = [AU1TAU1 AU1TAP; AU1TAP' APTAP];
                   Teuchos::SerialDenseMatrix<int,ScalarType> Gtmp( Teuchos::View, *G_, (numBlocks_+recycleBlocks_), (numBlocks_+recycleBlocks_) );
                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TAU1tmp( Teuchos::View, *AU1TAU1_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> G11( Teuchos::View, *G_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> G12( Teuchos::View, *G_, recycleBlocks_, numBlocks_, 0, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> G21( Teuchos::View, *G_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> G22( Teuchos::View, *G_, numBlocks_, numBlocks_, recycleBlocks_, recycleBlocks_ );
                   G11.assign(AU1TAU1tmp);
                   G12.assign(AU1TAPtmp);
                   // G21 = G12'; (no transpose operator exists for SerialDenseMatrix; Do copy manually)
                   for (int ii=0;ii<recycleBlocks_;++ii)
                     for (int jj=0;jj<numBlocks_;++jj)
                       G21(jj,ii) = G12(ii,jj);
                   G22.assign(APTAPtmp);

                   // compute harmonic Ritz vectors
                   Teuchos::SerialDenseMatrix<int,ScalarType> Ytmp( Teuchos::View, *Y_, (recycleBlocks_+numBlocks_), recycleBlocks_ );
                   getHarmonicVecs(Ftmp,Gtmp,Ytmp);

                   // U1 = [U1 P(:,2:end-1)]*Y;
                   index.resize( numBlocks_ );
                   for (int ii=0; ii<numBlocks_; ++ii) { index[ii] = ii+1; }
                   Teuchos::RCP<const MV> Ptmp = MVT::CloneView( *P_,  index );
                   index.resize( recycleBlocks_ );
                   for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
                   Teuchos::RCP<MV> PY2tmp = MVT::CloneViewNonConst( *PY2_,  index );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Y2( Teuchos::View, *Y_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   index.resize( recycleBlocks_ );
                   for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
                   Teuchos::RCP<MV> U1tmp  = MVT::CloneViewNonConst( *U1_,  index );
                   index.resize( recycleBlocks_ );
                   for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
                   Teuchos::RCP<MV> U1Y1tmp  = MVT::CloneViewNonConst( *U1Y1_,  index );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Y1( Teuchos::View, *Y_, recycleBlocks_, recycleBlocks_ );
                   MVT::MvTimesMatAddMv( one, *Ptmp, Y2, zero, *PY2tmp );
                   MVT::MvTimesMatAddMv( one, *U1tmp, Y1, zero, *U1Y1tmp );
                   MVT::MvAddMv(one,*U1Y1tmp, one, *PY2tmp, *U1tmp);

                   // Precompute some variables for next cycle

                   // AU1TAU1     = Y'*G*Y;
                   Teuchos::SerialDenseMatrix<int,ScalarType> GYtmp( Teuchos::View, *GY_, (numBlocks_+recycleBlocks_), recycleBlocks_ );
                   GYtmp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,Gtmp,Ytmp,zero);
                   AU1TAU1tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,Ytmp,GYtmp,zero);

                   // AU1TAP      = zeros(k,m);
                   // AU1TAP(:,1) = Y(end,:)' * (-1/Alpha(end));
                   AU1TAPtmp.putScalar(zero);
                   ScalarType alphatmp = -1.0 / Alphatmp(numBlocks_-1,0);
                   for (int ii=0; ii<recycleBlocks_; ++ii) {
                      AU1TAPtmp(ii,0) = Ytmp(numBlocks_+recycleBlocks_-1,ii) * alphatmp;
                   }

                   // AU1TU1      = Y'*F*Y;
                   Teuchos::SerialDenseMatrix<int,ScalarType> FYtmp( Teuchos::View, *FY_, (numBlocks_+recycleBlocks_), recycleBlocks_ );
                   FYtmp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,Ftmp,Ytmp,zero);
                   //Teuchos::SerialDenseMatrix<int,ScalarType> AU1TU1tmp( Teuchos::View, *AU1TU1_, recycleBlocks_, recycleBlocks_ );
                   AU1TU1tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,Ytmp,FYtmp,zero);

                   // Indicate the size of the P, Beta structures generated this cycle
                   lastp = numBlocks_+1;
                   lastBeta = numBlocks_;

                } // if (cycle != 1)
              } // if (!existU_)
              else { // U exists
                if (cycle == 0) { // No U1, but U exists
                   Teuchos::SerialDenseMatrix<int,ScalarType> Alphatmp( Teuchos::View, *Alpha_, numBlocks_, 1 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Betatmp( Teuchos::View, *Beta_, numBlocks_, 1 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Dtmp( Teuchos::View, *D_, numBlocks_, 1 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> APTAPtmp( Teuchos::View, *APTAP_, numBlocks_, numBlocks_ );
                   APTAPtmp.putScalar(zero);
                   for (int ii=0; ii<numBlocks_; ii++) {
                     APTAPtmp(ii,ii) = (Dtmp(ii,0) / Alphatmp(ii,0))*(1 + Betatmp(ii,0));
                     if (ii > 0) {
                       APTAPtmp(ii-1,ii) = -Dtmp(ii,0)/Alphatmp(ii-1,0);
                       APTAPtmp(ii,ii-1) = -Dtmp(ii,0)/Alphatmp(ii-1,0);
                     }
                   }

                   Teuchos::SerialDenseMatrix<int,ScalarType> L2tmp( Teuchos::View, *L2_, numBlocks_+1, numBlocks_ );
                   L2tmp.putScalar(zero);
                   for(int ii=0;ii<numBlocks_;ii++) {
                     L2tmp(ii,ii)   = 1./Alphatmp(ii,0);
                     L2tmp(ii+1,ii) = -1./Alphatmp(ii,0);
                   }

                   // AUTAP = UTAU*Delta*L2;
                   Teuchos::SerialDenseMatrix<int,ScalarType> AUTAPtmp( Teuchos::View, *AUTAP_, recycleBlocks_, numBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> UTAUtmp( Teuchos::View, *UTAU_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Deltatmp( Teuchos::View, *Delta_, recycleBlocks_, numBlocks_+1 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> DeltaL2tmp( Teuchos::View, *DeltaL2_, recycleBlocks_, numBlocks_ );
                   DeltaL2tmp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,Deltatmp,L2tmp,zero);
                   AUTAPtmp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,UTAUtmp,DeltaL2tmp,zero);

                   // F = [UTAU zeros(k,m); zeros(m,k) diag(D)];
                   Teuchos::SerialDenseMatrix<int,ScalarType> Ftmp( Teuchos::View, *F_, (numBlocks_+recycleBlocks_), (numBlocks_+recycleBlocks_) );
                   Teuchos::SerialDenseMatrix<int,ScalarType> F11( Teuchos::View, *F_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> F12( Teuchos::View, *F_, recycleBlocks_, numBlocks_, 0, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> F21( Teuchos::View, *F_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> F22( Teuchos::View, *F_, numBlocks_, numBlocks_, recycleBlocks_, recycleBlocks_ );
                   F11.assign(UTAUtmp);
                   F12.putScalar(zero);
                   F21.putScalar(zero);
                   F22.putScalar(zero);
                   for(int ii=0;ii<numBlocks_;ii++) {
                     F22(ii,ii) = Dtmp(ii,0);
                   }

                   // G = [AUTAU AUTAP; AUTAP' APTAP];
                   Teuchos::SerialDenseMatrix<int,ScalarType> Gtmp( Teuchos::View, *G_, (numBlocks_+recycleBlocks_), (numBlocks_+recycleBlocks_) );
                   Teuchos::SerialDenseMatrix<int,ScalarType> G11( Teuchos::View, *G_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> G12( Teuchos::View, *G_, recycleBlocks_, numBlocks_, 0, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> G21( Teuchos::View, *G_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> G22( Teuchos::View, *G_, numBlocks_, numBlocks_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> AUTAUtmp( Teuchos::View, *AUTAU_, recycleBlocks_, recycleBlocks_ );
                   G11.assign(AUTAUtmp);
                   G12.assign(AUTAPtmp);
                   // G21 = G12'; (no transpose operator exists for SerialDenseMatrix; Do copy manually)
                   for (int ii=0;ii<recycleBlocks_;++ii)
                     for (int jj=0;jj<numBlocks_;++jj)
                       G21(jj,ii) = G12(ii,jj);
                   G22.assign(APTAPtmp);

                   // compute harmonic Ritz vectors
                   Teuchos::SerialDenseMatrix<int,ScalarType> Ytmp( Teuchos::View, *Y_, (recycleBlocks_+numBlocks_), recycleBlocks_ );
                   getHarmonicVecs(Ftmp,Gtmp,Ytmp);

                   // U1 = [U P(:,1:end-1)]*Y;
                   index.resize( recycleBlocks_ );
                   for (int ii=0; ii<(recycleBlocks_); ++ii) { index[ii] = ii; }
                   Teuchos::RCP<const MV> Utmp = MVT::CloneView( *U_,  index );
                   index.resize( numBlocks_ );
                   for (int ii=0; ii<numBlocks_; ++ii) { index[ii] = ii; }
                   Teuchos::RCP<const MV> Ptmp = MVT::CloneView( *P_,  index );
                   index.resize( recycleBlocks_ );
                   for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
                   Teuchos::RCP<MV> PY2tmp = MVT::CloneViewNonConst( *PY2_,  index );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Y2( Teuchos::View, *Y_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   index.resize( recycleBlocks_ );
                   for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
                   Teuchos::RCP<MV> UY1tmp  = MVT::CloneViewNonConst( *U1Y1_,  index );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Y1( Teuchos::View, *Y_, recycleBlocks_, recycleBlocks_ );
                   index.resize( recycleBlocks_ );
                   for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
                   Teuchos::RCP<MV> U1tmp  = MVT::CloneViewNonConst( *U1_,  index );
                   MVT::MvTimesMatAddMv( one, *Ptmp, Y2, zero, *PY2tmp );
                   MVT::MvTimesMatAddMv( one, *Utmp, Y1, zero, *UY1tmp );
                   MVT::MvAddMv(one,*UY1tmp, one, *PY2tmp, *U1tmp);

                   // Precompute some variables for next cycle

                   // AU1TAU1     = Y'*G*Y;
                   Teuchos::SerialDenseMatrix<int,ScalarType> GYtmp( Teuchos::View, *GY_, (numBlocks_+recycleBlocks_), recycleBlocks_ );
                   GYtmp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,Gtmp,Ytmp,zero);
                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TAU1tmp( Teuchos::View, *AU1TAU1_, recycleBlocks_, recycleBlocks_ );
                   AU1TAU1tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,Ytmp,GYtmp,zero);

                   // AU1TU1      = Y'*F*Y;
                   Teuchos::SerialDenseMatrix<int,ScalarType> FYtmp( Teuchos::View, *FY_, (numBlocks_+recycleBlocks_), recycleBlocks_ );
                   FYtmp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,Ftmp,Ytmp,zero);
                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TU1tmp( Teuchos::View, *AU1TU1_, recycleBlocks_, recycleBlocks_ );
                   AU1TU1tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,Ytmp,FYtmp,zero);

                   // AU1TU   = UTAU;
                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TUtmp( Teuchos::View, *AU1TU_, recycleBlocks_, recycleBlocks_ );
                   AU1TUtmp.assign(UTAUtmp);

                   // dold    = D(end);
                   dold = Dtmp(numBlocks_-1,0);

                   // indicate that updated recycle space now defined
                   existU1_ = true;

                   // Indicate the size of the P, Beta structures generated this cycle
                   lastp = numBlocks_;
                   lastBeta = numBlocks_-1;
                }
                else { // Have U and U1
                   Teuchos::SerialDenseMatrix<int,ScalarType> Alphatmp( Teuchos::View, *Alpha_, numBlocks_, 1 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Betatmp( Teuchos::View, *Beta_, numBlocks_+1, 1 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Dtmp( Teuchos::View, *D_, numBlocks_, 1 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> APTAPtmp( Teuchos::View, *APTAP_, numBlocks_, numBlocks_ );
                   for (int ii=0; ii<numBlocks_; ii++) {
                     APTAPtmp(ii,ii) = (Dtmp(ii,0) / Alphatmp(ii,0))*(1 + Betatmp(ii+1,0));
                     if (ii > 0) {
                       APTAPtmp(ii-1,ii) = -Dtmp(ii,0)/Alphatmp(ii-1,0);
                       APTAPtmp(ii,ii-1) = -Dtmp(ii,0)/Alphatmp(ii-1,0);
                     }
                   }

                   Teuchos::SerialDenseMatrix<int,ScalarType> L2tmp( Teuchos::View, *L2_, numBlocks_+1, numBlocks_ );
                   for(int ii=0;ii<numBlocks_;ii++) {
                     L2tmp(ii,ii)   = 1./Alphatmp(ii,0);
                     L2tmp(ii+1,ii) = -1./Alphatmp(ii,0);
                   }

                   // M(end,1) = dold*(-Beta(1)/Alpha(1));
                   // AU1TAP = Y'*[AU1TU*Delta*L2; M];
                   Teuchos::SerialDenseMatrix<int,ScalarType> DeltaL2( Teuchos::View, *DeltaL2_, recycleBlocks_, numBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> Deltatmp( Teuchos::View, *Delta_, recycleBlocks_, numBlocks_+1 );
                   DeltaL2.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,Deltatmp,L2tmp,zero);
                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TUDeltaL2( Teuchos::View, *AU1TUDeltaL2_, recycleBlocks_, numBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TUtmp( Teuchos::View, *AU1TU_, recycleBlocks_, recycleBlocks_ );
                   AU1TUDeltaL2.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,AU1TUtmp,DeltaL2,zero);
                   Teuchos::SerialDenseMatrix<int,ScalarType> Y1( Teuchos::View, *Y_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TAPtmp( Teuchos::View, *AU1TAP_, recycleBlocks_, numBlocks_ );
                   AU1TAPtmp.putScalar(zero);
                   AU1TAPtmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,Y1,AU1TUDeltaL2,zero);
                   Teuchos::SerialDenseMatrix<int,ScalarType> Y2( Teuchos::View, *Y_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   ScalarType val = dold * (-Betatmp(0,0)/Alphatmp(0,0));
                   for(int ii=0;ii<recycleBlocks_;ii++) {
                     AU1TAPtmp(ii,0) += Y2(numBlocks_-1,ii)*val;
                   }

                   // AU1TU = Y1'*AU1TU
                   Teuchos::SerialDenseMatrix<int,ScalarType> Y1TAU1TU( Teuchos::View, *GY_, recycleBlocks_, recycleBlocks_ );
                   Y1TAU1TU.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,Y1,AU1TUtmp,zero);
                   AU1TUtmp.assign(Y1TAU1TU);

                   // F = [AU1TU1 zeros(k,m); zeros(m,k) diag(D)];
                   Teuchos::SerialDenseMatrix<int,ScalarType> Ftmp( Teuchos::View, *F_, (numBlocks_+recycleBlocks_), (numBlocks_+recycleBlocks_) );
                   Teuchos::SerialDenseMatrix<int,ScalarType> F11( Teuchos::View, *F_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> F12( Teuchos::View, *F_, recycleBlocks_, numBlocks_, 0, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> F21( Teuchos::View, *F_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> F22( Teuchos::View, *F_, numBlocks_, numBlocks_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TU1tmp( Teuchos::View, *AU1TU1_, recycleBlocks_, recycleBlocks_ );
                   F11.assign(AU1TU1tmp);
                   for(int ii=0;ii<numBlocks_;ii++) {
                     F22(ii,ii) = Dtmp(ii,0);
                   }

                   // G = [AU1TAU1 AU1TAP; AU1TAP' APTAP];
                   Teuchos::SerialDenseMatrix<int,ScalarType> Gtmp( Teuchos::View, *G_, (numBlocks_+recycleBlocks_), (numBlocks_+recycleBlocks_) );
                   Teuchos::SerialDenseMatrix<int,ScalarType> G11( Teuchos::View, *G_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> G12( Teuchos::View, *G_, recycleBlocks_, numBlocks_, 0, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> G21( Teuchos::View, *G_, numBlocks_, recycleBlocks_, recycleBlocks_, 0 );
                   Teuchos::SerialDenseMatrix<int,ScalarType> G22( Teuchos::View, *G_, numBlocks_, numBlocks_, recycleBlocks_, recycleBlocks_ );
                   Teuchos::SerialDenseMatrix<int,ScalarType> AU1TAU1tmp( Teuchos::View, *AU1TAU1_, recycleBlocks_, recycleBlocks_ );
                   G11.assign(AU1TAU1tmp);
                   G12.assign(AU1TAPtmp);
                   // G21 = G12'; (no transpose operator exists for SerialDenseMatrix; Do copy manually)
                   for (int ii=0;ii<recycleBlocks_;++ii)
                     for (int jj=0;jj<numBlocks_;++jj)
                       G21(jj,ii) = G12(ii,jj);
                   G22.assign(APTAPtmp);

                   // compute harmonic Ritz vectors
                   Teuchos::SerialDenseMatrix<int,ScalarType> Ytmp( Teuchos::View, *Y_, (recycleBlocks_+numBlocks_), recycleBlocks_ );
                   getHarmonicVecs(Ftmp,Gtmp,Ytmp);

                   // U1 = [U1 P(:,2:end-1)]*Y;
                   index.resize( numBlocks_ );
                   for (int ii=0; ii<numBlocks_; ++ii) { index[ii] = ii+1; }
                   Teuchos::RCP<const MV> Ptmp = MVT::CloneView( *P_,  index );
                   index.resize( recycleBlocks_ );
                   for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
                   Teuchos::RCP<MV> PY2tmp = MVT::CloneViewNonConst( *PY2_,  index );
                   index.resize( recycleBlocks_ );
                   for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
                   Teuchos::RCP<MV> U1tmp  = MVT::CloneViewNonConst( *U1_,  index );
                   index.resize( recycleBlocks_ );
                   for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
                   Teuchos::RCP<MV> U1Y1tmp  = MVT::CloneViewNonConst( *U1Y1_,  index );
                   MVT::MvTimesMatAddMv( one, *Ptmp, Y2, zero, *PY2tmp );
                   MVT::MvTimesMatAddMv( one, *U1tmp, Y1, zero, *U1Y1tmp );
                   MVT::MvAddMv(one,*U1Y1tmp, one, *PY2tmp, *U1tmp);

                   // Precompute some variables for next cycle

                   // AU1TAU1     = Y'*G*Y;
                   Teuchos::SerialDenseMatrix<int,ScalarType> GYtmp( Teuchos::View, *GY_, (numBlocks_+recycleBlocks_), recycleBlocks_ );
                   GYtmp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,Gtmp,Ytmp,zero);
                   AU1TAU1tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,Ytmp,GYtmp,zero);

                   // AU1TU1      = Y'*F*Y;
                   Teuchos::SerialDenseMatrix<int,ScalarType> FYtmp( Teuchos::View, *FY_, (numBlocks_+recycleBlocks_), recycleBlocks_ );
                   FYtmp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,Ftmp,Ytmp,zero);
                   AU1TU1tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,Ytmp,FYtmp,zero);

                   // dold    = D(end);
                   dold = Dtmp(numBlocks_-1,0);

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
            (*Beta_)(0,0) = (*Beta_)(lastBeta,0);

            // Delta = Delta(:,end);
            if (existU_) { // Delta only initialized if U exists
              Teuchos::SerialDenseMatrix<int,ScalarType> mu1( Teuchos::View, *Delta_, recycleBlocks_, 1, 0, 0 );
              Teuchos::SerialDenseMatrix<int,ScalarType> mu2( Teuchos::View, *Delta_, recycleBlocks_, 1, 0, numBlocks_ );
              mu1.assign(mu2);
            }

            // Now reinitialize state variables for next cycle
            newstate.P = Teuchos::null;
            index.resize( numBlocks_+1 );
            for (int ii=0; ii<(numBlocks_+1); ++ii) { index[ii] = ii+1; }
            newstate.P  = MVT::CloneViewNonConst( *P_,  index );

            newstate.Beta = Teuchos::null;
            newstate.Beta = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::View, *Beta_, numBlocks_, 1, 1, 0 ) );

            newstate.Delta = Teuchos::null;
            newstate.Delta = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::View, *Delta_, recycleBlocks_, numBlocks_, 0, 1 ) );

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
        index.resize(recycleBlocks_);
        for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
        MVT::SetBlock(*U1_,index,*U_);
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
          newstate.D = Teuchos::null;
          newstate.Delta = Teuchos::null;
          newstate.LUUTAU = Teuchos::null;
          newstate.ipiv = Teuchos::null;
          newstate.rTz_old = Teuchos::null;

          // Reinitialize AU, UTAU, AUTAU
          index.resize(recycleBlocks_);
          for (int i=0; i<recycleBlocks_; ++i) { index[i] = i; }
          Teuchos::RCP<const MV> Utmp  = MVT::CloneView( *U_,  index );
          index.resize(recycleBlocks_);
          for (int i=0; i<recycleBlocks_; ++i) { index[i] = i; }
          Teuchos::RCP<MV> AUtmp = MVT::CloneViewNonConst( *AU_, index );
          // Initialize AU
          problem_->applyOp( *Utmp, *AUtmp );
          // Initialize UTAU
          Teuchos::SerialDenseMatrix<int,ScalarType> UTAUtmp( Teuchos::View, *UTAU_, recycleBlocks_, recycleBlocks_ );
          MVT::MvTransMv( one, *Utmp, *AUtmp, UTAUtmp );
          // Initialize AUTAU  ( AUTAU = AU'*(M\AU) )
          Teuchos::SerialDenseMatrix<int,ScalarType> AUTAUtmp( Teuchos::View, *AUTAU_, recycleBlocks_, recycleBlocks_ );
          if ( precObj != Teuchos::null ) {
            index.resize(recycleBlocks_);
            for (int i=0; i<recycleBlocks_; ++i) { index[i] = i; }
            index.resize(recycleBlocks_);
            for (int ii=0; ii<recycleBlocks_; ++ii) { index[ii] = ii; }
            Teuchos::RCP<MV> LeftPCAU = MVT::CloneViewNonConst( *U1_, index ); // use U1 as temp storage
            OPT::Apply( *precObj, *AUtmp, *LeftPCAU );
            MVT::MvTransMv( one, *AUtmp, *LeftPCAU, AUTAUtmp );
          } else {
            MVT::MvTransMv( one, *AUtmp, *AUtmp, AUTAUtmp );
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
    typedef StatusTestGenResNorm<ScalarType,MV,OP> conv_test_type;
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
template<class ScalarType, class MV, class OP>
void RCGSolMgr<ScalarType,MV,OP,true>::getHarmonicVecs(const Teuchos::SerialDenseMatrix<int,ScalarType>& F,
                                                  const Teuchos::SerialDenseMatrix<int,ScalarType>& G,
                                                        Teuchos::SerialDenseMatrix<int,ScalarType>& Y) {
  // order of F,G
  int n = F.numCols();

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
  Teuchos::SerialDenseMatrix<int,ScalarType> F2( Teuchos::Copy, *F_ );
  Teuchos::SerialDenseMatrix<int,ScalarType> G2( Teuchos::Copy, *G_ );

  // query for optimal workspace size
  lapack.SYGV(itype, jobz, uplo, n, G2.values(), G2.stride(), F2.values(), F2.stride(), &w[0], &work[0], lwork, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, RCGSolMgrLAPACKFailure,
                     "Belos::RCGSolMgr::solve(): LAPACK SYGV failed to query optimal work size.");
  lwork = (int)work[0];
  work.resize(lwork);
  lapack.SYGV(itype, jobz, uplo, n, G2.values(), G2.stride(), F2.values(), F2.stride(), &w[0], &work[0], lwork, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, RCGSolMgrLAPACKFailure,
                     "Belos::RCGSolMgr::solve(): LAPACK SYGV failed to compute eigensolutions.");


  // Construct magnitude of each harmonic Ritz value
  this->sort(w,n,iperm);

  // Select recycledBlocks_ smallest eigenvectors
  for( int i=0; i<recycleBlocks_; i++ ) {
    for( int j=0; j<n; j++ ) {
      Y(j,i) = G2(j,iperm[i]);
    }
  }

}

// This method sorts list of n floating-point numbers and return permutation vector
template<class ScalarType, class MV, class OP>
void RCGSolMgr<ScalarType,MV,OP,true>::sort(std::vector<ScalarType>& dlist, int n, std::vector<int>& iperm)
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
template<class ScalarType, class MV, class OP>
std::string RCGSolMgr<ScalarType,MV,OP,true>::description() const
{
  std::ostringstream oss;
  oss << "Belos::RCGSolMgr<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
  return oss.str();
}

} // end Belos namespace

#endif /* BELOS_RCG_SOLMGR_HPP */
