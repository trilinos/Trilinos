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

#ifndef BELOS_PCPG_SOLMGR_HPP
#define BELOS_PCPG_SOLMGR_HPP

/// \file BelosPCPGSolMgr.hpp
/// \brief Declaration and definition of Belos::PCPGSolMgr
///   (PCPG iterative linear solver).

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"

#include "BelosPCPGIter.hpp"

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
#  include "Teuchos_TimeMonitor.hpp"
#endif
#if defined(HAVE_TEUCHOSCORE_CXX11)
#  include <type_traits>
#endif // defined(HAVE_TEUCHOSCORE_CXX11)
#include "Teuchos_TypeTraits.hpp"

namespace Belos {

  //! @name PCPGSolMgr Exceptions
  //@{

  /** \brief PCPGSolMgrLinearProblemFailure is thrown when the linear problem is
   * not setup (i.e. setProblem() was not called) when solve() is called.
   *
   * This exception is thrown from the PCPGSolMgr::solve() method.
   *
   */
  class PCPGSolMgrLinearProblemFailure : public BelosError {public:
    PCPGSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief PCPGSolMgrOrthoFailure is thrown when the orthogonalization manager is
   * unable to generate orthonormal columns from the initial basis vectors.
   * This exception is thrown from the PCPGSolMgr::solve() method.
   *
   */
  class PCPGSolMgrOrthoFailure : public BelosError {public:
    PCPGSolMgrOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief PCPGSolMgrLAPACKFailure is thrown when a nonzero value is retuned
   * from an LAPACK call.
   *
   * This exception is thrown from the PCPGSolMgr::solve() method.
   *
   */
  class PCPGSolMgrLAPACKFailure : public BelosError {public:
    PCPGSolMgrLAPACKFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief PCPGSolMgrRecyclingFailure is thrown when any problem occurs in using/creating
   * the recycling subspace.
   *
   * The PCPGSolMgr::solve() method throws the exception.
   *
   */
  class PCPGSolMgrRecyclingFailure : public BelosError {public:
    PCPGSolMgrRecyclingFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  //@}


  /// \class Belos::PCPGSolMgr
  /// \brief PCPG iterative linear solver.
  /// \author David Day
  /// \ingroup belos_solver_framework
  ///
  /// PCPG is a CG-based "seed solver."  This means that it does
  /// preconditioned CG to build up a matrix polynomial, then can
  /// reuse that polynomial to compute solutions of successive linear
  /// systems, possibly with different right-hand sides.  Belos also
  /// implements a Block GMRES - based seed solver,
  /// Belos::GmresPolySolMgr.
  ///
  /// Users must ensure that each linear system has the same coefficient
  /// matrix.  The seed space is invariant during an individual linear
  /// system solve.  Finally, due to finite precision arithmetic, the
  /// off-diaognal "P'AP" terms grow.
  ///
  /// One often sees PCPG in context with the FETI domain
  /// decomposition method.
  ///
  /// \example PCPG/PCPGEpetraExFile.cpp
  ///
  /// The provided example uses PCPGSolMgr with an ML preconditioner.

  // Partial specialization for complex ScalarType.
  // This contains a trivial implementation.
  // See discussion in the class documentation above.
  //
  // FIXME (mfh 09 Sep 2015) This also is a stub for types other than
  // float or double.
  template<class ScalarType, class MV, class OP,
           const bool supportsScalarType =
             Belos::Details::LapackSupportsScalar<ScalarType>::value &&
             ! Teuchos::ScalarTraits<ScalarType>::isComplex>
  class PCPGSolMgr :
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
    PCPGSolMgr () :
      base_type ()
    {}
    PCPGSolMgr (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                const Teuchos::RCP<Teuchos::ParameterList> &pl) :
      base_type ()
    {}
    virtual ~PCPGSolMgr () {}

    //! clone for Inverted Injection (DII)
    Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const override {
      return Teuchos::rcp(new PCPGSolMgr<ScalarType,MV,OP,supportsScalarType>);
    }
  };

  template<class ScalarType, class MV, class OP>
  class PCPGSolMgr<ScalarType, MV, OP, true> :
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

    /*! \brief Empty constructor for PCPGSolMgr.
     * This constructor takes no arguments and sets the default values for the solver.
     * The linear problem must be passed in using setProblem() before solve() is called on this object.
     * In most instances, LinearProblem setProblem(...) methods are used.
     * Solver values may be changed using setParameters().
     */
     PCPGSolMgr();

    /*! \brief Basic constructor for PCPGSolMgr.
     * The constructor accepts a LinearProblem to be solved and a parameter list of these options:
     *
     *   - "Num Deflated Blocks" - a \c int specifying the number of blocks deflated from the linear system. Default: 2
     *                     The parameter distinguishes PCPG from CG.
     *   - "Num Saved Blocks" - a \c int specifying the maximum number of blocks saved from old Krylov bases. Default: 16
     *                     The parameter distinguishes PCPG from CG.
     *   - "Block Size" - an \c int specifying the block size to be used by the underlying block
     *                    conjugate-gradient solver.  In PCPC block size = one.  Many parameters are
     *                    meaningless in the unit block size case.  Default: 1
     *   - "Adaptive Block Size" - a \c bool specifying whether the block size can be modified
     *                             throughout the solve. Default: true
     *                             Meaningless with unit block size
     *   - "Maximum Iterations" - an \c int specifying the maximum number of iterations the
     *                            underlying solver is allowed to perform. Default: 1000
     *   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms
     *                               must reach to decide convergence. Default: 1e-8.
     *   - "Orthogonalization" - a \c string specifying the desired orthogonalization:  DGKS, ICGS, IMGS. Default: "DGKS"
     *                           Meaningless with unit block size
     *   - "Orthogonalization Constant" - a \c MagnitudeType used by DGKS orthogonalization to
     *                                    determine whether another step of classical Gram-Schmidt
     *                                    is necessary.  Default: -1 (use DGKS default)
     *                                    Meaningless with unit block size
     *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Belos::Errors
     *   - "Output Style" - a OutputType specifying the style of output. Default: Belos::General
     *   - "Output Stream" - a reference-counted pointer to the output stream where all
     *                       solver output is sent.  Default: Teuchos::rcp(&std::cout,false)
     *   - "Output Frequency" - an \c int specifying how often convergence information should be
     *                          outputted.  Default: -1 (never)
     *   - "Show Maximum Residual Norm Only" - a \c bool specifying whether that only the maximum
     *                                         relative residual norm is printed if convergence
     *                                         information is printed. Default: false
     *                                         Meaningless with unit block size
     *   - "Timer Label" - a \c std::string to use as a prefix for the timer labels.  Default: "Belos"
     */
    PCPGSolMgr( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                      const Teuchos::RCP<Teuchos::ParameterList> &pl );

    //! Destructor.
    virtual ~PCPGSolMgr() {};

    //! clone for Inverted Injection (DII)
    virtual Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const {
      return Teuchos::rcp(new PCPGSolMgr<ScalarType,MV,OP>);
    }
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
      return Teuchos::tuple(timerSolve_);
    }

    /// \brief Tolerance achieved by the last \c solve() invocation.
    ///
    /// This is the maximum over all right-hand sides' achieved
    /// convergence tolerances, and is set whether or not the solve
    /// actually managed to achieve the desired convergence tolerance.
    MagnitudeType achievedTol() const {
      return achievedTol_;
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

    /*! \brief The method either solves the problem or decides to quit.  On each call, a (possibly null)
     * seed space is used to accelerate convergence.
     *
     * The method calls PCPGIter::iterate(), which will return either because a specially constructed status
     * test evaluates to ::Passed or an exception is thrown.  The first Krylov vectors are appended to the
     * seed space.
     *
     * A return from PCPGIter::iterate() signifies one of the following scenarios:
     *    - the maximum number of restarts has been exceeded. In this scenario, the current solutions to the linear system
     *      will be placed in the linear problem and return ::Unconverged.
     *    - global convergence has been met. In this case, the current solutions to the linear system will be
     *      placed in the linear problem and the solver manager will return ::Converged
     *
     * \returns ::ReturnType specifying:
     *     - ::Converged: the linear problem was solved to the specification required by the solver manager.
     *     - ::Unconverged: the linear problem was not solved to the specification desired by the solver manager.
     */
    ReturnType solve();

    //@}

    /** \name Overridden from Teuchos::Describable */
    //@{

    /** \brief Method to return description of the PCPG solver manager */
    std::string description() const;

    //@}

  private:

    // In the A-inner product, perform an RRQR decomposition without using A unless absolutely necessary.  Given
    // the seed space U and C = A U, find U1 and C1 with span(U1)=span(U) such that C1'U1 = I maintaining C=AU.
    int ARRQR(int numVecs, int numOrthVecs, const Teuchos::SerialDenseMatrix<int,ScalarType>& D);

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
    static constexpr int maxIters_default_ = 1000;
    static constexpr int deflatedBlocks_default_ = 2;
    static constexpr int savedBlocks_default_ = 16;
    static constexpr int verbosity_default_ = Belos::Errors;
    static constexpr int outputStyle_default_ = Belos::General;
    static constexpr int outputFreq_default_ = -1;
    static constexpr const char * label_default_ = "Belos";
    static constexpr const char * orthoType_default_ = "DGKS";
    static constexpr std::ostream * outputStream_default_ = &std::cout;

    //
    // Current solver values.
    //

    //! Convergence tolerance (read from parameter list).
    MagnitudeType convtol_;

    //! Orthogonalization parameter (read from parameter list).
    MagnitudeType orthoKappa_;

    //! Tolerance achieved by the last \c solve() invocation.
    MagnitudeType achievedTol_;

    //! Number of iterations taken by the last \c solve() invocation.
    int numIters_;

    //! Maximum iteration count (read from parameter list).
    int maxIters_;

    int deflatedBlocks_, savedBlocks_, verbosity_, outputStyle_, outputFreq_;
    std::string orthoType_;

    // Recycled subspace, its image and the residual
    Teuchos::RCP<MV> U_, C_, R_;

    // Actual dimension of current recycling subspace (<= savedBlocks_ )
    int dimU_;

    // Timers.
    std::string label_;
    Teuchos::RCP<Teuchos::Time> timerSolve_;

    // Internal state variables.
    bool isSet_;
  };


// Empty Constructor
template<class ScalarType, class MV, class OP>
PCPGSolMgr<ScalarType,MV,OP,true>::PCPGSolMgr() :
  outputStream_(Teuchos::rcp(outputStream_default_,false)),
  convtol_(DefaultSolverParameters::convTol),
  orthoKappa_(DefaultSolverParameters::orthoKappa),
  achievedTol_(Teuchos::ScalarTraits<MagnitudeType>::zero()),
  numIters_(0),
  maxIters_(maxIters_default_),
  deflatedBlocks_(deflatedBlocks_default_),
  savedBlocks_(savedBlocks_default_),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  orthoType_(orthoType_default_),
  dimU_(0),
  label_(label_default_),
  isSet_(false)
{}


// Basic Constructor
template<class ScalarType, class MV, class OP>
PCPGSolMgr<ScalarType,MV,OP,true>::PCPGSolMgr(
                                             const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                                             const Teuchos::RCP<Teuchos::ParameterList> &pl ) :
  problem_(problem),
  outputStream_(Teuchos::rcp(outputStream_default_,false)),

  convtol_(DefaultSolverParameters::convTol),
  orthoKappa_(DefaultSolverParameters::orthoKappa),
  achievedTol_(Teuchos::ScalarTraits<MagnitudeType>::zero()),
  numIters_(0),
  maxIters_(maxIters_default_),
  deflatedBlocks_(deflatedBlocks_default_),
  savedBlocks_(savedBlocks_default_),
  verbosity_(verbosity_default_),
  outputStyle_(outputStyle_default_),
  outputFreq_(outputFreq_default_),
  orthoType_(orthoType_default_),
  dimU_(0),
  label_(label_default_),
  isSet_(false)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    problem_.is_null (), std::invalid_argument,
    "Belos::PCPGSolMgr two-argument constructor: "
    "'problem' is null.  You must supply a non-null Belos::LinearProblem "
    "instance when calling this constructor.");

  if (! pl.is_null ()) {
    // Set the parameters using the list that was passed in.
    setParameters (pl);
  }
}


template<class ScalarType, class MV, class OP>
void PCPGSolMgr<ScalarType,MV,OP,true>::setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params )
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

  // Check for the maximum numbers of saved and deflated blocks.
  if (params->isParameter("Num Saved Blocks")) {
    savedBlocks_ = params->get("Num Saved Blocks",savedBlocks_default_);
    TEUCHOS_TEST_FOR_EXCEPTION(savedBlocks_ <= 0, std::invalid_argument,
                       "Belos::PCPGSolMgr: \"Num Saved Blocks\" must be strictly positive.");

    // savedBlocks > number of matrix rows and columns, not known in parameters.
    //TEUCHOS_TEST_FOR_EXCEPTION(savedBlocks_ >= maxIters_, std::invalid_argument,
    //"Belos::PCPGSolMgr: \"Num Saved Blocks\" must be less than \"Maximum Iterations\".");

    // Update parameter in our list.
    params_->set("Num Saved Blocks", savedBlocks_);
  }
  if (params->isParameter("Num Deflated Blocks")) {
    deflatedBlocks_ = params->get("Num Deflated Blocks",deflatedBlocks_default_);
    TEUCHOS_TEST_FOR_EXCEPTION(deflatedBlocks_ < 0, std::invalid_argument,
                       "Belos::PCPGSolMgr: \"Num Deflated Blocks\" must be positive.");

    TEUCHOS_TEST_FOR_EXCEPTION(deflatedBlocks_ > savedBlocks_, std::invalid_argument,
                       "Belos::PCPGSolMgr: \"Num Deflated Blocks\" must be <= \"Num Saved Blocks\".");

    // Update parameter in our list.
    // The static_cast is for clang link issues with the constexpr before c++17
    params_->set("Num Deflated Blocks", static_cast<int>(deflatedBlocks_));
  }

  // Check to see if the timer label changed.
  if (params->isParameter("Timer Label")) {
    std::string tempLabel = params->get("Timer Label", label_default_);

    // Update parameter in our list and solver timer
    if (tempLabel != label_) {
      label_ = tempLabel;
      params_->set("Timer Label", label_);
      std::string solveLabel = label_ + ": PCPGSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerSolve_ = Teuchos::TimeMonitor::getNewCounter(solveLabel);
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
                        "Belos::PCPGSolMgr: \"Orthogonalization\" must be either \"DGKS\", \"ICGS\", or \"IMGS\".");
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

  // Create status tests if we need to.

  // Basic test checks maximum iterations and native residual.
  if (maxIterTest_ == Teuchos::null)
    maxIterTest_ = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>( maxIters_ ) );

  if (convTest_ == Teuchos::null)
    convTest_ = Teuchos::rcp( new StatusTestResNorm_t( convtol_, 1 ) );

  sTest_ = Teuchos::rcp( new StatusTestCombo_t( StatusTestCombo_t::OR, maxIterTest_, convTest_ ) );

  // Create the status test output class.
  // This class manages and formats the output from the status test.
  StatusTestOutputFactory<ScalarType,MV,OP> stoFactory( outputStyle_ );
  outputTest_ = stoFactory.create( printer_, sTest_, outputFreq_, Passed+Failed+Undefined );

  // Set the solver string for the output test
  std::string solverDesc = " PCPG ";
  outputTest_->setSolverDesc( solverDesc );


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
                         "Belos::PCPGSolMgr(): Invalid orthogonalization type.");
    }
  }

  // Create the timer if we need to.
  if (timerSolve_ == Teuchos::null) {
    std::string solveLabel = label_ + ": PCPGSolMgr total solve time";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    timerSolve_ = Teuchos::TimeMonitor::getNewCounter(solveLabel);
#endif
  }

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}


template<class ScalarType, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList>
PCPGSolMgr<ScalarType,MV,OP,true>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    // Set all the valid parameters and their default values.
    pl->set("Convergence Tolerance", static_cast<MagnitudeType>(DefaultSolverParameters::convTol),
      "The relative residual tolerance that needs to be achieved by the\n"
      "iterative solver in order for the linear system to be declared converged.");
    pl->set("Maximum Iterations", static_cast<int>(maxIters_default_),
      "The maximum number of iterations allowed for each\n"
      "set of RHS solved.");
    pl->set("Num Deflated Blocks", static_cast<int>(deflatedBlocks_default_),
      "The maximum number of vectors in the seed subspace." );
    pl->set("Num Saved Blocks", static_cast<int>(savedBlocks_default_),
      "The maximum number of vectors saved from old Krylov subspaces." );
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
    pl->set("Timer Label", static_cast<const char *>(label_default_),
      "The string to use as a prefix for the timer labels.");
    pl->set("Orthogonalization", static_cast<const char *>(orthoType_default_),
      "The type of orthogonalization to use: DGKS, ICGS, IMGS");
    pl->set("Orthogonalization Constant",static_cast<MagnitudeType>(DefaultSolverParameters::orthoKappa),
      "The constant used by DGKS orthogonalization to determine\n"
      "whether another step of classical Gram-Schmidt is necessary.");
    validPL = pl;
  }
  return validPL;
}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType PCPGSolMgr<ScalarType,MV,OP,true>::solve() {

  // Set the current parameters if are not set already.
  if (!isSet_) { setParameters( params_ ); }

  Teuchos::BLAS<int,ScalarType> blas;
  Teuchos::LAPACK<int,ScalarType> lapack;
  ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

  TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null,PCPGSolMgrLinearProblemFailure,
                     "Belos::PCPGSolMgr::solve(): Linear problem is not a valid object.");

  TEUCHOS_TEST_FOR_EXCEPTION(!problem_->isProblemSet(),PCPGSolMgrLinearProblemFailure,
                     "Belos::PCPGSolMgr::solve(): Linear problem is not ready, setProblem() has not been called.");

  // Create indices for the linear systems to be solved.
  int numRHS2Solve = MVT::GetNumberVecs( *(problem_->getRHS()) );
  std::vector<int> currIdx(1);
  currIdx[0] = 0;

   bool debug = false;

  // Inform the linear problem of the current linear system to solve.
  problem_->setLSIndex( currIdx ); // block size == 1

  // Assume convergence is achieved, then let any failed convergence set this to false.
  bool isConverged = true;

  //////////////////////////////////////////////////////////////////////////////////////
  // PCPG iteration parameter list
  Teuchos::ParameterList plist;
  plist.set("Saved Blocks", savedBlocks_);
  plist.set("Block Size", 1);
  plist.set("Keep Diagonal", true);
  plist.set("Initialize Diagonal", true);

  //////////////////////////////////////////////////////////////////////////////////////
  // PCPG solver

  Teuchos::RCP<PCPGIter<ScalarType,MV,OP> > pcpg_iter;
  pcpg_iter = Teuchos::rcp( new PCPGIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,ortho_,plist) );
  // Number of iterations required to generate initial recycle space (if needed)

  // Enter solve() iterations
  {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*timerSolve_);
#endif
    while ( numRHS2Solve > 0 ) {  // test for quick return

      // Reset the status test.
      outputTest_->reset();

      // Create the first block in the current Krylov basis (residual).
      if (R_ == Teuchos::null)
        R_ = MVT::Clone( *(problem_->getRHS()), 1 );

      problem_->computeCurrResVec( &*R_ );


      // Hypothesis: if U_ is not null, then neither is C_ and furthermore U'C= I.
      // TODO: ensure hypothesis right here ... I have to think about use cases.

      if( U_ != Teuchos::null ){
        // Hypothesis: if U_ is not null, then neither is C_ and furthermore U'C= I.

        // possibly over solved equation ...  I want residual norms
        // relative to the initial residual, not what I am about to compute.
        Teuchos::RCP<MV> cur_soln_vec = problem_->getCurrLHSVec();
        std::vector<MagnitudeType> rnorm0(1);
        MVT::MvNorm( *R_, rnorm0 ); // rnorm0  = norm(R_);

        // Z := U_'*R_; xo += U_*Z ;R_ -= C_*Z
        std::cout  << "Solver Manager:  dimU_ = " << dimU_ << std::endl;
        Teuchos::SerialDenseMatrix<int,ScalarType> Z( dimU_, 1 );

        Teuchos::RCP<const MV> Uactive, Cactive;
        std::vector<int> active_columns( dimU_ );
        for (int i=0; i < dimU_; ++i) active_columns[i] = i;
        Uactive = MVT::CloneView(*U_, active_columns);
        Cactive = MVT::CloneView(*C_, active_columns);

        if( debug ){
          std::cout << " Solver Manager : check duality of seed basis " << std::endl;
          Teuchos::SerialDenseMatrix<int,ScalarType> H( dimU_, dimU_ );
          MVT::MvTransMv( one, *Uactive, *Cactive, H );
          H.print( std::cout );
        }

        MVT::MvTransMv( one, *Uactive, *R_, Z );
        Teuchos::RCP<MV> tempU = MVT::Clone( *R_, 1 );
        MVT::MvTimesMatAddMv( one, *Uactive, Z, zero, *tempU );  // UZ
        MVT::MvAddMv( one, *tempU, one, *cur_soln_vec, *cur_soln_vec );  // xo += tmp;
        MVT::MvTimesMatAddMv( one, *Cactive, Z, zero, *tempU );  // CZ
        MVT::MvAddMv( -one, *tempU, one, *R_, *R_ );  // R_ -= tmp;
        std::vector<MagnitudeType> rnorm(1);
        MVT::MvNorm( *R_, rnorm );
        if( rnorm[0] < rnorm0[0] * .001 ){  //reorthogonalize
          MVT::MvTransMv( one, *Uactive, *R_, Z );
          MVT::MvTimesMatAddMv( one, *Uactive, Z, zero, *tempU );
          MVT::MvAddMv( one, *tempU, one, *cur_soln_vec, *cur_soln_vec );  // xo += UZ;
          MVT::MvTimesMatAddMv( one, *Cactive, Z, zero, *tempU );
          MVT::MvAddMv( -one, *tempU, one, *R_, *R_ );  // R_ -= CZ;
        }
        Uactive = Teuchos::null;
        Cactive = Teuchos::null;
        tempU = Teuchos::null;
      }
      else {
        dimU_ = 0;
      }


      // Set the new state and initialize the solver.
      PCPGIterState<ScalarType,MV> pcpgState; // fails if R == null.

      pcpgState.R = R_;
      if( U_ != Teuchos::null ) pcpgState.U = U_;
      if( C_ != Teuchos::null ) pcpgState.C = C_;
      if( dimU_ > 0 ) pcpgState.curDim = dimU_;
      pcpg_iter->initialize(pcpgState);

      // treat initialize() exceptions here?  how to use try-catch-throw? DMD

      // Get the current number of deflated blocks with the PCPG iteration
      dimU_ = pcpgState.curDim;
      if( !dimU_ ) printer_->stream(Debug) << " No recycled subspace available for RHS index " << currIdx[0] << std::endl << std::endl;
      pcpg_iter->resetNumIters();

      if( dimU_ > savedBlocks_ )
        std::cout << "Error: dimU_  = " << dimU_ << " > savedBlocks_ = " << savedBlocks_ << std::endl;

      while(1) { // dummy loop for break

        // tell pcpg_iter to iterate
        try {
          if( debug ) printf("********** Calling iterate...\n");
          pcpg_iter->iterate();

          ////////////////////////////////////////////////////////////////////////////////////
          //
          // check convergence first
          //
          ////////////////////////////////////////////////////////////////////////////////////
          if ( convTest_->getStatus() == Passed ) {
            // we have convergence
            break;  // break from while(1){pcpg_iter->iterate()}
          }
          ////////////////////////////////////////////////////////////////////////////////////
          //
          // check for maximum iterations
          //
          ////////////////////////////////////////////////////////////////////////////////////
          else if ( maxIterTest_->getStatus() == Passed ) {
            // we don't have convergence
            isConverged = false;
            break;  // break from while(1){pcpg_iter->iterate()}
          }
          else {

          ////////////////////////////////////////////////////////////////////////////////////
          //
          // we returned from iterate(), but none of our status tests Passed.
          // Something is wrong, and it is probably the developers fault.
          //
          ////////////////////////////////////////////////////////////////////////////////////

            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                               "Belos::PCPGSolMgr::solve(): Invalid return from PCPGIter::iterate().");
          } // end if
        } // end try
        catch (const PCPGIterOrthoFailure &e) {

          // Check to see if the most recent solution yielded convergence.
          sTest_->checkStatus( &*pcpg_iter );
          if (convTest_->getStatus() != Passed)
            isConverged = false;
          break;
        }
        catch (const std::exception &e) {
          printer_->stream(Errors) << "Error! Caught exception in PCPGIter::iterate() at iteration "
                                   << pcpg_iter->getNumIters() << std::endl
                                   << e.what() << std::endl;
          throw;
        }
      } // end of while(1)

      // Update the linear problem.
      Teuchos::RCP<MV> update = pcpg_iter->getCurrentUpdate();
      problem_->updateSolution( update, true );

      // Inform the linear problem that we are finished with this block linear system.
      problem_->setCurrLS();

      // Get the state.   How did pcpgState die?
      PCPGIterState<ScalarType,MV> oldState = pcpg_iter->getState();

      dimU_ = oldState.curDim;
      int q = oldState.prevUdim;

      std::cout << "SolverManager: dimU_ " << dimU_ << "   prevUdim= " << q << std::endl;

      if( q > deflatedBlocks_ )
        std::cout << "SolverManager: Error deflatedBlocks = " << deflatedBlocks_ << std::endl;

      int rank;
      if( dimU_ > q ){ // Orthogonalize [U;C](:,prevUdim:dimU_)
        //Given the seed space U and C = A U for some symmetric positive definite A,
        //find U1 and C1 with span(U1)=span(U) such that C1'U1 = I maintaining C=AU

        //oldState.D->print( std::cout ); D = diag( C'*U )

        U_ = oldState.U; //MVT::MvPrint( *U_, std::cout );
        C_ = oldState.C; //MVT::MvPrint( *C_, std::cout );
        rank = ARRQR(dimU_,q, *oldState.D );
        if( rank < dimU_ ) {
                std::cout << " rank decreased in ARRQR, something to do? " << std::endl;
        }
        dimU_ = rank;

      } // Now U_ and C_ = AU are dual bases.

      if( dimU_ > deflatedBlocks_ ){

        if( !deflatedBlocks_ ){
           U_ = Teuchos::null;
           C_ = Teuchos::null;
           dimU_ = deflatedBlocks_;
           break;
        }

        bool Harmonic = false; // (Harmonic) Ritz vectors

        Teuchos::RCP<MV> Uorth;

        std::vector<int> active_cols( dimU_ );
        for (int i=0; i < dimU_; ++i) active_cols[i] = i;

        if( Harmonic ){
          Uorth = MVT::CloneCopy(*C_, active_cols);
        }
        else{
          Uorth = MVT::CloneCopy(*U_, active_cols);
        }

        // Explicitly construct Q and R factors
        Teuchos::SerialDenseMatrix<int,ScalarType> R(dimU_,dimU_);
        rank = ortho_->normalize(*Uorth, Teuchos::rcp(&R,false));
        Uorth = Teuchos::null;
        // TODO:  During the previous solve, the matrix that normalizes U(1:q) was computed and discarded.
        // One might save it, reuse it here, and just normalize columns U(q+1:dimU_) here.

        // throw an error if U is both A-orthonormal and rank deficient
        TEUCHOS_TEST_FOR_EXCEPTION(rank < dimU_,PCPGSolMgrOrthoFailure,
                           "Belos::PCPGSolMgr::solve(): Failed to compute orthonormal basis for initial recycled subspace.");


        // R VT' = Ur S,
        Teuchos::SerialDenseMatrix<int,ScalarType> VT; // Not referenced
        Teuchos::SerialDenseMatrix<int,ScalarType> Ur; // Not referenced
        int lwork = 5*dimU_;                           // minimal, extra computation < 67*dimU_
        int info = 0;  // Hermite
        int lrwork = 1;
        if( problem_->isHermitian() ) lrwork = dimU_;
        std::vector<ScalarType> work(lwork); //
        std::vector<ScalarType> Svec(dimU_); //
        std::vector<ScalarType> rwork(lrwork);
        lapack.GESVD('N', 'O',
                   R.numRows(),R.numCols(),R.values(), R.numRows(),
                   &Svec[0],
                   Ur.values(),1,
                   VT.values(),1, // Output: VT stored in R
                   &work[0], lwork,
                   &rwork[0], &info);

        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, PCPGSolMgrLAPACKFailure,
                             "Belos::PCPGSolMgr::solve(): LAPACK _GESVD failed to compute singular values.");

        if( work[0] !=  67. * dimU_ )
           std::cout << " SVD " << dimU_ <<  " lwork " << work[0]  << std::endl;
        for( int i=0; i< dimU_; i++)
           std::cout << i << " " << Svec[i] << std::endl;

        Teuchos::SerialDenseMatrix<int,ScalarType> wholeV( R, Teuchos::TRANS);

        int startRow = 0, startCol = 0;
        if( Harmonic )
          startCol = dimU_ - deflatedBlocks_;

        Teuchos::SerialDenseMatrix<int,ScalarType> V(Teuchos::Copy,
                                                     wholeV,
                                                     wholeV.numRows(),
                                                     deflatedBlocks_,
                                                     startRow,
                                                     startCol);
        std::vector<int> active_columns( dimU_ );
        std::vector<int> def_cols( deflatedBlocks_ );
        for (int i=0; i < dimU_; ++i) active_columns[i] = i;
        for (int i=0; i < deflatedBlocks_; ++i) def_cols[i] = i;

        Teuchos::RCP<MV> Uactive = MVT::CloneViewNonConst(*U_, def_cols);
        Teuchos::RCP<MV> Ucopy = MVT::CloneCopy( *U_, active_columns );
        MVT::MvTimesMatAddMv( one, *Ucopy, V, zero, *Uactive ); //  U:= U*V
        Ucopy   = Teuchos::null;
        Uactive = Teuchos::null;
        Teuchos::RCP<MV> Cactive = MVT::CloneViewNonConst(*C_, def_cols);
        Teuchos::RCP<MV> Ccopy = MVT::CloneCopy( *C_, active_columns );
        MVT::MvTimesMatAddMv( one, *Ccopy, V, zero, *Cactive ); //  C:= C*V
        Ccopy  = Teuchos::null;
        Cactive = Teuchos::null;
        dimU_ = deflatedBlocks_;
      }
      printer_->stream(Debug) << " Generated recycled subspace using RHS index " << currIdx[0] << " of dimension " << dimU_ << std::endl << std::endl;

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

  // Save the convergence test value ("achieved tolerance") for this solve.
  {
    using Teuchos::rcp_dynamic_cast;
    typedef StatusTestGenResNorm<ScalarType,MV,OP> conv_test_type;
    // testValues is nonnull and not persistent.
    const std::vector<MagnitudeType>* pTestValues =
      rcp_dynamic_cast<conv_test_type>(convTest_)->getTestValue();

    TEUCHOS_TEST_FOR_EXCEPTION(pTestValues == NULL, std::logic_error,
      "Belos::PCPGSolMgr::solve(): The convergence test's getTestValue() "
      "method returned NULL.  Please report this bug to the Belos developers.");

    TEUCHOS_TEST_FOR_EXCEPTION(pTestValues->size() < 1, std::logic_error,
      "Belos::PCPGSolMgr::solve(): The convergence test's getTestValue() "
      "method returned a vector of length zero.  Please report this bug to the "
      "Belos developers.");

    // FIXME (mfh 12 Dec 2011) Does pTestValues really contain the
    // achieved tolerances for all vectors in the current solve(), or
    // just for the vectors from the last deflation?
    achievedTol_ = *std::max_element (pTestValues->begin(), pTestValues->end());
  }

  // get iteration information for this solve
  numIters_ = maxIterTest_->getNumIters();

  if (!isConverged) {
    return Unconverged; // return from PCPGSolMgr::solve()
  }
  return Converged; // return from PCPGSolMgr::solve()
}

// A-orthogonalize the Seed Space
// Note that Anasazi::GenOrthoManager provides simplified versions of the algorithm,
// that are not rank revealing, and are not designed for PCPG in other ways too.
template<class ScalarType, class MV, class OP>
int PCPGSolMgr<ScalarType,MV,OP,true>::ARRQR(int p, int q, const Teuchos::SerialDenseMatrix<int,ScalarType>& D)
{
  using Teuchos::RCP;
  ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

  // Allocate memory for scalars.
  Teuchos::SerialDenseMatrix<int,ScalarType> alpha( 1, 1 );
  Teuchos::SerialDenseMatrix<int,ScalarType> gamma( 1, 1 );
  Teuchos::SerialDenseMatrix<int,ScalarType> anorm( 1, 1 );
  std::vector<int> curind(1);
  std::vector<int> ipiv(p - q); // RRQR Pivot indices
  std::vector<ScalarType> Pivots(p); // RRQR Pivots
  int i, imax, j, l;
  ScalarType rteps = 1.5e-8;

  // Scale such that diag( U'C) = I
  for( i = q ; i < p ; i++ ){
    ipiv[i-q] = i;
    curind[0] = i;
    RCP<MV> P = MVT::CloneViewNonConst(*U_,curind);
    RCP<MV> AP = MVT::CloneViewNonConst(*C_,curind);
    anorm(0,0) = one / Teuchos::ScalarTraits<ScalarType>::squareroot( D(i-q,i-q) ) ;
    MVT::MvAddMv( anorm(0,0), *P, zero, *AP, *P );
    MVT::MvAddMv( zero, *P, anorm(0,0), *AP, *AP );
    Pivots[i]  = one;
  }

  for( i = q ; i < p ; i++ ){
    if( q < i && i < p-1 ){ // Find the largest pivot
      imax = i;
      l = ipiv[imax-q];
      for( j = i+1 ; j < p ; j++ ){
         const int k = ipiv[j-q];
         if( Pivots[k] > Pivots[l] ){
           imax = j;
           l = k;
         }
      } // end for
      if( imax > i ){
          l = ipiv[imax-q]; // swap ipiv( imax ) and ipiv(i+1)
          ipiv[imax-q] = ipiv[i-q];
          ipiv[i-q] = l;
      }
    } // largest pivot found
    int k = ipiv[i-q];

    if( Pivots[k]  > 1.5625e-2 ){
      anorm(0,0) =  Pivots[k]; // A-norm of u
    }
    else{ // anorm(0,0) = sqrt( U(:,k)'*C(:,k) );
      curind[0] = k;
      RCP<const MV> P = MVT::CloneView(*U_,curind);
      RCP<const MV> AP = MVT::CloneView(*C_,curind);
      MVT::MvTransMv( one, *P, *AP, anorm );
      anorm(0,0) = Teuchos::ScalarTraits<ScalarType>::squareroot( anorm(0,0) ) ;
    }
    if( rteps <= anorm(0,0) && anorm(0,0) < 9.765625e-4){
       /*
       C(:,k) = A*U(:,k);  % Change C
       fixC = U(:, ipiv(1:i-1) )'*C(:,k);
       U(:,k) = U(:,k) - U(:, ipiv(1:i-1) )*fixC;
       C(:,k) = C(:,k) - C(:, ipiv(1:i-1) )*fixC;
       anorm = sqrt( U(:,k)'*C(:,k) );
       */
       std::cout << "ARRQR: Bad case not implemented" << std::endl;
    }
    if( anorm(0,0) < rteps ){ // rank [U;C] = i-1
       std::cout << "ARRQR : deficient case not implemented " << std::endl;
       //U = U(:, ipiv(1:i-1) );
       //C = C(:, ipiv(1:i-1) );
       p = q + i;
       // update curDim_ in State
       break;
    }
    curind[0] = k;
    RCP<MV> P = MVT::CloneViewNonConst(*U_,curind);
    RCP<MV> AP = MVT::CloneViewNonConst(*C_,curind);
    MVT::MvAddMv( anorm(0,0), *P, zero, *AP, *P ); // U(:,k) = U(:,k)/anorm;
    MVT::MvAddMv( zero, *P, anorm(0,0), *AP, *AP ); // C(:,k) = C(:,k)/anorm;
    P = Teuchos::null;
    AP = Teuchos::null;
    Pivots[k] = one;                 // delete,  for diagonostic purposes
    P = MVT::CloneViewNonConst(*U_,curind);  // U(:,k)
    AP = MVT::CloneViewNonConst(*C_,curind); // C(:,k)
    for( j = i+1 ; j < p ; j++ ){
      l = ipiv[j-q];   // ahhh
      curind[0] = l;
      RCP<MV> Q = MVT::CloneViewNonConst(*U_,curind); // segmentation fault,  j=i+1=5
      MVT::MvTransMv( one, *Q, *AP, alpha); // alpha(0,0) = U(:,l)'*C(:,k);
      MVT::MvAddMv( -alpha(0,0), *P, one, *Q, *Q ); // U(:,l) -= U(:,k) * alpha(0,0);
      Q = Teuchos::null;
      RCP<MV> AQ = MVT::CloneViewNonConst(*C_,curind);
      MVT::MvAddMv( -alpha(0,0), *AP, one, *AQ, *AQ ); // C(:,l) -= C(:,l) - C(:,k) * alpha(0,0);
      AQ = Teuchos::null;
      gamma(0,0) = ( Pivots[l] - alpha(0,0))*( Pivots[l] + alpha(0,0));
      if( gamma(0,0) > 0){
        Pivots[l] = Teuchos::ScalarTraits<ScalarType>::squareroot( gamma(0,0) );
      }
      else {
        Pivots[l] = zero; //rank deficiency revealed
      }
    }
  }
  return p;
}

//  The method returns a string describing the solver manager.
template<class ScalarType, class MV, class OP>
std::string PCPGSolMgr<ScalarType,MV,OP,true>::description() const
{
  std::ostringstream oss;
  oss << "Belos::PCPGSolMgr<...,"<<Teuchos::ScalarTraits<ScalarType>::name()<<">";
  oss << "{";
  oss << "Ortho Type='"<<orthoType_;
  oss << "}";
  return oss.str();
}

} // end Belos namespace

#endif /* BELOS_PCPG_SOLMGR_HPP */
