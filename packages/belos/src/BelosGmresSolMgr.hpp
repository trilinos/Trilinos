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

#ifndef __Belos_GmresSolMgr_hpp
#define __Belos_GmresSolMgr_hpp

/// \file BelosGmresSolMgr.hpp
/// \author Mark Hoemmen
/// \brief Solver manager for CA-GMRES and standard GMRES.

#include <BelosGmresBaseIteration.hpp>

#include <BelosOutputManager.hpp>
#include <BelosSolverManager.hpp>
#include <BelosStatusTestGenResNorm.hpp>
#include <BelosStatusTestOutputFactory.hpp>
#include <BelosOrthoManagerFactory.hpp>

#include <functional> // std::logical_and
#include <numeric> // std::accumulate

namespace Belos {

  /// \class GmresSolMgr
  /// \author Mark Hoemmen
  /// \brief Solver manager for CA-GMRES and standard GMRES.
  ///
  template<class Scalar, class MV, class OP>
  class GmresSolMgr : public SolverManager<Scalar, MV, OP> {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    typedef MV multivector_type;
    typedef OP operator_type;

  private:
    typedef MultiVecTraits<Scalar,MV> MVT;
    typedef OperatorTraits<Scalar,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;
    typedef GmresBaseIteration<Scalar, MV, OP> iteration_type;

  public:
    //! \name Constructors and destructor 
    //@{ 

    /// \brief Preferred constructor.
    ///
    /// \param problem [in/out] The linear problem to solve.
    /// \param params [in] Parameters for the solve.  If null, 
    ///   we use defaults, else we make a deep copy.
    GmresSolMgr (const Teuchos::RCP<LinearProblem<Scalar,MV,OP> >& problem,
		 const Teuchos::RCP<const Teuchos::ParameterList>& params,
		 const bool debug = false) : 
      problem_ (validatedProblem (problem)), 
      debug_ (debug)
    { 
      setParametersImpl (params);
    }

    //! Default constructor; defers setting parameters.
    GmresSolMgr() 
    {
      throw std::logic_error("TODO: implement deferring setting parameters until solve()");
    }

    //! Destructor, defined virtual for safe inheritance.
    virtual ~GmresSolMgr() {}

    //@}

    //! \name "Get" methods
    //@{ 

    //! Const reference to the linear problem being solved.
    const LinearProblem<Scalar,MV,OP>& getProblem() const {
      return *problem_;
    }

    /// \brief Valid default parameters for this solver manager.
    ///
    /// This class method is preferred to the instance method, because
    /// you can call it before you've instantiated the solver manager,
    /// and pass the parameters into the solver manager's constructor.
    ///
    /// \warning This routine is not reentrant.  It caches the
    ///   default parameter list for later reuse.
    static Teuchos::RCP<const Teuchos::ParameterList> getDefaultParameters();

    /// \brief Valid default parameters for this solver manager.
    ///
    /// The class method \c getDefaultParameters() is preferred, since
    /// it may be called before the solver manager has been
    /// instantiated.
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const {
      return getDefaultParameters();
    }

    //! Current parameters for this solver manager instance.
    Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const {
      return params_;
    }

    /// \brief Iteration count for the most recent call to \c solve().
    ///
    /// It's not entirely clear what the SolutionManager interface
    /// expects this to mean.  We've interpreted it to mean the
    /// maximum (over all right-hand side(s)) of the total number of
    /// iterations over all restart cycle(s).
    int getNumIters() const { 
      if (iter_.is_null())
	throw std::logic_error("The number of iterations is undefined, since "
			       "the iteration has not yet been initialized.");
      return iter_->getNumIters();
    }

    /// \brief Was a loss of accuracy detected?
    ///
    /// Some GMRES-type solvers have the capability to detect loss of
    /// accuracy.  If this solver does not, this method always returns
    /// false.  If the solver does have this capability, this flag
    /// starts out false, and will be (re)set the next time \c solve()
    /// is called.
    bool isLOADetected() const {
      return false;
    }

    /// \brief Last solve's number of restarts and total iterations.
    ///
    /// After calling solve(), you may call this method.  For each
    /// right-hand side solved, it returns (total number of restart
    /// cycles, total number of iterations over all restart cycles).
    /// Before calling solve(), the result of calling this method is
    /// undefined.
    std::vector<std::pair<int, int> >
    totalNumIters () const {
      return totalNumIters_;
    }

    //@}

    //! \name "Set" methods
    //@{

    /// \brief Set the linear problem to solve.
    ///
    /// This method restarts the current solve that might be in progress.
    void 
    setProblem (const Teuchos::RCP<LinearProblem<Scalar,MV,OP> >& problem)
    {
      problem_ = validatedProblem (problem);
      // Changing the problem requires (re)initializing possibly a lot
      // of different things, including the status tests,
      // orthogonalization method, and the GmresBaseIteration object.
      reset (Belos::Problem);
    }

    /// \brief Set parameters for solving the linear problem.
    ///
    /// If necessary, this method restarts (by calling
    /// reset(Belos::Problem)) the current solve that might be in
    /// progress.  Currently, it does so only if the parameters have
    /// not yet been set, or if the maximum number of iterations has
    /// changed (which affects GMRES storage).
    ///
    /// \param params [in] New parameters for the linear solve.  The
    ///   original parameter list is not modified.  Since
    ///   setParameters() implements a pure virtual method of
    ///   Belos::SolutionManager, the ParameterList has to be passed
    ///   in as non-const, even though we don't modify it.
    ///
    /// \note We don't actually keep a pointer to params.  This is
    ///   because we might fill in unsupplied parameters with their
    ///   default values.  Instead, we make a deep copy.
    void 
    setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params);

    /// \brief Set a user-defined convergence stopping criterion.
    ///
    /// This test will be applied in short-circuiting Boolean AND
    /// sequence after the other convergence tests.  All convergence
    /// tests are invoked (in short-circuiting Boolean OR fashion)
    /// after testing whether the maximum number of iterations has
    /// been exceeded.
    ///
    /// \param userConvTest [in] The user-defined convergence stopping
    ///   criterion.
    void 
    setUserConvStatusTest (const Teuchos::RCP<StatusTest<Scalar,MV,OP> >& userConvTest);
    //@}

    //! \name "Reset" methods
    //@{

    /// \brief Reset the solver manager.
    ///
    /// "Reset" with type = Belos::Problem means the following:
    /// - Tell the LinearProblem instance to recompute initial
    ///   residual(s) (both left-preconditioned, if applicable, and
    ///   unpreconditioned) with its stored current approximate
    ///   solution
    /// - Restart GMRES, using the newly (re)computed residual
    ///   vector(s).
    ///
    /// This method does not currently support values of type for
    /// which (type & Belos::RecycleSubspace) != 0.
    ///
    /// \param type [in] The type of reset to perform.  Only type =
    ///   Belos::Problem is currently supported.
    void reset (const ResetType type);
    //@}

    //! \name Solver application methods
    //@{ 

    /// \brief Attempt to solve the linear system \f$AX=B\f$.
    ///
    /// Do so by calling the underlying linear solver's iterate()
    /// routine zero or more times, until the problem has been solved
    /// (as decided by the solver manager) or the solver manager
    /// decides to quit.
    ///
    /// \return ReturnType enum specifying Converged (the linear
    ///   problem was solved to the desired tolerance) or Unconverged
    ///   (the linear problem was not thusly solved).
    ReturnType solve();

    /// \brief Change the right-hand side of the linear system to solve.
    ///
    /// This works like a special case of setProblem(), when the only
    /// part of the problem that has changed is the right-hand side,
    /// and the new right-hand side is in the same vector space as the
    /// previous right-hand side.  This avoids possibly expensive
    /// reinitialization of things like the orthogonalization manager.
    /// The initial guess will be left the same; in fact, if solve()
    /// has been called before, the initial guess will be the
    /// approximate solution from the last invocation of solve().
    ///
    /// \note This does _not_ create a new LinearProblem instance.
    ///   The LinearProblem reference returned by getProblem() will
    ///   still be valid and will still point to the original
    ///   LinearProblem, except that the original LinearProblem's
    ///   right-hand side will be different.  The original right-hand
    ///   side will not be overwritten, though, so it will not go away
    ///   if you retain an RCP to it.
    void 
    setRHS (const Teuchos::RCP<const MV>& B) 
    {
      TEST_FOR_EXCEPTION(B.is_null(), std::invalid_argument,
			 "Belos::GmresSolMgr::setRHS(): the given new right-"
			 "hand side B is null.");
      // This "unsets" the problem: it makes the initial residual (and
      // the left-preconditioned initial residual, if applicable)
      // invalid.  We don't have to call setProblem() to reset the
      // problem, since setProblem() will be called in solve().
      problem_->setRHS (B);
    }

    /// \brief Change the initial guess of the linear system to solve.
    ///
    /// This works like a special case of setProblem(), when the only
    /// part of the problem that has changed is the initial guess, and
    /// the new initial guess is in the same vector space as the
    /// previous initial guess.  This avoids possibly expensive
    /// reinitialization of things like the orthogonalization manager.
    ///
    /// \note This does _not_ create a new LinearProblem instance.
    ///   The LinearProblem reference returned by getProblem() will
    ///   still be valid and will still point to the original
    ///   LinearProblem, except that the original LinearProblem's
    ///   initial guess will be different.  The original initial guess
    ///   will not be overwritten, though, so it will not go away if
    ///   you retain an RCP to it.
    void 
    setLHS (const Teuchos::RCP<MV>& X) 
    {
      TEST_FOR_EXCEPTION(X.is_null(), std::invalid_argument,
			 "Belos::GmresSolMgr::setLHS(): the given new initial "
			 "guess X is null.");
      // This "unsets" the problem; it makes the initial residual (and
      // the left-preconditioned initial residual, if applicable)
      // invalid.  We don't have to call setProblem() to reset the
      // problem, since setProblem() will be called in solve().
      problem_->setLHS (X);
    }

    /// Change stopping criteria for next invocation of solve().
    ///
    /// \param convTol [in] The new convergence tolerance.
    ///   Interpretation of this depends on how the residual norm
    ///   test(s) were initially set up.
    /// \param maxItersPerRestart [in] Maximum number of iterations
    ///   per restart cycle.
    /// \param maxNumRestarts [in] Maximum number of restart cycle(s).
    void
    changeStoppingCriteria (const magnitude_type convTol,
			    const int maxItersPerRestart,
			    const int maxNumRestarts)
    {
      // Testing maxNumRestarts now, rather than below, ensures the
      // strong exception guarantee: either no state will be changed,
      // or no exceptions will be thrown.
      TEST_FOR_EXCEPTION(maxNumRestarts < 0, std::invalid_argument,
			 "The maximum number of restart cycle(s) must be "
			 "nonnegative, but a value of " << maxNumRestarts 
			 << " was specified.");
      rebuildStatusTests (convTol, maxItersPerRestart);
      params_->set ("Maximum Restarts", maxNumRestarts,
		    "Maximum number of restart cycle(s) allowed for each "
		    "right-hand side solved.");
    }
    //@}
    
  private:
    //! The linear problem to solve.
    Teuchos::RCP<LinearProblem<Scalar, MV, OP> > problem_;

    //! Output manager.
    Teuchos::RCP<OutputManager<Scalar> > outMan_;

    //! Current parameters for this solver manager instance.
    Teuchos::RCP<Teuchos::ParameterList> params_;

    //! Convergence stopping criterion for the current solve.
    Teuchos::RCP<StatusTest<Scalar,MV,OP> > convTest_;

    //! Cumulative stopping criterion for the current solve.
    Teuchos::RCP<StatusTest<Scalar,MV,OP> > statusTest_;

    /// \brief User-defined custom convergence test.
    ///
    /// When the caller sets this, the current status test
    /// (statusTest_) is modified to include the user-defined test in
    /// sequential AND fashion.  Any previous user-defined status test
    /// is discarded.
    Teuchos::RCP<StatusTest<Scalar,MV,OP> > userConvTest_;

    //! "Status test" that outputs intermediate iteration results.
    Teuchos::RCP<StatusTestOutput<Scalar,MV,OP> > outTest_;

    //! The orthogonalization method to use for GMRES.
    Teuchos::RCP<const OrthoManager<Scalar, MV> > orthoMan_;

    //! Instance of the Belos::Iteration subclass.
    Teuchos::RCP<iteration_type> iter_;

    /// After calling solve(): For each right-hand side, the total
    /// number of restart cycle(s) (first element in the pair), and
    /// the total number of iterations over all the restart cycle(s)
    /// (second element in the pair).  Undefined if solve() has not
    /// yet been called.
    std::vector<std::pair<int, int> > totalNumIters_;

    //! Whether or not to print debug output.
    bool debug_;

    /// \brief Set parameters for solving the linear problem.
    ///
    /// If necessary, this method restarts (by calling
    /// reset(Belos::Problem)) the current solve that might be in
    /// progress.  Currently, it does so only if the parameters have
    /// not yet been set, or if the maximum number of iterations has
    /// changed (which affects GMRES storage).
    ///
    /// \param params [in] New parameters for the linear solve.
    ///   The original parameter list is not modified.
    void 
    setParametersImpl (const Teuchos::RCP<const Teuchos::ParameterList>& params);

    /// \brief (Re)build all the iteration stopping criteria.
    ///
    /// \param plist [in] If supplied, read the parameters for
    ///   constructing the stopping criteria from the given parameter
    ///   list.  Otherwise, use the solver manager's current parameter
    ///   list (params_).
    void 
    rebuildStatusTests (Teuchos::RCP<const Teuchos::ParameterList> plist = 
			Teuchos::null);

    /// \brief (Re)build all the iteration stopping criteria.
    ///
    /// Variant of the one-argument rebuildStatusTests(), when you
    /// just want to change the convergence tolerance and maximum
    /// number of iterations per restart cycle.  (The status tests
    /// don't control the maximum number of restarts.)
    ///
    /// \param convTol [in] The new convergence tolerance.
    ///   Interpretation of this depends on how the residual norm
    ///   test(s) were initially set up.
    /// \param maxItersPerRestart [in] Maximum number of iterations
    ///   per restart cycle.
    void
    rebuildStatusTests (const magnitude_type convTol,
			const int maxItersPerRestart);

    /// \brief Initialize the OrthoManager (orthogonalization method).
    ///
    /// \param plist [in] If supplied, read the parameters for
    ///   constructing the orthogonalization method from the given
    ///   parameter list.  Otherwise, use the solver manager's current
    ///   parameter list (params_).
    void
    rebuildOrthoManager (Teuchos::RCP<const Teuchos::ParameterList> plist = 
			 Teuchos::null);

    /// \brief (Re)build the Iteration subclass.
    ///
    /// \param plist [in] If supplied, read the parameters for
    ///   constructing the Iteration subclass from the given parameter
    ///   list.  Otherwise, use the solver manager's current parameter
    ///   list (params_).
    void 
    rebuildIteration (Teuchos::RCP<const Teuchos::ParameterList> plist = 
		      Teuchos::null);

    /// \brief Make sure that GmresSolMgr can solve the given problem.
    ///
    /// If it can, return a validated version of the problem (the same
    /// instance, perhaps with setProblem() called.  Otherwise, if the
    /// problemis invalid, throw an std::invalid_argument exception.
    ///
    /// \warning This method does _not_ check whether the
    ///   LinearProblem instance's setLSIndex() method has been called
    ///   in order to set the current linear system(s) to solve.
    ///   The GmresBaseIteration object cannot be (re)initialized 
    ///   until that is the case.
    static Teuchos::RCP<LinearProblem<Scalar, MV, OP> >
    validatedProblem (const Teuchos::RCP<LinearProblem<Scalar, MV, OP> >& problem)
    {
      const char prefix[] = "Belos::GmresSolMgr: ";
      TEST_FOR_EXCEPTION(problem.is_null(), std::invalid_argument,
			 prefix << "The linear problem (Belos::LinearProblem "
			 "instance) that you want me to solve is null.");
      if (! problem->isProblemSet())
	{
	  problem->setProblem();
	  const bool notInitialized = problem->getOperator().is_null() || 
	    problem->getRHS().is_null() || problem->getLHS().is_null();
	  if (notInitialized)
	    {
	      std::ostringstream os;
	      os << prefix << "The given linear problem (Belos::Linear"
		"Problem) instance is not fully initialized: "
		 << (problem->getOperator().is_null() ? "the operator A is null; " : "")
		 << (problem->getRHS().is_null() ? "the right-hand side B is null; " : "")
		 << (problem->getLHS().is_null() ? "the initial guess X is null" : "")
		 << ".";
	      TEST_FOR_EXCEPTION(notInitialized, std::invalid_argument, os.str());
	    }
	  TEST_FOR_EXCEPTION(! problem->isProblemSet(), std::logic_error,
			     prefix << "Although the given LinearProblem instance "
			     "has non-null operator (matrix A), right-hand side B,"
			     " and initial guess X, and although its setProblem() "
			     "method has been called, its isProblemSet() method "
			     "returns false.  This should not happen.");
	}
      return problem;
    }

    //! Convert string to enumerated type for residual test.
    static ScaleType
    stringToScaleType (const std::string& scaleType) 
    {
      if (scaleType == "Norm of Initial Residual") 
        return Belos::NormOfInitRes;
      else if (scaleType == "Norm of Preconditioned Initial Residual") 
        return Belos::NormOfPrecInitRes;
      else if (scaleType == "Norm of RHS" || scaleType == "Norm of Right-Hand Side")
        return Belos::NormOfRHS;
      else if (scaleType == "None")
        return Belos::None;
      else {
        TEST_FOR_EXCEPTION (true, std::logic_error,
			    "Invalid residual scaling type \"" << scaleType 
			    << "\".");
      }
    }

    /// Create a new "output status test" object.
    ///
    /// The "output status test" knows how to generate "progress
    /// reports" of convergence as the iterations progress.
    /// Prerequisites for creating this are the "output manager" and
    /// the stopping criterion ("status test").
    ///
    /// \note The output test is cheap enough to create that we just
    ///   recreate it each time we change the settings, rather than
    ///   trying to keep around its constituent components and change
    ///   their settings.
    static Teuchos::RCP<StatusTestOutput<Scalar,MV,OP> > 
    initOutputTest (const Teuchos::RCP<OutputManager<Scalar> >& outMan,
		    const Teuchos::RCP<StatusTest<Scalar,MV,OP> >& statusTest,
		    const OutputType outStyle,
		    const int outFreq,
		    const std::string& solverDesc = "",  // e.g., " CA-GMRES "
		    const std::string& precondDesc = ""); // Preconditioner description

    /// \brief Create a new convergence test.
    ///
    /// The convergence test is a prerequisite for creating the final
    /// stopping criterion (StatusTest) object, which checks both the
    /// number of iterations and convergence.
    //
    /// \param convTol [in] Relative residual convergence tolerance.
    /// \param haveLeftPreconditioner [in] Whether or not a left
    ///   preconditioner has been provided.
    /// \param implicitScaleType [in]
    /// \param explicitScaleType [in]
    /// \param userConvTest [in] If not null, a user-defined
    ///   convergence test that will be performed sequentially after
    ///   the maximum iteration count test and the implicit (and
    ///   explicit) residual norm tests.
    ///
    /// \return Convergence test stopping criterion, which stops if
    ///   both the relative residual error has dropped below the given
    ///   level, and whether the user-supplied convergence test (if
    ///   any) has passed.
    static Teuchos::RCP<StatusTest<Scalar,MV,OP> >
    initConvTest (const magnitude_type convTol,
		  const bool haveLeftPreconditioner,
		  const std::string& implicitScaleType,
		  const std::string& explicitScaleType,
		  Teuchos::RCP<StatusTest<Scalar,MV,OP> > userConvTest);

    /// Create a new stopping criterion ("StatusTest") object.
    ///
    /// The stopping criterion is a prerequisite for creating the
    /// "output status test" object.
    ///
    /// \param convTest [in] Convergence test, previously created by 
    ///   initConvTest.
    /// \param maxIters [in] Maximum number of iterations to execute.
    ///
    /// \return Stopping criterion which succeeds if the given maximum
    ///   iteration count has been reached, OR if the given
    ///   convergence test has passed.
    static Teuchos::RCP<StatusTest<Scalar,MV,OP> >
    initStatusTest (Teuchos::RCP<StatusTest<Scalar,MV,OP> > convTest,
		    const int maxIters);

    /// \brief Create or reinitialize the given output manager.
    ///
    /// If outMan is null, create and return a new OutputManager with
    /// the given verbosity, that reports to the given output stream.
    /// If outMan is not null, set the verbosity and output stream of
    /// the OutputManager, and return the (pointer to the)
    /// OutputManager.
    ///
    /// \param outMan [in/out] Either the output manager to reset, or
    ///   null (if a new output manager is to be created).
    ///
    /// \param verbosity [in] The new verbosity level.
    ///
    /// \param outStream [in/out] The output stream to which the
    ///   output manager should report.  If null, std::cout is used.
    ///
    /// \return Output manager with verbosity and output stream set to
    ///   the specified values.
    static Teuchos::RCP<OutputManager<Scalar> >
    initOutputManager (const Teuchos::RCP<OutputManager<Scalar> >& outMan,
		       const MsgType verbosity,
		       const Teuchos::RCP<std::ostream>& outStream);

  };


  template<class Scalar, class MV, class OP>
  Teuchos::RCP<OutputManager<Scalar> >
  GmresSolMgr<Scalar,MV,OP>::
  initOutputManager (const Teuchos::RCP<OutputManager<Scalar> >& outMan,
		     const MsgType verbosity,
		     const Teuchos::RCP<std::ostream>& outStream)
  {
    // If outStream is null, use std::cout as the default.
    Teuchos::RCP<std::ostream> out = 
      outStream.is_null() ? Teuchos::rcpFromRef(std::cout) : outStream;
    if (outMan.is_null())
      return Teuchos::rcp (new OutputManager<Scalar> (verbosity, out));
    else
      {
	outMan->setVerbosity (verbosity);
	outMan->setOStream (out);
	return outMan;
      }
  }

  template<class Scalar, class MV, class OP>
  Teuchos::RCP<StatusTestOutput<Scalar,MV,OP> > 
  GmresSolMgr<Scalar,MV,OP>::
  initOutputTest (const Teuchos::RCP<OutputManager<Scalar> >& outMan,
		  const Teuchos::RCP<StatusTest<Scalar,MV,OP> >& statusTest,
		  const OutputType outStyle,
		  const int outFreq,
		  const std::string& solverDesc,
		  const std::string& precondDesc)
  {
    using Teuchos::RCP;

    TEST_FOR_EXCEPTION(outMan.is_null(), std::logic_error,
			"Construction / reinitialization of the output status "
			"test depends on the OutputManager being initialized, "
			"but it has not yet been initialized.");
    TEST_FOR_EXCEPTION(statusTest.is_null(), std::logic_error,
		       "Construction / reinitialization of the output status "
		       "test depends on the status test being initialized, "
		       "but it has not yet been initialized.");
    StatusTestOutputFactory<Scalar,MV,OP> factory (outStyle);

    RCP<StatusTestOutput<Scalar, MV, OP> > newOutTest = 
      factory.create (outMan, statusTest, outFreq, 
		      Passed+Failed+Undefined);

    // newOutTest->setOutputManager (outMan);
    // newOutTest->setOutputFrequency (outFreq);
    // newOutTest->setChild (statusTest);

    // The factory doesn't know how to set the description strings,
    // so we do so here.
    if (solverDesc != "")
      newOutTest->setSolverDesc (solverDesc);
    if (precondDesc != "")
      newOutTest->setPrecondDesc (precondDesc);

    return newOutTest;
  }


  template<class Scalar, class MV, class OP>
  Teuchos::RCP<StatusTest<Scalar,MV,OP> >
  GmresSolMgr<Scalar,MV,OP>::
  initConvTest (const magnitude_type convTol,
		const bool haveLeftPreconditioner,
		const std::string& implicitScaleType,
		const std::string& explicitScaleType,
		Teuchos::RCP<StatusTest<Scalar,MV,OP> > userConvTest)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef StatusTest<Scalar,MV,OP> base_test;
    typedef StatusTestGenResNorm<Scalar,MV,OP> res_norm_test;
    typedef StatusTestCombo<Scalar,MV,OP> combo_test;

    // "Deflation Quorum": number of converged systems before
    // deflation is allowed.  Cannot be larger than "Block Size",
    // which in our case is 1.  -1 is the default in
    // StatusTestGenResNorm.
    const int defQuorum = -1;

    // "Show Maximum Residual Norm Only" is only meaningful when the
    // "Block Size" parameter is > 1.  Our solver only supports a
    // Block Size of 1.
    const bool showMaxResNormOnly = false;

    // The "implicit" residual test checks the "native" residual
    // norm to determine if convergence was achieved.  It is less
    // expensive than the "explicit" residual test.
    RCP<res_norm_test> implicitTest (new res_norm_test (convTol, defQuorum));
    implicitTest->defineScaleForm (stringToScaleType (implicitScaleType), 
				   Belos::TwoNorm);
    implicitTest->setShowMaxResNormOnly (showMaxResNormOnly);

    // If there's a left preconditioner, create a combined status
    // test that check first the "explicit," then the "implicit"
    // residual norm, requiring that both have converged to within
    // the specified tolerance.  Otherwise, we only perform the
    // "implicit" test.
    RCP<res_norm_test> explicitTest;
    if (haveLeftPreconditioner) // ! problem_->getLeftPrec().is_null()
      {
	explicitTest = rcp (new res_norm_test (convTol, defQuorum));
	explicitTest->defineResForm (res_norm_test::Explicit, Belos::TwoNorm);
	explicitTest->defineScaleForm (stringToScaleType (explicitScaleType),
				       Belos::TwoNorm);
	explicitTest->setShowMaxResNormOnly (showMaxResNormOnly);
      }
    // The "final" convergence test: 
    //
    // First, the implicit residual norm test,
    // Followed by the explicit residual norm test if applicable,
    // Followed by the user-defined convergence test if supplied.
    RCP<base_test> convTest;
    if (explicitTest.is_null())
      convTest = implicitTest;
    else
      // The "explicit" residual test is only performed once the
      // native ("implicit") residual is below the convergence
      // tolerance.
      convTest = rcp (new combo_test (combo_test::SEQ, 
				      implicitTest, 
				      explicitTest));
    if (! userConvTest.is_null())
      convTest = rcp (new combo_test (combo_test::SEQ, 
				      convTest, 
				      userConvTest));
    return convTest;
  }


  template<class Scalar, class MV, class OP>
  Teuchos::RCP<StatusTest<Scalar,MV,OP> >
  GmresSolMgr<Scalar,MV,OP>::
  initStatusTest (Teuchos::RCP<StatusTest<Scalar,MV,OP> > convTest,
		  const int maxIters)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef StatusTest<Scalar,MV,OP> base_test;
    typedef StatusTestMaxIters<Scalar,MV,OP> max_iter_test;
    typedef StatusTestCombo<Scalar,MV,OP> combo_test;

    // Stopping criterion for maximum number of iterations.
    //
    // We could just keep this around and call setMaxIters() on it
    // when the maximum number of iterations changes, but why
    // bother?  It's not a heavyweight object.  Just make a new one.
    RCP<max_iter_test> maxIterTest (new max_iter_test (maxIters));

    // The "final" stopping criterion:
    //
    // Either we've run out of iterations, OR we've converged.
    return rcp (new combo_test (combo_test::OR, maxIterTest, convTest));
  }

  template<class Scalar, class MV, class OP>
  Teuchos::RCP<const Teuchos::ParameterList>
  GmresSolMgr<Scalar,MV,OP>::getDefaultParameters()
  {
    using Teuchos::RCP;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;

    // FIXME (mfh 15 Feb 2011) This makes the routine nonreentrant.
    static RCP<const ParameterList> defaultParams;
    if (defaultParams.is_null())
      {
	RCP<ParameterList> params = parameterList();
	params->setName ("GmresSolMgr");
	//
	// The OutputManager only depends on the verbosity and the
	// output stream.
	//
	// Only display errors by default.
	//
	// FIXME (mfh 16 Feb 2011) Annoyingly, I have to write this as
	// an int and not a MsgType, otherwise it's saved as a string.
	MsgType defaultVerbosity = Belos::Errors;
	params->set ("Verbosity", static_cast<int>(defaultVerbosity),
		     "What type(s) of solver information should be written "
		     "to the output stream.");
	// Output goes to std::cout by default.
	RCP<std::ostream> defaultOutStream = Teuchos::rcpFromRef (std::cout);
	params->set ("Output Stream", defaultOutStream, 
		     "A reference-counted pointer to the output stream where "
		     "all solver output is sent.");
	//
	// The status test output class (created via
	// StatusTestOutputFactory) depends on the output style and
	// frequency.
	//
	// FIXME (mfh 16 Feb 2011) Annoyingly, even when I insist (via
	// a cast to OutputType) that the value has type OutputType,
	// it's stored in the XML file as a Teuchos::any, and read
	// back in as a string.  Thus, I cast to int instead.
	//
	// Default output style.
	OutputType defaultOutputStyle = Belos::General;
	params->set ("Output Style", (int) defaultOutputStyle,
		     "What style is used for the solver information written "
		     "to the output stream.");
	// Output frequency level refers to the number of iterations
	// between status outputs.  The default is -1, which means
	// "never."
	const int defaultOutFreq = -1;
	params->set ("Output Frequency", defaultOutFreq,
		     "How often (in terms of number of iterations) "
		     "convergence information should be written to "
		     "the output stream.  (-1 means \"never.\")");

	// Default orthogonalization type is "Simple," which does
	// block Gram-Schmidt for project() and MGS for normalize().
	const std::string defaultOrthoType ("Simple");
	params->set ("Orthogonalization", defaultOrthoType);

	// Parameters for the orthogonalization.  These depend on
	// the orthogonalization type, so if you change that, you
	// should also change its parameters.
	OrthoManagerFactory<Scalar, MV, OP> orthoFactory;
	RCP<const ParameterList> defaultOrthoParams = 
	  orthoFactory.getDefaultParameters (defaultOrthoType);
	// Storing the ParameterList as a sublist, rather than as an
	// RCP, ensures correct input and output.
	params->set ("Orthogonalization Parameters", *defaultOrthoParams);

	// Maximum number of iterations per restart cycle.
	const int defaultMaxIters = 100;
	params->set ("Maximum Iterations", defaultMaxIters,
		     "Maximum number of iterations allowed per restart cycle, "
		     "for each right-hand side solved.");
	// Maximum number of restart cycle(s) per right-hand side.
	const int defaultMaxRestarts = 10;
	params->set ("Maximum Restarts", defaultMaxRestarts,
		     "Maximum number of restart cycle(s) allowed for each "
		     "right-hand side solved.");

	// Default convergence tolerance is the square root of
	// machine precision.
	const magnitude_type defaultConvTol = STM::squareroot (STM::eps());
	params->set ("Convergence Tolerance", defaultConvTol,
		     "Relative residual tolerance that must be attained by "
		     "the iterative solver in order for the linear system to "
		     "be declared Converged.");
	//
	// The implicit (and explicit, if applicable) relative
	// residual tests needs an initial scaling, which makes the
	// norms "relative."  The default for the implicit test is
	// to use the (left-)preconditioned initial residual (which
	// is just the initial residual in the case of no left
	// preconditioning).  The default for the explicit test is
	// to use the non-preconditioned initial residual.
	//
	const std::string defaultImplicitScaleType = 
	  "Norm of Preconditioned Initial Residual";
	params->set ("Implicit Residual Scaling", 
		     defaultImplicitScaleType,			
		     "The type of scaling used in the implicit residual "
		     "convergence test.");
	const std::string defaultExplicitScaleType =
	  "Norm of Initial Residual";
	params->set ("Explicit Residual Scaling", 
		     defaultExplicitScaleType,
		     "The type of scaling used in the explicit residual "
		     "convergence test.");

	// TODO (mfh 21 Feb 2011) Fill in default values for the GMRES
	// iteration parameters.  Right now, the parameter list exists
	// but is empty.
	ParameterList iterParams;
	params->set ("Iteration Parameters", iterParams);

	// TODO (mfh 15 Feb 2011) Set any the other default parameters.

	defaultParams = params;
      }
    return defaultParams;
  }

  template<class Scalar, class MV, class OP>
  void
  GmresSolMgr<Scalar,MV,OP>::
  setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using Teuchos::ParameterList;
    using Teuchos::rcp_const_cast;
    // const_cast is OK, because setParametersImpl doesn't modify its
    // input.  We just have to pass in a non-const ParameterList
    // because Belos::SolutionManager requires it for implementations
    // of setParameters().
    setParametersImpl (rcp_const_cast<const ParameterList> (params));
  }

  template<class Scalar, class MV, class OP>
  void
  GmresSolMgr<Scalar,MV,OP>::
  setParametersImpl (const Teuchos::RCP<const Teuchos::ParameterList>& params)
  {
    using Teuchos::Exceptions::InvalidParameter;
    using Teuchos::Exceptions::InvalidParameterType;
    using Teuchos::null;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;

    RCP<const ParameterList> defaultParams = getDefaultParameters();
    TEST_FOR_EXCEPTION(defaultParams.is_null(), std::logic_error, 
		       "Belos::GmresSolMgr::setParametersImpl: The default "
		       "parameter list is null; this should never happen.  "
		       "Please report this bug to the Belos developers.");
    RCP<ParameterList> actualParams;
    if (params.is_null())
      actualParams = parameterList (*defaultParams);
    else
      { // Make a deep copy of the given parameter list.  This ensures
	// that the solver's behavior won't change, even if users
	// modify params later on.  Users _must_ invoke
	// setParameters() in order to change the solver's behavior.
	actualParams = parameterList (*params);

	// Fill in default values for parameters that aren't provided,
	// and make sure that all the provided parameters' values are
	// correct.
	//
	// FIXME (mfh 16 Feb 2011) Reading the output stream (which
	// has type RCP<std::ostream>) from a ParameterList may be
	// impossible if the ParameterList was read in from a file.
	// We hackishly test for this by catching
	// InvalidParameterType, setting the output stream in
	// actualParams to its default value, and redoing the
	// validation.  This is a hack because we don't know whether
	// the "Output Stream" parameter really caused that exception
	// to be thrown.
	bool success = false;
	try {
	  actualParams->validateParametersAndSetDefaults (*defaultParams);
	  success = true;
	} catch (InvalidParameterType&) {
	  success = false;
	}
	if (! success)
	  {
	    RCP<std::ostream> outStream = 
	      defaultParams->get<RCP<std::ostream> > ("Output Stream");
	    actualParams->set ("Output Stream", outStream, 
			       "A reference-counted pointer to the output "
			       "stream where all solver output is sent.");
	    // Retry the validation.
	    actualParams->validateParametersAndSetDefaults (*defaultParams);
	    success = true;
	  }
      }

    // Use the given name if one was provided, otherwise name the
    // parameter list appropriately.
    if (params.is_null() || params->name() == "" || params->name() == "ANONYMOUS")
      actualParams->setName (defaultParams->name());
    else
      actualParams->setName (params->name());

    // If we don't have a problem to solve yet, we haven't previously
    // gotten parameters, or the maximum number of iterations per
    // restart cycle has changed, we will need to reset the solver.
    //
    // Changing the maximum number of iterations requires resetting
    // the solver, since GMRES usually preallocates this number of
    // vectors.
    //
    // FIXME (mfh 21 Feb 2011) Perhaps we should ask the specific
    // GMRES implementation whether it needs to be reset?
    const bool needToResetSolver = problem_.is_null() || params_.is_null() || 
      params_->get<int> ("Maximum Iterations") != 
      actualParams->get<int> ("Maximum Iterations");

    // Initialize the OutputManager.    
    {
      // FIXME (mfh 16 Feb 2011) Annoyingly, I have to read this back
      // in as an int rather than a Belos::MsgType.
      //
      // TODO (mfh 15 Feb 2011) Validate verbosity; MsgType is really
      // just an int, so it might have an invalid value.  (C++
      // typedefs don't have a nice type theory.)
      const MsgType verbosity = (MsgType) actualParams->get<int> ("Verbosity");

      // Reading the output stream from a ParameterList can be tricky
      // if the ParameterList was read in from a file.  Serialization
      // of pointers to arbitrary data is very, very difficult, and
      // ParameterList doesn't try.  If the XML file shows the type as
      // "any", the get() call below may throw InvalidParameterType.
      // In that case, we set the outStream to point to std::cout, as
      // a reasonable default.
      RCP<std::ostream> outStream;
      try {
	outStream = actualParams->get<RCP<std::ostream> > ("Output Stream");
      } catch (InvalidParameterType&) {
	outStream = Teuchos::rcpFromRef (std::cout);
      }
	
      // Sanity check
      //
      // FIXME (mfh 16 Feb 2011) Should this be a "black hole" stream
      // on procs other than Proc 0?
      if (outStream.is_null())
	outStream = Teuchos::rcpFromRef (std::cout);

      // This will do the right thing whether or not outMan_ has
      // previously been initialized (== is not null).
      outMan_ = initOutputManager (outMan_, verbosity, outStream);
    }
    TEST_FOR_EXCEPTION(outMan_.is_null(), std::logic_error,
		       "Belos::GmresSolMgr::setParametersImpl: OutputManager "
		       "instance is null after its initialization; this should"
		       " never happen.  Please report this bug to the Belos "
		       "developers.");

    // Now that we've initialized the output manager, we can print
    // debug output.
    std::ostream& dbg = outMan_->stream(Debug);
    dbg << "Belos::GmresSolMgr::setParametersImpl:" << endl
	<< "-- Initialized output manager; now we can print debug output." 
	<< endl;

    // (Re)initialize the stopping criteria.
    //
    // FIXME (mfh 21 Feb 2011) Currently, userConvTest_ remains unset
    // until setUserConvStatusTest() is called.  Should we look for a
    // custom test in the parameter list?
    rebuildStatusTests (actualParams);
    dbg << "-- Initialized status tests." << endl;

    // Initialize the "output status test."  This must be done after
    // initializing statusTest_.
    {
      // FIXME (mfh 16 Feb 2011) Annoyingly, I have to read this back
      // in as an int rather than a Belos::OutputType.
      //
      // TODO (mfh 15 Feb 2011) Validate.
      const OutputType outStyle = 
	(OutputType) actualParams->get<int> ("Output Style");

      // TODO (mfh 15 Feb 2011) Validate.
      const int outFreq = actualParams->get<int> ("Output Frequency");

      // TODO (mfh 15 Feb 2011) Set the solver description according
      // to the specific kind of GMRES (CA-GMRES, standard GMRES,
      // Flexible GMRES, Flexible CA-GMRES) being used.
      const std::string solverDesc (" GMRES ");
      // TODO (mfh 15 Feb 2011) Get a preconditioner description.
      const std::string precondDesc;
      outTest_ = initOutputTest (outMan_, statusTest_, outStyle, outFreq, 
				 solverDesc, precondDesc);
    }
    dbg << "-- Initialized output status test." << endl;

    // (Re)initialize the orthogonalization.
    rebuildOrthoManager (actualParams);
    dbg << "-- Initialized OrthoManager subclass instance." << endl;

    // We don't initialize the Iteration subclass (GmresBaseIteration)
    // here; we do that on demand.  This is because GmresBaseIteration
    // instantiates GmresBase, and GmresBase wants setLSIndex() to
    // have been called on the LinearProblem.  setLSIndex() won't be
    // called until the solve() method.

    // Note that the parameter list we store contains the actual
    // parameter values used by this solver manager, not necessarily
    // those supplied by the user.  We reserve the right to fill in
    // default values and silently correct bad values.
    params_ = actualParams;

    // If necessary, restart the current solve that might be in
    // progress.
    if (needToResetSolver)
      reset (Belos::Problem);

    dbg << "-- Done with setParametersImpl()." << endl;
  }

  template<class Scalar, class MV, class OP>
  void
  GmresSolMgr<Scalar,MV,OP>::
  setUserConvStatusTest (const Teuchos::RCP<StatusTest<Scalar,MV,OP> >& userConvTest)
  {
    userConvTest_ = userConvTest;

    // For now, we just rebuild the entire stopping criterion from
    // scratch.  We could instead keep the parts of the stopping
    // criterion that don't change, and just change the user-defined
    // convergence criterion.
    rebuildStatusTests ();
  }

  template<class Scalar, class MV, class OP>
  void
  GmresSolMgr<Scalar,MV,OP>::
  rebuildStatusTests (Teuchos::RCP<const Teuchos::ParameterList> plist)
  {
    using Teuchos::rcp_const_cast;
    using Teuchos::ParameterList;
    using Teuchos::RCP;

    // Default value for plist is null, in which case we use the
    // stored parameter list.  One of those two should be non-null.
    RCP<const ParameterList> theParams = plist.is_null() ? 
      rcp_const_cast<const ParameterList>(params_) : plist;
    TEST_FOR_EXCEPTION(theParams.is_null(),
		       std::logic_error,
		       "Belos::GmresSolMgr::rebuildStatusTests: We can't (re)"
		       "build the status tests without any parameters.");
    const magnitude_type convTol = 
      theParams->get<magnitude_type> ("Convergence Tolerance");
    TEST_FOR_EXCEPTION(convTol < STM::zero(), std::invalid_argument,
		       "Convergence tolerance " << convTol << " is negative.");
    const int maxIters = 
      theParams->get<int> ("Maximum Iterations");
    TEST_FOR_EXCEPTION(maxIters < 0, std::invalid_argument,
		       "Maximum number of iterations " << maxIters 
		       << " is negative.");
    // TODO (mfh 15 Feb 2011) Validate.
    const std::string implicitScaleType = 
      theParams->get<std::string> ("Implicit Residual Scaling");
    // TODO (mfh 15 Feb 2011) Validate.
    const std::string explicitScaleType = 
      theParams->get<std::string> ("Explicit Residual Scaling");

    // If we don't have a problem to solve yet (problem_ is null), do
    // both an implicit and an explicit convergence test by default.
    // Then, when the user later calls setProblem() with a problem to
    // solve, that method will call reset(Belos::Problem), which in
    // turn will call rebuildStatusTests() again.
    const bool haveLeftPrecond = problem_.is_null() ||
      ! problem_->getLeftPrec().is_null();

    convTest_ = initConvTest (convTol, haveLeftPrecond,
			      implicitScaleType, explicitScaleType,
			      userConvTest_);
    statusTest_ = initStatusTest (convTest_, maxIters);
  }


  template<class Scalar, class MV, class OP>
  void
  GmresSolMgr<Scalar,MV,OP>::
  rebuildStatusTests (const typename GmresSolMgr<Scalar,MV,OP>::magnitude_type convTol,
		      const int maxItersPerRestart)
  {
    using Teuchos::ParameterList;
    using Teuchos::RCP;

    TEST_FOR_EXCEPTION(convTol < STM::zero(), std::invalid_argument,
		       "Convergence tolerance " << convTol << " is negative.");
    params_->set ("Convergence Tolerance", convTol);
    TEST_FOR_EXCEPTION(maxItersPerRestart < 0, std::invalid_argument,
		       "Maximum number of iterations " << maxItersPerRestart
		       << " per restart cycle is negative.");
    params_->set ("Maximum Iterations", maxItersPerRestart,
		  "Maximum number of iterations allowed per restart cycle, "
		  "for each right-hand side solved.");
    // If we don't have a problem to solve yet (problem_ is null), do
    // both an implicit and an explicit convergence test by default.
    // Then, when the user later calls setProblem() with a problem to
    // solve, that method will call reset(Belos::Problem), which in
    // turn will call rebuildStatusTests() again.
    const bool haveLeftPrecond = problem_.is_null() ||
      ! problem_->getLeftPrec().is_null();
    const std::string implicitScaleType = 
      params_->get<std::string> ("Implicit Residual Scaling");
    const std::string explicitScaleType = 
      params_->get<std::string> ("Explicit Residual Scaling");
    convTest_ = initConvTest (convTol, haveLeftPrecond,
			      implicitScaleType, explicitScaleType,
			      userConvTest_);
    statusTest_ = initStatusTest (convTest_, maxItersPerRestart);
  }


  template<class Scalar, class MV, class OP>
  void
  GmresSolMgr<Scalar,MV,OP>::
  rebuildOrthoManager (Teuchos::RCP<const Teuchos::ParameterList> plist)
  {
    using Teuchos::rcp_const_cast;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::null;

    const char prefix[] = "Belos::GmresSolMgr::rebuildOrthoManager: ";
    TEST_FOR_EXCEPTION(outMan_.is_null(), std::logic_error,
		       prefix << "Cannot (re)build the orthogonalization "
		       "method unless the OutputManager has been "
		       "instantiated.");

    // Default value for plist is null, in which case we use the
    // stored parameter list.  One of those two should be non-null.
    RCP<const ParameterList> actualParams = plist.is_null() ? 
      rcp_const_cast<const ParameterList>(params_) : plist;
    TEST_FOR_EXCEPTION(actualParams.is_null(), std::logic_error,
		       prefix << "We can't (re)build the orthogonalization "
		       "method without any parameters.");

    // Get the orthogonalization method name.
    const std::string orthoType = 
      actualParams->get<std::string> ("Orthogonalization");
    
    // Validate the orthogonalization method name.
    //
    // TODO (mfh 16 Feb 2011) Encode the validator in the default
    // parameter list.
    OrthoManagerFactory<Scalar, MV, OP> factory;
    if (! factory.isValidName (orthoType))
      {
	std::ostringstream os;
	os << prefix << "Invalid orthogonalization type \"" << orthoType 
	   << "\".  Valid orthogonalization types: "
	   << factory.printValidNames(os) << ".";
	throw std::invalid_argument (os.str());
      }

    // Get the parameters for that orthogonalization method.
    //
    // FIXME (mfh 16 Feb 2011) Extraction via reference is legitimate
    // only if we know that the whole parameter list won't go away.
    // Some OrthoManager subclasses might not copy their input
    // parameter lists deeply.
    const ParameterList& orthoParams = 
      actualParams->sublist ("Orthogonalization Parameters");
      
    // (Re)instantiate the orthogonalization manager.  Don't bother
    // caching this, since it's too much of a pain to check whether
    // any of the parameters have changed.

    // Set the timer label for orthogonalization
    std::ostringstream os; 
    os << "Orthogonalization (method \"" << orthoType << "\")";
    orthoMan_ = factory.makeOrthoManager (orthoType, null, outMan_, os.str(),
					  Teuchos::rcpFromRef (orthoParams));
  }

  template<class Scalar, class MV, class OP>
  void
  GmresSolMgr<Scalar,MV,OP>::
  rebuildIteration (Teuchos::RCP<const Teuchos::ParameterList> plist)
  {
    using Teuchos::rcp_const_cast;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    const char prefix[] = "Belos::GmresSolMgr::rebuildIteration: ";

    std::ostream& dbg = outMan_->stream(Debug);
    dbg << prefix << endl;

    RCP<const ParameterList> theParams = plist.is_null() ? 
      rcp_const_cast<const ParameterList>(params_) : plist;
    TEST_FOR_EXCEPTION(theParams.is_null(), std::logic_error,
		       prefix << "We can't (re)build the Iteration subclass "
		       "instance without any parameters.");
    TEST_FOR_EXCEPTION(! theParams->isSublist ("Iteration Parameters"),
		       std::logic_error,
		       prefix << "The parameter list needs an \"Iteration "
		       "Parameters\" sublist.");

    // Extraction of a sublist via reference is legitimate only if we
    // know that the parent parameter list won't go away.  That's why
    // we make a deep copy.
    const ParameterList& theSubList = theParams->sublist ("Iteration Parameters");
    RCP<const ParameterList> iterParams = parameterList (theSubList);
    dbg << "-- Instantiating GmresBaseIteration instance:" << endl;
    TEST_FOR_EXCEPTION(problem_.is_null(), std::logic_error, 
		       prefix << "LinearProblem instance is null.");
    TEST_FOR_EXCEPTION(orthoMan_.is_null(), std::logic_error, 
		       prefix << "OrthoManager instance is null.");
    TEST_FOR_EXCEPTION(outMan_.is_null(), std::logic_error, 
		       prefix << "OutputManager instance is null.");
    TEST_FOR_EXCEPTION(statusTest_.is_null(), std::logic_error, 
		       prefix << "StatusTest instance is null.");
    TEST_FOR_EXCEPTION(iterParams.is_null(), std::logic_error, 
		       prefix << "Iteration parameters list is null.");
    iter_ = rcp (new iteration_type (problem_, orthoMan_, outMan_, 
				     statusTest_, iterParams));
    dbg << "-- Done instantiating GmresBaseIteration instance." << endl;
  }

  template<class Scalar, class MV, class OP>
  void
  GmresSolMgr<Scalar,MV,OP>::reset (const ResetType type) 
  {
    using Teuchos::rcp;
    using std::endl;

    const char prefix[] = "Belos::GmresSolMgr::reset: ";
    std::ostream& dbg = outMan_->stream(Debug);
    dbg << prefix << endl;

    if ((type & Belos::RecycleSubspace) != 0)
      // TODO (mfh 15 Feb 2011) Do we want to support recycling GMRES here?
      throw std::invalid_argument ("This version of the GMRES solution "
				   "manager does not currently support "
				   "Krylov subspace recycling.");
    else if ((type & Belos::Problem) != 0)
      { // Recompute the initial residual(s) (both left-
	// preconditioned (if applicable) and unpreconditioned).
	//
	// FIXME (mfh 21 Feb 2011) This recomputes _all_ the initial
	// residuals, since it calls setProblem() on the
	// LinearProblem instance.  Broken!  We should only compute
	// initial residuals for the unsolved right-hand sides.
	dbg << "-- Recomputing initial residual(s)" << endl;
	problem_ = validatedProblem (problem_);

	// Rebuild the status tests, using the solver manager's
	// current list of parameters.
	dbg << "-- Rebuilding status tests" << endl;
	rebuildStatusTests ();

	// Rebuild the orthogonalization method.  This is a bit of a
	// hack; TsqrOrthoManager specifically needs to be
	// reinitialized if the data layout of the vectors have
	// changed, but the MVT and OPT interfaces don't let me
	// check for that.
	dbg << "-- Rebuilding OrthoManager" << endl;
	rebuildOrthoManager ();
      }
    else
      {
	std::ostringstream os;
	os << prefix << "Invalid ResetType argument type = " 
	   << type << ".  Currently, this method only accepts accept type = "
	  "Belos::Problem (== " << Belos::Problem << ").";
	throw std::invalid_argument (os.str());
      }
  }

  template<class Scalar, class MV, class OP>
  ReturnType
  GmresSolMgr<Scalar,MV,OP>::solve ()
  {
    using Teuchos::RCP;
    using std::endl;
    using std::make_pair;
    using std::pair;
    using std::vector;

    const char prefix[] = "Belos::GmresSolMgr::solve(): ";
    std::ostream& dbg = outMan_->stream(Debug);
    dbg << prefix << endl;

    // Reset the status test output.
    dbg << "-- Resetting status test output" << endl;
    outTest_->reset();
    // Reset the number of calls about which the status test output knows.
    dbg << "-- Resetting status test number of calls" << endl;
    outTest_->resetNumCalls();

    TEST_FOR_EXCEPTION(problem_.is_null(), std::logic_error,
		       prefix << "The LinearProblem instance is null.");
    TEST_FOR_EXCEPTION(problem_->getOperator().is_null(), std::logic_error,
		       prefix << "The LinearProblem's operator (the matrix A "
		       "in the equation AX=B) is null.");
    TEST_FOR_EXCEPTION(problem_->getRHS().is_null(), std::logic_error,
		       prefix << "The LinearProblem's initial guess X is "
		       "null.");
    TEST_FOR_EXCEPTION(problem_->getRHS().is_null(), std::logic_error,
		       prefix << "The LinearProblem's right-hand side B is "
		       "null.");

    // Set up the linear problem by computing initial residual(s)
    // for all right-hand side(s).  setProblem() takes optional
    // argument(s) if we want to change the initial guess(es) or the
    // right-hand side(s).  At some point we might want to change
    // the initial guess(es) for future solve(s), if we know that
    // the right-hand side(s) are related.
    if (! problem_->isProblemSet())
      {
	dbg << "-- (Re)computing initial residual(s) for all right-hand "
	  "side(s)" << endl;
	problem_->setProblem ();
      }
    // Total number of linear systems (right-hand sides) to solve.
    const int numRHS = MVT::GetNumberVecs (*(problem_->getRHS()));
    dbg << "-- There are " << numRHS << " total right-hand side" 
	<< (numRHS != 1 ? "s" : "") << " to solve." << endl;

    // Keep track of which linear system(s) converged.  Initially,
    // none of them have yet converged, since we haven't tested for
    // their convergence yet.
    //
    // FIXME (mfh 21 Feb 2011) Oh dear, std::vector<bool>.... full
    // of hackish badness with respect to iterators.
    vector<bool> converged (numRHS, false);
    
    // For each right-hand side, keep track of the total number of
    // restart cycle(s) (first element in the pair), and the total
    // number of iterations over all the restart cycle(s) (second
    // element in the pair).
    vector<pair<int, int> > totalNumIters (numRHS, make_pair (0, 0));

    // Maximum number of restart cycles.
    const int maxNumRestarts = params_->get<int> ("Maximum Restarts");
    dbg << "-- Max number of restart cycles: " << maxNumRestarts << endl;
      
    // For each right-hand side, restart until GMRES converges or we
    // run out of restart cycles.
    for (int curRHS = 0; curRHS < numRHS; ++curRHS)
      {
	dbg << "-- Solving for right-hand side " << (curRHS+1) 
	    << " of " << numRHS << ":" << endl;

	// Tell the linear problem for which right-hand side we are
	// currently solving.  A block solver could solve for more
	// than one right-hand side at a time; the GMRES solvers that
	// this solution manager knows how to manage can only solve
	// for one right-hand side at a time.
	vector<int> curProbIndex (1);
	curProbIndex[0] = curRHS;
	// This method ensures that problem_->getCurr{LHS,RHS}Vec()
	// return the right vector, corresponding to the current
	// linear system to solve.
	problem_->setLSIndex (curProbIndex);

	// Sanity checks.
	{
	  // Make sure that setLSIndex() has been called on the linear
	  // problem.  Otherwise, getCurr{L,R}HSVec() (corresponding to
	  // the current linear system to solve) will return null.
	  RCP<const MV> X = problem_->getCurrLHSVec();
	  RCP<const MV> B = problem_->getCurrRHSVec();
	  TEST_FOR_EXCEPTION(X.is_null() || B.is_null(), std::invalid_argument,
			     prefix << "setLSIndex() must not yet have been "
			     "called on the given Belos::LinearProblem "
			     "instance, since getCurrLHSVec() and/or "
			     "getCurrRHSVec() return null.");
	  // Our GMRES implementations only know how to solve for one
	  // right-hand side at a time.  Make sure that X and B each have
	  // exactly one column.
	  const int nLHS = MVT::GetNumberVecs(*X);
	  const int nRHS = MVT::GetNumberVecs(*B);
	  TEST_FOR_EXCEPTION(nLHS != 1 || nRHS != 1, std::invalid_argument,
			     prefix << "The current linear system to solve has "
			     << nLHS << " initial guess" 
			     << (nLHS != 1 ? "es" : "")
			     << " and " << nRHS << " right-hand side"
			     << (nLHS != 1 ? "s" : "") 
			     << ", but our GMRES solvers only know how to "
			     "solve for one right-hand side at a time.");
	  if (debug_)
	    { // In debug mode, compute and print the initial
	      // residual independently of the solver framework, as
	      // a sanity check.
	      RCP<const OP> A = problem_->getOperator ();
	      RCP<MV> R = MVT::Clone (*B, MVT::GetNumberVecs (*B));
	      // R := A * X
	      OPT::Apply (*A, *X, *R);
	      // R := B - R
	      MVT::MvAddMv (STS::one(), *B, -STS::one(), *R, *R);
	      std::vector<magnitude_type> theNorm (MVT::GetNumberVecs (*R));
	      MVT::MvNorm (*R, theNorm);
	      dbg << "-- Initial residual norm: ||B - A*X||_2 = " 
		  << theNorm[0] << endl;
	    }
	}

	// (Re)initialize the GmresBaseIteration, and therefore the
	// GmresBase subclass instance, with the current linear
	// system to solve.  (We need to do this since we've moved
	// on to the next right-hand side, which means we're solving
	// a different linear system, even though the
	// "LinearProblem" object expresses the same set of linear
	// system(s) to solve.
	rebuildIteration ();

	// Restart until GMRES converges or we run out of restart
	// cycles.
	for (int restartCycle = 0; 
	     convTest_->getStatus() != Passed && restartCycle < maxNumRestarts;
	     ++restartCycle)
	  {
	    dbg << "-- Restart cycle " << (restartCycle+1) << " of " 
		<< maxNumRestarts << ":" << endl;
	    // reset() restarts the iteration, so we don't need to
	    // restart the first time.
	    if (restartCycle > 0)
	      iter_->restart ();
	    // Iterate to convergence or maximum number of
	    // iterations.
	    iter_->iterate ();
	    // For the current right-hand side, keep track of the
	    // number of restart cycles and the total number of
	    // iterations over all restart cycles.
	    totalNumIters[curRHS].first++;
	    totalNumIters[curRHS].second += iter_->getNumIters();

	    // Update the current approximate solution in the linear problem.
	    iter_->updateSolution ();
	  }
	// Remember whether the current linear system converged.
	converged[curRHS] = (convTest_->getStatus() == Passed);

	// Tell the linear problem that we have finished attempting
	// to solve the current linear system, and are ready to move
	// on to the next system (if there is one).
	//
	// FIXME (mfh 21 Feb 2011) Yeah, I know, this LinearProblem
	// method name doesn't make any sense.  Do _you_ have time to
	// rewrite most of Belos?
	problem_->setCurrLS();

	// Remember the total number of restarts and iterations.
	totalNumIters_ = totalNumIters;
      }

    // We've Converged if all of the right-hand sides converged.
    if (std::accumulate (converged.begin(), converged.end(), 
			 true, std::logical_and<bool>()))
      return Converged;
    else
      return Unconverged;
  }


} // namespace Belos

#endif // __Belos_GmresSolMgr_hpp
