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

#include <BelosOutputManager.hpp>
#include <BelosSolverManager.hpp>
#include <BelosStatusTestGenResNorm.hpp>
#include <BelosStatusTestOutputFactory.hpp>
#include <BelosOrthoManagerFactory.hpp>
//#include <Teuchos_ScalarTraits.hpp>

namespace Belos {

  /// \class GmresSolMgr
  /// \author MarkHoemmen
  /// \brief Solver manager for CA-GMRES and standard GMRES.
  ///
  template<class Scalar, class MV, class OP>
  class GmresSolMgr : public SolverManager<Scalar, MV, OP> {
  private:
    typedef MultiVecTraits<Scalar,MV> MVT;
    typedef OperatorTraits<Scalar,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    typedef typename SCT::magnitudeType magnitude_type;
    typedef Teuchos::ScalarTraits<magnitude_type> SMT;

  public:
    typedef Scalar scalar_type;
    typedef MV multivector_type;
    typedef OP operator_type;

    //! \name Constructors and destructor 
    //@{ 

    /// \brief Default constructor.
    ///
    /// \note I'm not a fan of default construction for heavyweight
    /// "solve this problem" classes like this one, since it tends to
    /// violate the "construction is initialization" preferred C++
    /// idiom.  Nevertheless, we have to support default construction
    /// in order to inherit from SolverManager.
    GmresSolMgr() {}

    /// \brief Preferred constructor.
    ///
    /// \param problem [in/out] The linear problem to solve.
    /// \param params [in] Parameters for the solve.  This is 
    ///   a pointer to nonconst only because the SolutionManager
    ///   interface demands it (more or less...).
    GmresSolMgr (const Teuchos::RCP<LinearProblem<Scalar,MV,OP> >& problem,
		 const Teuchos::RCP<Teuchos::ParameterList>& params) : 
      problem_ (validatedProblem (problem))
    { 
      setParameters (params);
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

    //! Iteration count for the most recent call to \c solve().
    int getNumIters() const { 
      // TODO (mfh 15 Feb 2011) Return the actual value.
      return 0; 
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
      // Restart the current solve that might be in progress.
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
    /// \param params [in] New parameters for the linear solve
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
    void 
    reset (const ResetType type) 
    {
      if ((type & Belos::RecycleSubspace) != 0)
	// TODO (mfh 15 Feb 2011) Do we want to support recycling GMRES here?
	throw std::invalid_argument ("This version of the GMRES solution "
				     "manager does not currently support "
				     "Krylov subspace recycling.");
      else if ((type & Belos::Problem) != 0)
	{ // Invalidate the current approximate solution, recheck
	  // validity of the linear poblem to solve, and recompute the
	  // initial residual(s) (both left-preconditioned (if
	  // applicable) and unpreconditioned).
	  if (! problem_.is_null())
	    (void) problem_->setProblem();
	    
	  // TODO (mfh 15 Feb 2011) Restart GMRES, etc.
	}
      else
	{
	  std::ostringstream os;
	  os << "GmresSolMgr::reset(): Invalid ResetType argument type = " 
	     << type << ".  Currently, this method only accepts accept type = "
	    "Belos::Problem (== " << Belos::Problem << ").";
	  throw std::invalid_argument (os.str());
	}
    }
    //@}

    //! \name Solver application methods
    //@{ 


    /// \brief Attempt to solve the linear system.
    ///
    /// Do so by calling the underlying linear solver's iterate()
    /// routine zero or more times, until the problem has been solved
    /// (as decided by the solver manager) or the solver manager
    /// decides to quit.
    ///
    /// \return ReturnType enum specifying Converged (the linear
    ///   problem was solved to the desired tolerance) or Unconverged
    ///   (the linear problem was not thusly solved).
    ReturnType solve() {
      // Reset the status test output.
      outTest_->reset();
      // Reset the number of calls about which the status test output knows.
      outTest_->resetNumCalls();


      // TODO (mfh 15 Feb 2011) This method should of course invoke
      // the Iteration and try to solve the problem.
      return Unconverged;
    }

    //@}
    
  private:
    //! The linear problem to solve.
    Teuchos::RCP<LinearProblem<Scalar, MV, OP> > problem_;

    //! Output manager.
    Teuchos::RCP<OutputManager<Scalar> > outMan_;

    //! Current parameters for this solver manager instance.
    Teuchos::RCP<const Teuchos::ParameterList> params_;

    //! Status test (stopping criterion) for the current solve.
    Teuchos::RCP<StatusTest<Scalar,MV,OP> > statusTest_;

    /// \brief User-defined status test (stopping criterion).
    ///
    /// This will be combined in AND fashion with the current status
    /// test (not including the current user-defined status test, if
    /// there is any).
    Teuchos::RCP<StatusTest<Scalar,MV,OP> > userStatusTest_;

    //! "Status test" that outputs intermediate iteration results.
    Teuchos::RCP<StatusTestOutput<Scalar,MV,OP> > outTest_;

    //! The orthogonalization method to use for GMRES.
    Teuchos::RCP<Belos::OrthoManager<Scalar, MV> > orthoMan_;

    /// \brief Make sure that GmresSolMgr can solve the given problem.
    ///
    /// If it can, return the problem untouched, else throw an
    /// std::invalid_argument.  Currently, GmresSolMgr only knows how
    /// to solve problems with one right-hand side; this is what we
    /// check.
    static Teuchos::RCP<LinearProblem<Scalar, MV, OP> >
    validatedProblem (const Teuchos::RCP<LinearProblem<Scalar, MV, OP> >& problem)
    {
      const int numRHS = MVT::GetNumberVecs (*(problem->getRHS()));
      TEST_FOR_EXCEPTION(! problem.is_null() && numRHS != 1,
			 std::invalid_argument,
			 "Currently, GmresSolMgr only knows how to solve linear "
			 "problems with one right-hand side, but the provided "
			 "linear problem has " << numRHS << "right-hand side"
			 << (numRHS != 1 ? "s" : "") << ".");
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

    /// Create a new stopping criterion ("StatusTest") object.
    ///
    /// The stopping criterion is a prerequisite for creating the
    /// "output status test" object.
    ///
    /// \param convTest [in] Previously initialized convergence test
    /// \param maxIters [in] Maximum number of iterations to execute
    /// \param haveLeftPreconditioner [in] Whether or not a left
    ///   preconditioner has been provided.
    /// \param implicitScaleType [in]
    /// \param explicitScaleType [in]
    /// \param userConvTest [in] If not null, a user-defined
    ///   convergence test that will be performed sequentially after
    ///   the maximum iteration count test and the implicit (and
    ///   explicit) residual norm tests.
    ///
    /// \return Stopping criterion which is the Boolean OR of the
    ///   given convergence test, and the given maximum iteration
    ///   count.
    ///
    /// \note We don't take a previously initialized status test here,
    ///   because it's easier just to create a new object, rather than
    ///   keep around the status test as well as all of its child
    ///   tests.
    static Teuchos::RCP<StatusTest<Scalar,MV,OP> >
    initStatusTest (const magnitude_type convTol,
		    const int maxIters,
		    const bool haveLeftPreconditioner,
		    const std::string& implicitScaleType,
		    const std::string& explicitScaleType,
		    Teuchos::RCP<StatusTest<Scalar,MV,OP> > userConvTest = Teuchos::null);

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
  initStatusTest (const magnitude_type convTol,
		  const int maxIters,
		  const bool haveLeftPreconditioner,
		  const std::string& implicitScaleType,
		  const std::string& explicitScaleType,
		  Teuchos::RCP<StatusTest<Scalar,MV,OP> > userConvTest)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef StatusTest<Scalar,MV,OP> base_test;
    typedef StatusTestMaxIters<Scalar,MV,OP> max_iter_test;
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

    // Stopping criterion for maximum number of iterations.
    //
    // We could just keep this around and call setMaxIters() on it
    // when the maximum number of iterations changes, but why
    // bother?  It's not a heavyweight object.  Just make a new one.
    RCP<max_iter_test> maxIterTest (new max_iter_test (maxIters));

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
    // Followed by any user-defined convergence test.
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
	MsgType defaultVerbosity = Belos::Errors;
	params->set ("Verbosity", static_cast<MsgType>(defaultVerbosity),
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
	// Default output style.
	OutputType defaultOutputStyle = Belos::General;
	params->set ("Output Style", (OutputType) defaultOutputStyle,
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
	params->set ("Orthogonalization Parameters", defaultOrthoParams);

	// Maximum number of iterations.
	const int defaultMaxIters = 1000;
	params->set ("Maximum Iterations", defaultMaxIters,
		     "Maximum number of iterations allowed for each right-"
		     "hand side solved.");

	// Default convergence tolerance is the square root of
	// machine precision.
	const magnitude_type defaultConvTol = SMT::squareroot (SMT::eps());
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
    using Teuchos::Exceptions::InvalidParameter;
    using Teuchos::null;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<const ParameterList> defaultParams = getDefaultParameters();
    RCP<ParameterList> actualParams;
    if (params.is_null())
      actualParams = rcp(*defaultParams);
    else
      { // Make a deep copy of the given parameter list.  This ensures
	// that the solver's behavior won't change, even if users
	// modify params later on.  Users _must_ invoke
	// setParameters() in order to change the solver's behavior.
	actualParams = rcp(*params);
	// Fill in default values for parameters that aren't provided,
	// and make sure that all the provided parameters' values are
	// correct.
	actualParams->validateParametersAndSetDefaults (*defaultParams);
      }

    // Use the given name if one was provided, otherwise name the
    // parameter list appropriately.
    if (params.is_null() || params->name() == "" || params->name() == "ANONYMOUS")
      actualParams->setName (defaultParams->name());
    else
      actualParams->setName (params->name());

    // Changing certain parameters may require resetting the solver.
    // Of course we have to (re)set the solver if there's no linear
    // problem to solve yet, or if the parameters have not yet been
    // set.
    bool needToResetSolver = problem_.is_null() || params_.is_null();

    // Initialize the stopping criteria.
    {
      magnitude_type convTol;
      try {
	convTol = params->get<magnitude_type> ("Convergence Tolerance");
      } catch (InvalidParameter&) {
	convTol = defaultParams->get<magnitude_type> ("Convergence Tolerance");
      }
      // Validate first, before storing the value.
      TEST_FOR_EXCEPTION(convTol < SMT::zero(), 
			 std::invalid_argument,
			 "Convergence tolerance " << convTol 
			 << " is negative.");
      actualParams->set ("Convergence Tolerance", convTol);

      int maxIters;
      try {
	maxIters = params->get<int> ("Maximum Iterations");
      } catch (InvalidParameter&) {
	maxIters = defaultParams->get<int> ("Maximum Iterations");
      }
      TEST_FOR_EXCEPTION(maxIters < 0,
			 std::invalid_argument,
			 "Maximum number of iterations " << maxIters 
			 << " is negative.");
      actualParams->set ("Maximum Iterations", maxIters);
      // Changing the maximum number of iterations requires
      // resetting the solver, since GMRES usually preallocates this
      // number of vectors.  (Yay for short-circuiting OR!)
      needToResetSolver = params_.is_null() ||
	params_->get<int> ("Maximum Iterations") != maxIters;

      std::string implicitScaleType;
      try {
	implicitScaleType = params->get<std::string> ("Implicit Residual Scaling");
      } catch (InvalidParameter&) {
	implicitScaleType = defaultParams->get<std::string> ("Implicit Residual Scaling");
      }
      // TODO (mfh 15 Feb 2011) Validate first.
      actualParams->set ("Implicit Residual Scaling", implicitScaleType);

      std::string explicitScaleType;
      try {
	explicitScaleType = params->get<std::string> ("Explicit Residual Scaling");
      } catch (InvalidParameter&) {
	explicitScaleType = defaultParams->get<std::string> ("Explicit Residual Scaling");
      }
      // TODO (mfh 15 Feb 2011) Validate first.
      actualParams->set ("Explicit Residual Scaling", explicitScaleType);

      // If we don't have a problem to solve yet, do both an
      // implicit and an explicit convergence test by default.  When
      // we later get a problem to solve, just rebuild the stopping
      // criterion.
      const bool haveLeftPrecond = problem_.is_null() || 
	! problem_->getLeftPrec().is_null();
      statusTest_ = initStatusTest (convTol, maxIters, haveLeftPrecond, 
				    implicitScaleType, explicitScaleType);
    }

    // Initialize the OutputManager.
    {
      MsgType verbosity;
      try {
	verbosity = params->get<MsgType> ("Verbosity");
	// TODO (mfh 15 Feb 2011) Validate verbosity; MsgType is
	// really just an int, so it might have an invalid value.
	// (C++ typedefs don't have a nice type theory.)
      } catch (InvalidParameter&) {
	verbosity = defaultParams->get<MsgType> ("Verbosity");
      }
      // TODO (mfh 15 Feb 2011) Validate first.
      actualParams->set ("Verbosity", (MsgType) verbosity);

      RCP<std::ostream> outStream;
      try {
	outStream = params->get<RCP<std::ostream> > ("Output Stream");
      } catch (InvalidParameter&) {
	outStream = defaultParams->get<RCP<std::ostream> > ("Output Stream");
      }
      actualParams->set ("OutputStream", outStream);

      // This will do the right thing whether or not outMan_ has
      // previously been initialized (== is not null).
      outMan_ = initOutputManager (outMan_, verbosity, outStream);
    }

    // Initialize the "output status test."
    {
      OutputType outStyle;
      try {
	outStyle = params->get<OutputType> ("Output Style");
      } catch (InvalidParameter&) {
	outStyle = defaultParams->get<OutputType> ("Output Style");
      }
      // TODO (mfh 15 Feb 2011) Validate first.
      actualParams->set ("Output Style", (OutputType) outStyle);

      int outFreq;
      try {
	outFreq = params->get<int> ("Output Frequency");
      } catch (InvalidParameter&) {
	outFreq = defaultParams->get<int> ("Output Frequency");
      }
      // TODO (mfh 15 Feb 2011) Validate first.
      actualParams->set ("Output Frequency", outFreq);

      // TODO (mfh 15 Feb 2011) Set the solver description according
      // to the specific kind of GMRES (CA-GMRES, standard GMRES,
      // Flexible GMRES, Flexible CA-GMRES) being used.
      const std::string solverDesc (" GMRES ");
      // TODO (mfh 15 Feb 2011) Get a preconditioner description.
      const std::string precondDesc;
      outTest_ = initOutputTest (outMan_, statusTest_, outStyle, outFreq, 
				 solverDesc, precondDesc);
    }

    // Initialize the orthogonalization.
    {
      OrthoManagerFactory<Scalar, MV, OP> factory;
      
      // Get the orthogonalization method name.
      std::string orthoType;
      try {
	orthoType = params->get<std::string> ("Orthogonalization");
      } catch (InvalidParameter&) {
	orthoType = defaultParams->get<std::string> ("Orthogonalization");
      }
      // Validate the orthogonalization type.
      if (! factory.isValidName (orthoType))
	{
	  std::ostringstream os;
	  os << "Invalid orthogonalization type \"" << orthoType 
	     << "\".  Valid orthogonalization types: ";
	  factory.printValidNames (os);
	  os << ".";
	  throw std::invalid_argument (os.str());
	}

      // Get the parameters for that orthogonalization method.
      RCP<const ParameterList> orthoParams;
      try {
	orthoParams = params->get<RCP<const ParameterList> > ("Orthogonalization Parameters");
      } catch (InvalidParameter&) {
	orthoParams = defaultParams->get<RCP<const ParameterList> > ("Orthogonalization Parameters");
      }
      
      // (Re)instantiate the orthogonalization manager.
      //
      // TODO (mfh 15 Feb 2011) Set the orthogonalization timer label
      // appropriately.
      {
	// Set the timer label for orthogonalization
	std::ostringstream os; 
	os << "Orthogonalization (method \"" << orthoType << "\")";
	orthoMan_ = factory.makeOrthoManager (orthoType, null, os.str(),
					      orthoParams);
      }
    }
    // TODO (mfh 15 Feb 2011) Validate the other new parameters.

    // Note that the parameter list we store contains the actual
    // parameter values used by this solver manager, not necessarily
    // those supplied by the user.  We reserve the right to fill in
    // default values and silently correct bad values.
    params_ = actualParams;

    // If necessary, restart the current solve that might be in
    // progress.
    if (needToResetSolver)
      reset (Belos::Problem);
  }

  template<class Scalar, class MV, class OP>
  void
  GmresSolMgr<Scalar,MV,OP>::
  setUserConvStatusTest (const Teuchos::RCP<StatusTest<Scalar,MV,OP> >& userConvTest)
  {
    // For now, we just rebuild the entire stopping criterion from
    // scratch.  We could instead keep the parts of the stopping
    // criterion that don't change, and just change the user-defined
    // convergence criterion.

    const magnitude_type convTol = 
      params_->get<magnitude_type> ("Convergence Tolerance");
    const int maxIters = params_->get<int> ("Maximum Iterations");
    const std::string implicitScaleType = 
      params_->get<std::string> ("Implicit Residual Scaling");
    const std::string explicitScaleType = 
      params_->get<std::string> ("Explicit Residual Scaling");

    // If we don't have a problem to solve yet, do both an
    // implicit and an explicit convergence test by default.  When
    // we later get a problem to solve, just rebuild the stopping
    // criterion.
    const bool haveLeftPrecond = problem_.is_null() || 
      ! problem_->getLeftPrec().is_null();
    statusTest_ = initStatusTest (convTol, maxIters, haveLeftPrecond, 
				  implicitScaleType, explicitScaleType,
				  userConvTest);
  }

} // namespace Belos

#endif // __Belos_GmresSolMgr_hpp
