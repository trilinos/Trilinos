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

#include <BelosStatusTestCombo.hpp>
#include <BelosStatusTestGenResNorm.hpp>
#include <BelosStatusTestImpResNorm.hpp>
#include <BelosStatusTestMaxIters.hpp>
#include <Teuchos_ParameterList.hpp>

namespace Belos {

  /// \class StatusTestFactory
  /// \author Mark Hoemmen
  /// \brief A factory for making common cases of stopping criteria.
  ///
  /// Belos::StatusTest is an abstract interface for different
  /// stopping criteria ("status tests") for the iterative methods
  /// implemented in Belos.  There are many different kinds of status
  /// tests implemented in Belos, and they can be combined in
  /// different ways to make aggregate status tests (e.g., "reached
  /// max number of iterations, or residual norm converged").  We
  /// include a few common cases of those aggregate status tests in
  /// this factory class so that people implementing subclasses of
  /// Belos::SolverManager can minimize code duplication.  
  ///
  /// This factory also provides ways to change the convergence
  /// tolerance and/or maximum number of iterations, in place, for
  /// existing and possibly aggregate status tests.  This is an
  /// experimental optimization that attempts to avoid the overhead of
  /// constructing new status tests, but it comes at the cost of
  /// requiring \f$O(1)\f$ dynamic casts for each test in an aggregate
  /// collection (StatusTestCombo) of tests.  It may be useful if you
  /// are not allowed to create new status test objects, for example
  /// if you only have access to them through a copied RCP and not
  /// through a reference to the RCP.
  /// 
  /// The general idea of factory classes like this one is to avoid
  /// duplicated code whenever possible.  Scientific programmers' time
  /// is limited, so it's often better to avoid duplicated code, even
  /// if it means adding a new class or making the nonduplicated code
  /// a bit more complicated.
  template<class Scalar, class MV, class OP>
  class StatusTestFactory {
  public:
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    typedef StatusTest<Scalar,MV,OP> base_test;
    typedef StatusTestGenResNorm<Scalar,MV,OP> res_norm_test;
    typedef StatusTestMaxIters<Scalar,MV,OP> max_iter_test;
    typedef StatusTestCombo<Scalar,MV,OP> combo_test;

    /// \brief Status test suitable for (preconditioned) GMRES.
    ///
    /// The returned status test stops iteration if the maximum number
    /// of iterations has been reached, OR if the residual norm has
    /// converged.  The latter first iterates until the implicit
    /// residual norm (the one computed as a side effect of solving
    /// the projected least-squares problem in GMRES) reaches the
    /// given tolerance.  Then, it switches to computing the explicit
    /// residual norm (\f$\|b - Ax\|_2\f$) and iterates until that has
    /// reached the tolerance.
    ///
    /// \param convTol [in] Convergence tolerance.  The meaning of
    ///   this depends on the scaling used.
    ///
    /// \param maxIterCount [in] Maximum number of iterations of GMRES
    ///   to execute before stopping.
    ///
    /// \param haveLeftPreconditioner [in] Whether we will be running
    ///   GMRES with a left preconditioner.  This affects the residual
    ///   norm computation.
    ///
    /// \param implicitScaleType [in] Type of scaling to use for the
    ///   implicit residual norm computation.  Default is to scale by
    ///   the norm of the preconditioned (if applicable) initial
    ///   residual vector.  The \c stringToScaleType() method might be
    ///   useful.
    ///
    /// \param explicitScaleType [in] Type of scaling to use for the
    ///   explicit residual norm computation.  Default is to scale by
    ///   the norm of the (nonpreconditioned) initial residual vector.
    ///   The \c stringToScaleType() method might be useful.
    ///
    /// \param blockSize [in] Number of linear system(s) that the
    ///   solver can solve at one time.  Block iterative methods
    ///   support blockSize > 1.  Non-block iterative methods require
    ///   blockSize = 1.  This parameter is only used to verify
    ///   defQuorum, if the latter is not -1.
    ///
    /// \param defQuorum [in] "Deflation Quorum": number of converged
    ///   systems before deflation is allowed.  Cannot be larger than
    ///   "Block Size".  -1 is the default in
    ///   Belos::StatusTestGenResNorm, and it means that all of the
    ///   systems must converge before deflation is allowed.
    ///
    /// \param showMaxResNormOnly [in] "Show Maximum Residual Norm
    ///   Only" is only meaningful when the "Block Size" parameter is
    ///   > 1.
    ///
    /// \return Aggregate status test for GMRES.
    static Teuchos::RCP<base_test>
    gmresTest (const magnitude_type convTol,
	       const int maxIterCount,
	       const bool haveLeftPreconditioner,
	       const ScaleType implicitScaleType = Belos::NormOfPrecInitRes,
	       const ScaleType explicitScaleType = Belos::NormOfInitRes,
	       const int blockSize = 1,
	       const int defQuorum = -1,
	       const bool showMaxResNormOnly = false);

    /// \brief Overloaded \c gmresTest() for ParameterList input.
    ///
    /// Does the same thing as the other \c gmresTest(), except it
    /// reads the values from the given parameter list.  The parameter
    /// list must be nonnull and valid.  We make only a modest effort
    /// to fill in default values for nonrequired parameters only.  We
    /// throw std::invalid_argument if any required parameter is
    /// missing or any provided value is invalid.
    static Teuchos::RCP<base_test>
    gmresTest (const bool haveLeftPreconditioner, 
	       const Teuchos::RCP<const Teuchos::ParameterList>& params);

    /// \brief Change convergence tolerance and max number of iterations.
    ///
    /// The changes take place in place in the given status test.  If
    /// the status test is a StatusTestCombo, we recurse on all
    /// children, changing them if they are of the appropriate type.
    /// If we don't find anything to change, we do nothing, except
    /// report back in the return value.
    ///
    /// \return First element of the pair is whether we found at least
    ///   one residual norm test for which to change the convergence
    ///   tolerance.  Second element is whether we found at least one
    ///   maximum number of iterations test for which to change the
    ///   maximum iteration count.
    ///
    /// \note This method does nothing to protect against infinite
    ///   recursion if the StatusTestCombo tests were constructed
    ///   incorrectly.  However, StatusTestCombo itself prevents you
    ///   from forming loops in the graph of status tests, so this
    ///   method should always terminate.
    ///
    /// \note The given status test does _not_ need to be one that
    ///   this factory created.  This method should work with any
    ///   Belos::StatusTest instance.
    ///
    static std::pair<bool, bool>
    changeConvTolAndMaxIters (const Teuchos::RCP<base_test>& test, 
			      const magnitude_type convTol,
			      const int maxIterCount);


    /// \brief Change max number of iterations in place.
    /// 
    /// See the notes for \c changeConvTolAndMaxIters().
    ///
    /// \return Whether we found at least one test to change.
    static bool
    changeMaxNumIters (const Teuchos::RCP<base_test>& test, 
		       const int maxIterCount);

    /// \brief Change convergence tolerance in place.
    /// 
    /// See the notes for \c changeConvTolAndMaxIters().
    ///
    /// \return Whether we found at least one test to change.
    static bool
    changeConvTol (const Teuchos::RCP<base_test>& test, 
		   const magnitude_type convTol);

    /// Convert string to enum that tells residual test how to scale.
    ///
    /// \param scaleType [in] A string describing the type of
    /// scaling to be used in a residual norm convergence test.
    /// Valid values include:
    /// - "Norm of Initial Residual"
    /// - "Norm of Preconditioned Initial Residual"
    /// - "Norm of RHS" / "Norm of Right-Hand Side"
    /// - "None" (no scaling)
    ///
    /// \return The ScaleType enum value that tells the residual test
    ///   how to perform scaling.
    static ScaleType
    stringToScaleType (const std::string& scaleType);
  };


  template<class Scalar, class MV, class OP>
  Teuchos::RCP<typename StatusTestFactory<Scalar, MV, OP>::base_test>
  StatusTestFactory<Scalar, MV, OP>::
  gmresTest (const typename StatusTestFactory<Scalar, MV, OP>::magnitude_type convTol,
	     const int maxIterCount,
	     const bool haveLeftPreconditioner,
	     const ScaleType implicitScaleType = NormOfPrecInitRes,
	     const ScaleType explicitScaleType = NormOfInitRes,
	     const int blockSize = 1,
	     const int defQuorum = -1,
	     const bool showMaxResNormOnly = false)
  {
    using Teuchos::null;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    TEST_FOR_EXCEPTION(blockSize < 1, std::invalid_argument,
		       "blockSize (= " << blockSize << ") must be >= 1.");
    TEST_FOR_EXCEPTION(defQuorum > blockSize, std::invalid_argument,
		       "defQuorum (= " << defQuorum << ") may be no larger "
		       "than blockSize (= " << blockSize << ").");

    // The "implicit" residual test checks the "native" residual
    // norm to determine if convergence was achieved.  It is less
    // expensive than the "explicit" residual test.
    RCP<res_norm_test> implicitTest 
      = rcp (new res_norm_test (convTol, defQuorum));
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
    else
      explicitTest = null;

    // Full convergence test:
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

    // Stopping criterion for maximum number of iterations.
    RCP<max_iter_test> maxIterTest = rcp (new max_iter_test (maxIterCount));

    // The "final" stopping criterion:
    //
    // Either we've run out of iterations, OR we've converged.
    return rcp (new combo_test (combo_test::OR, maxIterTest, convTest));
  }


  template<class Scalar, class MV, class OP>
  Teuchos::RCP<typename StatusTestFactory<Scalar, MV, OP>::base_test>
  StatusTestFactory<Scalar, MV, OP>::
  gmresTest (const bool haveLeftPreconditioner, 
	     const Teuchos::RCP<const Teuchos::ParameterList>& params)
  {
    using Teuchos::Exceptions::InvalidParameter;
    using std::string;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    const magnitude_type convTol = 
      params->get<magnitude_type> ("Convergence Tolerance");
    TEST_FOR_EXCEPTION(convTol < STM::zero(), std::invalid_argument,
		       "Convergence tolerance " << convTol 
		       << " is negative.");
    const int maxIterCount = params->get<int> ("Maximum Iterations");
    TEST_FOR_EXCEPTION(maxIterCount < 0, std::invalid_argument,
		       "Maximum number of iterations " << maxIterCount
		       << " is negative.");

    ScaleType implicitScaleType = NormOfPrecInitRes;
    try {
      implicitScaleType = 
	stringToScaleType (params->get<string> ("Implicit Residual Scaling"));
    } catch (InvalidParameter&) {
      // Do nothing; leave default value
    }
    ScaleType explicitScaleType = NormOfInitRes;
    try {
      explicitScaleType = 
	stringToScaleType (params->get<string> ("Explicit Residual Scaling"));
    } catch (InvalidParameter&) {
      // Do nothing; leave default value
    }
    int blockSize = 1;
    try {
      blockSize = params->get<int> ("Block Size");
    } catch (InvalidParameter&) {
      // Do nothing; leave default value
    }
    int defQuorum = -1;
    try {
      defQuorum = params->get<int> ("Deflation Quorum");
    } catch (InvalidParameter&) {
      // Do nothing; leave default value
    }
    bool showMaxResNormOnly = false;
    try {
      showMaxResNormOnly = 
	params->get<bool> ("Show Maximum Residual Norm Only");
    } catch (InvalidParameter&) {
      // Do nothing; leave default value
    }

    return gmresTest (convTol, maxIterCount, haveLeftPreconditioner,
		      implicitScaleType, explicitScaleType, blockSize,
		      defQuorum, showMaxResNormOnly);
  }


  template<class Scalar, class MV, class OP>
  bool
  StatusTestFactory<Scalar, MV, OP>::
  changeMaxNumIters (const Teuchos::RCP<typename StatusTestFactory<Scalar, MV, OP>::base_test>& test,
		     const int maxIterCount)
  {
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::nonnull;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // We declare "success" if the test or at least one of its
    // children (if it's a combo_test) is a max_iter_test.
    bool success = false;
    RCP<max_iter_test> maxIterTest = rcp_dynamic_cast<max_iter_test> (test);
    if (nonnull (maxIterTest))
      {
	test->setMaxIters (maxIterCount);
	success = true;
      }
    else
      {
	RCP<combo_test> comboTest = rcp_dynamic_cast<combo_test> (test);
	if (nonnull (comboTest))
	  {
	    typedef typename combo_test::st_vector st_vector;
	    typedef typename st_vector::size_type size_type;
	    st_vector tests = test->getStatusTests ();
	    // We could use boost lambda to remove this for loop.
	    // Standard STL doesn't have a version of mem_fun for
	    // member functions that take > 1 argument.
	    for (size_type k = 0; result || k < tests.end(); ++k)
	      { // Recurse on all children, since it's possible for
		// more than one child to be a max_iter_test.
		const bool result =
		  changeMaxNumIters (tests[k], maxIterCount);
		success = result || success;
	      }
	  }
      }
    return success;
  }

  template<class Scalar, class MV, class OP>
  std::pair<bool, bool>
  StatusTestFactory<Scalar, MV, OP>::
  changeConvTolAndMaxNumIters (const Teuchos::RCP<typename StatusTestFactory<Scalar, MV, OP>::base_test>& test,
			       const typename StatusTestFactory<Scalar, MV, OP>::magnitude_type convTol,
			       const int maxIterCount)
  {
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::nonnull;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef StatusTestResNorm<Scalar,MV,OP> res_norm_base_test;
    RCP<max_iter_test> maxIterTest = rcp_dynamic_cast<max_iter_test> (test);

    // We declare "success" if we encountered a res_norm_base_test
    // _and_ a max_iter_test somewhere along the recursion path.
    bool foundResNormTest = false;
    bool foundMaxIterTest = false;

    RCP<res_norm_base_test> normTest = 
      rcp_dynamic_cast<res_norm_base_test> (test);
    if (nonnull (normTest))
      { 
	// NOTE (mfh 03 Mar 2011) setTolerance() returns an int
	// result.  However, all subclasses' implementations return 0
	// here, and all of them always (re)set the tolerance, so I
	// think it's OK to ignore the result.
	(void) test->setTolerance (convTol);
	foundResNormTest = true;
      }
    else 
      {
	RCP<max_iter_test> maxIterTest = 
	  rcp_dynamic_cast<max_iter_test> (test);
	if (nonnull (maxIterTest))
	  {
	    test->setMaxIters (maxIterCount);
	    foundMaxIterTest = true;
	  }
      }
    if (! foundResNormTest && ! foundMaxIterTest)
      {
	RCP<combo_test> comboTest = rcp_dynamic_cast<combo_test> (test);
	if (nonnull (comboTest))
	  {
	    typedef typename combo_test::st_vector st_vector;
	    typedef typename st_vector::size_type size_type;
	    st_vector tests = test->getStatusTests ();
	    // We could use boost lambda to remove this for loop.
	    // Standard STL doesn't have a version of mem_fun for
	    // member functions that take > 1 argument.
	    for (size_type k = 0; result || k < tests.end(); ++k)
	      { // Recurse on all children.
		const std::pair<bool, bool> result = 
		  changeConvTolAndMaxIters (tests[k], convTol, maxIterCount);
		foundResNormTest = result.first || foundResNormTest;
		foundMaxIterTest = result.second || foundMaxIterTest;		  
	      }
	  }
      }
    return std::make_pair (foundResNormTest, foundMaxIterTest);
  }

  template<class Scalar, class MV, class OP>
  bool
  StatusTestFactory<Scalar, MV, OP>::
  changeConvTol (const Teuchos::RCP<typename StatusTestFactory<Scalar, MV, OP>::base_test>& test,
		 const typename StatusTestFactory<Scalar, MV, OP>::magnitude_type convTol)
  {
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::nonnull;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef StatusTestResNorm<Scalar,MV,OP> res_norm_base_test;

    // We declare "success" if the test or at least one of its
    // children (if it's a combo_test) is a res_norm_base_test.
    bool success = false;
    RCP<res_norm_base_test> normTest = 
      rcp_dynamic_cast<res_norm_base_test> (test);
    if (nonnull (normTest))
      { 
	// NOTE (mfh 03 Mar 2011) setTolerance() returns an int
	// result.  However, all subclasses' implementations return 0
	// here, and all of them always (re)set the tolerance, so I
	// think it's OK to ignore the result.
	(void) test->setTolerance (convTol);
	success = true;
      }
    else
      {
	RCP<combo_test> comboTest = rcp_dynamic_cast<combo_test> (test);
	if (nonnull (comboTest))
	  {
	    typedef typename combo_test::st_vector st_vector;
	    typedef typename st_vector::size_type size_type;
	    st_vector tests = test->getStatusTests ();
	    // We could use boost lambda to remove this for loop.
	    // Standard STL doesn't have a version of mem_fun for
	    // member functions that take > 1 argument.
	    for (size_type k = 0; result || k < tests.end(); ++k)
	      { // Recurse on all children, since it's possible for
		// more than one child to be a res_norm_base_test.
		const bool result = changeConvTol (tests[k], convTol);
		success = result || success;
	      }
	  }
      }
    return success;
  }


  template<class Scalar, class MV, class OP>
  static ScaleType
  StatusTestFactory<Scalar, MV, OP>::
  stringToScaleType (const std::string& scaleType) 
  {
    const char* validNames[] = {
      "Norm of Initial Residual", 
      "Norm of Preconditioned Initial Residual",
      "Norm of RHS",
      "Norm of Right-Hand Side",
      "None"
    };
    const int numValidNames = 5;
    const ScaleType correspondingOutputs[] = {
      Belos::NormOfInitRes,
      Belos::NormOfPrecInitRes,
      Belos::NormOfRHS,
      Belos::NormOfRHS,
      Belos::None
    };
    for (int k = 0; k < numValidNames; ++k)
      {
	if (scaleType == validNames[k])
	  return correspondingOutputs[k];
      }
    TEST_FOR_EXCEPTION (true, std::logic_error,
			"Invalid residual scaling type \"" << scaleType 
			<< "\".");
  }

} // namespace Belos
