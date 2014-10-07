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

#ifndef __Belos_GmresBaseIteration_hpp
#define __Belos_GmresBaseIteration_hpp

/// \file BelosGmresBaseIteration.hpp
/// \brief Arnoldi iteration interface
///

#include "BelosIteration.hpp"
#include "BelosGmresBaseFactory.hpp"

namespace Belos {

  /// \class GmresBaseIteration
  /// \brief Subclass of Belos::Iteration for GMRES-based solvers
  /// \author Mark Hoemmen
  ///
  /// \warning This is EXPERIMENTAL CODE.  DO NOT RELY ON THIS CODE.
  ///   The interface or implementation may change at any time.
  ///
  /// \note Belos already has a GmresIteration subclass of Iteration,
  ///   so we named this class GmresBaseIteration, to indicate that it
  ///   is implemented using GmresBase.
  template<class Scalar, class MV, class OP>
  class GmresBaseIteration : public Iteration<Scalar, MV, OP> {
  private:
    /// \typedef impl_type
    /// \brief The implementation of all functionality of this Iteration
    typedef GmresBase<Scalar, MV, OP> impl_type;
    
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    typedef MV multivector_type;
    typedef OP operator_type;

    /// \brief Iterate until error or status test triggers
    ///
    /// \note (mfh 30 Dec 2010) The convention in Belos, for some
    /// reason, is that iterate() does not call updateSolution() on
    /// the LinearProblem instance.  In fact, getCurrentUpdate()
    /// doesn't even do that.  Iteration subclasses that implement
    /// GMRES simply don't call updateSolution() on the LinearProblem
    /// instance.  Instead, the solution manager first calls
    /// getCurrentUpdate() on the Iteration subclass instance, and
    /// then calls updateSolution(update,true) on the LinearProblem
    /// instance.  This is a reasonable convention.  
    ///
    /// \note (mfh 30 Dec 2010) Nevertheless, Belos' (block) GMRES
    /// solution manager's solve() method is messy and invasive of the
    /// Iteration subclass' state.  For example, the Iteration should
    /// be responsible for the mechanics of restarting, not the
    /// solution manager.  The Iteration has the LinearProblem
    /// instance anyway, so it should call updateSolution() before the
    /// restart.
    void 
    iterate() 
    {
      // GmresBaseIteration is an Iteration subclass, and the status
      // checks (should) only call Iteration methods to get the
      // information they want.  We add another "status test," which
      // ensures that GMRES doesn't exceed the maximum number of
      // iterations allowed (for storage of the basis vectors).  The
      // status test will also perform intermediate output via the
      // OutputManager object.
      while (stest_->checkStatus(this) != Passed && impl_->canAdvance()) {
	impl_->advance();
      }
    }

    /// Initialize the solver
    ///
    /// \note This Iteration subclass assumes that the linear problem
    /// is ready to solve, i.e., that the matrix, any
    /// preconditioner(s), and right-hand side have been set.
    void initialize();

    /// \brief Return the current iteration count.
    ///
    /// Restarting resets the current iteration count.
    ///
    /// \note The first iteration (after the initial residual
    ///   computation) counts as 1.  Belos likes to measure the
    ///   iterations while they are in progress, not after they are
    ///   done (when we increment them.  That's why we add one in the
    ///   code below.
    int getNumIters() const { return impl_->getNumIters() + 1; }
    
    /// \brief Reset the iteration count to iter.
    ///
    /// Reset the iteration count to iter, typically zero (the default
    /// value).
    ///
    /// \param iter [in] New iteration count (>= 0)
    ///
    /// \note The implementation actually saves state up to and
    /// including iteration number iter, by "backing out" state beyond
    /// that iteration.  Iterations may be resumed from there.  Thus
    /// it promises more than the Iteration interface, which reserves
    /// the right to invalidate all iteration state (???).
    void resetNumIters (int iter = 0) { 
      TEUCHOS_TEST_FOR_EXCEPTION(iter < 0, std::invalid_argument, "Belos::GmresBaseIter"
			 "ation::resetNumIters: iter = " << iter << " is "
			 "invalid; only nonnegative values are allowed.");
      impl_->backOut (iter); 
    }

    //! Restart GMRES.
    void restart () { 
      TEUCHOS_TEST_FOR_EXCEPTION(impl_.is_null(), std::logic_error, 
			 "Belos::GmresBaseIteration::restart(): GmresBase "
			 "subclass instance is null.");
      impl_->restart ();
    }
    
    /// \brief Get the residuals native to the solver.
    ///
    /// \param norms [out] norms[j] is set to the "native" residual
    /// norm for the j-th residual vector.  (This is not block
    /// Arnoldi, so blockSize==1 always, and there should only be one
    /// residual vector.)
    ///
    /// \return Optionally, a multivector with blockSize vectors
    /// containing the native residuals.  Subclasses reserve the right
    /// to return Teuchos::null.
    ///
    /// \note The Iteration interface requires that this method be
    /// const.  However, we reserve the right to cache computed
    /// results.
    /// 
    /// \note (mfh 30 Dec 2010) Belos::StatusTestGenResNorm's
    /// checkStatus() method (re)computes the norms (using a norm of
    /// its own specification, like the 1-norm, 2-norm, etc. -- but
    /// not a general norm, since it calls MVT::MvNorm() to do the
    /// work) if this method returns non-null.  If it returns null,
    /// StatusTestGenResNorm's checkStatus() method uses the "norms"
    /// parameter.  It never calls getNativeResiduals(NULL).  This is
    /// a waste, since GmresBase is always able to return a native
    /// residual vector _and_ a native residual norm.  Iteration
    /// should really have two separate pure virtual methods:
    /// getNativeResidualVector() and getNativeResidualNorm().
    Teuchos::RCP<const MV> 
    getNativeResiduals (std::vector<magnitude_type>* norms) const
    {
      using Teuchos::RCP;

      RCP<const MV> nativeResVec = impl_->currentNativeResidualVector ();
      const magnitude_type nativeResNorm = impl_->currentNativeResidualNorm ();
      if (norms != NULL)
	{
	  // This should always be 1
	  const int blockSize = getBlockSize();
	  TEUCHOS_TEST_FOR_EXCEPTION(blockSize != 1, std::logic_error,
			     "This implementation of Arnoldi/GMRES only "
			     "supports a block size of 1, but the current "
			     "block size is " << blockSize << ".");
	  // Silence the compiler error about comparing signed and unsigned.
	  if (norms->size() < (unsigned int) blockSize)
	    norms->resize (blockSize);
	  for (int k = 0; k < blockSize; ++k)
	    (*norms)[k] = nativeResNorm;
	}
      return nativeResVec;
    }

    /// Force computation of the current solution update, if it has
    /// not already yet been computed.
    ///
    /// \return Current approximate solution update
    ///
    /// \note The Iteration interface requires that this method be
    /// const, but we reserve the right to cache the current solution
    /// update.
    Teuchos::RCP<MV> getCurrentUpdate() const { 
      return impl_->getCurrentUpdate();
    }

    //! Get a constant reference to the linear problem.
    const LinearProblem<Scalar, MV, OP>& getProblem() const { return *lp_; }

    /// The block size to be used by the iterative solver in solving
    /// this linear problem.  Our GMRES/Arnoldi implementations only
    /// solve for one right-hand side at a time, so the block size is
    /// always one.
    int getBlockSize() const { return 1; }

    /// \brief Set blocksize (not allowed unless blockSize==1)
    ///
    /// Set the block size to use when solving this linear problem.
    /// The Arnoldi/GMRES implementations used by this Iteration
    /// subclass can only solve for one right-hand side at a time, so
    /// the block size must always be one.  We deal with this
    /// hackishly by raising an std::invalid_argument if a non-1 input
    /// is given.  (It's hackish because it violates the "is-a"
    /// requirement of class hierarchies; GmresBaseIteration "is not
    /// an" Iteration, because it does not support arbitrary block
    /// sizes.)
    void setBlockSize (int blockSize) {
      TEUCHOS_TEST_FOR_EXCEPTION(blockSize != 1, std::invalid_argument,
			 "Belos::GmresBaseIteration::setBlockSize: blockSize = " 
			 << blockSize << " is invalid; only blockSize = 1 is "
			 "allowed for this iteration.");
    }

    //! Whether or not the solver has been initialized
    bool isInitialized() {
      return ! impl_.is_null();
    }

    /// \brief Constructor
    ///
    /// \param problem [in/out] The linear problem to solve, along
    ///   with any preconditioner(s)
    /// \param ortho [in] Orthogonalization manager
    /// \param printer [in/out] Output manager
    /// \param tester [in] Status test for stopping iteration
    /// \param params [in(/out???)] Parameters and options for GMRES.
    ///   An important parameter, required for correctness if the
    ///   preconditioner might change from iteration to iteration, is
    ///   the "Flexible" boolean parameter which tells us whether or
    ///   not to run Flexible GMRES.
    ///
    /// \note Other Belos solvers take a MatOrthoManager instead of an
    /// OrthoManager.  This iteration doesn't use the "Mat" features,
    /// so we only take an OrthoManager.
    GmresBaseIteration (const Teuchos::RCP<LinearProblem<Scalar, MV, OP> >& problem,
			const Teuchos::RCP<const OrthoManager<Scalar, MV> >& ortho,
			const Teuchos::RCP<OutputManager<Scalar> >& printer,
			const Teuchos::RCP<StatusTest<Scalar, MV, OP> >& tester,
			const Teuchos::RCP<const Teuchos::ParameterList >& params);
    
    //! Trivial destructor
    virtual ~GmresBaseIteration() {}

    /// \brief Update the current approximate solution.
    ///
    /// Modify the LinearProblem instance's current approximate
    /// solution (the multivector returned by getLHS()) by adding in
    /// the current solution update (xUpdate_).  Tell the
    /// LinearProblem to compute a new residual vector as well.
    void
    updateSolution ()
    {
      impl_->updateSolution ();
    }

  private:
    /// Current iteration state and implementation
    /// 
    /// \note This is mutable because of some odd design choices in
    /// Belos, in particular the convention that Iteration subclasses
    /// are "stateless," even though they have internal mutable state.
    /// For example, the "initialize()" member function enforces the
    /// maintenance of internal mutable state in Iteration subclasses
    /// (otherwise, how could it initialize anything?).
    mutable Teuchos::RCP<GmresBase<Scalar, MV, OP> > impl_;

    /// \brief The linear problem to solve.
    //
    /// \note We have to keep a pointer to this because the interface
    ///   that we must implement (in order for the iteration status
    ///   checks to work) requires that initialization be separate
    ///   from construction.
    Teuchos::RCP<LinearProblem<Scalar, MV, OP> > lp_;
    /// \brief Orthogonalization manager.
    ///
    /// We're not using MatOrthoManager, because GMRES doesn't need
    /// that expanded interface.
    Teuchos::RCP<const OrthoManager<Scalar, MV> > ortho_;
    //! Output manager, for intermediate output during iterations
    Teuchos::RCP<OutputManager<Scalar> > outMan_;
    //! Status test for stopping iteration
    Teuchos::RCP<StatusTest<Scalar, MV, OP> > stest_;
    //! Configuration parameters. May be null (signifies defaults).
    Teuchos::RCP<const Teuchos::ParameterList> params_;

    /// \brief Return a validated version of the given LinearProblem.
    ///
    /// "Validated" means that the returned LinearProblem is not null,
    /// and problem->isProblemSet() returns true.
    static Teuchos::RCP<LinearProblem<Scalar, MV, OP> >
    validatedProblem (const Teuchos::RCP<LinearProblem<Scalar, MV, OP> >& problem);
  };

  template<class Scalar, class MV, class OP>
  GmresBaseIteration<Scalar,MV,OP>::
  GmresBaseIteration (const Teuchos::RCP<LinearProblem<Scalar, MV, OP> >& problem,
		      const Teuchos::RCP<const OrthoManager<Scalar, MV> >& ortho,
		      const Teuchos::RCP<OutputManager<Scalar> >& printer,
		      const Teuchos::RCP<StatusTest<Scalar, MV, OP> >& tester,
		      const Teuchos::RCP<const Teuchos::ParameterList>& params) :
    lp_ (validatedProblem (problem)),
    ortho_ (ortho), 
    outMan_ (printer), 
    stest_ (tester), 
    params_ (params)
  {
    initialize();
  }

  template<class Scalar, class MV, class OP>
  void
  GmresBaseIteration<Scalar,MV,OP>::initialize ()
  {
    if (! isInitialized())
      {
	typedef GmresBaseFactory<Scalar, MV, OP> factory_type;
	impl_ = factory_type::create (lp_, ortho_, outMan_, params_);
      }
  }

  template<class Scalar, class MV, class OP>
  Teuchos::RCP<LinearProblem<Scalar, MV, OP> >
  GmresBaseIteration<Scalar,MV,OP>::
  validatedProblem (const Teuchos::RCP<LinearProblem<Scalar, MV, OP> >& problem)
  {
    const char prefix[] = "Belos::GmresBaseIteration constructor: ";
    TEUCHOS_TEST_FOR_EXCEPTION(problem.is_null(), std::invalid_argument,
		       prefix << "The linear problem (Belos::LinearProblem "
		       "instance) that you want me to solve is null.");
    if (! problem->isProblemSet())
      {
	const bool notInitialized = problem->getOperator().is_null() || 
	  problem->getRHS().is_null() || problem->getLHS().is_null();
	TEUCHOS_TEST_FOR_EXCEPTION(notInitialized, std::invalid_argument,
			   prefix << "The given linear problem (Belos::Linear"
			   "Problem) instance is not fully initialized: "
			   << (problem->getOperator().is_null() ? "the operator A is null; " : "")
			   << (problem->getRHS().is_null() ? "the right-hand side B is null; " : "")
			   << (problem->getLHS().is_null() ? "the initial guess X is null" : "")
			   << ".");
	if (! problem->isProblemSet())
	  problem->setProblem();
	TEUCHOS_TEST_FOR_EXCEPTION(! problem->isProblemSet(), std::logic_error,
			   prefix << "Although the given LinearProblem instance "
			   "has non-null operator (matrix A), right-hand side B,"
			   " and initial guess X, and although its setProblem() "
			   "method has been called, its isProblemSet() method "
			   "returns false.  This should not happen.");
      }
    return problem;
  }

} // namespace Belos

#endif // __Belos_GmresBaseIteration_hpp
