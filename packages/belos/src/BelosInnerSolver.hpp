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

#ifndef __Belos_InnerSolver_hpp
#define __Belos_InnerSolver_hpp

#include <BelosInnerSolveResult.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosMultiVecTraits.hpp>
#include <BelosOperatorTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace Belos {

  /// \brief New LinearProblem with different right-hand side.
  ///
  /// Given a LinearProblem instance, construct and return a new
  /// LinearProblem instance with the same operator,
  /// preconditioner(s), and current approximate solution, but a
  /// different right-hand side.  The current approximate solution is
  /// copied deeply, but the operator and preconditioner(s) are only
  /// copied shallowly.
  ///
  /// This is useful with InnerSolver subclasses.
  ///
  template<class Scalar, class MV, class OP>
  Teuchos::RCP<LinearProblem<Scalar, MV, OP> >
  problemWithNewRHS (const Teuchos::RCP<const LinearProblem<Scalar, MV, OP> >& problem,
		     const Teuchos::RCP<const MV>& B)
  {
    using Teuchos::is_null;
    using Teuchos::nonnull;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef LinearProblem<Scalar, MV, OP> lp_type;
    typedef MultiVecTraits<Scalar, MV> MVT;

    RCP<const OP> A = problem->getOperator (); 
    RCP<MV> X_orig = problem->getLHS ();
    TEST_FOR_EXCEPTION(is_null (X_orig), std::invalid_argument,
		       "problemWithNewRHS(): The original LinearProblem's "
		       "initial guess / current approximate solution (getLHS())"
		       " is null.  We need an initial guess or current approxim"
		       "ate solution in order to know the domain of the (right-"
		       "preconditioned, if applicable) operator.  This is "
		       "because Belos::MultiVecTraits does not include the idea"
		       " of the domain and range of an operator, or the space "
		       "to which a vector belongs.");
    TEST_FOR_EXCEPTION(is_null (B), std::invalid_argument,
		       "problemWithNewRHS(): the given new right-hand side B "
		       "is null.");
    RCP<MV> X = MVT::CloneCopy (problem->getLHS ());

    RCP<lp_type> lp (new lp_type (A, X, B));
    lp->setLeftPrec (problem->getLeftPrec ());
    lp->setRightPrec (problem->getRightPrec ());
    // Compute initial residual(s) and prepare the problem for solution.
    lp->setProblem ();
    return lp;
  }


  template<class Scalar, class MV, class OP>
  class InnerSolver {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    typedef MV multivector_type;
    typedef OP operator_type;
    typedef typename OperatorTraits<OP>::vector_space_type vector_space_type;

    virtual Teuchos::RCP<const vector_space_type> getDomain() const = 0;
    virtual Teuchos::RCP<const vector_space_type> getRange() const = 0;

    /// \brief Solve \f$AX=B\f$ for the given right-hand side(s) B.
    ///
    /// Implementations are free to interpret the convergence
    /// tolerance and maximum iteration count parameters in any way
    /// they please.  (They should interpret these parameters
    /// consistently between calls.)  They are also free to preserve
    /// more state between calls than just the matrix and
    /// preconditioner(s); for example, they may choose to perform
    /// Krylov subspace recycling.
    ///
    /// \param X [in/out] On input: The initial guess for the inner
    ///   solver, if the inner solver accepts an initial guess (it is
    ///   not required to do so).  On output: the approximate solution
    ///   to Ax=B as computed by the inner solver.  Whether or not the
    ///   solver accepts an initial guess, X must be allocated to hold
    ///   the output, and it must be in the correct vector space for
    ///   the solution vectors.
    ///
    /// \param B [in] Right-hand side(s) for which to solve
    ///
    /// \param convTol [in] "Convergence tolerance," the meaning of
    ///   which depends on the subclass
    /// \param maxItersPerRestart [in] Maximum number of iterations
    ///   per restart cycle in the inner solve.
    /// \param maxNumRestarts [in] Maximum number of restart cycle(s) 
    ///   in the inner solve.
    ///
    /// \return The result of the inner solve.  It is a single result,
    ///   aggregated over all right-hand side(s).
    ///
    virtual InnerSolveResult
    solve (const Teuchos::RCP<MV>& X,
	   const Teuchos::RCP<const MV>& B,
	   const magnitude_type convTol,
	   const int maxItersPerRestart,
	   const int maxNumRestarts) = 0;

    /// \brief Solve \f$AX=B\f$ for the given right-hand side(s) B.
    ///
    /// This should do the same thing as the five-argument version of
    /// solve(), except it should pick reasonable defaults for the
    /// convergence tolerance, maximum number of iterations, and
    /// maximum number of restart cycles.
    /// 
    /// \note Subclasses' implementations are welcome to change the
    ///   convergence tolerance, maximum number of iterations, maximum
    ///   number of restart cycles, and perhaps other solve parameters
    ///   after every invocation of solve().  For example, solve() may
    ///   implement an inexact Krylov method, in which the convergence
    ///   tolerance is gradually increased according to a computed
    ///   bound.  Callers who want to override this behavior of a
    ///   subclass' implementation should call the four-argument
    ///   version of solve().
    ///
    /// \param X [in/out] On input: The initial guess for the inner
    ///   solver, if the inner solver accepts an initial guess (it is
    ///   not required to do so).  On output: the approximate solution
    ///   to Ax=B as computed by the inner solver.  
    ///
    /// \param B [in] Right-hand side(s) for which to solve
    ///
    /// \return The result of the inner solve.  It is a single result,
    ///   aggregated over all right-hand side(s).
    ///
    virtual InnerSolveResult
    solve (const Teuchos::RCP<MV>& X,
	   const Teuchos::RCP<const MV>& B) = 0;
  };


  /// \brief Partial specialization of OperatorTraits for InnerSolver.
  ///
  /// This partial specialization lets you use InnerSolver (or any of
  /// its implementations) as an operator ("OP") in any of the Belos
  /// iterative methods.  Remember that in Belos, the preconditioner
  /// and the matrix operators must have the same type.  Thus, using
  /// InnerSolver as the OP template argument in Belos solvers would
  /// be best with unpreconditioned iterative methods.
  template<class Scalar, class MV, class OP>
  class OperatorTraits<Scalar, MV, InnerSolver<Scalar, MV, OP> > {
  public:
    static void
    Apply (const InnerSolver<Scalar, MV, OP>& Op,
	   const MV& x,
	   MV& y,
	   ETrans trans = NOTRANS)
    {
      using Teuchos::RCP;
      using Teuchos::rcpFromRef;

      TEST_FOR_EXCEPTION(trans != NOTRANS, std::invalid_argument,
			 "Belos::InnerSolver is not able to solve the "
			 "transposed system.");
      RCP<const MV> x_ptr = rcpFromRef (x);
      RCP<MV> y_ptr = rcpFromRef (y);
      (void) Op.solve (y_ptr, x_ptr);
    }

    typedef typename InnerSolver<Scalar, MV, OP>::vector_space_type vector_space_type;

    static Teuchos::RCP<const vector_space_type> 
    getDomain (const InnerSolver<Scalar, MV, OP>& A) {
      return A.getDomain ();
    }
    static Teuchos::RCP<const vector_space_type> 
    getRange (const InnerSolver<Scalar, MV, OP>& A) {
      return A.getRange ();
    }
  };

} // namespace Belos

#endif // __Belos_InnerSolver_hpp
