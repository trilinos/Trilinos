// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
    TEUCHOS_TEST_FOR_EXCEPTION(is_null (X_orig), std::invalid_argument,
		       "problemWithNewRHS(): The original LinearProblem's "
		       "initial guess / current approximate solution (getLHS())"
		       " is null.  We need an initial guess or current approxim"
		       "ate solution in order to know the domain of the (right-"
		       "preconditioned, if applicable) operator.  This is "
		       "because Belos::MultiVecTraits does not include the idea"
		       " of the domain and range of an operator, or the space "
		       "to which a vector belongs.");
    TEUCHOS_TEST_FOR_EXCEPTION(is_null (B), std::invalid_argument,
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


  /// \class InnerSolver
  /// \brief Inner solver interface.
  /// \author Mark Hoemmen
  ///
  /// \warning This is EXPERIMENTAL CODE.  DO NOT RELY ON THIS CODE.
  ///   The interface or implementation may change at any time.
  ///
  /// An "inner solver" wraps an existing linear solver, especially an
  /// iterative solver (such as any of the iterative solvers
  /// implemented in Belos).  InnerSolver is designed especially for
  /// implementing the inner solve in inner-outer iterations, such as
  /// the following:
  /// - Flexible GMRES (which was originally designed as an inner-outer
  ///   iteration)
  /// - Multiprecision algorithms (where e.g., the preconditioner uses a 
  ///   different floating-point precision than the matrix and vectors)
  /// - Inexact Krylov methods (where the matrix and/or preconditioner
  ///   are applied with accuracy that gradually decreases as the
  ///   outer iteration converges)
  ///
  /// InnerSolvers may be used directly as the "OP" template argument
  /// in \c Belos::OperatorTraits, but Belos' current architecture
  /// makes this hard to use in practice.  The more typical use is via
  /// the wrapper defined in \c Belos::InnerSolveTraits.  The \c
  /// InnerSolveTraits::makeInnerSolveOperator() method wraps an
  /// InnerSolve instance in an operator appropriate for your linear
  /// algebra library: \c Epetra_Operator for Epetra, \c
  /// Tpetra::Operator for Tpetra, or \c Thyra::LinearOpBase for
  /// Thyra.  Then, you can mix and match InnerSolver instances with
  /// ordinary matrices or other kinds of operators in your Belos
  /// iterative solvers.  The wrappers have an "envelope" capability,
  /// so your custom iterative solvers can receive the InnerSolver
  /// wrapped up in the wrapper, extract the InnerSolver (via \c
  /// Belos::InnerSolveTraits::getInnerSolver()), discard the wrapper
  /// if desired, and use the InnerSolver directly.
  template<class Scalar, class MV, class OP>
  class InnerSolver {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    typedef MV multivector_type;
    typedef OP operator_type;

    //! Virtual destructor, for correctness.
    virtual ~InnerSolver() {}

    /// \brief Current parameters for the inner solver implementation.
    ///
    /// These parameters may change values in place, if the
    /// five-argument version of the \c solve() method is called.  If
    /// you want to preserve the original parameter values, make a
    /// deep copy of the returned ParameterList.
    virtual Teuchos::RCP<const Teuchos::ParameterList> 
    getCurrentParameters() const = 0;

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
    ///   the output, and it must be in the range vector space of this
    ///   solver.
    ///
    /// \param B [in] Right-hand side(s) for which to solve.  It must
    ///   be in the domain vector space of this solver.
    ///
    /// \param convTol [in] "Convergence tolerance," the meaning of
    ///   which depends on the subclass
    ///
    /// \param maxItersPerRestart [in] Maximum number of iterations
    ///   per restart cycle in the inner solve.
    ///
    /// \param maxNumRestarts [in] Maximum number of restart cycle(s) 
    ///   in the inner solve.
    ///
    /// \return The result of the inner solve.  It is a single result,
    ///   aggregated over all right-hand side(s).
    ///
    /// \note Why isn't this method const?  First, inner solves almost
    ///   certainly involve side effects (other than just modifying
    ///   the output argument X), and likely those side effects affect
    ///   the state of the inner solver's implementation.  Second, \c
    ///   solve() may not be a pure function.  For example, if the
    ///   inner solve implementation is a recycling Krylov solver, it
    ///   may compute a different (hopefully better) approximate
    ///   solution if given the same inputs twice in a row, or at
    ///   least it may converge in fewer iterations the second time.
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
    ///   to Ax=B as computed by the inner solver.  Whether or not the
    ///   solver accepts an initial guess, X must be allocated to hold
    ///   the output, and it must be in the range vector space of this
    ///   solver.
    ///
    /// \param B [in] Right-hand side(s) for which to solve.  It must
    ///   be in the domain vector space of this solver.
    ///
    /// \return The result of the inner solve.  It is a single result,
    ///   aggregated over all right-hand side(s).
    ///
    /// \note Why isn't this method const?  First, inner solves almost
    ///   certainly involve side effects (other than just modifying
    ///   the output argument X), and likely those side effects affect
    ///   the state of the inner solver's implementation.  Second, \c
    ///   solve() may not be a pure function.  For example, if the
    ///   inner solve implementation is a recycling Krylov solver, it
    ///   may compute a different (hopefully better) approximate
    ///   solution if given the same inputs twice in a row, or at
    ///   least it may converge in fewer iterations the second time.
    ///   Third, the two-argument version of \c solve() reserves the
    ///   right to modify the stopping criteria on each call.
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

      TEUCHOS_TEST_FOR_EXCEPTION(trans != NOTRANS, std::invalid_argument,
			 "Belos::InnerSolver is not able to solve the "
			 "transposed system.");
      RCP<const MV> x_ptr = rcpFromRef (x);
      RCP<MV> y_ptr = rcpFromRef (y);
      (void) Op.solve (y_ptr, x_ptr);
    }

  };

  /// \class UndefinedWrapperType
  ///
  /// Undefined wrapper type, to check at compile time whether
  /// InnerSolverTraits has been specialized.  Specializations
  /// (partial or total) of InnerSolverTraits do not refer to this
  /// class.  (Deliberate) compiler errors will occur whenever the
  /// compiler attempts to instantiate the constructor of this class.
  template<class Scalar, class MV, class OP>
  class UndefinedWrapperType {
  public:
    UndefinedWrapperType() {
      (void) Scalar::this_specialization_is_not_defined();
      (void) MV::this_specialization_is_not_defined();
      (void) OP::this_specialization_is_not_defined();
    }
  };

  /// \brief Wrap an InnerSolver in an OP (operator).
  ///
  /// Take an InnerSolver instance, and wrap it in an implementation
  /// of OP, where OP is one of the base operator interfaces that
  /// Belos supports (currently Thyra::LinearOpBase, Tpetra::Operator,
  /// and Epetra_Operator).
  ///
  /// \note The reason for this class is that Belos solvers require
  ///   that the preconditioner(s) and the operator have the same type
  ///   (OP).  So, if we want to use the InnerSolver as e.g., a
  ///   preconditioner for a Tpetra::CrsMatrix (which implements
  ///   Tpetra::Operator), we have to wrap the InnerSolver in
  ///   something that implements the Tpetra::Operator interface.
  ///   Otherwise, we would just use the above definition of
  ///   OperatorTraits for InnerSolver above.
  /// 
  /// \note The reason this is a class and not just a template
  ///   function, is that C++ doesn't allow partial template
  ///   specialization of template functions.
  template<class Scalar, class MV, class OP>
  class InnerSolverTraits {
  public:
    typedef Scalar scalar_type;
    typedef MV multivector_type;
    typedef OP operator_type;
    typedef InnerSolver<scalar_type, multivector_type, operator_type> inner_solver_type;
    typedef UndefinedWrapperType<Scalar, MV, OP> wrapper_type;

    /// \brief Wrap the given inner solver in a wrapper_type.
    ///
    /// The wrapper_type class implements the operator_type interface,
    /// which can be used directly in Belos.
    static Teuchos::RCP<OP>
    makeInnerSolverOperator (const Teuchos::RCP<InnerSolver<Scalar, MV, OP> >& solver)
    {
      using Teuchos::rcp;
      using Teuchos::rcp_implicit_cast;
      // If this class is not specialized for the given combination of
      // (Scalar, MV, OP), the constructor of wrapper_type here will
      // (deliberately) raise a compiler error.
      return rcp_implicit_cast<operator_type> (rcp (new wrapper_type (solver)));
    }

    /// \brief Return the given wrapper's inner solver object.
    ///
    /// If op is an inner solver wrapper instance, return the inner
    /// solver object.  Otherwise, throw an std::bad_cast exception.  
    ///
    /// \note After calling this method, the inner solver object will
    ///   persist beyond the scope of op.  Thus, if you don't care
    ///   about the wrapper that implements the operator_type
    ///   interface, you can get rid of the wrapper (by setting the
    ///   RCP to null) and keep the inner solver.
    static Teuchos::RCP<inner_solver_type>
    getInnerSolver (const Teuchos::RCP<operator_type>& op)
    {
      using Teuchos::RCP;
      using Teuchos::rcp_dynamic_cast;
      // If this class is not specialized for the given combination of
      // (Scalar, MV, OP), the instantiation of the wrapper_type class
      // here will (deliberately) raise a compiler error.
      RCP<wrapper_type> wrapper = rcp_dynamic_cast<wrapper_type> (op, true);
      return wrapper->getInnerSolver();
    }
  };

} // namespace Belos

#endif // __Belos_InnerSolver_hpp
