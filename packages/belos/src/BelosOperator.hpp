// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_OPERATOR_HPP
#define BELOS_OPERATOR_HPP

/// \file BelosOperator.hpp
///
/// \brief Alternative run-time polymorphic interface for operators.

#include "BelosConfigDefs.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosMultiVec.hpp"
#ifdef HAVE_BELOS_EXPERIMENTAL
#  include "BelosInnerSolver.hpp"
#endif // HAVE_BELOS_EXPERIMENTAL


namespace Belos {

  /// \class Operator
  /// \brief Alternative run-time polymorphic interface for operators.
  /// \author Michael Heroux and Heidi Thornquist
  ///
  /// Belos' linear solvers are templated on the scalar (Scalar),
  /// multivector (MV), and operator (OP) types.  The term "operator"
  /// includes the matrix A in the linear system \f$AX = B\f$, any
  /// left or right preconditioners. If you have enabled the corresponding 
  /// Trilinos packages, you can use Belos' solvers directly with OP =
  /// Epetra_Operator, Tpetra::Operator, or Thyra::LinearOpBase.
  /// Alternately, you may wish to use some other object as an
  /// operator.  If so, you can make that object inherit from
  /// Belos::Operator<Scalar>, and make its corresponding multivector
  /// objects inherit from Belos::MultiVec<Scalar>.  Belos' solvers
  /// may also be instantiated with MV = Belos::MultiVec<Scalar> and
  /// OP = Belos::Operator<Scalar>.
  ///
  /// A concrete implementation of this class is necessary.  Users may
  /// create their own implementation if the supplied implementations
  /// are not suitable for their needs.
  template <class ScalarType>
  class Operator {
  public:
    //! @name Constructor/Destructor
    //@{ 
    
    //! Default constructor (does nothing).
    Operator() {};
    
    //! Virtual destructor, for memory safety of derived classes.
    virtual ~Operator() {};
    //@}
    
    //! @name Methods relating to applying the operator
    //@{ 

    /// \brief Apply the operator to x, putting the result in y.
    ///
    /// Take the Belos::MultiVec \c x and apply the operator (or its
    /// transpose or Hermitian transpose) to it, writing the result
    /// into the Belos::MultiVec \c y.
    ///
    /// \param x [in] The input multivector.
    ///
    /// \param y [out] The output multivector.  x and y may not alias
    ///   (i.e., be views of) one another.
    ///
    /// \param trans [in] Whether to apply the operator (NOTRANS), its
    ///   transpose (TRANS), or its Hermitian transpose (CONJTRANS).
    ///   The default is NOTRANS.
    ///
    /// Your Operator subclass' implementation of Apply() is not
    /// required to support applying the transpose (or Hermitian
    /// transpose, if applicable).  If the caller passes in a value of
    /// \c trans which your Apply() implementation does not support,
    /// it should throw a subclass of std::exception.  In general,
    /// subclasses' implementations should signal any problems
    /// applying the operator by throwing a subclass of
    /// std::exception.
    virtual void 
    Apply (const MultiVec<ScalarType>& x, 
	   MultiVec<ScalarType>& y, 
	   ETrans trans=NOTRANS) const = 0;

    /// \brief Whether this operator implements applying the transpose.
    ///
    /// Your Operator subclass' implementation of Apply() is not
    /// required to support applying the transpose (or Hermitian
    /// transpose, if applicable).  If it <i>does</i> support applying
    /// the transpose, this method should return true.  Otherwise, it
    /// should return false.  
    ///
    /// We assume that if an operator can apply its transpose, it can
    /// also apply its Hermitian transpose, if the operator is complex
    /// (otherwise the transpose and Hermitian transpose are the same
    /// operation).
    ///
    /// We provide a default implementation of this method that
    /// conservatively returns false.  If you want your Operator
    /// subclass to advertise that it implements applying the
    /// transpose, override the default implementation in your
    /// subclass.
    virtual bool HasApplyTranspose () const {
      return false;
    }
    //@}
  };
  
  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::OperatorTraits for Belos::Operator 
  //                                               and Belos::MultiVec.
  //
  ////////////////////////////////////////////////////////////////////  
  
  /// \brief Specialization of OperatorTraits for Operator and MultiVec.
  ///
  /// This is a partial template specialization of
  /// Belos::OperatorTraits class using the Belos::Operator and
  /// Belos::MultiVec abstract interfaces.  Any class that inherits
  /// from Belos::Operator will be accepted by the Belos templated
  /// solvers, due to this specialization of Belos::OperatorTraits.
  template<class ScalarType> 
  class OperatorTraits<ScalarType, MultiVec<ScalarType>, Operator<ScalarType> > 
  {
  public:
    //! Specialization of Apply() for Operator and MultiVec objects.
    static void 
    Apply (const Operator<ScalarType>& Op, 
	   const MultiVec<ScalarType>& x, 
	   MultiVec<ScalarType>& y,
	   ETrans trans=NOTRANS)
    { 
      Op.Apply (x, y, trans); 
    }

    //! Specialization of HasApplyTranspose() for Operator objects.
    static bool
    HasApplyTranspose (const Operator<ScalarType>& Op)
    {
      return Op.HasApplyTranspose ();
    }
  };

#ifdef HAVE_BELOS_EXPERIMENTAL

  /// \class OperatorInnerSolver
  /// \brief Adaptor between InnerSolver and Belos::Operator.
  /// 
  /// This wrapper lets you use as a Belos::Operator any
  /// implementation of Belos::InnerSolver<Scalar, MV, OP> with MV =
  /// Belos::MultiVec<Scalar> and OP = Belos::Operator<Scalar>.  You
  /// may also treat this wrapper as an "envelope" by extracting the
  /// underlying InnerSolver object (which has a richer interface) and
  /// discarding the Belos::Operator wrapper.
  ///
  /// \warning This interface is experimental and therefore subject to
  ///   change or removal at any time.  Do not rely on the stability
  ///   of this interface.
  template<class Scalar>
  class OperatorInnerSolver : public Operator<Scalar> {
  public:
    typedef Scalar scalar_type;
    typedef MultiVec<Scalar> multivector_type;
    typedef Operator<Scalar> operator_type;
    typedef InnerSolver<scalar_type, multivector_type, operator_type> inner_solver_type;

    /// \brief Constructor.
    ///
    /// \param solver [in/out] The actual inner solver implementation.
    OperatorInnerSolver (const Teuchos::RCP<inner_solver_type>& solver) :
      solver_ (solver)
    {}
    //! Virtual destructor implementation, for correctness.
    virtual ~OperatorInnerSolver() {}

    /// \brief Return the underlying inner solver object.
    ///
    /// This breach of encapsulation makes this class into an
    /// "envelope."  First, the inner solver hides inside the
    /// Belos::Operator until it gets inside a Belos outer solver that
    /// recognizes the Belos::Operator as an OperatorInnerSolver.
    /// Then, the Belos outer solver can take the InnerSolver out of
    /// the "envelope," destroy the envelope (by setting its RCP to
    /// null) if it wants, and work directly with the InnerSolver.
    /// This is useful because InnerSolver's interface has the
    /// features necessary for algorithms like inexact Krylov and
    /// other inner-outer iterations, where the outer iteration is
    /// responsible for controlling the convergence tolerance and/or
    /// the cost of the inner iteration.
    ///
    /// \note This method is declared const in order to cheat
    ///   Belos::LinearProblem into letting the operator act like an
    ///   envelope.  It's technically correct to call this method
    ///   const, since it doesn't let the caller assign to the pointer
    ///   (even though it lets the caller call nonconst methods on the
    ///   InnerSolver).
    Teuchos::RCP<inner_solver_type> getInnerSolver() const {
      return solver_;
    }

    /// \brief Compute Y := alpha*solver(Y,X) + beta*Y.
    ///
    /// \note The contents of Y on input may be relevant, depending on
    ///   the inner solver implementation.  For example, Y on input
    ///   may be treated as the initial guess of an iterative solver.
    ///
    /// This function is virtual, in case you would like to override
    /// its default behavior of not implementing the transpose
    /// operation.
    virtual void 
    Apply (const multivector_type& X,
	   multivector_type& Y,
	   ETrans mode = NOTRANS) const
    {
      using Teuchos::rcpFromRef;

      TEUCHOS_TEST_FOR_EXCEPTION(mode != NOTRANS, std::invalid_argument,
			 "Belos::OperatorInnerSolver only supports applying the"
			 " operator itself, not its transpose or conjugate "
			 "transpose.");
      solver_->solve (rcpFromRef (Y), rcpFromRef (X));
    }

    /// \brief Whether this operator implements applying the transpose.
    /// 
    /// By default, we assume that it doesn't.  You may override this
    /// behavior in derived classes.
    virtual bool HasApplyTranspose() const {
      return false;
    }

  private:
    //! Default construction is not allowed.
    OperatorInnerSolver ();

    //! The inner solver implementation.
    Teuchos::RCP<inner_solver_type> solver_;
  };

  /// Specialization of makeInnerSolverOperator() for Belos::Operator.
  ///
  /// This class knows how to take an InnerSolver instance and wrap it
  /// in an implementation of the Belos::Operator interface.  That way
  /// you can use it alongside any other implementation of the
  /// Belos::Operator interface in any of the the Belos solvers.
  ///
  /// \warning This interface is experimental and therefore subject to
  ///   change or removal at any time.  Do not rely on the stability
  ///   of this interface.
  ///
  template <class Scalar>
  class InnerSolverTraits<Scalar, MultiVec<Scalar>, Operator<Scalar> > {
  public:
    typedef Scalar scalar_type;
    typedef MultiVec<scalar_type> multivector_type;
    typedef Operator<scalar_type> operator_type;
    typedef InnerSolver<scalar_type, multivector_type, operator_type> inner_solver_type;
    typedef OperatorInnerSolver<scalar_type> wrapper_type;

    /// \brief Wrap the given inner solver in a wrapper_type.
    ///
    /// The wrapper_type class implements the operator_type interface,
    /// which can be used directly in Belos.
    static Teuchos::RCP<operator_type>
    makeInnerSolverOperator (const Teuchos::RCP<inner_solver_type>& solver)
    {
      using Teuchos::rcp;
      using Teuchos::rcp_implicit_cast;
      return rcp_implicit_cast<operator_type> (rcp (new wrapper_type (solver)));
    }

    /// \brief Return the given wrapper's inner solver object.
    ///
    /// If op is an inner solver wrapper instance, return the inner
    /// solver object.  Otherwise, throw an std::bad_cast exception.
    ///
    /// \note The returned inner solver object will persist beyond the
    ///   scope of the input object.  Thus, if you don't care about
    ///   the wrapper that implements the operator_type interface, you
    ///   can get rid of the wrapper and keep the inner solver.
    static Teuchos::RCP<inner_solver_type>
    getInnerSolver (const Teuchos::RCP<operator_type>& op)
    {
      using Teuchos::RCP;
      using Teuchos::rcp_dynamic_cast;
      RCP<wrapper_type> wrapper = rcp_dynamic_cast<wrapper_type> (op, true);
      return wrapper->getInnerSolver();
    }
  };
#endif // HAVE_BELOS_EXPERIMENTAL
  
} // end Belos namespace

#endif

// end of file BelosOperator.hpp
