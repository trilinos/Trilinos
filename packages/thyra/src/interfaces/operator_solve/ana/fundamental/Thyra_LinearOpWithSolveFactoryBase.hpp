// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DECL_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DECL_HPP

#include "Thyra_LinearOpWithSolveBaseDecl.hpp"

namespace Thyra {

/** \brief Factory interface for creating <tt>LinearOpWithSolveBase</tt>
 * objects from <tt>LinearOpBase</tt> objects.
 *
 * \ingroup Thyra_Op_Vec_Interoperability_Extended_Interfaces_grp
 *
 * \section LOWSFB_outline_sec Outline
 *
 * <ul>
 * <li>\ref LOWSFB_intro_sec
 * <li>\ref LOWSFB_use_cases_sec
 * <li>\ref LOWSFB_developer_notes_sec
 * </ul>
 *
 * \section LOWSFB_intro_sec Introduction
 *
 * This strategy interface allows a client to take one or more "compatible"
 * <tt>LinearOpBase</tt> objects and then create one or more
 * <tt>LinearOpWithSolveBase</tt> objects that can then be used to solve for
 * linear systems.  This interface carefully separates the construction from
 * the initialization of a <tt>LinearOpWithSolveBase</tt> object.
 *
 * Note that the non-member functions defined \ref
 * Thyra_LinearOpWithSolveFactoryBase_helper_grp "here" provide for simpler
 * use cases and are demonstrated in the example code below.
 *
 * This interface can be implemented by both direct and iterative linear
 * solvers.
 *
 * This interface supports the concept of an external preconditioner that can
 * be created by the client an passed into the function
 *
 * \section LOWSFB_use_cases_sec Use cases
 *
 * The following use cases demonstrate and specify the behavior of this
 * interface in several different situations.  These use cases don't cover
 * very possible variation but implementors and clients should get a pretty
 * good idea of the behavior specified here.
 *
 * The following use cases are described below:
 *
 * <ul>
 * <li>\ref LOWSFB_single_linear_solve_sec
 * <li>\ref LOWSFB_scaled_adjoint_sec
 * <li>\ref LOWSFB_updated_of_linear_op_sec
 * <li>\ref LOWSFB_reuse_of_factorization_sec
 * <li>\ref LOWSFB_major_changes_in_op_sec
 * <li>\ref LOWSFB_external_prec_sec
 *   <ul>
 *   <li>\ref LOWSFB_external_prec_op_sec
 *   <li>\ref LOWSFB_external_prec_mat_sec
 *   <li>\ref LOWSFB_external_prec_op_reuse_sec
 *   <li>\ref LOWSFB_external_prec_mat_reuse_sec
 *   </ul>
 * </ul>
 *
 * \subsection LOWSFB_single_linear_solve_sec Performing a single linear solve with the forward operator
 *
 * Performing a single linear solve is at minimum a two step process but can
 * be performed in a single line of code.  The following example function
 * shows how to take a <tt>LinearOpBase</tt> object, a compatible
 * <tt>LinearOpWithSolveFactoryBase</tt> object, a RHS and a LHS and then
 * solve the linear system using a default solve:

 \code

  template<class Scalar>
  void singleLinearSolve(
    const Thyra::LinearOpBase<Scalar>                    &A
    ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>   &opFactory
    ,const Thyra::VectorBase<Scalar>                     &b
    ,Thyra::VectorBase<Scalar>                           *x
    )
  {
    // Create the LOWSB object that will be used to solve the linear system
    Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<Scalar> >
      invertibleA = createAndInitializeLinearOpWithSolve(opFactory,Teuchos::rcp(&A,false));
    // Solve the system using a default solve criteria using a non-member helper function 
    assign(&*x,Teuchos::ScalarTraits<Scalar>::zero()); // Must initialize to a guess before solve!
    solve(*invertibleA,NOTRANS,b,x);
  }

 \endcode

 * See the documentation for the <tt>LinearOpWithSolveBase</tt> interface for
 * how to specify solve criteria.
 *
 * Note that the forward operator <tt>A</tt> is passed in as a raw object
 * reference and not as a <tt>Teuchos::RefCountPtr</tt> wrapped object since
 * no persisting relationship with this object will be last after this
 * function exists, even if an exception is thrown.
 *
 * Also note that once the above function exits that all memory of the
 * invertible operator created as part of the <tt>invertibleA</tt> object will
 * be erased.  The following use cases show more sophisticated uses of these
 * interfaces.
 *
 * \subsection LOWSFB_scaled_adjoint_sec Creating invertible operators for scaled and/or adjoint forward operators
 *
 * This interface requires that all good implementations support implicitly
 * scaled and adjoint (or transposed) forward operators.  The following example
 * function shows how this looks:

 \code

  template<class Scalar>
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> 
  createScaledAdjoint(
    Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >     &A          // A persisting relationship!
    ,const Scalar                                                scalar
    ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>           &opFactory
    )
  {
    return createAndInitializeLinearOpWithSolve(opFactory,scale(scalar,adjoint(A)));
  }

 \endcode
 
 * In the above example, the functions <tt>adjoint()</tt> and <tt>scale()</tt>
 * create an implicitly scaled and adjoint operator of type
 * <tt>ScaledAdjointLinearOp</tt> which is then unwrapped by the
 * <tt>opFactory</tt> implementation.  The idea is that any operation that
 * works with a particular forward operator <tt>A</tt> should automatically
 * work with an implicitly scaled and/or adjoint view of that forward
 * operator.  The specification of this interface actually says that all
 * <tt>LinearOpBase</tt> objects that support the abstract mix-in interface
 * <tt>ScaledAdjointLinearOpBase</tt> should be supported.
 *
 * Above also note that the forward operator <tt>A</tt> is passed in as a
 * <tt>Teuchos::RefCountPtr</tt> wrapped object since it will be used to
 * create a persisting relationship with the returned
 * <tt>Thyra::LinearOpWithSolveBase</tt> object.  Also note that the
 * <tt>opFactory</tt> object is still passed in as a raw object reference
 * since no persisting relationship with <tt>opFactory</tt> is created as a
 * side effect of calling this function.  Remember, the specification of this
 * interface requires that the returned <tt>LinearOpWithSolveBase</tt> be
 * independent from the the <tt>opFactory</tt> object that creates it.
 *
 * \subsection LOWSFB_updated_of_linear_op_sec Updates of linear operator between linear solves
 *
 * In this use case, we show how to:
 * <ul>
 * <li> Create a <tt>LinearOpWithSolveBase</tt> object from a <tt>LinearOpBase</tt> objects
 * <li> Use the <tt>LinearOpWithSolveBase</tt> object to solve one or more linear systems
 * <li> Modify the <tt>LinearOpBase</tt> object and reinitialize the <tt>LinearOpWithSolveBase</tt> object
 * <li> Use the updated <tt>LinearOpWithSolveBase</tt> object to solve one or more linear systems
 * </ul>
 *
 * The below code fragment shows how this looks.

 \code

  template<class Scalar>
  void solveChangeSolve(
    const Thyra::LinearOpBase<Scalar>                    &A
    ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>   &opFactory
    )
  {
    // Get a local non-owned RCP to A to be used by opFactory
    Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
      rcpA = Teuchos::rcp(&A,false);
    // Create the LOWSB object that will be used to solve the linear system
    Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<Scalar> >
      invertibleA = opFactory.createOp();
    // Initialize the invertible linear operator given the forward operator
    opFactory.initializeOp(rcpA,&*invertibleA);
    // Solve the system using a default solve criteria using a non-member helper function
    solve(*invertibleA,...);
    // Change the operator and reinitialize the invertible operator
    opFactory.uninitializeOp(&*invertibleA);
      // Note that the above uninitialization is not required but
      // is recommended before a change of A since the object invertibleA may
      // become invalid and this avoids accidental use until it is reinitialized
      // below.
    someChange(&A);
    opFactory.initializeOp(rcpA,&*invertibleA);
      // Note that above will reuse any factorization structures that may have been
      // created in the first call to initializeOp(...)
    // Solve another linear system with new values of A
    solve(*invertibleA,...);
  }

 \endcode

 * In the above code fragment the call to the function
 * <tt>opFactory.uninitializeOp(&*invertibleA)</tt> may not fully uninitialize
 * the <tt>*invertibleA</tt> object as it may contain memory of a
 * factorization structure, or a communication pattern, or whatever may have
 * been computed in the first call to
 * <tt>opFactory.initializeOp(rcpA,&*invertibleA)</tt>.  This allows for a
 * more efficient reinitialization on the second call of
 * <tt>opFactory.initializeOp(rcpA,&*invertibleA)</tt>.
 *
 * \subsection LOWSFB_reuse_of_factorization_sec Reuse of factorizations for small changes in the forward operator
 *
 * This interface supports the notion of the reuse of all factorizations or
 * other expensive preprocessing used in the last initialization of the
 * <tt>LinearOpWithSolveBase</tt> object.  The primary use for this is the
 * reuse of a preconditioner for small changes in the matrix values.  While
 * this approach generally is not very effective in many cases, there are some
 * cases, such as in transient solvers, where this has been shown to be
 * effective in some problems.
 *
 * The below example function shows what this looks like:

 \code

  template<class Scalar>
  void solveSmallChangeSolve(
    const Thyra::LinearOpBase<Scalar>                    &A
    ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>   &opFactory
    )
  {
    // Get a local non-owned RCP to A to be used by opFactory
    Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
      rcpA = Teuchos::rcp(&A,false);
    // Create the LOWSB object that will be used to solve the linear system
    Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<Scalar> >
      invertibleA = opFactory.createOp();
    // Initialize the invertible linear operator given the forward operator
    opFactory.initializeOp(rcpA,&*invertibleA);
    // Solve the system using a default solve criteria using a non-member helper function
    solve(*invertibleA,...);
    // Change the operator and reinitialize the invertible operator
    opFactory.uninitializeOp(&*invertibleA);
      // Note that the above uninitialization is not required but
      // is recommended before a change of A since the object invertibleA may
      // become invalid and this avoids accidental use until it is reinitialized
      // below.
    someSmallChange(&A);
    opFactory.initializeAndReuseOp(rcpA,&*invertibleA);
      // Note that above will reuse any factorization values (i.e. preconditioners)
      // that may have been created in the first call to initializeOp(...)
    // Solve another linear system with new values of A
    solve(*invertibleA,...);
  }

 \endcode
 
 *
 * \subsection LOWSFB_major_changes_in_op_sec Major changes in the structure of the forward operator
 *
 * A major change in the shape or structure of a forward operator generally
 * eliminates the possibility of reusing any computations between different
 * calls to <tt>initializeOp()</tt>.  In these instances, the client might as
 * well just recreate the <tt>LinearOpWithSolveBase</tt> using <tt>createOp()</tt>
 * before the reinitialization.
 *
 * The below example function shows what this looks like:

 \code

  template<class Scalar>
  void solveMajorChangeSolve(
    const Thyra::LinearOpBase<Scalar>                    &A
    ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>   &opFactory
    )
  {
    // Get a local non-owned RCP to A to be used by opFactory
    Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
      rcpA = Teuchos::rcp(&A,false);
    // Create the LOWSB object that will be used to solve the linear system
    Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<Scalar> >
      invertibleA = opFactory.createOp();
    // Initialize the invertible linear operator given the forward operator
    opFactory.initializeOp(rcpA,&*invertibleA);
    // Solve the system using a default solve criteria using a non-member helper function
    solve(*invertibleA,...);
    // Change the operator and initialize the invertible operator anew
    invertibleA = Teuchos::null;
      // Note that the above whip out of the invertibleA object is not required but
      // is recommended before a change of A since the object invertibleA may
      // become invalid and this avoids accidental use until it is recreated below.
    someMajorChange(&A);
    invertibleA = opFactory.createOp(); // Starting from scratch
    opFactory.initializeOp(rcpA,&*invertibleA);
      // Note that above will reuse any factorization values (i.e. preconditioners)
      // that may have been created in the first call to initializeOp(...)
    // Solve another linear system with new values of A
    solve(*invertibleA,...);
  }

 \endcode

 *
 * \subsection LOWSFB_external_prec_sec The use of externally defined preconditioners
 *
 * This interface also supports the use of externally defined preconditioners
 * that are created and controlled by the client.  Client-created and
 * client-controlled preconditioners can be passed along with the forward
 * operator through the function <tt>initializePreconditionedOp()</tt>.  Only
 * objects that return <tt>supportsPreconditionerType()==true</tt> will
 * support externally defined preconditioners.  In general, iterative linear
 * solver implementations will support externally defined preconditioners and
 * direct linear solver implementations will not.
 *
 * Externally defined preconditioners can be passed in as operators or as matrices
 * and these two sub use cases are described next.
 *
 * \subsubsection LOWSFB_external_prec_op_sec The use of externally defined preconditioner operators
 *
 * In this use case, the preconditioner is passed in as an operator to be
 * directly applied to as an operator.
 *
 * The following example function shows how a <tt>LinearOpWithSolveBase</tt>
 * object is created and initialized using an externally defined
 * preconditioner operator:

 \code

  template<class Scalar>
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> 
  createOperatorPreconditioned(
    Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >     &A          // A persisting relationship!
    Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >     &Pop        // A persisting relationship!
    ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>           &opFactory
    )
  {
    Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<Scalar> >
      invertibleA = opFactory.createOp();
    opFactory.initializePreconditionedOp(A,Pop,Thyra::PRECONDITIONER_INPUT_TYPE_AS_OPERATOR,&*invertibleA);
    return invertibleA;
  }

 \endcode

 * As show above, the type of the preconditioner is flagged using the enum
 * value of <tt>Thyra::PRECONDITIONER_INPUT_TYPE_AS_OPERATOR</tt>.  This allows the
 * <tt>opFactory</tt> implementation to determine how the interpret the input
 * preconditioner object <tt>P</tt>.  Once this <tt>invertibleA</tt> object
 * is returned, it can be used just like any other
 * <tt>LinearOpWithSolveBase</tt> object.
 *
 * \subsubsection LOWSFB_external_prec_mat_sec The use of externally defined preconditioner matrices
 *
 * In this use case, the preconditioner is passed in as a matrix object to be
 * used to generate (using an algebraic method such as ILU or some multi-level
 * method) preconditioner operator that then is used as the preconditioner.
 *
 * The following example function shows how a <tt>LinearOpWithSolveBase</tt>
 * object is created and initialized using an externally defined
 * preconditioner matrix:

 \code

  template<class Scalar>
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> 
  createOperatorPreconditioned(
    Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >     &A          // A persisting relationship!
    Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >     &Pmat       // A persisting relationship!
    ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>           &opFactory
    )
  {
    Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<Scalar> >
      invertibleA = opFactory.createOp();
    opFactory.initializePreconditionedOp(A,Pmat,Thyra::PRECONDITIONER_INPUT_TYPE_AS_MATRIX,&*invertibleA);
    return invertibleA;
  }

 \endcode

 * As show above, the type of the preconditioner is flagged using the enum
 * value of <tt>Thyra::PRECONDITIONER_INPUT_TYPE_AS_MATRIX</tt>.  This allows the
 * <tt>opFactory</tt> implementation to determine how the interpret the input
 * preconditioner object <tt>P</tt>.  This mode is very similar to the default
 * mode as supported by <tt>initializeOp()</tt> where the preconditioner is
 * generated from the forward matrix itself.  Except here a different matrix
 * is used to generate the preconditioner. This would be the case where the
 * forward operator is not a matrix (i.e. it can only be applied) and a
 * different, perhaps simplified physics-based, matrix approximation is used to
 * generate the preconditioner.
 *
 * \subsubsection LOWSFB_external_prec_op_reuse_sec Reuse of externally defined preconditioner operators
 *
 * Reusing an externally defined preconditioner operator is a very
 * straightforward matter.  Since the client controls the creation and
 * initialization of the preconditioner operator, this can be passed into
 *
 * The below example function shows what it looks like to reuse a
 * preconditioner operator:

 \code

  template<class Scalar>
  void preconditionedSolveChangeSolve(
    const Thyra::LinearOpBase<Scalar>                    &A
    const Thyra::LinearOpBase<Scalar>                    &P
    ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>   &opFactory
    )
  {
    // Get a local non-owned RCP  objects to A and P to be used by opFactory
    Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
      rcpA = Teuchos::rcp(&A,false),
      rcpP = Teuchos::rcp(&A,false);
    // Create the LOWSB object that will be used to solve the linear system
    Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<Scalar> >
      invertibleA = opFactory.createOp();
    // Initialize the preconditioned invertible linear operator given the
    // forward and preconditioner operators
    opFactory.initializePreconditionedOp(rcpA,rcpP,Thyra::PRECONDITIONER_INPUT_TYPE_AS_OPERATOR,&*invertibleA);
    // Solve the system using a default solve criteria using a non-member helper function
    solve(*invertibleA,...);
    // Change the operator and reinitialize the invertible operator
    opFactory.uninitializeOp(&*invertibleA);
      // Note that the above uninitialization is not required but
      // is recommended before a change of A since the object invertibleA may
      // become invalid and this avoids accidental use until it is reinitialized
      // below.
    someSmallChange(&A);
    opFactory.initializePreconditionedOp(rcpA,rcpP,Thyra::PRECONDITIONER_INPUT_TYPE_AS_OPERATOR,&*invertibleA);
      // Note that above we assume that P is independent from A so that when A
      // is changed it will have not impact on P.  If P was dependent on A,
      // then the client would need to know what the side effects were and if
      // this was reasonable or not.
    // Solve another linear system with new values of A but same preconditioner P
    solve(*invertibleA,...);
  }

 \endcode

 *
 * \subsubsection LOWSFB_external_prec_mat_reuse_sec Reuse of externally defined preconditioner matrices
 *
 * This interface does not guarantee the correct reuse of externally defined
 * preconditioner matrices.  The problem is that the internally generated
 * preconditioner may be dependent on the input preconditioner matrix and it
 * is not clear for the current interface design how to cleanly handle this
 * use case.  However, when the client knows that the internally generated
 * preconditioner operator is independent from the externally defined
 * preconditioner input matrix, then the function
 * <tt>initializeAndReuseOp()</tt> can be called to reuse the internal
 * preconditioner but this is implementation defined.
 * 
 * \section LOWSFB_developer_notes_sec Notes to subclass developers
 *
 * This interface assumes a minimal default set of functionality that is
 * appropriate for direct and simple itertive linear solver implementations.
 * The pure virtual functions that must be overridden by a subclass are
 * <tt>isCompatible()</tt>, <tt>createOp()</tt>, <tt>initializeOp()</tt>, and
 * <tt>uninitializeOp()</tt>.  By far the most complex function to implement
 * is <tt>initializeOp()</tt> which is where the real guts of the factory
 * method is defined.
 *
 * If the concrete subclass can support some significant type of preprocessing
 * reuse then it may override to function <tt>initializeAndReuseOp()</tt>.
 * The most common use of this function it to support the reuse of
 * preconditioners between small changes in forward operator matrix values.
 *
 * If the concrete subclass can utilize a preconditioner, then it should
 * override the functions <tt>supportsPreconditionerInputType()</tt> and
 * <tt>initializePreconditionedOp()</tt>.  The subclass implementation can
 * decide if preconditioner operators and/or matrices are supported or not and
 * this is determined by the return value from
 * <tt>supportsPreconditionerInputType()</tt>.
 */
template <class RangeScalar, class DomainScalar = RangeScalar>
class LinearOpWithSolveFactoryBase : virtual public Teuchos::Describable {
public:

  /** @name Pure virtual public functions that must be overridden in subclasses */
  //@{

  /** \brief Check that a <tt>LinearOpBase</tt> object is compatible with
   * <tt>*this</tt> factory object.
   */
  virtual bool isCompatible( const LinearOpBase<RangeScalar,DomainScalar> &fwdOp ) const = 0;

  /** \brief Create an (uninitialized) <tt>LinearOpWithSolveBase</tt> object
   * to be initialized later in <tt>this->initializeOp()</tt>.
   *
   * Note that on output <tt>return->domain().get()==NULL</tt> may be true
   * which means that the operator is not fully initialized.  In fact, the
   * output operator object is not guaranteed to be fully initialized until
   * after it is passed through <tt>this->initializeOp()</tt>.
   */
  virtual Teuchos::RefCountPtr<LinearOpWithSolveBase<RangeScalar,DomainScalar> > createOp() const = 0;

  /** \brief Initialize a pre-created <tt>LinearOpWithSolveBase</tt> object
   * given a "compatible" <tt>LinearOpBase</tt> object.
   *
   * \param  fwdOp  [in] The forward linear operator that will be used to create
   *                the output <tt>LinearOpWithSolveBase</tt> object.  Note that this object is remembered
   *                by the <tt>*Op</tt> object on output.
   * \param  Op     [in/out] The output <tt>LinearOpWithSolveBase</tt> object.  This object must have
   *                be created first by <tt>this->createOp()</tt>.  The object may have also
   *                already been passed through this function several times.  Note that subclasses
   *                should always first strip off the transpose and scaling by calling <tt>unwrap()</tt>
   *                before attempting to dynamic cast the object.
   * \param  supportSolveUse
   *                [in] Determines if <tt>Op->solve(...)</tt> or <tt>Op->solveTranspose(...)</tt> will
   *                be called.  This allows <tt>*this</tt> factory object determine how to best initialae
   *                the <tt>*Op</tt> object.  Default <tt>supportSolveUse=SUPPORT_SOLVE_UNSPECIFIED</tt>
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>fwdOp.get()!=NULL</tt>
   * <li><tt>this->isCompatible(*fwdOp)==true</tt>
   * <li><tt>Op!=NULL</tt>
   * <li><tt>*Op</tt> must have been created by <tt>this->createOp()</tt> prior to calling
   *     this function.
   * <li>[<tt>supportSolveUse==SUPPORT_SOLVE_FORWARD_ONLY<tt>]
   *     <tt>this->solveSupportsConj(conj)==true</tt> for any value of <tt>conj</tt>
   * <li>[<tt>supportSolveUse==SUPPORT_SOLVE_TRANSPOSE_ONLY<tt>]
   *     <tt>this->solveTransposeSupportsConj(conj)==true</tt> for any value of <tt>conj</tt>
   * <li>[<tt>supportSolveUse==SUPPORT_SOLVE_FORWARD_AND_TRANSPOSE<tt>]
   *     <tt>this->solveSupportsConj(conj)==true && this->solveTransposeSupportsConj(conj)==true</tt>
   *     for any value of <tt>conj</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>Throws <tt>CatastrophicSolveFailure</tt> if the underlying linear solver could
   *     not be created successfully (e.g. due to a factorization failure or some other cause).
   * <li><tt>Op->range()->isCompatible(*fwdOp->range())==true</tt>
   * <li><tt>Op->domain()->isCompatible(*fwdOp->domain())==true</tt>
   * <li><tt>Op->apply()</tt> and <tt>Op->applyTranspose()</tt> must behave
   *     exactly the same as <tt>fwdOp->apply()</tt> and <tt>fwdOp->applyTranspose()</tt>
   * <li><t>Op->solveSupportsConj(conj)==this->solveSupportsConj(conj)</tt>
   * <li><t>Op->solveTransposeSupportsConj(conj)==this->solveTransposeSupportsConj(conj)</tt>
   * <li><t>Op->solveSupportsSolveTolType(conj,solveTolType)==this->solveSupportsSolveTolType(conj,solveTolType)</tt>
   * <li><t>Op->solveTransposeSupportsSolveTolType(conj,solveTolType)==this->solveTransposeSupportsSolveTolType(conj,solveTolType)</tt>
   * <li><tt>fwdOp.count()</tt> after output is greater than <tt>fwdOp.count()</tt>
   *     just before this call and therefore the client can assume that the <tt>*fwdOp</tt> object will 
   *     be remembered by the <tt>*Op</tt> object.  The client must be careful
   *     not to modify the <tt>*fwdOp</tt> object or else the <tt>*Op</tt> object may also
   *     be modified and become invalid.
   * </ul>
   */
  virtual void initializeOp(
    const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &fwdOp
    ,LinearOpWithSolveBase<RangeScalar,DomainScalar>                             *Op
    ,const ESupportSolveUse                                                      supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
    ) const = 0;

  /** \brief Uninitialize a <tt>LinearOpWithSolveBase</tt> object and return
   * its remembered forward linear operator and potentially also its
   * externally generated preconditioner.
   *
   * \param  Op     [in/out] On input, <tt>*Op</tt> is an initialized or uninitialized
   *                object and on output is uninitialized.  Note that "uninitialized"
   *                does not mean that <tt>Op</tt> is completely stateless.  It may still
   *                remember some aspect of the matrix <tt>fwdOp</tt> that will allow
   *                for a more efficient initialization next time through
   *                <tt>this->initializeOp()</tt>.
   * \param  fwdOp  [in/out] If <tt>fwdOp!=NULL</tt> on input, then on output this is set to the
   *                same forward operator passed into <tt>this->initializeOp()</tt>.
   * \param  precOp [in/out] If <tt>precOp!=NULL</tt> on input, then on output, this is set to
   *                same preconditioner operator passed into <tt>this->initializePreconditionedOp()</tt> will be
   *                returned if the preconditioner was not ignored.
   * \param  precOpType
   *                [in/out] If <tt>precOpType!=NULL</tt> on input, then on output this is set to
   *                same option value passed to <tt>this->initializePreconditionedOp()</tt> will be
   *                returned if the preconditioner was not ignored.
   * \param  ESupportSolveUse
   *                [in/out] If <tt>fwdOp!=NULL</tt> on input, then on output this is set to
   *                same option value passed to <tt>this->initializeOp()</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>*Op</tt> must have been created by <tt>this->createOp()</tt> prior to calling
   *     this function.
   * <li><tt>Op</tt> may or may not have been passed through a call to
   *     <tt>this->initializeOp()</tt> or <tt>this->initializePreconditionedOp()</tt>.
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>If <tt>*Op</tt> on input was initialized through a call to <tt>this->initializeOp()</tt>
   *     and if <tt>fwdOp!=NULL</tt> then <tt>(*fwdOp).get()!=NULL</tt>.
   * <li>If <tt>*Op</tt> was uninitialized on input and <tt>fwdOp!=NULL</tt> then <tt>fwdOp->get()==NULL</tt> out output.
   * <li>On output, <tt>*Op</tt> can be considered to be uninitialized and
   *     it is safe to modify the forward operator object <tt>*(*fwdOp)</tt> returned in <tt>fwdOp</tt>.
   *     The default is <tt>fwdOp==NULL</tt> in which case the forward operator will not be returned in <tt>*fwdOp</tt>.
   * </ul>
   *
   * This function should be called before the forward operator passed in to
   * <tt>this->initializeOp()</tt> is modified.  Otherwise, <tt>*this</tt>
   * could be left in an inconsistent state.  However, this is not required.
   */
  virtual void uninitializeOp(
    LinearOpWithSolveBase<RangeScalar,DomainScalar>                       *Op
    ,Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >  *fwdOp       = NULL
    ,Teuchos::RefCountPtr<const LinearOpBase<DomainScalar,RangeScalar> >  *precOp      = NULL
    ,EPreconditionerInputType                                             *precOpType  = NULL
    ,ESupportSolveUse                                                     *supportSolveUse = NULL
    ) const = 0;
  
  //@}

  /** @name Virtual public functions with default implementations */
  //@{

  /** \brief Return if <tt>solve()</tt> supports the argument <tt>conj</tt>.
   *
   * The default implementation returns <tt>true</tt> for real valued scalar types
   * or when <tt>conj==NONCONJ_ELE</tt> for complex valued types.
   */
  virtual bool solveSupportsConj(EConj conj) const;

  /** \brief Return if <tt>solveTranspose()</tt> supports the argument <tt>conj</tt>.
   *
   * The default implementation returns <tt>false</tt>.
   */
  virtual bool solveTransposeSupportsConj(EConj conj) const;

  /** \brief Return if <tt>solve()</tt> supports the solve tolerance type
   * <tt>ESolveTolType</tt>.
   *
   * The default implementation returns <tt>true</tt> for
   * <tt>solveTolType==SOLVE_TOL_REL_RESIDUAL_NORM</tt>.  Therefore, it is
   * assumed by default that the solver implementation will only be able to
   * check for and enforce a tolerance on the residual.
   */
  virtual bool solveSupportsSolveTolType(EConj conj, ESolveTolType solveTolType) const;

  /** \brief Return if <tt>solveTranspose()</tt> supports the solve tolerance
   * type <tt>ESolveTolType</tt>.
   *
   * The default implementation returns <tt>true</tt> for
   * <tt>solveTolType==SOLVE_TOL_REL_RESIDUAL_NORM</tt>.  Therefore, it is
   * assumed by default that the solver implementation will only be able to
   * check for and enforce a tolerance on the residual.
   */
  virtual bool solveTransposeSupportsSolveTolType(EConj conj, ESolveTolType solveTolType) const;

  /** \brief Initialize a pre-created <tt>LinearOpWithSolveBase</tt> object
   * given a "compatible" <tt>LinearOpBase</tt> object but allow for reuse of
   * any preprocessing that is in <tt>*Op</tt>..
   *
   * \param  fwdOp  [in] The forward linear operator that will be used to create
   *                the output <tt>LinearOpWithSolveBase</tt> object.
   * \param  Op     [in/out] The output <tt>LinearOpWithSolveBase</tt> object.  This object must have
   *                be created first by <tt>this->createOp()</tt> and may have already been through
   *                at least one previous set of calls to <tt>this->initializeOp()</tt> and
   *                <tt>this->uninitializeOp()</tt>.  Note that subclasses
   *                should always first strip off the transpose and scaling by calling <tt>unwrap()</tt>
   *                before attempting to dynamic cast the object.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>fwdOp.get()!=NULL</tt>
   * <li><tt>this->isCompatible(*fwdOp)==true</tt>
   * <li><tt>Op!=NULL</tt>
   * <li><tt>*Op</tt> must have been created by <tt>this->createOp()</tt> prior to calling
   *     this function.
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>Throws <tt>CatastrophicSolveFailure</tt> if the underlying linear solver could
   *     not be created successfully (do to a factorization failure or some other cause).
   * <li><tt>Op->range()->isCompatible(*fwdOp->range())==true</tt>
   * <li><tt>Op->domain()->isCompatible(*fwdOp->domain())==true</tt>
   * <li><tt>Op->apply()</tt> and <tt>Op->applyTranspose()</tt> must behave
   *     exactly the same as <tt>fwdOp->apply()</tt> and <tt>fwdOp->applyTranspose()</tt>
   * <li><t>Op->solveSupportsConj(conj)==this->solveSupportsConj(conj)</tt>
   * <li><t>Op->solveTransposeSupportsConj(conj)==this->solveTransposeSupportsConj(conj)</tt>
   * <li><t>Op->solveSupportsSolveTolType(conj,solveTolType)==this->solveSupportsSolveTolType(conj,solveTolType)</tt>
   * <li><t>Op->solveTransposeSupportsSolveTolType(conj,solveTolType)==this->solveTransposeSupportsSolveTolType(conj,solveTolType)</tt>
   * <li><tt>fwdOp.count()</tt> after output is greater than <tt>fwdOp.count()</tt>
   *     just before this call and therefore the client can assume that the <tt>*fwdOp</tt> object will 
   *     be remembered by the <tt>*Op</tt> object.  The client must be careful
   *     not to modify the <tt>*fwdOp</tt> object or else the <tt>*Op</tt> object may also
   *     be modified.
   * </ul>
   *
   * The purpose of this function is to allow the reuse of old factorizations
   * and/or preconditioners that may go into the initialization of the
   * <tt>*Op</tt> objects.  Note that by calling this function, the
   * performance <tt>Op->solve(...)</tt> may not be as good as when calling
   * the function <tt>this->initializeOp(...,Op)</tt> to initialize
   * <tt>*Op</tt>.
   *
   * The default implementation of this function just calls
   * <tt>this->initializeOp(fwdOp,Op)</tt> which does the default
   * implementation.
   */
  virtual void initializeAndReuseOp(
    const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &fwdOp
    ,LinearOpWithSolveBase<RangeScalar,DomainScalar>                             *Op
    ) const;

  /** \brief Determines if <tt>*this</tt> supports given preconditioner type.
   *
   * The default implementation returns <tt>false</tt>.
   */
  virtual bool supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const;

  /** \brief Initialize a pre-created <tt>LinearOpWithSolveBase</tt> object
   * given a "compatible" <tt>LinearOpBase</tt> object and an optional
   * preconditioner.
   *
   * \param  fwdOp  [in] The forward linear operator that will be used to create
   *                the output <tt>LinearOpWithSolveBase</tt> object.
   * \param  precOp [in] The preconditioner linear operator that will be used to create
   *                the output <tt>LinearOpWithSolveBase</tt> object if preconditioners are supported.
   * \param  precOpType
   *                [in] The type of the preconditioner being passed in.
   * \param  Op     [in/out] The output <tt>LinearOpWithSolveBase</tt> object.  This object must have
   *                be created first by <tt>this->createOp()</tt>.  The object may have also
   *                already been passed through this function several times.  Note that subclasses
   *                should always first strip off the transpose and scaling by calling <tt>unwrap()</tt>
   *                before attempting to dynamic cast the object.
   * \param  supportSolveUse
   *                [in] Determines if <tt>Op->solve(...)</tt> or <tt>Op->solveTranspose(...)</tt> will
   *                be called.  This allows <tt>*this</tt> factory object determine how to best initialae
   *                the <tt>*Op</tt> object.  Default <tt>supportSolveUse=SUPPORT_SOLVE_UNSPECIFIED</tt>
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>fwdOp.get()!=NULL</tt>
   * <li><tt>precOp.get()!=NULL</tt>
   * <li><tt>this->isCompatible(*fwdOp)==true</tt>
   * <li><tt>Op!=NULL</tt>
   * <li><tt>*Op</tt> must have been created by <tt>this->createOp()</tt> prior to calling
   *     this function.
   * <li>It is allowed for an implementation to throw an exception if 
   *     <tt>this->supportsPreconditionerInputType(precOpType)==false</tt> but this is
   *     not required.
   * <li>[<tt>supportSolveUse==SUPPORT_SOLVE_FORWARD_ONLY<tt>]
   *     <tt>this->solveSupportsConj(conj)==true</tt> for any value of <tt>conj</tt>
   * <li>[<tt>supportSolveUse==SUPPORT_SOLVE_TRANSPOSE_ONLY<tt>]
   *     <tt>this->solveTransposeSupportsConj(conj)==true</tt> for any value of <tt>conj</tt>
   * <li>[<tt>supportSolveUse==SUPPORT_SOLVE_FORWARD_AND_TRANSPOSE<tt>]
   *     <tt>this->solveSupportsConj(conj)==true && this->solveTransposeSupportsConj(conj)==true</tt>
   *     for any value of <tt>conj</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>Throws <tt>CatastrophicSolveFailure</tt> if the underlying linear solver could
   *     not be created successfully (do to a factorization failure or some other cause).
   * <li><tt>Op->range()->isCompatible(*fwdOp->range())==true</tt>
   * <li><tt>Op->domain()->isCompatible(*fwdOp->domain())==true</tt>
   * <li><tt>Op->apply()</tt> and <tt>Op->applyTranspose()</tt> must behave
   *     exactly the same as <tt>fwdOp->apply()</tt> and <tt>fwdOp->applyTranspose()</tt>
   * <li><t>Op->solveSupportsConj(conj)==this->solveSupportsConj(conj)</tt>
   * <li><t>Op->solveTransposeSupportsConj(conj)==this->solveTransposeSupportsConj(conj)</tt>
   * <li><t>Op->solveSupportsSolveTolType(conj,solveTolType)==this->solveSupportsSolveTolType(conj,solveTolType)</tt>
   * <li><t>Op->solveTransposeSupportsSolveTolType(conj,solveTolType)==this->solveTransposeSupportsSolveTolType(conj,solveTolType)</tt>
   * <li><tt>fwdOp.count()</tt> after output is greater than <tt>fwdOp.count()</tt>
   *     just before this call and therefore the client can assume that the <tt>*fwdOp</tt> object will 
   *     be remembered by the <tt>*Op</tt> object.  The client must be careful
   *     not to modify the <tt>*fwdOp</tt> object or else the <tt>*Op</tt> object may also
   *     be modified.
   * <li>If <tt>this->supportsPreconditionerInputType(precOpType)==true</tt> then
   *     <ul>
   *       <li><tt>precOp.count()</tt> after output is greater than <tt>precOp.count()</tt>
   *       just before this call and therefore the client can assume that the <tt>*precOp</tt> object will 
   *       be remembered by the <tt>*Op</tt> object.  The client must be careful
   *       not to modify the <tt>*precOp</tt> object or else the <tt>*Op</tt> object may also
   *       be modified.
   *     </ul>
   * <li>else if an exception is not thrown then
   *     <ul>
   *       <li><tt>precOp.count()</tt> after output is equal to <tt>precOp.count()</tt>
   *       just before this call and therefore the <tt>*precOp</tt> is ignored and is not remembered.
   *     </ul>
   * </ul>
   *
   * <b>Warning!</b> It is allowed for an implementation to throw an exception
   * if <tt>this->supportsPreconditionerInputType(precOpType)==false</tt> so
   * therefore a client should not call this function if preconditioners are
   * not supported!  The mode of silently ignoring preconditioners is
   * acceptable in some cases and is therefore allowed behavior.
   *
   * The default implementation throws an exception which is consistent with
   * the default implementation of
   * <tt>this->supportsPreconditionerInputType()</tt> which assumes by default
   * that preconditioners can not be supported..
   */
  virtual void initializePreconditionedOp(
    const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >     &fwdOp
    ,const Teuchos::RefCountPtr<const LinearOpBase<DomainScalar,RangeScalar> >    &precOp
    ,const EPreconditionerInputType                                               precOpType
    ,LinearOpWithSolveBase<RangeScalar,DomainScalar>                              *Op
    ,const ESupportSolveUse                                                       supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
    ) const;

  //@}

};

//@}

/** \defgroup Thyra_LinearOpWithSolveFactoryBase_helper_grp Non-member LinearOpWithSolveFactoryBase helper functions.
 *
 * These functions provide for simpler use cases for the use of
 * <tt>LinearOpWithSolveFactoryBase</tt> objects. and provide some
 * documentation for the various use cases.
 *
 * \ingroup Thyra_Op_Vec_Interoperability_Extended_Interfaces_grp
 */

/** \brief Create and initialize a <tt>LinearOpWithSolveBase</tt> object from
 * a <tt>LinearOpBase</tt> object using a
 * <tt>LinearOpWithSolveFactoryBase</tt> strategy object.
 *
 * \ingroup Thyra_LinearOpWithSolveFactoryBase_helper_grp
 */
template<class RangeScalar, class DomainScalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<RangeScalar,DomainScalar> >
createAndInitializeLinearOpWithSolve(
  const LinearOpWithSolveFactoryBase<RangeScalar,DomainScalar>                  &factory
  ,const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &fwdOp
  ,const ESupportSolveUse                                                       supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  Teuchos::RefCountPtr<LinearOpWithSolveBase<RangeScalar,DomainScalar> >
    Op = factory.createOp();
  factory.initializeOp(fwdOp,&*Op,supportSolveUse);
  return Op;
}

// /////////////////////////
// Implementations

template <class RangeScalar, class DomainScalar>
bool LinearOpWithSolveFactoryBase<RangeScalar,DomainScalar>::solveSupportsConj(EConj conj) const
{
  return true;
}

template <class RangeScalar, class DomainScalar>
bool LinearOpWithSolveFactoryBase<RangeScalar,DomainScalar>::solveTransposeSupportsConj(EConj conj) const
{
  return false;
}

template <class RangeScalar, class DomainScalar>
bool LinearOpWithSolveFactoryBase<RangeScalar,DomainScalar>::solveSupportsSolveTolType(EConj conj, ESolveTolType solveTolType) const
{
  return ( solveTolType == SOLVE_TOL_REL_RESIDUAL_NORM );
}

template <class RangeScalar, class DomainScalar>
bool LinearOpWithSolveFactoryBase<RangeScalar,DomainScalar>::solveTransposeSupportsSolveTolType(EConj conj, ESolveTolType solveTolType) const
{
  return ( solveTolType == SOLVE_TOL_REL_RESIDUAL_NORM );
}

template<class RangeScalar, class DomainScalar>
void LinearOpWithSolveFactoryBase<RangeScalar,DomainScalar>::initializeAndReuseOp(
  const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &fwdOp
  ,LinearOpWithSolveBase<RangeScalar,DomainScalar>                             *Op
  ) const
{
  this->initializeOp(fwdOp,Op);
}

template<class RangeScalar, class DomainScalar>
bool LinearOpWithSolveFactoryBase<RangeScalar,DomainScalar>::supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const
{
  return false;
}

template<class RangeScalar, class DomainScalar>
void LinearOpWithSolveFactoryBase<RangeScalar,DomainScalar>::initializePreconditionedOp(
  const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >     &fwdOp
  ,const Teuchos::RefCountPtr<const LinearOpBase<DomainScalar,RangeScalar> >    &precOp
  ,const EPreconditionerInputType                                               precOpType
  ,LinearOpWithSolveBase<RangeScalar,DomainScalar>                              *Op
  ,const ESupportSolveUse                                                       supportSolveUse
  ) const
{
  TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"Error, the concrete implementation described as \'"<<this->description()<<"\' did not override the "
    "initializePreconditionedOp(...) function and the default implementation throws this exception!"
    );
}

} // namespace Thyra

#endif // THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DECL_HPP
