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

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObject.hpp"


namespace Thyra {


/** \brief Factory interface for creating <tt>LinearOpWithSolveBase</tt>
 * objects from compatible <tt>LinearOpBase</tt> objects.
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 *
 * \section LOWSFB_outline_sec Outline
 *
 * <ul>
 * <li>\ref LOWSFB_intro_sec
 * <li>\ref LOWSFB_use_cases_sec
 * <li>\ref LOWSFB_verbosity_level_sec
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
 * <tt>initializePreconditionedOp()</tt>
 *
 * \section LOWSFB_use_cases_sec Use cases
 *
 * The following use cases demonstrate and specify the behavior of this
 * interface in several different situations.  These use cases don't cover
 * very possible variation but implementors and clients should get a pretty
 * good idea of the desired behavior from what is specified below.
 *
 * <b>Note:</b> All of the code fragments shown in the below use cases is code
 * that is compiled and run automatically by the test harness and the code
 * fragments are automatically pulled into this HTML documentation whenever
 * the documentation is rebuilt.  Therefore, one can have a very high degree
 * of confidence that the code shown in these examples is correct.
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
 * \subsection LOWSFB_single_linear_solve_sec Performing a single linear solve given a forward operator
 *
 * Performing a single linear solve is at minimum a two step process.  The
 * following example function shows how to take a <tt>LinearOpBase</tt>
 * object, a compatible <tt>LinearOpWithSolveFactoryBase</tt> object, a RHS
 * and a LHS and then solve the linear system using the default solve
 * criteria:
 *
 * \dontinclude Thyra_LinearOpWithSolveFactoryExamples.hpp
 * \skip begin singleLinearSolve
 * \skip template
 * \until end singleLinearSolve
 *
 * See the documentation for the <tt>LinearOpWithSolveBase</tt> interface for
 * how to specify solve criteria, how to perform adjoint (i.e. transpose)
 * solve, and how to solve systems with multiple RHSs.
 *
 * Note that the forward operator <tt>A</tt> is passed into this function
 * <tt>singleLinearSolve(...)</tt> as a raw object reference and not as a
 * <tt>Teuchos::RCP</tt> wrapped object since no persisting
 * relationship with this object will be last after this function exists, even
 * if an exception is thrown.
 *
 * Also note that once the above function exits that all memory of the
 * invertible operator created as part of the <tt>invertibleA</tt> object will
 * be erased.  The following use cases show more sophisticated uses of these
 * interfaces.
 *
 * \subsection LOWSFB_scaled_adjoint_sec Creating invertible operators for scaled and/or adjoint forward operators
 *
 * This interface requires that all good implementations support implicitly
 * scaled and adjoint (or transposed) forward operators.  The following
 * example function shows how this looks:
 *
 * \dontinclude Thyra_LinearOpWithSolveFactoryExamples.hpp
 * \skip begin createScaledAdjointLinearOpWithSolve
 * \skip template
 * \until end createScaledAdjointLinearOpWithSolve
 *
 * In the above example, the functions <tt>adjoint()</tt> and <tt>scale()</tt>
 * create an implicitly scaled adjoint operator of type
 * <tt>DefaultScaledAdjointLinearOp</tt> which is then unwrapped by the
 * <tt>lowsFactory</tt> implementation.  The idea is that any operation that
 * works with a particular forward operator <tt>A</tt> should automatically
 * work with an implicitly scaled and/or adjoint view of that forward
 * operator.  The specification of this interface actually says that all
 * <tt>LinearOpBase</tt> objects that support the abstract mix-in interface
 * <tt>ScaledAdjointLinearOpBase</tt> will allow this feature.
 *
 * Above also note that the forward operator <tt>A</tt> is passed in as a
 * <tt>Teuchos::RCP</tt> wrapped object since it will be used to
 * create a persisting relationship with the returned
 * <tt>Thyra::LinearOpWithSolveBase</tt> object.  Also note that the
 * <tt>lowsFactory</tt> object is still passed in as a raw object reference
 * since no persisting relationship with <tt>lowsFactory</tt> is created as a
 * side effect of calling this function.  Remember, the specification of this
 * interface requires that the returned <tt>LinearOpWithSolveBase</tt> objects
 * be independent from the <tt>lowsFactory</tt> object that creates them.  In
 * other words, once a <tt>LinearOpWithSolveFactoryBase</tt> object creates a
 * <tt>LinearOpWithSolveBase</tt> object, then the
 * <tt>LinearOpWithSolveFactoryBase</tt> object can be deleted and the
 * <tt>LinearOpWithSolveBase</tt> it created must still be valid.
 *
 * \subsection LOWSFB_updated_of_linear_op_sec Updates of linear operator between linear solves
 *
 * In this use case is comprised of:
 * <ul>
 * <li> Creating a <tt>LinearOpWithSolveBase</tt> object from a <tt>LinearOpBase</tt> object
 * <li> Using the <tt>LinearOpWithSolveBase</tt> object to solve one or more linear systems
 * <li> Modifying the <tt>LinearOpBase</tt> object and reinitializing the <tt>LinearOpWithSolveBase</tt> object
 * <li> Using the updated <tt>LinearOpWithSolveBase</tt> object to solve one or more linear systems
 * </ul>
 *
 * The below code fragment shows how this looks.
 *
 * \dontinclude Thyra_LinearOpWithSolveFactoryExamples.hpp
 * \skip begin solveNumericalChangeSolve
 * \skip template
 * \until end solveNumericalChangeSolve
 *
 * In the above code fragment the call to the function
 * <tt>lowsFactory.uninitializeOp(&*invertibleA)</tt> may not fully
 * uninitialize the <tt>*invertibleA</tt> object as it may contain memory of a
 * factorization structure, or a communication pattern, or whatever may have
 * been computed in the first call to
 * <tt>lowsFactory.initializeOp(rcpA,&*invertibleA)</tt>.  This allows for a
 * more efficient reinitialization on the second call of
 * <tt>lowsFactory.initializeOp(rcpA,&*invertibleA)</tt>.
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
 *
 * \dontinclude Thyra_LinearOpWithSolveFactoryExamples.hpp
 * \skip begin solveSmallNumericalChangeSolve
 * \skip template
 * \until end solveSmallNumericalChangeSolve
 *
 * Note that the <tt>*invertibleA</tt> object reinitialized in the second
 * call to <tt>lowsFactory.initializeAndReuseOp(rcpA,&*invertibleA)</tt> must
 * be able to correctly solve the linear systems to an appropriate accuracy.
 * For example, a preconditioned iterative linear solver might keep the same
 * preconditioner but would use the updated forward operator to define the
 * linear system.  A direct solver can not reuse the factorization unless it
 * has an iterative refinement feature (which uses the updated forward
 * operator) or some other device.  Or a direct solver might, for instance,
 * reuse the same factorization structure and not re-pivot where it might
 * re-pivot generally.
 *
 * \subsection LOWSFB_major_changes_in_op_sec Updating the linear solver for major changes in the structure of the forward operator
 *
 * A major change in the shape or structure of a forward operator generally
 * eliminates the possibility of reusing any computations between different
 * calls to <tt>initializeOp()</tt>.  In these instances, the client might as
 * well just recreate the <tt>LinearOpWithSolveBase</tt> using <tt>createOp()</tt>
 * before the reinitialization.
 *
 * The below example function shows what this looks like:
 *
 * \dontinclude Thyra_LinearOpWithSolveFactoryExamples.hpp
 * \skip begin solveMajorChangeSolve
 * \skip template
 * \until end solveMajorChangeSolve
 *
 * \subsection LOWSFB_external_prec_sec Externally defined preconditioners
 *
 * This interface also supports the use of externally defined preconditioners
 * that are created and controlled by the client.  Client-created and
 * client-controlled preconditioners can be passed along with the forward
 * operator through the function initializePreconditionedOp().  Only
 * <tt>LinearOpWithSolveFactoryBase</tt> objects that return
 * <tt>supportsPreconditionerType()==true</tt> will support externally defined
 * preconditioners.  In general, iterative linear solver implementations will
 * support externally defined preconditioners and direct linear solver
 * implementations will not.
 *
 * Externally defined preconditioners can be passed in as operators or as
 * matrices and these two sub use cases are described below.
 *
 * \subsubsection LOWSFB_external_prec_op_sec Use of externally defined preconditioner operators
 *
 * The client can use externally defined preconditioner linear operator(s) for
 * <tt>LinearOpWithSolveFactoryBase</tt> objects that accept preconditioners.
 * Preconditioner operator(s) are specified through the
 * <tt>PreconditionerBase</tt> interface.  A <tt>PreconditionerBase</tt>
 * object is essentially an encapsulation of a left, right, or split left/right
 * preconditioner.
 *
 * <b>Using an externally defined <tt>PreconditionerBase</tt> object</b>
 *
 * Given an externally defined <tt>PreconditionerBase</tt> object, the client
 * can pass it allow with a forward linear operator to initialize a
 * <tt>LinearOpWithSolveBase</tt> object as shown in the following code
 * fragment:
 *
 * \dontinclude Thyra_LinearOpWithSolveFactoryExamples.hpp
 * \skip begin createGeneralPreconditionedLinearOpWithSolve
 * \skip template
 * \until end createGeneralPreconditionedLinearOpWithSolve
 *
 * <b>Using a single externally defined <tt>LinearOpBase</tt> object as an unspecified preconditioner</b>
 *
 * A client can easily use any compatible appropriately defined
 * <tt>LinearOpBase</tt> object as a preconditioner operator.  To do so, it can
 * just be wrapped by the default <tt>PreconditionerBase</tt> implementation
 * class <tt>DefaultPreconditioner</tt>.
 *
 * The following code fragment shows show to use an externally defined
 * <tt>LinearOpBase</tt> object as a preconditioner operator without
 * specifying whether it should be applied on the right or on the left.
 *
 * \dontinclude Thyra_LinearOpWithSolveFactoryExamples.hpp
 * \skip begin createUnspecifiedPreconditionedLinearOpWithSolve
 * \skip template
 * \until end createUnspecifiedPreconditionedLinearOpWithSolve
 *
 * <b>Using a single externally defined <tt>LinearOpBase</tt> object as a left preconditioner</b>
 *
 * An externally defined <tt>LinearOpBase</tt> object can be used as a
 * preconditioner operator targeted as a left preconditioner as shown in the
 * following code fragment:
 *
 * \dontinclude Thyra_LinearOpWithSolveFactoryExamples.hpp
 * \skip begin createLeftPreconditionedLinearOpWithSolve
 * \skip template
 * \until end createLeftPreconditionedLinearOpWithSolve
 *
 * <b>Using a single externally defined <tt>LinearOpBase</tt> object as a right preconditioner</b>
 *
 * An externally defined <tt>LinearOpBase</tt> object can be used as a
 * preconditioner operator targeted as a right preconditioner as shown in the
 * following code fragement:
 *
 * \dontinclude Thyra_LinearOpWithSolveFactoryExamples.hpp
 * \skip begin createRightPreconditionedLinearOpWithSolve
 * \skip template
 * \until end createRightPreconditionedLinearOpWithSolve
 *
 * <b>Using two single externally defined <tt>LinearOpBase</tt> objects as a split left/right preconditioner</b>
 *
 * Two externally defined <tt>LinearOpBase</tt> objects can be used as
 * preconditioner operators for form a split preconditioner targeted to the
 * left and right as shown in the following code fragement:
 *
 * \dontinclude Thyra_LinearOpWithSolveFactoryExamples.hpp
 * \skip begin createLeftRightPreconditionedLinearOpWithSolve
 * \skip template
 * \until end createLeftRightPreconditionedLinearOpWithSolve
 *
 * \subsubsection LOWSFB_external_prec_mat_sec The use of externally defined preconditioner matrices
 *
 * A different linear operator can be used to build a preconditioner
 * internally when initializing a <tt>LinearOpWithSolveBase</tt> object.  The
 * following code fragement shows how to use an approximate forward linear
 * from which to form a preconditioner internally:
 *
 * \dontinclude Thyra_LinearOpWithSolveFactoryExamples.hpp
 * \skip begin createMatrixPreconditionedLinearOpWithSolve
 * \skip template
 * \until end createMatrixPreconditionedLinearOpWithSolve
 *
 * \subsubsection LOWSFB_external_prec_op_reuse_sec Reuse of externally defined preconditioner operators
 *
 * Reusing an externally defined preconditioner is a very straightforward
 * matter conceptually.  Since the client controls the creation and
 * initialization of the preconditioner, the client can control which
 * preconditioners are paired with which forward operators in order to form
 * <tt>LinearOpWithSolveBase</tt> objects.  However, when an abstract
 * <tt>PreconditionerFactoryBase</tt> object is used to create and maintain
 * the preconditioner, the exact nature of the preconditioner after the
 * operator is updated is undefined.  The following code fragment show a use
 * case where the operator is reused without explicitly updating the
 * preconditioner between solves:
 *
 * \dontinclude Thyra_LinearOpWithSolveFactoryExamples.hpp
 * \skip begin externalPreconditionerReuseWithSolves
 * \skip template
 * \until end externalPreconditionerReuseWithSolves
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
 * \section LOWSFB_verbosity_level_sec Output and verbosity
 *
 * <b>TODO:</b> Provide some guidance on how clients and subclasses should
 * interpret Teuchos::EVerbLevel.
 * 
 * \section LOWSFB_developer_notes_sec Notes to subclass developers
 *
 * This interface assumes a minimal default set of functionality that is
 * appropriate for direct and simple iterative linear solver implementations.
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
template<class Scalar>
class LinearOpWithSolveFactoryBase
  : virtual public Teuchos::Describable, 
    virtual public Teuchos::VerboseObject<LinearOpWithSolveFactoryBase<Scalar> >,
    virtual public Teuchos::ParameterListAcceptor
{
public:

  /** @name Preconditioner Factory Management */
  //@{

  /** \brief Determines if <tt>*this</tt> accepts external preconditioner factories.
   *
   * The default implementation returns <tt>false</tt>.
   */
  virtual bool acceptsPreconditionerFactory() const;

  /** \brief Set a preconditioner factory object.
   *
   * \param precFactory [in] The preconditioner factory to be used internally
   * to create preconditioners.
   *
   * \param precFactoryName [in] The name to give to the preconditioner
   * factory internally.  This name is used when setting parameters in the
   * parameter list.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->acceptsPreconditionerFactory()==true</tt>
   * <li><tt>precFactory.get()!=NULL</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->getPreconditionerFactory().get()==precFactory.get()</tt>
   * </ul>
   *
   * The default implementation thrown an exception which is consistent with
   * the default implementation <tt>acceptsPreconditionerFactory()==false</tt>.
   */
  virtual void setPreconditionerFactory(
    const RCP<PreconditionerFactoryBase<Scalar> >  &precFactory,
    const std::string &precFactoryName
    );

  /** \brief Get a preconditioner factory object.
   *
   * The default implementation returns <tt>Teuchos::null</tt>.
   */
  virtual RCP<PreconditionerFactoryBase<Scalar> >
  getPreconditionerFactory() const;

  /** \brief Unset the preconditioner factory (if one is set).
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->getPreconditionerFactory().get()==NULL</tt>
   * </ul>
   *
   * The default implementation returns <tt>Teuchos::null</tt>.
   */
  virtual void unsetPreconditionerFactory(
    RCP<PreconditionerFactoryBase<Scalar> >  *precFactory = NULL,
    std::string *precFactoryName = NULL
    );

  //@}

  /** @name Creation/Initialization of basic LinearOpWithSolveBase objects */
  //@{

  /** \brief Check that a <tt>LinearOpBase</tt> object is compatible with
   * <tt>*this</tt> factory object.
   */
  virtual bool isCompatible(
    const LinearOpSourceBase<Scalar> &fwdOpSrc ) const = 0;

  /** \brief Create an (uninitialized) <tt>LinearOpWithSolveBase</tt> object
   * to be initialized later in <tt>this->initializeOp()</tt>.
   *
   * Note that on output <tt>return->domain().get()==NULL</tt> may be true
   * which means that the operator is not fully initialized.  In fact, the
   * output operator object is not guaranteed to be fully initialized until
   * after it is passed through <tt>this->initializeOp()</tt>.
   */
  virtual RCP<LinearOpWithSolveBase<Scalar> >
  createOp() const = 0;

  /** \brief Initialize a pre-created <tt>LinearOpWithSolveBase</tt> object
   * given a "compatible" <tt>LinearOpBase</tt> object.
   *
   * \param fwdOpSrc [in] The forward linear operator that will be used to
   * create the output <tt>LinearOpWithSolveBase</tt> object.  Note that this
   * object is remembered by the <tt>*Op</tt> object on output.
   *
   * \param Op [in/out] The output <tt>LinearOpWithSolveBase</tt> object.
   * This object must have be created first by <tt>this->createOp()</tt>.  The
   * object may have also already been passed through this function several
   * times.  Note that subclasses should always first strip off the transpose
   * and scaling by calling <tt>unwrap()</tt> before attempting to dynamic
   * cast the object.
   *
   * \param supportSolveUse [in] Determines if <tt>Op->solve(...)</tt> or
   * <tt>Op->solveTranspose(...)</tt> will be called.  This allows
   * <tt>*this</tt> factory object to determine how to best initialize the
   * <tt>*Op</tt> object.  Default
   * <tt>supportSolveUse=SUPPORT_SOLVE_UNSPECIFIED</tt>
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>fwdOpSrc.get()!=NULL</tt>
   *
   * <li><tt>this->isCompatible(*fwdOpSrc)==true</tt>
   *
   * <li><tt>Op!=NULL</tt>
   *
   * <li><tt>*Op</tt> must have been created by <tt>this->createOp()</tt>
   * prior to calling this function.
   *
   * <li>[<tt>supportSolveUse==SUPPORT_SOLVE_FORWARD_ONLY<tt>]
   * <tt>this->solveSupportsConj(conj)==true</tt> for any value of
   * <tt>conj</tt>
   *
   * <li>[<tt>supportSolveUse==SUPPORT_SOLVE_TRANSPOSE_ONLY<tt>]
   * <tt>this->solveTransposeSupportsConj(conj)==true</tt> for any value of
   * <tt>conj</tt>
   *
   * <li>[<tt>supportSolveUse==SUPPORT_SOLVE_FORWARD_AND_TRANSPOSE<tt>]
   * <tt>this->solveSupportsConj(conj)==true &&
   * this->solveTransposeSupportsConj(conj)==true</tt> for any value of
   * <tt>conj</tt>
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li>Throws <tt>CatastrophicSolveFailure</tt> if the underlying linear
   * solver could not be created successfully (e.g. due to a factorization
   * failure or some other cause).
   *
   * <li><tt>Op->range()->isCompatible(*fwdOpSrc->range())==true</tt>
   *
   * <li><tt>Op->domain()->isCompatible(*fwdOpSrc->domain())==true</tt>
   *
   * <li><tt>Op->apply()</tt> and <tt>Op->applyTranspose()</tt> must behave
   * exactly the same as <tt>fwdOpSrc->apply()</tt> and
   * <tt>fwdOpSrc->applyTranspose()</tt>
   *
   * <li><tt>Op->solveSupportsConj(conj)==this->solveSupportsConj(conj)</tt>
   *
   * <li><tt>Op->solveTransposeSupportsConj(conj)==this->solveTransposeSupportsConj(conj)</tt>
   *
   * <li><tt>fwdOpSrc.count()</tt> after output is greater than
   * <tt>fwdOpSrc.count()</tt> just before this call and therefore the client
   * can assume that the <tt>*fwdOpSrc</tt> object will be remembered by the
   * <tt>*Op</tt> object.  The client must be careful not to modify the
   * <tt>*fwdOpSrc</tt> object or else the <tt>*Op</tt> object may also be
   * modified and become invalid.
   *
   * </ul>
   */
  virtual void initializeOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
    ) const = 0;

  /** \brief Initialize a pre-created <tt>LinearOpWithSolveBase</tt> object
   * given a "compatible" <tt>LinearOpBase</tt> object but allow for reuse of
   * any preprocessing that is in <tt>*Op</tt>.
   *
   * \param fwdOpSrc [in] The forward linear operator that will be used to
   * create the output <tt>LinearOpWithSolveBase</tt> object.
   *
   * \param Op [in/out] The output <tt>LinearOpWithSolveBase</tt> object.
   * This object must have be created first by <tt>this->createOp()</tt> and
   * may have already been through at least one previous set of calls to
   * <tt>this->initializeOp()</tt> and <tt>this->uninitializeOp()</tt>.  Note
   * that subclasses should always first strip off the transpose and scaling
   * by calling <tt>unwrap()</tt> before attempting to dynamic cast the
   * object.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>fwdOpSrc.get()!=NULL</tt>
   *
   * <li><tt>this->isCompatible(*fwdOpSrc)==true</tt>
   *
   * <li><tt>Op!=NULL</tt>
   *
   * <li><tt>*Op</tt> must have been created by <tt>this->createOp()</tt>
   * prior to calling this function.
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li>Throws <tt>CatastrophicSolveFailure</tt> if the underlying linear
   * solver could not be created successfully (e.g. due to a factorization
   * failure or some other cause).
   *
   * <li><tt>Op->range()->isCompatible(*fwdOpSrc->range())==true</tt>
   *
   * <li><tt>Op->domain()->isCompatible(*fwdOpSrc->domain())==true</tt>
   *
   * <li><tt>Op->apply()</tt> and <tt>Op->applyTranspose()</tt> must behave
   * exactly the same as <tt>fwdOpSrc->apply()</tt> and
   * <tt>fwdOpSrc->applyTranspose()</tt>
   *
   * <li><tt>Op->solveSupportsConj(conj)==this->solveSupportsConj(conj)</tt>
   *
   * <li><tt>Op->solveTransposeSupportsConj(conj)==this->solveTransposeSupportsConj(conj)</tt>
   *
   * <li><tt>fwdOpSrc.count()</tt> after output is greater than
   * <tt>fwdOpSrc.count()</tt> just before this call and therefore the client
   * can assume that the <tt>*fwdOpSrc</tt> object will be remembered by the
   * <tt>*Op</tt> object.  The client must be careful not to modify the
   * <tt>*fwdOpSrc</tt> object or else the <tt>*Op</tt> object may also be
   * modified.
   *
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
   * <tt>this->initializeOp(fwdOpSrc,Op)</tt> which does the default
   * implementation.
   */
  virtual void initializeAndReuseOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op
    ) const;

  /** \brief Uninitialize a <tt>LinearOpWithSolveBase</tt> object and return
   * its remembered forward linear operator and potentially also its
   * externally generated preconditioner.
   *
   * \param Op [in/out] On input, <tt>*Op</tt> is an initialized or
   * uninitialized object and on output is uninitialized.  Note that
   * "uninitialized" does not mean that <tt>Op</tt> is completely stateless.
   * It may still remember some aspect of the matrix <tt>fwdOpSrc</tt> that
   * will allow for a more efficient initialization next time through
   * <tt>this->initializeOp()</tt>.
   *
   * \param fwdOpSrc [in/out] If <tt>fwdOpSrc!=NULL</tt> on input, then on
   * output this is set to the same forward operator passed into
   * <tt>this->initializeOp()</tt>.
   *
   * \param prec [in/out] If <tt>prep!=NULL</tt> on input, then on output,
   * this this is set to same preconditioner that was passed into
   * <tt>this->initializePreconditionedOp()</tt>.
   *
   * \param approxFwdOpSrc [in/out] If <tt>approxFwdOpSrc!=NULL</tt> on input,
   * then on output, this is set to same approximate forward operator that was
   * passed into <tt>this->initializePreconditionedOp()</tt>.
   *
   * \param ESupportSolveUse [in/out] If <tt>fwdOpSrc!=NULL</tt> on input,
   * then on output this is set to same option value passed to
   * <tt>this->initializeOp()</tt>.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>*Op</tt> must have been created by <tt>this->createOp()</tt>
   * prior to calling this function.
   *
   * <li><tt>Op</tt> may or may not have been passed through a call to
   * <tt>this->initializeOp()</tt> or
   * <tt>this->initializePreconditionedOp()</tt>.
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li>If <tt>*Op</tt> on input was initialized through a call to
   * <tt>this->initializeOp()</tt> and if <tt>fwdOpSrc!=NULL</tt> then
   * <tt>(*fwdOpSrc).get()!=NULL</tt>.
   *
   * <li>If <tt>*Op</tt> was uninitialized on input and
   * <tt>fwdOpSrc!=NULL</tt> then <tt>fwdOpSrc->get()==NULL</tt> out output.
   *
   * <li>On output, <tt>*Op</tt> can be considered to be uninitialized and it
   * is safe to modify the forward operator object <tt>*(*fwdOpSrc)</tt>
   * returned in <tt>fwdOpSrc</tt>.  The default is <tt>fwdOpSrc==NULL</tt> in
   * which case the forward operator will not be returned in
   * <tt>*fwdOpSrc</tt>.
   *
   * </ul>
   *
   * This function should be called before the forward operator passed in to
   * <tt>this->initializeOp()</tt> is modified.  Otherwise, <tt>*this</tt>
   * could be left in an inconsistent state.  However, this is not required.
   */
  virtual void uninitializeOp(
    LinearOpWithSolveBase<Scalar> *Op,
    RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc = NULL,
    RCP<const PreconditionerBase<Scalar> > *prec = NULL,
    RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc = NULL,
    ESupportSolveUse *supportSolveUse = NULL
    ) const = 0;
  
  //@}

  /** @name Creation/Initialization of Preconditioned LinearOpWithSolveBase objects */
  //@{

  /** \brief Determines if <tt>*this</tt> supports given preconditioner type.
   *
   * The default implementation returns <tt>false</tt>.
   */
  virtual bool supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const;

  /** \brief Initialize a pre-created <tt>LinearOpWithSolveBase</tt> object
   * given a "compatible" <tt>LinearOpBase</tt> object and an optional
   * <tt>PreconditionerBase</tt> object.
   *
   * \param fwdOpSrc [in] The forward linear operator that will be used to
   * create the output <tt>LinearOpWithSolveBase</tt> object.
   *
   * \param prec [in] The preconditioner that will be used to create the
   * output <tt>LinearOpWithSolveBase</tt> object if preconditioners are
   * supported.
   *
   * \param Op [in/out] The output <tt>LinearOpWithSolveBase</tt> object.
   * This object must have be created first by <tt>this->createOp()</tt>.  The
   * object may have also already been passed through this function several
   * times.  Note that subclasses should always first strip off the transpose
   * and scaling by calling <tt>unwrap()</tt> before attempting to dynamic
   * cast the object.
   *
   * \param supportSolveUse [in] Determines if <tt>Op->solve(...)</tt> or
   * <tt>Op->solveTranspose(...)</tt> will be called.  This allows
   * <tt>*this</tt> factory object determine how to best initialize the
   * <tt>*Op</tt> object.  Default
   * <tt>supportSolveUse=SUPPORT_SOLVE_UNSPECIFIED</tt>
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>fwdOpSrc.get()!=NULL</tt>
   *
   * <li><tt>prec.get()!=NULL</tt>
   *
   * <li><tt>this->isCompatible(*fwdOpSrc)==true</tt>
   *
   * <li><tt>Op!=NULL</tt>
   *
   * <li><tt>*Op</tt> must have been created by <tt>this->createOp()</tt>
   * prior to calling this function.
   *
   * <li>It is allowed for an implementation to throw an exception if
   * <tt>this->supportsPreconditionerInputType(PRECONDITIONER_INPUT_TYPE_AS_OPERATOR)==false</tt>
   * but this is not required.
   *
   * <li>[<tt>supportSolveUse==SUPPORT_SOLVE_FORWARD_ONLY<tt>]
   * <tt>this->solveSupportsConj(conj)==true</tt> for any value of
   * <tt>conj</tt>
   *
   * <li>[<tt>supportSolveUse==SUPPORT_SOLVE_TRANSPOSE_ONLY<tt>]
   * <tt>this->solveTransposeSupportsConj(conj)==true</tt> for any value of
   * <tt>conj</tt>
   *
   * <li>[<tt>supportSolveUse==SUPPORT_SOLVE_FORWARD_AND_TRANSPOSE<tt>]
   * <tt>this->solveSupportsConj(conj)==true &&
   * this->solveTransposeSupportsConj(conj)==true</tt> for any value of
   * <tt>conj</tt>
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li>Throws <tt>CatastrophicSolveFailure</tt> if the underlying linear solver could
   *     not be created successfully (e.g. due to a factorization failure or some other cause).
   *
   * <li><tt>Op->range()->isCompatible(*fwdOpSrc->range())==true</tt>
   *
   * <li><tt>Op->domain()->isCompatible(*fwdOpSrc->domain())==true</tt>
   *
   * <li><tt>Op->apply()</tt> and <tt>Op->applyTranspose()</tt> must behave
   *     exactly the same as <tt>fwdOpSrc->apply()</tt> and <tt>fwdOpSrc->applyTranspose()</tt>
   *
   * <li><t>Op->solveSupportsConj(conj)==this->solveSupportsConj(conj)</tt>
   *
   * <li><t>Op->solveTransposeSupportsConj(conj)==this->solveTransposeSupportsConj(conj)</tt>
   *
   * <li><tt>fwdOpSrc.count()</tt> after output is greater than <tt>fwdOpSrc.count()</tt>
   *     just before this call and therefore the client can assume that the <tt>*fwdOpSrc</tt> object will 
   *     be remembered by the <tt>*Op</tt> object.  The client must be careful
   *     not to modify the <tt>*fwdOpSrc</tt> object or else the <tt>*Op</tt> object may also
   *     be modified.
   *
   * <li>If
   * <tt>this->supportsPreconditionerInputType(PRECONDITIONER_INPUT_TYPE_AS_OPERATOR)==true</tt>
   * then
   *
   *     <ul>
   *
   *       <li><tt>prec.count()</tt> after output is greater than
   *       <tt>prec.count()</tt> just before this call and therefore the
   *       client can assume that the <tt>*prec</tt> object will be remembered
   *       by the <tt>*Op</tt> object.  The client must be careful not to
   *       modify the <tt>*prec</tt> object or else the <tt>*Op</tt> object
   *       may also be modified.
   *
   *     </ul>
   *
   * <li>else if an exception is not thrown then
   *
   *     <ul>
   *
   *       <li><tt>prec.count()</tt> after output is equal to
   *       <tt>prec.count()</tt> just before this call and therefore the
   *       <tt>*prec</tt> is ignored and is not remembered.
   *
   *     </ul>
   *
   * </ul>
   *
   * <b>Warning!</b> It is allowed for an implementation to throw an exception
   * if
   * <tt>this->supportsPreconditionerInputType(PRECONDITIONER_INPUT_TYPE_AS_OPERATOR)==false</tt>
   * so therefore a client should not call this function if preconditioners
   * are not supported!  The mode of silently ignoring preconditioners is
   * acceptable in some cases and is therefore allowed behavior.
   *
   * The default implementation throws an exception which is consistent with
   * the default implementation of
   * <tt>this->supportsPreconditionerInputType()</tt> which assumes by default
   * that preconditioners can not be supported..
   */
  virtual void initializePreconditionedOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const RCP<const PreconditionerBase<Scalar> > &prec,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
    ) const;

  /** \brief Initialize a pre-created <tt>LinearOpWithSolveBase</tt> object
   * given a "compatible" forward <tt>LinearOpBase</tt> object and an
   * approximate forward <tt>LinearOpBase</tt> object.
   *
   * \param fwdOpSrc [in] The forward linear operator that will be used to
   * create the output <tt>LinearOpWithSolveBase</tt> object.
   *
   * \param approxFwdOpSrc [in] Approximation to <tt>fwdOpSrc</tt> from which
   * a preconditioner will be created for.
   *
   * \param Op [in/out] The output <tt>LinearOpWithSolveBase</tt> object.
   * This object must have be created first by <tt>this->createOp()</tt>.  The
   * object may have also already been passed through this function several
   * times.  Note that subclasses should always first strip off the transpose
   * and scaling by calling <tt>unwrap()</tt> before attempting to dynamic
   * cast the object.
   *
   * \param supportSolveUse [in] Determines if <tt>Op->solve(...)</tt> or
   * <tt>Op->solveTranspose(...)</tt> will be called.  This allows
   * <tt>*this</tt> factory object determine how to best initialize the
   * <tt>*Op</tt> object.  Default
   * <tt>supportSolveUse=SUPPORT_SOLVE_UNSPECIFIED</tt>
   *
   * ToDo: finish documentation!
   */
  virtual void initializeApproxPreconditionedOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
    ) const;
  
  //@}

private:
  
  // Not defined and not to be called
  LinearOpWithSolveFactoryBase<Scalar>&
  operator=(const LinearOpWithSolveFactoryBase<Scalar>&);

};


} // namespace Thyra


#endif // THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DECL_HPP
