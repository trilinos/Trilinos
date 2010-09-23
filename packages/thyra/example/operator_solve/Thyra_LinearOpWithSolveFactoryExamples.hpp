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

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_EXAMPLES_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_EXAMPLES_HPP


#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorBase.hpp"


namespace Thyra {


//
// Helper code
//

/** \brief Silly abstract strategy interface for changing
 * <tt>Thyra::LinearOpBase</tt> objects.
 *
 * The only purpose of this class is to allow me to write the use cases
 * involving <tt>LinearOpWithSolveFactoryBase</tt> and
 * <tt>PreconditionerFactoryBase</tt> in a general way.  This class is not to
 * be used to do anything really useful!
 *
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
class LinearOpChanger {
public:
  /** \brief . */
  virtual ~LinearOpChanger() {}
  /** \brief . */
  virtual void changeOp( const Teuchos::Ptr<LinearOpBase<Scalar> > &op ) const = 0;
};

/** \brief No-op changer.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
class NullLinearOpChanger : public LinearOpChanger<Scalar> {
public:
  /** \brief . */
  void changeOp( const Teuchos::Ptr<LinearOpBase<Scalar> > &op ) const {}
};


} // namespace Thyra


//
// Individual non-externally preconditioned use cases
//


// begin singleLinearSolve
/** \brief Performing a single linear solve given a forward operator.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
void singleLinearSolve(
  const Thyra::LinearOpBase<Scalar> &A,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const Thyra::VectorBase<Scalar> &b,
  const Teuchos::Ptr<Thyra::VectorBase<Scalar> > &x,
  Teuchos::FancyOStream &out
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::rcpFromRef;
  Teuchos::OSTab tab(out);
  out << "\nPerforming a single linear solve ...\n";
  // Create the LOWSB object that will be used to solve the linear system
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > invertibleA =
    Thyra::linearOpWithSolve(lowsFactory, rcpFromRef(A));
  // Solve the system using a default solve criteria using a non-member helper function 
  assign(x, ST::zero()); // Must initialize to a guess before solve!
  Thyra::SolveStatus<Scalar> 
    status = Thyra::solve<Scalar>(*invertibleA, Thyra::NOTRANS, b, x);
  out << "\nSolve status:\n" << status;
} // end singleLinearSolve


// begin createScaledAdjointLinearOpWithSolve
/** \brief Create a <tt>LinearOpWithSolveBase</tt> object from an implicitly
 * scaled adjoint <tt>LinearOpBase</tt> object.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
createScaledAdjointLinearOpWithSolve(
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > &A,
  const Scalar &scalar,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  Teuchos::FancyOStream &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating a scaled adjoint LinearOpWithSolveBase object ...\n";
  // Create the LOWSB object that will be used to solve the linear system
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > invertibleAdjointA =
    Thyra::linearOpWithSolve(lowsFactory,scale(scalar,adjoint(A)));
  out << "\nCreated LOWSB object:\n" << describe(*invertibleAdjointA,
    Teuchos::VERB_MEDIUM);
  return invertibleAdjointA;
} // end createScaledAdjointLinearOpWithSolve


// begin solveNumericalChangeSolve
/** \brief Create a <tt>LinearOpWithSolveBase</tt> object and perform a solve,
 * then change the <tt>LinearOpBase</tt> object and solve again.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
void solveNumericalChangeSolve(
  const Teuchos::Ptr<Thyra::LinearOpBase<Scalar> > &A,
  const Thyra::LinearOpChanger<Scalar> &opChanger,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const Thyra::VectorBase<Scalar> &b1,
  const Teuchos::Ptr<Thyra::VectorBase<Scalar> > &x1,
  const Thyra::VectorBase<Scalar> &b2,
  const Teuchos::Ptr<Thyra::VectorBase<Scalar> > &x2,
  Teuchos::FancyOStream &out
  )
{
  using Teuchos::as; using Teuchos::ptr; using Teuchos::rcpFromPtr;
  Teuchos::OSTab tab(out);
  out << "\nPerforming a solve, changing the operator, then performing another"
      << " solve ...\n";
  // Get a local non-owned RCP to A to be used by lowsFactory
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > rcpA = rcpFromPtr(A);
  // Create the LOWSB object that will be used to solve the linear system
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > invertibleA =
    lowsFactory.createOp();
  // Initialize the invertible linear operator given the forward operator
  Thyra::initializeOp<Scalar>(lowsFactory, rcpA, invertibleA.ptr());
  // Solve the system using a default solve criteria using a non-member helper function
  Thyra::assign(x1, as<Scalar>(0.0));
  Thyra::solve<Scalar>(*invertibleA, Thyra::NOTRANS, b1, x1);
  // Before the forward operator A is changed it is recommended that you
  // uninitialize *invertibleA first to avoid accidental use of *invertiableA
  // while it may be in an inconsistent state from the time between *A changes
  // and *invertibleA is explicitly updated. However, this step is not
  // required!
  Thyra::uninitializeOp<Scalar>(lowsFactory, invertibleA.ptr());
  // Change the operator and reinitialize the invertible operator
  opChanger.changeOp(A);
  Thyra::initializeOp<Scalar>(lowsFactory, rcpA, invertibleA.ptr());
  // Note that above will reuse any factorization structures that may have been
  // created in the first call to initializeOp(...).
  // Finally, solve another linear system with new values of A
  Thyra::assign<Scalar>(x2, as<Scalar>(0.0));
  Thyra::solve<Scalar>(*invertibleA, Thyra::NOTRANS, b2, x2);
} // end solveNumericalChangeSolve


// begin solveSmallNumericalChangeSolve
/** \brief Create a <tt>LinearOpWithSolveBase</tt> object and perform a solve,
 * then change the <tt>LinearOpBase</tt> object in a very small way and solve
 * again.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
void solveSmallNumericalChangeSolve(
  const Teuchos::Ptr<Thyra::LinearOpBase<Scalar> > &A,
  const Thyra::LinearOpChanger<Scalar> &opSmallChanger,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const Thyra::VectorBase<Scalar> &b1,
  const Teuchos::Ptr<Thyra::VectorBase<Scalar> > &x1,
  const Thyra::VectorBase<Scalar> &b2,
  const Teuchos::Ptr<Thyra::VectorBase<Scalar> > &x2,
  Teuchos::FancyOStream &out
  )
{
  using Teuchos::ptr; using Teuchos::as; using Teuchos::rcpFromPtr;
  Teuchos::OSTab tab(out);
  out << "\nPerforming a solve, changing the operator in a very small way,"
      << " then performing another solve ...\n";
  // Get a local non-owned RCP to A to be used by lowsFactory
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > rcpA = rcpFromPtr(A);
  // Create the LOWSB object that will be used to solve the linear system
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > invertibleA =
    lowsFactory.createOp();
  // Initialize the invertible linear operator given the forward operator
  Thyra::initializeOp<Scalar>(lowsFactory, rcpA, invertibleA.ptr());
  // Solve the system using a default solve criteria using a non-member helper function
  Thyra::assign(x1, as<Scalar>(0.0));
  Thyra::solve<Scalar>(*invertibleA, Thyra::NOTRANS, b1, x1);
  // Before the forward operator A is changed it is recommended that you
  // uninitialize *invertibleA first to avoid accidental use of *invertiableA
  // while it may be in an inconsistent state from the time between *A changes
  // and *invertibleA is explicitly updated. However, this step is not
  // required!
  Thyra::uninitializeOp<Scalar>(lowsFactory, invertibleA.ptr());
  // Change the operator and reinitialize the invertible operator
  opSmallChanger.changeOp(A);
  Thyra::initializeAndReuseOp<Scalar>(lowsFactory, rcpA, invertibleA.ptr());
  // Note that above a maximum amount of reuse will be achieved, such as
  // keeping the same precondtioner.
  Thyra::assign(x2, as<Scalar>(0.0));
  Thyra::solve<Scalar>(*invertibleA, Thyra::NOTRANS, b2, x2);
} // end solveSmallNumericalChangeSolve


// begin solveMajorChangeSolve
/** \brief Create a <tt>LinearOpWithSolveBase</tt> object and perform a solve,
 * then change the <tt>LinearOpBase</tt> object in a very small way and solve
 * again.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
void solveMajorChangeSolve(
  const Teuchos::Ptr<Thyra::LinearOpBase<Scalar> > &A,
  const Thyra::LinearOpChanger<Scalar> &opMajorChanger,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const Thyra::VectorBase<Scalar> &b1,
  const Teuchos::Ptr<Thyra::VectorBase<Scalar> > &x1,
  const Thyra::VectorBase<Scalar> &b2,
  const Teuchos::Ptr<Thyra::VectorBase<Scalar> > &x2,
  Teuchos::FancyOStream &out
  )
{
  using Teuchos::as; using Teuchos::rcpFromPtr;
  Teuchos::OSTab tab(out);
  out << "\nPerforming a solve, changing the operator in a major way, then performing"
      << " another solve ...\n";
  // Get a local non-owned RCP to A to be used by lowsFactory
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > rcpA = rcpFromPtr(A);
  // Create the LOWSB object that will be used to solve the linear system
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > invertibleA =
    lowsFactory.createOp();
  // Initialize the invertible linear operator given the forward operator
  Thyra::initializeOp<Scalar>(lowsFactory, rcpA, invertibleA.ptr());
  // Solve the system using a default solve criteria using a non-member helper function
  Thyra::assign(x1, as<Scalar>(0.0));
  Thyra::solve<Scalar>(*invertibleA, Thyra::NOTRANS, b1, x1);
  // Before the forward operator A is changed it is recommended that you
  // uninitialize *invertibleA first to avoid accidental use of *invertiableA
  // while it may be in an inconsistent state from the time between *A changes
  // and *invertibleA is explicitly updated. However, this step is not
  // required!
  Thyra::uninitializeOp<Scalar>(lowsFactory, invertibleA.ptr());
  // Change the operator in some major way (perhaps even changing its structure)
  opMajorChanger.changeOp(A);
  // Recreate the LOWSB object and initialize it from scratch
  invertibleA = lowsFactory.createOp();
  Thyra::initializeOp<Scalar>(lowsFactory, rcpA, invertibleA.ptr());
  // Solve another set of linear systems
  Thyra::assign(x2, as<Scalar>(0.0));
  Thyra::solve<Scalar>(*invertibleA, Thyra::NOTRANS, b2, x2);
} // end solveMajorChangeSolve


//
// Individual externally preconditioned use cases
//


// begin createGeneralPreconditionedLinearOpWithSolve
/** \brief Create a <tt>Thyra::LinearOpWithSolveBase</tt> object given an
 * externally created preconditioner specified as a
 * <tt>Thyra::PreconditionerBase</tt> object.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
createGeneralPreconditionedLinearOpWithSolve(
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > &A,
  const Teuchos::RCP<const Thyra::PreconditionerBase<Scalar> > &P,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  Teuchos::FancyOStream &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating an externally preconditioned LinearOpWithSolveBase object ...\n";
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > invertibleA =
    lowsFactory.createOp();
  Thyra::initializePreconditionedOp<Scalar>(lowsFactory, A, P, invertibleA.ptr());
  out << "\nCreated LOWSB object:\n" << describe(*invertibleA, Teuchos::VERB_MEDIUM);
  return invertibleA;
} // end createGeneralPreconditionedLinearOpWithSolve


// begin createUnspecifiedPreconditionedLinearOpWithSolve
/** \brief Create a <tt>Thyra::LinearOpWithSolveBase</tt> object given an
 * externally created preconditioner specified as a single
 * <tt>Thyra::LinearOpBase</tt> object not targeted for the left or right.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
createUnspecifiedPreconditionedLinearOpWithSolve(
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > &A,
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > &P_op,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  Teuchos::FancyOStream &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating an LinearOpWithSolveBase object given a preconditioner operator"
      << " not targeted to the left or right ...\n";
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > invertibleA =
    lowsFactory.createOp();
  Thyra::initializePreconditionedOp<Scalar>(lowsFactory, A,
    Thyra::unspecifiedPrec<Scalar>(P_op), invertibleA.ptr());
  // Above, the lowsFactory object will decide whether to apply the single
  // preconditioner operator on the left or on the right.
  out << "\nCreated LOWSB object:\n" << describe(*invertibleA, Teuchos::VERB_MEDIUM);
  return invertibleA;
} // end createUnspecifiedPreconditionedLinearOpWithSolve


// begin createLeftPreconditionedLinearOpWithSolve
/** \brief Create a <tt>Thyra::LinearOpWithSolveBase</tt> object given an
 * externally created preconditioner specified as <tt>Thyra::LinearOpBase</tt>
 * object to be applied on the left.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
createLeftPreconditionedLinearOpWithSolve(
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > &A,
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > &P_op_left,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  Teuchos::FancyOStream &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating an LinearOpWithSolveBase object given a left preconditioner"
      << " operator ...\n";
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > invertibleA =
    lowsFactory.createOp();
  Thyra::initializePreconditionedOp<Scalar>(lowsFactory, A,
    Thyra::leftPrec<Scalar>(P_op_left), invertibleA.ptr());
  out << "\nCreated LOWSB object:\n" << describe(*invertibleA, Teuchos::VERB_MEDIUM);
  return invertibleA;
} // end createLeftPreconditionedLinearOpWithSolve


// begin createRightPreconditionedLinearOpWithSolve
/** \brief Create a <tt>Thyra::LinearOpWithSolveBase</tt> object given an
 * externally created preconditioner specified as <tt>Thyra::LinearOpBase</tt>
 * object to be applied on the right.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
createRightPreconditionedLinearOpWithSolve(
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > &A,
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > &P_op_right,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  Teuchos::FancyOStream &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating an LinearOpWithSolveBase object given a right"
      << " preconditioner operator ...\n";
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > invertibleA =
    lowsFactory.createOp();
  Thyra::initializePreconditionedOp<Scalar>(lowsFactory, A,
    Thyra::rightPrec<Scalar>(P_op_right), invertibleA.ptr());
  out << "\nCreated LOWSB object:\n" << describe(*invertibleA, Teuchos::VERB_MEDIUM);
  return invertibleA;
} // end createRightPreconditionedLinearOpWithSolve


// begin createLeftRightPreconditionedLinearOpWithSolve
/** \brief Create a <tt>Thyra::LinearOpWithSolveBase</tt> object given an
 * externally created preconditioner specified as left and right
 * <tt>Thyra::LinearOpBase</tt> objects.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
createLeftRightPreconditionedLinearOpWithSolve(
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > &A,
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > &P_op_left,
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > &P_op_right,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  Teuchos::FancyOStream &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating an LinearOpWithSolveBase object given a left and"
      << "right preconditioner operator ...\n";
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > invertibleA =
    lowsFactory.createOp();
  Thyra::initializePreconditionedOp<Scalar>(lowsFactory, A,
    Thyra::splitPrec<Scalar>(P_op_left, P_op_right), invertibleA.ptr());
  out << "\nCreated LOWSB object:\n" << describe(*invertibleA, Teuchos::VERB_MEDIUM);
  return invertibleA;
} // end createLeftRightPreconditionedLinearOpWithSolve


// begin createMatrixPreconditionedLinearOpWithSolve
/** \brief Create a <tt>Thyra::LinearOpWithSolveBase</tt> object given an
 * approximate forward operator that will be used to create the preconditioner
 * internally.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
createMatrixPreconditionedLinearOpWithSolve(
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > &A,
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > &A_approx,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  Teuchos::FancyOStream &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating a LinearOpWithSolveBase object given an approximate forward"
      << " operator to define the preconditioner ...\n";
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > invertibleA =
    lowsFactory.createOp();
  Thyra::initializeApproxPreconditionedOp<Scalar>(lowsFactory, A, A_approx,
    invertibleA.ptr());
  out << "\nCreated LOWSB object:\n" << describe(*invertibleA, Teuchos::VERB_MEDIUM);
  return invertibleA;
} // end createMatrixPreconditionedLinearOpWithSolve


// begin externalPreconditionerReuseWithSolves
/** \brief Example use case for preconditioner reuse
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
void externalPreconditionerReuseWithSolves(
  const Teuchos::Ptr<Thyra::LinearOpBase<Scalar> > &A_inout,
  const Thyra::LinearOpChanger<Scalar> &opChanger,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const Thyra::PreconditionerFactoryBase<Scalar> &precFactory,
  const Thyra::VectorBase<Scalar> &b1,
  const Teuchos::Ptr<Thyra::VectorBase<Scalar> > &x1,
  const Thyra::VectorBase<Scalar> &b2,
  const Teuchos::Ptr<Thyra::VectorBase<Scalar> > &x2,
  Teuchos::FancyOStream &out
  )
{
  using Teuchos::tab; using Teuchos::rcpFromPtr;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Teuchos::OSTab tab2(out);
  out << "\nShowing resuse of the preconditioner ...\n";
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > A = rcpFromPtr(A_inout);
  // Create the initial preconditioner for the input forward operator
  Teuchos::RCP<Thyra::PreconditionerBase<Scalar> > P =
    precFactory.createPrec();
  Thyra::initializePrec<Scalar>(precFactory, A, P.ptr());
  // Create the invertible LOWS object given the preconditioner
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > invertibleA =
    lowsFactory.createOp();
  Thyra::initializePreconditionedOp<Scalar>(lowsFactory, A, P, invertibleA.ptr());
  // Solve the first linear system
  assign(x1, ST::zero());
  Thyra::SolveStatus<Scalar> status1 = Thyra::solve<Scalar>(*invertibleA,
    Thyra::NOTRANS, b1, x1);
  out << "\nSolve status:\n" << status1;
  // Change the forward linear operator without changing the preconditioner
  opChanger.changeOp(A.ptr());
  // Warning! After the above change the integrity of the preconditioner
  // linear operators in P is undefined. For some implementations of the
  // preconditioner, its behavior will remain unchanged (e.g. ILU) which in
  // other cases the behavior will change but the preconditioner will still
  // work (e.g. Jacobi). However, there may be valid implementations where
  // the preconditioner will simply break if the forward operator that it is
  // based on breaks.
  //
  // Reinitialize the LOWS object given the updated forward operator A and the
  // old preconditioner P.
  Thyra::initializePreconditionedOp<Scalar>(lowsFactory, A, P, invertibleA.ptr());
  // Solve the second linear system
  assign(x2, ST::zero());
  Thyra::SolveStatus<Scalar>status2 = Thyra::solve<Scalar>(*invertibleA,
    Thyra::NOTRANS, b2, x2);
  out << "\nSolve status:\n" << status2;
} // end externalPreconditionerReuseWithSolves


//
// Combined use cases
//


/** \brief Combined use cases for the use of a
 * <tt>Thyra::LinearOpWithSolveFactoryBase</tt> object without an externally
 * set preconditioner.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
void nonExternallyPreconditionedLinearSolveUseCases(
  const Thyra::LinearOpBase<Scalar> &A,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  bool supportsAdjoints,
  Teuchos::FancyOStream &out
  )
{
  using Teuchos::as;
  Teuchos::OSTab tab(out);
  out << "\nRunning example use cases for a LinearOpWithSolveFactoryBase object ...\n";
  // Create a non-const A object (don't worry, it will not be changed)
  const Teuchos::Ptr<Thyra::LinearOpBase<Scalar> > A_nonconst =
    Teuchos::ptrFromRef(const_cast<Thyra::LinearOpBase<Scalar>&>(A));
  // Create the RHS (which is just a random set of coefficients)
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
    b1 = Thyra::createMember(A.range()),
    b2 = Thyra::createMember(A.range());
  Thyra::randomize( as<Scalar>(-1.0), as<Scalar>(+1.0), b1.ptr() );
  Thyra::randomize( as<Scalar>(-1.0), as<Scalar>(+1.0), b2.ptr() );
  // Create the LHS for the linear solve
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
    x1 = Thyra::createMember(A.domain()),
    x2 = Thyra::createMember(A.domain());
  // Perform a single, non-adjoint, linear solve
  singleLinearSolve(A, lowsFactory, *b1, x1.ptr(), out);
  // Creating a scaled adjoint LinearOpWithSolveBase object
  if(supportsAdjoints) {
    Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
      invertibleAdjointA = createScaledAdjointLinearOpWithSolve(
        Teuchos::rcp(&A,false),as<Scalar>(2.0),lowsFactory,out);
  }
  // Perform a solve, change the operator, and then solve again.
  solveNumericalChangeSolve<Scalar>(
    A_nonconst,                            // Don't worry, it will not be changed!
    Thyra::NullLinearOpChanger<Scalar>(),  // This object will not really change A!
    lowsFactory, *b1, x1.ptr(), *b2, x1.ptr(), out );
  // Perform a solve, change the operator in a very small way, and then solve again.
  solveSmallNumericalChangeSolve<Scalar>(
    A_nonconst,                            // Don't worry, it will not be changed!
    Thyra::NullLinearOpChanger<Scalar>(),  // This object will not really change A!
    lowsFactory, *b1, x1.ptr(), *b2, x1.ptr(), out );
  // Perform a solve, change the operator in a major way, and then solve again.
  solveMajorChangeSolve<Scalar>(
    A_nonconst,                            // Don't worry, it will not be changed!
    Thyra::NullLinearOpChanger<Scalar>(),  // This object will not really change A!
    lowsFactory, *b1, x1.ptr(), *b2, x1.ptr(), out );
}


/** \brief Combined use cases for the use of a
 * <tt>Thyra::LinearOpWithSolveFactoryBase</tt> object which use an externally
 * set preconditioner.
 * \ingroup thyra_operator_solve_support_LOWSF_examples
 */
template<class Scalar>
void externallyPreconditionedLinearSolveUseCases(
  const Thyra::LinearOpBase<Scalar> &A,
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const Thyra::PreconditionerFactoryBase<Scalar> &precFactory,
  const bool supportsLeftPrec,
  const bool supportsRightPrec,
  Teuchos::FancyOStream &out
  )
{
  using Teuchos::rcpFromRef; using Teuchos::as;
  Teuchos::OSTab tab(out);
  out << "\nRunning example use cases with an externally defined"
      << " preconditioner with a LinearOpWithSolveFactoryBase object ...\n";
  // Create a non-const A object (don't worry, it will not be changed)
  const Teuchos::Ptr<Thyra::LinearOpBase<Scalar> > A_nonconst =
    Teuchos::ptrFromRef(const_cast<Thyra::LinearOpBase<Scalar>&>(A));
  // Create the RHS (which is just a random set of coefficients)
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
    b1 = Thyra::createMember(A.range()),
    b2 = Thyra::createMember(A.range());
  Thyra::randomize( as<Scalar>(-1.0), as<Scalar>(+1.0), b1.ptr() );
  Thyra::randomize( as<Scalar>(-1.0), as<Scalar>(+1.0), b2.ptr() );
  // Create the LHS for the linear solve
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
    x1 = Thyra::createMember(A.domain()),
    x2 = Thyra::createMember(A.domain());
  // Create a preconditioner for the input forward operator
  Teuchos::RCP<Thyra::PreconditionerBase<Scalar> >
    P = precFactory.createPrec();
  Thyra::initializePrec<Scalar>(precFactory, rcpFromRef(A), P.ptr());
  // Above, we don't really know the nature of the preconditioner. It could a
  // single linear operator to be applied on the left or the right or it could
  // be a split preconditioner with different linear operators to be applied
  // on the right or left. Or, it could be a single linear operator that is
  // not targeted to the left or the right.
  //
  // Create a LOWSB object given the created preconditioner
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleA = createGeneralPreconditionedLinearOpWithSolve<Scalar>(
      Teuchos::rcp(&A, false), P.getConst(), lowsFactory, out);
  // Grab a preconditioner operator out of the preconditioner object
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > P_op;
  if (nonnull(P_op=P->getUnspecifiedPrecOp()));
  else if (nonnull(P_op=P->getLeftPrecOp()));
  else if (nonnull(P_op=P->getRightPrecOp()));
  // Create a LOWSB object given an unspecified preconditioner operator
  invertibleA = createUnspecifiedPreconditionedLinearOpWithSolve(
    rcpFromRef(A), P_op, lowsFactory, out);
  // Create a LOWSB object given a left preconditioner operator
  if(supportsLeftPrec) {
    invertibleA = createLeftPreconditionedLinearOpWithSolve(
      rcpFromRef(A), P_op, lowsFactory,out);
  }
  // Create a LOWSB object given a right preconditioner operator
  if(supportsRightPrec) {
    invertibleA = createRightPreconditionedLinearOpWithSolve(
      rcpFromRef(A), P_op, lowsFactory, out);
  }
  // Create a LOWSB object given (bad set of) left and right preconditioner
  // operators
  if( supportsLeftPrec && supportsRightPrec ) {
    invertibleA = createLeftRightPreconditionedLinearOpWithSolve(
      rcpFromRef(A), P_op, P_op, lowsFactory, out);
  }
  // Create a LOWSB object given a (very good) approximate forward linear
  // operator to construct the preconditoner from..
  invertibleA = createMatrixPreconditionedLinearOpWithSolve<Scalar>(
    rcpFromRef(A), rcpFromRef(A), lowsFactory,out);
  // Preconditioner reuse example
  externalPreconditionerReuseWithSolves<Scalar>(
    A_nonconst,                            // Don't worry, it will not be changed!
    Thyra::NullLinearOpChanger<Scalar>(),  // This object will not really change A!
    lowsFactory, precFactory,
    *b1, x1.ptr(), *b2, x2.ptr(), out );
}


#endif // THYRA_LINEAR_OP_WITH_SOLVE_EXAMPLES_HPP
