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
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"

namespace Thyra {

//
// Helper code
//

/** \brief Abstract strategy interface for changing
 * <tt>Thyra::LinearOpBase</tt> objects.
 */
template<class Scalar>
class LinearOpChanger {
public:
  /** \brief . */
  virtual ~LinearOpChanger() {}
  /** \brief . */
  virtual void changeOp( LinearOpBase<Scalar> *op ) const = 0;
};

/** \brief No-op changer. */
template<class Scalar>
class NullLinearOpChanger : public LinearOpChanger<Scalar> {
public:
  /** \brief . */
  void changeOp( LinearOpBase<Scalar> *op ) const {}
};

} // namespace Thyra

//
// Individual non-externally preconditioned use cases
//
  
// begin singleLinearSolve
/** \brief Performing a single linear solve given a forward operator. */
template<class Scalar>
void singleLinearSolve(
  const Thyra::LinearOpBase<Scalar>                    &A
  ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>   &lowsFactory
  ,const Thyra::VectorBase<Scalar>                     &b
  ,Thyra::VectorBase<Scalar>                           *x
  ,Teuchos::FancyOStream                               &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nPerforming a single linear solve ...\n";
  // Create the LOWSB object that will be used to solve the linear system
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleA = Thyra::createAndInitializeLinearOpWithSolve(lowsFactory,Teuchos::rcp(&A,false));
  // Solve the system using a default solve criteria using a non-member helper function 
  assign(&*x,Teuchos::ScalarTraits<Scalar>::zero()); // Must initialize to a guess before solve!
  Thyra::SolveStatus<Scalar>
    status = Thyra::solve(*invertibleA,Thyra::NOTRANS,b,x);
  out << "\nSolve status:\n" << status;
} // end singleLinearSolve

// begin createScaledAdjointLinearOpWithSolve
/** \brief Create a <tt>LinearOpWithSolveBase</tt> object from an implicitly
 * scaled adjoint <tt>LinearOpBase</tt> object.
 */
template<class Scalar>
Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
createScaledAdjointLinearOpWithSolve(
  const Thyra::LinearOpBase<Scalar>                    &A
  ,const Scalar                                        &scalar
  ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>   &lowsFactory
  ,Teuchos::FancyOStream                               &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating a scaled adjoint LinearOpWithSolveBase object ...\n";
  // Create the LOWSB object that will be used to solve the linear system
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleAdjointA
    = Thyra::createAndInitializeLinearOpWithSolve(
      lowsFactory
      ,Thyra::scale(scalar,Thyra::adjoint(Teuchos::rcp(&A,false)))
      );
  out << "\nCreated LOWSB object:\n" << describe(*invertibleAdjointA,Teuchos::VERB_MEDIUM);
  return invertibleAdjointA;
} // end createScaledAdjointLinearOpWithSolve

// begin solveNumericalChangeSolve
/** \brief Create a <tt>LinearOpWithSolveBase</tt> object and perform a solve,
 * then change the <tt>LinearOpBase</tt> object and solve again.
 */
template<class Scalar>
void solveNumericalChangeSolve(
  Thyra::LinearOpBase<Scalar>                          *A
  ,const Thyra::LinearOpChanger<Scalar>                &opChanger
  ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>   &lowsFactory
  ,const Thyra::VectorBase<Scalar>                     &b1
  ,Thyra::VectorBase<Scalar>                           *x1
  ,const Thyra::VectorBase<Scalar>                     &b2
  ,Thyra::VectorBase<Scalar>                           *x2
  ,Teuchos::FancyOStream                               &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nPerforming a solve, changing the operator, then performing another solve ...\n";
  // Get a local non-owned RCP to A to be used by lowsFactory
  Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
    rcpA = Teuchos::rcp(A,false); // Note: This is the right way to do this!
  // Create the LOWSB object that will be used to solve the linear system
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleA = lowsFactory.createOp();
  // Initialize the invertible linear operator given the forward operator
  lowsFactory.initializeOp(rcpA,&*invertibleA);
  // Solve the system using a default solve criteria using a non-member helper function
  Thyra::assign(x1,Scalar(0.0));
  solve(*invertibleA,Thyra::NOTRANS,b1,x1);
  // Before the forward operator A is changed it is recommended that you
  // uninitialize *invertibleA first to avoid accidental use of *invertiableA
  // while it may be in an inconsistent state from the time between *A changes
  // and *invertibleA is explicitly updated.  However, this step is not
  // required!
  lowsFactory.uninitializeOp(&*invertibleA);
  // Change the operator and reinitialize the invertible operator
  opChanger.changeOp(A);
  lowsFactory.initializeOp(rcpA,&*invertibleA);
  // Note that above will reuse any factorization structures that may have been
  // created in the first call to initializeOp(...)
  // Solve another linear system with new values of A
  Thyra::assign(x2,Scalar(0.0));
  solve(*invertibleA,Thyra::NOTRANS,b2,x2);
} // end solveNumericalChangeSolve

// begin solveSmallNumericalChangeSolve
/** \brief Create a <tt>LinearOpWithSolveBase</tt> object and perform a solve,
 * then change the <tt>LinearOpBase</tt> object in a very small way and solve
 * again.
 */
template<class Scalar>
void solveSmallNumericalChangeSolve(
  Thyra::LinearOpBase<Scalar>                          *A
  ,const Thyra::LinearOpChanger<Scalar>                &opSmallChanger
  ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>   &lowsFactory
  ,const Thyra::VectorBase<Scalar>                     &b1
  ,Thyra::VectorBase<Scalar>                           *x1
  ,const Thyra::VectorBase<Scalar>                     &b2
  ,Thyra::VectorBase<Scalar>                           *x2
  ,Teuchos::FancyOStream                               &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nPerforming a solve, changing the operator in a very small way, then performing another solve ...\n";
  // Get a local non-owned RCP to A to be used by lowsFactory
  Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
    rcpA = Teuchos::rcp(A,false); // Note: This is the right way to do this!
  // Create the LOWSB object that will be used to solve the linear system
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleA = lowsFactory.createOp();
  // Initialize the invertible linear operator given the forward operator
  lowsFactory.initializeOp(rcpA,&*invertibleA);
  // Solve the system using a default solve criteria using a non-member helper function
  Thyra::assign(x1,Scalar(0.0));
  solve(*invertibleA,Thyra::NOTRANS,b1,x1);
  // Before the forward operator A is changed it is recommended that you
  // uninitialize *invertibleA first to avoid accidental use of *invertiableA
  // while it may be in an inconsistent state from the time between *A changes
  // and *invertibleA is explicitly updated.  However, this step is not
  // required!
  lowsFactory.uninitializeOp(&*invertibleA);
  // Change the operator and reinitialize the invertible operator
  opSmallChanger.changeOp(A);
  lowsFactory.initializeAndReuseOp(rcpA,&*invertibleA);
  // Note that above a maximum amount of reuse will be achieved, such as
  // keeping the same precondtioner.
  Thyra::assign(x2,Scalar(0.0));
  solve(*invertibleA,Thyra::NOTRANS,b2,x2);
} // end solveSmallNumericalChangeSolve

// begin solveMajorChangeSolve
/** \brief Create a <tt>LinearOpWithSolveBase</tt> object and perform a solve,
 * then change the <tt>LinearOpBase</tt> object in a very small way and solve
 * again.
 */
template<class Scalar>
void solveMajorChangeSolve(
  Thyra::LinearOpBase<Scalar>                          *A
  ,const Thyra::LinearOpChanger<Scalar>                &opSmallChanger
  ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>   &lowsFactory
  ,const Thyra::VectorBase<Scalar>                     &b1
  ,Thyra::VectorBase<Scalar>                           *x1
  ,const Thyra::VectorBase<Scalar>                     &b2
  ,Thyra::VectorBase<Scalar>                           *x2
  ,Teuchos::FancyOStream                               &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nPerforming a solve, changing the operator in a major way, then performing another solve ...\n";
  // Get a local non-owned RCP to A to be used by lowsFactory
  Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >
    rcpA = Teuchos::rcp(A,false); // Note: This is the right way to do this!
  // Create the LOWSB object that will be used to solve the linear system
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleA = lowsFactory.createOp();
  // Initialize the invertible linear operator given the forward operator
  lowsFactory.initializeOp(rcpA,&*invertibleA);
  // Solve the system using a default solve criteria using a non-member helper function
  Thyra::assign(x1,Scalar(0.0));
  solve(*invertibleA,Thyra::NOTRANS,b1,x1);
  // Before the forward operator A is changed it is recommended that you
  // uninitialize *invertibleA first to avoid accidental use of *invertiableA
  // while it may be in an inconsistent state from the time between *A changes
  // and *invertibleA is explicitly updated.  However, this step is not
  // required!
  lowsFactory.uninitializeOp(&*invertibleA);
  // Change the operator in some major way (perhaps even changing its structure)
  opSmallChanger.changeOp(A);
  // Recreate the LOWSB object and initialize it from scratch
  invertibleA = lowsFactory.createOp();
  lowsFactory.initializeOp(rcpA,&*invertibleA);
  // Solve another set of linear systems
  Thyra::assign(x2,Scalar(0.0));
  solve(*invertibleA,Thyra::NOTRANS,b2,x2);
} // end solveMajorChangeSolve

//
// Individual externally preconditioned use cases
//

/** \brief Create a <tt>Thyra::LinearOpWithSolveBase</tt> object given an
 * externally created preconditioner specified as a
 * <tt>Thyra::PreconditionerBase</tt> object.
 */
template<class Scalar>
Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
createGeneralPreconditionedLinearOpWithSolve(
  const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >        &A
  ,const Teuchos::RefCountPtr<const Thyra::PreconditionerBase<Scalar> > &P
  ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>                    &lowsFactory
  ,Teuchos::FancyOStream                                                &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating an externally preconditioned LinearOpWithSolveBase object ...\n";
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleA = lowsFactory.createOp();
  lowsFactory.initializePreconditionedOp(A,P,&*invertibleA);
  out << "\nCreated LOWSB object:\n" << describe(*invertibleA,Teuchos::VERB_MEDIUM);
  return invertibleA;
}

/** \brief Create a <tt>Thyra::LinearOpWithSolveBase</tt> object given an
 * externally created preconditioner specified as a single
 * <tt>Thyra::LinearOpBase</tt> object not targeted for the left or right.
 */
template<class Scalar>
Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
createUnspecifiedPreconditionedLinearOpWithSolve(
  const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >      &A
  ,const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >     &P_op
  ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>                  &lowsFactory
  ,Teuchos::FancyOStream                                              &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating an LinearOpWithSolveBase object given a preconditioner operator ...\n";
  Teuchos::RefCountPtr<Thyra::PreconditionerBase<Scalar> >
    P = Teuchos::rcp(new Thyra::DefaultPreconditioner<Scalar>(P_op));
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleA = lowsFactory.createOp();
  lowsFactory.initializePreconditionedOp(A,P,&*invertibleA);
  out << "\nCreated LOWSB object:\n" << describe(*invertibleA,Teuchos::VERB_MEDIUM);
  return invertibleA;
}

/** \brief Create a <tt>Thyra::LinearOpWithSolveBase</tt> object given an
 * externally created preconditioner specified as <tt>Thyra::LinearOpBase</tt>
 * object to be applied on the left.
 */
template<class Scalar>
Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
createLeftPreconditionedLinearOpWithSolve(
  const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >      &A
  ,const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >     &P_op_left
  ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>                  &lowsFactory
  ,Teuchos::FancyOStream                                              &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating an LinearOpWithSolveBase object given a left preconditioner operator ...\n";
  Teuchos::RefCountPtr<Thyra::PreconditionerBase<Scalar> >
    P = Teuchos::rcp(new Thyra::DefaultPreconditioner<Scalar>(P_op_left,Teuchos::null));
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleA = lowsFactory.createOp();
  lowsFactory.initializePreconditionedOp(A,P,&*invertibleA);
  out << "\nCreated LOWSB object:\n" << describe(*invertibleA,Teuchos::VERB_MEDIUM);
  return invertibleA;
}

/** \brief Create a <tt>Thyra::LinearOpWithSolveBase</tt> object given an
 * externally created preconditioner specified as <tt>Thyra::LinearOpBase</tt>
 * object to be applied on the right.
 */
template<class Scalar>
Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
createRightPreconditionedLinearOpWithSolve(
  const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >      &A
  ,const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >     &P_op_right
  ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>                  &lowsFactory
  ,Teuchos::FancyOStream                                              &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating an LinearOpWithSolveBase object given a right preconditioner operator ...\n";
  Teuchos::RefCountPtr<Thyra::PreconditionerBase<Scalar> >
    P = Teuchos::rcp(new Thyra::DefaultPreconditioner<Scalar>(Teuchos::null,P_op_right));
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleA = lowsFactory.createOp();
  lowsFactory.initializePreconditionedOp(A,P,&*invertibleA);
  out << "\nCreated LOWSB object:\n" << describe(*invertibleA,Teuchos::VERB_MEDIUM);
  return invertibleA;
}

/** \brief Create a <tt>Thyra::LinearOpWithSolveBase</tt> object given an
 * externally created preconditioner specified as left and right
 * <tt>Thyra::LinearOpBase</tt> objects.
 */
template<class Scalar>
Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
createLeftRightPreconditionedLinearOpWithSolve(
  const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >      &A
  ,const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >     &P_op_left
  ,const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >     &P_op_right
  ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>                  &lowsFactory
  ,Teuchos::FancyOStream                                              &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating an LinearOpWithSolveBase object given a left and right preconditioner operator ...\n";
  Teuchos::RefCountPtr<Thyra::PreconditionerBase<Scalar> >
    P = Teuchos::rcp(new Thyra::DefaultPreconditioner<Scalar>(P_op_left,Teuchos::null));
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleA = lowsFactory.createOp();
  lowsFactory.initializePreconditionedOp(A,P,&*invertibleA);
  out << "\nCreated LOWSB object:\n" << describe(*invertibleA,Teuchos::VERB_MEDIUM);
  return invertibleA;
}

/** \brief Create a <tt>Thyra::LinearOpWithSolveBase</tt> object given an
 * approximate forward operator that will be used to create the preconditioner
 * internally.
 */
template<class Scalar>
Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
createMatrixPreconditionedLinearOpWithSolve(
  const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >      &A
  ,const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> >     &A_approx
  ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>                  &lowsFactory
  ,Teuchos::FancyOStream                                              &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nCreating a LinearOpWithSolveBase object given an approximate forward"
      << " operator to define the preconditioner ...\n";
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleA = lowsFactory.createOp();
  lowsFactory.initializeApproxPreconditionedOp(A,A_approx,&*invertibleA);
  out << "\nCreated LOWSB object:\n" << describe(*invertibleA,Teuchos::VERB_MEDIUM);
  return invertibleA;
}

//
// Combined use cases
//

/** \brief Combined use cases for the use of a
 * <tt>Thyra::LinearOpWithSolveFactoryBase</tt> object without an externally
 * set preconditioner.
 */
template<class Scalar>
void nonExternallyPreconditionedLinearSolveUseCases(
  const Thyra::LinearOpBase<Scalar>                    &A
  ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>   &lowsFactory
  ,Teuchos::FancyOStream                               &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nRunning example use cases for a LinearOpWithSolveFactoryBase object ...\n";
  // Create the RHS (which is just a random set of coefficients)
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
    b1 = Thyra::createMember(A.range()),
    b2 = Thyra::createMember(A.range());
  Thyra::randomize( Scalar(-1.0), Scalar(+1.0), &*b1 );
  Thyra::randomize( Scalar(-1.0), Scalar(+1.0), &*b2 );
  // Create the LHS for the linear solve
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
    x1 = Thyra::createMember(A.domain()),
    x2 = Thyra::createMember(A.domain());
  // Perform a single, non-adjoint, linear solve
  singleLinearSolve(A,lowsFactory,*b1,&*x1,out);
  // Creating a scaled adjoint LinearOpWithSolveBase object
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleAdjointA = createScaledAdjointLinearOpWithSolve(
      A,Scalar(2.0),lowsFactory,out);
  // Perform a solve, change the operator, and then solve again.
  solveNumericalChangeSolve(
    const_cast<Thyra::LinearOpBase<Scalar>*>(&A) // Don't worry, it will not be changed!
    ,Thyra::NullLinearOpChanger<Scalar>()        // This object will not really change A!
    ,lowsFactory,*b1,&*x1,*b2,&*x1,out
    );
  // Perform a solve, change the operator in a very small way, and then solve again.
  solveSmallNumericalChangeSolve(
    const_cast<Thyra::LinearOpBase<Scalar>*>(&A) // Don't worry, it will not be changed!
    ,Thyra::NullLinearOpChanger<Scalar>()        // This object will not really change A!
    ,lowsFactory,*b1,&*x1,*b2,&*x1,out
    );
  // Perform a solve, change the operator in a major way, and then solve again.
  solveMajorChangeSolve(
    const_cast<Thyra::LinearOpBase<Scalar>*>(&A) // Don't worry, it will not be changed!
    ,Thyra::NullLinearOpChanger<Scalar>()        // This object will not really change A!
    ,lowsFactory,*b1,&*x1,*b2,&*x1,out
    );
}

/** \brief Combined use cases for the use of a
 * <tt>Thyra::LinearOpWithSolveFactoryBase</tt> object which use an externally
 * set preconditioner.
 */
template<class Scalar>
void externallyPreconditionedLinearSolveUseCases(
  const Thyra::LinearOpBase<Scalar>                    &A
  ,const Thyra::LinearOpWithSolveFactoryBase<Scalar>   &lowsFactory
  ,const Thyra::PreconditionerFactoryBase<Scalar>      &precFactory
  ,Teuchos::FancyOStream                               &out
  )
{
  Teuchos::OSTab tab(out);
  out << "\nRunning example use cases with an externally defined"
      << " preconditioner with a LinearOpWithSolveFactoryBase object ...\n";
  // Create the RHS (which is just a random set of coefficients)
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
    b1 = Thyra::createMember(A.range()),
    b2 = Thyra::createMember(A.range());
  Thyra::randomize( Scalar(-1.0), Scalar(+1.0), &*b1 );
  Thyra::randomize( Scalar(-1.0), Scalar(+1.0), &*b2 );
  // Create the LHS for the linear solve
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
    x1 = Thyra::createMember(A.domain()),
    x2 = Thyra::createMember(A.domain());
  // Create a preconditioner for the input forward operator
  Teuchos::RefCountPtr<Thyra::PreconditionerBase<Scalar> >
    P = precFactory.createPrec();
  precFactory.initializePrec(Teuchos::rcp(&A,false),&*P);
  // Above, we don't really know the nature of the preconditioner.  It could a
  // single linear operator to be applied on the left or the right or it could
  // be a split preconditioner with different linear operators to be applied
  // on the right or left.  Or, it could be a single linear operator that is
  // not targeted to the left or the right.
  //
  // Create a LOWSB object given the created preconditioner
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
    invertibleA = createGeneralPreconditionedLinearOpWithSolve<Scalar>(
      Teuchos::rcp(&A,false),P,lowsFactory,out);
  // Grab a preconditioner operator out of the preconditioner object
  Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> > P_op;
  if((P_op=P->getUnspecifiedPrecOp()).get());
  else if((P_op=P->getLeftPrecOp()).get());
  else if((P_op=P->getRightPrecOp()).get());
  // Create a LOWSB object given an unspecified preconditioner operator
  invertibleA = createUnspecifiedPreconditionedLinearOpWithSolve(
    Teuchos::rcp(&A,false),P_op,lowsFactory,out);
  // Create a LOWSB object given a left preconditioner operator
  invertibleA = createLeftPreconditionedLinearOpWithSolve(
    Teuchos::rcp(&A,false),P_op,lowsFactory,out);
  // Create a LOWSB object given a right preconditioner operator
  invertibleA = createRightPreconditionedLinearOpWithSolve(
    Teuchos::rcp(&A,false),P_op,lowsFactory,out);
  // Create a LOWSB object given (bad set of) left and right preconditioner operators
  invertibleA = createLeftRightPreconditionedLinearOpWithSolve(
    Teuchos::rcp(&A,false),P_op,P_op,lowsFactory,out);
  // Create a LOWSB object given a (very good) approximate forward linear
  // operator to construct the preconditoner out of.
  invertibleA = createMatrixPreconditionedLinearOpWithSolve<Scalar>(
    Teuchos::rcp(&A,false),Teuchos::rcp(&A,false),lowsFactory,out);
}

#endif // THYRA_LINEAR_OP_WITH_SOLVE_EXAMPLES_HPP
