// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX.H"
#include "Teuchos_ParameterList.hpp"

#include "NOX_Petsc_Group.H"    // class definition
#include "NOX_Petsc_Interface.H"
#include "NOX_Petsc_Vector.H"
#include "NOX_Petsc_SharedJacobian.H"
#include "Teuchos_ParameterList.hpp"

// External include files - linking to Petsc
#include "petscversion.h"
#include "petscksp.h"

using namespace NOX;
using namespace NOX::Petsc;

Group::Group(Interface& i, Vec& x, Mat& J) :
  xVector(x, "Solution"), // deep copy x
  RHSVector(x, "RHS", ShapeCopy), // new vector of same size
  gradVector(x, "Grad", ShapeCopy), // new vector of same size
  NewtonVector(x, "Newton", ShapeCopy), // new vector of same size
  sharedJacobianPtr(new SharedJacobian(J)), // pass J to SharedJacobian
  sharedJacobian(*sharedJacobianPtr), // pass J to SharedJacobian
  jacType("User Supplied"), // the only option for now
  userInterface(i)
{
  resetIsValid();
}

Group::Group(const Group& source, CopyType type) :
  xVector(source.xVector, type),
  RHSVector(source.RHSVector, type),
  gradVector(source.gradVector, type),
  NewtonVector(source.NewtonVector, type),
  sharedJacobianPtr(NULL),
  sharedJacobian(source.sharedJacobian),
  jacType(source.jacType),
  userInterface(source.userInterface)
{
  switch (type) {

  case DeepCopy:

    isValidRHS = source.isValidRHS;
    isValidGrad = source.isValidGrad;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    isValidPreconditioner = source.isValidPreconditioner;
    normRHS = source.normRHS;

    // New copy takes ownership of the shared Jacobian
    if (isValidJacobian)
      sharedJacobian.getJacobian(this);

    break;

  case ShapeCopy:
    resetIsValid();
    break;

  default:
    std::cerr << "ERROR: Invalid ConstructorType for group copy constructor." << std::endl;
    throw std::runtime_error("NOX Error");
  }

}

Group::~Group()
{
  delete sharedJacobianPtr;
}

void Group::resetIsValid() //private
{
  isValidRHS = false;
  isValidJacobian = false;
  isValidGrad = false;
  isValidNewton = false;
  isValidPreconditioner = false;
}

Teuchos::RCP<NOX::Abstract::Group>
Group::clone(CopyType type) const
{
  Teuchos::RCP<NOX::Abstract::Group> newgrp =
    Teuchos::rcp(new NOX::Petsc::Group(*this, type));
  return newgrp;
}

Abstract::Group&
Group::operator=(const Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

Abstract::Group&
Group::operator=(const Group& source)
{
  // Copy the xVector
  xVector = source.xVector;

  // Copy reference to sharedJacobian
  sharedJacobian = source.sharedJacobian;

  // Update the isValidVectors
  isValidRHS = source.isValidRHS;
  isValidGrad = source.isValidGrad;
  isValidNewton = source.isValidNewton;
  isValidJacobian = source.isValidJacobian;
  isValidPreconditioner = source.isValidPreconditioner;

  // Only copy vectors that are valid
  if (isValidRHS) {
    RHSVector = source.RHSVector;
    normRHS = source.normRHS;
  }

  if (isValidGrad)
    gradVector = source.gradVector;

  if (isValidNewton)
    NewtonVector = source.NewtonVector;

  // If valid, this takes ownership of the shared Jacobian
  if (isValidJacobian)
    sharedJacobian.getJacobian(this);

  return *this;
}

void
Group::setX(const Abstract::Vector& y)
{
  setX(dynamic_cast<const Vector&> (y));
  return;
}

void
Group::setX(const Vector& y)
{
  resetIsValid();
  xVector = y;
  return;
}

void
Group::computeX(const Abstract::Group& grp, const Abstract::Vector& d, double step)
{
  // Cast to appropriate type, then call the "native" computeX
  const Group& petscgrp = dynamic_cast<const Group&> (grp);
  const Vector& petscd = dynamic_cast<const Vector&> (d);
  computeX(petscgrp, petscd, step);
  return;
}

void
Group::computeX(const Group& grp, const Vector& d, double step)
{
  resetIsValid();
  xVector.update(1.0, grp.xVector, step, d);
  return;
}

Abstract::Group::ReturnType
Group::computeF()
{
  if (isF())
    return Abstract::Group::Ok;

  bool status = false;

  status = userInterface.computeF(xVector.getPetscVector(), RHSVector.getPetscVector());

  if(status == false) {
    std::cout << "ERROR: Petsc::Group::computeF() - fill failed!!!"
         << std::endl;
    throw std::runtime_error("NOX Error:) RHS Fill Failed";
  }

  normRHS = RHSVector.norm();

  isValidRHS = true;

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType
Group::computeJacobian()
{
  // Skip if the Jacobian is already valid
  if (isJacobian())
    return Abstract::Group::Ok;

  // Take ownership of the Jacobian and get a reference
  Mat& Jacobian = sharedJacobian.getJacobian(this);

  // Fill the Jacobian

  bool status = false;

  if(jacType == "User Supplied") {

    status = userInterface.computeJacobian(xVector.getPetscVector(), Jacobian);

    if (status == false) {
      std::cout << "ERROR: Petsc::Group::computeJacobian() - fill failed!!!"
           << std::endl;
      throw std::runtime_error("NOX Error:) Jacobian Fill Failed";
    }
  }
  else if(jacType == "Finite Difference") {
    std::cout << "Finite Difference evaluation not yet supported !!\n\n";
    throw std::runtime_error("NOX Error");
  }

  // Update status
  isValidJacobian = true;

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType
Group::computeGradient()
{
  if (isGradient())
    return Abstract::Group::Ok;

  if (!isF()) {
    std::cerr << "ERROR: NOX::Petsc::Group::computeGradient() - RHS is out of date wrt X!" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  if (!isJacobian()) {
    std::cerr << "ERROR: NOX::Petsc::Group::computeGradient() - Jacobian is out of date wrt X!" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  // Get a reference to the Jacobian (it's validity was checked above)
  const Mat& Jacobian = sharedJacobian.getJacobian();

  // Compute grad = Jacobian^T * RHS.
  // Need to add a check on ierr
  PetscErrorCode ierr = MatMultTranspose(Jacobian, RHSVector.getPetscVector(), gradVector.getPetscVector());
  TEUCHOS_ASSERT(ierr == 0);

  // Update state
  isValidGrad = true;

  // Return result
  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType
Group::computeNewton(Teuchos::ParameterList& p)
{
  if (isNewton())
    return Abstract::Group::Ok;

  if (!isF()) {
    std::cerr << "ERROR: NOX::Petsc::Group::computeNewton() - invalid RHS" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  if (!isJacobian()) {
    std::cerr << "ERROR: NOX::Petsc::Group::computeNewton() - invalid Jacobian" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  // Get the Jacobian
  Mat& Jacobian = sharedJacobian.getJacobian(this);

  // Create Petsc KSP problem for the linear solve

  KSP ksp;
  KSPConvergedReason reason;

  int ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
#if  (PETSC_VERSION_MAJOR >= 3) || (PETSC_VERSION_MINOR >= 5)
  ierr = KSPSetOperators(ksp,Jacobian,Jacobian);
#else
  ierr = KSPSetOperators(ksp,Jacobian,Jacobian, DIFFERENT_NONZERO_PATTERN);
#endif

  /*
     Set runtime options (e.g., -ksp_type <type> -pc_type <type>)
  */

  ierr = KSPSetFromOptions(ksp);

  /*
     Solve linear system.  Here we explicitly call KSPSetUp() for more
     detailed performance monitoring of certain preconditioners, such
     as ICC and ILU.  This call is optional, as KSPSetUp() will
     automatically be called within KSPSolve() if it hasn't been
     called already.
  */
  ierr = KSPSetUp(ksp);

  int its = 0;
  // Need to put a check on ierr.
  ierr = KSPSolve( ksp, RHSVector.getPetscVector(), NewtonVector.getPetscVector() );

  // Ascertain convergence status
  ierr = KSPGetConvergedReason(ksp,&reason);

  if( reason == KSP_DIVERGED_INDEFINITE_PC )
  {
    std::cout << "\nDivergence because of indefinite preconditioner;\n";
    std::cout << "Run the executable again but with -pc_ilu_shift option.\n";
  }
  else if( reason < 0 )
  {
    std::cout << "\nOther kind of divergence: this should not happen.\n";
  }
  else
  {
    ierr = KSPGetIterationNumber( ksp, &its );
    std::cout << "\nConvergence in " << its << " iterations.\n";
  }

  TEUCHOS_ASSERT(ierr == 0);

  // Scale soln by -1
  NewtonVector.scale(-1.0);

  // Update state
  isValidNewton = true;

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType
Group::applyJacobian(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& petscinput = dynamic_cast<const Vector&> (input);
  Vector& petscresult = dynamic_cast<Vector&> (result);
  return applyJacobian(petscinput, petscresult);
}

Abstract::Group::ReturnType
Group::applyJacobian(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian())
    return Abstract::Group::BadDependency;

  // Get a reference to the Jacobian (it's validity was checked above)
  const Mat& Jacobian = sharedJacobian.getJacobian();

  // Apply the Jacobian
  MatMult(Jacobian, input.getPetscVector(), result.getPetscVector());

  return Abstract::Group::Ok;
}


Abstract::Group::ReturnType
Group::applyRightPreconditioning(Teuchos::ParameterList& params,
                                 const Abstract::Vector& input,
                                 Abstract::Vector& result) const
{
  const Vector& petscinput = dynamic_cast<const Vector&> (input);
  Vector& petscresult = dynamic_cast<Vector&> (result);
  return applyRightPreconditioning(petscinput, petscresult);
}

Abstract::Group::ReturnType
Group::applyRightPreconditioning(const Vector& input, Vector& result) const
{
  if (!isJacobian())
    return Abstract::Group::BadDependency;

  // Get a reference to the Jacobian
  const Mat& Jacobian = sharedJacobian.getJacobian();

  // Get petsc reference to the result vector
  Vec& r = result.getPetscVector();

  // Set up preconditioner context
  PC pc;

  int ierr = PCCreate(PETSC_COMM_WORLD,&pc);

  // Here a default to jacobi (jacobian-diagonal-inverse) is established
  // but can be overridden via specification of pc_type in .petscrc

  // This allows more general preconditioning via specification of -pc_type
  // in .petscrc
  ierr = PCSetFromOptions(pc);

  /*
     Set operators and vector. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
#if  (PETSC_VERSION_MAJOR >= 3) || (PETSC_VERSION_MINOR >= 5)
  ierr = PCSetOperators(pc,Jacobian,Jacobian);
#else
  ierr = PCSetOperators(pc,Jacobian,Jacobian, DIFFERENT_NONZERO_PATTERN);
#endif
  ierr = PCSetUp(pc);

  // Apply the preconditioner
  ierr = PCApply(pc,input.getPetscVector(),r);

  // Cleanup
#if  (PETSC_VERSION_MAJOR >= 3) || (PETSC_VERSION_MINOR >= 5)
  ierr = PCDestroy(&pc);
#else
  ierr = PCDestroy(pc);
#endif

  TEUCHOS_ASSERT(ierr == 0);

  return Abstract::Group::Ok;
}


Abstract::Group::ReturnType
Group::applyJacobianTranspose(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& petscinput = dynamic_cast<const Vector&> (input);
  Vector& petscresult = dynamic_cast<Vector&> (result);
  return applyJacobianTranspose(petscinput, petscresult);
}

Abstract::Group::ReturnType
Group::applyJacobianTranspose(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian())
    return Abstract::Group::BadDependency;

  // Get a reference to the Jacobian (it's validity was check above)
  const Mat& Jacobian = sharedJacobian.getJacobian();

  // Apply the Jacobian
  int ierr = MatMultTranspose(Jacobian, input.getPetscVector(), result.getPetscVector());

  TEUCHOS_ASSERT(ierr == 0);

  return Abstract::Group::Ok;
}


bool
Group::isF() const
{
  return isValidRHS;
}

bool
Group::isJacobian() const
{
  return ((sharedJacobian.isOwner(this)) && (isValidJacobian));
}

bool
Group::isGradient() const
{
  return isValidGrad;
}

bool
Group::isNewton() const
{
  return isValidNewton;
}

bool
Group::isPreconditioner() const
{
  return isValidPreconditioner;
}

const Abstract::Vector&
Group::getX() const
{
  return xVector;
}

const Abstract::Vector&
Group::getF() const
{
  if (!isF()) {
    std::cerr << "ERROR: NOX::Petsc::Group::getF() - invalid RHS" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return RHSVector;
}

double
Group::getNormF() const
{
  if (!isF()) {
    std::cerr << "ERROR: NOX::Petsc::Group::getNormF() - invalid RHS" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return normRHS;
}

const Abstract::Vector&
Group::getGradient() const
{
  if (!isGradient()) {
    std::cerr << "ERROR: NOX::Petsc::Group::getGradient() - invalid gradient" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return gradVector;
}

const Abstract::Vector&
Group::getNewton() const
{
  if (!isNewton()) {
    std::cerr << "ERROR: NOX::Petsc::Group::getNewton() - invalid Newton vector" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return NewtonVector;
}

Teuchos::RCP< const Abstract::Vector >
Group::getXPtr() const
{
  return Teuchos::rcp< const NOX::Abstract::Vector >(&xVector, false);
}

Teuchos::RCP< const Abstract::Vector >
Group::getFPtr() const
{
  if (!isF()) {
    std::cerr << "ERROR: NOX::Petsc::Group::getFPtr() - invalid RHS" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return Teuchos::rcp< const NOX::Abstract::Vector >(&RHSVector, false);
}

Teuchos::RCP< const Abstract::Vector >
Group::getGradientPtr() const
{
  if (!isGradient()) {
    std::cerr << "ERROR: NOX::Petsc::Group::getGradientPtr() - invalid gradient" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return Teuchos::rcp< const NOX::Abstract::Vector >(&gradVector, false);
}

Teuchos::RCP< const Abstract::Vector >
Group::getNewtonPtr() const
{
  if (!isNewton()) {
    std::cerr << "ERROR: NOX::Petsc::Group::getNewtonPtr() - invalid Newton vector" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return Teuchos::rcp< const NOX::Abstract::Vector >(&NewtonVector, false);
}


