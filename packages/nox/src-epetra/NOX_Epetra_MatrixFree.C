// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_CompilerCodeTweakMacros.hpp"
#include "NOX_Abstract_Group.H"
#include "NOX_Epetra_Vector.H"
#include "NOX_Epetra_VectorSpace.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Utils.H"

#include "NOX_Epetra_MatrixFree.H"

using namespace NOX;
using namespace NOX::Epetra;

MatrixFree::MatrixFree(Teuchos::ParameterList& printParams,
               const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,
               const NOX::Epetra::Vector& x, bool p) :
  label("NOX::Matrix-Free"),
  interface(i),
  currentX(x),
  perturbX(x),
  fo(x),
  fp(x),
  diffType(Forward),
  lambda(1.0e-6),
  eta(0.0),
  userEta(1.0e-6),
  computeEta(true),
  useGroupForComputeF(false),
  useSolverForComputeJacobian(false),
  useNewPerturbation(p),
  utils(printParams)
{
  // Zero out Vectors
  perturbX.init(0.0);
  fo.init(0.0);
  fp.init(0.0);

  // Epetra_Operators require Epetra_Maps, so anyone using block maps
  // (Epetra_BlockMap) won't be able to directly use the AztecOO solver.
  // We get around this by creating an Epetra_Map from the Epetra_BlockMap.
  const Epetra_Map* testMap = 0;
  testMap = dynamic_cast<const Epetra_Map*>(&currentX.getEpetraVector().Map());
  if (testMap != 0) {
    epetraMap = Teuchos::rcp(new Epetra_Map(*testMap));
  }
  else {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    int size = currentX.getEpetraVector().Map().NumGlobalPoints();
    int indexBase = currentX.getEpetraVector().Map().IndexBase();
#else
    long long size = currentX.getEpetraVector().Map().NumGlobalPoints64();
    long long indexBase = currentX.getEpetraVector().Map().IndexBase64();
#endif
    int mySize = currentX.getEpetraVector().Map().NumMyPoints();
    const Epetra_Comm& comm = currentX.getEpetraVector().Map().Comm();
    epetraMap = Teuchos::rcp(new Epetra_Map(size, mySize, indexBase, comm));
  }

}

MatrixFree::~MatrixFree()
{

}

int MatrixFree::SetUseTranspose(bool use_transpose)
{
  if (use_transpose == true) {
    utils.out() << "ERROR: NOX::Epetra::MatrixFree::SetUseTranspose() - Transpose is "
     << "unavailable in Matrix-Free mode!" << std::endl;
    throw std::runtime_error("NOX Error");
  }
  return (-1);
}

int MatrixFree::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  // Use a directional derivative to compute y = Jx
  /*
   * eta = scalar perturbation
   * u = solution vector used to evaluate f
   * f = function evaluation (RHS)
   * x = vector that J is applied to
   *
   *        f(u+eta*x) - f(u)
   * Jx =   -----------------
   *               eta
   */

  // Convert X and Y from an Epetra_MultiVector to a Epetra_Vectors
  // and NOX::epetra::Vectors.  This is done so we use a consistent
  // vector space for norms and inner products.
  Teuchos::RCP<Epetra_Vector> wrappedX =
    Teuchos::rcp(new Epetra_Vector(View, X, 0));
  Teuchos::RCP<Epetra_Vector> wrappedY =
    Teuchos::rcp(new Epetra_Vector(View, Y, 0));
  NOX::Epetra::Vector nevX(wrappedX, NOX::Epetra::Vector::CreateView);
  NOX::Epetra::Vector nevY(wrappedY, NOX::Epetra::Vector::CreateView);

  // Compute perturbation constant, eta
  // Taken from LOCA v1.0 manual SAND2002-0396 p. 28 eqn. 2.43
  // eta = lambda*(lambda + 2norm(u)/2norm(x))
  double solutionNorm = 1.0;
  double vectorNorm = 1.0;

  solutionNorm = currentX.norm();
  vectorNorm = currentX.getVectorSpace()->norm(*wrappedX);

  // Make sure the norm is not zero, otherwise we can get an inf perturbation
  if (vectorNorm == 0.0) {
    //utils.out(Utils::Warning) << "Warning: NOX::Epetra::MatrixFree::Apply() - vectorNorm is zero" << std::endl;
    vectorNorm = 1.0;
    wrappedY->PutScalar(0.0);
    return 0;
  }

  // Create an extra perturbed residual vector pointer if needed
  if ( diffType == Centered )
    if ( Teuchos::is_null(fmPtr) )
      fmPtr = Teuchos::rcp(new NOX::Epetra::Vector(fo));

  const double scaleFactor = (diffType == Backward ? -1.0 : 1.0);

  if (computeEta) {
    if (useNewPerturbation) {
      double dotprod = currentX.getVectorSpace()->
    innerProduct(currentX.getEpetraVector(), *wrappedX);
      if (dotprod==0.0)
    dotprod = 1.0e-12;
      eta = lambda*(1.0e-12/lambda + fabs(dotprod)/(vectorNorm * vectorNorm))
    * dotprod/fabs(dotprod);
    }
    else
      eta = lambda*(lambda + solutionNorm/vectorNorm);
  }
  else
    eta = userEta;

  // Compute the perturbed RHS
  perturbX.update(1.0,currentX,eta*scaleFactor,nevX,0.0);

  if (!useGroupForComputeF)
      interface->computeF(perturbX.getEpetraVector(), fp.getEpetraVector(),
              NOX::Epetra::Interface::Required::MF_Res);
  else{
    groupPtr->setX(perturbX);
    groupPtr->computeF();
    fp = dynamic_cast<const NOX::Epetra::Vector&>
      (groupPtr->getF());
  }

  if ( diffType == Centered ) {
    perturbX.update(-2.0*eta*scaleFactor,nevX,1.0);
    if (!useGroupForComputeF)
      interface->computeF(perturbX.getEpetraVector(), fmPtr->getEpetraVector(),
              NOX::Epetra::Interface::Required::MF_Res);
    else{
      groupPtr->setX(perturbX);
      groupPtr->computeF();
      *fmPtr = dynamic_cast<const NOX::Epetra::Vector&>
        (groupPtr->getF());
    }
  }

  // Compute the directional derivative
  if ( diffType != Centered ) {
    nevY.update(1.0, fp, -1.0, fo, 0.0);
    nevY.scale( 1.0/(scaleFactor * eta) );
  }
  else {
    nevY.update(1.0, fp, -1.0, *fmPtr, 0.0);
    nevY.scale( 1.0/(2.0 * eta * scaleFactor) );
  }

  return 0;
}

int MatrixFree::ApplyInverse(const Epetra_MultiVector& /* X */, Epetra_MultiVector& /* Y */) const
{
  utils.out() << "ERROR: NOX::MatrixFree::ApplyInverse - Not available for Matrix Free!"
       << std::endl;
  throw std::runtime_error("NOX Error");
  TEUCHOS_UNREACHABLE_RETURN(-1);
}

double MatrixFree::NormInf() const
{
  utils.out() << "ERROR: NOX::Epetra::MatrixFree::NormInf() - Not Available for "
       << "Matrix-Free mode!" << std::endl;
  throw std::runtime_error("NOX Error");
  TEUCHOS_UNREACHABLE_RETURN(1.0);
}


const char* MatrixFree::Label () const
{
  return label.c_str();
}

bool MatrixFree::UseTranspose() const
{
  return false;
}

bool MatrixFree::HasNormInf() const
{
  return false;
}

const Epetra_Comm & MatrixFree::Comm() const
{
  return currentX.getEpetraVector().Map().Comm();
}
const Epetra_Map& MatrixFree::OperatorDomainMap() const
{
  return *epetraMap;
}

const Epetra_Map& MatrixFree::OperatorRangeMap() const
{
  return *epetraMap;
}

bool MatrixFree::computeJacobian(const Epetra_Vector& x, Epetra_Operator& /* Jac */)
{
  // Since we have no explicit Jacobian we set our currentX to the
  // incoming value and evaluate the RHS.  When the Jacobian is applied,
  // we compute the perturbed residuals and the directional
  // derivative.  Ignore Jac.
  currentX = x;

  bool ok = false;
  if (!useSolverForComputeJacobian) {
    if (!useGroupForComputeF)
      ok = interface->computeF(x, fo.getEpetraVector(),
                   NOX::Epetra::Interface::Required::MF_Jac);
    else {
      groupPtr->setX(currentX);
      groupPtr->computeF();
      fo = dynamic_cast<const NOX::Epetra::Vector&>
        (groupPtr->getF());
      ok = true;
    }
  }
  else {
    // If this is throwing an invalid resid, switch to
    // Full Step (no line search) so last Resid is current.
    fo  = dynamic_cast<const NOX::Epetra::Vector&>
        (slvrPtr->getSolutionGroup().getF());
    ok = true;
  }
  return ok;
}

void MatrixFree::setDifferenceMethod(DifferenceType diffType_)
{
  diffType = diffType_;
}

void MatrixFree::setLambda(double lambda_)
{
  lambda = lambda_;
}

void MatrixFree::setComputePerturbation(bool bVal)
{
  computeEta = bVal;
}

void MatrixFree::setPerturbation(double eta_)
{
  userEta = eta_;
  computeEta = false;
}

double MatrixFree::getPerturbation() const
{
  return eta;
}

void MatrixFree::setGroupForComputeF(const NOX::Abstract::Group& group)
{
  useGroupForComputeF = true;
  groupPtr = group.clone();
  return;
}

void MatrixFree::setSolverForComputeJacobian(const Teuchos::RCP<NOX::Solver::Generic>& slvr)
{
  useSolverForComputeJacobian = true;
  slvrPtr = slvr.create_weak();
  return;
}
