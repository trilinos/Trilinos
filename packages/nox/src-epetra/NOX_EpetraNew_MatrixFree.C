//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER
                                                                                
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "NOX_Abstract_Group.H"
#include "NOX_Epetra_Vector.H"
#include "NOX_EpetraNew_Interface_Required.H"
#include "NOX_Utils.H"

#include "NOX_EpetraNew_MatrixFree.H"

using namespace NOX;
using namespace NOX::EpetraNew;

MatrixFree::MatrixFree(Interface::Required& i, const Epetra_Vector& x) :
  label("NOX::Matrix-Free"),
  interface(i),
  currentX(x),
  perturbX(x),
  fo(x),
  fp(x),
  fmPtr(0),
  epetraMap(0),
  ownsMap(false),
  lambda(1.0e-6),
  diffType(Forward),
  eta(0.0),
  userEta(1.0e-6),
  computeEta(true),
  useGroupForComputeF(false),
  groupPtr(0)
{
  // Zero out Vectors
  perturbX.PutScalar(0.0);
  fo.PutScalar(0.0);
  fp.PutScalar(0.0);

  // Create an EpetraMap if needed (First try and get it from the 
  // solution vector currentX)
  const Epetra_BlockMap* bmap = 0;
  bmap = dynamic_cast<const Epetra_Map*>(&currentX.Map());
  if (bmap != 0)
    epetraMap = dynamic_cast<const Epetra_Map*>(bmap);
  else {
    
    int size = currentX.Map().NumGlobalElements();
    int mySize = currentX.Map().NumMyElements();
    int indexBase = currentX.Map().IndexBase();
    int* globalElementList = currentX.Map().MyGlobalElements();
    const Epetra_Comm& comm = currentX.Map().Comm();
    epetraMap = new Epetra_Map(size, mySize, globalElementList, indexBase, comm);
    ownsMap = true;
  }

}

MatrixFree::~MatrixFree()
{
  delete groupPtr;
  if (ownsMap)
    delete epetraMap;
}

int MatrixFree::SetUseTranspose(bool UseTranspose) 
{
  if (UseTranspose == true) {
    cout << "ERROR: NOX::EpetraNew::MatrixFree::SetUseTranspose() - Transpose is "
	 << "unavailable in Matrix-Free mode!" << endl;
    throw "NOX Error";
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

  // Compute perturbation constant, eta
  // Taken from LOCA v1.0 manual SAND2002-0396 p. 28 eqn. 2.43
  // eta = lambda*(lambda + 2norm(u)/2norm(x))
  double solutionNorm = 1.0;
  double vectorNorm = 1.0;

  int test = currentX.Norm2(&solutionNorm);

  // Make sure the norm computed correctly
  if (test != 0) {
    if (NOX::Utils::doPrint(Utils::Warning)) 
      cout << "Warning: NOX::EpetraNew::MatrixFree::Apply() - solutionNorm "
	   << "failed!" << endl;
    solutionNorm = 1.0;
  }

  test = X.Norm2(&vectorNorm);

  // Make sure the norm computed correctly
  if (test != 0) {
    if (NOX::Utils::doPrint(Utils::Warning)) 
      cout << "Warning: NOX::EpetraNew::MatrixFree::Apply() - vectorNorm failed!" 
	   << endl;
    vectorNorm = 1.0;
  }

  // Make sure the norm is not zero, otherwise we can get an inf perturbation
  if (vectorNorm == 0.0) {
    //if (NOX::Utils::doPrint(Utils::Warning)) 
    //cout << "Warning: NOX::EpetraNew::MatrixFree::Apply() - vectorNorm is zero" 
    //<< endl;
    vectorNorm = 1.0;
  }

  // Create an extra perturbed residual vector pointer if needed
  if ( diffType == Centered )
    if ( !fmPtr )
      fmPtr = new Epetra_Vector(fo);
  
  // Create a reference to the extra perturbed residual vector
  Epetra_Vector& fm = *fmPtr;
  
  double scaleFactor = 1.0;
  if ( diffType == Backward )
  scaleFactor = -1.0;

  if (computeEta)
    eta = lambda*(lambda + solutionNorm/vectorNorm);
  else
    eta = userEta;

  // Compute the perturbed RHS
  perturbX = currentX;
  Y = X;
  Y.Scale(eta);
  perturbX.Update(1.0,Y,1.0);

  if (!useGroupForComputeF)
    interface.computeF(perturbX, fp, NOX::EpetraNew::Interface::Required::MF_Res);
  else{
    NOX::Epetra::Vector noxX(perturbX, NOX::DeepCopy, true);
    groupPtr->setX(noxX);
    groupPtr->computeF();
    fp = dynamic_cast<const NOX::Epetra::Vector&>
      (groupPtr->getF()).getEpetraVector();
  } 

  if ( diffType == Centered ) {
    Y.Scale(-2.0);
    perturbX.Update(scaleFactor,Y,1.0);
    if (!useGroupForComputeF)
      interface.computeF(perturbX, fm, NOX::EpetraNew::Interface::Required::MF_Res);
    else{
      NOX::Epetra::Vector noxX(perturbX, NOX::DeepCopy, true);
      groupPtr->setX(noxX);
      groupPtr->computeF();
      fm = dynamic_cast<const NOX::Epetra::Vector&>
        (groupPtr->getF()).getEpetraVector();
    } 
  }

  // Compute the directional derivative
  if ( diffType != Centered ) {
    Y.Update(1.0, fp, -1.0, fo, 0.0);
    Y.Scale( 1.0/(scaleFactor * eta) );
  }
  else {
    Y.Update(1.0, fp, -1.0, fm, 0.0);
    Y.Scale( 1.0/(2.0 * eta) );
  }
  
  return 0;
}

int MatrixFree::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  cout << "ERROR: NOX::MatrixFree::ApplyInverse - Not available for Matrix Free!"
       << endl;
  throw "NOX Error";
  return (-1);
}

double MatrixFree::NormInf() const
{
  cout << "ERROR: NOX::EpetraNew::MatrixFree::NormInf() - Not Available for "
       << "Matrix-Free mode!" << endl;
  throw "NOX Error";
  return 1.0;
}


char* MatrixFree::Label () const
{
  return const_cast<char*>(label.c_str());
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
  return currentX.Map().Comm();
}
const Epetra_Map& MatrixFree::OperatorDomainMap() const
{
  return *epetraMap;
}

const Epetra_Map& MatrixFree::OperatorRangeMap() const
{
  return *epetraMap;
}

bool MatrixFree::computeJacobian(const Epetra_Vector& x)
{
  // Since we have no explicit Jacobian we set our currentX to the 
  // incoming value and evaluate the RHS.  When the Jacobian is applied,
  // we compute the perturbed residuals and the directional 
  // derivative.
  currentX = x;

  bool ok = false;
  if (!useGroupForComputeF)
    ok = interface.computeF(x, fo, NOX::EpetraNew::Interface::Required::MF_Res);
  else {
    NOX::Epetra::Vector noxX(currentX, NOX::DeepCopy, true);
    groupPtr->setX(noxX);
    groupPtr->computeF();
    fo = dynamic_cast<const NOX::Epetra::Vector&>
      (groupPtr->getF()).getEpetraVector();
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

void MatrixFree::setGroupForComputeF(NOX::Abstract::Group& group)
{
  useGroupForComputeF = true;
  delete groupPtr;
  groupPtr = 0;
  groupPtr = group.clone();
  return;
}
