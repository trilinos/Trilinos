//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
                                                                                
#include "NOX_Epetra_SchurCouplingOperator.H"

using namespace NOX;
using namespace NOX::Epetra;

SchurCouplingOp::SchurCouplingOp( int probId_, int depId_, Problem_Manager & probManager) :
  problemManager(probManager),
  probId(probId_),
  depId(depId_),
  label("Schur Coupling Operator")
{
}

SchurCouplingOp::~SchurCouplingOp()
{
}

int
SchurCouplingOp::SetUseTranspose( bool UseTranspose )
{
  if (UseTranspose == true) 
  {
    cout << "ERROR: NOX::Epetra::SchurCouplingOp::SetUseTranspose() - Transpose is "
	 << "unavailable for this operator!" << endl;
    throw "NOX Error";
  }
  return (-1);
}

int
SchurCouplingOp::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // Convert X and Y from an Epetra_MultiVector to a Epetra_Vectors
  // and NOX::epetra::Vectors.  This is done so we use a consistent
  // vector space for norms and inner products.
  Teuchos::RefCountPtr<Epetra_Vector> wrappedX = Teuchos::rcp(new Epetra_Vector(View, X, 0));
  Teuchos::RefCountPtr<Epetra_Vector> wrappedY = Teuchos::rcp(new Epetra_Vector(View, Y, 0));
  Teuchos::RefCountPtr<Epetra_Vector> tempX    = Teuchos::rcp(new Epetra_Vector(*wrappedX) );
  Teuchos::RefCountPtr<Epetra_Vector> tempY    = Teuchos::rcp(new Epetra_Vector(problemManager.getSolutionVec(depId)) );

  // The substance of this operator ----- to do RWH 10/24/2006

  cout << "Incoming X :\n" << *wrappedX << endl;
  problemManager.applyBlockAction( depId, probId, *wrappedX, *tempY);
  cout << "After applyBlockAction(" << depId << ", " << probId << ") :\n" << *tempY << endl;
  problemManager.getBlockInverseOperator(depId)->ApplyInverse(*tempY, *tempY);
  cout << "After ApplyInverse :\n" << *tempY << endl;
  problemManager.applyBlockAction( probId, depId, *tempY, *tempX);
  cout << "After applyBlockAction(" << probId << ", " << depId << ") :\n" << *tempX << endl;

  problemManager.getBlockJacobianMatrix(probId)->Apply(*wrappedX, *wrappedY);
  cout << "After Apply of diagonal block :\n" << *wrappedY << endl;
  wrappedY->Update(-1.0, *tempX, 1.0);
  cout << "After combingin; final result :\n" << *wrappedY << endl;

  //Teuchos::RefCountPtr<Epetra_Vector> tmpA = Teuchos::rcp( new Epetra_Vector( *problemManager.getProblem(1).getSolution()) );
  //Teuchos::RefCountPtr<Epetra_Vector> tmpB = Teuchos::rcp( new Epetra_Vector( *problemManager.getProblem(1).getSolution()) );

  //Teuchos::ParameterList& lsParams = (*problemManager.nlParams).sublist("Direction").sublist("Newton").sublist("Linear Solver");

  //problemManager.setGroupX(1);
  //problemManager.setGroupX(2);
  //problemManager.getGroup(1).computeJacobian();
  //problemManager.getGroup(2).computeJacobian();
  //problemManager.createBlockInverseOperator(1, lsParams);
  //problemManager.createBlockInverseOperator(2, lsParams);

  ////cout << "\nBefore doing apply(1,1) :\n" << *tmpA << *tmpB << endl;
  ////problemManager.getBlockJacobianMatrix(1,1)->Apply(*tmpA, *tmpB);
  ////cout << "\nAfter doing apply(1,1) :\n" << *tmpA << *tmpB << endl;

  ////problemManager.getBlockInverseOperator(1)->ApplyInverse(*tmpB, *tmpA);
  ////cout << "\nAfter doing applyinverse(1) :\n" << *tmpA << *tmpB << endl;

  //cout << "Problem 1 solution :\n" << problemManager.getSolutionVec(1) << endl;
  //cout << "Problem 2 solution :\n" << problemManager.getSolutionVec(2) << endl;

  //problemManager.getGroup(1).computeF();
  //problemManager.getGroup(2).computeF();

  //cout << "Problem 1 residual :\n" << problemManager.getResidualVec(1) << endl;
  //cout << "Problem 2 residual :\n" << problemManager.getResidualVec(2) << endl;

  //*tmpB = problemManager.getResidualVec(2);
  //problemManager.getBlockInverseOperator(2)->ApplyInverse(*tmpB, *tmpB);
  //cout << "After applying blockInverse(2) to Problem 2 residual :\n" << *tmpB << endl;

  //problemManager.applyBlockAction( 1, 2, *tmpB, *tmpA);
  //cout << "After applying blockAction(1, 2) to modified Problem 2 residual :\n" << *tmpA << endl;

  return (0);
}

int
SchurCouplingOp::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  cout << "ERROR: NOX::Epetra::SchurCouplingOp::ApplyInverse() - Not valid "
       << "for this operator!" << endl;
  throw "NOX Error";

  return (-1);
}

double 
SchurCouplingOp::NormInf() const
{
  cout << "ERROR: NOX::Epetra::SchurCouplingOp::NormInf() - Not Available for "
       << "this operator!" << endl;
  throw "NOX Error";

  return 1.0;
}

const char* 
SchurCouplingOp::Label () const
{
  return label.c_str();
}

bool 
SchurCouplingOp::UseTranspose() const
{
  return false;
}

bool 
SchurCouplingOp::HasNormInf() const
{
  return false;
}

const Epetra_Comm & 
SchurCouplingOp::Comm() const
{
  return problemManager.getJacobianOperator()->Comm();
}

const Epetra_Map &
SchurCouplingOp::OperatorDomainMap() const
{
  return problemManager.getBlockJacobianMatrix(probId)->OperatorDomainMap();
}

const Epetra_Map & 
SchurCouplingOp::OperatorRangeMap() const
{
  return problemManager.getBlockJacobianMatrix(probId)->OperatorRangeMap();
}

void
SchurCouplingOp::modifyRHS( Teuchos::RefCountPtr<Epetra_Vector> & probVec, Teuchos::RefCountPtr<Epetra_Vector> & depVec )
{
  Teuchos::RefCountPtr<Epetra_Vector> tempX = Teuchos::rcp( new Epetra_Vector(*probVec) );
  Teuchos::RefCountPtr<Epetra_Vector> tempY = Teuchos::rcp( new Epetra_Vector(*depVec ) );

  problemManager.getBlockInverseOperator(depId)->ApplyInverse(*depVec, *tempY);
  problemManager.applyBlockAction( probId, depId, *tempY, *tempX);

  probVec->Update(-1.0, *tempX, 1.0);

  return;
}
