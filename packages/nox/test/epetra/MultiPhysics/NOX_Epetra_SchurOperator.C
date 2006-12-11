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
                                                                                
#include "NOX_Epetra_SchurOperator.H"
#include "Problem_Manager.H"

using namespace NOX;
using namespace NOX::Epetra;

SchurOp::SchurOp( int probId_, int depId_, SchurInterface & interface_) :
  schurInterface(interface_),
  probId(probId_),
  depId(depId_),
  label("Schur-Coupling Operator")
{
}

//-----------------------------------------------------------------------------

SchurOp::~SchurOp()
{
}

//-----------------------------------------------------------------------------

int
SchurOp::SetUseTranspose( bool UseTranspose )
{
  if (UseTranspose == true) 
  {
    string msg = "ERROR: NOX::Epetra::SchurOp::SetUseTranspose() - Transpose is unavailable for this operator!";
    throw msg;
  }
  return (-1);
}

//-----------------------------------------------------------------------------

int
SchurOp::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  //Problem_Manager & problemManager =  dynamic_cast<Problem_Manager&>(schurInterface);

  // Convert X and Y from an Epetra_MultiVectors to Epetra_Vectors
  // and NOX::Epetra::Vectors.  This is done so we use a consistent
  // vector space for norms and inner products.
  Teuchos::RefCountPtr<Epetra_Vector> wrappedX = Teuchos::rcp(new Epetra_Vector(View, X, 0));
  Teuchos::RefCountPtr<Epetra_Vector> wrappedY = Teuchos::rcp(new Epetra_Vector(View, Y, 0));
  Teuchos::RefCountPtr<Epetra_Vector> tempX    = Teuchos::rcp(new Epetra_Vector(*wrappedX) );
  Teuchos::RefCountPtr<Epetra_Vector> tempY    = Teuchos::rcp(new Epetra_Vector(*wrappedY) );

  // The substance of this operator ----- to do RWH 10/24/2006

  //cout << "Incoming X :\n" << *wrappedX << endl;
  //problemManager.applyBlockAction( depId, probId, *wrappedX, *tempY);
  schurInterface.applyBlockAction( depId, probId, *wrappedX, *tempY );
  //cout << "After applyBlockAction(" << depId << ", " << probId << ") :\n" << *tempY << endl;
  //problemManager.getBlockInverseOperator(depId)->ApplyInverse(*tempY, *tempY);
  schurInterface.applyBlockInverseAction( depId, depId, *tempY, *tempY );
  //cout << "After ApplyInverse :\n" << *tempY << endl;
  //problemManager.applyBlockAction( probId, depId, *tempY, *tempX);
  schurInterface.applyBlockAction( probId, depId, *tempY, *tempX );
  //cout << "After applyBlockAction(" << probId << ", " << depId << ") :\n" << *tempX << endl;

  //problemManager.getBlockJacobianMatrix(probId)->Apply(*wrappedX, *wrappedY);
  schurInterface.applyBlockAction( probId, probId, *wrappedX, *wrappedY );
  //cout << "After Apply of diagonal block :\n" << *wrappedY << endl;
  wrappedY->Update(-1.0, *tempX, 1.0);
  //cout << "After combingin; final result :\n" << *wrappedY << endl;

  return (0);
}

//-----------------------------------------------------------------------------

int
SchurOp::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  cout << "ERROR: NOX::Epetra::SchurOp::ApplyInverse() - Not valid "
       << "for this operator!" << endl;
  throw "NOX Error";

  return (-1);
}

//-----------------------------------------------------------------------------

double 
SchurOp::NormInf() const
{
  cout << "ERROR: NOX::Epetra::SchurOp::NormInf() - Not Available for "
       << "this operator!" << endl;
  throw "NOX Error";

  return 1.0;
}

//-----------------------------------------------------------------------------

const char* 
SchurOp::Label () const
{
  return label.c_str();
}

//-----------------------------------------------------------------------------

bool 
SchurOp::UseTranspose() const
{
  return false;
}

//-----------------------------------------------------------------------------

bool 
SchurOp::HasNormInf() const
{
  return false;
}

//-----------------------------------------------------------------------------

const Epetra_Comm & 
SchurOp::Comm() const
{

  return schurInterface.getComm();
}

//-----------------------------------------------------------------------------

const Epetra_Map &
SchurOp::OperatorDomainMap() const
{
  //Problem_Manager & problemManager =  dynamic_cast<Problem_Manager&>(schurInterface);

  //return problemManager.getBlockJacobianMatrix(probId)->OperatorDomainMap();

  if( !schurInterface.hasExplicitOperator( probId, probId ) )
    throw "ERROR: NOX::Epetra::SchurOp::OperatorDomainMap() - No available explicit operator for this block.";
  
  return schurInterface.getExplicitOperator(probId, probId)->OperatorDomainMap();
}

//-----------------------------------------------------------------------------

const Epetra_Map & 
SchurOp::OperatorRangeMap() const
{
  //Problem_Manager & problemManager =  dynamic_cast<Problem_Manager&>(schurInterface);

  //return problemManager.getBlockJacobianMatrix(probId)->OperatorRangeMap();
  if( !schurInterface.hasExplicitOperator( probId, probId ) )
    throw "ERROR: NOX::Epetra::SchurOp::OperatorRangeMap() - No available explicit operator for this block.";

  return schurInterface.getExplicitOperator(probId, probId)->OperatorRangeMap();
}

//-----------------------------------------------------------------------------

void
SchurOp::modifyRHS( Teuchos::RefCountPtr<Epetra_Vector> & probVec, Teuchos::RefCountPtr<Epetra_Vector> & depVec )
{
  //Problem_Manager & problemManager =  dynamic_cast<Problem_Manager&>(schurInterface);

  Teuchos::RefCountPtr<Epetra_Vector> tempX = Teuchos::rcp( new Epetra_Vector(*probVec) );
  Teuchos::RefCountPtr<Epetra_Vector> tempY = Teuchos::rcp( new Epetra_Vector(*depVec ) );

  //problemManager.getBlockInverseOperator(depId)->ApplyInverse(*depVec, *tempY);
  //problemManager.applyBlockAction( probId, depId, *tempY, *tempX);
  schurInterface.applyBlockInverseAction( depId, depId, *depVec, *tempY );
  schurInterface.applyBlockAction( probId, depId, *tempY, *tempX);

  probVec->Update(-1.0, *tempX, 1.0);

  return;
}

//-----------------------------------------------------------------------------

