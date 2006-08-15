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

#include "NOX_Epetra_Scaling.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_LinearProblem.h"

#include "NOX_Utils.H"

NOX::Epetra::Scaling::Scaling()
{
  
}

NOX::Epetra::Scaling::~Scaling()
{

}

void NOX::Epetra::Scaling::addUserScaling(ScaleType type, const Teuchos::RefCountPtr<Epetra_Vector>& D)
{
  if ( Teuchos::is_null(tmpVectorPtr) ) 
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(*D));

  scaleType.push_back(type);
  sourceType.push_back(UserDefined);
  scaleVector.push_back(D);
}

void NOX::Epetra::Scaling::addRowSumScaling(ScaleType type, const Teuchos::RefCountPtr<Epetra_Vector>& D)
{
  if ( Teuchos::is_null(tmpVectorPtr) )
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(*D));

  scaleType.push_back(type);
  sourceType.push_back(RowSum);
  scaleVector.push_back(D);
}

void NOX::Epetra::Scaling::addColSumScaling(ScaleType type, const Teuchos::RefCountPtr<Epetra_Vector>& D)
{
  if ( Teuchos::is_null(tmpVectorPtr) )
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(*D));

  scaleType.push_back(type);
  sourceType.push_back(ColSum);
  scaleVector.push_back(D);
}

void NOX::Epetra::Scaling::computeScaling(const Epetra_LinearProblem& problem)
{
  
  Epetra_Vector* diagonal = 0;
  for (unsigned int i = 0; i < scaleVector.size(); i ++) {
 
    if (sourceType[i] == RowSum) {
      
      diagonal = scaleVector[i].get();

      // Make sure the Jacobian is an Epetra_RowMatrix, otherwise we can't 
      // perform a row sum scale!
      const Epetra_RowMatrix* test = 0;
      test = dynamic_cast<const Epetra_RowMatrix*>(problem.GetOperator());
      if (test == 0) {
	cout << "ERROR: NOX::Epetra::Scaling::scaleLinearSystem() - "
	     << "For \"Row Sum\" scaling, the Matrix must be an "
	     << "Epetra_RowMatrix derived object!" << endl;
	throw "NOX Error";
      }
      
      test->InvRowSums(*diagonal);
      diagonal->Reciprocal(*diagonal);

    }

    else if (sourceType[i] == ColSum) {
      
      diagonal = scaleVector[i].get();

      // Make sure the Jacobian is an Epetra_RowMatrix, otherwise we can't 
      // perform a row sum scale!
      const Epetra_RowMatrix* test = 0;
      test = dynamic_cast<const Epetra_RowMatrix*>(problem.GetOperator());
      if (test == 0) {
	cout << "ERROR: NOX::Epetra::Scaling::scaleLinearSystem() - "
	     << "For \"Column Sum\" scaling, the Matrix must be an "
	     << "Epetra_RowMatrix derived object!" << endl;
	throw "NOX Error";
      }
      
      test->InvColSums(*diagonal);
      diagonal->Reciprocal(*diagonal);

    }

  }

}

void NOX::Epetra::Scaling::scaleLinearSystem(Epetra_LinearProblem& problem)
{
  Epetra_Vector* diagonal = 0;
  for (unsigned int i = 0; i < scaleVector.size(); i ++) {
 
    diagonal = scaleVector[i].get();

    if (scaleType[i] == Left) {
 
      tmpVectorPtr->Reciprocal(*diagonal);
      problem.LeftScale(*tmpVectorPtr);

    }
    else if (scaleType[i] == Right) {

      tmpVectorPtr->Reciprocal(*diagonal);
      problem.RightScale(*tmpVectorPtr);
    }

  }
}

void NOX::Epetra::Scaling::unscaleLinearSystem(Epetra_LinearProblem& problem)
{
  Epetra_Vector* diagonal = 0;
  for (unsigned int i = 0; i < scaleVector.size(); i ++) {
    
    diagonal = scaleVector[i].get();
    
    if (scaleType[i] == Left) { 
      problem.LeftScale(*diagonal);
    }
    else if (scaleType[i] == Right) {
      problem.RightScale(*diagonal);
      
    }
  }
}

void NOX::Epetra::Scaling::applyRightScaling(const Epetra_Vector& input, 
					     Epetra_Vector& result)
{
  if (scaleVector.size() == 0) {
    result = input;
  }
  else {
    Epetra_Vector* diagonal = 0;
    for (unsigned int i = 0; i < scaleVector.size(); i ++) {
      
      if (scaleType[i] == Right) {
	diagonal = scaleVector[i].get();
	
	tmpVectorPtr->Reciprocal(*diagonal);
	
	result.Multiply(1.0, input, *tmpVectorPtr, 0.0);
      }
    }
  }
}

void NOX::Epetra::Scaling::applyLeftScaling(const Epetra_Vector& input, 
					    Epetra_Vector& result)
{ 
  if (scaleVector.size() == 0) {
    result = input;
  }
  else {
    Epetra_Vector* diagonal = 0;
    for (unsigned int i = 0; i < scaleVector.size(); i ++) {
      
      if (scaleType[i] == Left) {
	diagonal = scaleVector[i].get();
	
	tmpVectorPtr->Reciprocal(*diagonal);
	
	result.Multiply(1.0, input, *tmpVectorPtr, 0.0);
      }
    }
  }
}

void NOX::Epetra::Scaling::print(ostream& os)
{

  os << "\n       LINEAR SOLVER SCALING:" << endl;

  for (unsigned int i = 0; i < scaleVector.size(); i ++) {

    string source = " ";
    if (sourceType[i] == UserDefined)
      source = "User Defined Vector";
    else if (sourceType[i] == RowSum)
      source = "Row Sum";
    else if (sourceType[i] == ColSum)
      source = "Col Sum";

    if (scaleType[i] == Left) {
      os << "       " << (i+1) << ".  Left Scaled with " << source << endl;

    }
    else if (scaleType[i] == Right)
      os << "       " << (i+1) << ".  Right Scaled with " << source << endl;
  }

  return;
}

ostream& operator<<(ostream& os, NOX::Epetra::Scaling& scalingObject)
{
  scalingObject.print(os);
  return os;
}
