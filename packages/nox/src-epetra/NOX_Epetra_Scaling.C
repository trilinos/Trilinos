//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

void NOX::Epetra::Scaling::addUserScaling(ScaleType type, const Teuchos::RCP<Epetra_Vector>& D)
{
  if ( Teuchos::is_null(tmpVectorPtr) ) 
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(*D));

  scaleType.push_back(type);
  sourceType.push_back(UserDefined);
  scaleVector.push_back(D);
}

void NOX::Epetra::Scaling::addRowSumScaling(ScaleType type, const Teuchos::RCP<Epetra_Vector>& D)
{
  if ( Teuchos::is_null(tmpVectorPtr) )
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(*D));

  scaleType.push_back(type);
  sourceType.push_back(RowSum);
  scaleVector.push_back(D);
}

void NOX::Epetra::Scaling::addColSumScaling(ScaleType type, const Teuchos::RCP<Epetra_Vector>& D)
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
	     << "Epetra_RowMatrix derived object!" << std::endl;
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
	     << "Epetra_RowMatrix derived object!" << std::endl;
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

  os << "\n       LINEAR SOLVER SCALING:" << std::endl;

  for (unsigned int i = 0; i < scaleVector.size(); i ++) {

    std::string source = " ";
    if (sourceType[i] == UserDefined)
      source = "User Defined Vector";
    else if (sourceType[i] == RowSum)
      source = "Row Sum";
    else if (sourceType[i] == ColSum)
      source = "Col Sum";

    if (scaleType[i] == Left) {
      os << "       " << (i+1) << ".  Left Scaled with " << source << std::endl;

    }
    else if (scaleType[i] == Right)
      os << "       " << (i+1) << ".  Right Scaled with " << source << std::endl;
  }

  return;
}

ostream& 
NOX::Epetra::operator<<(ostream& os, NOX::Epetra::Scaling& scalingObject)
{
  scalingObject.print(os);
  return os;
}
