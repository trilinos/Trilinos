
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#include "NOX_Epetra_Scaling.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_LinearProblem.h"

#include "NOX_Utils.H"

NOX::Epetra::Scaling::Scaling() :
  tmpVectorPtr(0)
{
  
}

NOX::Epetra::Scaling::~Scaling()
{
  delete tmpVectorPtr;
}

void NOX::Epetra::Scaling::addUserScaling(ScaleType type, Epetra_Vector& D)
{
  scaleType.push_back(type);
  sourceType.push_back(UserDefined);
  scaleVector.push_back(&D);
}

void NOX::Epetra::Scaling::addRowSumScaling(ScaleType type, Epetra_Vector& D)
{
  scaleType.push_back(type);
  sourceType.push_back(RowSum);
  scaleVector.push_back(&D);
}

void NOX::Epetra::Scaling::computeScaling(const Epetra_LinearProblem& problem)
{
  
  Epetra_Vector* diagonal = 0;
  for (unsigned int i = 0; i < scaleVector.size(); i ++) {
 
    if (sourceType[i] == RowSum) {
      
      diagonal = scaleVector[i];

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
  }

}

void NOX::Epetra::Scaling::scaleLinearSystem(Epetra_LinearProblem& problem)
{
  Epetra_Vector* diagonal = 0;
  for (unsigned int i = 0; i < scaleVector.size(); i ++) {
 
    diagonal = scaleVector[i];

    if (scaleType[i] == Left) {
 
      if (tmpVectorPtr == 0)
	tmpVectorPtr = new Epetra_Vector(*diagonal);
     
      tmpVectorPtr->Reciprocal(*diagonal);
      problem.LeftScale(*tmpVectorPtr);

    }
    else if (scaleType[i] == Right)
      problem.RightScale(*diagonal);

  }
}

void NOX::Epetra::Scaling::unscaleLinearSystem(Epetra_LinearProblem& problem)
{
  Epetra_Vector* diagonal = 0;
  for (unsigned int i = 0; i < scaleVector.size(); i ++) {
    
    diagonal = scaleVector[i];
    
    if (scaleType[i] == Left) { 
      problem.LeftScale(*diagonal);
    }
    else if (scaleType[i] == Right) {

      if (tmpVectorPtr == 0)
	tmpVectorPtr = new Epetra_Vector(*diagonal);

      tmpVectorPtr->Reciprocal(*diagonal);
      problem.RightScale(*tmpVectorPtr);
      
    }
  }
}

void NOX::Epetra::Scaling::applyRightScaling(const Epetra_Vector& input, 
					     Epetra_Vector& result)
{
  Epetra_Vector* diagonal = 0;
  for (unsigned int i = 0; i < scaleVector.size(); i ++) {

    if (scaleType[i] == Right) {
      diagonal = scaleVector[i];

      if (tmpVectorPtr == 0)
	tmpVectorPtr = new Epetra_Vector(*diagonal);
     
      tmpVectorPtr->Reciprocal(*diagonal);

      result.Multiply(1.0, input, *tmpVectorPtr, 0.0);
    }
  }
}

void NOX::Epetra::Scaling::applyLeftScaling(const Epetra_Vector& input, 
					    Epetra_Vector& result)
{ 
  Epetra_Vector* diagonal = 0;
  for (unsigned int i = 0; i < scaleVector.size(); i ++) {
    
    if (scaleType[i] == Left) {
      diagonal = scaleVector[i];

      if (tmpVectorPtr == 0)
	tmpVectorPtr = new Epetra_Vector(*diagonal);
     
      tmpVectorPtr->Reciprocal(*diagonal);

      result.Multiply(1.0, input, *tmpVectorPtr, 0.0);
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
    else
      source = "Row Sum";

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
