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

#include "NOX_Common.H"
#include "NOX_EpetraNew_BroydenOperator.H"

#include "Epetra_Map.h"

using namespace NOX;
using namespace NOX::EpetraNew;

BroydenOperator::BroydenOperator(NOX::Parameter::List & nlParams, 
                                 Epetra_Vector & solnVec,
                                 Epetra_CrsMatrix & mat,
                                 bool verbose_ ) :
  verbose(verbose_),
  updateVectorPtr( new NOX::Epetra::Vector(solnVec)),
  updateVector(*updateVectorPtr),
  broydenVecPtr(0),
  residualVecPtr(0),
  tempVecPtr(0),
  crsMatrix(mat),
  jacIntPtr(0),
  jacMatrixPtr(0),
  precIntPtr(0),
  precMatrixPtr(0),
  myType("Broyden Operator")
{
  // Set ourself as the Pre/Post Operator so we can update our data during
  // runPostIterate.  Note: we need to make provision for storing and calling 
  // any existing Pre/Post Operator.

  nlParams.sublist("Solver Options").setParameter(
                     "User Defined Pre/Post Operator", *this);
}

BroydenOperator::BroydenOperator(NOX::Parameter::List & nlParams, 
                  Epetra_Vector & solnVec,
                  Epetra_CrsMatrix & mat,
		  NOX::EpetraNew::Interface::Jacobian & jacInt,
                  Epetra_CrsMatrix & jacMatrix, 
                  bool verbose_ ) :
  verbose(verbose_),
  updateVectorPtr( new NOX::Epetra::Vector(solnVec)),
  updateVector(*updateVectorPtr),
  broydenVecPtr(0),
  residualVecPtr(0),
  tempVecPtr(0),
  crsMatrix(mat),
  jacIntPtr(&jacInt),
  jacMatrixPtr(&jacMatrix),
  precIntPtr(0),
  precMatrixPtr(0),
  myType("Broyden Operator")
{
  // Set ourself as the Pre/Post Operator so we can update our data during
  // runPostIterate.  Note: we need to make provision for storing and calling 
  // any existing Pre/Post Operator.

  nlParams.sublist("Solver Options").setParameter(
                     "User Defined Pre/Post Operator", *this);
}

BroydenOperator::BroydenOperator(NOX::Parameter::List & nlParams, 
                  Epetra_Vector & solnVec,
                  Epetra_CrsMatrix & mat,
		  NOX::EpetraNew::Interface::Preconditioner & precInt,
                  Epetra_CrsMatrix & precMatrix, 
                  bool verbose_ ) :
  verbose(verbose_),
  updateVectorPtr( new NOX::Epetra::Vector(solnVec)),
  updateVector(*updateVectorPtr),
  broydenVecPtr(0),
  residualVecPtr(0),
  tempVecPtr(0),
  crsMatrix(mat),
  jacIntPtr(0),
  jacMatrixPtr(0),
  precIntPtr(&precInt),
  precMatrixPtr(&precMatrix),
  myType("Broyden Operator")
{
  // Set ourself as the Pre/Post Operator so we can update our data during
  // runPostIterate.  Note: we need to make provision for storing and calling 
  // any existing Pre/Post Operator.

  nlParams.sublist("Solver Options").setParameter(
                     "User Defined Pre/Post Operator", *this);
}

BroydenOperator::BroydenOperator(const BroydenOperator & bOp) :
  verbose(bOp.verbose),
  updateVectorPtr( new NOX::Epetra::Vector(bOp.updateVector) ),
  updateVector(*updateVectorPtr),
  broydenVecPtr(0),
  residualVecPtr(0),
  tempVecPtr(0),
  crsMatrix(bOp.crsMatrix),
  jacIntPtr(bOp.jacIntPtr),
  jacMatrixPtr(bOp.jacMatrixPtr),
  precIntPtr(bOp.precIntPtr),
  precMatrixPtr(bOp.precMatrixPtr),
  myType("Broyden Operator")
{ 
  if( bOp.broydenVecPtr )
    broydenVecPtr = new NOX::Epetra::Vector(*bOp.broydenVecPtr);
  if( bOp.residualVecPtr )
    residualVecPtr = new NOX::Epetra::Vector(*bOp.residualVecPtr);
  if( bOp.tempVecPtr )
    tempVecPtr = new NOX::Epetra::Vector(*bOp.tempVecPtr);
}

BroydenOperator::~BroydenOperator()
{
  delete updateVectorPtr; updateVectorPtr = 0;
  delete broydenVecPtr  ; broydenVecPtr   = 0;
  delete residualVecPtr ; residualVecPtr  = 0;
  delete tempVecPtr     ; tempVecPtr      = 0;
}

bool BroydenOperator::computeJacobian(const Epetra_Vector & x)
{
  if( jacIntPtr ) {
    jacIntPtr->computeJacobian(x);
    replaceBroydenMatrixValues(*jacMatrixPtr);
  }

  return true;
}

bool BroydenOperator::computePreconditioner(const Epetra_Vector & x,
           NOX::Parameter::List * pList)
{
  if( precIntPtr ) {
    precIntPtr->computePreconditioner(x, pList);
    replaceBroydenMatrixValues(*precMatrixPtr);
  }

  return true;
}

void BroydenOperator::runPostIterate( const NOX::Solver::Generic & solver)
{ 
  // Get and update using the solver object.

  NOX::Abstract::Group::ReturnType status;
  int ierr;

  if( solver.getNumIterations() > 0 )
  {
    if( !broydenVecPtr )
      broydenVecPtr = new NOX::Epetra::Vector(updateVector);
    if( !residualVecPtr )
      residualVecPtr = new NOX::Epetra::Vector(updateVector);
    if( verbose )
      if( !tempVecPtr )
        tempVecPtr = new NOX::Epetra::Vector(updateVector);
      
    // Store previous Newton vector as the update, s
    const Abstract::Group & oldSolnGrp = solver.getPreviousSolutionGroup();
    if( !oldSolnGrp.isNewton() )
    {
      cout << "ERROR: NOX::EpetraNew::BroydenOperator::runPostIterate(...) "
           << "- getNewton() failed!!!" << endl;
      throw "NOX Error: Broyden Update Failed";
    }
    updateVector.update(1.0, oldSolnGrp.getNewton(), 0.0);
    if( verbose ) {
      oldSolnGrp.applyJacobian(updateVector, *tempVecPtr);
       cout << "Js vector : " << endl << tempVecPtr->getEpetraVector()
            << "\nOld residual vector : " << endl 
            << dynamic_cast<const NOX::Epetra::Vector&>(oldSolnGrp.getF()).getEpetraVector() << endl;
    }
      

    // Do the Broyden update to our matrix
    ierr = crsMatrix.Multiply( false, updateVector.getEpetraVector(), 
                                residualVecPtr->getEpetraVector() );
    if( ierr )
    {
      cout << "ERROR: NOX::EpetraNew::BroydenOperator::runPostIterate(...) "
           << "- crsMatrix.Multiply() failed!!!" << endl;
      throw "NOX Error: Broyden Update Failed";
    }
    status = oldSolnGrp.applyJacobian(updateVector, *broydenVecPtr);
    if( status != NOX::Abstract::Group::Ok )
    {
      cout << "ERROR: NOX::EpetraNew::BroydenOperator::runPostIterate(...) "
           << "- applyJacobian() failed!!!" << endl;
      throw "NOX Error: Broyden Update Failed";
    }

    const Abstract::Group & solnGrp = solver.getSolutionGroup();
    // Form the difference needed for the outer product with the update vec
    broydenVecPtr->update(1.0, solnGrp.getF(), -1.0, *residualVecPtr, 1.0);

    double invUpdateNorm = 1.0 / updateVector.dot(updateVector);

    double * values = 0;
    int * indices = 0;
    int numEntries;
    for( int row = 0; row < crsMatrix.NumMyRows(); ++row)
    {
      ierr = crsMatrix.ExtractMyRowView(row, numEntries, values, indices);
      if( ierr )
      {
        cout << "ERROR (" << ierr << ") : "
             << "NOX::EpetraNew::BroydenOperator::runPostIterate(...) "
             << "- crsMatrix.ExtractGlobalRowView(...) failed for row --> "
             << row << endl;
        throw "NOX Error: Broyden Update Failed";
      }
      for( int col = 0; col < numEntries; ++col )
        // Could threshhold values here.
        (*values++) += broydenVecPtr->getEpetraVector()[row] * updateVector.getEpetraVector()[(*indices++)] * invUpdateNorm;

    }

    // Our %Broyden matrix has been updated and is now ready to use as a 
    // preconditioning matrix or as the Jacobian.
 
    if( verbose ) {
      NOX::Epetra::Vector * tempVecPtr2 = new NOX::Epetra::Vector(*tempVecPtr);
      ierr = crsMatrix.Multiply( false, updateVector.getEpetraVector(), 
                                 tempVecPtr2->getEpetraVector() );
      cout << "New Bs product vector : " << endl 
           << tempVecPtr2->getEpetraVector() << endl;
      tempVecPtr2->update(1.0, *tempVecPtr, -1.0);
      double maxDiff = tempVecPtr2->norm(NOX::Abstract::Vector::MaxNorm);
      cout << "Max difference applied to old update vector --> "
           << maxDiff << endl
           << "... and L-Inf norm of new residual -->"
           << solnGrp.getF().norm(NOX::Abstract::Vector::MaxNorm) << endl;
      delete tempVecPtr2; tempVecPtr2 = 0;
    }
  }

}

NOX::Parameter::Arbitrary * BroydenOperator::clone() const
{
  return new BroydenOperator(*this);
}

const string & BroydenOperator::getType() const
{
  return myType;
}

void BroydenOperator::replaceBroydenMatrixValues( const Epetra_CrsMatrix & mat)
{
  double * values = 0;
  int * indices = 0;
  int numEntries;
  int ierr;
  for( int row = 0; row < mat.NumMyRows(); ++row) {
    ierr = mat.ExtractMyRowView(row, numEntries, values, indices);
    ierr += crsMatrix.ReplaceGlobalValues(row, numEntries, values, indices);
    if( ierr )
    {
      cout << "ERROR (" << ierr << ") : "
           << "NOX::EpetraNew::BroydenOperator::replaceBroydenMatrixValues(...)"
           << " - Extract or Replace values error for row --> "
           << row << endl;
      throw "NOX Broyden Operator Error";
    }
  }
}
