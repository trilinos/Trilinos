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
                                 Epetra_CrsMatrix & mat) :
  updateVectorPtr( new NOX::Epetra::Vector(solnVec)),
  updateVector(*updateVectorPtr),
  broydenVecPtr(0),
  residualVecPtr(0),
  crsMatrix(mat),
  myType("Broyden Operator")
{
  // Set ourself as the Pre/Post Operator so we can update our data during
  // runPostIterate.  Note: we need to make provision for storing and calling 
  // any existing Pre/Post Operator.

  nlParams.sublist("Solver Options").setParameter(
                     "User Defined Pre/Post Operator", *this);
}

BroydenOperator::BroydenOperator(const BroydenOperator & bOp) :
  updateVectorPtr( new NOX::Epetra::Vector(bOp.updateVector) ),
  updateVector(*updateVectorPtr),
  broydenVecPtr(0),
  residualVecPtr(0),
  crsMatrix(bOp.crsMatrix),
  myType("Broyden Operator")
{ 
  if( bOp.broydenVecPtr )
    broydenVecPtr = new NOX::Epetra::Vector(*bOp.broydenVecPtr);
  if( bOp.residualVecPtr )
    residualVecPtr = new NOX::Epetra::Vector(*bOp.residualVecPtr);
}

//! Pure virtual destructor
BroydenOperator::~BroydenOperator()
{
  delete updateVectorPtr; updateVectorPtr = 0;
  delete broydenVecPtr  ; broydenVecPtr   = 0;
  delete residualVecPtr ; residualVecPtr  = 0;
}

bool BroydenOperator::computePreconditioner(const Epetra_Vector & x,
           NOX::Parameter::List * pList)
{
  return true;
}

void BroydenOperator::runPostIterate( const NOX::Solver::Generic & solver)
{ 
  // Get and update using the solver object.

  NOX::Abstract::Group::ReturnType status;
  int ierr;

  if( solver.getNumIterations() > 1 )
  {
    if( !broydenVecPtr )
      broydenVecPtr = new NOX::Epetra::Vector(updateVector);
    if( !residualVecPtr )
      residualVecPtr = new NOX::Epetra::Vector(updateVector);
      
    // Store previous Newton vector as the update, s
    const Abstract::Group & oldSolnGrp = solver.getPreviousSolutionGroup();
    if( !oldSolnGrp.isNewton() )
    {
      cout << "ERROR: NOX::EpetraNew::BroydenOperator::runPostIterate(...) "
           << "- getNewton() failed!!!" << endl;
      throw "NOX Error: Broyden Update Failed";
    }
    updateVector.update(1.0, oldSolnGrp.getNewton(), 0.0);

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

    // Form the difference needed for the outer product with the update vec
    broydenVecPtr->update(1.0, oldSolnGrp.getF(), -1.0, *residualVecPtr, 1.0);

    //cout << "\n\tBroyden vector :\n" << endl;
    //broydenVecPtr->print();
    
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

    cout << "Before updating Preconditioner matrix :" << endl;
    crsMatrix.Print(cout);

    // Our crsMatrix has been updated and is now ready to use as a 
    // preconditioner.
    cout << "After updating Preconditioner matrix :" << endl;
    crsMatrix.Print(cout);
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
