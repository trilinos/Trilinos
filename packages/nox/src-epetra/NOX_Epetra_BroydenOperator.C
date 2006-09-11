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

#include "NOX_Common.H"
#include "NOX_Epetra_BroydenOperator.H"

#include "Epetra_Map.h"

using namespace NOX;
using namespace NOX::Epetra;

BroydenOperator::BroydenOperator(
       Teuchos::ParameterList & nlParams, 
       Epetra_Vector & solnVec,
       const Teuchos::RefCountPtr<Epetra_CrsMatrix>& mat,
       bool verbose_ ) :
  verbose(verbose_),
  updateVectorPtr( Teuchos::rcp(new NOX::Epetra::Vector(solnVec)) ),
  updateVector(*updateVectorPtr),
  crsMatrix(mat),
  myType("Broyden Operator")
{
  initialize( nlParams, solnVec );
}

//-----------------------------------------------------------------------------

BroydenOperator::BroydenOperator(
      Teuchos::ParameterList & nlParams, 
      Epetra_Vector & solnVec,
      const Teuchos::RefCountPtr<Epetra_CrsMatrix> & mat,
      const Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian>& jacInt,
      const Teuchos::RefCountPtr<Epetra_CrsMatrix>& jacMatrix, 
      bool verbose_ ) :
  verbose(verbose_),
  updateVectorPtr( Teuchos::rcp(new NOX::Epetra::Vector(solnVec)) ),
  updateVector(*updateVectorPtr),
  crsMatrix(mat),
  jacIntPtr(jacInt),
  jacMatrixPtr(jacMatrix),
  myType("Broyden Operator")
{
  // Set ourself as the Pre/Post Operator so we can update our data during
  // runPostIterate.  Note: we need to make provision for storing and calling 
  // any existing Pre/Post Operator.

  // RPP 9/20/2005: This is a very bad idea!  It breaks the rcp and
  // user's expectations.  For now we will have to create rcp without
  // ownership.  What happens if a user write their own PPO?
  Teuchos::RefCountPtr<NOX::Abstract::PrePostOperator> me = Teuchos::rcp(this, false);
  nlParams.sublist("Solver Options").set("User Defined Pre/Post Operator", me);

  // Create out working vectors
}

//-----------------------------------------------------------------------------

BroydenOperator::BroydenOperator(
 Teuchos::ParameterList & nlParams, 
 Epetra_Vector & solnVec,
 const Teuchos::RefCountPtr<Epetra_CrsMatrix>& mat,
 const Teuchos::RefCountPtr<NOX::Epetra::Interface::Preconditioner>& precInt,
 const Teuchos::RefCountPtr<Epetra_CrsMatrix>& precMatrix, 
 bool verbose_ ) :
  verbose(verbose_),
  updateVectorPtr( Teuchos::rcp(new NOX::Epetra::Vector(solnVec)) ),
  updateVector(*updateVectorPtr),
  crsMatrix(mat),
  precIntPtr(precInt),
  precMatrixPtr(precMatrix),
  myType("Broyden Operator")
{
  // Set ourself as the Pre/Post Operator so we can update our data during
  // runPostIterate.  Note: we need to make provision for storing and calling 
  // any existing Pre/Post Operator.

  // RPP 9/20/2005: This is a very bad idea!  It breaks the rcp and
  // user's expectations.  For now we will have to create rcp without
  // ownership.  What happens if a user write their own PPO?
  Teuchos::RefCountPtr<NOX::Abstract::PrePostOperator> me = Teuchos::rcp(this, false);
  nlParams.sublist("Solver Options").set("User Defined Pre/Post Operator", me);
}

//-----------------------------------------------------------------------------

BroydenOperator::BroydenOperator(const BroydenOperator & bOp) :
  verbose(bOp.verbose),
  updateVectorPtr(Teuchos::rcp(new NOX::Epetra::Vector(*bOp.updateVectorPtr))),
  updateVector(*updateVectorPtr),
  crsMatrix(bOp.crsMatrix),
  jacIntPtr(bOp.jacIntPtr),
  jacMatrixPtr(bOp.jacMatrixPtr),
  precIntPtr(bOp.precIntPtr),
  precMatrixPtr(bOp.precMatrixPtr),
  myType("Broyden Operator")
{ 
  if( !Teuchos::is_null(bOp.broydenVecPtr) )
    broydenVecPtr = Teuchos::rcp(new NOX::Epetra::Vector(*bOp.broydenVecPtr));
  if( !Teuchos::is_null(bOp.residualVecPtr) )
    residualVecPtr = Teuchos::rcp(new NOX::Epetra::Vector(*bOp.residualVecPtr));
  if( !Teuchos::is_null(bOp.tempVecPtr) )
    tempVecPtr = Teuchos::rcp(new NOX::Epetra::Vector(*bOp.tempVecPtr));
}

//-----------------------------------------------------------------------------

bool
BroydenOperator::initialize( Teuchos::ParameterList & nlParams, const Epetra_Vector & vec )
{
  stepVec  = Teuchos::rcp( new NOX::Epetra::Vector(vec) );
  yieldVec = Teuchos::rcp( new NOX::Epetra::Vector(vec) );
  workVec  = Teuchos::rcp( new NOX::Epetra::Vector(vec) );

  // Set ourself as the Pre/Post Operator so we can update our data during
  // runPostIterate.  Note: we need to make provision for storing and calling 
  // any existing Pre/Post Operator.

  // RPP 9/20/2005: This is a very bad idea!  It breaks the rcp and
  // user's expectations.  For now we will have to create rcp without
  // ownership.  What happens if a user write their own PPO?
  Teuchos::RefCountPtr<NOX::Abstract::PrePostOperator> me = Teuchos::rcp(this, false);
  nlParams.sublist("Solver Options").set("User Defined Pre/Post Operator", me);

  return true;
}

//-----------------------------------------------------------------------------

BroydenOperator::~BroydenOperator()
{
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::computeJacobian( const Epetra_Vector & x, Epetra_Operator& Jac )
{
  if( !Teuchos::is_null(jacIntPtr) ) 
  {
    jacIntPtr->computeJacobian(x, Jac);
    replaceBroydenMatrixValues(*jacMatrixPtr);
  }

  return true;
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::computePreconditioner( const Epetra_Vector & x,
					     Epetra_Operator& Prec,
                                             Teuchos::ParameterList * pList )
{
  if( !Teuchos::is_null(precIntPtr) ) 
  {
    precIntPtr->computePreconditioner(x, Prec, pList);
    replaceBroydenMatrixValues(*precMatrixPtr);
  }

  return true;
}

//-----------------------------------------------------------------------------

void 
BroydenOperator::runPostIterate( const NOX::Solver::Generic & solver)
{ 
  // Get and update using the solver object.

  NOX::Abstract::Group::ReturnType status;
  int ierr;

  if( solver.getNumIterations() > 0 )
  {
    if( Teuchos::is_null(broydenVecPtr) )
      broydenVecPtr = Teuchos::rcp(new NOX::Epetra::Vector(*updateVectorPtr));
    if( Teuchos::is_null(residualVecPtr) )
      residualVecPtr = Teuchos::rcp(new NOX::Epetra::Vector(*updateVectorPtr));
    if( verbose )
      if( Teuchos::is_null(tempVecPtr) )
        tempVecPtr = Teuchos::rcp(new NOX::Epetra::Vector(*updateVectorPtr));
      
    // Store previous Newton vector as the update, s
    const Abstract::Group & oldSolnGrp = solver.getPreviousSolutionGroup();
    if( !oldSolnGrp.isNewton() )
    {
      cout << "ERROR: NOX::Epetra::BroydenOperator::runPostIterate(...) "
           << "- getNewton() failed!!!" << endl;
      throw "NOX Error: Broyden Update Failed";
    }
    updateVector.update(1.0, oldSolnGrp.getNewton(), 0.0);
    if( verbose ) 
    {
      oldSolnGrp.applyJacobian(updateVector, *tempVecPtr);
       cout << "Js vector : " << endl << tempVecPtr->getEpetraVector()
            << "\nOld residual vector : " << endl 
            << dynamic_cast<const NOX::Epetra::Vector&>(oldSolnGrp.getF()).getEpetraVector() << endl;
    }
      

    // Do the Broyden update to our matrix
    ierr = crsMatrix->Multiply( false, updateVector.getEpetraVector(), 
                                residualVecPtr->getEpetraVector() );
    if( ierr )
    {
      cout << "ERROR: NOX::Epetra::BroydenOperator::runPostIterate(...) "
           << "- crsMatrix->Multiply() failed!!!" << endl;
      throw "NOX Error: Broyden Update Failed";
    }

    status = oldSolnGrp.applyJacobian(updateVector, *broydenVecPtr);
    if( status != NOX::Abstract::Group::Ok )
    {
      cout << "ERROR: NOX::Epetra::BroydenOperator::runPostIterate(...) "
           << "- applyJacobian() failed!!!" << endl;
      throw "NOX Error: Broyden Update Failed";
    }

    const Abstract::Group & solnGrp = solver.getSolutionGroup();
    // Form the difference needed for the outer product with the update vec
    broydenVecPtr->update(1.0, solnGrp.getF(), -1.0, *residualVecPtr, 1.0);

    double invUpdateNorm = 1.0 / updateVector.innerProduct(updateVector);

    double * values = 0;
    int * indices = 0;
    int numEntries;
    for( int row = 0; row < crsMatrix->NumMyRows(); ++row)
    {
      ierr = crsMatrix->ExtractMyRowView(row, numEntries, values, indices);
      if( ierr )
      {
        cout << "ERROR (" << ierr << ") : "
             << "NOX::Epetra::BroydenOperator::runPostIterate(...) "
             << "- crsMatrix->ExtractGlobalRowView(...) failed for row --> "
             << row << endl;
        throw "NOX Error: Broyden Update Failed";
      }
      for( int col = 0; col < numEntries; ++col )
        // Could threshhold values here.
        (*values++) += broydenVecPtr->getEpetraVector()[row] * updateVector.getEpetraVector()[(*indices++)] * invUpdateNorm;

    }

    // Our %Broyden matrix has been updated and is now ready to use as a 
    // preconditioning matrix or as the Jacobian.
    if( verbose ) 
    {
      NOX::Epetra::Vector * tempVecPtr2 = new NOX::Epetra::Vector(*tempVecPtr);
      ierr = crsMatrix->Multiply( false, updateVector.getEpetraVector(), 
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

//-----------------------------------------------------------------------------

const string & 
BroydenOperator::getType() const
{
  return myType;
}

//-----------------------------------------------------------------------------

void
BroydenOperator::setStepVector( Epetra_Vector & vec )
{
  stepVec->getEpetraVector() = vec;
}

//-----------------------------------------------------------------------------

void
BroydenOperator::setStepVector( NOX::Epetra::Vector & vec )
{
  *stepVec = vec;
}

//-----------------------------------------------------------------------------

void
BroydenOperator::setYieldVector( NOX::Epetra::Vector & vec )
{
  *yieldVec = vec;
}

//-----------------------------------------------------------------------------

void
BroydenOperator::setYieldVector( Epetra_Vector & vec )
{
  yieldVec->getEpetraVector() = vec;
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::computeSparseBroydenUpdate()
{
  //NOX::Abstract::Group::ReturnType status;
  //int ierr;

  //if( Teuchos::is_null(broydenVecPtr) )
  //  broydenVecPtr = Teuchos::rcp(new NOX::Epetra::Vector(*updateVectorPtr));
  //if( Teuchos::is_null(residualVecPtr) )
  //  residualVecPtr = Teuchos::rcp(new NOX::Epetra::Vector(*updateVectorPtr));
  //if( verbose )
  //  if( Teuchos::is_null(tempVecPtr) )
  //    tempVecPtr = Teuchos::rcp(new NOX::Epetra::Vector(*updateVectorPtr));
    
  // Do the Broyden update to our matrix
  int ierr = crsMatrix->Multiply( false, stepVec->getEpetraVector(), workVec->getEpetraVector() );
  if( ierr )
  {
    cout << "ERROR: NOX::Epetra::BroydenOperator::runPostIterate(...) "
         << "- crsMatrix->Multiply() failed!!!" << endl;
    throw "NOX Error: Broyden Update Failed";
  }

  double * values  = 0 ;
  int    * indices = 0 ;
  int      numEntries  ;

  for( int row = 0; row < crsMatrix->NumMyRows(); ++row )
  {
    ierr = crsMatrix->ExtractMyRowView( row, numEntries, values, indices );
    if( ierr )
    {
      cout << "ERROR (" << ierr << ") : NOX::Epetra::BroydenOperator::computeSparseBroydenUdate() "
           << "- crsMatrix->ExtractGlobalRowView(...) failed for row --> " << row << endl;
      throw "NOX Error: Broyden Update Failed";
    }

    double diffVal = yieldVec->getEpetraVector()[row] - workVec->getEpetraVector()[row];

    double rowStepInnerProduct = 0.0;

    for( int col = 0; col < numEntries; ++col )
      rowStepInnerProduct += stepVec->getEpetraVector()[indices[col]] * stepVec->getEpetraVector()[indices[col]];

    for( int col = 0; col < numEntries; ++col )
      (*values++) += diffVal * stepVec->getEpetraVector()[(*indices++)] / rowStepInnerProduct;
  }

  // Our %Broyden matrix has been updated and is now ready to use as a 
  // preconditioning matrix or as the Jacobian.

  return true;
}

//-----------------------------------------------------------------------------

void 
BroydenOperator::replaceBroydenMatrixValues( const Epetra_CrsMatrix & mat)
{
  double * values = 0;
  int * indices = 0;
  int numEntries;
  int ierr;
  for( int row = 0; row < mat.NumMyRows(); ++row) {
    ierr = mat.ExtractMyRowView(row, numEntries, values, indices);
    ierr += crsMatrix->ReplaceGlobalValues(row, numEntries, values, indices);
    if( ierr )
    {
      cout << "ERROR (" << ierr << ") : "
           << "NOX::Epetra::BroydenOperator::replaceBroydenMatrixValues(...)"
           << " - Extract or Replace values error for row --> "
           << row << endl;
      throw "NOX Broyden Operator Error";
    }
  }
}

//-----------------------------------------------------------------------------
