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
       Teuchos::ParameterList & nlParams_                ,
       const Teuchos::RefCountPtr<NOX::Utils>& utils_    ,
       Epetra_Vector & solnVec                           ,
       const Teuchos::RefCountPtr<Epetra_CrsMatrix>& mat ,
       bool verbose_ ) :
  verbose(verbose_),
  crsMatrix(mat),
  nlParams(nlParams_),
  utils(utils_),
  prePostOperator( utils, nlParams.sublist("Solver Options") ),
  label("NOX::Epetra::BroydenOperator"),
  isValidStep(false),
  isValidYield(false),
  isValidBroyden(false)
{
  initialize( nlParams, solnVec );
}

//-----------------------------------------------------------------------------

BroydenOperator::BroydenOperator(const BroydenOperator & bOp) :
  verbose        ( bOp.verbose        ) ,
  stepVec        ( bOp.stepVec        ) ,
  yieldVec       ( bOp.yieldVec       ) ,
  workVec        ( bOp.workVec        ) ,
  oldX           ( bOp.oldX           ) ,
  oldF           ( bOp.oldF           ) ,
  crsMatrix      ( bOp.crsMatrix      ) ,
  jacIntPtr      ( bOp.jacIntPtr      ) ,
  jacMatrixPtr   ( bOp.jacMatrixPtr   ) ,
  precIntPtr     ( bOp.precIntPtr     ) ,
  precMatrixPtr  ( bOp.precMatrixPtr  ) ,
  nlParams       ( bOp.nlParams       ) ,
  utils          ( bOp.utils          ) ,
  prePostOperator( utils, nlParams.sublist("Solver Options") ),
  label          ( "NOX::Epetra::BroydenOperator" ),
  isValidStep    ( bOp.isValidStep    ) ,
  isValidYield   ( bOp.isValidYield   ) ,
  isValidBroyden ( bOp.isValidBroyden )
{ 
}

//-----------------------------------------------------------------------------

bool
BroydenOperator::initialize( Teuchos::ParameterList & nlParams, const Epetra_Vector & vec )
{
  stepVec  = Teuchos::rcp( new NOX::Epetra::Vector(vec) );
  yieldVec = Teuchos::rcp( new NOX::Epetra::Vector(vec) );
  workVec  = Teuchos::rcp( new NOX::Epetra::Vector(vec) );
  oldX     = Teuchos::rcp( new NOX::Epetra::Vector(vec) );
  oldF     = Teuchos::rcp( new NOX::Epetra::Vector(vec) );

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

const char* BroydenOperator::Label () const
{
  return label.c_str();
}

int BroydenOperator::SetUseTranspose(bool UseTranspose)
{
  return crsMatrix->SetUseTranspose(UseTranspose);
}

int BroydenOperator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return crsMatrix->Apply(X, Y);
}

int BroydenOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return crsMatrix->ApplyInverse(X, Y);
}

bool BroydenOperator::UseTranspose() const
{
  return crsMatrix->UseTranspose();
}

bool BroydenOperator::HasNormInf() const
{
  return crsMatrix->HasNormInf();
}

const Epetra_Map& BroydenOperator::OperatorDomainMap() const
{
  return crsMatrix->OperatorDomainMap();
}

const Epetra_Map& BroydenOperator::OperatorRangeMap() const
{
  return crsMatrix->OperatorRangeMap();
}

bool BroydenOperator::Filled() const
{
  return crsMatrix->Filled();
}

int BroydenOperator::NumMyRowEntries(int MyRow, int & NumEntries) const
{
  return crsMatrix->NumMyRowEntries(MyRow, NumEntries);
}

int BroydenOperator::MaxNumEntries() const
{
  return crsMatrix->MaxNumEntries();
}

inline int BroydenOperator::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const
{
  return crsMatrix->ExtractMyRowCopy(MyRow, Length, NumEntries, Values, Indices);
}

int BroydenOperator::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  return crsMatrix->ExtractDiagonalCopy(Diagonal);
}

int BroydenOperator::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return crsMatrix->Multiply(TransA, X, Y);
}

int BroydenOperator::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X,  Epetra_MultiVector& Y) const
{
  return crsMatrix->Solve(Upper, Trans, UnitDiagonal, X, Y);
}

int BroydenOperator::InvRowSums(Epetra_Vector& x) const
{
  return crsMatrix->InvRowSums(x);
}

int BroydenOperator::LeftScale(const Epetra_Vector& x)
{
  return crsMatrix->LeftScale(x);
}

int BroydenOperator::InvColSums(Epetra_Vector& x) const
{
  return crsMatrix->InvColSums(x);
}

int BroydenOperator::RightScale(const Epetra_Vector& x)
{
  return crsMatrix->RightScale(x);
}

double BroydenOperator::NormInf() const
{
  return crsMatrix->NormInf();
}

double BroydenOperator::NormOne() const
{
  return crsMatrix->NormOne();
}

int BroydenOperator::NumGlobalNonzeros() const
{
  return crsMatrix->NumGlobalNonzeros();
}

int BroydenOperator::NumGlobalRows() const
{
  return crsMatrix->NumGlobalRows();
}

int BroydenOperator::NumGlobalCols() const
{
  return crsMatrix->NumGlobalCols();
}

int BroydenOperator::NumGlobalDiagonals() const
{
  return crsMatrix->NumGlobalDiagonals();
}

int BroydenOperator::NumMyNonzeros() const
{
  return crsMatrix->NumMyNonzeros();
}

int BroydenOperator::NumMyRows() const
{
  return crsMatrix->NumMyRows();
}

int BroydenOperator::NumMyCols() const
{
  return crsMatrix->NumMyCols();
}

int BroydenOperator::NumMyDiagonals() const
{
  return crsMatrix->NumMyDiagonals();
}

bool BroydenOperator::LowerTriangular() const
{
  return crsMatrix->LowerTriangular();
}

bool BroydenOperator::UpperTriangular() const
{
  return crsMatrix->UpperTriangular();
}

const Epetra_Comm& BroydenOperator::Comm() const
{
  return crsMatrix->Comm();
}

const Epetra_Map& BroydenOperator::RowMatrixRowMap() const
{
  return crsMatrix->RowMatrixRowMap();
}

const Epetra_Map& BroydenOperator::RowMatrixColMap() const
{
  return crsMatrix->RowMatrixColMap();
}

const Epetra_Import* BroydenOperator::RowMatrixImporter() const
{
  return crsMatrix->RowMatrixImporter();
}

const Epetra_BlockMap& BroydenOperator::Map() const
{
  return crsMatrix->Map();
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::computeJacobian( const Epetra_Vector & x, Epetra_Operator& Jac )
{
  return computeSparseBroydenUpdate();
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::computePreconditioner( const Epetra_Vector & x,
					Epetra_Operator& Prec,
                                        Teuchos::ParameterList * pList )
{
  return computeJacobian( x, Prec );
}

//-----------------------------------------------------------------------------

void 
BroydenOperator::runPreSolve( const NOX::Solver::Generic & solver)
{
  // Reset valid step and yield flags. NOTE: we leave the Broyden valid flag
  // alone as it may be valid from a previous solve
  isValidStep  = false;
  isValidYield = false;

  prePostOperator.runPreSolve( solver );

  return;
}

//-----------------------------------------------------------------------------

void 
BroydenOperator::runPreIterate( const NOX::Solver::Generic & solver)
{
  // We need to store the "old" solution and corresponding residual vectors here
  // in a manner consistent with infrequent updates of the preconditioner
  const NOX::Epetra::Group & grp = 
    dynamic_cast<const NOX::Epetra::Group &>(solver.getSolutionGroup());

  NOX::Epetra::LinearSystem::PreconditionerReusePolicyType reuse = 
    Teuchos::rcp_const_cast<NOX::Epetra::LinearSystem>(grp.getLinearSystem())->getPreconditionerPolicy(false);

  if( NOX::Epetra::LinearSystem::PRPT_REBUILD   == reuse ||
      NOX::Epetra::LinearSystem::PRPT_RECOMPUTE == reuse    )
  {
    *oldX = solver.getSolutionGroup().getX();
    *oldF = solver.getSolutionGroup().getF();
    if( verbose )
      cout << "\t...Updated oldX and oldF..." << endl;
  }

  prePostOperator.runPreIterate( solver );

  return;
}

//-----------------------------------------------------------------------------

void 
BroydenOperator::runPostIterate( const NOX::Solver::Generic & solver)
{
  // Get and update using the solver object.

  if( solver.getNumIterations() > 0 )
  {
    // To consistently compute changes to the step and the yield, eg effects of
    // linesearch or other globalizations, we need to get and subtract the old and
    // new solution states and corresponding residuals 

    const Abstract::Group & oldSolnGrp = solver.getPreviousSolutionGroup();
    const Abstract::Group & solnGrp    = solver.getSolutionGroup();

    // Set the Step vector
#ifdef HAVE_NOX_DEBUG
    if( verbose )
    {
      *workVec = *oldX;
      workVec->update(-1.0, oldSolnGrp.getX(), 1.0);
      double norm = workVec->norm();
      cout << "BroydenOperator::runPostIterate(...) : oldX diff norm = " << norm << endl;
    }
#endif
    *workVec = solnGrp.getX();
    workVec->update(-1.0, oldSolnGrp.getX(), 1.0);
    //workVec->update(-1.0, *oldX, 1.0);

    setStepVector( *workVec );

    // Set the Yield vector
#ifdef HAVE_NOX_DEBUG
    if( verbose )
    {
      *workVec = *oldF;
      workVec->update(-1.0, oldSolnGrp.getF(), 1.0);
      double norm = workVec->norm();
      cout << "BroydenOperator::runPostIterate(...) : oldF diff norm = " << norm << endl;
    }
#endif
    *workVec = solnGrp.getF();
    workVec->update(-1.0, oldSolnGrp.getF(), 1.0);
    //workVec->update(-1.0, *oldF, 1.0);

    setYieldVector( *workVec );

    // Use of these updated vectors occurs when computeJacobian or computePreconditioner 
    // gets called
  }

  prePostOperator.runPostIterate( solver );

  return;
}

//-----------------------------------------------------------------------------

void 
BroydenOperator::runPostSolve( const NOX::Solver::Generic & solver)
{
  prePostOperator.runPostSolve( solver );

  return;
}

//-----------------------------------------------------------------------------

void
BroydenOperator::setStepVector( Epetra_Vector & vec )
{
  stepVec->getEpetraVector() = vec;
  isValidStep    = true  ;
  isValidBroyden = false ;

  return;
}

//-----------------------------------------------------------------------------

void
BroydenOperator::setStepVector( NOX::Epetra::Vector & vec )
{
  *stepVec = vec;
  isValidStep    = true  ;
  isValidBroyden = false ;

  return;
}

//-----------------------------------------------------------------------------

void
BroydenOperator::setYieldVector( NOX::Epetra::Vector & vec )
{
  *yieldVec = vec;
  isValidYield   = true  ;
  isValidBroyden = false ;

  return;
}

//-----------------------------------------------------------------------------

void
BroydenOperator::setYieldVector( Epetra_Vector & vec )
{
  yieldVec->getEpetraVector() = vec;
  isValidYield   = true  ;
  isValidBroyden = false ;

  return;
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::computeSparseBroydenUpdate()
{

  if( isValidBroyden )
    return true;
  else if( isValidStep && isValidYield )
  {
    // Do the Broyden update to our matrix
    int ierr = crsMatrix->Multiply( false, stepVec->getEpetraVector(), workVec->getEpetraVector() );
    if( ierr )
    {
      cout << "ERROR: NOX::Epetra::BroydenOperator::computeSparseBroydenUpdate(...) "
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
        cout << "ERROR (" << ierr << ") : NOX::Epetra::BroydenOperator::computeSparseBroydenUpdate() "
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

    if( verbose )
      cout << "\t... BroydenOperator::computeSparseBroydenUpdate()..." << endl;

  }
  else
  {
    cout << "Warning: NOX::Epetra::BroydenOperator::computeSparseBroydenUpdate(...) "
         << "- either the step vector or the yield vector or both is not valid." << endl;
    cout <<  "Leaving existing matrix unchanged." << endl;
  }

  return true;
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::isStep() const
{
  return isValidStep;
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::isYield() const
{
  return isValidYield;
}

//-----------------------------------------------------------------------------

bool 
BroydenOperator::isBroyden() const
{
  return isValidBroyden;
}

//-----------------------------------------------------------------------------

void 
BroydenOperator::replaceBroydenMatrixValues( const Epetra_CrsMatrix & mat)
{
  double * values = 0;
  int * indices = 0;
  int numEntries;
  int ierr;
  for( int row = 0; row < mat.NumMyRows(); ++row) 
  {
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
