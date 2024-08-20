// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
                                                                               #include "Epetra_Map.h"
#include "NOX_Epetra_SchurCoupler.H"

using namespace NOX;
using namespace NOX::Epetra;

SchurCoupler::SchurCoupler( Problem_Manager & probManager ) :
  problemManager(probManager),
  label("Schur_Coupler")
{
  int numProblems = probManager.getProblemCount();

  // For now, just create diagonal Schur operators
  for( int i = 0; i < numProblems - 1; ++i )
    schurOperators[i] = Teuchos::rcp( new SchurOp(i, i+1, *this ) );
}

SchurCoupler::~SchurCoupler()
{
}

int
SchurCoupler::SetUseTranspose( bool UseTranspose )
{
  if (UseTranspose == true)
  {
    std::cout << "ERROR: NOX::Epetra::SchurCoupler::SetUseTranspose() - Transpose is "
     << "unavailable for this operator!" << std::endl;
    throw std::runtime_error("NOX Error");
  }
  return (-1);
}

int
SchurCoupler::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
//  // Convert X and Y from an Epetra_MultiVectors to Epetra_Vectors
//  // and NOX::Epetra::Vectors.  This is done so we use a consistent
//  // vector space for norms and inner products.
//  Teuchos::RCP<Epetra_Vector> wrappedX = Teuchos::rcp(new Epetra_Vector(View, X, 0));
//  Teuchos::RCP<Epetra_Vector> wrappedY = Teuchos::rcp(new Epetra_Vector(View, Y, 0));
//  Teuchos::RCP<Epetra_Vector> tempX    = Teuchos::rcp(new Epetra_Vector(*wrappedX) );
//  Teuchos::RCP<Epetra_Vector> tempY    = Teuchos::rcp(new Epetra_Vector(*wrappedY) );
//
//  // The substance of this operator ----- to do RWH 10/24/2006
//
//  //cout << "Incoming X :\n" << *wrappedX << std::endl;
//  problemManager.applyBlockAction( depId, probId, *wrappedX, *tempY);
//  //cout << "After applyBlockAction(" << depId << ", " << probId << ") :\n" << *tempY << std::endl;
//  problemManager.getBlockInverseOperator(depId)->ApplyInverse(*tempY, *tempY);
//  //cout << "After ApplyInverse :\n" << *tempY << std::endl;
//  problemManager.applyBlockAction( probId, depId, *tempY, *tempX);
//  //cout << "After applyBlockAction(" << probId << ", " << depId << ") :\n" << *tempX << std::endl;
//
//  problemManager.getBlockJacobianMatrix(probId)->Apply(*wrappedX, *wrappedY);
//  //cout << "After Apply of diagonal block :\n" << *wrappedY << std::endl;
//  wrappedY->Update(-1.0, *tempX, 1.0);
//  //cout << "After combingin; final result :\n" << *wrappedY << std::endl;
//
  return (0);
}

int
SchurCoupler::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  std::cout << "ERROR: NOX::Epetra::SchurCoupler::ApplyInverse() - Not valid "
       << "for this operator!" << std::endl;
  throw std::runtime_error("NOX Error");

  return (-1);
}

double
SchurCoupler::NormInf() const
{
  std::cout << "ERROR: NOX::Epetra::SchurCoupler::NormInf() - Not Available for "
       << "this operator!" << std::endl;
  throw std::runtime_error("NOX Error");

  return 1.0;
}

const char*
SchurCoupler::Label () const
{
  return label.c_str();
}

bool
SchurCoupler::UseTranspose() const
{
  return false;
}

bool
SchurCoupler::HasNormInf() const
{
  return false;
}

const Epetra_Comm &
SchurCoupler::Comm() const
{
  return problemManager.getJacobianOperator()->Comm();
}

const Epetra_Map &
SchurCoupler::OperatorDomainMap() const
{
  return *(problemManager.getCompositeMap());
}

const Epetra_Map &
SchurCoupler::OperatorRangeMap() const
{
  return *(problemManager.getCompositeMap());
}


bool
SchurCoupler::applyBlockAction( int rowBlock, int colBlock, const Epetra_Vector & x, Epetra_Vector & result )
{
  // Note that we need to convert to 1-based when interacting with Problem_Manager
  return problemManager.applyBlockAction( rowBlock+1, colBlock+1, x, result );
}

bool
SchurCoupler::applyBlockInverseAction( int rowBlock, int colBlock, const Epetra_Vector & x, Epetra_Vector & result )
{
  return problemManager.applyBlockInverseAction( rowBlock+1, colBlock+1, x, result );
}

bool
SchurCoupler::hasExplicitOperator( int rowBlock, int colBlock )
{
  return ( rowBlock == colBlock );
}

Teuchos::RCP<Epetra_Operator>
SchurCoupler::getExplicitOperator( int rowBlock, int colBlock )
{
  return problemManager.getBlockJacobianMatrix(rowBlock+1);
}

Teuchos::RCP<SchurOp>
SchurCoupler::getSchurOperator( int rowBlock, int colBlock )
{

  if( rowBlock != colBlock )
  {
    std::cout << "ERROR: NOX::Epetra::SchurCoupler::getSchurOperator() - Off operators not supported."
         << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return schurOperators[rowBlock];
}

