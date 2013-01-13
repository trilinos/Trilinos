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
    throw "NOX Error";
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
  throw "NOX Error";

  return (-1);
}

double 
SchurCoupler::NormInf() const
{
  std::cout << "ERROR: NOX::Epetra::SchurCoupler::NormInf() - Not Available for "
       << "this operator!" << std::endl;
  throw "NOX Error";

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
    throw "NOX Error";
  }

  return schurOperators[rowBlock];
}

