
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

#include "NOX_Epetra_Preconditioner_Jacobi.H" 

#include "Epetra_RowMatrix.h"
#include "Epetra_BlockMap.h"

namespace NOX {

namespace Epetra {

//! Constructor 
JacobiPreconditioner::JacobiPreconditioner(const Epetra_Vector& shape, 
					   const double value) :
  label("NOX::Epetra::JacobiPreconditioner"),
  diagonalVectorPtr(new Epetra_Vector(shape)),
  diagonalVector(*diagonalVectorPtr),
  minValue(value)
{}

//! Destructor
JacobiPreconditioner::~JacobiPreconditioner()
{
  delete diagonalVectorPtr;
}

int JacobiPreconditioner::SetUseTranspose(bool UseTranspose) 
{
  // Disable this option
  return false;
}

int JacobiPreconditioner::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
{
  // Not implemented: Throw an error!
  cout << "ERROR: NOX::Epetra::JacobiPreconditioner::Apply() - "
       << "method is NOT implemented!!  " << endl;
  throw "NOX Error";
  return false;
}

int JacobiPreconditioner::ApplyInverse(const Epetra_MultiVector& input, 
				       Epetra_MultiVector& result) const 
{

  // Calculate r = input ./ diagonalVector (./ is element-by-element divide)
  int retcode = result.ReciprocalMultiply(1.0, *diagonalVectorPtr, input, 0.0);

  // Check if this worked
  if (retcode != 0)
    return false;

  return true;
}
  
double JacobiPreconditioner::NormInf() const 
{
  // Not implemented: Throw an error!
  cout << "ERROR: NOX::Epetra::JacobiPreconditioner::NormInf() - "
       << "method is NOT implemented!!  " << endl;
  throw "NOX Error";
  return 0.0;
}

char* JacobiPreconditioner::Label () const 
{
  return const_cast<char*>(label.c_str());
}
  
bool JacobiPreconditioner::UseTranspose() const 
{
  // Not implemented: Throw an error!
  cout << "ERROR: NOX::Epetra::JacobiPreconditioner::UseTranspose() - "
       << "method is NOT implemented!!  " << endl;
  throw "NOX Error";
  return false;
}

bool JacobiPreconditioner::HasNormInf() const 
{
  // NormInf is not implemented
  return false;
}

const Epetra_Comm & JacobiPreconditioner::Comm() const 
{
  return diagonalVector.Map().Comm();
}
const Epetra_BlockMap& JacobiPreconditioner::DomainMap () const 
{
  return diagonalVector.Map();
}

const Epetra_BlockMap& JacobiPreconditioner::RangeMap () const 
{
  return diagonalVector.Map();
}
bool JacobiPreconditioner::computePreconditioner(const Epetra_Vector& x, 
					     const Epetra_Operator* jacobianPtr) 
{
  // We need a Jacobian pointer
  if (jacobianPtr == NULL) {
    cout << "ERROR: NOX::Epetra::JacobiPreconditioner::computePreconditioner() - "
	 << "the Jacobian point passed in is NULL.  This means the Jacobian is "
	 << "NOT valid with respect to the Group's curretn solution!" << endl;
    throw "NOX Error";
  }

  /* To extract the diagonal inverse, we must have a real matrix 
   * (we can NOT be using matrix-free operators).  Thus it must 
   * be an Epetra_RowMatrix.  Check for this and if not, throw an 
   * error.
   */
  // Try and cast it to an Epetra_RowMatrix
  const Epetra_RowMatrix* testRowMatrix = dynamic_cast<const Epetra_RowMatrix*>(jacobianPtr);
  if (testRowMatrix == NULL) {
    cout << "ERROR: NOX::Epetra::JacobiPreconditioner::computePreconditioner - "
	 << "the Jacobian operator is NOT a matrix!" << endl;
    throw "NOX Error";
  }

  // Convert the Epetra_RowMatrix to a reference
  const Epetra_RowMatrix& jacobian = *testRowMatrix;

  // Put a copy of the diagonal of the Jacobian into tmpVector
  int retcode = jacobian.ExtractDiagonalCopy(diagonalVector);
  
  // Check if ExtractDiagonalCopy is supported
  if (retcode != 0) {
    cout << "ERROR: NOX::Epetra::JacobiPreconditioner::computePreconditioner - "
	 << "ExtractDiagonalCopy() on Jacobian Epetra_Operator failed!" << endl;
    throw "NOX Error";
  }

  // Take element-wise absolute value of diagonal vector
  retcode = diagonalVector.Abs(diagonalVector);
  
  // Check minimum absolute value of diagonal vector
  double minAbsValue = 0.0;
  retcode = diagonalVector.MinValue(&minAbsValue);

  if(minAbsValue <= minValue) // This minimum threshold can be adjusted
  {
    cout << "Poor scaling on Jacobian diagonal (min abs value: " 
	 << minAbsValue << " ) --> No nonlinear preconditioning will be used!!" 
	 << endl;
    return false; 
  }
  
  return true;
}

} // namespace Epetra
} // namespace NOX
