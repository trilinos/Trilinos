
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

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "NOX_Epetra_Interface.H"

#include "NOX_Epetra_MatrixFree.H"

using namespace NOX;
using namespace NOX::Epetra;

MatrixFree::MatrixFree(Interface& i, const Epetra_Vector& x) :
  interface(i),
  currentX(x),
  perturbX(x),
  fo(x),
  fp(x),
  importer(x.Map(),x.Map())
{
  // Zero out Vectors
  perturbX.PutScalar(0.0);
  fo.PutScalar(0.0);
  fp.PutScalar(0.0);
}

MatrixFree::~MatrixFree()
{
}

bool MatrixFree::Filled() const
{
  return true;
}

int MatrixFree::NumMyRowEntries(int MyRow, int & NumEntries) const
{
  cout << "throw NumMyRowEntires" << endl;
  throw;
}
  
int MatrixFree::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const
{
  cout << "throw ExtractMyRowCopy" << endl;
  throw;
}
  
int MatrixFree::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  cout << "throw ExtractDiagonalCopy" << endl;
  throw;
}

int MatrixFree::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // Use a directional derivative to compute y = Jx
  /*
   * eta = scalar perturbation
   * f = function evaluation (RHS)
   *
   *        f(x+eta*x) - f(x)
   * Jx =   -----------------
   *               eta
   */

  // Check if the transpose matrix is required.  Throw for now until we 
  // figure out how to handle this in a matrix free mode. 
  if (TransA == true) {
    cout << "ERROR: Matrix-Free can not compute a Transpose!" << endl;
    throw;  
  }

  // Compute eta
  double eta = 1.0e-4;

  // Compute the perturbed RHS
  perturbX = currentX;
  Y = X;
  Y.Scale(eta);
  perturbX.Update(1.0,Y,1.0);
  interface.computeRHS(perturbX,fp);
  
  // Compute the directional derivative
  Y.Update(1.0, fp, -1.0, fo, 0.0);
  Y.Scale(1.0/eta);

  return 0;
}

int MatrixFree::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X,  Epetra_MultiVector& Y) const
{
  // Check that flags we can't handle are not set.
  if (Upper == true) {
    cout << "ERROR: Matrix-Free can not compute an Upper Multiply!" << endl;
    throw;  
  }
  else if (UnitDiagonal == true) {
    cout << "ERROR: Matrix-Free can not compute a Diagonal Multiply!" << endl;
    throw;  
  }

  Multiply(Trans, X, Y);
  return 0;
}

int MatrixFree::InvRowSums(Epetra_Vector& x) const{}
  
int MatrixFree::LeftScale(const Epetra_Vector& x){}
  
int MatrixFree::InvColSums(Epetra_Vector& x) const{}
  
int MatrixFree::RightScale(const Epetra_Vector& x){}
  
double MatrixFree::NormInf() const{}

double MatrixFree::NormOne() const{}
  
int MatrixFree::NumGlobalNonzeros() const{}
  
int MatrixFree::NumGlobalRows() const{}
  
int MatrixFree::NumGlobalCols() const{}
  
int MatrixFree::NumGlobalDiagonals() const{}
  
int MatrixFree::NumMyNonzeros() const{}
  
int MatrixFree::NumMyRows() const{}
  
int MatrixFree::NumMyCols() const{}
  
int MatrixFree::NumMyDiagonals() const{}
  
bool MatrixFree::LowerTriangular() const{}

bool MatrixFree::UpperTriangular() const{}

const Epetra_Comm& MatrixFree::Comm() const
{
  return currentX.Comm();
}

const Epetra_BlockMap& MatrixFree::BlockRowMap() const
{
  return currentX.Map();
}

const Epetra_BlockMap& MatrixFree::BlockImportMap() const
{
  return currentX.Map();
}
  
const Epetra_Import* MatrixFree::Importer() const
{
  return &importer;
}

bool MatrixFree::computeJacobian(const Epetra_Vector& x, Epetra_RowMatrix& Jac)
{
  // Since we have no explicit Jacobian we set our currentX to the 
  // incoming value and evaluate the RHS.  When the Jacobian is used,
  // we compute the perturbed residuals and the directional 
  // derivative.
  currentX = x;

  return interface.computeRHS(x,fo);
}
