
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
#include "Epetra_RowMatrix.h"
#include "NOX_Epetra_Interface.H"
#include "NOX_Utils.H"

#include "NOX_Epetra_MatrixFree.H"

using namespace NOX;
using namespace NOX::Epetra;

MatrixFree::MatrixFree(Interface& i, const Epetra_Vector& x) :
  label("NOX::Matrix-Free"),
  interface(i),
  currentX(x),
  perturbX(x),
  fo(x),
  fp(x)
{
  // Zero out Vectors
  perturbX.PutScalar(0.0);
  fo.PutScalar(0.0);
  fp.PutScalar(0.0);
}

MatrixFree::~MatrixFree()
{
}

int MatrixFree::SetUseTranspose(bool UseTranspose) 
{
  if (UseTranspose == true) {
    cout << "ERROR: NOX::Epetra::MatrixFree::SetUseTranspose() - Transpose is "
	 << "unavailable in Matrix-Free mode!" << endl;
    throw "NOX Error";
  }
  return (-1);
}

int MatrixFree::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  // Use a directional derivative to compute y = Jx
  /*
   * eta = scalar perturbation
   * u = solution vector used to evaluate f
   * f = function evaluation (RHS)
   * x = vector that J is applied to
   *
   *        f(u+eta*x) - f(u)
   * Jx =   -----------------
   *               eta
   */

  // Compute perturbation constant, eta
  // Taken from LOCA v1.0 manual SAND2002-0396 p. 28 eqn. 2.43
  // eta = lambda*(lambda + 2norm(u)/2norm(x))
  double lambda = 1.0e-6;
  double solutionNorm = 1.0;
  double vectorNorm = 1.0;

  int test = currentX.Norm2(&solutionNorm);

  // Make sure the norm computed correctly
  if (test != 0) {
    if (NOX::Utils::doPrint(Utils::Warning)) 
      cout << "Warning: NOX::Epetra::MatrixFree::Apply() - solutionNorm "
	   << "failed!" << endl;
    solutionNorm = 1.0;
  }

  test = X.Norm2(&vectorNorm);

  // Make sure the norm computed correctly
  if (test != 0) {
    if (NOX::Utils::doPrint(Utils::Warning)) 
      cout << "Warning: NOX::Epetra::MatrixFree::Apply() - vectorNorm failed!" 
	   << endl;
    vectorNorm = 1.0;
  }

  // Make sure the norm is not zero, otherwise we can get an inf perturbation
  if (vectorNorm == 0.0) {
    if (NOX::Utils::doPrint(Utils::Warning)) 
      cout << "Warning: NOX::Epetra::MatrixFree::Apply() - vectorNorm is zero" 
	   << endl;
    vectorNorm = 1.0;
  }

  double eta = lambda*(lambda + solutionNorm/vectorNorm);

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

int MatrixFree::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  cout << "ERROR: NOX::MatrixFree::ApplyInverse - Not available for Matrix Free!"
       << endl;
  throw "NOX Error";
  return (-1);
}

double MatrixFree::NormInf() const
{
  cout << "ERROR: NOX::Epetra::MatrixFree::NormInf() - Not Available for "
       << "Matrix-Free mode!" << endl;
  throw "NOX Error";
  return 1.0;
}


char* MatrixFree::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool MatrixFree::UseTranspose() const
{
  return false;
}

bool MatrixFree::HasNormInf() const
{
  return false;
}

const Epetra_Comm & MatrixFree::Comm() const
{
  return currentX.Map().Comm();
}
const Epetra_BlockMap& MatrixFree::DomainMap() const
{
  return currentX.Map();
}

const Epetra_BlockMap& MatrixFree::RangeMap() const
{
  return currentX.Map();
}

bool MatrixFree::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  // Since we have no explicit Jacobian we set our currentX to the 
  // incoming value and evaluate the RHS.  When the Jacobian is applied,
  // we compute the perturbed residuals and the directional 
  // derivative.
  currentX = x;

  return interface.computeRHS(x,fo);
}
