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

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

#include "NOX_Epetra_Interface.H"
#include "NOX_Utils.H"

#include "NOX_Epetra_FiniteDifference.H"

using namespace NOX;
using namespace NOX::Epetra;

FiniteDifference::FiniteDifference(Interface& i, const Epetra_Vector& x, double beta_, double alpha_) :
  map(x.Map()),
  graph(0),
  jacobian(0),
  interface(i),
  x_perturb(x),
  fo(x),
  fp(x),
  fmPtr(0),
  Jc(x),
  alpha(alpha_),
  beta(beta_),
  betaVector(0),
  betaType(Scalar),
  diffType(Forward),
  label("NOX::FiniteDifference Jacobian")
{
  // Create the finite difference Jacobian matrix
  jacobian = createGraphAndJacobian(i, x);
}

FiniteDifference::FiniteDifference(Interface& i, const Epetra_Vector& x, const Epetra_Vector& beta_, double alpha_) :
  map(x.Map()),
  graph(0),
  jacobian(0),
  interface(i),
  x_perturb(x),
  fo(x),
  fp(x),
  fmPtr(0),
  Jc(x),
  alpha(alpha_),
  beta(0),
  betaVector(&beta_),
  betaType(Vector),
  diffType(Forward),
  label("NOX::FiniteDifference Jacobian")
{
  // Create the finite difference Jacobian matrix
  jacobian = createGraphAndJacobian(i, x);
}

FiniteDifference::FiniteDifference(Interface& i, const Epetra_Vector& x, const Epetra_CrsGraph& userGraph, double beta_, double alpha_) :
  map(x.Map()),
  graph(0),
  jacobian(0),
  interface(i),
  x_perturb(x),
  fo(x),
  fp(x),
  fmPtr(0),
  Jc(x),
  alpha(alpha_),
  beta(beta_),
  betaVector(0),
  betaType(Scalar),
  diffType(Forward),
  label("NOX::FiniteDifference Jacobian")
{
  // Create the finite difference Jacobian matrix directly using a
  // user supplied graph.
  jacobian = new Epetra_CrsMatrix(Copy, userGraph);
  jacobian->TransformToLocal();
}

FiniteDifference::FiniteDifference(Interface& i, const Epetra_Vector& x, const Epetra_CrsGraph& userGraph, const Epetra_Vector& beta_, double alpha_) :
  map(x.Map()),
  graph(0),
  jacobian(0),
  interface(i),
  x_perturb(x),
  fo(x),
  fp(x),
  fmPtr(0),
  Jc(x),
  alpha(alpha_),
  beta(0),
  betaVector(&beta_),
  betaType(Vector),
  diffType(Forward),
  label("NOX::FiniteDifference Jacobian")
{
  // Create the finite difference Jacobian matrix directly using a
  // user supplied graph.
  jacobian = new Epetra_CrsMatrix(Copy, userGraph);
  jacobian->TransformToLocal();
}

FiniteDifference::~FiniteDifference()
{
  delete fmPtr;
  delete jacobian;
  delete graph;
}

const char* FiniteDifference::Label () const
{
  return label.c_str();
}

int FiniteDifference::SetUseTranspose(bool UseTranspose)
{
  return jacobian->SetUseTranspose(UseTranspose);
}

int FiniteDifference::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return jacobian->Apply(X, Y);
}

int FiniteDifference::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return jacobian->ApplyInverse(X, Y);
}

bool FiniteDifference::UseTranspose() const
{
  return jacobian->UseTranspose();
}

bool FiniteDifference::HasNormInf() const
{
  return jacobian->HasNormInf();
}

const Epetra_Map& FiniteDifference::OperatorDomainMap() const
{
  return jacobian->OperatorDomainMap();
}

const Epetra_Map& FiniteDifference::OperatorRangeMap() const
{
  return jacobian->OperatorRangeMap();
}

bool FiniteDifference::Filled() const
{
  return jacobian->Filled();
}

int FiniteDifference::NumMyRowEntries(int MyRow, int & NumEntries) const
{
  return jacobian->NumMyRowEntries(MyRow, NumEntries);
}

int FiniteDifference::MaxNumEntries() const
{
  return jacobian->MaxNumEntries();
}

int FiniteDifference::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const
{
  return jacobian->ExtractMyRowCopy(MyRow, Length, NumEntries, Values, Indices);
}

int FiniteDifference::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  return jacobian->ExtractDiagonalCopy(Diagonal);
}

int FiniteDifference::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return jacobian->Multiply(TransA, X, Y);
}

int FiniteDifference::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X,  Epetra_MultiVector& Y) const
{
  return jacobian->Solve(Upper, Trans, UnitDiagonal, X, Y);
}

int FiniteDifference::InvRowSums(Epetra_Vector& x) const
{
  return jacobian->InvRowSums(x);
}

int FiniteDifference::LeftScale(const Epetra_Vector& x)
{
  return jacobian->LeftScale(x);
}

int FiniteDifference::InvColSums(Epetra_Vector& x) const
{
  return jacobian->InvColSums(x);
}

int FiniteDifference::RightScale(const Epetra_Vector& x)
{
  return jacobian->RightScale(x);
}

double FiniteDifference::NormInf() const
{
  return jacobian->NormInf();
}

double FiniteDifference::NormOne() const
{
  return jacobian->NormOne();
}

int FiniteDifference::NumGlobalNonzeros() const
{
  return jacobian->NumGlobalNonzeros();
}

int FiniteDifference::NumGlobalRows() const
{
  return jacobian->NumGlobalRows();
}

int FiniteDifference::NumGlobalCols() const
{
  return jacobian->NumGlobalCols();
}

int FiniteDifference::NumGlobalDiagonals() const
{
  return jacobian->NumGlobalDiagonals();
}

int FiniteDifference::NumMyNonzeros() const
{
  return jacobian->NumMyNonzeros();
}

int FiniteDifference::NumMyRows() const
{
  return jacobian->NumMyRows();
}

int FiniteDifference::NumMyCols() const
{
  return jacobian->NumMyCols();
}

int FiniteDifference::NumMyDiagonals() const
{
  return jacobian->NumMyDiagonals();
}

bool FiniteDifference::LowerTriangular() const
{
  return jacobian->LowerTriangular();
}

bool FiniteDifference::UpperTriangular() const
{
  return jacobian->UpperTriangular();
}

const Epetra_Comm& FiniteDifference::Comm() const
{
  return jacobian->Comm();
}

const Epetra_Map& FiniteDifference::RowMatrixRowMap() const
{
  return jacobian->RowMatrixRowMap();
}

const Epetra_Map& FiniteDifference::RowMatrixColMap() const
{
  return jacobian->RowMatrixColMap();
}

const Epetra_Import* FiniteDifference::RowMatrixImporter() const
{
  return jacobian->RowMatrixImporter();
}

const Epetra_BlockMap& FiniteDifference::Map() const
{
  return jacobian->Map();
}

bool FiniteDifference::computeJacobian(const Epetra_Vector& x)
{
  return( computeJacobian(x, *this));
}

bool FiniteDifference::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  // First check to make sure Jac is a NOX::Epetra::FiniteDifference object
  FiniteDifference* testMatrix = dynamic_cast<FiniteDifference*>(&Jac);
  if (testMatrix == 0) {
    cout << "ERROR: NOX::Epetra::FiniteDifference::computeJacobian() - "
	 << "Jacobian to evaluate is not a FiniteDifference object!" << endl;
    throw "NOX Error";
  }

  // We need the Epetra_CrsMatrix inside the FiniteDifference object
  // for the correct insertion commands.
  Epetra_CrsMatrix& jac = *testMatrix->jacobian;

  // Zero out Jacobian
  jac.PutScalar(0.0);

  // Create an extra perturbed residual vector pointer if needed
  if ( diffType == Centered )
    if ( !fmPtr )
      fmPtr = new Epetra_Vector(x);

  // Create a reference to the extra perturbed residual vector
  Epetra_Vector& fm = *fmPtr;

  double scaleFactor = 1.0;
  if ( diffType == Backward )
    scaleFactor = -1.0;

  double eta = 0.0;  // Value to perturb the solution vector

  int min = map.MinAllGID();  // Minimum Global ID value
  int max = map.MaxAllGID();  // Maximum Global ID value
  int myMin = map.MinMyGID(); // Minimum Local ID value
  int myMax = map.MaxMyGID(); // Maximum Local ID value

  // Compute the RHS at the initial solution
  interface.computeF(x, fo, Interface::Jacobian);

  x_perturb = x;

  // loop over each global unknown
  for (int k = min; k < max+1; k++) {

    // Perturb the solution vector via local IDs only if it is owned by
    // this processor
    int proc = 0;
    if (map.MyGID(k)) {

      if (betaType == Scalar)
	eta = alpha*x[map.LID(k)] + beta;
      else
	eta = alpha*x[map.LID(k)] + (*betaVector)[map.LID(k)];

      x_perturb[map.LID(k)] += scaleFactor * eta;
      proc = map.Comm().MyPID();
    }

    // Find what proc eta is on
    int broadcastProc = 0;
    map.Comm().SumAll(&proc ,&broadcastProc , 1);
    // Send the perturbation variable, eta, to all processors
    map.Comm().Broadcast(&eta, 1, broadcastProc);

    // Compute the perturbed RHS
    interface.computeF(x_perturb,fp, Interface::Jacobian);

    if ( diffType == Centered ) {
      if (map.MyGID(k))
        x_perturb[map.LID(k)] -= 2.0 * eta;
      interface.computeF(x_perturb,fm, Interface::Jacobian);
    }

    // Compute the column k of the Jacobian
    if ( diffType != Centered ) {
      Jc.Update(1.0, fp, -1.0, fo, 0.0);
      Jc.Scale( 1.0/(scaleFactor * eta) );
    }
    else {
      Jc.Update(1.0, fp, -1.0, fm, 0.0);
      Jc.Scale( 1.0/(2.0 * eta) );
    }

    // Insert nonzero column entries into the jacobian
    for (int j = myMin; j < myMax+1; j++) {
      if (Jc[map.LID(j)] != 0.0) {
	jac.ReplaceGlobalValues(j,1,&Jc[map.LID(j)],&k);
      }
    }

    // Unperturb the solution vector
    x_perturb = x;

  }

  jac.TransformToLocal();

  return true;
}

bool FiniteDifference::computePreconditioner(const Epetra_Vector& x, Epetra_RowMatrix& M)
{
  return computeJacobian(x, M);
}

Epetra_CrsMatrix*  FiniteDifference::createGraphAndJacobian(Interface& i,
						  const Epetra_Vector& x)
{

  double eta = 0.0;  // Value to perturb the solution vector

  int min = map.MinAllGID();  // Minimum Global ID value
  int max = map.MaxAllGID();  // Maximum Global ID value
  int myMin = map.MinMyGID(); // Minimum Local ID value
  int myMax = map.MaxMyGID(); // Maximum Local ID value

  // Create the graph
  graph = new Epetra_CrsGraph(Copy,map,10);

  // Compute the RHS at the initial solution
  i.computeF(x,fo);

  // loop over each global unknown
  for (int k = min; k < max+1; k++) {

    // Perturb the solution vector via local IDs only if it is owned by
    // this processor
    if (map.MyGID(k)) {

      if (betaType == Scalar)
	eta = alpha*x[map.LID(k)] + beta;
      else
	eta = alpha*x[map.LID(k)] + (*betaVector)[map.LID(k)];

      x_perturb[map.LID(k)] += eta;
    }

    // Compute the perturbed RHS
    i.computeF(x_perturb,fp);

    // Compute the column k of the Jacobian
    Jc.Update(1.0, fp, -1.0, fo, 0.0);
    //Jc.Scale(1.0/eta);

    // Insert column entries into the graph
    for (int j = myMin; j < myMax+1; j++) {
      if (Jc[map.LID(j)] != 0.0) {
	graph->InsertGlobalIndices(j,1,&k);
      }
    }

    // Unperturb the solution vector
    x_perturb = x;

  }

  graph->TransformToLocal();
   // No-ops because graph is already sorted (maherou) graph->SortIndices();
   // No-ops because graph is already sorted (maherou) graph->RemoveRedundantIndices();
  jacobian = new Epetra_CrsMatrix(Copy, *graph);
  jacobian->TransformToLocal();

  return jacobian;

}

void FiniteDifference::setDifferenceMethod(DifferenceType diffType_)
{
  diffType = diffType_;
}

Epetra_CrsMatrix& FiniteDifference::getUnderlyingMatrix() const
{
  return *jacobian;
}

void FiniteDifference::Print(ostream& strm) const
{
  jacobian->Print(strm);
}
