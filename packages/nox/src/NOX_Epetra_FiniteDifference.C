
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
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "NOX_Epetra_Interface.H"
#include "NOX_Utils.H"

#include "NOX_Epetra_FiniteDifference.H"

using namespace NOX;
using namespace NOX::Epetra;

FiniteDifference::FiniteDifference(Interface& i, const Epetra_Vector& x) :
  map(x.Map()),
  interface(i),
  x_perturb(x),
  fo(x),
  fp(x),
  Jc(x)
{
  // Create the finite difference Jacobian matrix
  jacobian = createJacobian(i, x);
}

FiniteDifference::~FiniteDifference()
{
  delete jacobian;
  delete graph;
}

char* FiniteDifference::Label () const
{
  return jacobian->Label();
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

const Epetra_BlockMap& FiniteDifference::DomainMap() const
{
  return jacobian->DomainMap();
}

const Epetra_BlockMap& FiniteDifference::RangeMap() const
{
  return jacobian->RangeMap();
}

bool FiniteDifference::Filled() const
{
  return jacobian->Filled();
}

int FiniteDifference::NumMyRowEntries(int MyRow, int & NumEntries) const
{
  return jacobian->NumMyRowEntries(MyRow, NumEntries);
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

const Epetra_BlockMap& FiniteDifference::BlockRowMap() const
{
  return jacobian->BlockRowMap();
}

const Epetra_BlockMap& FiniteDifference::BlockImportMap() const
{
  return jacobian->BlockImportMap();
}
  
const Epetra_Import* FiniteDifference::Importer() const
{
  return jacobian->Importer();
}

bool FiniteDifference::computeJacobian(const Epetra_Vector& x, Epetra_RowMatrix& Jac)
{
  // First check to make sure Jac is a NOX::Epetra::FiniteDifference object
  FiniteDifference* testMatrix = dynamic_cast<FiniteDifference*>(&Jac);
  if (testMatrix == 0) {
    cout << "ERROR: NOX::Epetra::FiniteDifference::computeJacobian() - "
	 << "Jacobian to evaluate is not a FiniteDifference object!" << endl;
    throw 1;
  } 

  // We need the Epetra_CrsMatrix inside the FiniteDifference object 
  // for the correct insertion commands.
  Epetra_CrsMatrix& jac = *testMatrix->jacobian;

  // Zero out Jacobian
  jacobian->PutScalar(0.0);

  double eta = 0.0;  // Value to perturb the solution vector 
  
  int min = map.MinAllGID();  // Minimum Global ID value
  int max = map.MaxAllGID();  // Maximum Global ID value
  int myMin = map.MinMyGID(); // Minimum Local ID value
  int myMax = map.MaxMyGID(); // Maximum Local ID value

  // Compute the RHS at the initial solution
  interface.computeRHS(x,fo);

  x_perturb = x;

  // loop over each global unknown
  for (int k = min; k < max+1; k++) {
  
    // Perturb the solution vector via local IDs only if it is owned by 
    // this processor
    int proc = 0;
    if (map.MyGID(k)) {
      eta = 1.0e-4*x[map.LID(k)] + 1.0e-4;
      x_perturb[map.LID(k)] += eta;
      proc = map.Comm().MyPID();
    }  

    // Find what proc eta is on
    int broadcastProc = 0;
    map.Comm().SumAll(&proc ,&broadcastProc , 1); 
    // Send the perturbation variable, eta, to all processors
    map.Comm().Broadcast(&eta, 1, broadcastProc); 

    // Compute the perturbed RHS
    interface.computeRHS(x_perturb,fp);
    
    // Compute the column k of the Jacobian
    Jc.Update(1.0, fp, -1.0, fo, 0.0);
    Jc.Scale(1.0/eta);

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

Epetra_CrsMatrix*  FiniteDifference::createJacobian(Interface& i, 
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
  i.computeRHS(x,fo);

  // loop over each global unknown
  for (int k = min; k < max+1; k++) {
  
    // Perturb the solution vector via local IDs only if it is owned by 
    // this processor
    if (map.MyGID(k)) {
      eta = 1.0e-4*x[map.LID(k)] + 1.0e-6;
      x_perturb[map.LID(k)] += eta;
    }  

    // Compute the perturbed RHS
    i.computeRHS(x_perturb,fp);
    
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
  graph->SortIndices();
  graph->RemoveRedundantIndices();
  jacobian = new Epetra_CrsMatrix(Copy, *graph);
  jacobian->TransformToLocal();

  return jacobian;

}
