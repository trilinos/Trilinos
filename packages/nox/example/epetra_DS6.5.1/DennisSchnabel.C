#include <iostream>
#include <cmath>
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

#include "DennisSchnabel.H"

// Constructor - creates the Epetra objects (maps and vectors) 
DennisSchnabel::DennisSchnabel(int numGlobalElements, Epetra_Comm& comm) :
  NumGlobalElements(numGlobalElements),
  Comm(&comm)
{

  // Commonly used variables
  int i;
  MyPID = Comm->MyPID();      // Process ID
  NumProc = Comm->NumProc();  // Total number of processes

  // Construct a Source Map that puts approximately the same 
  // Number of equations on each processor in uniform global ordering
  StandardMap = new Epetra_Map(NumGlobalElements, 0, *Comm);

  // Get the number of elements owned by this processor
  NumMyElements = StandardMap->NumMyElements();

  // Construct an overlaped map for the finite element fill *************
  // For single processor jobs, the overlap and standard map are the same
  if (NumProc == 1) {
    OverlapMap = new Epetra_Map(*StandardMap);
  } else {

    int OverlapNumMyElements;
    int OverlapMinMyGID;
    
    OverlapNumMyElements = 2;    
    int OverlapMyGlobalElements[OverlapNumMyElements];
    
    for (i = 0; i < OverlapNumMyElements; i ++) 
      OverlapMyGlobalElements[i] = i;
    
    OverlapMap = new Epetra_Map(-1, OverlapNumMyElements, 
			    OverlapMyGlobalElements, 0, *Comm);
  } // End Overlap map construction *************************************

  // Construct Linear Objects  
  Importer = new Epetra_Import(*OverlapMap, *StandardMap);
  initialSolution = new Epetra_Vector(*StandardMap);
  AA = new Epetra_CrsGraph(Copy, *StandardMap, 5);

  // Allocate the memory for a matrix dynamically (i.e. the graph is dynamic).
  generateGraph(*AA);

  // Create a second matrix using graph of first matrix - this creates a 
  // static graph so we can refill the new matirx after TransformToLocal()
  // is called.
  A = new Epetra_CrsMatrix (Copy, *AA);
  A->TransformToLocal();
}

// Destructor
DennisSchnabel::~DennisSchnabel()
{
  delete AA;
  delete A;
  delete initialSolution;
  delete Importer;
  delete OverlapMap;
  delete StandardMap;
}

// Matrix and Residual Fills
void DennisSchnabel::evaluate(FillType f, 
				    const Epetra_Vector* tmp_soln, 
				    Epetra_Vector* tmp_rhs, 
				    Epetra_RowMatrix* tmp_matrix)
{
  flag = f;

  // Set the incoming linear objects
  if (flag == RHS_ONLY) {
    soln = const_cast<Epetra_Vector*>(tmp_soln);
    rhs = tmp_rhs;
  } else if (flag == MATRIX_ONLY) {
    soln = const_cast <Epetra_Vector*> (tmp_soln);
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
  } else if (flag == ALL) { 
    soln = const_cast<Epetra_Vector*>(tmp_soln);
    rhs = tmp_rhs;
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
  } else {
    cout << "ERROR: DennisSchnabel::fillMatrix() - FillType flag is broken" << endl;
    throw;
  }

  // Create the overlapped solution and position vectors
  Epetra_Vector u(*OverlapMap);

  // Export Solution to Overlap vector
  u.Import(*soln, *Importer, Insert)==0;

  // Declare required variables
  int i,j,ierr;
  int OverlapNumMyElements = OverlapMap->NumMyElements();
  int OverlapMinMyGID;
  if (MyPID==0) OverlapMinMyGID = StandardMap->MinMyGID();
  else OverlapMinMyGID = StandardMap->MinMyGID()-1;

  // Zero out the objects that will be filled
  if ((flag == MATRIX_ONLY) || (flag == ALL)) i=A->PutScalar(0.0);
  if ((flag == RHS_ONLY)    || (flag == ALL)) i=rhs->PutScalar(0.0);

 

  if((flag==RHS_ONLY) || (flag==ALL)) {

    if (MyPID==0) { 
      (*rhs)[0]=(u[0]*u[0] + u[1]*u[1] - 2.);
      if (NumProc==1) (*rhs)[1]=(exp(u[0]-1.) + u[1]*u[1]*u[1] - 2.);
    } else { 
      (*rhs)[0]=(exp(u[0]-1.) + u[1]*u[1]*u[1] - 2.);
    }
  }

  
  int row; 
  int* column = new int[2];
  double* jac = new double[2];

  if((flag==MATRIX_ONLY) || (flag==ALL)) {

    if (MyPID==0) {
      // fill global row 0 on proc 0
      column[0] = 0;  jac[0] = 2.*u[0];
      column[1] = 1;  jac[1] = 2.*u[1];
      ierr=A->ReplaceGlobalValues(0, 2, jac, column);
      if (NumProc==1) {
	// Fill global row 1 on proc 0 if single processor job 
	column[0] = 0;  jac[0] = exp(u[0]-1.);
	column[1] = 1;  jac[1] = 3.*u[1]*u[1];
	ierr=A->ReplaceGlobalValues(1, 2, jac, column);
      }
    } else {
      // Fill global row 1 on proc 2 if 2 processor job
      column[0] = 0;  jac[0] = exp(u[0]-1.);
      column[1] = 1;  jac[1] = 3.*u[1]*u[1];
      ierr=A->ReplaceGlobalValues(1, 2, jac, column);
    }

    delete [] column;
    delete [] jac;
  } 

  // Sync up processors to be safe
  Comm->Barrier();
 
  A->TransformToLocal();

  return ;
}

Epetra_Vector& DennisSchnabel::getSolution()
{
  return *initialSolution;
}
  
Epetra_CrsMatrix& DennisSchnabel::getJacobian()
{
  return *A;
}

Epetra_CrsGraph& DennisSchnabel::generateGraph(Epetra_CrsGraph& AA)
{
  
  int* index = new int[2];

  if (MyPID==0) {
    index[0]=0;
    index[1]=1;
    AA.InsertGlobalIndices(0, 2, index);
  
    if (NumProc==1) {
      index[0]=0;
      index[1]=1;
      AA.InsertGlobalIndices(1, 2, index);
    }
  } else {
    index[0]=0;
    index[1]=1;
    AA.InsertGlobalIndices(1, 2, index);
  }
  
  delete [] index;
  
  AA.TransformToLocal();
  AA.SortIndices();
  AA.RemoveRedundantIndices();
  return AA;
}
