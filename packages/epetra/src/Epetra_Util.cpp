
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


#include "Epetra_Util.h"
#include "Epetra_Object.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"

const double Epetra_Util::chopVal_ = 1.0e-15;

//=============================================================================
  void Epetra_Util::Sort(bool SortAscending, int NumKeys, int * Keys, 
			int NumDoubleCompanions,double ** DoubleCompanions, 
			int NumIntCompanions, int ** IntCompanions) const {

  int i;

  int n = NumKeys;
  int * const list = Keys;
  int m = n/2;
  
  while (m > 0) {
    int max = n - m;
    for (int j=0; j<max; j++)
      {
	for (int k=j; k>=0; k-=m)
	  {
	    if ((SortAscending && list[k+m] >= list[k]) || 
		( !SortAscending && list[k+m] >= list[k]))
	      break;
	    int temp = list[k+m];
	    list[k+m] = list[k];
	    list[k] = temp;
	    for (i=0; i<NumDoubleCompanions; i++) {
	      double dtemp = DoubleCompanions[i][k+m];
	    DoubleCompanions[i][k+m] = DoubleCompanions[i][k];
	    DoubleCompanions[i][k] = dtemp;
	    }
	    for (i=0; i<NumIntCompanions; i++) {
	      int itemp = IntCompanions[i][k+m];
	    IntCompanions[i][k+m] = IntCompanions[i][k];
	    IntCompanions[i][k] = itemp;
	    }
	  }
      }
    m = m/2;
  }

}

//----------------------------------------------------------------------------
int Epetra_Util_binary_search(int item,
                              const int* list,
                              int len,
                              int& insertPoint)
{
  if (len < 2) {
    if (len < 1) {
      insertPoint = 0;
      return(-1);
    }

    if (list[0] == item) return(0);
    else {
      insertPoint = list[0] < item ? 1 : 0;
      return(-1);
    }
  }

  unsigned start = 0, end = len - 1;

  while(end - start > 1) {
    unsigned mid = (start + end) >> 1;
    if (list[mid] < item) start = mid;
    else end = mid;
  }

  if (list[end] < item) {
    insertPoint = (int)end+1;
    return(-1);
  }

  if (list[start] == item) return((int)start);
  if (list[end] == item) return((int)end);

  if (list[start] < item) insertPoint = (int)end;
  else insertPoint = (int)start;

  return(-1);
}

//=========================================================================
int Epetra_Util_ExtractHbData(Epetra_CrsMatrix * A, Epetra_MultiVector * LHS,
			      Epetra_MultiVector * RHS,
			      int & M, int & N, int & nz, int * & ptr,
			      int * & ind, double * & val, int & Nrhs,
			      double * & rhs, int & ldrhs,
			      double * & lhs, int & ldlhs) {

  int ierr = 0;
  if (A==0) EPETRA_CHK_ERR(-1); // This matrix is defined
  if (!A->IndicesAreContiguous()) { // Data must be contiguous for this to work
    EPETRA_CHK_ERR(A->MakeDataContiguous()); // Call MakeDataContiguous() method on the matrix
    ierr = 1; // Warn User that we changed the matrix
  }
  
  M = A->NumMyRows();
  N = A->NumMyCols();
  nz = A->NumMyNonzeros();
  val = (*A)[0];        // Dangerous, but cheap and effective way to access first element in 
  
  const Epetra_CrsGraph & Graph = A->Graph();
  ind = Graph[0];  // list of values and indices
  
  Nrhs = 0; // Assume no rhs, lhs

  if (RHS!=0) {
    Nrhs = RHS->NumVectors();
    if (Nrhs>1)
    if (!RHS->ConstantStride()) {EPETRA_CHK_ERR(-2)}; // Must have strided vectors
    ldrhs = RHS->Stride();
    rhs = (*RHS)[0]; // Dangerous but effective (again)
  }
  if (LHS!=0) {
    int Nlhs = LHS->NumVectors();
    if (Nlhs!=Nrhs) {EPETRA_CHK_ERR(-3)}; // Must have same number of rhs and lhs
    if (Nlhs>1)
    if (!LHS->ConstantStride()) {EPETRA_CHK_ERR(-4)}; // Must have strided vectors
  ldlhs = LHS->Stride();
  lhs = (*LHS)[0];
  }
  
  // Finally build ptr vector
  
  if (ptr==0) {
    ptr = new int[M+1];
    ptr[0] = 0;
    for (int i=0; i<M; i++) ptr[i+1] = ptr[i] + Graph.NumMyIndices(i);
  }
  EPETRA_CHK_ERR(ierr);
  return(0);
}
