
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_Util.h"
#include "Epetra_Object.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"

const double Epetra_Util::chopVal_ = 1.0e-15;

//=========================================================================
unsigned int Epetra_Util::RandomInt() {

  const int a = 16807;
	const int m = 2147483647;
	const int q = 127773;
	const int r = 2836;

	int hi = Seed_ / q;
	int lo = Seed_ % q;
	int test = a * lo - r * hi;
	if (test > 0)
		Seed_ = test;
	else
		Seed_ = test + m;
	
	return(Seed_);
}

//=========================================================================
double Epetra_Util::RandomDouble() {
	const double Modulus = 2147483647.0;
	const double DbleOne = 1.0;
	const double DbleTwo = 2.0;

	double randdouble = RandomInt(); // implicit conversion from int to double
	randdouble = DbleTwo * (randdouble / Modulus) - DbleOne; // scale to (-1.0, 1.0)

	return(randdouble);
}

//=========================================================================
unsigned int Epetra_Util::Seed() const {
	return(Seed_);
}

//=========================================================================
int Epetra_Util::SetSeed(unsigned int Seed) {
	Seed_ = Seed;
	return(0);
}

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

  if (list[start] == item) return(start);
  if (list[end] == item) return(end);

  if (list[end] < item) {
    insertPoint = end+1;
    return(-1);
  }

  if (list[start] < item) insertPoint = end;
  else insertPoint = start;

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
