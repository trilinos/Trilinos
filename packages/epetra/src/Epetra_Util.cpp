
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_Util.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_Directory.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"

const double Epetra_Util::chopVal_ = 1.0e-15;

//=========================================================================
double Epetra_Util::Chop(const double & Value){
  if (std::abs(Value) < chopVal_) return 0;
  return Value;
}

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
int Epetra_Util::SetSeed(unsigned int Seed_in) {
	Seed_ = Seed_in;
	return(0);
}

//=============================================================================
void Epetra_Util::Sort(bool SortAscending, int NumKeys, int * Keys, 
		       int NumDoubleCompanions,double ** DoubleCompanions, 
		       int NumIntCompanions, int ** IntCompanions)
{
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
		( !SortAscending && list[k+m] <= list[k]))
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
Epetra_Map
Epetra_Util::Create_OneToOne_Map(const Epetra_Map& usermap,
				 bool high_rank_proc_owns_shared)
{
  //if usermap is already 1-to-1 then we'll just return a copy of it.
  if (usermap.IsOneToOne()) {
    Epetra_Map newmap(usermap);
    return(newmap);
  }

  int myPID = usermap.Comm().MyPID();
  Epetra_Directory* directory = usermap.Comm().CreateDirectory(usermap);

  int numMyElems = usermap.NumMyElements();
  const int* myElems = usermap.MyGlobalElements();

  int* owner_procs = new int[numMyElems];

  directory->GetDirectoryEntries(usermap, numMyElems, myElems, owner_procs,
				 0, 0, high_rank_proc_owns_shared);

  //we'll fill a list of map-elements which belong on this processor

  int* myOwnedElems = new int[numMyElems];
  int numMyOwnedElems = 0;

  for(int i=0; i<numMyElems; ++i) {
    int GID = myElems[i];
    int owner = owner_procs[i];

    if (myPID == owner) {
      myOwnedElems[numMyOwnedElems++] = GID;
    }
  }

  Epetra_Map one_to_one_map(-1, numMyOwnedElems, myOwnedElems,
			 usermap.IndexBase(), usermap.Comm());

  delete [] myOwnedElems;
  delete [] owner_procs;
  delete directory;

  return(one_to_one_map);
}

//----------------------------------------------------------------------------
Epetra_Map
Epetra_Util::Create_Root_Map(const Epetra_Map& usermap,
				 int root)
{

  int numProc = usermap.Comm().NumProc();
  if (numProc==1) {
    Epetra_Map newmap(usermap);
    return(newmap);
  }

  const Epetra_Comm & comm = usermap.Comm();
  bool isRoot = usermap.Comm().MyPID()==root;

  //if usermap is already completely owned by root then we'll just return a copy of it.
  int quickreturn = 0;
  int globalquickreturn = 0;

  if (isRoot) {
    if (usermap.NumMyElements()==usermap.NumGlobalElements()) quickreturn = 1;
  }
  else {
    if (usermap.NumMyElements()==0) quickreturn = 1;
  }
  usermap.Comm().MinAll(&quickreturn, &globalquickreturn, 1);
  
  if (globalquickreturn==1) {
    Epetra_Map newmap(usermap);
    return(newmap);
  }
  
  // Linear map: Simple case, just put all GIDs linearly on root processor
  if (usermap.LinearMap() && root!=-1) {
    int numMyElements = 0;
    if (isRoot) numMyElements = usermap.MaxAllGID()+1;
    Epetra_Map newmap(-1, numMyElements, usermap.IndexBase(), comm);
    return(newmap);
  }

  if (!usermap.UniqueGIDs()) 
    throw usermap.ReportError("usermap must have unique GIDs",-1);

  // General map

  // Build IntVector of the GIDs, then ship them to root processor
  int numMyElements = usermap.NumMyElements();
  Epetra_Map allGidsMap(-1, numMyElements, 0, comm);
  Epetra_IntVector allGids(allGidsMap);
  for (int i=0; i<numMyElements; i++) allGids[i] = usermap.GID(i);
  
  int numGlobalElements = usermap.NumGlobalElements();
  if (root!=-1) {
    int n1 = 0; if (isRoot) n1 = numGlobalElements;
    Epetra_Map allGidsOnRootMap(-1, n1, 0, comm);
    Epetra_Import importer(allGidsOnRootMap, allGidsMap);
    Epetra_IntVector allGidsOnRoot(allGidsOnRootMap);
    allGidsOnRoot.Import(allGids, importer, Insert);
    
    Epetra_Map rootMap(-1, allGidsOnRoot.MyLength(), allGidsOnRoot.Values(), usermap.IndexBase(), comm);
    return(rootMap);
  }
  else {
    int n1 = numGlobalElements;
    Epetra_LocalMap allGidsOnRootMap(n1, 0, comm);
    Epetra_Import importer(allGidsOnRootMap, allGidsMap);
    Epetra_IntVector allGidsOnRoot(allGidsOnRootMap);
    allGidsOnRoot.Import(allGids, importer, Insert);
    
    Epetra_Map rootMap(-1, allGidsOnRoot.MyLength(), allGidsOnRoot.Values(), usermap.IndexBase(), comm);

    return(rootMap);
  }
}

//----------------------------------------------------------------------------
Epetra_BlockMap
Epetra_Util::Create_OneToOne_BlockMap(const Epetra_BlockMap& usermap,
				      bool high_rank_proc_owns_shared)
{
  //if usermap is already 1-to-1 then we'll just return a copy of it.
  if (usermap.IsOneToOne()) {
    Epetra_BlockMap newmap(usermap);
    return(newmap);
  }

  int myPID = usermap.Comm().MyPID();
  Epetra_Directory* directory = usermap.Comm().CreateDirectory(usermap);

  int numMyElems = usermap.NumMyElements();
  const int* myElems = usermap.MyGlobalElements();

  int* owner_procs = new int[numMyElems*2];
  int* sizes = owner_procs+numMyElems;

  directory->GetDirectoryEntries(usermap, numMyElems, myElems, owner_procs,
				 0, sizes, high_rank_proc_owns_shared);

  //we'll fill a list of map-elements which belong on this processor

  int* myOwnedElems = new int[numMyElems*2];
  int* ownedSizes = myOwnedElems+numMyElems;
  int numMyOwnedElems = 0;

  for(int i=0; i<numMyElems; ++i) {
    int GID = myElems[i];
    int owner = owner_procs[i];

    if (myPID == owner) {
      ownedSizes[numMyOwnedElems] = sizes[i];
      myOwnedElems[numMyOwnedElems++] = GID;
    }
  }

  Epetra_BlockMap one_to_one_map(-1, numMyOwnedElems, myOwnedElems,
				 sizes, usermap.IndexBase(), usermap.Comm());

  delete [] myOwnedElems;
  delete [] owner_procs;
  delete directory;

  return(one_to_one_map);
}

//----------------------------------------------------------------------------
int Epetra_Util_binary_search(int item,
                              const int* list,
                              int len,
                              int& insertPoint)
{
  if (len < 1) {
    insertPoint = 0;
    return(-1);
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
