/*
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
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#ifdef PETRA_MPI
#include "mpi.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Petra_Comm.h"
#include "Petra_Map.h"
#include "Petra_Time.h"
#include "Petra_RDP_CRS_Matrix.h"
 
int main(int argc, char *argv[])
{
  int ierr = 0, i, j;
  bool debug = false;

#ifdef PETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;



#ifdef PETRA_MPI
  Petra_Comm & Comm = *new Petra_Comm( MPI_COMM_WORLD );
#else
  Petra_Comm & Comm = *new Petra_Comm();
#endif


  //char tmp;
  //if (rank==0) cout << "Press any key to continue..."<< endl;
  //if (rank==0) cin >> tmp;
  //Comm.Barrier();

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
              << " is alive."<<endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && rank!=0) verbose = false;

  int NumMyPoints = 10000;
  int NumGlobalPoints = NumMyPoints*NumProc+minfn(NumProc,3);
  if (MyPID < 3) NumMyPoints++;
  int IndexBase = 0;
  bool DistributedGlobal = (NumGlobalPoints>NumMyPoints);

  // Construct a Source Map that puts approximately the same Number of equations on each processor in 
  // uniform global ordering

  Petra_Map& SourceMap = *new Petra_Map(NumGlobalPoints, NumMyPoints, 0, Comm);
  
  // Get update list and number of local equations from newly created Map
  int NumMyElements = SourceMap.NumMyElements();
  int * SourceMyGlobalElements = new int[NumMyElements];
  SourceMap.MyGlobalElements(SourceMyGlobalElements);


  // Construct a Target Map that will contain:
  //  some unchanged elements (relative to the soure map),
  //  some permuted elements
  //  some off-processor elements
  Petra_RDP_Vector & RandVec = *new Petra_RDP_Vector(SourceMap);
  RandVec.Random(); // This creates a vector of random numbers between negative one and one.

  int *TargetMyGlobalElements = new int[NumMyElements];

  for (i=0; i< NumMyPoints/2; i++) TargetMyGlobalElements[i] = i; // Half will be the same...
  for (i=NumMyPoints/2; i<NumMyPoints; i++) {
    int index = abs((int)(((double) (NumGlobalPoints-1) ) * RandVec[i]));
    TargetMyGlobalElements[i] = minfn(NumGlobalPoints-1,maxfn(0,index));
  }

  int NumSameIDs = 0;
  int NumPermutedIDs = 0;
  int NumRemoteIDs = 0;
  bool StillContiguous = true;
  for (i=0; i < NumMyPoints; i++) {
    if (SourceMyGlobalElements[i]==TargetMyGlobalElements[i] && StillContiguous)
      NumSameIDs++;
    else if (SourceMap.MyGID(TargetMyGlobalElements[i])) {
      StillContiguous = false;
      NumPermutedIDs++;
    }
    else {
      StillContiguous = false;
      NumRemoteIDs++;
    }
  }
  assert(NumMyPoints==NumSameIDs+NumPermutedIDs+NumRemoteIDs);

  Petra_Map & TargetMap = *new Petra_Map(-1, NumMyElements, TargetMyGlobalElements, 0, Comm);

  // Create a multivector whose elements are GlobalID * (column number +1)

  int NumVectors = 3;
  Petra_RDP_MultiVector & SourceMultiVector = *new Petra_RDP_MultiVector(SourceMap, NumVectors);
  for (j=0; j < NumVectors; j++)
    for (i=0; i < NumMyElements; i++)
      SourceMultiVector[j][i] = (double) SourceMyGlobalElements[i]*(j+1);

  // Create a target multivector that we will fill using an Import

  Petra_RDP_MultiVector & TargetMultiVector = *new Petra_RDP_MultiVector(TargetMap, NumVectors);

  Petra_Import & Importer = *new Petra_Import(TargetMap, SourceMap);

  assert(TargetMultiVector.Import(SourceMultiVector, Importer, Insert)==0);

  // Test Target against expected values

  for (j=0; j < NumVectors; j++)
    for (i=0; i < NumMyElements; i++) {
      if (TargetMultiVector[j][i]!= (double) TargetMyGlobalElements[i]*(j+1))
	cout << "TargetMultiVector["<<i<<"]["<<j<<"] = " << TargetMultiVector[j][i] 
	     <<  "  TargetMyGlobalElements[i]*(j+1) = " <<  TargetMyGlobalElements[i]*(j+1) << endl;
      assert(TargetMultiVector[j][i]== (double) TargetMyGlobalElements[i]*(j+1));
    }

  if (verbose) cout << "MultiVector Import using Importer Check OK" << endl << endl;


  //////////////////////////////////////////////////////////////////////////////

  // Now use Importer to do an export

  Petra_RDP_Vector & TargetVector = *new  Petra_RDP_Vector(SourceMap);
  Petra_RDP_Vector & ExpectedTarget = *new  Petra_RDP_Vector(SourceMap);
  Petra_RDP_Vector & SourceVector = *new  Petra_RDP_Vector(TargetMap);

  NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumExportIDs = Importer.NumExportIDs();
  int *PermuteToLIDs = Importer.PermuteToLIDs();
  int *PermuteFromLIDs = Importer.PermuteFromLIDs();
  int *ExportLIDs = Importer.ExportLIDs();
  int *ExportPIDs = Importer.ExportPIDs();

  for (i=0; i < NumSameIDs; i++) ExpectedTarget[i] = (double) (MyPID+1);
  for (i=0; i < NumPermuteIDs; i++) ExpectedTarget[PermuteFromLIDs[i]] = 
				      (double) (MyPID+1);
  for (i=0; i < NumExportIDs; i++) ExpectedTarget[ExportLIDs[i]] += 
				     (double) (ExportPIDs[i]+1);

  for (i=0; i < NumMyElements; i++) SourceVector[i] =  (double) (MyPID+1);

  assert(TargetVector.Export(SourceVector, Importer, Add)==0);

    for (i=0; i < NumMyElements; i++) {
      if (TargetVector[i]!= ExpectedTarget[i])
	cout <<  "     TargetVector["<<i<<"] = " << TargetVector[i] 
	     <<  "   ExpectedTarget["<<i<<"] = " <<  ExpectedTarget[i] << " on PE " << MyPID << endl;
      assert(TargetVector[i]== ExpectedTarget[i]);
    }

  if (verbose) cout << "Vector Export using Importer Check OK" << endl << endl;



  //////////////////////////////////////////////////////////////////////////////////////////
  //  Build a tridiagonal system two ways:
  //  1) From "standard" matrix view where equations are uniquely owned.
  //  2) From 1D PDE view where nodes (equations) between processors are shared and partial contributions are done
  //     in parallel, then merged together at the end of the construction process.
  //
  //////////////////////////////////////////////////////////////////////////////////////////



  // Construct a Standard Map that puts approximately the same number of equations on each processor in 
  // uniform global ordering

  Petra_Map& StandardMap = *new Petra_Map(NumGlobalPoints, NumMyPoints, 0, Comm);
  
  // Get update list and number of local equations from newly created Map
  NumMyElements = StandardMap.NumMyElements();
  int * StandardMyGlobalElements = new int[NumMyElements];
  StandardMap.MyGlobalElements(StandardMyGlobalElements);


  // Create a standard Petra_CRS_Graph

  Petra_CRS_Graph& StandardGraph = *new Petra_CRS_Graph(Copy, StandardMap, 3);
  assert(!StandardGraph.IndicesAreGlobal());
  assert(!StandardGraph.IndicesAreLocal());
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  int *Indices = new int[2];
  int NumEntries;
  
  for (i=0; i<NumMyPoints; i++)
    {
    if (StandardMyGlobalElements[i]==0)
      {
	Indices[0] = 1;
	NumEntries = 1;
      }
    else if (StandardMyGlobalElements[i] == NumGlobalPoints-1)
      {
	Indices[0] = NumGlobalPoints-2;
	NumEntries = 1;
      }
    else
      {
	Indices[0] = StandardMyGlobalElements[i]-1;
	Indices[1] = StandardMyGlobalElements[i]+1;
	NumEntries = 2;
      }
     assert(StandardGraph.InsertGlobalIndices(StandardMyGlobalElements[i], NumEntries, Indices)==0);
     assert(StandardGraph.InsertGlobalIndices(StandardMyGlobalElements[i], 1, StandardMyGlobalElements+i)==0); // Put in the diagonal entry
    }
  
  // Finish up
  assert(StandardGraph.IndicesAreGlobal());
  assert(StandardGraph.FillComplete()==0);
  assert(StandardGraph.IndicesAreLocal());
  assert(!StandardGraph.StorageOptimized());
  StandardGraph.OptimizeStorage();
  assert(StandardGraph.StorageOptimized());
  assert(!StandardGraph.UpperTriangular());
  assert(!StandardGraph.LowerTriangular());


  // Create Petra_RDP_CRS_Matrix using the just-built graph

  Petra_RDP_CRS_Matrix& StandardMatrix = *new Petra_RDP_CRS_Matrix(Copy, StandardGraph);
  assert(!StandardMatrix.IndicesAreGlobal());
  assert(StandardMatrix.IndicesAreLocal());
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  double *Values = new double[2];
  Values[0] = -1.0; Values[1] = -1.0;
  double two = 2.0;
  
  for (i=0; i<NumMyPoints; i++)
    {
    if (StandardMyGlobalElements[i]==0)
      {
	Indices[0] = 1;
	NumEntries = 1;
      }
    else if (StandardMyGlobalElements[i] == NumGlobalPoints-1)
      {
	Indices[0] = NumGlobalPoints-2;
	NumEntries = 1;
      }
    else
      {
	Indices[0] = StandardMyGlobalElements[i]-1;
	Indices[1] = StandardMyGlobalElements[i]+1;
	NumEntries = 2;
      }
     assert(StandardMatrix.ReplaceGlobalValues(StandardMyGlobalElements[i], NumEntries, Values, Indices)==0);
     assert(StandardMatrix.ReplaceGlobalValues(StandardMyGlobalElements[i], 1, &two, StandardMyGlobalElements+i)==0); // Put in the diagonal entry
    }
  
  // Finish up
  assert(StandardMatrix.IndicesAreLocal());
  assert(StandardMatrix.FillComplete()==0);
  assert(StandardMatrix.IndicesAreLocal());
  assert(StandardMatrix.StorageOptimized());
  assert(!StandardMatrix.UpperTriangular());
  assert(!StandardMatrix.LowerTriangular());

  // Construct an Overlapped Map of StandardMap that include the endpoints from two neighboring processors.

  int OverlapNumMyElements;
  int OverlapMinMyGID;

  OverlapNumMyElements = NumMyElements + 1;
  if (MyPID==0) OverlapNumMyElements--;

  if (MyPID==0) OverlapMinMyGID = StandardMap.MinMyGID();
  else OverlapMinMyGID = StandardMap.MinMyGID()-1;

  int * OverlapMyGlobalElements = new int[OverlapNumMyElements];

  for (i=0; i< OverlapNumMyElements; i++) OverlapMyGlobalElements[i] = OverlapMinMyGID + i;

  Petra_Map& OverlapMap = *new Petra_Map(-1, OverlapNumMyElements, OverlapMyGlobalElements, 0, Comm);

  // Create the Overlap Petra_Matrix

  Petra_RDP_CRS_Matrix& OverlapMatrix = *new Petra_RDP_CRS_Matrix(Copy, OverlapMap, 4);
  assert(!OverlapMatrix.IndicesAreGlobal());
  assert(!OverlapMatrix.IndicesAreLocal());
  
  // Add  matrix element one cell at a time.
  // Each cell does an incoming and outgoing flux calculation


  double pos_one = 1.0;
  double neg_one = -1.0;

  for (i=0; i<OverlapNumMyElements; i++)
    {
      int node_left = OverlapMyGlobalElements[i]-1;
      int node_center = node_left + 1;
      int node_right = node_left + 2;
      if (i>0) {
	if (node_left>-1)
	  assert(OverlapMatrix.InsertGlobalValues(node_center, 1, &neg_one, &node_left)==0);
	assert(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0);
      }
      if (i<OverlapNumMyElements-1) {
	assert(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0);
	if (node_right<NumGlobalPoints) 
	  assert(OverlapMatrix.InsertGlobalValues(node_center, 1, &neg_one, &node_right)==0);
      }
    }

  // Handle endpoints
  if (MyPID==0) {
    int node_center = 0;
    assert(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0);
  }
  if (MyPID==NumProc-1) {
    int node_center = OverlapMyGlobalElements[OverlapNumMyElements-1];
    assert(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0);
  }
    
  assert(OverlapMatrix.FillComplete()==0);

  // Make a gathered matrix from OverlapMatrix.  It should be identical to StandardMatrix

  Petra_RDP_CRS_Matrix& GatheredMatrix = *new Petra_RDP_CRS_Matrix(Copy, StandardGraph);
  Petra_Export & Exporter = *new Petra_Export(OverlapMap, StandardMap);
  assert(GatheredMatrix.Export(OverlapMatrix, Exporter, Add)==0);
  assert(GatheredMatrix.FillComplete()==0);

  // Check if entries of StandardMatrix and GatheredMatrix are identical

  int StandardNumEntries, GatheredNumEntries;
  int * StandardIndices, * GatheredIndices;
  double * StandardValues, * GatheredValues;

  int StandardNumMyNonzeros = StandardMatrix.NumMyNonzeros();
  int GatheredNumMyNonzeros = GatheredMatrix.NumMyNonzeros();
  assert(StandardNumMyNonzeros==GatheredNumMyNonzeros);

  int StandardNumMyRows = StandardMatrix.NumMyRows();
  int GatheredNumMyRows = GatheredMatrix.NumMyRows();
  assert(StandardNumMyRows==GatheredNumMyRows);

  for (i=0; i< StandardNumMyRows; i++)
    {
      assert(StandardMatrix.ExtractMyRowView(i, StandardNumEntries, StandardValues, StandardIndices)==0);
      assert(GatheredMatrix.ExtractMyRowView(i, GatheredNumEntries, GatheredValues, GatheredIndices)==0);
      assert(StandardNumEntries==GatheredNumEntries);
      for (j=0; j < StandardNumEntries; j++) {
	//if (StandardIndices[j]!=GatheredIndices[j])
	// cout << "MyPID = " << MyPID << " i = " << i << "   StandardIndices[" << j << "] = " << StandardIndices[j] 
	//      << "   GatheredIndices[" << j << "] = " << GatheredIndices[j] << endl;
	//if (StandardValues[j]!=GatheredValues[j])
	//cout << "MyPID = " << MyPID << " i = " << i << "    StandardValues[" << j << "] = " <<  StandardValues[j] 
	//     << "    GatheredValues[" << j << "] = " <<  GatheredValues[j] << endl;
	assert(StandardIndices[j]==GatheredIndices[j]);
	assert(StandardValues[j]==GatheredValues[j]);
      }
    }

  if (verbose) cout << "Matrix Export Check OK" << endl;
  // Release all objects

  delete &SourceVector;
  delete &TargetVector;
  delete &ExpectedTarget;


  delete &Importer;
  delete &SourceMap;
  delete &TargetMap;

  delete [] SourceMyGlobalElements;
  delete [] TargetMyGlobalElements;

  delete &SourceMultiVector;
  delete &TargetMultiVector;
  delete &RandVec;

  delete &Exporter;
  delete &GatheredMatrix;
  delete &OverlapMatrix;
  delete &OverlapMap;
  delete [] OverlapMyGlobalElements;

  delete &StandardMatrix;
  delete &StandardGraph;
  delete &StandardMap;
  delete [] StandardMyGlobalElements;

  delete [] Values;
  delete [] Indices;
  delete &Comm;

#ifdef PETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
