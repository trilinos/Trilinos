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

#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "../epetra_test_err.h"

int main(int argc, char *argv[])
{
  int ierr = 0, i, j, forierr = 0;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;




  //char tmp;
  //if (Comm.MyPID()==0) cout << "Press any key to continue..."<< endl;
  //if (Comm.MyPID()==0) cin >> tmp;
  //Comm.Barrier();

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
              << " is alive."<<endl;

  // Redefine verbose to only print on PE 0
  if (verbose && Comm.MyPID()!=0) verbose = false;

  int NumMyEquations = 20;
  int NumGlobalEquations = NumMyEquations*NumProc+EPETRA_MIN(NumProc,3);
  if (MyPID < 3) NumMyEquations++;
  // Construct a Source Map that puts approximately the same Number of equations on each processor in 
  // uniform global ordering

  Epetra_Map SourceMap(NumGlobalEquations, NumMyEquations, 0, Comm);
  
  // Get update list and number of local equations from newly created Map
  int NumMyElements = SourceMap.NumMyElements();
  int * SourceMyGlobalElements = new int[NumMyElements];
  SourceMap.MyGlobalElements(SourceMyGlobalElements);

  // Construct a Target Map that will contain:
  //  some unchanged elements (relative to the soure map),
  //  some permuted elements
  //  some off-processor elements
  Epetra_Vector RandVec(SourceMap);
  RandVec.Random(); // This creates a vector of random numbers between negative one and one.

  int *TargetMyGlobalElements = new int[NumMyElements];

  int MinGID = SourceMap.MinMyGID();
  for (i=0; i< NumMyEquations/2; i++) TargetMyGlobalElements[i] = i + MinGID; // Half will be the same...
  for (i=NumMyEquations/2; i<NumMyEquations; i++) {
    int index = abs((int)(((double) (NumGlobalEquations-1) ) * RandVec[i]));
    TargetMyGlobalElements[i] = EPETRA_MIN(NumGlobalEquations-1,EPETRA_MAX(0,index));
  }

  int NumSameIDs = 0;
  int NumPermutedIDs = 0;
  int NumRemoteIDs = 0;
  bool StillContiguous = true;
  for (i=0; i < NumMyEquations; i++) {
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
  EPETRA_TEST_ERR(!(NumMyEquations==NumSameIDs+NumPermutedIDs+NumRemoteIDs),ierr);

  Epetra_Map TargetMap(-1, NumMyElements, TargetMyGlobalElements, 0, Comm);

  // Create a multivector whose elements are GlobalID * (column number +1)

  int NumVectors = 3;
  Epetra_MultiVector SourceMultiVector(SourceMap, NumVectors);
  for (j=0; j < NumVectors; j++)
    for (i=0; i < NumMyElements; i++)
      SourceMultiVector[j][i] = (double) SourceMyGlobalElements[i]*(j+1);

  // Create a target multivector that we will fill using an Import

  Epetra_MultiVector TargetMultiVector(TargetMap, NumVectors);

  Epetra_Import Importer(TargetMap, SourceMap);

  EPETRA_TEST_ERR(!(TargetMultiVector.Import(SourceMultiVector, Importer, Insert)==0),ierr);

  // Test Target against expected values
  forierr = 0;
  for (j=0; j < NumVectors; j++)
    for (i=0; i < NumMyElements; i++) {
      if (TargetMultiVector[j][i]!= (double) TargetMyGlobalElements[i]*(j+1))
	cout << "TargetMultiVector["<<i<<"]["<<j<<"] = " << TargetMultiVector[j][i] 
	     <<  "  TargetMyGlobalElements[i]*(j+1) = " <<  TargetMyGlobalElements[i]*(j+1) << endl;
      forierr += !(TargetMultiVector[j][i]== (double) TargetMyGlobalElements[i]*(j+1));
    }
  EPETRA_TEST_ERR(forierr,ierr);

  if (verbose) cout << "MultiVector Import using Importer Check OK" << endl << endl;


  //////////////////////////////////////////////////////////////////////////////

  // Now use Importer to do an export

  Epetra_Vector TargetVector(SourceMap);
  Epetra_Vector ExpectedTarget(SourceMap);
  Epetra_Vector SourceVector(TargetMap);

  NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumExportIDs = Importer.NumExportIDs();
  int *PermuteFromLIDs = Importer.PermuteFromLIDs();
  int *ExportLIDs = Importer.ExportLIDs();
  int *ExportPIDs = Importer.ExportPIDs();

  for (i=0; i < NumSameIDs; i++) ExpectedTarget[i] = (double) (MyPID+1);
  for (i=0; i < NumPermuteIDs; i++) ExpectedTarget[PermuteFromLIDs[i]] = 
				      (double) (MyPID+1);
  for (i=0; i < NumExportIDs; i++) ExpectedTarget[ExportLIDs[i]] += 
				     (double) (ExportPIDs[i]+1);

  for (i=0; i < NumMyElements; i++) SourceVector[i] =  (double) (MyPID+1);

  EPETRA_TEST_ERR(!(TargetVector.Export(SourceVector, Importer, Add)==0),ierr);

  forierr = 0;
  for (i=0; i < NumMyElements; i++) {
    if (TargetVector[i]!= ExpectedTarget[i])
      cout <<  "     TargetVector["<<i<<"] = " << TargetVector[i] 
	   <<  "   ExpectedTarget["<<i<<"] = " <<  ExpectedTarget[i] << " on PE " << MyPID << endl;
    forierr += !(TargetVector[i]== ExpectedTarget[i]);
  }
  EPETRA_TEST_ERR(forierr,ierr);

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

  Epetra_Map StandardMap(NumGlobalEquations, NumMyEquations, 0, Comm);
  
  // Get update list and number of local equations from newly created Map
  NumMyElements = StandardMap.NumMyElements();
  int * StandardMyGlobalElements = new int[NumMyElements];
  StandardMap.MyGlobalElements(StandardMyGlobalElements);


  // Create a standard Epetra_CrsGraph

  Epetra_CrsGraph StandardGraph(Copy, StandardMap, 3);
  EPETRA_TEST_ERR(StandardGraph.IndicesAreGlobal(),ierr);
  EPETRA_TEST_ERR(StandardGraph.IndicesAreLocal(),ierr);
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  int *Indices = new int[2];
  int NumEntries;
  
  forierr = 0;
  for (i=0; i<NumMyEquations; i++)
    {
    if (StandardMyGlobalElements[i]==0)
      {
	Indices[0] = 1;
	NumEntries = 1;
      }
    else if (StandardMyGlobalElements[i] == NumGlobalEquations-1)
      {
	Indices[0] = NumGlobalEquations-2;
	NumEntries = 1;
      }
    else
      {
	Indices[0] = StandardMyGlobalElements[i]-1;
	Indices[1] = StandardMyGlobalElements[i]+1;
	NumEntries = 2;
      }
    forierr += !(StandardGraph.InsertGlobalIndices(StandardMyGlobalElements[i], NumEntries, Indices)==0);
    forierr += !(StandardGraph.InsertGlobalIndices(StandardMyGlobalElements[i], 1, StandardMyGlobalElements+i)==0); // Put in the diagonal entry
    }
  EPETRA_TEST_ERR(forierr,ierr);

  // Finish up
  EPETRA_TEST_ERR(!(StandardGraph.IndicesAreGlobal()),ierr);
  EPETRA_TEST_ERR(!(StandardGraph.TransformToLocal()==0),ierr);
  EPETRA_TEST_ERR(!(StandardGraph.IndicesAreLocal()),ierr);
  EPETRA_TEST_ERR(StandardGraph.StorageOptimized(),ierr);
  StandardGraph.OptimizeStorage();
  EPETRA_TEST_ERR(!(StandardGraph.StorageOptimized()),ierr);
  EPETRA_TEST_ERR(StandardGraph.UpperTriangular(),ierr);
  EPETRA_TEST_ERR(StandardGraph.LowerTriangular(),ierr);

  // Create Epetra_CrsMatrix using the just-built graph

  Epetra_CrsMatrix StandardMatrix(Copy, StandardGraph);
  EPETRA_TEST_ERR(StandardMatrix.IndicesAreGlobal(),ierr);
  EPETRA_TEST_ERR(!(StandardMatrix.IndicesAreLocal()),ierr);
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  double *Values = new double[2];
  Values[0] = -1.0; Values[1] = -1.0;
  double two = 2.0;
  
  forierr = 0;
  for (i=0; i<NumMyEquations; i++)
    {
    if (StandardMyGlobalElements[i]==0)
      {
	Indices[0] = 1;
	NumEntries = 1;
      }
    else if (StandardMyGlobalElements[i] == NumGlobalEquations-1)
      {
	Indices[0] = NumGlobalEquations-2;
	NumEntries = 1;
      }
    else
      {
	Indices[0] = StandardMyGlobalElements[i]-1;
	Indices[1] = StandardMyGlobalElements[i]+1;
	NumEntries = 2;
      }
    forierr += !(StandardMatrix.ReplaceGlobalValues(StandardMyGlobalElements[i], NumEntries, Values, Indices)==0);
    // Put in the diagonal entry
    forierr += !(StandardMatrix.ReplaceGlobalValues(StandardMyGlobalElements[i], 1, &two, StandardMyGlobalElements+i)==0); 
    }
  EPETRA_TEST_ERR(forierr,ierr);

  // Finish up
  EPETRA_TEST_ERR(!(StandardMatrix.IndicesAreLocal()),ierr);
  EPETRA_TEST_ERR(!(StandardMatrix.TransformToLocal()==0),ierr);
  EPETRA_TEST_ERR(!(StandardMatrix.IndicesAreLocal()),ierr);
  EPETRA_TEST_ERR(!(StandardMatrix.StorageOptimized()),ierr);
  EPETRA_TEST_ERR(StandardMatrix.UpperTriangular(),ierr);
  EPETRA_TEST_ERR(StandardMatrix.LowerTriangular(),ierr);

  // Construct an Overlapped Map of StandardMap that include the endpoints from two neighboring processors.

  int OverlapNumMyElements;
  int OverlapMinMyGID;

  OverlapNumMyElements = NumMyElements + 1;
  if (MyPID==0) OverlapNumMyElements--;

  if (MyPID==0) OverlapMinMyGID = StandardMap.MinMyGID();
  else OverlapMinMyGID = StandardMap.MinMyGID()-1;

  int * OverlapMyGlobalElements = new int[OverlapNumMyElements];

  for (i=0; i< OverlapNumMyElements; i++) OverlapMyGlobalElements[i] = OverlapMinMyGID + i;

  Epetra_Map OverlapMap(-1, OverlapNumMyElements, OverlapMyGlobalElements, 0, Comm);

  // Create the Overlap Epetra_Matrix

  Epetra_CrsMatrix OverlapMatrix(Copy, OverlapMap, 4);
  EPETRA_TEST_ERR(OverlapMatrix.IndicesAreGlobal(),ierr);
  EPETRA_TEST_ERR(OverlapMatrix.IndicesAreLocal(),ierr);
  
  // Add  matrix element one cell at a time.
  // Each cell does an incoming and outgoing flux calculation


  double pos_one = 1.0;
  double neg_one = -1.0;

  forierr = 0;
  for (i=0; i<OverlapNumMyElements; i++)
    {
      int node_left = OverlapMyGlobalElements[i]-1;
      int node_center = node_left + 1;
      int node_right = node_left + 2;
      if (i>0) {
	if (node_left>-1)
	  forierr += !(OverlapMatrix.InsertGlobalValues(node_center, 1, &neg_one, &node_left)==0);
	forierr += !(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0);
      }
      if (i<OverlapNumMyElements-1) {
	forierr += !(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0);
	if (node_right<NumGlobalEquations) 
	  forierr += !(OverlapMatrix.InsertGlobalValues(node_center, 1, &neg_one, &node_right)==0);
      }
    }
  EPETRA_TEST_ERR(forierr,ierr);

  // Handle endpoints
  if (MyPID==0) {
    int node_center = 0;
    EPETRA_TEST_ERR(!(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0),ierr);
  }
  if (MyPID==NumProc-1) {
    int node_center = OverlapMyGlobalElements[OverlapNumMyElements-1];
    EPETRA_TEST_ERR(!(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0),ierr);
  }
    
  EPETRA_TEST_ERR(!(OverlapMatrix.TransformToLocal()==0),ierr);

  // Make a gathered matrix from OverlapMatrix.  It should be identical to StandardMatrix

  Epetra_CrsMatrix GatheredMatrix(Copy, StandardGraph);
  Epetra_Export Exporter(OverlapMap, StandardMap);
  EPETRA_TEST_ERR(!(GatheredMatrix.Export(OverlapMatrix, Exporter, Add)==0),ierr);
  EPETRA_TEST_ERR(!(GatheredMatrix.TransformToLocal()==0),ierr);

  // Check if entries of StandardMatrix and GatheredMatrix are identical

  int StandardNumEntries, GatheredNumEntries;
  int * StandardIndices, * GatheredIndices;
  double * StandardValues, * GatheredValues;

  int StandardNumMyNonzeros = StandardMatrix.NumMyNonzeros();
  int GatheredNumMyNonzeros = GatheredMatrix.NumMyNonzeros();
  EPETRA_TEST_ERR(!(StandardNumMyNonzeros==GatheredNumMyNonzeros),ierr);

  int StandardNumMyRows = StandardMatrix.NumMyRows();
  int GatheredNumMyRows = GatheredMatrix.NumMyRows();
  EPETRA_TEST_ERR(!(StandardNumMyRows==GatheredNumMyRows),ierr);

  forierr = 0;
  for (i=0; i< StandardNumMyRows; i++)
    {
      forierr += !(StandardMatrix.ExtractMyRowView(i, StandardNumEntries, StandardValues, StandardIndices)==0);
      forierr += !(GatheredMatrix.ExtractMyRowView(i, GatheredNumEntries, GatheredValues, GatheredIndices)==0);
      forierr += !(StandardNumEntries==GatheredNumEntries);
      for (j=0; j < StandardNumEntries; j++) {
	//if (StandardIndices[j]!=GatheredIndices[j])
	// cout << "MyPID = " << MyPID << " i = " << i << "   StandardIndices[" << j << "] = " << StandardIndices[j] 
	//      << "   GatheredIndices[" << j << "] = " << GatheredIndices[j] << endl;
	//if (StandardValues[j]!=GatheredValues[j])
	//cout << "MyPID = " << MyPID << " i = " << i << "    StandardValues[" << j << "] = " <<  StandardValues[j] 
	//     << "    GatheredValues[" << j << "] = " <<  GatheredValues[j] << endl;
	forierr += !(StandardIndices[j]==GatheredIndices[j]);
	forierr += !(StandardValues[j]==GatheredValues[j]);
      }
    }
  EPETRA_TEST_ERR(forierr,ierr);

  if (verbose) cout << "Matrix Export Check OK" << endl << endl;

  bool passed;
  Epetra_IntVector v1(StandardMap); v1.PutValue(2);
  Epetra_IntVector v2(StandardMap); v2.PutValue(3);

  Epetra_Export identExporter(StandardMap,StandardMap); // Identity exporter
  EPETRA_TEST_ERR(!(v2.Export(v1, identExporter, Insert)==0),ierr);
  passed = (v2.MinValue()==2);
  EPETRA_TEST_ERR(!passed,ierr);

  v1.PutValue(1);
  Epetra_Import identImporter(StandardMap,StandardMap); // Identity importer
  EPETRA_TEST_ERR(!(v2.Import(v1, identExporter, Insert)==0),ierr);
  passed = passed && (v2.MaxValue()==1);
  EPETRA_TEST_ERR(!passed,ierr);

  if (verbose) {
    if (passed) cout << "Identity Import/Export Check OK" << endl << endl;
    else cout << "Identity Import/Export Check Failed" << endl << endl;
  }

  int NumSubMapElements = StandardMap.NumMyElements()/2;
  int SubStart = Comm.MyPID();
  NumSubMapElements = EPETRA_MIN(NumSubMapElements,StandardMap.NumMyElements()-SubStart);
  Epetra_Map SubMap(-1, NumSubMapElements, StandardMyGlobalElements+SubStart, 0, Comm);

  Epetra_IntVector v3(View, SubMap, SubMap.MyGlobalElements()); // Fill v3 with GID values for variety
  Epetra_Export subExporter(SubMap, StandardMap); // Export to a subset of indices of standard map
  EPETRA_TEST_ERR(!(v2.Export(v3,subExporter,Insert)==0),ierr);

  forierr = 0;
  for (i=0; i<SubMap.NumMyElements(); i++) {
    int i1 = StandardMap.LID(SubMap.GID(i));
    forierr += !(v3[i]==v2[i1]);
  }
  EPETRA_TEST_ERR(forierr,ierr);

  Epetra_Import subImporter(StandardMap, SubMap); // Import to a subset of indices of standard map
  EPETRA_TEST_ERR(!(v1.Import(v3,subImporter,Insert)==0),ierr);

  for (i=0; i<SubMap.NumMyElements(); i++) {
    int i1 = StandardMap.LID(SubMap.GID(i));
    forierr += !(v3[i]==v1[i1]);
  }
  EPETRA_TEST_ERR(forierr,ierr);

  if (verbose) {
    if (forierr==0) cout << "SubMap Import/Export Check OK" << endl << endl;
    else cout << "SubMap Import/Export Check Failed" << endl << endl;
  }

  // Release all objects

  delete [] SourceMyGlobalElements;
  delete [] TargetMyGlobalElements;
  delete [] OverlapMyGlobalElements;
  delete [] StandardMyGlobalElements;

  delete [] Values;
  delete [] Indices;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
