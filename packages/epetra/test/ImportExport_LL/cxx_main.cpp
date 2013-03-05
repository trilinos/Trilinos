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


#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LongLongVector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_OffsetIndex.h"
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
#include "Epetra_Version.h"

int special_submap_import_test(Epetra_Comm& Comm);
int combine_mode_test(Epetra_Comm& Comm);
int alternate_import_constructor_test(Epetra_Comm& Comm);

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

  if (verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
              << " is alive."<<endl;

  // Redefine verbose to only print on PE 0
  if (verbose && Comm.MyPID()!=0) verbose = false;

  int NumMyEquations = 20;
  long long NumGlobalEquations = NumMyEquations*NumProc+EPETRA_MIN(NumProc,3);
  if (MyPID < 3) NumMyEquations++;
  // Construct a Source Map that puts approximately the same Number of equations on each processor in 
  // uniform global ordering

  Epetra_Map SourceMap(NumGlobalEquations, NumMyEquations, 0LL, Comm);
  
  // Get update list and number of local equations from newly created Map
  int NumMyElements = SourceMap.NumMyElements();
  long long * SourceMyGlobalElements = new long long[NumMyElements];
  SourceMap.MyGlobalElements(SourceMyGlobalElements);

  // Construct a Target Map that will contain:
  //  some unchanged elements (relative to the soure map),
  //  some permuted elements
  //  some off-processor elements
  Epetra_Vector RandVec(SourceMap);
  RandVec.Random(); // This creates a vector of random numbers between negative one and one.

  long long *TargetMyGlobalElements = new long long[NumMyElements];

  long long MinGID = SourceMap.MinMyGID64();
  for (i=0; i< NumMyEquations/2; i++) TargetMyGlobalElements[i] = i + MinGID; // Half will be the same...
  for (i=NumMyEquations/2; i<NumMyEquations; i++) {
    int index = abs((int)(((double) (NumGlobalEquations-1) ) * RandVec[i]));
    TargetMyGlobalElements[i] = EPETRA_MIN(NumGlobalEquations-1,(long long) EPETRA_MAX(0,index));
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

  Epetra_Map TargetMap((long long) -1, NumMyElements, TargetMyGlobalElements, 0LL, Comm);

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

  Epetra_Map StandardMap(NumGlobalEquations, NumMyEquations, 0LL, Comm);
  
  // Get update list and number of local equations from newly created Map
  NumMyElements = StandardMap.NumMyElements();
  long long * StandardMyGlobalElements = new long long[NumMyElements];
  StandardMap.MyGlobalElements(StandardMyGlobalElements);


  // Create a standard Epetra_CrsGraph

  Epetra_CrsGraph StandardGraph(Copy, StandardMap, 3);
  EPETRA_TEST_ERR(StandardGraph.IndicesAreGlobal(),ierr);
  EPETRA_TEST_ERR(StandardGraph.IndicesAreLocal(),ierr);
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  long long *Indices = new long long[2];
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
  EPETRA_TEST_ERR(!(StandardGraph.FillComplete()==0),ierr);
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
  EPETRA_TEST_ERR(!(StandardMatrix.FillComplete()==0),ierr);
  EPETRA_TEST_ERR(!(StandardMatrix.IndicesAreLocal()),ierr);
  //  EPETRA_TEST_ERR((StandardMatrix.StorageOptimized()),ierr);
  EPETRA_TEST_ERR((StandardMatrix.OptimizeStorage()),ierr);
  EPETRA_TEST_ERR(!(StandardMatrix.StorageOptimized()),ierr);
  EPETRA_TEST_ERR(StandardMatrix.UpperTriangular(),ierr);
  EPETRA_TEST_ERR(StandardMatrix.LowerTriangular(),ierr);

  // Construct an Overlapped Map of StandardMap that include the endpoints from two neighboring processors.

  int OverlapNumMyElements;
  long long OverlapMinMyGID;

  OverlapNumMyElements = NumMyElements + 1;
  if (MyPID==0) OverlapNumMyElements--;

  if (MyPID==0) OverlapMinMyGID = StandardMap.MinMyGID64();
  else OverlapMinMyGID = StandardMap.MinMyGID64()-1;

  long long * OverlapMyGlobalElements = new long long[OverlapNumMyElements];

  for (i=0; i< OverlapNumMyElements; i++) OverlapMyGlobalElements[i] = OverlapMinMyGID + i;

  Epetra_Map OverlapMap((long long) -1, OverlapNumMyElements, OverlapMyGlobalElements, 0LL, Comm);

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
      long long node_left = OverlapMyGlobalElements[i]-1;
      long long node_center = node_left + 1;
      long long node_right = node_left + 2;
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
    long long node_center = 0;
    EPETRA_TEST_ERR(!(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0),ierr);
  }
  if (MyPID==NumProc-1) {
    long long node_center = OverlapMyGlobalElements[OverlapNumMyElements-1];
    EPETRA_TEST_ERR(!(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0),ierr);
  }
    
  EPETRA_TEST_ERR(!(OverlapMatrix.FillComplete()==0),ierr);

  // Make a gathered matrix from OverlapMatrix.  It should be identical to StandardMatrix

  Epetra_CrsMatrix GatheredMatrix(Copy, StandardGraph);
  Epetra_Export Exporter(OverlapMap, StandardMap);
  EPETRA_TEST_ERR(!(GatheredMatrix.Export(OverlapMatrix, Exporter, Add)==0),ierr);
  EPETRA_TEST_ERR(!(GatheredMatrix.FillComplete()==0),ierr);

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

  //Do Again with use of Epetra_OffsetIndex object for speed
  Epetra_OffsetIndex OffsetIndex( OverlapMatrix.Graph(), GatheredMatrix.Graph(), Exporter );
  EPETRA_TEST_ERR(!(GatheredMatrix.Export(OverlapMatrix, Exporter, Add)==0),ierr);

  if (verbose) cout << "Optimized Matrix Export Check OK" << endl << endl;

  bool passed;
  Epetra_LongLongVector v1(StandardMap); v1.PutValue(2);
  Epetra_LongLongVector v2(StandardMap); v2.PutValue(3);

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
  Epetra_Map SubMap((long long) -1, NumSubMapElements, StandardMyGlobalElements+SubStart, 0LL, Comm);

  Epetra_LongLongVector v3(View, SubMap, SubMap.MyGlobalElements64()); // Fill v3 with GID values for variety
  Epetra_Export subExporter(SubMap, StandardMap); // Export to a subset of indices of standard map
  EPETRA_TEST_ERR(!(v2.Export(v3,subExporter,Insert)==0),ierr);

  forierr = 0;
  for (i=0; i<SubMap.NumMyElements(); i++) {
    int i1 = StandardMap.LID(SubMap.GID64(i));
    forierr += !(v3[i]==v2[i1]);
  }
  EPETRA_TEST_ERR(forierr,ierr);

  Epetra_Import subImporter(StandardMap, SubMap); // Import to a subset of indices of standard map
  EPETRA_TEST_ERR(!(v1.Import(v3,subImporter,Insert)==0),ierr);

  for (i=0; i<SubMap.NumMyElements(); i++) {
    int i1 = StandardMap.LID(SubMap.GID64(i));
    forierr += !(v3[i]==v1[i1]);
  }
  EPETRA_TEST_ERR(forierr,ierr);

  if (verbose) {
    if (forierr==0) cout << "SubMap Import/Export Check OK" << endl << endl;
    else cout << "SubMap Import/Export Check Failed" << endl << endl;
  }

#ifdef DOESNT_WORK_IN_PARALLEL
  forierr = special_submap_import_test(Comm);
  EPETRA_TEST_ERR(forierr, ierr);

  if (verbose) {
    if (forierr==0) cout << "Special SubMap Import Check OK" << endl << endl;
    else cout << "Special SubMap Import Check Failed" << endl << endl;
  }
#endif

  forierr =  alternate_import_constructor_test(Comm);
  EPETRA_TEST_ERR(forierr, ierr);

  if (verbose) {
    if (forierr==0) cout << "Alternative Import Constructor Check OK" << endl << endl;
    else cout << "Alternative Import Constructor Check Failed" << endl << endl;
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

int special_submap_import_test(Epetra_Comm& Comm)
{
  int localProc = Comm.MyPID();

  //set up ids_source and ids_target such that ids_source are only
  //a subset of ids_target, and furthermore that ids_target are ordered
  //such that the LIDs don't match up. In other words, even if gid 2 does
  //exist in both ids_source and ids_target, it will correspond to different
  //LIDs on at least 1 proc.
  //
  //This is to test a certain bug-fix in Epetra_Import where the 'RemoteLIDs'
  //array wasn't being calculated correctly on all procs.

  long long ids_source[1];
  ids_source[0] = localProc*2+2;

  long long ids_target[3];
  ids_target[0] = localProc*2+2;
  ids_target[1] = localProc*2+1;
  ids_target[2] = localProc*2+0;

  Epetra_Map map_source((long long) -1, 1, &ids_source[0], 0LL, Comm);
  Epetra_Map map_target((long long) -1, 3, &ids_target[0], 0LL, Comm);

  Epetra_Import importer(map_target, map_source);

  Epetra_LongLongVector vec_source(map_source);
  Epetra_LongLongVector vec_target(map_target);

  vec_target.PutValue(0);

  //set vec_source's contents so that entry[i] == GID[i].
  long long* GIDs = map_source.MyGlobalElements64();
  for(int i=0; i<map_source.NumMyElements(); ++i) {
    vec_source[i] = GIDs[i];
  }

  //Import vec_source into vec_target. This should result in the contents
  //of vec_target remaining 0 for the entries that don't exist in vec_source,
  //and other entries should be equal to the corresponding GID in the map.

  vec_target.Import(vec_source, importer, Insert);

  GIDs = map_target.MyGlobalElements64();
  int test_failed = 0;

  //the test passes if the i-th entry in vec_target equals either 0 or
  //GIDs[i].
  for(int i=0; i<vec_target.MyLength(); ++i) {
    if (vec_target[i] != GIDs[i] && vec_target[i] != 0) test_failed = 1;
  }

  int global_result;
  Comm.MaxAll(&test_failed, &global_result, 1);

  //If test didn't fail on any procs, global_result should be 0.
  //If test failed on any proc, global_result should be 1.
  return global_result;
}

int combine_mode_test(Epetra_Comm& Comm)
{
  int localProc = Comm.MyPID();


  long long ids_source[1];
  ids_source[0] = localProc*2+2;

  long long ids_target[3];
  ids_target[0] = localProc*2+2;
  ids_target[1] = localProc*2+1;
  ids_target[2] = localProc*2+0;

  Epetra_Map map_source((long long) -1, 1, &ids_source[0], 0LL, Comm);
  Epetra_Map map_target((long long) -1, 3, &ids_target[0], 0LL, Comm);

  Epetra_Import importer(map_target, map_source);

  Epetra_LongLongVector vec_source(map_source);
  Epetra_LongLongVector vec_target(map_target);

  vec_target.PutValue(0);

  //set vec_source's contents so that entry[i] == GID[i].
  long long* GIDs = map_source.MyGlobalElements64();
  for(int i=0; i<map_source.NumMyElements(); ++i) {
    vec_source[i] = GIDs[i];
  }

  //Import vec_source into vec_target. This should result in the contents
  //of vec_target remaining 0 for the entries that don't exist in vec_source,
  //and other entries should be equal to the corresponding GID in the map.

  vec_target.Import(vec_source, importer, Insert);

  GIDs = map_target.MyGlobalElements64();
  int test_failed = 0;

  //the test passes if the i-th entry in vec_target equals either 0 or
  //GIDs[i].
  for(int i=0; i<vec_target.MyLength(); ++i) {
    if (vec_target[i] != GIDs[i] && vec_target[i] != 0) test_failed = 1;
  }

  int global_result;
  Comm.MaxAll(&test_failed, &global_result, 1);

  //If test didn't fail on any procs, global_result should be 0.
  //If test failed on any proc, global_result should be 1.
  return global_result;
}


int test_import_gid(const char * name,Epetra_LongLongVector & Source, Epetra_LongLongVector & Target, const Epetra_Import & Import){
  int i;
  bool test_passed=true;

  // Setup
  for(i=0; i<Source.MyLength(); i++)
    Source[i] = Source.Map().GID64(i);
  Target.PutValue(0);

  // Import
  Target.Import(Source,Import,Add);

  // Test
  for(i=0; i<Target.MyLength(); i++){
    if(Target[i] != Target.Map().GID64(i)) test_passed=false;
  }

  if(!test_passed){
    printf("[%d] test_import_gid %s failed: ",Source.Map().Comm().MyPID(),name);
    for(i=0; i<Target.MyLength(); i++)
      printf("%2lld(%2lld) ",Target[i],Target.Map().GID64(i));
    printf("\n");
    fflush(stdout);
  }

  return !test_passed;
}



int alternate_import_constructor_test(Epetra_Comm& Comm) {
  int rv=0;
  int nodes_per_proc=10;
  int numprocs = Comm.NumProc();
  int mypid    = Comm.MyPID();

  // Only run if we have multiple procs & MPI
  if(numprocs==0) return 0;
#ifndef HAVE_MPI
  return 0;
#endif

  // Build Map 1 - linear
  Epetra_Map Map1((long long)-1,nodes_per_proc,(long long)0,Comm);
  
  // Build Map 2 - mod striped
  std::vector<long long> MyGIDs(nodes_per_proc);
  for(int i=0; i<nodes_per_proc; i++)
    MyGIDs[i] = (mypid*nodes_per_proc + i) % numprocs;
  Epetra_Map Map2((long long)-1,nodes_per_proc,&MyGIDs[0],(long long)0,Comm);

  // For testing
  Epetra_LongLongVector Source(Map1), Target(Map2);


  // Build Import 1 - normal
  Epetra_Import Import1(Map2,Map1);
  rv = rv|| test_import_gid("Alt test: 2 map constructor",Source,Target, Import1);

  // Build Import 2 - no-comm constructor
  int Nremote=Import1.NumRemoteIDs();
  const int * RemoteLIDs = Import1.RemoteLIDs();
  std::vector<int> RemotePIDs(Nremote+1); // I hate you, stl vector....
  std::vector<int> AllPIDs;  
  Epetra_Util::GetPids(Import1,AllPIDs,true);

  for(int i=0; i<Nremote; i++) {
    RemotePIDs[i]=AllPIDs[RemoteLIDs[i]];
  }
  Epetra_Import Import2(Import1.TargetMap(),Import1.SourceMap(),Nremote,&RemotePIDs[0],Import1.NumExportIDs(),Import1.ExportLIDs(),Import1.ExportPIDs());

  rv = rv || test_import_gid("Alt test: no comm constructor",Source,Target,Import2);


  // Build Import 3 - Remotes only
  Epetra_Import Import3(Import1.TargetMap(),Import1.SourceMap(),Nremote,&RemotePIDs[0]);
  rv = rv || test_import_gid("Alt test: remote only constructor",Source,Target, Import3);


  return rv;
}
