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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h" 
#include "Epetra_Import.h"
#include "Epetra_Export.h"
int main(int argc, char *argv[])
{
  using namespace std;
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

  cout << Comm <<endl;

  if (Comm.NumProc()!=3) cout << "Must be run on 3 MPI processes." << endl;

  char tmp;
  if (Comm.MyPID()==0) cout << "Press any key to continue..."<< endl;
  if (Comm.MyPID()==0) cin >> tmp;
  Comm.Barrier();
  int NumMyElements = 2;
  Epetra_Map SourceMap((long long) -1, NumMyElements, 0, Comm);
  long long NumGlobalElements = SourceMap.NumGlobalElements();
  
  Epetra_Vector SourceX(SourceMap);

  // Define SourceX = [0, 1, 2, 3, 4, 5]' with 2 elements per process

  SourceX[0] = (double) SourceMap.GID(0);
  SourceX[1] = (double) SourceMap.GID(1);

  cout << SourceX; // Print it to stdout

  // Construct first target map: PE0 gets all elements of SourceX

  long long TargetOneNumMyElements;
  if (Comm.MyPID()==0) TargetOneNumMyElements = NumGlobalElements;
  else TargetOneNumMyElements = 0;

  Epetra_Map TargetOneMap((long long) -1, TargetOneNumMyElements, 0, Comm);

  Epetra_Import ImporterOne(TargetOneMap, SourceMap);

  Epetra_Vector TargetOneX(TargetOneMap);

  TargetOneX.Import(SourceX, ImporterOne, Insert);

  cout << TargetOneX << endl;

  long long * GlobalElementList = new long long[2];
  GlobalElementList[0] = 5 - SourceMap.GID(0);
  GlobalElementList[1] = 5 - SourceMap.GID(1);
  
  Epetra_Map TargetTwoMap((long long) -1, 2, GlobalElementList, 0, Comm);

  Epetra_Import ImporterTwo(TargetTwoMap, TargetOneMap);

  Epetra_Vector TargetTwoX(TargetTwoMap);

  TargetTwoX.Import(TargetOneX, ImporterTwo, Insert);

  cout << TargetTwoX << endl;
  
  Epetra_LocalMap TargetThreeMap((long long) 6, 0, Comm);

  Epetra_Import ImporterThree(TargetThreeMap, TargetOneMap);

  Epetra_Vector TargetThreeX(TargetThreeMap);

  TargetThreeX.Import(TargetOneX, ImporterThree, Insert);

  cout << TargetThreeX << endl;
  
  Epetra_Export ExporterOne(TargetThreeMap, SourceMap);

  Epetra_Vector TargetFourX(SourceMap);

  TargetFourX.Import(TargetThreeX, ExporterOne, Add);

  cout << TargetFourX << endl;
  


  MPI_Finalize();
  return 0;
}
