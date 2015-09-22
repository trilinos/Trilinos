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


// Epetra_BlockMap Test routine

#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_IntVector.h"
#include "Epetra_LongLongVector.h"
#include "Epetra_Import.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif

#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"
#include "Epetra_Version.h"

int main(int argc, char *argv[]) {

  int i, returnierr=0;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Uncomment to debug in parallel int tmp; if (Comm.MyPID()==0) cin >> tmp; Comm.Barrier();

  bool verbose = false;
  bool veryVerbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  // Check if we should print lots of results to standard out
  if (argc>2) if (argv[2][0]=='-' && argv[2][1]=='v') veryVerbose = true;

  if (verbose && Comm.MyPID()==0)
    std::cout << Epetra_Version() << std::endl << std::endl;

  if (!verbose) Comm.SetTracebackMode(0); // This should shut down any error traceback reporting

  if (verbose) std::cout << Comm << std::endl << std::flush;

  bool verbose1 = verbose;
  if (verbose) verbose = (Comm.MyPID()==0);

  bool veryVerbose1 = veryVerbose;
  if (veryVerbose) veryVerbose = (Comm.MyPID()==0);

  int NumMyElements = 100;
  if (veryVerbose1) NumMyElements = 10;
  NumMyElements += Comm.MyPID();
  int MaxNumMyElements = NumMyElements+Comm.NumProc()-1;
  int * ElementSizeList = new int[NumMyElements];
  long long * MyGlobalElements = new long long[NumMyElements];

  for (i = 0; i<NumMyElements; i++) {
    MyGlobalElements[i] = (long long)(Comm.MyPID()*MaxNumMyElements+i)*2;
    ElementSizeList[i] = i%6 + 2; // elementsizes go from 2 to 7
  }

  Epetra_BlockMap Map(
      -1LL, NumMyElements, MyGlobalElements, ElementSizeList,
      0, Comm
      );

  delete [] ElementSizeList;
  delete [] MyGlobalElements;

  Epetra_MapColoring C0(Map);

  int * elementColors = new int[NumMyElements];

  int maxcolor = 24;
  int * colorCount = new int[maxcolor];
  int ** colorLIDs = new int*[maxcolor];
  for (i=0; i<maxcolor; i++) colorCount[i] = 0;
  for (i=0; i<maxcolor; i++) colorLIDs[i] = 0;

  int defaultColor = C0.DefaultColor();
  for (i=0; i<Map.NumMyElements(); i++) {
    assert(C0[i]==defaultColor);
    assert(C0(Map.GID64(i))==defaultColor);
    if (i%2==0) C0[i] = i%6+5+i%14; // cycle through 5...23 on even elements
    else C0(Map.GID64(i)) = i%5+1; // cycle through 1...5 on odd elements
    elementColors[i] = C0[i]; // Record color of ith element for use below
    colorCount[C0[i]]++; // Count how many of each color for checking below
  }

  if (veryVerbose)
    std::cout << "Original Map Coloring using element-by-element definitions" << std::endl;
  if (veryVerbose1)
    std::cout <<  C0 << std::endl;

  int numColors = 0;
  for (i=0; i<maxcolor; i++)
    if (colorCount[i]>0) {
      numColors++;
      colorLIDs[i] = new int[colorCount[i]];
    }
  for (i=0; i<maxcolor; i++) colorCount[i] = 0;
  for (i=0; i<Map.NumMyElements(); i++) colorLIDs[C0[i]][colorCount[C0[i]]++] = i;



  int newDefaultColor = -1;
  Epetra_MapColoring C1(Map, elementColors, newDefaultColor);
  if (veryVerbose)
    std::cout << "Same Map Coloring using one-time construction" << std::endl;
  if (veryVerbose1)
    std::cout <<  C1 << std::endl;
  assert(C1.DefaultColor()==newDefaultColor);
  for (i=0; i<Map.NumMyElements(); i++) assert(C1[i]==C0[i]);

  Epetra_MapColoring C2(C1);
  if (veryVerbose)
    std::cout << "Same Map Coloring using copy constructor" << std::endl;
  if (veryVerbose1)
    std::cout <<  C1 << std::endl;
  for (i=0; i<Map.NumMyElements(); i++) assert(C2[i]==C0[i]);
  assert(C2.DefaultColor()==newDefaultColor);

  assert(numColors==C2.NumColors());

  for (i=0; i<maxcolor; i++) {
    int curNumElementsWithColor = C2.NumElementsWithColor(i);
    assert(colorCount[i]==curNumElementsWithColor);
    int * curColorLIDList = C2.ColorLIDList(i);
    if (curNumElementsWithColor==0) {
      assert(curColorLIDList==0);
    }
    else
      for (int j=0; j<curNumElementsWithColor; j++) assert(curColorLIDList[j]==colorLIDs[i][j]);
  }
  int curColor = 1;
  Epetra_Map * Map1 = C2.GenerateMap(curColor);
  Epetra_BlockMap * Map2 = C2.GenerateBlockMap(curColor);

  assert(Map1->NumMyElements()==colorCount[curColor]);
  assert(Map2->NumMyElements()==colorCount[curColor]);

  for (i=0; i<Map1->NumMyElements(); i++) {
    assert(Map1->GID64(i)==Map.GID64(colorLIDs[curColor][i]));
    assert(Map2->GID64(i)==Map.GID64(colorLIDs[curColor][i]));
    assert(Map2->ElementSize(i)==Map.ElementSize(colorLIDs[curColor][i]));
  }

  // Now test data redistribution capabilities


  Epetra_Map ContiguousMap(-1LL, Map.NumMyElements(), Map.IndexBase64(), Comm);
  // This vector contains the element sizes for the original map.
  Epetra_IntVector elementSizes(Copy, ContiguousMap, Map.ElementSizeList());
  Epetra_LongLongVector elementIDs(Copy, ContiguousMap, Map.MyGlobalElements64());
  Epetra_IntVector elementColorValues(Copy, ContiguousMap, C2.ElementColors());


  long long NumMyElements0 = 0;
  if (Comm.MyPID()==0) NumMyElements0 = Map.NumGlobalElements64();
  Epetra_Map CMap0(-1LL, NumMyElements0, Map.IndexBase64(), Comm);
  Epetra_Import importer(CMap0, ContiguousMap);
  Epetra_IntVector elementSizes0(CMap0);
  Epetra_LongLongVector elementIDs0(CMap0);
  Epetra_IntVector elementColorValues0(CMap0);
  elementSizes0.Import(elementSizes, importer, Insert);
  elementIDs0.Import(elementIDs, importer, Insert);
  elementColorValues0.Import(elementColorValues, importer, Insert);

  Epetra_BlockMap MapOnPE0(
      -1LL,NumMyElements0, elementIDs0.Values(),
      elementSizes0.Values(), Map.IndexBase64(), Comm
      );

  Epetra_Import importer1(MapOnPE0, Map);
  Epetra_MapColoring ColoringOnPE0(MapOnPE0);
  ColoringOnPE0.Import(C2, importer1, Insert);

  for (i=0; i<MapOnPE0.NumMyElements(); i++)
    assert(ColoringOnPE0[i]==elementColorValues0[i]);

  if (veryVerbose)
    std::cout << "Same Map Coloring on PE 0 only" << std::endl;
  if (veryVerbose1)
    std::cout <<  ColoringOnPE0 << std::endl;
  Epetra_MapColoring C3(Map);
  C3.Export(ColoringOnPE0, importer1, Insert);
  for (i=0; i<Map.NumMyElements(); i++) assert(C3[i]==C2[i]);
  if (veryVerbose)
    std::cout << "Same Map Coloring after Import/Export exercise" << std::endl;
  if (veryVerbose1)
    std::cout <<  ColoringOnPE0 << std::endl;


  if (verbose) std::cout << "Checked OK\n\n" << std::endl;

  if (verbose1) {
    if (verbose) std::cout << "Test ostream << operator" << std::endl << std::flush;
    std::cout << C0 << std::endl;
  }


  delete [] elementColors;
  for (i=0; i<maxcolor; i++) if (colorLIDs[i]!=0) delete [] colorLIDs[i];
  delete [] colorLIDs;
  delete [] colorCount;

  delete Map1;
  delete Map2;


#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return returnierr;
}

