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

// Epetra_BlockMap Test routine

#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_IntVector.h"
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
    cout << Epetra_Version() << endl << endl;

  if (!verbose) Comm.SetTracebackMode(0); // This should shut down any error traceback reporting

  if (verbose) cout << Comm << endl << flush;

  bool verbose1 = verbose;
  if (verbose) verbose = (Comm.MyPID()==0);

  bool veryVerbose1 = veryVerbose;
  if (veryVerbose) veryVerbose = (Comm.MyPID()==0);

  int NumMyElements = 100;
  if (veryVerbose1) NumMyElements = 10;
  NumMyElements += Comm.MyPID();
  int MaxNumMyElements = NumMyElements+Comm.NumProc()-1;
  int * ElementSizeList = new int[NumMyElements];
  int * MyGlobalElements = new int[NumMyElements];

  for (i = 0; i<NumMyElements; i++) {
    MyGlobalElements[i] = (Comm.MyPID()*MaxNumMyElements+i)*2;
    ElementSizeList[i] = i%6 + 2; // elementsizes go from 2 to 7
  }

  Epetra_BlockMap Map(-1, NumMyElements, MyGlobalElements, ElementSizeList,
		      0, Comm);

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
    assert(C0(Map.GID(i))==defaultColor);
    if (i%2==0) C0[i] = i%6+5+i%14; // cycle through 5...23 on even elements
    else C0(Map.GID(i)) = i%5+1; // cycle through 1...5 on odd elements
    elementColors[i] = C0[i]; // Record color of ith element for use below
    colorCount[C0[i]]++; // Count how many of each color for checking below
  }
  
  if (veryVerbose)
    cout << "Original Map Coloring using element-by-element definitions" << endl;
  if (veryVerbose1)
    cout <<  C0 << endl;

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
    cout << "Same Map Coloring using one-time construction" << endl;
  if (veryVerbose1)
    cout <<  C1 << endl;
  assert(C1.DefaultColor()==newDefaultColor);
  for (i=0; i<Map.NumMyElements(); i++) assert(C1[i]==C0[i]);

  Epetra_MapColoring C2(C1);
  if (veryVerbose)
    cout << "Same Map Coloring using copy constructor" << endl;
  if (veryVerbose1)
    cout <<  C1 << endl;
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
    assert(Map1->GID(i)==Map.GID(colorLIDs[curColor][i]));
    assert(Map2->GID(i)==Map.GID(colorLIDs[curColor][i]));
    assert(Map2->ElementSize(i)==Map.ElementSize(colorLIDs[curColor][i]));
  }

  // Now test data redistribution capabilities


  Epetra_Map ContiguousMap(-1, Map.NumMyElements(), Map.IndexBase(), Comm);
  // This vector contains the element sizes for the original map.
  Epetra_IntVector elementSizes(Copy, ContiguousMap, Map.ElementSizeList());
  Epetra_IntVector elementIDs(Copy, ContiguousMap, Map.MyGlobalElements());
  Epetra_IntVector elementColorValues(Copy, ContiguousMap, C2.ElementColors());


  int NumMyElements0 = 0;
  if (Comm.MyPID()==0) NumMyElements0 = Map.NumGlobalElements();
  Epetra_Map CMap0(-1, NumMyElements0, Map.IndexBase(), Comm);
  Epetra_Import importer(CMap0, ContiguousMap);
  Epetra_IntVector elementSizes0(CMap0);
  Epetra_IntVector elementIDs0(CMap0);
  Epetra_IntVector elementColorValues0(CMap0);
  elementSizes0.Import(elementSizes, importer, Insert);
  elementIDs0.Import(elementIDs, importer, Insert);
  elementColorValues0.Import(elementColorValues, importer, Insert);

  Epetra_BlockMap MapOnPE0(-1,NumMyElements0, elementIDs0.Values(), 
			   elementSizes0.Values(), Map.IndexBase(), Comm);

  Epetra_Import importer1(MapOnPE0, Map);
  Epetra_MapColoring ColoringOnPE0(MapOnPE0);
  ColoringOnPE0.Import(C2, importer1, Insert);

  for (i=0; i<MapOnPE0.NumMyElements(); i++)
    assert(ColoringOnPE0[i]==elementColorValues0[i]);

  if (veryVerbose)
    cout << "Same Map Coloring on PE 0 only" << endl;
  if (veryVerbose1)
    cout <<  ColoringOnPE0 << endl;
  Epetra_MapColoring C3(Map);
  C3.Export(ColoringOnPE0, importer1, Insert);
  for (i=0; i<Map.NumMyElements(); i++) assert(C3[i]==C2[i]);
  if (veryVerbose)
    cout << "Same Map Coloring after Import/Export exercise" << endl;
  if (veryVerbose1)
    cout <<  ColoringOnPE0 << endl;
   
  
  if (verbose) cout << "Checked OK\n\n" <<endl;

  if (verbose1) {
    if (verbose) cout << "Test ostream << operator" << endl << flush;
    cout << C0 << endl;
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

