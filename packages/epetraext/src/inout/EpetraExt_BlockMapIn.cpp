//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
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
// ***********************************************************************
//@HEADER
#include "EpetraExt_BlockMapIn.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"

using namespace EpetraExt;
namespace EpetraExt {

int MatrixMarketFileToMap( const char *filename, const Epetra_Comm & comm, Epetra_Map * & map) {

  Epetra_BlockMap * bmap;
  if (MatrixMarketFileToBlockMap(filename, comm, bmap)) return(-1);
  map = dynamic_cast<Epetra_Map *>(bmap);
  return(0);
}

int MatrixMarketFileToBlockMap( const char *filename, const Epetra_Comm & comm, Epetra_BlockMap * & map) {

  const int lineLength = 1025;
  char line[lineLength];
  char token[lineLength];
  int M, N, numProc, MaxElementSize, MinElementSize, NumMyElements, IndexBase, NumGlobalElements, firstGid;

  FILE * handle = 0;

  bool inHeader = true;

  handle = fopen(filename,"r");
  while (inHeader) {
    if(fgets(line, lineLength, handle)==0) return(-1);
    if(sscanf(line, "%s", token)==0) return(-1);
    if (!strcmp(token, "%NumProc:")) inHeader = false;
  }

  if(fgets(line, lineLength, handle)==0) return(-1); // numProc value
  if(sscanf(line, "%s %d", token, &numProc)==0) return(-1);

  if(fgets(line, lineLength, handle)==0) return(-1); // MaxElementSize header line
  if(fgets(line, lineLength, handle)==0) return(-1); // MaxElementSize value
  if(sscanf(line, "%s %d", token, &MaxElementSize)==0) return(-1);

  if(fgets(line, lineLength, handle)==0) return(-1); // MinElementSize header line
  if(fgets(line, lineLength, handle)==0) return(-1); // MinElementSize value
  if(sscanf(line, "%s %d", token, &MinElementSize)==0) return(-1);

  if(fgets(line, lineLength, handle)==0) return(-1); // IndexBase header line
  if(fgets(line, lineLength, handle)==0) return(-1); // IndexBase value
  if(sscanf(line, "%s %d", token, &IndexBase)==0) return(-1);

  if(fgets(line, lineLength, handle)==0) return(-1); // NumGlobalElements header line
  if(fgets(line, lineLength, handle)==0) return(-1); // NumGlobalElements value
  if(sscanf(line, "%s %d", token, &NumGlobalElements)==0) return(-1);

  int ierr = 0;
  if (comm.NumProc()==numProc) {
    if(fgets(line, lineLength, handle)==0) return(-1); // NumMyElements header line
    firstGid = 0;
    for (int i=0; i<comm.MyPID(); i++) {
      if(fgets(line, lineLength, handle)==0) return(-1); // ith NumMyElements value
      if(sscanf(line, "%s %d", token, &NumMyElements)==0) return(-1);
      firstGid += NumMyElements;
    }
 
    if(fgets(line, lineLength, handle)==0) return(-1); // This PE's NumMyElements value
    if(sscanf(line, "%s %d", token, &NumMyElements)==0) return(-1);

    for (int i=comm.MyPID()+1; i<numProc; i++) {
      if(fgets(line, lineLength, handle)==0) return(-1); // ith NumMyElements value (dump these)
    }
  }
  else {
    ierr = 1; // Warning error, different number of processors.

    if(fgets(line, lineLength, handle)==0) return(-1); // NumMyElements header line
    for (int i=0; i<numProc; i++) {
      if(fgets(line, lineLength, handle)==0) return(-1); // ith NumMyElements value (dump these)
    }

    NumMyElements = NumGlobalElements/comm.NumProc();
    firstGid = comm.MyPID()*NumMyElements;
    int remainder = NumGlobalElements%comm.NumProc();
    if (comm.MyPID()<remainder) NumMyElements++;
    int extra = remainder;
    if (comm.MyPID()<remainder) extra = comm.MyPID();
    firstGid += extra;
  }
  if(fgets(line, lineLength, handle)==0) return(-1); // Number of rows, columns
  if(sscanf(line, "%d %d", &M, &N)==0) return(-1);

  bool doSizes = (N>1);
  Epetra_IntSerialDenseVector v1(NumMyElements);
  Epetra_IntSerialDenseVector v2(NumMyElements);
  for (int i=0; i<firstGid; i++) {
    if(fgets(line, lineLength, handle)==0) return(-1); // dump these
  }

  if (doSizes) {
    for (int i=0; i<NumMyElements; i++) {
      if(fgets(line, lineLength, handle)==0) return(-1);
      if(sscanf(line, "%d %d", &v1[i], &v2[i])==0) return(-1); // load v1, v2
    }
  }
  else {
    for (int i=0; i<NumMyElements; i++) {
      if(fgets(line, lineLength, handle)==0) return(-1);
      if(sscanf(line, "%d", &v1[i])==0) return(-1); // load v1
      v2[i] = MinElementSize; // Fill with constant size
    }
  }
  if (fclose(handle)) return(-1);

  comm.Barrier();

  if (MinElementSize==1 && MaxElementSize==1)
    map = new Epetra_Map(-1, NumMyElements, v1.Values(), IndexBase, comm);
  else
    map = new Epetra_BlockMap(-1, NumMyElements, v1.Values(), v2.Values(), IndexBase, comm);
  return(0);
}
} // namespace EpetraExt
