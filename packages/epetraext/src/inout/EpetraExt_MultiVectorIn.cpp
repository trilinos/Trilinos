//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER
#include "EpetraExt_MultiVectorIn.h"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_BlockMap.h"

using namespace EpetraExt;
namespace EpetraExt {

int MatrixMarketFileToMultiVector( const char *filename, const Epetra_BlockMap & map, Epetra_MultiVector * & A) {

  const int lineLength = 1025;
  const int tokenLength = 35;
  char line[lineLength];
  char token1[tokenLength];
  char token2[tokenLength];
  char token3[tokenLength];
  char token4[tokenLength];
  char token5[tokenLength];
  int M, N;

  FILE * handle = 0;

  handle = fopen(filename,"r");  // Open file
  if (handle == 0)
    EPETRA_CHK_ERR(-1); // file not found

  // Check first line, which should be "%%MatrixMarket matrix coordinate real general" (without quotes)
  if(fgets(line, lineLength, handle)==0) return(-1);
  if(sscanf(line, "%s %s %s %s %s", token1, token2, token3, token4, token5 )==0) return(-1);
  if (strcmp(token1, "%%MatrixMarket") ||
      strcmp(token2, "matrix") ||
      strcmp(token3, "array") ||
      strcmp(token4, "real") ||
      strcmp(token5, "general")) return(-1);

  // Next, strip off header lines (which start with "%")
  do {
    if(fgets(line, lineLength, handle)==0) return(-1);
  } while (line[0] == '%');

  // Next get problem dimensions: M, N
  if(sscanf(line, "%d %d", &M, &N)==0) return(-1);

  // Compute the offset for each processor for when it should start storing values
  int numMyPoints = map.NumMyPoints();
  int offset;
  map.Comm().ScanSum(&numMyPoints, &offset, 1); // ScanSum will compute offsets for us
  offset -= numMyPoints; // readjust for my PE

  // Now construct vector/multivector
  if (N==1)
    A = new Epetra_Vector(map);
  else
    A = new Epetra_MultiVector(map, N);

  double ** Ap = A->Pointers();

  for (int j=0; j<N; j++) {
    double * v = Ap[j];

    // Now read in lines that we will discard
    for (int i=0; i<offset; i++)
      if(fgets(line, lineLength, handle)==0) return(-1);
    
    // Now read in each value and store to the local portion of the the  if the row is owned.
    double V;
    for (int i=0; i<numMyPoints; i++) {
      if(fgets(line, lineLength, handle)==0) return(-1);
      if(sscanf(line, "%lg\n", &V)==0) return(-1);
      v[i] = V;
    }
    // Now read in the rest of the lines to discard
    for (int i=0; i < M-numMyPoints-offset; i++) {
      if(fgets(line, lineLength, handle)==0) return(-1);
    }
  }

  if (fclose(handle)) return(-1);
  
  return(0);
}

} // namespace EpetraExt
