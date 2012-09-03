//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#include "Poisson2dOperator.h"
#include <stdlib.h>
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"

//==============================================================================
Poisson2dOperator::Poisson2dOperator(int nx, int ny, const Epetra_Comm & comm) 
  : nx_(nx),
    ny_(ny),
    useTranspose_(false),
    comm_(comm),
    map_(0),
    numImports_(0),
    importIDs_(0),
    importMap_(0),
    importer_(0),
    importX_(0),
    Label_("2D Poisson Operator") {

  int numProc = comm.NumProc(); // Get number of processors
  int myPID = comm.MyPID(); // My rank
  if (2*numProc > ny) { // ny must be >= 2*numProc (to avoid degenerate cases)
    ny = 2*numProc; ny_ = ny;
    cout << " Increasing ny to " << ny << " to avoid degenerate distribution on " << numProc << " processors." << endl;
  }
  
  int chunkSize = ny/numProc;
  int remainder = ny%numProc;

  if (myPID+1 <= remainder) chunkSize++; // add on remainder

  myny_ = chunkSize;

  map_ = new Epetra_Map((long long)-1, nx*chunkSize, 0, comm_);

  if (numProc>1) {
    // Build import GID list to build import map and importer
    if (myPID>0) numImports_ += nx;
    if (myPID+1<numProc) numImports_ += nx;
    
    if (numImports_>0) importIDs_ = new long long[numImports_];
    long long * ptr = importIDs_;
    long long minGID = map_->MinMyGID64();
    long long maxGID = map_->MaxMyGID64();
    
    if (myPID>0) for (int i=0; i< nx; i++) *ptr++ = minGID - nx + i;
    if (myPID+1<numProc) for (int i=0; i< nx; i++) *ptr++ = maxGID + i +1;
    
    // At the end of the above step importIDs_ will have a list of global IDs that are needed
    // to compute the matrix multiplication operation on this processor.  Now build import map 
    // and importer
    
    
    importMap_ = new Epetra_Map((long long)-1, numImports_, importIDs_, 0, comm_);
    
    importer_ = new Epetra_Import(*importMap_, *map_);

  }
}
//==============================================================================
Poisson2dOperator::~Poisson2dOperator() {
  if (importX_!=0) delete importX_;
  if (importer_!=0) delete importer_;
  if (importMap_!=0) delete importMap_;
  if (importIDs_!=0) delete [] importIDs_;
  if (map_!=0) delete map_;
}
//==============================================================================
int Poisson2dOperator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {


  // This is a very brain-dead implementation of a 5-point finite difference stencil, but
  // it should illustrate the basic process for implementing the Epetra_Operator interface.

  if (!X.Map().SameAs(OperatorDomainMap())) abort();  // These aborts should be handled as int return codes.
  if (!Y.Map().SameAs(OperatorRangeMap())) abort();
  if (Y.NumVectors()!=X.NumVectors()) abort();

  if (comm_.NumProc()>1) {
    if (importX_==0)
      importX_ = new Epetra_MultiVector(*importMap_, X.NumVectors());
    else if (importX_->NumVectors()!=X.NumVectors()) {
      delete importX_;
      importX_ = new Epetra_MultiVector(*importMap_, X.NumVectors());
    }
    importX_->Import(X, *importer_, Insert); // Get x values we need
  }

  double * importx1 = 0;
  double * importx2 = 0;
  int nx = nx_;
  //int ny = ny_;

  for (int j=0; j < X.NumVectors(); j++) {

    const double * x = X[j];
    if (comm_.NumProc()>1) {
      importx1 = (*importX_)[j];
      importx2 = importx1+nx;
      if (comm_.MyPID()==0) importx2=importx1;
    }
    double * y = Y[j];
    if (comm_.MyPID()==0) {
      y[0] = 4.0*x[0]-x[nx]-x[1];
      y[nx-1] = 4.0*x[nx-1]-x[nx-2]-x[nx+nx-1];
      for (int ix=1; ix< nx-1; ix++)
	y[ix] = 4.0*x[ix]-x[ix-1]-x[ix+1]-x[ix+nx];
    }
    else {
      y[0] = 4.0*x[0]-x[nx]-x[1]-importx1[0];
      y[nx-1] = 4.0*x[nx-1]-x[nx-2]-x[nx+nx-1]-importx1[nx-1];
      for (int ix=1; ix< nx-1; ix++)
	y[ix] = 4.0*x[ix]-x[ix-1]-x[ix+1]-x[ix+nx]-importx1[ix];
    }
    if (comm_.MyPID()+1==comm_.NumProc()) {
      int curxy = nx*myny_-1;
      y[curxy] = 4.0*x[curxy]-x[curxy-nx]-x[curxy-1];
      curxy -= (nx-1);
      y[curxy] = 4.0*x[curxy]-x[curxy-nx]-x[curxy+1];
      for (int ix=1; ix< nx-1; ix++) {
	curxy++;
	y[curxy] = 4.0*x[curxy]-x[curxy-1]-x[curxy+1]-x[curxy-nx];
      }
    }
    else {
      int curxy = nx*myny_-1;
      y[curxy] = 4.0*x[curxy]-x[curxy-nx]-x[curxy-1]-importx2[nx-1];
      curxy -= (nx-1);
      y[curxy] = 4.0*x[curxy]-x[curxy-nx]-x[curxy+1]-importx2[0];
      for (int ix=1; ix< nx-1; ix++) {
	curxy++;
	y[curxy] = 4.0*x[curxy]-x[curxy-1]-x[curxy+1]-x[curxy-nx]-importx2[ix];
      }
    }
    for (int iy=1; iy< myny_-1; iy++) {
      int curxy = nx*(iy+1)-1;
      y[curxy] = 4.0*x[curxy]-x[curxy-nx]-x[curxy-1]-x[curxy+nx];
      curxy -= (nx-1);
      y[curxy] = 4.0*x[curxy]-x[curxy-nx]-x[curxy+1]-x[curxy+nx];
      for (int ix=1; ix< nx-1; ix++) {
	curxy++;
	y[curxy] = 4.0*x[curxy]-x[curxy-1]-x[curxy+1]-x[curxy-nx]-x[curxy+nx];
      }
    }
  }
  return(0);
}
//==============================================================================
Epetra_CrsMatrix * Poisson2dOperator::GeneratePrecMatrix() const {

  // Generate a tridiagonal matrix as an Epetra_CrsMatrix
  // This method illustrates how to generate a matrix that is an approximation to the true
  // operator.  Given this matrix, we can use any of the Aztec or IFPACK preconditioners.


  // Create a Epetra_Matrix
  Epetra_CrsMatrix * A = new Epetra_CrsMatrix(Copy, *map_, 3);

  int NumMyElements = map_->NumMyElements();
  long long NumGlobalElements = map_->NumGlobalElements64();

  // Add  rows one-at-a-time
  double negOne = -1.0;
  double posTwo = 4.0;
  for (int i=0; i<NumMyElements; i++) {
    long long GlobalRow = A->GRID64(i); long long RowLess1 = GlobalRow - 1; long long RowPlus1 = GlobalRow + 1;

    if (RowLess1!=-1) A->InsertGlobalValues(GlobalRow, 1, &negOne, &RowLess1);
    if (RowPlus1!=NumGlobalElements) A->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus1);
    A->InsertGlobalValues(GlobalRow, 1, &posTwo, &GlobalRow);
  }

  // Finish up
  A->FillComplete();

  return(A);
}
