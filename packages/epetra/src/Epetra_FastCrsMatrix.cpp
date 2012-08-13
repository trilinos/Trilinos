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

#include "Epetra_FastCrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"
#include <limits.h>

//==============================================================================
Epetra_FastCrsMatrix::Epetra_FastCrsMatrix(const Epetra_CrsMatrix & Matrix, bool UseFloats) 
  : CrsMatrix_(Matrix),
    Values_(0),
    NumMyRows_(Matrix.NumMyRows()),
    NumMyNonzeros(Matrix.NumMyNonzeros()),
    ImportVector_(0),
    ExportVector_(0),
    CV_(Copy)
{
  if (!CrsMatrix_.Filled()) throw CrsMatrix_.ReportError("Input matrix must have called FillComplete()", -1);
  Allocate(UseFloats);
}

//==============================================================================
Epetra_FastCrsMatrix::~Epetra_FastCrsMatrix(){

  if (Values_!=0 && ValuesAllocated_) delete [] Values_;
  if (FloatValues_!=0) delete [] FloatValues_;
  if (Indices_!=0) delete [] Indices_;
  if (ShortIndices_!=0) delete [] ShortIndices_;

  if (ImportVector_!=0) delete ImportVector_;
  ImportVector_=0;
  if (ExportVector_!=0) delete ExportVector_;
  ExportVector_=0;
}

//==============================================================================
int Epetra_FastCrsMatrix::UpdateValues(const Epetra_CrsMatrix & Matrix) {
}

//==============================================================================
int Epetra_FastCrsMatrix::Allocate(bool UseFloats) {

  int i, j;
  
  // First determine how to handle values.  Three possibilities:
  // 1) Values are contiguously stored in user matrix, UseFloats is false so we point to the values in
  //    the user matrix.
  // 2) Values are not contiguously stored, UseFloats is false, so we copy values into a contiguous double array.
  // 3) UseFloats is true so we create a single precision copy of matrix values.

  if (CrsMatrix.IndicesAreContiguous() && !UseFloats) {
    ValuesAllocated_ = false;
    Values_ = CrsMatrix_[0]; // First value in the user's matrix
  }
  else if (!UseFloats) {
    // Allocate Values array
    Values_ = new double[NumMyNonzeros_];
    double * ptr = Values_;
    for (i=0; i<NumMyRows_; i++) {
      int NumEntries = CrsMatrix_.NumMyEntries(i);
      double * rowi = CrsMatrix_[i];
      for (j=0; j< NumEntries; j++) *ptr++ = *rowi++; // Copy values to contiguous storage
    }
  }
  else { // UseFloats == true
    FloatValues_ = new float[NumMyNonzeros_];
    float * ptr = FloatValues_;
    for (i=0; i<NumMyRows_; i++) {
      int NumEntries = CrsMatrix_.NumMyEntries(i);
      double * rowi = CrsMatrix_[i];
      for (j=0; j< NumEntries; j++) *ptr++ = (float) *rowi++; // convert and copy values
    }
  }

  // Next determine how to handle integers.  Two possibilities:
  // 1) Local column range is within the range of unsigned short ints, so we copy the indices to an array of unsigned shorts.
  // 2) Local column range is outside range of unsigned shorts and we copy to array of ints.
  // In both cases we embed the nonzero count per row into the index array.

  if (CrsMatrix_.NumMyCols()>USHRT_MAX) {
    // Allocate Values array
    Indices_ = new int[NumMyNonzeros_+NumMyRows_];
    int * ptr = Indices_;
    for (i=0; i<NumMyRows_; i++) {
      int NumEntries = CrsMatrix_.NumMyEntries(i);
      int * rowi = CrsMatrix_.Graph()[i];
      *ptr++ = NumEntries; // Pack this value first
      for (j=0; j< NumEntries; j++) *ptr++ = *rowi++; // Copy values to contiguous storage
    }
  }
  else { // CrsMatrix_.NumMyCols()<=USHRT_MAX
    ShortIndices_ = new unsigned short[NumMyNonzeros_+NumMyRows_];
    unsigned short * ptr = ShortIndices_;
    for (i=0; i<NumMyRows_; i++) {
      int NumEntries = CrsMatrix_.NumMyEntries(i);
      unsigned short * rowi = CrsMatrix_Graph()[i];
      *ptr++ = NumEntries; // Pack this value first
      for (j=0; j< NumEntries; j++) *ptr++ = (unsigned short) *rowi++; // convert and copy indices
    }
  }

  return(0);
}
//=============================================================================
int Epetra_FastCrsMatrix::Multiply(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const {
  //
  // This function forms the product y = A * x or y = A' * x
  //

  int i, j;
  double * xp = (double*)x.Values();
  double *yp = (double*)y.Values();
  int NumMyCols_ = NumMyCols();


  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()!=1) { delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),1); // Create import vector if needed
      ImportVector_->Import(x, *Importer(), Insert);
      xp = (double*)ImportVector_->Values();
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()!=1) { delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),1); // Create Export vector if needed
      yp = (double*)ExportVector_->Values();
    }

    // Do actual computation

    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *Values++;
      double sum = 0.0;
      for (j=0; j < NumEntries; j++) sum += RowValues[j] * xp[RowIndices[j]];

      yp[i] = sum;

    }
    if (Exporter()!=0) y.Export(*ExportVector_, *Exporter(), Add); // Fill y with Values from export vector
  }

  else { // Transpose operation


    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors

    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()!=1) { delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),1); // Create Export vector if needed
      ExportVector_->Import(x, *Exporter(), Insert);
      xp = (double*)ExportVector_->Values();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()!=1) { delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),1); // Create import vector if needed
      yp = (double*)ImportVector_->Values();
    }

    // Do actual computation

    for (i=0; i < NumMyCols_; i++) yp[i] = 0.0; // Initialize y for transpose multiply
        
    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *Values++;
      for (j=0; j < NumEntries; j++) yp[RowIndices[j]] += RowValues[j] * xp[i];
    }
    if (Importer()!=0) y.Export(*ImportVector_, *Importer(), Add); // Fill y with Values from export vector
  }

  UpdateFlops(2*NumGlobalNonzeros64());
  return(0);
}

//=============================================================================
int Epetra_FastCrsMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  //
  // This function forms the product Y = A * Y or Y = A' * X
  //
  if (X.NumVectors()==1 && Y.NumVectors()==1) {
    double * xp = (double *) X[0];
    double * yp = (double *) Y[0];
    Epetra_Vector x(View, X.Map(), xp);
    Epetra_Vector y(View, Y.Map(), yp);
    return(Multiply(TransA, x, y));
  }
  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.

  int i, j, k;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;

  double **Xp = (double**)X.Pointers();
  double **Yp = (double**)Y.Pointers();

  int NumVectors = X.NumVectors();
  int NumMyCols_ = NumMyCols();


  // Need to better manage the Import and Export vectors:
  // - Need accessor functions
  // - Need to make the NumVector match (use a View to do this)
  // - Need to look at RightScale and ColSum routines too.

  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()!=NumVectors) { delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),NumVectors); // Create import vector if needed
      ImportVector_->Import(X, *Importer(), Insert);
      Xp = (double**)ImportVector_->Pointers();
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()!=NumVectors) { delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),NumVectors); // Create Export vector if needed
      Yp = (double**)ExportVector_->Pointers();
    }

    // Do actual computation

    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *Values++;
      for (k=0; k<NumVectors; k++) {
	double sum = 0.0;
	for (j=0; j < NumEntries; j++) sum += RowValues[j] * Xp[k][RowIndices[j]];
	Yp[k][i] = sum;
      }
    }
    if (Exporter()!=0) Y.Export(*ExportVector_, *Exporter(), Add); // Fill Y with Values from export vector
  }
  else { // Transpose operation


    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors

    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()!=NumVectors) { delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),NumVectors); // Create Export vector if needed
      ExportVector_->Import(X, *Exporter(), Insert);
      Xp = (double**)ExportVector_->Pointers();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()!=NumVectors) { delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),NumVectors); // Create import vector if needed
      Yp = (double**)ImportVector_->Pointers();
    }

    // Do actual computation



    for (k=0; k<NumVectors; k++) 
      for (i=0; i < NumMyCols_; i++) Yp[k][i] = 0.0; // Initialize y for transpose multiply
    
    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *Values++;
      for (k=0; k<NumVectors; k++) {
	for (j=0; j < NumEntries; j++) Yp[k][RowIndices[j]] += RowValues[j] * Xp[k][i];
      }
    }
    if (Importer()!=0) Y.Export(*ImportVector_, *Importer(), Add); // Fill Y with Values from export vector
  }

  UpdateFlops(2*NumVectors*NumGlobalNonzeros64());
  return(0);
}

//=============================================================================
int Epetra_FastCrsMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_Vector& x, Epetra_Vector& y) const {
  //
  // This function find y such that Ly = x or Uy = x or the transpose cases.
  //

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.

  if ((Upper) && (!UpperTriangular())) EPETRA_CHK_ERR(-2);
  if ((!Upper) && (!LowerTriangular())) EPETRA_CHK_ERR(-3);
  if ((!UnitDiagonal) && (NoDiagonal())) EPETRA_CHK_ERR(-4); // If matrix has no diagonal, we must use UnitDiagonal
  if ((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_)) EPETRA_CHK_ERR(-5); // Need each row to have a diagonal
      

  int i, j, j0;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;
  int NumMyCols_ = NumMyCols();

  // If upper, point to last row
  if ((Upper && !Trans) || (!Upper && Trans)) {
    NumEntriesPerRow += NumMyRows_-1;
    Indices += NumMyRows_-1;
    Values += NumMyRows_-1;
  }
    
  double *xp = (double*)x.Values();
  double *yp = (double*)y.Values();

  if (!Trans) {

    if (Upper) {

      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=NumMyRows_-1; i >=0; i--) {
	int      NumEntries = *NumEntriesPerRow--;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Values--;
	double sum = 0.0;
	for (j=j0; j < NumEntries; j++) sum += RowValues[j] * yp[RowIndices[j]];
	
	if (UnitDiagonal) yp[i] = xp[i] - sum;
	else yp[i] = (xp[i] - sum)/RowValues[0];

      }
    }
    else {
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++ - j0;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Values++;
	double sum = 0.0;
	for (j=0; j < NumEntries; j++) sum += RowValues[j] * yp[RowIndices[j]];
	
	if (UnitDiagonal) yp[i] = xp[i] - sum;
	else yp[i] = (xp[i] - sum)/RowValues[NumEntries];

      }
    }
  }

  // ***********  Transpose case *******************************

  else {

    if (xp!=yp) for (i=0; i < NumMyCols_; i++) yp[i] = xp[i]; // Initialize y for transpose solve
    
    if (Upper) {

      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
    
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Values++;
	if (!UnitDiagonal) yp[i] = yp[i]/RowValues[0];
	for (j=j0; j < NumEntries; j++) yp[RowIndices[j]] -= RowValues[j] * yp[i];
      }
    }
    else {

      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
    
      for (i=NumMyRows_-1; i >= 0; i--) {
	int      NumEntries = *NumEntriesPerRow-- - j0;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Values--;
	if (!UnitDiagonal)  yp[i] = yp[i]/RowValues[NumEntries];
	for (j=0; j < NumEntries; j++) yp[RowIndices[j]] -= RowValues[j] * yp[i];
      }
    }

  }
  UpdateFlops(2*NumGlobalNonzeros64());
  return(0);
}

//=============================================================================
int Epetra_FastCrsMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  //
  // This function find Y such that LY = X or UY = X or the transpose cases.
  //
  if (X.NumVectors()==1 && Y.NumVectors()==1) {
    double * xp = (double *) X[0];
    double * yp = (double *) Y[0];
    Epetra_Vector x(View, X.Map(), xp);
    Epetra_Vector y(View, Y.Map(), yp);
    return(Solve(Upper, Trans, UnitDiagonal, x, y));
  }
  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.

  if ((Upper) && (!UpperTriangular())) EPETRA_CHK_ERR(-2);
  if ((!Upper) && (!LowerTriangular())) EPETRA_CHK_ERR(-3);
  if ((!UnitDiagonal) && (NoDiagonal())) EPETRA_CHK_ERR(-4); // If matrix has no diagonal, we must use UnitDiagonal
  if ((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_)) EPETRA_CHK_ERR(-5); // Need each row to have a diagonal

  int i, j, j0, k;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;
  double diag;

  // If upper, point to last row
  if ((Upper && !Trans) || (!Upper && Trans)) {
    NumEntriesPerRow += NumMyRows_-1;
    Indices += NumMyRows_-1;
    Values += NumMyRows_-1;
  }

  double **Xp = (double**)X.Pointers();
  double **Yp = (double**)Y.Pointers();

  int NumVectors = X.NumVectors();

  if (!Trans) {
    

    if (Upper) {
      
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=NumMyRows_-1; i >=0; i--) {
	int      NumEntries = *NumEntriesPerRow--;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Values--;
	if (!UnitDiagonal) diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
	for (k=0; k<NumVectors; k++) {
	  double sum = 0.0;
	  for (j=j0; j < NumEntries; j++) sum += RowValues[j] * Yp[k][RowIndices[j]];
	  
	  if (UnitDiagonal) Yp[k][i] = Xp[k][i] - sum;
	  else Yp[k][i] = (Xp[k][i] - sum)*diag;
	}
      }
    }
    else {
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++ - j0;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Values++;
	if (!UnitDiagonal) diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
	for (k=0; k<NumVectors; k++) {
	  double sum = 0.0;
	  for (j=0; j < NumEntries; j++) sum += RowValues[j] * Yp[k][RowIndices[j]];
	  
	  if (UnitDiagonal) Yp[k][i] = Xp[k][i] - sum;
	  else Yp[k][i] = (Xp[k][i] - sum)*diag;
	}
      }
    }
  }
  // ***********  Transpose case *******************************

  else {

    for (k=0; k<NumVectors; k++) 
      if (Yp[k]!=Xp[k]) for (i=0; i < NumMyRows_; i++) Yp[k][i] = Xp[k][i]; // Initialize y for transpose multiply
    
    if (Upper) {
      
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Values++;
	if (!UnitDiagonal) diag = 1.0/RowValues[j0]; // Take inverse of diagonal once for later use
	for (k=0; k<NumVectors; k++) {
	  if (!UnitDiagonal) Yp[k][i] = Yp[k][i]*diag;
	  for (j=j0; j < NumEntries; j++) Yp[k][RowIndices[j]] -= RowValues[j] * Yp[k][i];
	}
      }
    }
    else {
      
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      
      for (i=NumMyRows_-1; i>=0; i--) {
	int      NumEntries = *NumEntriesPerRow-- - j0;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Values--;
	for (k=0; k<NumVectors; k++) {
	  if (!UnitDiagonal)  Yp[k][i] = Yp[k][i]/Xp[k][i];
	  for (j=0; j < NumEntries; j++) Yp[k][RowIndices[j]] -= RowValues[j] * Yp[k][i];
        }
      }
    }
  }
  
  UpdateFlops(2*NumVectors*NumGlobalNonzeros64());
  return(0);
}
