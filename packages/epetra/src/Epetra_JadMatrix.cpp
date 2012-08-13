
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

#include "Epetra_JadMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Util.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

//==============================================================================
Epetra_JadMatrix::Epetra_JadMatrix(const Epetra_RowMatrix & Matrix) 
  : Epetra_BasicRowMatrix(Matrix.RowMatrixRowMap().Comm()),
    Values_(0),
    Indices_(0),
    IndexOffset_(0),
    Profile_(0),
    RowPerm_(0),
    InvRowPerm_(0),
    NumJaggedDiagonals_(Matrix.MaxNumEntries())
{
  SetMaps(Matrix.RowMatrixRowMap(), Matrix.RowMatrixColMap(), Matrix.OperatorDomainMap(), Matrix.OperatorRangeMap());
  if (!Matrix.Filled()) throw Matrix.RowMatrixRowMap().ReportError("Input matrix must have called FillComplete()", -1);
  Allocate(Matrix);
  SetLabel("Epetra::JadMatrix");
}

//==============================================================================
Epetra_JadMatrix::~Epetra_JadMatrix(){}

//==============================================================================
int Epetra_JadMatrix::UpdateValues(const Epetra_RowMatrix & Matrix, bool CheckStructure) {

  int NumEntries;
  int * Indices = 0;
  double * Values =0;

  int ierr = 0;

  try { // If matrix is an Epetra_CrsMatrix, we can get date much more cheaply

    const Epetra_CrsMatrix & A = dynamic_cast<const Epetra_CrsMatrix &>(Matrix);

    for (int i1=0; i1<NumMyRows_; i1++) {
      
      EPETRA_CHK_ERR(A.ExtractMyRowView(i1, NumEntries, Values, Indices)); // Get the current row based on the permutation
      int i = InvRowPerm_[i1]; // Determine permuted row location
      for (int j=0; j< NumEntries; j++) Values_[IndexOffset_[j]+i] = Values[j];
      if (CheckStructure)
	for (int j=0; j< NumEntries; j++) if (Indices_[IndexOffset_[j]+i] != Indices[j]) ierr = - 1;
    }
  }
  catch (...) { // Otherwise just live with RowMatrix interface
    
    Epetra_SerialDenseVector curValues(NumJaggedDiagonals_);
    Epetra_IntSerialDenseVector curIndices(NumJaggedDiagonals_);
    Indices = curIndices.Values();
    Values = curValues.Values();
    for (int i1=0; i1<NumMyRows_; i1++) {
      EPETRA_CHK_ERR(Matrix.ExtractMyRowCopy(i1, NumJaggedDiagonals_, NumEntries, Values, Indices)); // Get current row based on the permutation
      int i = InvRowPerm_[i1]; // Determine permuted row location
      for (int j=0; j< NumEntries; j++) Values_[IndexOffset_[j]+i] = Values[j];
      if (CheckStructure)
	for (int j=0; j< NumEntries; j++) if (Indices_[IndexOffset_[j]+i] != Indices[j]) ierr = - 1;
    }
  }

  HaveNumericConstants_ = false;
  EPETRA_CHK_ERR(ierr);
  return(ierr);
}

//==============================================================================
void Epetra_JadMatrix::Allocate(const Epetra_RowMatrix & Matrix) {

  // Allocate IndexOffset storage
  int numMyRows = Matrix.NumMyRows();
  int numMyNonzeros = Matrix.NumMyNonzeros();

  IndexOffset_.Resize(NumJaggedDiagonals_+1);

  // Next compute permutation of rows
  RowPerm_.Resize(numMyRows);
  InvRowPerm_.Resize(numMyRows);
  Profile_.Resize(numMyRows);
  for (int i=0; i<numMyRows; i++) {
    int NumEntries;
    Matrix.NumMyRowEntries(i, NumEntries);
    Profile_[i] = NumEntries;
    RowPerm_[i] = i;
  }

  Epetra_Util sorter;
  int * RowPerm = RowPerm_.Values();
  sorter.Sort(false, numMyRows, Profile_.Values(), 0, 0, 1, &RowPerm, 0, 0);
  //cout << "Profile = " << Profile_ << endl;
  //cout << "RowPerm = " << RowPerm_ << endl;
  for (int i=0; i<numMyRows; i++) InvRowPerm_[RowPerm[i]] = i; // Compute inverse row permutation
  //cout << "InvRowPerm = " << InvRowPerm_ << endl;

  // Now build IndexOffsets:  These contain the lengths of the jagged diagonals

  for (int i=0; i<NumJaggedDiagonals_; i++) IndexOffset_[i] = 0;

  int curOffset = numMyRows;
  int * curIndex = IndexOffset_.Values(); // point to first index (will be incremented immediately below)
  for (int i=1; i<NumJaggedDiagonals_+1; i++) {
    curIndex++;
    while (*curIndex==0) {
      if (Profile_[curOffset-1]<i) curOffset--;
      else *curIndex = *(curIndex-1) + curOffset; // Set the length of the current jagged diagonal (plus scan sum)
    }
  }

  Values_.Resize(numMyNonzeros);
  Indices_.Resize(numMyNonzeros);

  int NumEntries;
  int * Indices = 0;
  double * Values =0;

  try { // If matrix is an Epetra_CrsMatrix, we can get data much more cheaply

    const Epetra_CrsMatrix & A = dynamic_cast<const Epetra_CrsMatrix &>(Matrix);

    for (int i1=0; i1<numMyRows; i1++) {

      A.ExtractMyRowView(i1, NumEntries, Values, Indices); // Get the current row
      int i = InvRowPerm_[i1]; // Determine permuted row location
      //cout << "i1, i, NumEntries = " << i1 <<" "<< i <<" "<< NumEntries << endl;
      for (int j=0; j< NumEntries; j++) {
	Values_[IndexOffset_[j]+i] = Values[j];
	Indices_[IndexOffset_[j]+i] = Indices[j];
      }
    }
  }
  catch (...) { // Otherwise just live with RowMatrix interface

  Epetra_SerialDenseVector curValues(NumJaggedDiagonals_);
  Epetra_IntSerialDenseVector curIndices(NumJaggedDiagonals_);
  Indices = curIndices.Values();
  Values = curValues.Values();
    for (int i1=0; i1<numMyRows; i1++) {
      Matrix.ExtractMyRowCopy(i1, NumJaggedDiagonals_, NumEntries, Values, Indices); // Get  current row based on the permutation
      int i = InvRowPerm_[i1]; // Determine permuted row location
      for (int j=0; j< NumEntries; j++) {
	Values_[IndexOffset_[j]+i] = Values[j];
	Indices_[IndexOffset_[j]+i] = Indices[j];
      }
    }
  }
}
//=============================================================================
int Epetra_JadMatrix::NumMyRowEntries(int MyRow, int & NumEntries) const {
  int i = InvRowPerm_[MyRow]; // Determine permuted row location
  NumEntries = Profile_[i]; // NNZ in current row
  return(0);
}
//=============================================================================
int Epetra_JadMatrix::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const {

  if(MyRow < 0 || MyRow >= NumMyRows_)
    EPETRA_CHK_ERR(-1); // Not in Row range

  int i = InvRowPerm_[MyRow]; // Determine permuted row location
  NumEntries = Profile_[i]; // NNZ in current row
  if(NumEntries > Length)
    EPETRA_CHK_ERR(-2); // Not enough space for copy. Needed size is passed back in NumEntries

  for (int j=0; j< NumEntries; j++) Values[j] = Values_[IndexOffset_[j]+i];
  for (int j=0; j< NumEntries; j++) Indices[j] = Indices_[IndexOffset_[j]+i];
  return(0);
}
//=============================================================================
int Epetra_JadMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  //
  // This function forms the product Y = A * Y or Y = A' * X
  //

  int NumVectors = X.NumVectors();
  if (NumVectors!=Y.NumVectors()) {
    EPETRA_CHK_ERR(-1); // Need same number of vectors in each MV
  }

  double** Xp = (double**) X.Pointers();
  double** Yp = (double**) Y.Pointers();
  int LDX = X.ConstantStride() ? X.Stride() : 0;
  int LDY = Y.ConstantStride() ? Y.Stride() : 0;
  UpdateImportVector(NumVectors); // Make sure Import and Export Vectors are compatible
  UpdateExportVector(NumVectors);

  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      EPETRA_CHK_ERR(ImportVector_->Import(X, *Importer(), Insert));
      Xp = (double**)ImportVector_->Pointers();
      LDX = ImportVector_->ConstantStride() ? ImportVector_->Stride() : 0;
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      Yp = (double**)ExportVector_->Pointers();
      LDY = ExportVector_->ConstantStride() ? ExportVector_->Stride() : 0;
    }

    // Do actual computation
    if (NumVectors==1)
      GeneralMV(TransA, *Xp, *Yp);
    else
      GeneralMM(TransA, Xp, LDX, Yp, LDY, NumVectors);
    if (Exporter()!=0) {
      Y.PutScalar(0.0);  // Make sure target is zero
      Y.Export(*ExportVector_, *Exporter(), Add); // Fill Y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!OperatorRangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Y.Reduce());
  }
  else { // Transpose operation
		

    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors

    if (Exporter()!=0) {
      EPETRA_CHK_ERR(ExportVector_->Import(X, *Exporter(), Insert));
      Xp = (double**)ExportVector_->Pointers();
      LDX = ExportVector_->ConstantStride() ? ExportVector_->Stride() : 0;
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      Yp = (double**)ImportVector_->Pointers();
      LDY = ImportVector_->ConstantStride() ? ImportVector_->Stride() : 0;
    }

    // Do actual computation
    if (NumVectors==1)
      GeneralMV(TransA, *Xp, *Yp);
    else
      GeneralMM(TransA, Xp, LDX, Yp, LDY, NumVectors);
    if (Importer()!=0) {
      Y.PutScalar(0.0);  // Make sure target is zero
      EPETRA_CHK_ERR(Y.Export(*ImportVector_, *Importer(), Add)); // Fill Y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!OperatorDomainMap().DistributedGlobal() && Comm().NumProc()>1)  EPETRA_CHK_ERR(Y.Reduce());
  }

  UpdateFlops(2*NumVectors*NumGlobalNonzeros64());
  return(0);
}
//=======================================================================================================
void Epetra_JadMatrix::GeneralMM(bool TransA, double ** X, int LDX, double ** Y, int LDY, int NumVectors) const {

  if (LDX==0 || LDY==0 || NumVectors==1) {// Can't unroll RHS if X or Y not strided
    for (int k=0; k<NumVectors; k++) GeneralMV(TransA, X[k], Y[k]);
  }
  else if (NumVectors==2) // Special 2 RHS case (does unrolling in both NumVectors and NumJaggedDiagonals)
    GeneralMM2RHS(TransA, X[0], LDX, Y[0], LDY);
  // Otherwise unroll RHS only
  else
    GeneralMM3RHS(TransA, X, LDX, Y, LDY, NumVectors);

  return;
}
//=======================================================================================================
void Epetra_JadMatrix::GeneralMM3RHS(bool TransA, double ** X, int ldx, double ** Y, int ldy, int NumVectors) const {

#ifdef _CRAY 
#define Pragma(S) _Pragma(S) 
#else 
#define Pragma(S) 
#endif

  // Routine for 3 or more RHS

  const double * Values = Values_.Values();
  const int * Indices = Indices_.Values();
  const int * IndexOffset = IndexOffset_.Values();
  const int * RowPerm = RowPerm_.Values();
  for (int j=0; j<NumVectors; j++) {
    double * y = Y[j];
    if (!TransA)
      for (int i=0; i<NumMyRows_; i++) y[i] = 0.0;
    else
      for (int i=0; i<NumMyCols_; i++) y[i] = 0.0;
  }

  int nv = NumVectors%5; if (nv==0) nv=5;
    double * x = X[0];
    double * y = Y[0];
 

  for (int k=0; k<NumVectors; k+=5) {
    
    for (int j=0; j<NumJaggedDiagonals_; j++) {
      const int * curIndices = Indices+IndexOffset[j];
      const double * curValues = Values+IndexOffset[j];
      int jaggedDiagonalLength = IndexOffset[j+1]-IndexOffset[j];
      switch (nv){
      case 1:
	{
	  if (!TransA) {
Pragma("_CRI ivdep")
	    for (int i=0; i<jaggedDiagonalLength; i++) {
	      int ix = curIndices[i];
	      int iy = RowPerm[i];
	      double val = curValues[i];
	      y[iy] += val*x[ix];
	    }
	  }
	  else {
Pragma("_CRI ivdep")
	    for (int i=0; i<jaggedDiagonalLength; i++) {
	      int iy = curIndices[i];
	      int ix = RowPerm[i];
	      double val = curValues[i];
	      y[iy] += val*x[ix];
	    }
	  }
	  break;
	}
      case 2:
	{
	  if (!TransA) {
Pragma("_CRI ivdep")
	    for (int i=0; i<jaggedDiagonalLength; i++) {
	      int ix = curIndices[i];
	      int iy = RowPerm[i];
	      double val = curValues[i];
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	    }
	  }
	  else {
Pragma("_CRI ivdep")
	    for (int i=0; i<jaggedDiagonalLength; i++) {
	      int iy = curIndices[i];
	      int ix = RowPerm[i];
	      double val = curValues[i];
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	    }
	  }
	  break;
	}
      case 3:
	{
	  if (!TransA) {
Pragma("_CRI ivdep")
	    for (int i=0; i<jaggedDiagonalLength; i++) {
	      int ix = curIndices[i];
	      int iy = RowPerm[i];
	      double val = curValues[i];
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	    }
	  }
	  else {
Pragma("_CRI ivdep")
	    for (int i=0; i<jaggedDiagonalLength; i++) {
	      int iy = curIndices[i];
	      int ix = RowPerm[i];
	      double val = curValues[i];
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	    }
	  }
	  break;
	}
      case 4:
	{
	  if (!TransA) {
Pragma("_CRI ivdep")
	    for (int i=0; i<jaggedDiagonalLength; i++) {
	      int ix = curIndices[i];
	      int iy = RowPerm[i];
	      double val = curValues[i];
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	    }
	  }
	  else {
Pragma("_CRI ivdep")
	    for (int i=0; i<jaggedDiagonalLength; i++) {
	      int iy = curIndices[i];
	      int ix = RowPerm[i];
	      double val = curValues[i];
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	    }
	  }
	  break;
	}
      case 5:
	{
	  if (!TransA) {
Pragma("_CRI ivdep")
	    for (int i=0; i<jaggedDiagonalLength; i++) {
	      int ix = curIndices[i];
	      int iy = RowPerm[i];
	      double val = curValues[i];
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	    }
	  }
	  else {
Pragma("_CRI ivdep")
	    for (int i=0; i<jaggedDiagonalLength; i++) {
	      int iy = curIndices[i];
	      int ix = RowPerm[i];
	      double val = curValues[i];
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	      iy+=ldy; ix+=ldx;
	      y[iy] += val*x[ix];
	    }
	  }
	  break;
	}
      }
    }
    x += nv*ldx;
    y += nv*ldy;
    nv = 5; // After initial remainder, we will always do 5 RHS
  }
  return;
}
//=======================================================================================================
void Epetra_JadMatrix::GeneralMM2RHS(bool TransA, double * x, int ldx, double * y, int ldy) const {

  // special 2 rhs case

  const double * Values = Values_.Values();
  const int * Indices = Indices_.Values();
  const int * IndexOffset = IndexOffset_.Values();
  const int * RowPerm = RowPerm_.Values();
  if (!TransA) 
    for (int i=0; i<NumMyRows_; i++) {
      y[i] = 0.0;
      y[i+ldy] = 0.0;
    }
  else
    for (int i=0; i<NumMyCols_; i++) {
      y[i] = 0.0;
      y[i+ldy] = 0.0;
    }

  int j = 0;
  while (j<NumJaggedDiagonals_) {
  int j0 = j;
  int jaggedDiagonalLength = IndexOffset[j+1]-IndexOffset[j];
    j++;
    // check if other diagonals have same length up to a max of 2
    while ((j<NumJaggedDiagonals_-1) && (IndexOffset[j+1]-IndexOffset[j]==jaggedDiagonalLength) && (j-j0<2)) j++;
    
    int numDiags = j-j0;
    assert(numDiags<3 && numDiags>0);
    assert(j<NumJaggedDiagonals_+1);
    
    switch (numDiags){
    case 1:
      {
	const int * curIndices = Indices+IndexOffset[j0];
	const double * curValues = Values+IndexOffset[j0];
	if (!TransA) {
Pragma("_CRI ivdep")
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    int ix = curIndices[i];
	    int iy = RowPerm[i];
	    y[iy] += curValues[i]*x[ix];
	    iy+=ldy; ix+=ldx;
	    y[iy] += curValues[i]*x[ix];
	  }
	}
	else {
Pragma("_CRI ivdep")
	  for (int i=0; i<jaggedDiagonalLength; i++){
	    int iy = curIndices[i];
	    int ix = RowPerm[i];
	    y[iy] += curValues[i]*x[ix];
	    iy+=ldy; ix+=ldx;
	    y[iy] += curValues[i]*x[ix];
	  }
	}
	break;
      }
    case 2:
      {
	const int * curIndices0 = Indices+IndexOffset[j0];
	const double * curValues0 = Values+IndexOffset[j0++];
	const int * curIndices1 = Indices+IndexOffset[j0];
	const double * curValues1 = Values+IndexOffset[j0];
	if (!TransA) {
Pragma("_CRI ivdep")
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    int ix0 = curIndices0[i];
	    int ix1 = curIndices1[i];
	    int iy = RowPerm[i];
	    y[iy] += 
	      curValues0[i]*x[ix0] +
	      curValues1[i]*x[ix1];
	    iy+=ldy; ix0+=ldx; ix1+=ldx;
	    y[iy] += 
	      curValues0[i]*x[ix0] +
	      curValues1[i]*x[ix1];
	  }
	}
	else {
Pragma("_CRI ivdep")
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    int iy0 = curIndices0[i];
	    int iy1 = curIndices1[i];
	    int ix = RowPerm[i];
	    double xval = x[ix];
	    y[iy0] += curValues0[i]*xval;
	    y[iy1] += curValues1[i]*xval;
	    ix+=ldx; iy0+=ldy; iy1+=ldy;
	    xval = x[ix];
	    y[iy0] += curValues0[i]*xval;
	    y[iy1] += curValues1[i]*xval;
	  }
	}
      }
      break;
    }
  }
  return;
}
//=======================================================================================================
void Epetra_JadMatrix::GeneralMV(bool TransA, double * x, double * y)  const {
  
  const double * Values = Values_.Values();
  const int * Indices = Indices_.Values();
  const int * IndexOffset = IndexOffset_.Values();
  const int * RowPerm = RowPerm_.Values();
  if (!TransA)
    for (int i=0; i<NumMyRows_; i++) y[i] = 0.0;
  else
    for (int i=0; i<NumMyCols_; i++) y[i] = 0.0;

  int j = 0;
  while (j<NumJaggedDiagonals_) {
  int j0 = j;
  int jaggedDiagonalLength = IndexOffset[j+1]-IndexOffset[j];
    j++;
    // check if other diagonals have same length up to a max of 5
    while ((j<NumJaggedDiagonals_-1) && (IndexOffset[j+1]-IndexOffset[j]==jaggedDiagonalLength) && (j-j0<5)) j++;
    
    int numDiags = j-j0;
    assert(numDiags<6 && numDiags>0);
    assert(j<NumJaggedDiagonals_+1);
    
    switch (numDiags){
    case 1:
      {
	const int * curIndices = Indices+IndexOffset[j0];
	const double * curValues = Values+IndexOffset[j0];
	if (!TransA) {
Pragma("_CRI ivdep")
	  for (int i=0; i<jaggedDiagonalLength; i++)
	    y[RowPerm[i]] += curValues[i]*x[curIndices[i]];
	}
	else {
Pragma("_CRI ivdep")
	  for (int i=0; i<jaggedDiagonalLength; i++)
	    y[curIndices[i]] += curValues[i]*x[RowPerm[i]];
	}
	break;
      }
    case 2:
      {
	const int * curIndices0 = Indices+IndexOffset[j0];
	const double * curValues0 = Values+IndexOffset[j0++];
	const int * curIndices1 = Indices+IndexOffset[j0];
	const double * curValues1 = Values+IndexOffset[j0];
	if (!TransA) {
Pragma("_CRI ivdep")
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    y[RowPerm[i]] += 
	      curValues0[i]*x[curIndices0[i]] +
	      curValues1[i]*x[curIndices1[i]];
	  }
	}
	else {
	  //Pragma("_CRI ivdep")  (Collisions possible)
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    double xval = x[RowPerm[i]];
	    y[curIndices0[i]] += curValues0[i]*xval;
	    y[curIndices1[i]] += curValues1[i]*xval;
	  }
	}
      }
      break;
    case 3:
      {
	const int * curIndices0 = Indices+IndexOffset[j0];
	const double * curValues0 = Values+IndexOffset[j0++];
	const int * curIndices1 = Indices+IndexOffset[j0];
	const double * curValues1 = Values+IndexOffset[j0++];
	const int * curIndices2 = Indices+IndexOffset[j0];
	const double * curValues2 = Values+IndexOffset[j0];
	if (!TransA) {
Pragma("_CRI ivdep")
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    y[RowPerm[i]] += 
	      curValues0[i]*x[curIndices0[i]] +
	      curValues1[i]*x[curIndices1[i]] +
	      curValues2[i]*x[curIndices2[i]];
	  }
	}
	else {
	  //Pragma("_CRI ivdep")  (Collisions possible)
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    double xval = x[RowPerm[i]];
	    y[curIndices0[i]] += curValues0[i]*xval;
	    y[curIndices1[i]] += curValues1[i]*xval;
	    y[curIndices2[i]] += curValues2[i]*xval;
	  }
	}
      }
      break;
    case 4:
      {
	const int * curIndices0 = Indices+IndexOffset[j0];
	const double * curValues0 = Values+IndexOffset[j0++];
	const int * curIndices1 = Indices+IndexOffset[j0];
	const double * curValues1 = Values+IndexOffset[j0++];
	const int * curIndices2 = Indices+IndexOffset[j0];
	const double * curValues2 = Values+IndexOffset[j0++];
	const int * curIndices3 = Indices+IndexOffset[j0];
	const double * curValues3 = Values+IndexOffset[j0];
	if (!TransA) {
Pragma("_CRI ivdep")
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    y[RowPerm[i]] += 
	      curValues0[i]*x[curIndices0[i]] +
	      curValues1[i]*x[curIndices1[i]] +
	      curValues2[i]*x[curIndices2[i]] +
	      curValues3[i]*x[curIndices3[i]];
	  }
	}
	else {
	  //Pragma("_CRI ivdep")  (Collisions possible)
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    double xval = x[RowPerm[i]];
	    y[curIndices0[i]] += curValues0[i]*xval;
	    y[curIndices1[i]] += curValues1[i]*xval;
	    y[curIndices2[i]] += curValues2[i]*xval;
	    y[curIndices3[i]] += curValues3[i]*xval;
	  }
	}
      }
      break;
    case 5:
      {
	const int * curIndices0 = Indices+IndexOffset[j0];
	const double * curValues0 = Values+IndexOffset[j0++];
	const int * curIndices1 = Indices+IndexOffset[j0];
	const double * curValues1 = Values+IndexOffset[j0++];
	const int * curIndices2 = Indices+IndexOffset[j0];
	const double * curValues2 = Values+IndexOffset[j0++];
	const int * curIndices3 = Indices+IndexOffset[j0];
	const double * curValues3 = Values+IndexOffset[j0++];
	const int * curIndices4 = Indices+IndexOffset[j0];
	const double * curValues4 = Values+IndexOffset[j0];
	if (!TransA) {
Pragma("_CRI ivdep")
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    y[RowPerm[i]] += 
	      curValues0[i]*x[curIndices0[i]] +
	      curValues1[i]*x[curIndices1[i]] +
	      curValues2[i]*x[curIndices2[i]] +
	      curValues3[i]*x[curIndices3[i]] +
	      curValues4[i]*x[curIndices4[i]];
	  }
	}
	else {
	  // Pragma("_CRI ivdep") (Collisions possible)
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    double xval = x[RowPerm[i]];
	    y[curIndices0[i]] += curValues0[i]*xval;
	    y[curIndices1[i]] += curValues1[i]*xval;
	    y[curIndices2[i]] += curValues2[i]*xval;
	    y[curIndices3[i]] += curValues3[i]*xval;
	    y[curIndices4[i]] += curValues4[i]*xval;
	  }
	}
      }
      break;
    }
  }
  return;
}
