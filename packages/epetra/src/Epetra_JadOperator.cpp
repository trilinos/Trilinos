
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

#include "Epetra_JadOperator.h"
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
#include <climits>

//==============================================================================
Epetra_JadOperator::Epetra_JadOperator(const Epetra_RowMatrix & Matrix, bool UseFloats, bool UseShorts) 
  : Comm_(Matrix.RowMatrixRowMap().Comm().Clone()),
    OperatorDomainMap_(Matrix.OperatorDomainMap()),
    OperatorRangeMap_(Matrix.OperatorRangeMap()),
    NumMyRows_(Matrix.NumMyRows()),
    NumMyCols_(Matrix.NumMyCols()),
    NumMyNonzeros_(Matrix.NumMyNonzeros()),
    NumGlobalNonzeros_(Matrix.NumGlobalNonzeros()),
    Values_(0),
    FloatValues_(0),
    Indices_(0),
    ShortIndices_(0),
    IndexOffset_(0),
    RowPerm_(0),
    UseTranspose_(Matrix.UseTranspose()),
    HasNormInf_(Matrix.HasNormInf()),
    UsingFloats_(UseFloats),
    UsingShorts_(UseShorts),
    NumJaggedDiagonals_(Matrix.MaxNumEntries()),
    ImportVector_(0),
    ExportVector_(0),
    Importer_(0),
    Exporter_(0)
{
  if (!Matrix.Filled()) throw ReportError("Input matrix must have called FillComplete()", -1);
  Allocate(Matrix, UseFloats);
  SetLabel("Epetra::JadOperator");
}

//==============================================================================
Epetra_JadOperator::~Epetra_JadOperator(){

  if (FloatValues_!=0) delete [] FloatValues_;
  if (ShortIndices_!=0) delete [] ShortIndices_;

  if (ImportVector_!=0) delete ImportVector_;
  ImportVector_=0;
  if (ExportVector_!=0) delete ExportVector_;
  ExportVector_=0;
  if (Importer_!=0) delete Importer_;
  Importer_=0;
  if (Exporter_!=0) delete Exporter_;
  Exporter_=0;
  delete Comm_;
}

//==============================================================================
int Epetra_JadOperator::UpdateValues(const Epetra_RowMatrix & Matrix, bool CheckStructure) {

  int NumEntries;
  int * Indices = 0;
  double * Values =0;

  int ierr = 0;

  try { // If matrix is an Epetra_CrsMatrix, we can get date much more cheaply

    const Epetra_CrsMatrix & A = dynamic_cast<const Epetra_CrsMatrix &>(Matrix);

    for (int i=0; i<NumMyRows_; i++) {

      EPETRA_CHK_ERR(A.ExtractMyRowView(RowPerm_[i], NumEntries, Values, Indices)); // Get the current row based on the permutation
      if (!UsingFloats_) 
	for (int j=0; j< NumEntries; j++) Values_[IndexOffset_[j]+i] = Values[i];
      else
	for (int j=0; j< NumEntries; j++) FloatValues_[IndexOffset_[j]+i] = (float) Values[i];
      if (CheckStructure) {
	if (!UsingShorts_) 
	  for (int j=0; j< NumEntries; j++) if (Indices_[IndexOffset_[j]+i] != Indices[i]) ierr = - 1;
	else
	  for (int j=0; j< NumEntries; j++) if (ShortIndices_[IndexOffset_[j]+i] != (unsigned short) Indices[i]) ierr = - 1;
      }
    }
  }
  catch (...) { // Otherwise just live with RowMatrix interface

  Epetra_SerialDenseVector curValues(NumJaggedDiagonals_);
  Epetra_IntSerialDenseVector curIndices(NumJaggedDiagonals_);
  Indices = curIndices.Values();
  Values = curValues.Values();
    for (int i=0; i<NumMyRows_; i++) {
      EPETRA_CHK_ERR(Matrix.ExtractMyRowCopy(RowPerm_[i], NumJaggedDiagonals_, NumEntries, Values, Indices)); // Get current row based on the permutation
      if (!UsingFloats_) 
	for (int j=0; j< NumEntries; j++) Values_[IndexOffset_[j]+i] = Values[i];
      else
	for (int j=0; j< NumEntries; j++) FloatValues_[IndexOffset_[j]+i] = (float) Values[i];
      if (CheckStructure) {
	if (!UsingShorts_) 
	  for (int j=0; j< NumEntries; j++) if (Indices_[IndexOffset_[j]+i] != Indices[i]) ierr = - 1;
	else
	  for (int j=0; j< NumEntries; j++) if (ShortIndices_[IndexOffset_[j]+i] != (unsigned short) Indices[i]) ierr = - 1;
      }
    }
  }
}

//==============================================================================
int Epetra_JadOperator::Allocate(const Epetra_RowMatrix & Matrix, bool UseFloats) {

  // Check if non-trivial import/export operators
  if (!(Matrix.RowMatrixRowMap().SameAs(Matrix.OperatorRangeMap()))) 
    Exporter_ = new Epetra_Export(Matrix.RowMatrixRowMap(), Matrix.OperatorRangeMap());
  
  if (Matrix.RowMatrixImporter()!=0) 
    Importer_ = new Epetra_Import(Matrix.RowMatrixColMap(), Matrix.OperatorDomainMap());

  // Allocate IndexOffset storage

  IndexOffset_.Resize(NumJaggedDiagonals_+1);

  // Next compute permutation of rows
  RowPerm_.Resize(NumMyRows_);
  Epetra_IntSerialDenseVector profile(NumMyRows_);
  for (int i=0; i<NumMyRows_; i++) {
    int NumEntries;
    Matrix.NumMyRowEntries(i, NumEntries);
    profile[i] = NumEntries;
    RowPerm_[i] = i;
  }

  Epetra_Util sorter;
  int * RowPerm = RowPerm_.Values();
  sorter.Sort(false, NumMyRows_, profile.Values(), 0, 0, 1, &RowPerm);

  // Now build IndexOffsets:  These contain the lengths of the jagged diagonals

  for (int i=0; i<NumJaggedDiagonals_; i++) IndexOffset_[i] = 0;

  int curOffset = NumMyRows_;
  int * curIndex = IndexOffset_.Values(); // point to first index (will be incremented immediately below)
  for (int i=1; i<NumJaggedDiagonals_+1; i++) {
    curIndex++;
    while (*curIndex==0) {
      if (profile[curOffset-1]<i) curOffset--;
      else *curIndex = *(curIndex-1) + curOffset; // Set the length of the current jagged diagonal (plus scan sum)
    }
  }

  // Next determine how to handle values.  Two possibilities:
  // 1) UseFloats is false, so we copy values into a contiguous double array.
  // 3) UseFloats is true so we create a single precision copy of matrix values.

  if (!UsingFloats_)  // Allocate storage in Values_
    Values_.Resize(NumMyNonzeros_);
  else // UseFloats == true
    FloatValues_ = new float[NumMyNonzeros_];

  // Next determine how to handle integers.  Two possibilities:
  // 1) Local column range is within the range of unsigned short ints, so we copy the indices to an array of unsigned shorts.
  // 2) Local column range is outside range of unsigned shorts and we copy to array of ints.
  // In both cases we embed the nonzero count per row into the index array.

  if (Matrix.NumMyCols()<=USHRT_MAX && UsingShorts_) UsingShorts_ = true;

  if (!UsingShorts_)
    Indices_.Resize(NumMyNonzeros_);
  else // Matrix.NumMyCols()<=USHRT_MAX
    ShortIndices_ = new unsigned short[NumMyNonzeros_];

  int NumEntries;
  int * Indices = 0;
  double * Values =0;

  try { // If matrix is an Epetra_CrsMatrix, we can get date much more cheaply

    const Epetra_CrsMatrix & A = dynamic_cast<const Epetra_CrsMatrix &>(Matrix);

    for (int i=0; i<NumMyRows_; i++) {

      EPETRA_CHK_ERR(A.ExtractMyRowView(RowPerm_[i], NumEntries, Values, Indices)); // Get the current row based on the permutation
      if (!UsingFloats_) 
	for (int j=0; j< NumEntries; j++) Values_[IndexOffset_[j]+i] = Values[j];
      else
	for (int j=0; j< NumEntries; j++) FloatValues_[IndexOffset_[j]+i] = (float) Values[j];
      if (!UsingShorts_) 
	for (int j=0; j< NumEntries; j++) Indices_[IndexOffset_[j]+i] = Indices[j];
      else
	for (int j=0; j< NumEntries; j++) ShortIndices_[IndexOffset_[j]+i] = (unsigned short) Indices[j];
    }
  }
  catch (...) { // Otherwise just live with RowMatrix interface

  Epetra_SerialDenseVector curValues(NumJaggedDiagonals_);
  Epetra_IntSerialDenseVector curIndices(NumJaggedDiagonals_);
  Indices = curIndices.Values();
  Values = curValues.Values();
    for (int i=0; i<NumMyRows_; i++) {
      EPETRA_CHK_ERR(Matrix.ExtractMyRowCopy(RowPerm_[i], NumJaggedDiagonals_, NumEntries, Values, Indices)); // Get  current row based on the permutation
      if (!UsingFloats_) 
	for (int j=0; j< NumEntries; j++) Values_[IndexOffset_[j]+i] = Values[j];
      else
	for (int j=0; j< NumEntries; j++) FloatValues_[IndexOffset_[j]+i] = (float) Values[j];
      if (!UsingShorts_) 
	for (int j=0; j< NumEntries; j++) Indices_[IndexOffset_[j]+i] = Indices[j];
      else
	for (int j=0; j< NumEntries; j++) ShortIndices_[IndexOffset_[j]+i] = (unsigned short) Indices[j];
    }
  }
  return(0);
}
//=============================================================================
int Epetra_JadOperator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
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

  bool TransA = UseTranspose_;
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

  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  return(0);
}
//=======================================================================================================
void Epetra_JadOperator::UpdateImportVector(int NumVectors) const {
  if(Importer() != 0) {
    if(ImportVector_ != 0) {
      if(ImportVector_->NumVectors() != NumVectors) {
     delete ImportVector_;
     ImportVector_= 0;
      }
    }
    if(ImportVector_ == 0)
      ImportVector_ = new Epetra_MultiVector(Importer_->TargetMap(),NumVectors); // Create import vector if needed
  }
  return;
}
//=======================================================================================================
void Epetra_JadOperator::UpdateExportVector(int NumVectors) const {
  if(Exporter() != 0) {
    if(ExportVector_ != 0) {
      if(ExportVector_->NumVectors() != NumVectors) {
     delete ExportVector_;
     ExportVector_= 0;
      }
    }
    if(ExportVector_ == 0)
      ExportVector_ = new Epetra_MultiVector(Exporter_->SourceMap(),NumVectors); // Create Export vector if needed
  }
  return;
}
//=======================================================================================================
void Epetra_JadOperator::GeneralMM(bool TransA, double ** X, int LDX, double ** Y, int LDY, int NumVectors) const {

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
void Epetra_JadOperator::GeneralMM3RHS(bool TransA, double ** X, int ldx, double ** Y, int ldy, int NumVectors) const {

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
	    for (int i=0; i<jaggedDiagonalLength; i++) {
	      int ix = curIndices[i];
	      int iy = RowPerm[i];
	      double val = curValues[i];
	      y[iy] += val*x[ix];
	    }
	  }
	  else {
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
void Epetra_JadOperator::GeneralMM2RHS(bool TransA, double * x, int ldx, double * y, int ldy) const {

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
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    int ix = curIndices[i];
	    int iy = RowPerm[i];
	    y[iy] += curValues[i]*x[ix];
	    iy+=ldy; ix+=ldx;
	    y[iy] += curValues[i]*x[ix];
	  }
	}
	else {
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
void Epetra_JadOperator::GeneralMV(bool TransA, double * x, double * y)  const {
  
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
	  for (int i=0; i<jaggedDiagonalLength; i++)
	    y[RowPerm[i]] += curValues[i]*x[curIndices[i]];
	}
	else {
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
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    y[RowPerm[i]] += 
	      curValues0[i]*x[curIndices0[i]] +
	      curValues1[i]*x[curIndices1[i]];
	  }
	}
	else {
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
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    y[RowPerm[i]] += 
	      curValues0[i]*x[curIndices0[i]] +
	      curValues1[i]*x[curIndices1[i]] +
	      curValues2[i]*x[curIndices2[i]];
	  }
	}
	else {
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
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    y[RowPerm[i]] += 
	      curValues0[i]*x[curIndices0[i]] +
	      curValues1[i]*x[curIndices1[i]] +
	      curValues2[i]*x[curIndices2[i]] +
	      curValues3[i]*x[curIndices3[i]];
	  }
	}
	else {
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
//=======================================================================================================
void Epetra_JadOperator::Print(ostream& os) const {

  const Epetra_BlockMap * RowMap;
  const Epetra_BlockMap * ColMap;
  if (Importer_==0) 
    ColMap = &(OperatorDomainMap());
  else
    ColMap = &(Importer_->TargetMap());
  if (Exporter_==0) 
    RowMap = &(OperatorRangeMap());
  else
    RowMap = &(Exporter_->SourceMap());

  int MyPID = RowMap->Comm().MyPID();
  int NumProc = RowMap->Comm().NumProc();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      if (MyPID==0) {
	os <<    "Number of Global Nonzeros     = "; os << NumGlobalNonzeros(); os << endl;
      }
			
      os <<  "\nNumber of My Rows               = "; os << NumMyRows_; os << endl;
      os <<    "Number of My Jagged Diagonals   = "; os << NumJaggedDiagonals_; os << endl;
      os <<    "Number of My Nonzeros           = "; os << NumMyNonzeros_; os << endl; os << endl;

      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    Comm().Barrier();
    Comm().Barrier();
    Comm().Barrier();
  }
	
  {for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumMyRows1 = NumMyRows_;
			
      if (MyPID==0) {
	os.width(8);
	os <<  "   Processor ";
	os.width(10);
	os <<  "   Row Index ";
	os.width(10);
	os <<  "   Col Index ";
	os.width(20);
	os <<  "   Value     ";
	os << endl;
      }
      for (int i=0; i<NumMyRows1; i++) {
	int Row = RowMap->GID(RowPerm_[i]);; // Get global row number
	
	for (int j = 0; j < NumJaggedDiagonals_ ; j++) {   
	  if (IndexOffset_[j+1]-IndexOffset_[j]>i) {
	    int Index = ColMap->GID(Indices_[IndexOffset_[j]+i]);
	    double Value = Values_[IndexOffset_[j]+i];
	    os.width(8);
	    os <<  MyPID ; os << "    ";	
	    os.width(10);
	    os <<  Row ; os << "    ";	
	    os.width(10);
	    os <<  Index; os << "    ";
	    os.width(20);
	    os <<  Value; os << "    ";
	    os << endl;
	  }
	}
      }
		             
      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    Comm().Barrier();
    Comm().Barrier();
    Comm().Barrier();
  }}
	
  return;
}
