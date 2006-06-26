
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

// At this point, we are not support reduced storage capabilities
#undef REDUCED_STORAGE_SUPPORT

// When we do, we will need to add a check for climits vs limits.h because some machine, esp. IRIX 
// platforms do not have climits.
#ifdef REDUCED_STORAGE_SUPPORT
#include <climits>
#endif

//==============================================================================
Epetra_JadMatrix::Epetra_JadMatrix(const Epetra_RowMatrix & Matrix) 
  : Comm_(Matrix.RowMatrixRowMap().Comm().Clone()),
    OperatorDomainMap_(Matrix.OperatorDomainMap()),
    OperatorRangeMap_(Matrix.OperatorRangeMap()),
    RowMatrixRowMap_(Matrix.RowMatrixRowMap()),
    RowMatrixColMap_(Matrix.RowMatrixColMap()),
    NumMyRows_(Matrix.NumMyRows()),
    NumMyCols_(Matrix.NumMyCols()),
    NumMyNonzeros_(Matrix.NumMyNonzeros()),
    NumGlobalNonzeros_(Matrix.NumGlobalNonzeros()),
    Values_(0),
    Indices_(0),
    IndexOffset_(0),
    Profile_(0),
    RowPerm_(0),
    InvRowPerm_(0),
    UseTranspose_(Matrix.UseTranspose()),
    HasNormInf_(Matrix.HasNormInf()),
    LowerTriangular_(Matrix.LowerTriangular()),
    UpperTriangular_(Matrix.UpperTriangular()),
    NumJaggedDiagonals_(Matrix.MaxNumEntries()),
    ImportVector_(0),
    ExportVector_(0),
    Importer_(0),
    Exporter_(0)
{
  if (!Matrix.Filled()) throw ReportError("Input matrix must have called FillComplete()", -1);
  Allocate(Matrix);
  SetLabel("Epetra::JadMatrix");
}

//==============================================================================
Epetra_JadMatrix::~Epetra_JadMatrix(){

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
      for (int j=0; j< NumEntries; j++) Values_[IndexOffset_[j]+i] = Values[i];
      if (CheckStructure)
	for (int j=0; j< NumEntries; j++) if (Indices_[IndexOffset_[j]+i] != Indices[i]) ierr = - 1;
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
      for (int j=0; j< NumEntries; j++) Values_[IndexOffset_[j]+i] = Values[i];
      if (CheckStructure)
	for (int j=0; j< NumEntries; j++) if (Indices_[IndexOffset_[j]+i] != Indices[i]) ierr = - 1;
    }
  }
  EPETRA_CHK_ERR(ierr);
  return(ierr);
}

//==============================================================================
int Epetra_JadMatrix::Allocate(const Epetra_RowMatrix & Matrix) {

  // Check if non-trivial import/export operators
  if (!(Matrix.RowMatrixRowMap().SameAs(Matrix.OperatorRangeMap()))) 
    Exporter_ = new Epetra_Export(Matrix.RowMatrixRowMap(), Matrix.OperatorRangeMap());
  
  if (Matrix.RowMatrixImporter()!=0) 
    Importer_ = new Epetra_Import(Matrix.RowMatrixColMap(), Matrix.OperatorDomainMap());

  // Allocate IndexOffset storage

  IndexOffset_.Resize(NumJaggedDiagonals_+1);

  // Next compute permutation of rows
  RowPerm_.Resize(NumMyRows_);
  InvRowPerm_.Resize(NumMyRows_);
  Profile_.Resize(NumMyRows_);
  for (int i=0; i<NumMyRows_; i++) {
    int NumEntries;
    Matrix.NumMyRowEntries(i, NumEntries);
    Profile_[i] = NumEntries;
    RowPerm_[i] = i;
  }

  Epetra_Util sorter;
  int * RowPerm = RowPerm_.Values();
  sorter.Sort(false, NumMyRows_, Profile_.Values(), 0, 0, 1, &RowPerm);
  for (int i=0; i<NumMyRows_; i++) InvRowPerm_[RowPerm[i]] = i; // Compute inverse row permutation

  // Now build IndexOffsets:  These contain the lengths of the jagged diagonals

  for (int i=0; i<NumJaggedDiagonals_; i++) IndexOffset_[i] = 0;

  int curOffset = NumMyRows_;
  int * curIndex = IndexOffset_.Values(); // point to first index (will be incremented immediately below)
  for (int i=1; i<NumJaggedDiagonals_+1; i++) {
    curIndex++;
    while (*curIndex==0) {
      if (Profile_[curOffset-1]<i) curOffset--;
      else *curIndex = *(curIndex-1) + curOffset; // Set the length of the current jagged diagonal (plus scan sum)
    }
  }

  Values_.Resize(NumMyNonzeros_);
  Indices_.Resize(NumMyNonzeros_);

  int NumEntries;
  int * Indices = 0;
  double * Values =0;

  try { // If matrix is an Epetra_CrsMatrix, we can get data much more cheaply

    const Epetra_CrsMatrix & A = dynamic_cast<const Epetra_CrsMatrix &>(Matrix);

    for (int i1=0; i1<NumMyRows_; i1++) {

      EPETRA_CHK_ERR(A.ExtractMyRowView(i1, NumEntries, Values, Indices)); // Get the current row
      int i = InvRowPerm_[i1]; // Determine permuted row location
      for (int j=0; j< NumEntries; j++) Values_[IndexOffset_[j]+i] = Values[j];
      for (int j=0; j< NumEntries; j++) Indices_[IndexOffset_[j]+i] = Indices[j];
    }
  }
  catch (...) { // Otherwise just live with RowMatrix interface

  Epetra_SerialDenseVector curValues(NumJaggedDiagonals_);
  Epetra_IntSerialDenseVector curIndices(NumJaggedDiagonals_);
  Indices = curIndices.Values();
  Values = curValues.Values();
    for (int i1=0; i1<NumMyRows_; i1++) {
      EPETRA_CHK_ERR(Matrix.ExtractMyRowCopy(i1, NumJaggedDiagonals_, NumEntries, Values, Indices)); // Get  current row based on the permutation
      int i = InvRowPerm_[i1]; // Determine permuted row location
	for (int j=0; j< NumEntries; j++) Values_[IndexOffset_[j]+i] = Values[j];
	for (int j=0; j< NumEntries; j++) Indices_[IndexOffset_[j]+i] = Indices[j];
    }
  }
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
int Epetra_JadMatrix::NumMyRowEntries(int MyRow, int & NumEntries) const {

  // Crude implementation in terms of ExtractMyRowCopy

  Epetra_SerialDenseVector Values(MaxNumEntries());
  Epetra_IntSerialDenseVector Indices(MaxNumEntries());
  EPETRA_CHK_ERR(ExtractMyRowCopy(MyRow, MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
  
  
  return(1); // Return 1 because this is an expensive implementation (in runtime)
}
//=============================================================================
int Epetra_JadMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const {

  if(!RowMatrixRowMap().SameAs(Diagonal.Map())) 
    EPETRA_CHK_ERR(-2); // Maps must be the same

  // Crude implementation in terms of ExtractMyRowCopy

  Epetra_SerialDenseVector Values(MaxNumEntries());
  Epetra_IntSerialDenseVector Indices(MaxNumEntries());
  int NumEntries;

  for(int i = 0; i < NumMyRows_; i++) {
    EPETRA_CHK_ERR(ExtractMyRowCopy(RowPerm_[i], MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
    int ii = RowMatrixRowMap().GID(i);
    
    Diagonal[i] = 0.0;
    for(int j = 0; j < NumEntries; j++) {
      if(ii == RowMatrixColMap().GID(Indices[j])) {
	Diagonal[i] = Values[j];
	break;
      }
    }
  }
  return(0);
}
//=============================================================================
int Epetra_JadMatrix::InvRowSums(Epetra_Vector & x) const {
  int ierr = 0;
  int i, j;
  Epetra_SerialDenseVector Values(MaxNumEntries());
  Epetra_IntSerialDenseVector Indices(MaxNumEntries());
  int NumEntries;
  x.PutScalar(0.0); // Make sure we sum into a vector of zeros.
  double * xp = (double*)x.Values();
  if (OperatorRangeMap().SameAs(x.Map()) && Exporter() != 0) {
    Epetra_Vector x_tmp(RowMatrixRowMap());
    x_tmp.PutScalar(0.0);
    double * x_tmp_p = (double*)x_tmp.Values();
    for (i=0; i < NumMyRows_; i++) {
      EPETRA_CHK_ERR(ExtractMyRowCopy(RowPerm_[i], MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
      for (j=0; j < NumEntries; j++)  x_tmp_p[i] += fabs(Values[j]);
    }
    EPETRA_CHK_ERR(x.Export(x_tmp, *Exporter(), Add)); //Export partial row sums to x.
    int myLength = x.MyLength();
    for (i=0; i<myLength; i++) { 
      if (xp[i]<Epetra_MinDouble) {
        if (xp[i]==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
        else if (ierr!=1) ierr = 2;
        xp[i] = Epetra_MaxDouble;
      }
      else
        xp[i] = 1.0/xp[i];
    }
  }
  else if (RowMatrixRowMap().SameAs(x.Map())) {
    for (i=0; i < NumMyRows_; i++) {
      EPETRA_CHK_ERR(ExtractMyRowCopy(RowPerm_[i], MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
      double scale = 0.0;
      for (j=0; j < NumEntries; j++) scale += fabs(Values[j]);
      if (scale<Epetra_MinDouble) {
        if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
        else if (ierr!=1) ierr = 2;
        xp[i] = Epetra_MaxDouble;
      }
      else
        xp[i] = 1.0/scale;
    }
  }
  else { // x.Map different than both RowMatrixRowMap() and OperatorRangeMap()
    EPETRA_CHK_ERR(-2); // The map of x must be the RowMap or RangeMap of A.
  }
  EPETRA_CHK_ERR(ierr);  
  UpdateFlops(NumGlobalNonzeros());
  return(0);
}
//=============================================================================
int Epetra_JadMatrix::LeftScale(const Epetra_Vector & x) {
  double* xp = 0;
  double curValue;
  int curRowIndex, curColIndex;
  if(OperatorRangeMap().SameAs(x.Map()) && Exporter() != 0) {
    Epetra_Vector xtmp(RowMatrixRowMap());
    xtmp.Import(x,*Exporter(),Insert);
    for (int i=0; i<NumMyNonzeros_; i++) {
      ExtractMyEntryView(i, &curValue, curRowIndex, curColIndex);
      curValue *= xtmp[curRowIndex];
    }
  }
  else if (RowMatrixRowMap().SameAs(x.Map()))
    for (int i=0; i<NumMyNonzeros_; i++) {
      ExtractMyEntryView(i, &curValue, curRowIndex, curColIndex);
      curValue *= x[curRowIndex];
    }
  else {
    EPETRA_CHK_ERR(-2); // The Map of x must be the RowMap or RangeMap of A.
  }
  UpdateFlops(NumGlobalNonzeros());
  return(0);
}
//=============================================================================
int Epetra_JadMatrix::InvColSums(Epetra_Vector & x) const {
  int ierr = 0;
  int i, j;
  Epetra_SerialDenseVector Values(MaxNumEntries());
  Epetra_IntSerialDenseVector Indices(MaxNumEntries());
  int NumEntries;
  int MapNumMyElements = x.Map().NumMyElements();
  x.PutScalar(0.0); // Make sure we sum into a vector of zeros.
  double* xp = (double*)x.Values();
  if(OperatorDomainMap().SameAs(x.Map()) && Importer() != 0) {
    Epetra_Vector x_tmp(RowMatrixColMap());
    x_tmp.PutScalar(0.0);
    double * x_tmp_p = (double*)x_tmp.Values();
    for(i = 0; i < NumMyRows_; i++) {
      EPETRA_CHK_ERR(ExtractMyRowCopy(RowPerm_[i], MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
      for(j = 0; j < NumEntries; j++) 
        x_tmp_p[Indices[j]] += fabs(Values[j]);
    }
    EPETRA_CHK_ERR(x.Export(x_tmp, *Importer(), Add)); // Fill x with partial column sums
  }
  else if(RowMatrixColMap().SameAs(x.Map())) {
    for(i = 0; i < NumMyRows_; i++) {
      EPETRA_CHK_ERR(ExtractMyRowCopy(RowPerm_[i], MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
      for(j = 0; j < NumEntries; j++) 
        xp[Indices[j]] += fabs(Values[j]);
    }
  }
  else { //x.Map different than both RowMatrixColMap() and OperatorDomainMap()
    EPETRA_CHK_ERR(-2); // x must have the same distribution as the domain of A
  }

  // Invert values, don't allow them to get too large
  for(i = 0; i < MapNumMyElements; i++) {
    double scale = xp[i];
    if(scale < Epetra_MinDouble) {
      if(scale == 0.0) 
	ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
      else if(ierr != 1) 
	ierr = 2;
      xp[i] = Epetra_MaxDouble;
    }
    else
      xp[i] = 1.0 / scale;
  }

  EPETRA_CHK_ERR(ierr);
  UpdateFlops(NumGlobalNonzeros());
  return(0);
}
//=============================================================================
int Epetra_JadMatrix::RightScale(const Epetra_Vector & x) {
  double* xp = 0;
  double curValue;
  int curRowIndex, curColIndex;
  if(OperatorDomainMap().SameAs(x.Map()) && Importer() != 0) {
    Epetra_Vector xtmp(RowMatrixColMap());
    xtmp.Import(x,*Importer(),Insert);
    for (int i=0; i<NumMyNonzeros_; i++) {
      ExtractMyEntryView(i, &curValue, curRowIndex, curColIndex);
      curValue *= xtmp[curColIndex];
    }
  }
  else if (RowMatrixColMap().SameAs(x.Map()))
    for (int i=0; i<NumMyNonzeros_; i++) {
      ExtractMyEntryView(i, &curValue, curRowIndex, curColIndex);
      curValue *= x[curColIndex];
    }
  else {
    EPETRA_CHK_ERR(-2); // The Map of x must be the RowMap or RangeMap of A.
  }
  UpdateFlops(NumGlobalNonzeros());
  return(0);
}
//=============================================================================
double Epetra_JadMatrix::NormInf() const {
  double NormInf = -1.0;
  Epetra_SerialDenseVector Values(MaxNumEntries());
  Epetra_IntSerialDenseVector Indices(MaxNumEntries());
  int NumEntries;
  Epetra_Vector x(RowMatrixRowMap()); // Need temp vector for row sums
  for(int i = 0; i < NumMyRows_; i++) {
    x[i] = 0.0;
    EPETRA_CHK_ERR(ExtractMyRowCopy(RowPerm_[i], MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
    for(int j = 0; j < NumEntries; j++) 
      x[i] += fabs(Values[j]);
  }

  // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
  if(Exporter() != 0) {
    Epetra_Vector xtmp(OperatorRangeMap()); // Create temporary import vector if needed
    xtmp.Export(x,*Exporter(),Add);
    xtmp.MaxValue(&NormInf); // This is the NormInf
  }
  else
    x.MaxValue(&NormInf); // Find max
  UpdateFlops(NumGlobalNonzeros());
  return(NormInf);
}
//=============================================================================
double Epetra_JadMatrix::NormOne() const {
  double NormOne = -1.0;
  Epetra_SerialDenseVector Values(MaxNumEntries());
  Epetra_IntSerialDenseVector Indices(MaxNumEntries());
  int NumEntries;
  Epetra_Vector x(RowMatrixColMap()); // Need temp vector for column sums
  for(int i = 0; i < NumMyRows_; i++) {
    x[i] = 0.0;
    EPETRA_CHK_ERR(ExtractMyRowCopy(RowPerm_[i], MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
    for(int j = 0; j < NumEntries; j++) 
      x[Indices[i]] += fabs(Values[j]);
  }
  
  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  if(Importer() != 0) {
    Epetra_Vector xtmp(OperatorDomainMap()); // Create temporary import vector if needed
    xtmp.Export(x,*Importer(),Add);
    xtmp.MaxValue(&NormOne); // This is the NormOne
  }
  else
    x.MaxValue(&NormOne); // Find max

  UpdateFlops(NumGlobalNonzeros());
  return(NormOne);
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

  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  return(0);
}
//=======================================================================================================
void Epetra_JadMatrix::UpdateImportVector(int NumVectors) const {
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
void Epetra_JadMatrix::UpdateExportVector(int NumVectors) const {
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
#pragma _CRI ivdep
	    for (int i=0; i<jaggedDiagonalLength; i++) {
	      int ix = curIndices[i];
	      int iy = RowPerm[i];
	      double val = curValues[i];
	      y[iy] += val*x[ix];
	    }
	  }
	  else {
#pragma _CRI ivdep
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
#pragma _CRI ivdep
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
#pragma _CRI ivdep
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
#pragma _CRI ivdep
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
#pragma _CRI ivdep
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
#pragma _CRI ivdep
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
#pragma _CRI ivdep
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
#pragma _CRI ivdep
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
#pragma _CRI ivdep
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
#pragma _CRI ivdep
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    int ix = curIndices[i];
	    int iy = RowPerm[i];
	    y[iy] += curValues[i]*x[ix];
	    iy+=ldy; ix+=ldx;
	    y[iy] += curValues[i]*x[ix];
	  }
	}
	else {
#pragma _CRI ivdep
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
#pragma _CRI ivdep
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
#pragma _CRI ivdep
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
#pragma _CRI ivdep
	  for (int i=0; i<jaggedDiagonalLength; i++)
	    y[RowPerm[i]] += curValues[i]*x[curIndices[i]];
	}
	else {
#pragma _CRI ivdep
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
#pragma _CRI ivdep
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    y[RowPerm[i]] += 
	      curValues0[i]*x[curIndices0[i]] +
	      curValues1[i]*x[curIndices1[i]];
	  }
	}
	else {
	  //#pragma _CRI ivdep  (Collisions possible)
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
#pragma _CRI ivdep
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    y[RowPerm[i]] += 
	      curValues0[i]*x[curIndices0[i]] +
	      curValues1[i]*x[curIndices1[i]] +
	      curValues2[i]*x[curIndices2[i]];
	  }
	}
	else {
	  //#pragma _CRI ivdep  (Collisions possible)
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
#pragma _CRI ivdep
	  for (int i=0; i<jaggedDiagonalLength; i++) {
	    y[RowPerm[i]] += 
	      curValues0[i]*x[curIndices0[i]] +
	      curValues1[i]*x[curIndices1[i]] +
	      curValues2[i]*x[curIndices2[i]] +
	      curValues3[i]*x[curIndices3[i]];
	  }
	}
	else {
	  //#pragma _CRI ivdep  (Collisions possible)
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
#pragma _CRI ivdep
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
	  // #pragma _CRI ivdep (Collisions possible)
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
void Epetra_JadMatrix::Print(ostream& os) const {

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
      for (int i1=0; i1<NumMyRows1; i1++) {
	int i = InvRowPerm_[i1];
	int Row = RowMap->GID(i);; // Get global row number
	
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
