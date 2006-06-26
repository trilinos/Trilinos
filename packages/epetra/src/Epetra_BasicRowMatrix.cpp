
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

#include "Epetra_BasicRowMatrix.h"
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
Epetra_BasicRowMatrix::Epetra_BasicRowMatrix(const Epetra_Map & RowMap, const Epetra_Map & ColMap) 
  : Comm_(RowMap.Comm().Clone()),
    OperatorDomainMap_(RowMap()),
    OperatorRangeMap_(RowMap()),
    RowMatrixRowMap_(RowMap()),
    RowMatrixColMap_(ColMap()),
    NumMyNonzeros_(0),
    NumGlobalNonzeros_(0),
    MaxNumEntries_(0),
    UseTranspose_(false),
    HasNormInf_(true),
    LowerTriangular_(false),
    UpperTriangular_(false),
    ImportVector_(0),
    ExportVector_(0),
    Importer_(0),
    Exporter_(0)
{
  Setup();
  SetLabel("Epetra::BasicRowMatrix");
}

//==============================================================================
Epetra_BasicRowMatrix::Epetra_BasicRowMatrix(const Epetra_Map & RowMap, const Epetra_Map & ,  Epetra_Map & DomainMap, const Epetra_Map & RangeMap) 
  : Comm_(RowMap.Comm().Clone()),
    OperatorDomainMap_(DomainMap()),
    OperatorRangeMap_(RangeMap()),
    RowMatrixRowMap_(RowMap()),
    RowMatrixColMap_(ColMap()),
    NumMyNonzeros_(0),
    NumGlobalNonzeros_(0),
    MaxNumEntries_(0),
    UseTranspose_(false),
    HasNormInf_(true),
    LowerTriangular_(true),
    UpperTriangular_(true),
    ImportVector_(0),
    ExportVector_(0),
    Importer_(0),
    Exporter_(0)
{
  Setup();
  SetLabel("Epetra::BasicRowMatrix");
}

//==============================================================================
Epetra_BasicRowMatrix::~Epetra_BasicRowMatrix(){

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
int Epetra_BasicRowMatrix::Setup() {

  // Check if non-trivial import/export operators
  if (!(RowMatrixRowMap().SameAs(OperatorRangeMap()))) 
    Exporter_ = new Epetra_Export(RowMatrixRowMap(), OperatorRangeMap());
  
  if (!(RowMatrixColMap().SamesAs(OperatorDomainMap())))
    Importer_ = new Epetra_Import(RowMatrixColMap(), OperatorDomainMap());

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
int Epetra_BasicRowMatrix::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const {

  return(-1);
}
//=============================================================================
int Epetra_BasicRowMatrix::NumMyRowEntries(int MyRow, int & NumEntries) const {

  return(-1);
}
//=============================================================================
int Epetra_BasicRowMatrix::MaxNumEntries() const {

  int NumMyRows = RowMatrixRowMap().NumMyPoints();
  MaxNumEntries_ = 0;
  int NumEntries = 0;
  for (int i=0; i<NumMyRows; i++) {
    EPETRA_CHK_ERR(ExtractMyRowCopy(MyRow, NumEntries));
    if (NumEntries>MaxNumEntries_) MaxNumEntries_ = NumEntries;
  }
  
  return(0); // Return 1 because this is an expensive implementation (in runtime)
}
//=============================================================================
int Epetra_BasicRowMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const {

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
int Epetra_BasicRowMatrix::InvRowSums(Epetra_Vector & x) const {
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
int Epetra_BasicRowMatrix::LeftScale(const Epetra_Vector & x) {
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
int Epetra_BasicRowMatrix::InvColSums(Epetra_Vector & x) const {
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
int Epetra_BasicRowMatrix::RightScale(const Epetra_Vector & x) {
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
double Epetra_BasicRowMatrix::NormInf() const {
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
double Epetra_BasicRowMatrix::NormOne() const {
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
int Epetra_BasicRowMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
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
void Epetra_BasicRowMatrix::UpdateImportVector(int NumVectors) const {
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
void Epetra_BasicRowMatrix::UpdateExportVector(int NumVectors) const {
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
void Epetra_BasicRowMatrix::Print(ostream& os) const {

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
