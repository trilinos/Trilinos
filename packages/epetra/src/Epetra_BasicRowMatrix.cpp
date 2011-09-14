
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
Epetra_BasicRowMatrix::Epetra_BasicRowMatrix(const Epetra_Comm & comm) 
  : Comm_(comm.Clone()),
    OperatorDomainMap_(Epetra_Map(0,0,comm)),
    OperatorRangeMap_(Epetra_Map(0,0,comm)),
    RowMatrixRowMap_(Epetra_Map(0,0,comm)),
    RowMatrixColMap_(Epetra_Map(0,0,comm)),
    NumMyNonzeros_(0),
    NumGlobalNonzeros_(0),
    MaxNumEntries_(0),
    NormInf_(0.0),
    NormOne_(0.0),
    UseTranspose_(false),
    HasNormInf_(true),
    LowerTriangular_(true),
    UpperTriangular_(true),
    HaveStructureConstants_(false),
    HaveNumericConstants_(false),
    HaveMaps_(false),
    ImportVector_(0),
    ExportVector_(0),
    Importer_(0),
    Exporter_(0)
{
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
void Epetra_BasicRowMatrix::SetMaps(const Epetra_Map & RowMap, const Epetra_Map & ColMap) {

  SetMaps(RowMap, ColMap, RowMap, RowMap);
}

//==============================================================================
void Epetra_BasicRowMatrix::SetMaps(const Epetra_Map & RowMap, const Epetra_Map & ColMap, 
			const Epetra_Map & DomainMap, const Epetra_Map & RangeMap) {

  RowMatrixRowMap_ = RowMap;
  RowMatrixColMap_ = ColMap;
  OperatorDomainMap_ = DomainMap;
  OperatorRangeMap_ = RangeMap;
  HaveMaps_ = true;
  HaveStructureConstants_ = false;
  HaveNumericConstants_ = false;

  if (!OperatorDomainMap().UniqueGIDs()) throw RowMatrixRowMap().ReportError("At least one GID is repeated in domain map. Domain and range maps must have unique GIDs", -1);
  if (!OperatorRangeMap().UniqueGIDs()) throw RowMatrixRowMap().ReportError("At least one GID is repeated in range map. Domain and range maps must have unique GIDs", -2);
  SetImportExport();
}

//==============================================================================
void Epetra_BasicRowMatrix::SetImportExport() {

  // Check if non-trivial import/export operators
  if (!(RowMatrixRowMap().SameAs(OperatorRangeMap()))) 
    Exporter_ = new Epetra_Export(RowMatrixRowMap(), OperatorRangeMap());
  
  if (!(RowMatrixColMap().SameAs(OperatorDomainMap())))
    Importer_ = new Epetra_Import(RowMatrixColMap(), OperatorDomainMap());

  NumMyRows_ = RowMatrixRowMap().NumMyPoints();
  NumMyCols_ = RowMatrixColMap().NumMyPoints();
}

//==============================================================================
void Epetra_BasicRowMatrix::ComputeStructureConstants() const {
  MaxNumEntries_ = 0;
  NumMyNonzeros_ = 0;
  NumGlobalNonzeros_ = 0;
  int NumEntries = 0;
  for (int i=0; i<NumMyRows_; i++) {
    NumMyRowEntries(i, NumEntries);
    NumMyNonzeros_ += NumEntries;
    if (NumEntries>MaxNumEntries_) MaxNumEntries_ = NumEntries;
  }

  RowMatrixRowMap().Comm().SumAll(&NumMyNonzeros_, &NumGlobalNonzeros_, 1);
  HaveStructureConstants_ = true;
}
//=============================================================================
void Epetra_BasicRowMatrix::ComputeNumericConstants() const {
  Epetra_SerialDenseVector Values(MaxNumEntries());
  Epetra_IntSerialDenseVector Indices(MaxNumEntries());
  int NumEntries;
  Epetra_Vector x1(RowMatrixRowMap()); // Need temp vector for row sums
  Epetra_Vector x2(RowMatrixColMap()); // Need temp vector for column sums
  for(int i = 0; i < NumMyRows_; i++) {
    ExtractMyRowCopy(i, MaxNumEntries(), NumEntries, Values.Values(), Indices.Values());
    for(int j = 0; j < NumEntries; j++) {
      x1[i] += std::abs(Values[j]);
      x2[Indices[j]] += std::abs(Values[j]);
      if (Indices[j]<i) UpperTriangular_ = false;
      if (Indices[j]>i) LowerTriangular_ = false;
    }
  }

  // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
  if(Exporter() != 0) {
    Epetra_Vector xtmp(OperatorRangeMap()); // Create temporary import vector if needed
    xtmp.Export(x1,*Exporter(),Add);
    xtmp.MaxValue(&NormInf_); // This is the NormInf
  }
  else
    x1.MaxValue(&NormInf_); // Find max

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  if(Importer() != 0) {
    Epetra_Vector xtmp(OperatorDomainMap()); // Create temporary import vector if needed
    xtmp.Export(x2,*Importer(),Add);
    xtmp.MaxValue(&NormOne_); // This is the NormOne
  }
  else
    x2.MaxValue(&NormOne_); // Find max

  UpdateFlops(2*NumGlobalNonzeros());
  HaveNumericConstants_ = true;
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
    EPETRA_CHK_ERR(ExtractMyRowCopy(i, MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
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
      EPETRA_CHK_ERR(ExtractMyRowCopy(i, MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
      for (j=0; j < NumEntries; j++)  x_tmp_p[i] += std::abs(Values[j]);
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
      EPETRA_CHK_ERR(ExtractMyRowCopy(i, MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
      double scale = 0.0;
      for (j=0; j < NumEntries; j++) scale += std::abs(Values[j]);
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
  double *curValue;
  int curRowIndex, curColIndex;
  if(OperatorRangeMap().SameAs(x.Map()) && Exporter() != 0) {
    Epetra_Vector xtmp(RowMatrixRowMap());
    xtmp.Import(x,*Exporter(),Insert);
    for (int i=0; i<NumMyNonzeros_; i++) {
      EPETRA_CHK_ERR(ExtractMyEntryView(i, curValue, curRowIndex, curColIndex));
      *curValue *= xtmp[curRowIndex];
    }
  }
  else if (RowMatrixRowMap().SameAs(x.Map()))
    for (int i=0; i<NumMyNonzeros_; i++) {
      EPETRA_CHK_ERR(ExtractMyEntryView(i, curValue, curRowIndex, curColIndex));
      *curValue *= x[curRowIndex];
    }
  else {
    EPETRA_CHK_ERR(-2); // The Map of x must be the RowMap or RangeMap of A.
  }
  HaveNumericConstants_ = false;
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
      EPETRA_CHK_ERR(ExtractMyRowCopy(i, MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
      for(j = 0; j < NumEntries; j++) 
        x_tmp_p[Indices[j]] += std::abs(Values[j]);
    }
    EPETRA_CHK_ERR(x.Export(x_tmp, *Importer(), Add)); // Fill x with partial column sums
  }
  else if(RowMatrixColMap().SameAs(x.Map())) {
    for(i = 0; i < NumMyRows_; i++) {
      EPETRA_CHK_ERR(ExtractMyRowCopy(i, MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
      for(j = 0; j < NumEntries; j++) 
        xp[Indices[j]] += std::abs(Values[j]);
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
  double *curValue;
  int curRowIndex, curColIndex;
  if(OperatorDomainMap().SameAs(x.Map()) && Importer() != 0) {
    Epetra_Vector xtmp(RowMatrixColMap());
    xtmp.Import(x,*Importer(),Insert);
    for (int i=0; i<NumMyNonzeros_; i++) {
      EPETRA_CHK_ERR(ExtractMyEntryView(i, curValue, curRowIndex, curColIndex));
      *curValue *= xtmp[curColIndex];
    }
  }
  else if (RowMatrixColMap().SameAs(x.Map()))
    for (int i=0; i<NumMyNonzeros_; i++) {
      EPETRA_CHK_ERR(ExtractMyEntryView(i, curValue, curRowIndex, curColIndex));
      *curValue *= x[curColIndex];
    }
  else {
    EPETRA_CHK_ERR(-2); // The Map of x must be the RowMap or RangeMap of A.
  }
  HaveNumericConstants_ = false;
  UpdateFlops(NumGlobalNonzeros());
  return(0);
}
//=============================================================================
int Epetra_BasicRowMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  //
  // This function forms the product Y = A * Y or Y = A' * X
  //

  Epetra_SerialDenseVector Values(MaxNumEntries());
  Epetra_IntSerialDenseVector Indices(MaxNumEntries());
  int NumEntries;

  int NumVectors = X.NumVectors();
  if (NumVectors!=Y.NumVectors()) {
    EPETRA_CHK_ERR(-1); // Need same number of vectors in each MV
  }

  UpdateImportVector(NumVectors); // Make sure Import and Export Vectors are compatible
  UpdateExportVector(NumVectors);

  double ** Xp = (double**) X.Pointers();
  double ** Yp = (double**) Y.Pointers();

  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      EPETRA_CHK_ERR(ImportVector_->Import(X, *Importer(), Insert));
      Xp = (double**)ImportVector_->Pointers();
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      Yp = (double**)ExportVector_->Pointers();
    }

    // Do actual computation
    for(int i = 0; i < NumMyRows_; i++) {
      EPETRA_CHK_ERR(ExtractMyRowCopy(i, MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
      for (int k=0; k<NumVectors; k++) {
	double sum = 0.0;
	for(int j = 0; j < NumEntries; j++)
	  sum += Values[j]*Xp[k][Indices[j]];
	Yp[k][i] = sum;
      }
    }
    
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
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      Yp = (double**)ImportVector_->Pointers();
      ImportVector_->PutScalar(0.0);  // Make sure target is zero
    }
    else Y.PutScalar(0.0); // Make sure target is zero

    // Do actual computation
    for(int i = 0; i < NumMyRows_; i++) {
      EPETRA_CHK_ERR(ExtractMyRowCopy(i, MaxNumEntries(), NumEntries, Values.Values(), Indices.Values()));
      for (int k=0; k<NumVectors; k++) {
	double xtmp = Xp[k][i];
	for(int j = 0; j < NumEntries; j++)
	  Yp[k][Indices[j]] += Values[j]*xtmp;
      }
    }
    
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

  int MyPID = RowMatrixRowMap().Comm().MyPID();
  int NumProc = RowMatrixRowMap().Comm().NumProc();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      if (MyPID==0) {
    os <<    "\nNumber of Global Rows         = "; os << NumGlobalRows();    os << endl;
    os <<    "Number of Global Cols         = "; os << NumGlobalCols();    os << endl;
    os <<    "Number of Global Diagonals    = "; os << NumGlobalDiagonals(); os << endl;
	os <<    "Number of Global Nonzeros     = "; os << NumGlobalNonzeros_; os << endl;
      }
      
      os <<  "\nNumber of My Rows               = "; os << NumMyRows_; os << endl;
      os <<    "Number of My Cols               = "; os << NumMyCols_; os << endl;
      os <<    "Number of My Diagonals          = "; os << NumMyDiagonals(); os << endl;
      os <<    "Number of My Nonzeros           = "; os << NumMyNonzeros_; os << endl;
      os <<    "My Maximum Num Entries          = "; os << MaxNumEntries_; os << endl; os << endl;
      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    Comm().Barrier();
    Comm().Barrier();
    Comm().Barrier();
  }
  
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
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
      Epetra_SerialDenseVector Values(MaxNumEntries());
      Epetra_IntSerialDenseVector Indices(MaxNumEntries());
      int NumEntries;
      
      for(int i = 0; i < NumMyRows_; i++) {
	ExtractMyRowCopy(i, MaxNumEntries(), NumEntries, Values.Values(), Indices.Values());
	int Row = RowMatrixRowMap().GID(i);; // Get global row number
	
	for (int j = 0; j < NumEntries ; j++) {   
	  int Index = RowMatrixColMap().GID(Indices[j]);
	  os.width(8);
	  os <<  MyPID ; os << "    ";	
	  os.width(10);
	  os <<  Row ; os << "    ";	
	  os.width(10);
	  os <<  Index; os << "    ";
	  os.width(20);
	  os <<  Values[j]; os << "    ";
	  os << endl;
	}
      }
    
      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    Comm().Barrier();
    Comm().Barrier();
    Comm().Barrier();
  }
	
  return;
}
