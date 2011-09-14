
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

#include "Epetra_ConfigDefs.h"
#include "Epetra_Map.h"

#ifdef HAVE_OSKI
#ifdef HAVE_EPETRA_TEUCHOS
#include "Epetra_OskiMatrix.h"
#include "Epetra_Import.h"

//=============================================================================

Epetra_OskiMatrix::Epetra_OskiMatrix(const Epetra_OskiMatrix& Source) 
  : Epetra_CrsMatrix(Source), 
  Epetra_View_(Source.Epetra_View_),
  Copy_Created_(true) {
    A_tunable_ = oski_CopyMat(Source.A_tunable_);
}

Epetra_OskiMatrix::Epetra_OskiMatrix(const Epetra_CrsMatrix& Source, 
				     const Teuchos::ParameterList& List) 
  : Epetra_CrsMatrix(Source), 
  Epetra_View_(&Source) {
    bool AutoTune = false;
    bool DeepCopy = false;
    string Matrix = "general\0";
    bool IsDiagNotStored = false;
    bool IsArrayZeroBased = false;
    int MyIndexBase = 1;
    bool AreIndicesSorted = false;
    bool AreIndicesRepeated = false;
    oski_inmatprop_t MatrixType = MAT_GENERAL;
    oski_inmatprop_t Diagonal = MAT_DIAG_EXPLICIT;
    oski_inmatprop_t ArrayBasis = INDEX_ONE_BASED;
    oski_inmatprop_t SortedIndices = INDEX_UNSORTED;
    oski_inmatprop_t RepeatedIndices = INDEX_REPEATED;
    int* RowPtr = NULL;
    int* IndPtr = NULL;
    double* ValPtr = NULL;
    AutoTune = const_cast <Teuchos::ParameterList &>(List).get("autotune", false);
    if(List.isParameter("deepcopy")) 
      DeepCopy = Teuchos::getParameter<bool>(List, "deepcopy");
    if(AutoTune){  //Use parameters from the Epetra matrix to set as many fields as possible
      if(LowerTriangular())
        MatrixType = MAT_TRI_LOWER;
      if(UpperTriangular())
        MatrixType = MAT_TRI_UPPER;
      if(Sorted())
        SortedIndices = INDEX_SORTED;
      MyIndexBase = IndexBase();
      if(MyIndexBase == 0)
        ArrayBasis = INDEX_ZERO_BASED;
      else if(MyIndexBase == 1)
        ArrayBasis = INDEX_ONE_BASED;
      else if(!List.isParameter("zerobased")) {
        std::cerr << "An OskiMatrix must be either zero or one based.\n";
        return;
      }
      if(NoRedundancies())
        RepeatedIndices = INDEX_UNIQUE;
    }
    if(List.isParameter("matrixtype")) {
      Matrix = Teuchos::getParameter<string>(List, "matrixtype");
      if(Matrix == "general")
        MatrixType = MAT_GENERAL;
      else if(Matrix == "uppertri")
        MatrixType = MAT_TRI_UPPER;
      else if(Matrix == "lowertri")
        MatrixType = MAT_TRI_LOWER;
      else if(Matrix == "uppersymm")
        MatrixType = MAT_SYMM_UPPER;
      else if(Matrix == "lowersymm")
        MatrixType = MAT_SYMM_LOWER;
      else if(Matrix == "fullsymm")
        MatrixType = MAT_SYMM_FULL;
      else if(Matrix == "upperherm")
        MatrixType = MAT_HERM_UPPER;
      else if(Matrix == "lowerherm")
        MatrixType = MAT_HERM_LOWER;
      else if(Matrix == "fullherm")
        MatrixType = MAT_HERM_FULL;
    }
    if(List.isParameter("diagstored")) 
      IsDiagNotStored = Teuchos::getParameter<bool>(List, "diagstored");
    if(List.isParameter("zerobased")) 
      IsArrayZeroBased = Teuchos::getParameter<bool>(List, "zerobased");
    if(List.isParameter("sorted")) 
      AreIndicesSorted = Teuchos::getParameter<bool>(List, "sorted");
    if(List.isParameter("unique")) 
      DeepCopy = Teuchos::getParameter<bool>(List, "unique");
    if(IsDiagNotStored)
      Diagonal = MAT_UNIT_DIAG_IMPLICIT;
    if(IsArrayZeroBased)
      ArrayBasis = INDEX_ZERO_BASED;
    if(AreIndicesSorted)
      SortedIndices = INDEX_SORTED;
    if(AreIndicesRepeated)
      RepeatedIndices = INDEX_UNIQUE;
    if(ExtractCrsDataPointers(RowPtr, IndPtr, ValPtr)) {
      std::cerr << "Cannot wrap matrix as an Oski Matrix because at least one of FillComplete and Optimize Storage has not been called\n";
      return;
    }
    if(DeepCopy) {  
      Copy_Created_ = true; 
      A_tunable_ = oski_CreateMatCSR(RowPtr, IndPtr, ValPtr, Source.NumMyRows(), Source.NumMyCols(), COPY_INPUTMAT, 5, MatrixType, Diagonal, ArrayBasis, SortedIndices, RepeatedIndices);
    }
    else {
      Copy_Created_ = false;
      A_tunable_ = oski_CreateMatCSR(RowPtr, IndPtr, ValPtr, Source.NumMyRows(), Source.NumMyCols(), SHARE_INPUTMAT, 5, MatrixType, Diagonal, ArrayBasis, SortedIndices, RepeatedIndices);
    }
}

Epetra_OskiMatrix::~Epetra_OskiMatrix() {
  if(oski_DestroyMat(A_tunable_))
    std::cerr << "Destroy Matrix failed.\n";
}

int Epetra_OskiMatrix::ReplaceMyValues(int MyRow, 
				       int NumEntries, 
				       double* Values, 
				       int* Indices) {
  int ReturnVal = 0;
  if (Copy_Created_) {
    for(int i = 0; i < NumEntries; i++) {
      ReturnVal = oski_SetMatEntry(A_tunable_, MyRow, Indices[i], Values[i]);
      if(ReturnVal)
	break;
    }
    if(!ReturnVal)
      ReturnVal = const_cast <Epetra_CrsMatrix*> (Epetra_View_)->ReplaceMyValues(MyRow, NumEntries, Values, Indices);
  }
  else
    ReturnVal = const_cast <Epetra_CrsMatrix*> (Epetra_View_)->ReplaceMyValues(MyRow, NumEntries, Values, Indices);
  if(ReturnVal)
    std::cerr << "Error in ReplaceMyValues\n";
  return ReturnVal;
}

int Epetra_OskiMatrix::SumIntoMyValues(int MyRow, 
				       int NumEntries, 
				       double* Values, 
				       int* Indices) {
  int ReturnVal = 0;
  if (Copy_Created_) {
    for(int i = 0; i < NumEntries; i++) {
      ReturnVal = oski_SetMatEntry(A_tunable_, MyRow, Indices[i], oski_GetMatEntry(A_tunable_, MyRow, Indices[i]));
      if(ReturnVal)
	break;
    }
    if(!ReturnVal)
      ReturnVal = const_cast <Epetra_CrsMatrix*> (Epetra_View_)->SumIntoMyValues(MyRow, NumEntries, Values, Indices);
  }
  else
    ReturnVal = const_cast <Epetra_CrsMatrix*> (Epetra_View_)->SumIntoMyValues(MyRow, NumEntries, Values, Indices);
  if(ReturnVal)
    std::cerr << "Error in SumIntoMyValues\n";
  return ReturnVal;
}

int Epetra_OskiMatrix::ExtractDiagonalCopy(Epetra_Vector& Diagonal) const {
  int ReturnValue = 0;
  ReturnValue = Epetra_View_->ExtractDiagonalCopy(Diagonal);
  if (ReturnValue)
    std::cerr << "Error in ExtractDiagonalCopy\n";
  return ReturnValue;  
}

int Epetra_OskiMatrix::ReplaceDiagonalValues(const Epetra_OskiVector& Diagonal) {
  int ReturnVal = 0;
  if (Copy_Created_) {
    ReturnVal = oski_SetMatDiagValues(A_tunable_, 0, Diagonal.Oski_View());
    if(!ReturnVal)
      ReturnVal = const_cast <Epetra_CrsMatrix*> (Epetra_View_)->ReplaceDiagonalValues(*Diagonal.Epetra_View());
  }
  else
    ReturnVal = const_cast <Epetra_CrsMatrix*> (Epetra_View_)->ReplaceDiagonalValues(*Diagonal.Epetra_View());
  if(ReturnVal)
    std::cerr << "Error in ReplaceDiagonalValues\n";
  return ReturnVal;
}

int Epetra_OskiMatrix::Multiply(bool TransA, 
				const Epetra_Vector& x, 
				Epetra_Vector& y) const {
  int ReturnVal;
  ReturnVal = this->Multiply(TransA, x, y, 1.0, 0.0);
  return ReturnVal;
}

int Epetra_OskiMatrix::Multiply(bool TransA, 
				const Epetra_Vector& x, 
				Epetra_Vector& y, 
				double Alpha, 
				double Beta) const {
  int ReturnVal;
 
  if(!Filled())
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  double* xp = (double*) x.Values();
  double* yp = (double*) y.Values();

  Epetra_Vector * xcopy = 0;
  if (&x==&y && Importer()==0 && Exporter()==0) {
    xcopy = new Epetra_Vector(x);
    xp = (double *) xcopy->Values();
  }
  UpdateImportVector(1); // Refresh import and output vectors if needed
  UpdateExportVector(1);

  if(!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if(Importer() != 0) {
      EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer(), Insert));
      xp = (double*) ImportVector_->Values();
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if(Exporter() != 0) 
      yp = (double*) ExportVector_->Values();
    

    oski_vecview_t oskiX;
    oski_vecview_t oskiY;
    if(Importer() != 0) 
      oskiX = oski_CreateVecView(xp,ImportVector_->MyLength(),1);
    else               
      oskiX = oski_CreateVecView(xp,x.MyLength(),1);
    if(Exporter() != 0) 
      oskiY = oski_CreateVecView(yp,ExportVector_->MyLength(),1);
    else               
      oskiY = oski_CreateVecView(yp,y.MyLength(),1);

    //Do actual computation
    ReturnVal = oski_MatMult(A_tunable_, OP_NORMAL, Alpha, oskiX, Beta, oskiY);

    if(Exporter() != 0) {
      y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(y.Export(*ExportVector_, *Exporter(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());
  } //if(!TransA)
  else { // Transpose operation

    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
    if(Exporter() != 0) {
      EPETRA_CHK_ERR(ExportVector_->Import(x, *Exporter(), Insert));
      xp = (double*) ExportVector_->Values();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if(Importer() != 0) 
      yp = (double*) ImportVector_->Values();

    oski_vecview_t oskiX;
    oski_vecview_t oskiY;
    if(Importer() != 0) 
      oskiY = oski_CreateVecView(yp,ImportVector_->MyLength(),1);
    else               
      oskiY = oski_CreateVecView(yp,y.MyLength(),1);
    if(Exporter() != 0) 
      oskiX = oski_CreateVecView(xp,ExportVector_->MyLength(),1);
    else               
      oskiX = oski_CreateVecView(xp,x.MyLength(),1);
    
    // Do actual computation
    ReturnVal = oski_MatMult(A_tunable_, OP_TRANS, Alpha, oskiX, Beta, oskiY);

    if(Importer() != 0) {
      y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(y.Export(*ImportVector_, *Importer(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());
  }
  if(ReturnVal)
    std::cerr << "OskiVector multiply error\n";
  UpdateFlops(2 * NumGlobalNonzeros());
  return ReturnVal;
}

int Epetra_OskiMatrix::Multiply(bool TransA, 
				const Epetra_MultiVector& X, 
				Epetra_MultiVector& Y) const {
  int ReturnVal;
  ReturnVal = this->Multiply(TransA, X, Y, 1.0, 0.0);
  return ReturnVal;
}

int Epetra_OskiMatrix::Multiply(bool TransA, 
				const Epetra_MultiVector& X, 
				Epetra_MultiVector& Y, 
				double Alpha, 
				double Beta) const {
  int ReturnVal;

  if(!Filled())
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  int NumVectors = X.NumVectors();
  if (NumVectors!=Y.NumVectors()) {
    EPETRA_CHK_ERR(-2); // Need same number of vectors in each MV
  }

  double** Xp = (double**) X.Pointers();
  double** Yp = (double**) Y.Pointers();

  int LDX = X.ConstantStride() ? X.Stride() : 0;
  int LDY = Y.ConstantStride() ? Y.Stride() : 0;

  Epetra_MultiVector* Xcopy = 0;
  if (&X==&Y && Importer()==0 && Exporter()==0) {
    Xcopy = new Epetra_MultiVector(X);
    Xp = (double **) Xcopy->Pointers();
    LDX = Xcopy->ConstantStride() ? Xcopy->Stride() : 0;
  }
  UpdateImportVector(NumVectors); // Refresh import and output vectors if needed
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

    oski_vecview_t oskiX;
    oski_vecview_t oskiY;
    if (Importer()!=0) 
      oskiX = oski_CreateMultiVecView(*Xp,ImportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDX);
    else               
      oskiX = oski_CreateMultiVecView(*Xp,X.MyLength(),NumVectors,LAYOUT_COLMAJ,LDX);
    if (Exporter()!=0) 
      oskiY = oski_CreateMultiVecView(*Yp,ExportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDY);
    else               
      oskiY = oski_CreateMultiVecView(*Yp,Y.MyLength(),NumVectors,LAYOUT_COLMAJ,LDY);
    // Do actual computation
    ReturnVal = oski_MatMult(A_tunable_, OP_NORMAL, Alpha, oskiX, Beta, oskiY);
    if(ReturnVal)
      std::cerr << "OskiMultiVector multiply error\n";
    if (Exporter()!=0) {
      Y.PutScalar(0.0);  // Make sure target is zero
      Y.Export(*ExportVector_, *Exporter(), Add); // Fill Y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Y.Reduce());
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
    
    oski_vecview_t oskiX;
    oski_vecview_t oskiY;
    if (Importer()!=0) 
      oskiY = oski_CreateMultiVecView(*Yp,ImportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDY);
    else               
      oskiY = oski_CreateMultiVecView(*Yp,Y.MyLength(),NumVectors,LAYOUT_COLMAJ,LDY);
    if (Exporter()!=0) 
      oskiX = oski_CreateMultiVecView(*Xp,ExportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDX);
    else               
      oskiX = oski_CreateMultiVecView(*Xp,X.MyLength(),NumVectors,LAYOUT_COLMAJ,LDX);

    // Do actual computation
    ReturnVal = oski_MatMult(A_tunable_, OP_TRANS, Alpha, oskiX, Beta, oskiY);
    if(ReturnVal)
      std::cerr << "OskiMultiVector multiply error\n";
    if (Importer()!=0) {
      Y.PutScalar(0.0);  // Make sure target is zero
      EPETRA_CHK_ERR(Y.Export(*ImportVector_, *Importer(), Add)); // Fill Y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1)  EPETRA_CHK_ERR(Y.Reduce());
  }
  UpdateFlops(2 * NumGlobalNonzeros());
  //Y.ResetView(Yp);
  return ReturnVal;
}

int Epetra_OskiMatrix::Solve(bool Upper, bool TransA, bool UnitDiagonal, const Epetra_Vector& x, Epetra_Vector& y) const {
  int ReturnVal;
  ReturnVal = this->Solve(TransA, x, y, 1.0);
  return ReturnVal;
}

int Epetra_OskiMatrix::Solve(bool TransA, const Epetra_Vector& x, Epetra_Vector& y, double Alpha) const {
  std::cout << "This function Epetra_OskiMatrix::Solve probably works in serial but has not been tested.\n  It will not work in parallel.\n  If you wish to use it feel free to comment out this line and the next return statement.\n  However, correctness and performance are not guaranteed.\n";
  return(-1);
  Epetra_OskiVector* xCast = NULL;
  Epetra_OskiVector* yCast = NULL;
  Epetra_OskiVector* tCast = NULL;
  bool xCreate = false;
  bool yCreate = false;
  int ReturnVal;
  xCast = dynamic_cast<Epetra_OskiVector*>(const_cast<Epetra_Vector*>(&x));
  yCast = dynamic_cast<Epetra_OskiVector*>(&y);
  if (xCast == NULL) {
    xCast = new Epetra_OskiVector(x);
    xCreate = true;
  }
  if (yCast == NULL) {
    yCast = new Epetra_OskiVector(y);
    yCreate = true;
  }
  tCast = new Epetra_OskiVector(x);
  if(TransA)
    ReturnVal = oski_MatTrisolve(A_tunable_, OP_TRANS, Alpha, (*tCast).Oski_View());
  else
    ReturnVal = oski_MatTrisolve(A_tunable_, OP_NORMAL, Alpha, (*tCast).Oski_View());
  if(ReturnVal)
    std::cerr << "OskiVector Solve error\n";
  if(xCreate)
    delete xCast;
  yCast = tCast;
  if(yCreate)
    delete yCast;
  delete tCast;
  return ReturnVal;
}

int Epetra_OskiMatrix::Solve(bool Upper, bool TransA, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  int ReturnVal;
  ReturnVal = this->Solve(TransA, X, Y, 1.0);
  return ReturnVal;
}

int Epetra_OskiMatrix::Solve(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y, double Alpha) const {
  std::cout << "This function Epetra_OskiMatrix::Solve probably works in serial but has not been tested.\n  It will not work in parallel.\n  If you wish to use it feel free to comment out this line and the next return statement.\n  However, correctness and performance are not guaranteed.\n";
  return(-1);
  Epetra_OskiMultiVector* XCast = NULL;
  Epetra_OskiMultiVector* YCast = NULL;
  Epetra_OskiMultiVector* TCast = NULL;
  bool XCreate = false;
  bool YCreate = false;
  int ReturnVal;
  XCast = dynamic_cast<Epetra_OskiMultiVector*>(const_cast<Epetra_MultiVector*>(&X));
  YCast = dynamic_cast<Epetra_OskiMultiVector*>(&Y);
  if (XCast == NULL) {
    XCast = new Epetra_OskiMultiVector(X);
    XCreate = true;
  }
  if (YCast == NULL) {
    YCast = new Epetra_OskiMultiVector(Y);
    YCreate = true;
  }
  TCast = new Epetra_OskiMultiVector(X);
  if(TransA)
    ReturnVal = oski_MatTrisolve(A_tunable_, OP_TRANS, Alpha, (*TCast).Oski_View());
  else
    ReturnVal = oski_MatTrisolve(A_tunable_, OP_NORMAL, Alpha, (*TCast).Oski_View());
  if(ReturnVal)
    std::cerr << "OskiMultiVector Solve error\n";
  if(XCreate)
    delete XCast;
  YCast = TCast;
  if(YCreate)
    delete YCast;
  delete TCast;
  return ReturnVal;
}

int Epetra_OskiMatrix::MatTransMatMultiply(bool ATA, 
					   const Epetra_Vector& x, 
					   Epetra_Vector& y, 
					   Epetra_Vector* t, 
					   double Alpha, 
					   double Beta) const {
  int ReturnVal;
  
  if(!Filled())
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  double* xp = (double*) x.Values();
  double* xp2 = (double*) x.Values();
  double* yp = (double*) y.Values();
  double* tp = 0;
  if(t != NULL)
    tp = (double*) t->Values();

  Epetra_Vector* xcopy = 0;
  if (&x==&y && Importer()==0 && Exporter()==0) {
    xcopy = new Epetra_Vector(x);
    xp = (double *) xcopy->Values();
  }
  UpdateImportVector(1); // Refresh import and output vectors if needed
  UpdateExportVector(1);
  

  if(ATA) {
    if(Importer() != 0) {
      EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer(), Insert));
      xp = (double*) ImportVector_->Values();
      xp2 = new double[ImportVector_->MyLength()];
      for(int i = 0; i < ImportVector_->MyLength(); i++) 
        xp2[i] = xp[i];
      EPETRA_CHK_ERR(ImportVector_->Import(y, *Importer(), Insert));
      yp = (double*) ImportVector_->Values();
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if(Exporter() != 0 && (t) != NULL) {
      tp = (double*) ExportVector_->Values();
    }

    oski_vecview_t oskiX=0;
    oski_vecview_t oskiY=0;
    oski_vecview_t oskiT=0;
    if(Importer() != 0) 
      oskiX = oski_CreateVecView(xp2,ImportVector_->MyLength(),1);
    else
      oskiX = oski_CreateVecView(xp2,x.MyLength(),1);
    if(Importer() != 0) 
      oskiY = oski_CreateVecView(yp,ImportVector_->MyLength(),1);
    else
      oskiY = oski_CreateVecView(yp,y.MyLength(),1);

    if (t != NULL) {
      if(Exporter() != 0) 
          oskiT = oski_CreateVecView(tp,ExportVector_->MyLength(),1);
      else
        oskiT = oski_CreateVecView(tp,t->MyLength(),1);

    }
    else
      oskiT = INVALID_VEC;
    ReturnVal = oski_MatTransMatMult(A_tunable_, OP_AT_A, Alpha, oskiX, Beta, oskiY, oskiT);

    if(Importer() != 0) {
      y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(y.Export(*ImportVector_, *Importer(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());
    
    if(Exporter() != 0 && (t != NULL)) {
      t->PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(t->Export(*ExportVector_, *Exporter(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(t->Reduce());
    UpdateFlops(4 * NumGlobalNonzeros());

  }
  else {
    if(this->Comm().NumProc() == 1) {
      oski_vecview_t oskiX=0;
      oski_vecview_t oskiY=0;
      oski_vecview_t oskiT=0;
      if (t != NULL)
        oskiT = oski_CreateVecView(tp,t->MyLength(),1);
      oskiX = oski_CreateVecView(xp,x.MyLength(),1);
      oskiY = oski_CreateVecView(yp,y.MyLength(),1);
      ReturnVal = oski_MatTransMatMult(A_tunable_, OP_A_AT, Alpha, oskiX, Beta, oskiY, oskiT);
      UpdateFlops(4 * NumGlobalNonzeros());
    }
    else {

      if(t == NULL) {
        Epetra_Vector tempResult(this->DomainMap());
        ReturnVal = this->Multiply(false, x, tempResult, 1.0, 0.0);
        int ReturnVal2 = this->Multiply(true, tempResult, y, Alpha, Beta);
        if(ReturnVal < ReturnVal2)
          ReturnVal = ReturnVal2;
      }
      else {
        ReturnVal = this->Multiply(false, x, *t, 1.0, 0.0);
        int ReturnVal2 = this->Multiply(true,*t, y, Alpha, Beta);
        if(ReturnVal < ReturnVal2)
          ReturnVal = ReturnVal2;
      }
    } 
  }
  if (xcopy!=0) {
    delete xcopy;
    EPETRA_CHK_ERR(1); // Return positive code to alert the user about needing extra copy of x
    return(1);
  }
  return ReturnVal;
}

int Epetra_OskiMatrix::MatTransMatMultiply(bool ATA, 
					   const Epetra_MultiVector& X, 
					   Epetra_MultiVector& Y, 
					   Epetra_MultiVector* T, 
					   double Alpha, 
				           double Beta) const {
  int ReturnVal;
  
  if(!Filled())
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  int NumVectors = X.NumVectors();
  if (NumVectors!=Y.NumVectors()) {
    EPETRA_CHK_ERR(-2); // Need same number of vectors in each MV
  }

  

  double** Xp = (double**) X.Pointers();
  double** Xp2 = (double**) X.Pointers();
  double** Yp = (double**) Y.Pointers();
  double** Tp = 0;
  if(T != NULL)
    Tp = (double**) T->Pointers();
     
  int LDX = X.ConstantStride() ? X.Stride() : 0;
  int LDY = Y.ConstantStride() ? Y.Stride() : 0;
  int LDT = 0;
  if(T != NULL)
    LDT = T->ConstantStride() ? T->Stride() : 0;

  Epetra_MultiVector* Xcopy = 0;
  Epetra_MultiVector* X2 = 0;
  if (&X==&Y && Importer()==0 && Exporter()==0) {
    Xcopy = new Epetra_MultiVector(X);
    Xp = (double **) Xcopy->Pointers();
    LDX = Xcopy->ConstantStride() ? Xcopy->Stride() : 0;
  }
  UpdateImportVector(NumVectors); // Refresh import and output vectors if needed
  UpdateExportVector(NumVectors);
  

  if(ATA) {
    if(Importer() != 0) {
      EPETRA_CHK_ERR(ImportVector_->Import(X, *Importer(), Insert));
      Xp = (double**) ImportVector_->Pointers();
      LDX = ImportVector_->ConstantStride() ? ImportVector_->Stride() : 0;
//      need to think about this
      X2 = new Epetra_MultiVector(X);
      double** Xp2 = (double**) X2->Pointers();
      Xp2 = (double **) X2->Pointers();
      EPETRA_CHK_ERR(ImportVector_->Import(Y, *Importer(), Insert));
      Yp = (double**) ImportVector_->Pointers();
      LDY = ImportVector_->ConstantStride() ? ImportVector_->Stride() : 0;
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if(Exporter() != 0 && T != NULL) {
      Tp = (double**) ExportVector_->Pointers();
      LDT = ImportVector_->ConstantStride() ? ImportVector_->Stride() : 0;
    }
    oski_vecview_t oskiX=0;
    oski_vecview_t oskiY=0;
    oski_vecview_t oskiT=0;
    if(Importer() != 0) 
      oskiX = oski_CreateMultiVecView(*Xp2,ImportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDX);
    else
      oskiX = oski_CreateMultiVecView(*Xp2,X.MyLength(),NumVectors,LAYOUT_COLMAJ,LDX);
    if(Importer() != 0) 
      oskiY = oski_CreateMultiVecView(*Yp,ImportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDY);
    else
      oskiY = oski_CreateMultiVecView(*Yp,Y.MyLength(),NumVectors,LAYOUT_COLMAJ,LDY);

    if (T != NULL) {
      if(Exporter() != 0) 
        oskiT = oski_CreateMultiVecView(*Tp,ExportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDT);
      else
        oskiT = oski_CreateMultiVecView(*Tp,T->MyLength(),NumVectors,LAYOUT_COLMAJ,LDT);

    }
    else
      oskiT = INVALID_VEC;
    ReturnVal = oski_MatTransMatMult(A_tunable_, OP_AT_A, Alpha, oskiX, Beta, oskiY, oskiT);

    if(Importer() != 0) {
      Y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(Y.Export(*ImportVector_, *Importer(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Y.Reduce());
    
    if(Exporter() != 0 && (T != NULL)) {
      T->PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(T->Export(*ExportVector_, *Exporter(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(T->Reduce());
    UpdateFlops(4 * NumGlobalNonzeros());

  }
  else {
    if(this->Comm().NumProc() == 1) {
      oski_vecview_t oskiX=0;
      oski_vecview_t oskiY=0;
      oski_vecview_t oskiT=0;
      if (T != NULL)
        oskiT = oski_CreateMultiVecView(*Tp,T->MyLength(),NumVectors, LAYOUT_COLMAJ,LDT);
      oskiX = oski_CreateMultiVecView(*Xp,X.MyLength(),NumVectors, LAYOUT_COLMAJ,LDX);
      oskiY = oski_CreateMultiVecView(*Yp,Y.MyLength(),NumVectors, LAYOUT_COLMAJ,LDY);
      ReturnVal = oski_MatTransMatMult(A_tunable_, OP_A_AT, Alpha, oskiX, Beta, oskiY, oskiT);
      UpdateFlops(4 * NumGlobalNonzeros() *NumVectors);
    }
    else {
      if(T == NULL) {
        Epetra_MultiVector TempResult(this->DomainMap(), NumVectors);
        ReturnVal = this->Multiply(false, X, TempResult, 1.0, 0.0);
        int ReturnVal2 = this->Multiply(true, TempResult, Y, Alpha, Beta);
        if(ReturnVal < ReturnVal2)
          ReturnVal = ReturnVal2;
      }
      else {
        ReturnVal = this->Multiply(false, X, *T, 1.0, 0.0);
        int ReturnVal2 = this->Multiply(true, *T, Y, Alpha, Beta);
        if(ReturnVal < ReturnVal2)
          ReturnVal = ReturnVal2;
      }
    } 
  }
  if (Xcopy!=0) {
    delete Xcopy;
    EPETRA_CHK_ERR(1); // Return positive code to alert the user about needing extra copy of x
    return(1);
  }
  return ReturnVal;
}

int Epetra_OskiMatrix::MultiplyAndMatTransMultiply(bool TransA, 
						   const Epetra_Vector& x, 
						   Epetra_Vector& y, 
						   const Epetra_Vector& w, 
						   Epetra_Vector& z, 
						   double Alpha, 
						   double Beta, 
						   double Omega, 
						   double Zeta) const {

  int ReturnVal;

  if(!Filled())
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  double* xp = (double*) x.Values();
  double* xp2 = (double*) x.Values();
  double* wp = (double*) w.Values();
  double* yp = (double*) y.Values();
//  double* yp2 = (double*) y.Values();
  double* zp = (double*) z.Values();
  Epetra_MultiVector* yp2 = 0;
  Epetra_Vector* xcopy = 0;
  if (&x==&y && Importer()==0 && Exporter()==0) {
    xcopy = new Epetra_Vector(x);
    xp = (double *) xcopy->Values();
  }
  Epetra_Vector* wcopy = 0;
  if (&w==&z && Importer()==0 && Exporter()==0) {
    wcopy = new Epetra_Vector(w);
    wp = (double *) wcopy->Values();
  }
  UpdateImportVector(1); // Refresh import and output vectors if needed
  UpdateExportVector(1);

  if(TransA) {

    if(Importer() != 0) {
      EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer(), Insert));
      xp = (double*) ImportVector_->Values();
      xp2 = new double[ImportVector_->MyLength()];
      for(int i = 0; i < ImportVector_->MyLength(); i++)
        xp2[i] = xp[i];
      EPETRA_CHK_ERR(ImportVector_->Import(z, *Importer(), Insert));
      zp = (double*) ImportVector_->Values();
    }
    if(Exporter() != 0) {
      yp2 = new Epetra_MultiVector(*ExportVector_);
      yp = (double*) yp2->Values();
      
      //for(int i = 0; i < ExportVector_->MyLength(); i++)
        //yp2[i] = yp[i];
      wp = (double*) ExportVector_->Values();
    }

    oski_vecview_t oskiX=0;
    oski_vecview_t oskiY=0;
    oski_vecview_t oskiW=0;
    oski_vecview_t oskiZ=0;
    if(Importer() != 0) {
      oskiX = oski_CreateVecView(xp2,ImportVector_->MyLength(),1);
      oskiZ = oski_CreateVecView(zp,ImportVector_->MyLength(),1);
    }
    else {
      oskiX = oski_CreateVecView(xp2,x.MyLength(),1);
      oskiZ = oski_CreateVecView(zp,z.MyLength(),1);
    }
    if(Exporter() != 0) {
      oskiY = oski_CreateVecView(yp,ExportVector_->MyLength(),1);
      oskiW = oski_CreateVecView(wp,ExportVector_->MyLength(),1);
    }
    else {
      oskiY = oski_CreateVecView(yp,y.MyLength(),1);
      oskiW = oski_CreateVecView(wp,w.MyLength(),1);
    }
   
    ReturnVal = oski_MatMultAndMatTransMult(A_tunable_, Alpha, oskiX, Beta, oskiY, OP_TRANS, Omega, oskiW, Zeta, oskiZ);
    if(Exporter() != 0) {
      y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(y.Export(*yp2, *Exporter(), Add)); // Fill y with Values from export vector
    }
    if(Importer() != 0) {
      z.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(z.Export(*ImportVector_, *Importer(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(z.Reduce());

    UpdateFlops(4 * NumGlobalNonzeros());
  }
  //  ReturnVal = oski_MatMultAndMatTransMult(A_tunable_, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), OP_TRANS, Omega, (*wCast).Oski_View(), Zeta, (*zCast).Oski_View());
  else {

   if(Importer() != 0) {
      EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer(), Insert));
      xp = (double*) ImportVector_->Values();
      xp2 = new double[ImportVector_->MyLength()];
      for(int i = 0; i < ImportVector_->MyLength(); i++)
        xp2[i] = xp[i];
      EPETRA_CHK_ERR(ImportVector_->Import(w, *Importer(), Insert));
      wp = (double*) ImportVector_->Values();
    }
    if(Exporter() != 0) {
      yp2 = new Epetra_MultiVector(*ExportVector_);
      yp = (double*) yp2->Values();
      
      //for(int i = 0; i < ExportVector_->MyLength(); i++)
        //yp2[i] = yp[i];
      zp = (double*) ExportVector_->Values();
    }

    oski_vecview_t oskiX=0;
    oski_vecview_t oskiY=0;
    oski_vecview_t oskiW=0;
    oski_vecview_t oskiZ=0;
    if(Importer() != 0) {
      oskiX = oski_CreateVecView(xp2,ImportVector_->MyLength(),1);
      oskiW = oski_CreateVecView(wp,ImportVector_->MyLength(),1);
    }
    else {
      oskiX = oski_CreateVecView(xp2,x.MyLength(),1);
      oskiW = oski_CreateVecView(wp,w.MyLength(),1);
    }
    if(Exporter() != 0) {
      oskiY = oski_CreateVecView(yp,ExportVector_->MyLength(),1);
      oskiZ = oski_CreateVecView(zp,ExportVector_->MyLength(),1);
    }
    else {
      oskiY = oski_CreateVecView(yp,y.MyLength(),1);
      oskiZ = oski_CreateVecView(zp,z.MyLength(),1);
    }
   
    ReturnVal = oski_MatMultAndMatTransMult(A_tunable_, Alpha, oskiX, Beta, oskiY, OP_NORMAL, Omega, oskiW, Zeta, oskiZ);
    if(Exporter() != 0) {
      y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(y.Export(*yp2, *Exporter(), Add)); // Fill y with Values from export vector
      z.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(z.Export(*ExportVector_, *Exporter(), Add)); // Fill y with Values from export vector*/
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(z.Reduce());

    UpdateFlops(4 * NumGlobalNonzeros());

  }
  if(ReturnVal)
    std::cerr << "OskiVector multiply error\n";
  return ReturnVal;
}

int Epetra_OskiMatrix::MultiplyAndMatTransMultiply(bool TransA, 
						   const Epetra_MultiVector& X, 
						   Epetra_MultiVector& Y, 
						   const Epetra_MultiVector& W, 
						   Epetra_MultiVector& Z, 
						   double Alpha, 
						   double Beta, 
						   double Omega, 
						   double Zeta) const {
  int ReturnVal;
  if(!Filled())
    EPETRA_CHK_ERR(-1); // Matrix must be filled.
  int NumVectors = X.NumVectors();
  if (NumVectors!=Y.NumVectors()) {
    EPETRA_CHK_ERR(-2); // Need same number of vectors in each MV
  }


  double** Xp = (double**) X.Pointers();
  double** Xp2 = (double**) X.Pointers();
  double** Wp = (double**) W.Pointers();
  double** Yp = (double**) Y.Pointers();
  double** Zp = (double**) Z.Pointers();
  int LDX = X.ConstantStride() ? X.Stride() : 0;
  int LDY = Y.ConstantStride() ? Y.Stride() : 0;
  int LDW = W.ConstantStride() ? W.Stride() : 0;
  int LDZ = Z.ConstantStride() ? Z.Stride() : 0;

  Epetra_MultiVector* Yp2 = 0;
  Epetra_MultiVector* X2 = 0;
  Epetra_MultiVector* Xcopy = 0;
  if (&X==&Y && Importer()==0 && Exporter()==0) {
    Xcopy = new Epetra_MultiVector(X);
    Xp = (double **) Xcopy->Pointers();
    LDX = Xcopy->ConstantStride() ? Xcopy->Stride() : 0;
  }
  Epetra_MultiVector* Wcopy = 0;
  if (&W==&Z && Importer()==0 && Exporter()==0) {
    Wcopy = new Epetra_MultiVector(W);
    Wp = (double **) Wcopy->Values();
    LDW = Wcopy->ConstantStride() ? Wcopy->Stride() : 0;
  }
  UpdateImportVector(NumVectors); // Refresh import and output vectors if needed
  UpdateExportVector(NumVectors);

  if(TransA) {

    if(Importer() != 0) {
      EPETRA_CHK_ERR(ImportVector_->Import(X, *Importer(), Insert));
      Xp = (double**) ImportVector_->Pointers();
      LDX = ImportVector_->ConstantStride() ? ImportVector_->Stride() : 0;
      X2 = new Epetra_MultiVector(X);
      double** Xp2 = (double**) X2->Pointers();
      Xp2 = (double **) X2->Pointers();
      EPETRA_CHK_ERR(ImportVector_->Import(Y, *Importer(), Insert));
      Zp = (double**) ImportVector_->Pointers();
      LDZ = ImportVector_->ConstantStride() ? ImportVector_->Stride() : 0;
    }
    if(Exporter() != 0) {
      Yp2 = new Epetra_MultiVector(*ExportVector_);
      Yp = (double**) Yp2->Pointers();
      Wp = (double**) ExportVector_->Pointers();
    }

    oski_vecview_t oskiX=0;
    oski_vecview_t oskiY=0;
    oski_vecview_t oskiW=0;
    oski_vecview_t oskiZ=0;
    if(Importer() != 0) {
      oskiX = oski_CreateMultiVecView(*Xp2,ImportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDX);
      oskiZ = oski_CreateMultiVecView(*Zp,ImportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDZ);
    }
    else {
      oskiX = oski_CreateMultiVecView(*Xp2,X.MyLength(),NumVectors,LAYOUT_COLMAJ,LDX);
      oskiZ = oski_CreateMultiVecView(*Zp,Z.MyLength(),NumVectors,LAYOUT_COLMAJ,LDZ);
    }
    if(Exporter() != 0) {
      oskiY = oski_CreateMultiVecView(*Yp,ExportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDY);
      oskiW = oski_CreateMultiVecView(*Wp,ExportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDW);
    }
    else {
      oskiY = oski_CreateMultiVecView(*Yp,Y.MyLength(),NumVectors,LAYOUT_COLMAJ,LDY);
      oskiW = oski_CreateMultiVecView(*Wp,W.MyLength(),NumVectors,LAYOUT_COLMAJ,LDW);
    }
   
    ReturnVal = oski_MatMultAndMatTransMult(A_tunable_, Alpha, oskiX, Beta, oskiY, OP_TRANS, Omega, oskiW, Zeta, oskiZ);
    if(Exporter() != 0) {
      Y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(Y.Export(*Yp2, *Exporter(), Add)); // Fill y with Values from export vector
    }
    if(Importer() != 0) {
      Z.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(Z.Export(*ImportVector_, *Importer(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Y.Reduce());
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Z.Reduce());

    UpdateFlops(4 * NumGlobalNonzeros());
  }
  //  ReturnVal = oski_MatMultAndMatTransMult(A_tunable_, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), OP_TRANS, Omega, (*wCast).Oski_View(), Zeta, (*zCast).Oski_View());
  else {

   if(Importer() != 0) {
      EPETRA_CHK_ERR(ImportVector_->Import(X, *Importer(), Insert));
      Xp = (double**) ImportVector_->Pointers();
      LDX = ImportVector_->ConstantStride() ? ImportVector_->Stride() : 0;
      X2 = new Epetra_MultiVector(X);
      Xp2 = (double**) X2->Pointers();
      EPETRA_CHK_ERR(ImportVector_->Import(W, *Importer(), Insert));
      Wp = (double**) ImportVector_->Pointers();
      LDW = ImportVector_->ConstantStride() ? ImportVector_->Stride() : 0;
    }
    if(Exporter() != 0) {
      Yp2 = new Epetra_MultiVector(*ExportVector_);
      Yp = (double**) Yp2->Pointers();
      Zp = (double**) ExportVector_->Pointers();
    }

    oski_vecview_t oskiX=0;
    oski_vecview_t oskiY=0;
    oski_vecview_t oskiW=0;
    oski_vecview_t oskiZ=0;
    if(Importer() != 0) {
      oskiX = oski_CreateMultiVecView(*Xp2,ImportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDX);
      oskiW = oski_CreateMultiVecView(*Wp,ImportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDW);
    }
    else {
      oskiX = oski_CreateMultiVecView(*Xp2,X.MyLength(),NumVectors,LAYOUT_COLMAJ,LDX);
      oskiW = oski_CreateMultiVecView(*Wp,W.MyLength(),NumVectors,LAYOUT_COLMAJ,LDW);
    }
    if(Exporter() != 0) {
      oskiY = oski_CreateMultiVecView(*Yp,ExportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDY);
      oskiZ = oski_CreateMultiVecView(*Zp,ExportVector_->MyLength(),NumVectors,LAYOUT_COLMAJ,LDZ);
    }
    else {
      oskiY = oski_CreateMultiVecView(*Yp,Y.MyLength(),NumVectors,LAYOUT_COLMAJ,LDY);
      oskiZ = oski_CreateMultiVecView(*Zp,Z.MyLength(),NumVectors,LAYOUT_COLMAJ,LDZ);
    }
   
    ReturnVal = oski_MatMultAndMatTransMult(A_tunable_, Alpha, oskiX, Beta, oskiY, OP_NORMAL, Omega, oskiW, Zeta, oskiZ);
    if(Exporter() != 0) {
      Y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(Y.Export(*Yp2, *Exporter(), Add)); // Fill y with Values from export vector
      Z.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(Z.Export(*ExportVector_, *Exporter(), Add)); // Fill y with Values from export vector*/
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Y.Reduce());
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Z.Reduce());

    UpdateFlops(4 * NumGlobalNonzeros());

  }
  return ReturnVal;
}

int Epetra_OskiMatrix::MatPowMultiply(bool TransA, 
				      const Epetra_Vector& x, 
				      Epetra_Vector& y, 
 				      Epetra_MultiVector& T, 
				      int Power,
				      double Alpha, 
				      double Beta) const {
  //The code has not been tested.  It should work in serial but not in parallel.
  std::cerr << "MatPowMultiply not implemented in oski-1.01h release.\n";
  return -1;

  int ReturnVal; 

  if(!Filled())
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  double* xp = (double*) x.Values();
  double* yp = (double*) y.Values();
  double** Tp = (double**) T.Pointers();

  Epetra_MultiVector *Tptr;

  int LDT = T.ConstantStride() ? T.Stride() : 0;


  if(this->Comm().NumProc() == 1) {
    oski_vecview_t oskiX=0;
    oski_vecview_t oskiY=0;
    oski_vecview_t oskiT=0;
    if (&T != NULL)
      oskiT = oski_CreateMultiVecView(*Tp,T.MyLength(),Power,LAYOUT_COLMAJ,LDT);
    oskiX = oski_CreateVecView(xp,x.MyLength(),1);
    oskiY = oski_CreateVecView(yp,y.MyLength(),1);

    if(TransA)
      ReturnVal = oski_MatPowMult(A_tunable_, OP_TRANS, Power, Alpha, oskiX, Beta, oskiY, oskiT);
    else
      ReturnVal = oski_MatPowMult(A_tunable_, OP_NORMAL, Power, Alpha, oskiX, Beta, oskiY, oskiT);
    if(ReturnVal)
      std::cerr << "OskiVector matpow multiply error\n";
  }
  else {
    if(&T == NULL) 
      Tptr = new Epetra_MultiVector(x.Map(), Power);
    else
      Tptr = &T;
    if(TransA) {
      ReturnVal = this->Multiply(true, x, Tptr[1], 1.0, 0.0);
      for(int i = 1; i < Power-1; i++)
        ReturnVal = this->Multiply(true, Tptr[i], Tptr[i+1], 1.0, 0.0);
      ReturnVal = this->Multiply(true, Tptr[Power-2], y, Alpha, Beta);
    }
    else {
      ReturnVal = this->Multiply(false, x, Tptr[1], 1.0, 0.0);
      for(int i = 1; i < Power-1; i++)
        ReturnVal = this->Multiply(false, Tptr[i], Tptr[i+1], 1.0, 0.0);
      ReturnVal = this->Multiply(false, Tptr[Power-2], y, Alpha, Beta);
    }
    if(ReturnVal)
      std::cerr << "OskiVector matpow multiply error\n";
    if(&T == NULL)
      delete(Tptr);
  }    
  UpdateFlops(2 * Power * NumGlobalNonzeros());
  return ReturnVal;
}

int Epetra_OskiMatrix::MatPowMultiply(bool TransA, 
				      const Epetra_Vector& x, 
				      Epetra_Vector& y, 
				      int Power,
				      double Alpha, 
				      double Beta) const {

  //The code has not been tested.  It should work in serial but not in parallel.
  std::cerr << "MatPowMultiply not implemented in oski-1.01h release.\n";
  return -1;
  std::cerr << "mult\n";
  Epetra_OskiVector* xCast = NULL;
  Epetra_OskiVector* yCast = NULL;
  bool xCreate = false;
  bool yCreate = false;
  int ReturnVal;
  xCast = dynamic_cast<Epetra_OskiVector*>(const_cast <Epetra_Vector*>(&x));
  yCast = dynamic_cast<Epetra_OskiVector*>(&y);
  if (xCast == NULL) {
    xCast = new Epetra_OskiVector(x);
    xCreate = true;
  }
  if (yCast == NULL) {
    yCast = new Epetra_OskiVector(y);
    yCreate = true;
  }
  std::cerr << "mult\n";
  if(TransA)
    ReturnVal = oski_MatPowMult(A_tunable_, OP_TRANS, Power, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), NULL);
  else
    ReturnVal = oski_MatPowMult(A_tunable_, OP_NORMAL, Power, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), NULL);
  std::cerr << "done\n";
  if(ReturnVal)
    std::cerr << "OskiVector matpow multiply error\n";
  if(xCreate)
    delete xCast;
  if(yCreate)
    delete yCast;
  return ReturnVal;
}

int Epetra_OskiMatrix::SetHint(const Teuchos::ParameterList& List) {
  int* ArgArray = NULL;
  int Diags;
  int Blocks;
  char Number[10];
  char Row[20];
  char Col[20];
  char Diag[20];
  int ReturnVal = 0;
  if(List.isParameter("randompattern"))
    if(Teuchos::getParameter<bool>(List, "randompattern"))
      if(ReturnVal = oski_SetHint(A_tunable_, HINT_RANDOM_PATTERN))
        std::cerr << "Error setting hint random pattern.\n";
  if(List.isParameter("correlatedpattern"))
    if(Teuchos::getParameter<bool>(List, "correlatedpattern"))
      if(ReturnVal = oski_SetHint(A_tunable_, HINT_CORREL_PATTERN))
        std::cerr << "Error setting hint correlated pattern.\n";
  if(List.isParameter("symmetricpattern"))
    if(Teuchos::getParameter<bool>(List, "symmetricpattern"))
      if(ReturnVal = oski_SetHint(A_tunable_, HINT_SYMM_PATTERN))
        std::cerr << "Error setting hint symmetric pattern.\n";
  if(List.isParameter("nonsymmetricpattern"))
    if(Teuchos::getParameter<bool>(List, "nonsymmetricpattern"))
      if(ReturnVal = oski_SetHint(A_tunable_, HINT_NONSYMM_PATTERN))
        std::cerr << "Error setting hint nonsymmetric pattern.\n";
  if(List.isParameter("alignedblocks"))
    if(Teuchos::getParameter<bool>(List, "alignedblocks"))
    {
      if(ReturnVal = oski_SetHint(A_tunable_, HINT_ALIGNED_BLOCKS))
        std::cerr << "Error setting hint aligned blocks.\n";
    }
  if(List.isParameter("unalignedblocks"))
    if(Teuchos::getParameter<bool>(List, "unalignedblocks"))
      if(ReturnVal = oski_SetHint(A_tunable_, HINT_UNALIGNED_BLOCKS))
        std::cerr << "Error setting hint unaligned blocks.\n";
  if(List.isParameter("nodiags"))
    if(Teuchos::getParameter<bool>(List, "nodiags"))
      if(ReturnVal = oski_SetHint(A_tunable_, HINT_NO_DIAGS))
        std::cerr << "Error setting hint no diags.\n";
  if(List.isParameter("noblocks"))
    if(Teuchos::getParameter<bool>(List, "noblocks"))
      if(ReturnVal = oski_SetHint(A_tunable_, HINT_NO_BLOCKS))
        std::cerr << "Error setting hint no blocks.\n";
  if(List.isParameter("singleblocksize"))
    if(Teuchos::getParameter<bool>(List, "singleblocksize")) {
      ArgArray = new int[2];
      if(List.isParameter("row"))
      {
        ArgArray[0] = Teuchos::getParameter<int>(List, "row");
        if(List.isParameter("col"))
          ArgArray[1] = Teuchos::getParameter<int>(List, "col");
        if(ReturnVal = oski_SetHint(A_tunable_, HINT_SINGLE_BLOCKSIZE, ArgArray[0], ArgArray[1]))
          std::cerr << "Error setting hint single block size.\n";
      }
      else
        if(ReturnVal = oski_SetHint(A_tunable_, HINT_SINGLE_BLOCKSIZE))
          std::cerr << "Error setting hint single block size.\n";
      delete [] ArgArray;
      ArgArray = NULL;
    }
  if(List.isParameter("multipleblocksize"))
    if(Teuchos::getParameter<bool>(List, "multipleblocksize"))
      if(List.isParameter("blocks")) {
        Blocks = Teuchos::getParameter<int>(List, "blocks");
        ArgArray = new int[2*Blocks+1];
	ArgArray[0] = Blocks;
        for(int i = 0; i < Blocks; i++) {
	  sprintf(Number, "%d", i+1);
	  strcpy(Row, "row");
	  strcpy(Col, "col");
	  strcat(Row, Number);
	  strcat(Col, Number);
          if(List.isParameter(Row))
            ArgArray[i*2 + 1] = Teuchos::getParameter<int>(List, Row);
          if(List.isParameter(Col))
            ArgArray[i*2 + 2] = Teuchos::getParameter<int>(List, Col);
        }
        switch(Blocks) {
          case 1 :  if(ReturnVal = oski_SetHint(A_tunable_, HINT_MULTIPLE_BLOCKSIZES, ArgArray[0], ArgArray[1], ArgArray[2]))
          	      std::cerr << "Error setting hint multiple blocks.\n";
		    break;
          case 2 :  if(ReturnVal = oski_SetHint(A_tunable_, HINT_MULTIPLE_BLOCKSIZES, ArgArray[0], ArgArray[1], ArgArray[2], ArgArray[3], ArgArray[4]))
          	      std::cerr << "Error setting hint multiple blocks.\n";
		    break;
          case 3 :  if(ReturnVal = oski_SetHint(A_tunable_, HINT_MULTIPLE_BLOCKSIZES, ArgArray[0], ArgArray[1], ArgArray[2], ArgArray[3], ArgArray[4], ArgArray[5], ArgArray[6]))
          	      std::cerr << "Error setting hint multiple blocks.\n";
		    break;
          case 4 :  if(ReturnVal = oski_SetHint(A_tunable_, HINT_MULTIPLE_BLOCKSIZES, ArgArray[0], ArgArray[1], ArgArray[2], ArgArray[3], ArgArray[4], ArgArray[5], ArgArray[6], ArgArray[7], ArgArray[8]))
          	      std::cerr << "Error setting hint multiple blocks.\n";
		    break;
          case 5 :  if(ReturnVal = oski_SetHint(A_tunable_, HINT_MULTIPLE_BLOCKSIZES, ArgArray[0], ArgArray[1], ArgArray[2], ArgArray[3], ArgArray[4], ArgArray[5], ArgArray[6], ArgArray[7], ArgArray[8], ArgArray[9], ArgArray[10]))
          	      std::cerr << "Error setting hint multiple blocks.\n";
		    break;
          default : if(ReturnVal = oski_SetHint(A_tunable_, HINT_MULTIPLE_BLOCKSIZES))
                      std::cerr << "Error setting hint multiple blocks.\n";
                    break;
        }
        delete [] ArgArray;
        ArgArray = NULL;
      }
      else
        if(ReturnVal = oski_SetHint(A_tunable_, HINT_MULTIPLE_BLOCKSIZES))
          std::cerr << "Error setting hint multiple blocks.\n";
  if(List.isParameter("diags"))
    if(Teuchos::getParameter<bool>(List, "diags"))
      if(List.isParameter("numdiags")) {
        Diags = Teuchos::getParameter<int>(List, "numdiags");
        ArgArray = new int[Diags+1];
	ArgArray[0] = Diags;
        for(int i = 0; i < Diags; i++) {
	  sprintf(Number, "%d", i + 1);
	  strcpy(Diag, "diag");
	  strcat(Diag, Number);
          if(List.isParameter(Diag))
            ArgArray[i + 1] = Teuchos::getParameter<int>(List, Diag);
        }
        switch(Diags) {
          case 1 : if(ReturnVal = oski_SetHint(A_tunable_, HINT_DIAGS, ArgArray[0], ArgArray[1]))
                     std::cerr << "Error setting hint diags\n";
                   break;
          case 2 : if(ReturnVal = oski_SetHint(A_tunable_, HINT_DIAGS, ArgArray[0], ArgArray[1], ArgArray[2]))
                     std::cerr << "Error setting hint diags\n";
                   break;
          case 3 : if(ReturnVal = oski_SetHint(A_tunable_, HINT_DIAGS, ArgArray[0], ArgArray[1], ArgArray[2], ArgArray[3]))
                     std::cerr << "Error setting hint diags\n";
                   break;
          case 4 : if(ReturnVal = oski_SetHint(A_tunable_, HINT_DIAGS, ArgArray[0], ArgArray[1], ArgArray[2], ArgArray[3], ArgArray[4]))
                     std::cerr << "Error setting hint diags\n";
                   break;
          case 5 : if(ReturnVal = oski_SetHint(A_tunable_, HINT_DIAGS, ArgArray[0], ArgArray[1], ArgArray[2], ArgArray[3], ArgArray[4], ArgArray[5]))
                     std::cerr << "Error setting hint diags\n";
                   break;
          default : if(ReturnVal = oski_SetHint(A_tunable_, HINT_DIAGS, ArgArray[0]))
                      std::cerr << "Error setting hint diags\n";
                    break;
        }
        delete [] ArgArray;
      }
      else
      {
        if(ReturnVal = oski_SetHint(A_tunable_, HINT_DIAGS))
          std::cerr << "Error setting hint digs.\n";
      }
  return ReturnVal;
}

int Epetra_OskiMatrix::SetHintMultiply(bool TransA, 
				       double Alpha, 
				       const Epetra_OskiMultiVector& InVec, 
				       double Beta, 
				       const Epetra_OskiMultiVector& OutVec, 
				       int NumCalls, 
				       const Teuchos::ParameterList& List) {
  int ReturnVal;
  oski_vecview_t InView = NULL;
  oski_vecview_t OutView = NULL;
  InView = InVec.Oski_View();
  OutView = OutVec.Oski_View();
  if(List.isParameter("tune"))
    if(Teuchos::getParameter<bool>(List, "tune"))
      NumCalls = ALWAYS_TUNE;
  if(List.isParameter("tuneaggressive"))
    if(Teuchos::getParameter<bool>(List, "tuneaggressive"))
      NumCalls = ALWAYS_TUNE_AGGRESSIVELY;
  if(List.isParameter("symminvec"))
    if(Teuchos::getParameter<bool>(List, "symminvec"))
      InView = SYMBOLIC_VEC;
  if(List.isParameter("symminmultivec"))
    if(Teuchos::getParameter<bool>(List, "symminmultivec"))
      InView = SYMBOLIC_MULTIVEC;
  if(List.isParameter("symmoutvec"))
    if(Teuchos::getParameter<bool>(List, "symmoutvec"))
      OutView = SYMBOLIC_VEC;
  if(List.isParameter("symmoutmultivec"))
    if(Teuchos::getParameter<bool>(List, "symmoutmultivec"))
      OutView = SYMBOLIC_MULTIVEC;
  if(this->Comm().NumProc() > 1) {
    if(InVec.NumVectors() == 1)
      InView = SYMBOLIC_VEC;
    else
      InView = SYMBOLIC_MULTIVEC;
    if(OutVec.NumVectors() == 1)
      OutView = SYMBOLIC_VEC;
    else
      OutView = SYMBOLIC_MULTIVEC;
  }
  if(TransA)
    ReturnVal = oski_SetHintMatMult(A_tunable_, OP_TRANS, Alpha, InView, Beta, OutView, NumCalls);
  else
    ReturnVal = oski_SetHintMatMult(A_tunable_, OP_NORMAL, Alpha, InView, Beta, OutView, NumCalls);
  if(ReturnVal)
    std::cerr << "Set hint multivector multiply error\n";
  return ReturnVal;
}

int Epetra_OskiMatrix::SetHintSolve(bool TransA,
				    double Alpha, 
				    const Epetra_OskiMultiVector& Vector, 
				    int NumCalls, 
				    const Teuchos::ParameterList& List) {
  int ReturnVal;
  oski_vecview_t VecView = NULL;
  VecView = Vector.Oski_View();
  if(List.isParameter("tune"))
    if(Teuchos::getParameter<bool>(List, "tune"))
      NumCalls = ALWAYS_TUNE;
  if(List.isParameter("tuneaggressive"))
    if(Teuchos::getParameter<bool>(List, "tuneaggressive"))
      NumCalls = ALWAYS_TUNE_AGGRESSIVELY;
  if(List.isParameter("symmvec"))
    if(Teuchos::getParameter<bool>(List, "symmvec"))
      VecView = SYMBOLIC_VEC;
  if(List.isParameter("symmmultivec"))
    if(Teuchos::getParameter<bool>(List, "symmmultivec"))
      VecView = SYMBOLIC_MULTIVEC;
  if(this->Comm().NumProc() > 1) {
    if(Vector.NumVectors() == 1)
      VecView = SYMBOLIC_VEC;
    else
      VecView = SYMBOLIC_MULTIVEC;
  }
  if(TransA)
    ReturnVal = oski_SetHintMatTrisolve(A_tunable_, OP_TRANS, Alpha, VecView, NumCalls);
  else
    ReturnVal = oski_SetHintMatTrisolve(A_tunable_, OP_NORMAL, Alpha, VecView, NumCalls);
  if(ReturnVal)
    std::cerr << "Set hint solve error\n";
  return ReturnVal;
}

int Epetra_OskiMatrix::SetHintMatTransMatMultiply (bool ATA, 
						   double Alpha, 
						   const Epetra_OskiMultiVector& InVec, 
						   double Beta, 
						   const Epetra_OskiMultiVector& OutVec, 
						   const Epetra_OskiMultiVector& Intermediate, 
						   int NumCalls, 
						   const Teuchos::ParameterList& List) {
  int ReturnVal;
  oski_vecview_t InView = NULL;
  oski_vecview_t OutView = NULL;
  oski_vecview_t IntermediateView = NULL;
  InView = InVec.Oski_View();
  OutView = OutVec.Oski_View();
  if(&Intermediate != NULL)
    IntermediateView = Intermediate.Oski_View();
  if(List.isParameter("tune"))
    if(Teuchos::getParameter<bool>(List, "tune"))
      NumCalls = ALWAYS_TUNE;
  if(List.isParameter("tuneaggressive"))
    if(Teuchos::getParameter<bool>(List, "tuneaggressive"))
      NumCalls = ALWAYS_TUNE_AGGRESSIVELY;
  if(List.isParameter("symminvec"))
    if(Teuchos::getParameter<bool>(List, "symminvec"))
      InView = SYMBOLIC_VEC;
  if(List.isParameter("symminmultivec"))
    if(Teuchos::getParameter<bool>(List, "symminmultivec"))
      InView = SYMBOLIC_MULTIVEC;
  if(List.isParameter("symmoutvec"))
    if(Teuchos::getParameter<bool>(List, "symmoutvec"))
      OutView = SYMBOLIC_VEC;
  if(List.isParameter("symmoutmultivec"))
    if(Teuchos::getParameter<bool>(List, "symmoutmultivec"))
      OutView = SYMBOLIC_MULTIVEC;
  if(List.isParameter("symmintervec"))
    if(Teuchos::getParameter<bool>(List, "symmintervec"))
      IntermediateView = SYMBOLIC_VEC;
  if(List.isParameter("symmintermultivec"))
    if(Teuchos::getParameter<bool>(List, "symmintermultivec"))
      IntermediateView = SYMBOLIC_MULTIVEC;
  if(List.isParameter("invalidinter"))
    if(Teuchos::getParameter<bool>(List, "invalidinter"))
      IntermediateView = NULL;
  if(this->Comm().NumProc() > 1) {
    if(InVec.NumVectors() == 1)
      InView = SYMBOLIC_VEC;
    else
      InView = SYMBOLIC_MULTIVEC;
    if(OutVec.NumVectors() == 1)
      OutView = SYMBOLIC_VEC;
    else
      OutView = SYMBOLIC_MULTIVEC;
    if(Intermediate.NumVectors() == 1)
      IntermediateView = SYMBOLIC_VEC;
    else
      IntermediateView = SYMBOLIC_MULTIVEC;
  }
  if(ATA)
    ReturnVal = oski_SetHintMatTransMatMult(A_tunable_, OP_AT_A, Alpha, InView, Beta, OutView, IntermediateView, NumCalls);
  else
    ReturnVal = oski_SetHintMatTransMatMult(A_tunable_, OP_A_AT, Alpha, InView, Beta, OutView, IntermediateView, NumCalls);
  if(ReturnVal)
    std::cerr << "Set hint mattransmat multiply error\n";
  return ReturnVal;
}

int Epetra_OskiMatrix::SetHintMultiplyAndMatTransMultiply(bool TransA, 
							  double Alpha, 
							  const Epetra_OskiMultiVector& InVec, 
							  double Beta, 
							  const Epetra_OskiMultiVector& OutVec, 
							  double Omega, 
							  const Epetra_OskiMultiVector& InVec2, 
							  double Zeta, 
							  const Epetra_OskiMultiVector& OutVec2, 
							  int NumCalls, 
							  const Teuchos::ParameterList& List) {
  int ReturnVal;
  oski_vecview_t InView = NULL;
  oski_vecview_t OutView = NULL;
  oski_vecview_t InView2 = NULL;
  oski_vecview_t OutView2 = NULL;
  InView = InVec.Oski_View();
  OutView = OutVec.Oski_View();
  InView2 = InVec2.Oski_View();
  OutView2 = OutVec2.Oski_View();
  if(List.isParameter("tune"))
    if(Teuchos::getParameter<bool>(List, "tune"))
      NumCalls = ALWAYS_TUNE;
  if(List.isParameter("tuneaggressive"))
    if(Teuchos::getParameter<bool>(List, "tuneaggressive"))
      NumCalls = ALWAYS_TUNE_AGGRESSIVELY;
  if(List.isParameter("symminvec"))
    if(Teuchos::getParameter<bool>(List, "symminvec"))
      InView = SYMBOLIC_VEC;
  if(List.isParameter("symminmultivec"))
    if(Teuchos::getParameter<bool>(List, "symminmultivec"))
      InView = SYMBOLIC_MULTIVEC;
  if(List.isParameter("symmoutvec"))
    if(Teuchos::getParameter<bool>(List, "symmoutvec"))
      OutView = SYMBOLIC_VEC;
  if(List.isParameter("symmoutmultivec"))
    if(Teuchos::getParameter<bool>(List, "symmoutmultivec"))
      OutView = SYMBOLIC_MULTIVEC;
  if(List.isParameter("symminvec2"))
    if(Teuchos::getParameter<bool>(List, "symminvec2"))
      InView2 = SYMBOLIC_VEC;
  if(List.isParameter("symminmultivec2"))
    if(Teuchos::getParameter<bool>(List, "symminmultivec2"))
      InView2 = SYMBOLIC_MULTIVEC;
  if(List.isParameter("symmoutvec2"))
    if(Teuchos::getParameter<bool>(List, "symmoutvec2"))
      OutView2 = SYMBOLIC_VEC;
  if(List.isParameter("symmoutmultivec2"))
    if(Teuchos::getParameter<bool>(List, "symmoutmultivec2"))
      OutView2 = SYMBOLIC_MULTIVEC;
  if(this->Comm().NumProc() > 1) {
    if(InVec.NumVectors() == 1)
      InView = SYMBOLIC_VEC;
    else
      InView = SYMBOLIC_MULTIVEC;
    if(OutVec.NumVectors() == 1)
      OutView = SYMBOLIC_VEC;
    else
      OutView = SYMBOLIC_MULTIVEC;
    if(InVec2.NumVectors() == 1)
      InView2 = SYMBOLIC_VEC;
    else
      InView2 = SYMBOLIC_MULTIVEC;
    if(OutVec2.NumVectors() == 1)
      OutView2 = SYMBOLIC_VEC;
    else
      OutView2 = SYMBOLIC_MULTIVEC;
  }
  if(TransA)
    ReturnVal = oski_SetHintMatMultAndMatTransMult(A_tunable_, Alpha, InView, Beta, OutView, OP_TRANS, Omega, InView2, Zeta, OutView2, NumCalls);
  else
    ReturnVal = oski_SetHintMatMultAndMatTransMult(A_tunable_, Alpha, InView, Beta, OutView, OP_NORMAL, Omega, InView2, Zeta, OutView2, NumCalls);
  if(ReturnVal)
    std::cerr << "Set hint multivector multiply and mattransmultiply error\n";
  return ReturnVal;
}

int Epetra_OskiMatrix::SetHintPowMultiply(bool TransA, 
					  double Alpha, 
					  const Epetra_OskiMultiVector& InVec, 
					  double Beta, 
					  const Epetra_OskiMultiVector& OutVec, 
					  const Epetra_OskiMultiVector& Intermediate, 
					  int Power, 
					  int NumCalls, 
					  const Teuchos::ParameterList& List) {
  int ReturnVal;
  oski_vecview_t InView = NULL;
  oski_vecview_t OutView = NULL;
  oski_vecview_t IntermediateView = NULL;
  InView = InVec.Oski_View();
  OutView = OutVec.Oski_View();
  if(&Intermediate != NULL)
    IntermediateView = Intermediate.Oski_View();
  if(List.isParameter("tune"))
    if(Teuchos::getParameter<bool>(List, "tune"))
      NumCalls = ALWAYS_TUNE;
  if(List.isParameter("tuneaggressive"))
    if(Teuchos::getParameter<bool>(List, "tuneaggressive"))
      NumCalls = ALWAYS_TUNE_AGGRESSIVELY;
  if(List.isParameter("symminvec"))
    if(Teuchos::getParameter<bool>(List, "symminvec"))
      InView = SYMBOLIC_VEC;
  if(List.isParameter("symminmultivec"))
    if(Teuchos::getParameter<bool>(List, "symminmultivec"))
      InView = SYMBOLIC_MULTIVEC;
  if(List.isParameter("symmoutvec"))
    if(Teuchos::getParameter<bool>(List, "symmoutvec"))
      OutView = SYMBOLIC_VEC;
  if(List.isParameter("symmoutmultivec"))
    if(Teuchos::getParameter<bool>(List, "symmoutmultivec"))
      OutView = SYMBOLIC_MULTIVEC;
  if(List.isParameter("symmintervec"))
    if(Teuchos::getParameter<bool>(List, "symmintermultivec"))
      IntermediateView = SYMBOLIC_MULTIVEC;
  if(List.isParameter("invalidinter"))
    if(Teuchos::getParameter<bool>(List, "invalidinter"))
      IntermediateView = NULL;
  if(this->Comm().NumProc() > 1) {
    if(InVec.NumVectors() == 1)
      InView = SYMBOLIC_VEC;
    else
      InView = SYMBOLIC_MULTIVEC;
    if(OutVec.NumVectors() == 1)
      OutView = SYMBOLIC_VEC;
    else
      OutView = SYMBOLIC_MULTIVEC;
    if(Intermediate.NumVectors() == 1)
      IntermediateView = SYMBOLIC_VEC;
    else
      IntermediateView = SYMBOLIC_MULTIVEC;
  }
  if(TransA)
    ReturnVal = oski_SetHintMatPowMult(A_tunable_, OP_TRANS, Power, Alpha, InView, Beta, OutView, IntermediateView, NumCalls);
  else
    ReturnVal = oski_SetHintMatPowMult(A_tunable_, OP_NORMAL, Power, Alpha, InView, Beta, OutView, IntermediateView, NumCalls);
  if(ReturnVal)
    std::cerr << "Set hint matpow multiply error\n";
  return ReturnVal;
}


int Epetra_OskiMatrix::TuneMatrix() {
  int ReturnVal;
  ReturnVal = oski_TuneMat(A_tunable_);
  if(IsMatrixTransformed())
    Copy_Created_ = true;
  return ReturnVal;
}

int Epetra_OskiMatrix::IsMatrixTransformed() const {
  return oski_IsMatPermuted(A_tunable_);
}

const Epetra_OskiMatrix& Epetra_OskiMatrix::ViewTransformedMat() const {
  //might need a more efficient way to do this
  Epetra_OskiMatrix* Transformed = NULL;
  Epetra_OskiMatrix Temp(*this); //should be some lightweight copy
  Transformed = &Temp;
  Transformed->A_tunable_ = const_cast <oski_matrix_t> (oski_ViewPermutedMat(A_tunable_));
  return *Transformed;
}

const Epetra_OskiPermutation& Epetra_OskiMatrix::ViewRowPermutation() const {
  Epetra_OskiPermutation* RowPerm = NULL;
  RowPerm = new Epetra_OskiPermutation[1];
  (*RowPerm).ReplacePermutation(const_cast <oski_perm_t> (oski_ViewPermutedMatRowPerm(A_tunable_)));
  return *RowPerm;
}

const Epetra_OskiPermutation& Epetra_OskiMatrix::ViewColumnPermutation() const {
  Epetra_OskiPermutation* ColPerm;
  ColPerm = new Epetra_OskiPermutation[1];
  (*ColPerm).ReplacePermutation(const_cast <oski_perm_t> (oski_ViewPermutedMatColPerm(A_tunable_)));
  return *ColPerm;
}

char* Epetra_OskiMatrix::GetMatrixTransforms() const {
  char* ReturnVal = NULL;
  ReturnVal = oski_GetMatTransforms(A_tunable_);
  if(ReturnVal == NULL)
    std::cerr << "Error in GetMatrixTransforms\n";
  return ReturnVal;
}

int Epetra_OskiMatrix::ApplyMatrixTransforms(const char* Transforms) {
  int ReturnVal;
  ReturnVal = oski_ApplyMatTransforms(A_tunable_, Transforms);
  if(ReturnVal)
    std::cerr << "Error in ApplyMatrixTransforms\n";
  return ReturnVal;
}
#endif
#endif
