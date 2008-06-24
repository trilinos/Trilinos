
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

#include "Epetra_ConfigDefs.h"

#ifdef WITH_EPETRA_PRERELEASE
#ifdef HAVE_OSKI
#include "Epetra_OskiMatrix.h"

//=============================================================================

Epetra_OskiMatrix::Epetra_OskiMatrix(const Epetra_OskiMatrix& Source) 
  : Epetra_CrsMatrix(Source), 
  Copy_Created_(true), 
  Epetra_View_(Source.Epetra_View_) {
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
  if(TransA)
    ReturnVal = oski_MatMult(A_tunable_, OP_TRANS, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View());
  else
    ReturnVal = oski_MatMult(A_tunable_, OP_NORMAL, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View());
  if(ReturnVal)
    std::cerr << "OskiVector multiply error\n";
  if(xCreate)
    delete xCast;
  if(yCreate)
    delete yCast;
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
  Epetra_OskiMultiVector* XCast = NULL;
  Epetra_OskiMultiVector* YCast = NULL;
  bool XCreate = false;
  bool YCreate = false;
  int ReturnVal;
  XCast = dynamic_cast<Epetra_OskiMultiVector*>(const_cast <Epetra_MultiVector*>(&X));
  YCast = dynamic_cast<Epetra_OskiMultiVector*>(&Y);
  if (XCast == NULL) {
    XCast = new Epetra_OskiMultiVector(X);
    XCreate = true;
  }
  if (YCast == NULL) {
    YCast = new Epetra_OskiMultiVector(Y);
    YCreate = true;
  }
  if(TransA)
    ReturnVal = oski_MatMult(A_tunable_, OP_TRANS, Alpha, (*XCast).Oski_View(), Beta, (*YCast).Oski_View());
  else
    ReturnVal = oski_MatMult(A_tunable_, OP_NORMAL, Alpha, (*XCast).Oski_View(), Beta, (*YCast).Oski_View());
  if(ReturnVal)
    std::cerr << "OskiMultiVector multiply error\n";
  if(XCreate)
    delete XCast;
  if(YCreate)
    delete YCast;
  return ReturnVal;
}

int Epetra_OskiMatrix::Solve(bool Upper, bool TransA, bool UnitDiagonal, const Epetra_Vector& x, Epetra_Vector& y) const {
  int ReturnVal;
  ReturnVal = this->Solve(TransA, x, y, 1.0);
  return ReturnVal;
}

int Epetra_OskiMatrix::Solve(bool TransA, const Epetra_Vector& x, Epetra_Vector& y, double Alpha) const {
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
					   Epetra_Vector& t, 
					   double Alpha, 
					   double Beta) const {
  Epetra_OskiVector* xCast = NULL;
  Epetra_OskiVector* yCast = NULL;
  Epetra_OskiVector* tCast = NULL;
  bool xCreate = false;
  bool yCreate = false;
  bool tCreate = false;
  int ReturnVal;
  tCast = dynamic_cast<Epetra_OskiVector*>(&t);
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
  if (tCast == NULL) {
    tCast = new Epetra_OskiVector(t);
    tCreate = true;
  }
  if(ATA)
    ReturnVal = oski_MatTransMatMult(A_tunable_, OP_AT_A, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), (*tCast).Oski_View());
  else
    ReturnVal = oski_MatTransMatMult(A_tunable_, OP_A_AT, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), (*tCast).Oski_View());
  if(ReturnVal)
    std::cerr << "OskiVector MatTransMatMultiply error\n";
  if(xCreate)
    delete xCast;
  if(yCreate)
    delete yCast;
  if(tCreate)
    delete tCast;
  return ReturnVal;
}

int Epetra_OskiMatrix::MatTransMatMultiply(bool ATA, 
					   const Epetra_MultiVector& X, 
					   Epetra_MultiVector& Y, 
					   Epetra_MultiVector& T, 
					   double Alpha, 
				           double Beta) const {
  Epetra_OskiMultiVector* XCast = NULL;
  Epetra_OskiMultiVector* YCast = NULL;
  Epetra_OskiMultiVector* TCast = NULL;
  bool XCreate = false;
  bool YCreate = false;
  bool TCreate = false;
  int ReturnVal;
  TCast = dynamic_cast<Epetra_OskiMultiVector*>(&T);
  XCast = dynamic_cast<Epetra_OskiMultiVector*>(const_cast <Epetra_MultiVector*>(&X));
  YCast = dynamic_cast<Epetra_OskiMultiVector*>(&Y);
  if (XCast == NULL) {
    XCast = new Epetra_OskiMultiVector(X);
    XCreate = true;
  }
  if (YCast == NULL) {
    YCast = new Epetra_OskiMultiVector(Y);
    YCreate = true;
  }
  if (TCast == NULL && (&T) != NULL) {
    TCast = new Epetra_OskiMultiVector(T);
    TCreate = true;
  }
  if(ATA)
    ReturnVal = oski_MatTransMatMult(A_tunable_, OP_AT_A, Alpha, (*XCast).Oski_View(), Beta, (*YCast).Oski_View(), (*TCast).Oski_View());
  else
    ReturnVal = oski_MatTransMatMult(A_tunable_, OP_A_AT, Alpha, (*XCast).Oski_View(), Beta, (*YCast).Oski_View(), (*TCast).Oski_View());
  if(ReturnVal)
    std::cerr << "OskiMultiVector MatTransMatMultiply error\n";
  if(TCreate)
    delete TCast;
  if(XCreate)
    delete XCast;
  if(YCreate)
    delete YCast;
  return ReturnVal;
}

int Epetra_OskiMatrix::MatTransMatMultiply(bool ATA, 
					   const Epetra_Vector& x, 
					   Epetra_Vector& y, 
					   double Alpha, 
					   double Beta) const {
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
  if(ATA)
    ReturnVal = oski_MatTransMatMult(A_tunable_, OP_AT_A, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), INVALID_VEC);
  else
    ReturnVal = oski_MatTransMatMult(A_tunable_, OP_A_AT, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), INVALID_VEC);
  if(ReturnVal)
    std::cerr << "OskiVector MatTransMatMultiply error\n";
  if(xCreate)
    delete xCast;
  if(yCreate)
    delete yCast;
  return ReturnVal;
}

int Epetra_OskiMatrix::MatTransMatMultiply(bool ATA, 
					   const Epetra_MultiVector& X, 
					   Epetra_MultiVector& Y, 
					   double Alpha, 
				           double Beta) const {
  Epetra_OskiMultiVector* XCast = NULL;
  Epetra_OskiMultiVector* YCast = NULL;
  bool XCreate = false;
  bool YCreate = false;
  int ReturnVal;
  XCast = dynamic_cast<Epetra_OskiMultiVector*>(const_cast <Epetra_MultiVector*>(&X));
  YCast = dynamic_cast<Epetra_OskiMultiVector*>(&Y);
  if (XCast == NULL) {
    XCast = new Epetra_OskiMultiVector(X);
    XCreate = true;
  }
  if (YCast == NULL) {
    YCast = new Epetra_OskiMultiVector(Y);
    YCreate = true;
  }
  if(ATA)
    ReturnVal = oski_MatTransMatMult(A_tunable_, OP_AT_A, Alpha, (*XCast).Oski_View(), Beta, (*YCast).Oski_View(), INVALID_VEC);
  else
    ReturnVal = oski_MatTransMatMult(A_tunable_, OP_A_AT, Alpha, (*XCast).Oski_View(), Beta, (*YCast).Oski_View(), INVALID_VEC);
  if(ReturnVal)
    std::cerr << "OskiMultiVector MatTransMatMultiply error\n";
  if(XCreate)
    delete XCast;
  if(YCreate)
    delete YCast;
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
  Epetra_OskiVector* xCast = NULL;
  Epetra_OskiVector* yCast = NULL;
  Epetra_OskiVector* wCast = NULL;
  Epetra_OskiVector* zCast = NULL;
  bool xCreate = false;
  bool yCreate = false;
  bool wCreate = false;
  bool zCreate = false;
  int ReturnVal;
  xCast = dynamic_cast<Epetra_OskiVector*>(const_cast <Epetra_Vector*>(&x));
  yCast = dynamic_cast<Epetra_OskiVector*>(&y);
  wCast = dynamic_cast<Epetra_OskiVector*>(const_cast <Epetra_Vector*>(&w));
  zCast = dynamic_cast<Epetra_OskiVector*>(&z);
  if (xCast == NULL) {
    xCast = new Epetra_OskiVector(x);
    xCreate = true;
  }
  if (yCast == NULL) {
    yCast = new Epetra_OskiVector(y);
    yCreate = true;
  }
  if (wCast == NULL) {
    wCast = new Epetra_OskiVector(w);
    wCreate = true;
  }
  if (zCast == NULL) {
    zCast = new Epetra_OskiVector(z);
    zCreate = true;
  }
  if(TransA)
    ReturnVal = oski_MatMultAndMatTransMult(A_tunable_, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), OP_TRANS, Omega, (*wCast).Oski_View(), Zeta, (*zCast).Oski_View());
  else
    ReturnVal = oski_MatMultAndMatTransMult(A_tunable_, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), OP_NORMAL, Omega, (*wCast).Oski_View(), Zeta, (*zCast).Oski_View());
  if(ReturnVal)
    std::cerr << "OskiVector multiply error\n";
  if(xCreate)
    delete xCast;
  if(yCreate)
    delete yCast;
  if(wCreate)
    delete wCast;
  if(zCreate)
    delete zCast;
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
  Epetra_OskiMultiVector* XCast = NULL;
  Epetra_OskiMultiVector* YCast = NULL;
  Epetra_OskiMultiVector* WCast = NULL;
  Epetra_OskiMultiVector* ZCast = NULL;
  bool XCreate = false;
  bool YCreate = false;
  bool WCreate = false;
  bool ZCreate = false;
  int ReturnVal;
  XCast = dynamic_cast<Epetra_OskiMultiVector*>(const_cast <Epetra_MultiVector*>(&X));
  YCast = dynamic_cast<Epetra_OskiMultiVector*>(&Y);
  WCast = dynamic_cast<Epetra_OskiMultiVector*>(const_cast <Epetra_MultiVector*>(&W));
  ZCast = dynamic_cast<Epetra_OskiMultiVector*>(&Z);
  if (XCast == NULL) {
    XCast = new Epetra_OskiMultiVector(X);
    XCreate = true;
  }
  if (YCast == NULL) {
    YCast = new Epetra_OskiMultiVector(Y);
    YCreate = true;
  }
  if (WCast == NULL) {
    WCast = new Epetra_OskiMultiVector(W);
    WCreate = true;
  }
  if (ZCast == NULL) {
    ZCast = new Epetra_OskiMultiVector(Z);
    ZCreate = true;
  }
  if(TransA)
    ReturnVal = oski_MatMultAndMatTransMult(A_tunable_, Alpha, (*XCast).Oski_View(), Beta, (*YCast).Oski_View(), OP_TRANS, Omega, (*WCast).Oski_View(), Zeta, (*ZCast).Oski_View());
  else
    ReturnVal = oski_MatMultAndMatTransMult(A_tunable_, Alpha, (*XCast).Oski_View(), Beta, (*YCast).Oski_View(), OP_NORMAL, Omega, (*WCast).Oski_View(), Zeta, (*ZCast).Oski_View());
  if(ReturnVal)
    std::cerr << "OskiVector multiply error\n";
  if(XCreate)
    delete XCast;
  if(YCreate)
    delete YCast;
  if(WCreate)
    delete WCast;
  if(ZCreate)
    delete ZCast;
  return ReturnVal;
}

int Epetra_OskiMatrix::MatPowMultiply(bool TransA, 
				      const Epetra_Vector& x, 
				      Epetra_Vector& y, 
 				      Epetra_MultiVector& T, 
				      int Power,
				      double Alpha, 
				      double Beta) const {
  Epetra_OskiVector* xCast = NULL;
  Epetra_OskiVector* yCast = NULL;
  Epetra_OskiMultiVector* TCast = NULL;
  bool xCreate = false;
  bool yCreate = false;
  bool TCreate = false;
  int ReturnVal;
  TCast = dynamic_cast<Epetra_OskiMultiVector*>(&T);
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
  if (TCast == NULL) {
    TCast = new Epetra_OskiMultiVector(T);
    TCreate = true;
  }
  if(TransA)
    ReturnVal = oski_MatPowMult(A_tunable_, OP_TRANS, Power, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), (*TCast).Oski_View());
  else
    ReturnVal = oski_MatPowMult(A_tunable_, OP_NORMAL, Power, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), (*TCast).Oski_View());
  if(ReturnVal)
    std::cerr << "OskiVector matpow multiply error\n";
  if(xCreate)
    delete xCast;
  if(yCreate)
    delete yCast;
  if(TCreate)
    delete TCast;
  return ReturnVal;
}

int Epetra_OskiMatrix::MatPowMultiply(bool TransA,
				      const Epetra_MultiVector& X,
				      Epetra_MultiVector& Y, 
				      Epetra_MultiVector& T, 
				      int Power, 
				      double Alpha, 
				      double Beta) const {
  Epetra_OskiMultiVector* XCast = NULL;
  Epetra_OskiMultiVector* YCast = NULL;
  Epetra_OskiMultiVector* TCast = NULL;
  bool XCreate = false;
  bool YCreate = false;
  bool TCreate = false;
  int ReturnVal;
  TCast = dynamic_cast<Epetra_OskiMultiVector*>(&T);
  XCast = dynamic_cast<Epetra_OskiMultiVector*>(const_cast <Epetra_MultiVector*>(&X));
  YCast = dynamic_cast<Epetra_OskiMultiVector*>(&Y);
  if (XCast == NULL) {
    XCast = new Epetra_OskiMultiVector(X);
    XCreate = true;
  }
  if (YCast == NULL) {
    YCast = new Epetra_OskiMultiVector(Y);
    YCreate = true;
  }
  if (TCast == NULL) {
    TCast = new Epetra_OskiMultiVector(T);
    TCreate = true;
  }
  if(TransA)
    ReturnVal = oski_MatPowMult(A_tunable_, OP_TRANS, Power, Alpha, (*XCast).Oski_View(), Beta, (*YCast).Oski_View(), (*TCast).Oski_View());
  else
    ReturnVal = oski_MatPowMult(A_tunable_, OP_NORMAL, Power, Alpha, (*XCast).Oski_View(), Beta, (*YCast).Oski_View(), (*TCast).Oski_View());
  if(ReturnVal)
    std::cerr << "OskiVector matpow multiply error\n";
  if(XCreate)
    delete XCast;
  if(YCreate)
    delete YCast;
  if(TCreate)
    delete TCast;
  return ReturnVal;
}

int Epetra_OskiMatrix::MatPowMultiply(bool TransA, 
				      const Epetra_Vector& x, 
				      Epetra_Vector& y, 
				      int Power,
				      double Alpha, 
				      double Beta) const {
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
  if(TransA)
    ReturnVal = oski_MatPowMult(A_tunable_, OP_TRANS, Power, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), NULL);
  else
    ReturnVal = oski_MatPowMult(A_tunable_, OP_NORMAL, Power, Alpha, (*xCast).Oski_View(), Beta, (*yCast).Oski_View(), NULL);
  if(ReturnVal)
    std::cerr << "OskiVector matpow multiply error\n";
  if(xCreate)
    delete xCast;
  if(yCreate)
    delete yCast;
  return ReturnVal;
}

int Epetra_OskiMatrix::MatPowMultiply(bool TransA,
				      const Epetra_MultiVector& X,
				      Epetra_MultiVector& Y, 
				      int Power, 
				      double Alpha, 
				      double Beta) const {
  Epetra_OskiMultiVector* XCast = NULL;
  Epetra_OskiMultiVector* YCast = NULL;
  bool XCreate = false;
  bool YCreate = false;
  int ReturnVal;
  XCast = dynamic_cast<Epetra_OskiMultiVector*>(const_cast <Epetra_MultiVector*>(&X));
  YCast = dynamic_cast<Epetra_OskiMultiVector*>(&Y);
  if (XCast == NULL) {
    XCast = new Epetra_OskiMultiVector(X);
    XCreate = true;
  }
  if (YCast == NULL) {
    YCast = new Epetra_OskiMultiVector(Y);
    YCreate = true;
  }
  if(TransA)
    ReturnVal = oski_MatPowMult(A_tunable_, OP_TRANS, Power, Alpha, (*XCast).Oski_View(), Beta, (*YCast).Oski_View(), NULL);
  else
    ReturnVal = oski_MatPowMult(A_tunable_, OP_NORMAL, Power, Alpha, (*XCast).Oski_View(), Beta, (*YCast).Oski_View(), NULL);
  if(ReturnVal)
    std::cerr << "OskiVector matpow multiply error\n";
  if(XCreate)
    delete XCast;
  if(YCreate)
    delete YCast;
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
      if(ReturnVal = oski_SetHint(A_tunable_, HINT_ALIGNED_BLOCKS))
        std::cerr << "Error setting hint aligned blocks.\n";
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
        if(ReturnVal = oski_SetHint(A_tunable_, HINT_SINGLE_BLOCKSIZE, ArgArray))
          std::cerr << "Error setting hint no blocks.\n";
      }
      else
        if(ReturnVal = oski_SetHint(A_tunable_, HINT_SINGLE_BLOCKSIZE, ArgArray))
          std::cerr << "Error setting hint no blocks.\n";
      delete ArgArray;
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
        if(ReturnVal = oski_SetHint(A_tunable_, HINT_MULTIPLE_BLOCKSIZES, ArgArray))
          std::cerr << "Error setting hint multiple blocks.\n";
        delete ArgArray;
        ArgArray = NULL;
      }
      else
        if(ReturnVal = oski_SetHint(A_tunable_, HINT_MULTIPLE_BLOCKSIZES, ArgArray))
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
        if(ReturnVal = oski_SetHint(A_tunable_, HINT_DIAGS, ArgArray))
          std::cerr << "Error setting hint no blocks.\n";
        delete ArgArray;
      }
      else
      {
        if(ReturnVal = oski_SetHint(A_tunable_, HINT_DIAGS, ArgArray))
          std::cerr << "Error setting hint no blocks.\n";
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
