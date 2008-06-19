
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

// Epetra_BlockMap Constructor

/*Epetra_OskiMultiVector::Epetra_Matrix(const Epetra_Matrix& Source) {

}*/

//need to read in autotune or not.  Input mode.  Any overridden params.
Epetra_OskiMatrix::Epetra_OskiMatrix(const Epetra_CrsMatrix& Source, const Teuchos::ParameterList& List) : Epetra_CrsMatrix(Source), Epetra_View_(&Source) {
 
  bool AutoTune = false;
  bool DeepCopy = false;
  char Matrix[20] = "general\0";
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
  int* RowPtr;
  int* IndPtr;
  double* ValPtr;
  if(List.isParameter("autotune")) 
    AutoTune = List.get<bool>("autotune");
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
    strcpy(Matrix, Teuchos::getParameter<char*>(List, "matrixtype"));
    if(!strcmp(Matrix, "general"))
      MatrixType = MAT_GENERAL;
    else if(!strcmp(Matrix, "uppertri"))
      MatrixType = MAT_TRI_UPPER;
    else if(!strcmp(Matrix, "lowertri"))
      MatrixType = MAT_TRI_LOWER;
    else if(!strcmp(Matrix, "uppersymm"))
      MatrixType = MAT_SYMM_UPPER;
    else if(!strcmp(Matrix, "lowersymm"))
      MatrixType = MAT_SYMM_LOWER;
    else if(!strcmp(Matrix, "fullsymm"))
      MatrixType = MAT_SYMM_FULL;
    else if(!strcmp(Matrix, "upperherm"))
      MatrixType = MAT_HERM_UPPER;
    else if(!strcmp(Matrix, "lowerherm"))
      MatrixType = MAT_HERM_LOWER;
    else if(!strcmp(Matrix, "fullherm"))
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

Epetra_OskiMatrix::~Epetra_OskiMatrix (){
}

int Epetra_OskiMatrix::ReplaceMyValues(int MyRow, int NumEntries, double* Values, int* Indices) {

}

int Epetra_OskiMatrix::SumIntoMyValues(int MyRow, int NumEntries, double* Values, int* Indices) {

}

int Epetra_OskiMatrix::ExtractDiagonalCopy(Epetra_OskiVector& Diagonal) const {

}

int Epetra_OskiMatrix::ReplaceDiagonalValues(const Epetra_OskiVector& Diagonal) {

}

int Epetra_OskiMatrix::Multiply(bool TransA, const Epetra_Vector& x, Epetra_Vector& y, double Alpha, double Beta) const {

}

int Epetra_OskiMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y, double Alpha, double Beta) const {
  const Epetra_OskiMultiVector* XCast;
  Epetra_OskiMultiVector* YCast;
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
  if(XCreate)
    delete XCast;
  if(YCreate)
    delete YCast;
  if(ReturnVal)
    std::cerr << "MultiVector multiply error\n";
  return ReturnVal;
}

int Epetra_OskiMatrix::SetHint(const Teuchos::ParameterList& List){
}

int Epetra_OskiMatrix::SetHintMatrixMultiply(bool TransA, double Alpha, const Epetra_OskiMultiVector InVec, double Beta, const Epetra_OskiMultiVector OutVec, int NumCalls, const Teuchos::ParameterList& List){
/*  int ReturnVal;
  
  if(Trans)
    ReturnVal = oski_SetHintMatMult(A_tunable_, OP_TRANS, Alpha, some invec might be symbolic, Beta, some outvec might be symbolic, numcalls might be symbolic);
  else
    ReturnVal = oski_SetHintMatMult(A_tunable_, OP_NORMAL, Alpha, some invec might be symbolic, Beta, some outvec might be symbolic, numcalls might be symbolic);
  if(ReturnVal)
    std::cerr << "Set hint multivector multiply error\n";
  return ReturnVal;*/
}

int Epetra_OskiMatrix::TuneMatrix() {
  return oski_TuneMat(A_tunable_);
}

int Epetra_OskiMatrix::IsMatrixTransformed() const {
  return oski_IsMatPermuted(A_tunable_);
}

const Epetra_OskiMatrix& Epetra_OskiMatrix::ViewTransformedMat() const {

}

const Epetra_OskiPermutation& Epetra_OskiMatrix::ViewRowPermutation() const {

}

const Epetra_OskiPermutation& Epetra_OskiMatrix::ViewColumnPermutation() const {

}

char* Epetra_OskiMatrix::GetMatrixTransforms() const {

}

int Epetra_OskiMatrix::ApplyMatrixTransforms(const char* Transforms) {

}
#endif
#endif
