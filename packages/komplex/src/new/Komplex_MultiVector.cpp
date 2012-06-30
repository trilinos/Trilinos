//@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
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

#include "Komplex_MultiVector.hpp" 
#include "Epetra_DataAccess.h"
#include "Komplex_Ordering.hpp"
#include "Komplex_KForms.hpp"  
#include "Epetra_MultiVector.h" 
#include "Epetra_BlockMap.h" 
#include "Epetra_Vector.h"
#include "Epetra_Comm.h" 

//==========================================================================
Komplex_MultiVector::Komplex_MultiVector(const Epetra_BlockMap& Map, 
						     int NumVectors, 
						     bool RHS,
						     bool zeroOut,
						     Komplex_KForms KForm)
  : Imag_(0),
    IsOneObject_(true),
    RHS_(RHS),
    OtherMap_(0),
    TempDouble_(0)
{
  Real_ = new Epetra_MultiVector(Map, NumVectors, zeroOut);
  Ordering_ = new Komplex_Ordering(Map, KForm, true);
  CreateOtherMap();
}

//==========================================================================
Komplex_MultiVector::Komplex_MultiVector(const Komplex_MultiVector& Source)
  : Real_(Source.Real_),
    Imag_(Source.Imag_),
    Ordering_(Source.Ordering_),
    IsOneObject_(Source.IsOneObject_),
    RHS_(Source.RHS_),
    OtherMap_(Source.OtherMap_),
    TempDouble_(Source.TempDouble_)
{
}

//==========================================================================  
Komplex_MultiVector::Komplex_MultiVector(Epetra_DataAccess CV, 
						     const Epetra_BlockMap& Map, 
						     double* A, 
                                         int MyLDA, 
						     int NumVectors,
						     bool RHS, 
						     Komplex_KForms KForm)
  : Imag_(0),
    IsOneObject_(true),
    RHS_(RHS),
    OtherMap_(0),
    TempDouble_(0)
{
  Real_ = new Epetra_MultiVector(CV, Map, A, MyLDA, NumVectors);
  Ordering_ = new Komplex_Ordering(Map, KForm, true);
  CreateOtherMap();
}

//==========================================================================
Komplex_MultiVector::Komplex_MultiVector(Epetra_DataAccess CV, 
						     const Epetra_BlockMap& Map,
						     double* Real,
		                             double* Imag, 
						     int MyLDA, 
						     int NumVectors,
						     bool RHS,
						     Komplex_KForms KForm)
  : IsOneObject_(false),
    RHS_(RHS),
    OtherMap_(0),
    TempDouble_(0)
{
  Real_ = new Epetra_MultiVector(CV, Map, Real, MyLDA, NumVectors);
  Imag_ = new Epetra_MultiVector(CV, Map, Imag, MyLDA, NumVectors);
  Imag_->SetSeed(Real_->Seed());
  Ordering_ = new Komplex_Ordering(Map, KForm, false);
  CreateOtherMap();
}

//==========================================================================
Komplex_MultiVector::Komplex_MultiVector(Epetra_DataAccess CV,
						     const Epetra_BlockMap& Map,
					           double** ArrayOfPointers, 
                                         int NumVectors, 
						     bool RHS, 
						     Komplex_KForms KForm)
  : Imag_(0),
    IsOneObject_(true),
    RHS_(RHS), 
    OtherMap_(0),
    TempDouble_(0)
{
  Real_ = new Epetra_MultiVector(CV, Map, ArrayOfPointers, NumVectors);
  Ordering_ = new Komplex_Ordering(Map, KForm, true);
  CreateOtherMap();
}

//==========================================================================
Komplex_MultiVector::Komplex_MultiVector(Epetra_DataAccess CV,
						     const Epetra_BlockMap& Map,
						     double** AOPReal, 
						     double** AOPImag, 
                                         int NumVectors, 
						     bool RHS, 
						     Komplex_KForms KForm)
  : IsOneObject_(false),
    RHS_(RHS),
    OtherMap_(0),
    TempDouble_(0)
{
  Real_ = new Epetra_MultiVector(CV, Map, AOPReal, NumVectors);
  Imag_ = new Epetra_MultiVector(CV, Map, AOPImag, NumVectors);
  Imag_->SetSeed(Real_->Seed());
  Ordering_ = new Komplex_Ordering(Map, KForm, true);
  CreateOtherMap();
}

//==========================================================================
Komplex_MultiVector::Komplex_MultiVector(Epetra_DataAccess CV,
						     const Epetra_MultiVector& Source, 
					           int* Indices, 
						     int NumVectors, 
						     bool RHS,
						     Komplex_KForms KForm)
  : Imag_(0),
    IsOneObject_(true),
    RHS_(RHS), 
    OtherMap_(0),
    TempDouble_(0)
{
  Real_ = new Epetra_MultiVector(CV, Source, Indices, NumVectors);
  Ordering_ = new Komplex_Ordering(Source.Map(), KForm, true);
  CreateOtherMap();
}

//==========================================================================
Komplex_MultiVector::Komplex_MultiVector(Epetra_DataAccess CV,
						     const Epetra_MultiVector& Source, 
                                         int StartIndex, 
						     int NumVectors, 
						     bool RHS, 
						     Komplex_KForms KForm)
  : Imag_(0),
    IsOneObject_(true),
    RHS_(RHS),
    OtherMap_(0),
    TempDouble_(0)
{
  Real_ = new Epetra_MultiVector(CV, Source, StartIndex, NumVectors);
  Ordering_ = new Komplex_Ordering(Source.Map(), KForm, true);
  CreateOtherMap();
}  

//==========================================================================
Komplex_MultiVector::Komplex_MultiVector(Epetra_DataAccess CV,
						     const Epetra_MultiVector& Source, 
                                         bool RHS, 
						     Komplex_KForms KForm)
  : Imag_(0),
    IsOneObject_(true),
    RHS_(RHS),
    OtherMap_(0),
    TempDouble_(0)
{ 
  Real_ = new Epetra_MultiVector(CV, Source, 0, Source.NumVectors());
  Ordering_ = new Komplex_Ordering(Source.Map(), KForm, true);
  CreateOtherMap();
}

//==========================================================================
Komplex_MultiVector::Komplex_MultiVector(Epetra_DataAccess CV, 
						     const Epetra_MultiVector& Real, 
                                         const Epetra_MultiVector& Imag, 
						     bool RHS, 
						     Komplex_KForms KForm)
  : IsOneObject_(false),
    RHS_(RHS),
    OtherMap_(0),
    TempDouble_(0)
{ 
  Real_ = new Epetra_MultiVector(CV, Real, 0, Real.NumVectors());
  Imag_ = new Epetra_MultiVector(CV, Imag, 0, Imag.NumVectors());
  Imag_->SetSeed(Real_->Seed());
  Ordering_ = new Komplex_Ordering(Real.Map(), KForm, true);
  CreateOtherMap();
}

//==========================================================================  
Komplex_MultiVector::~Komplex_MultiVector() {
  delete Real_;
  if (Imag_ != 0) {
    delete Imag_;
  }
  delete Ordering_;
  delete OtherMap_;
  if (TempDouble_ != 0) {
    delete [] TempDouble_;
  }
}

//==========================================================================
int Komplex_MultiVector::ReplaceGlobalValue(int GlobalRow,
						        int VectorIndex, 
							  double ScalarValue) {
  if (IsOneObject_) {
    if (RHS_) { //it's a right-hand side multivector
      int Index = 0;
      int error = Ordering_->GlobalIndex(GlobalRow, Index);
      if (error != 0) { //probably want to keep this error checking in because user could have typed
                        //an index that is too big or negative or not belonging to "me" or something
        return error;
      } 
      if (Index == -1) {
        if (GlobalRow % 2 == 0) {
          return Real_->ReplaceGlobalValue(GlobalRow+1, VectorIndex, ScalarValue);
        } else {
          return Real_->ReplaceGlobalValue(GlobalRow-1, VectorIndex, ScalarValue);
        }
      } else {
	  return Real_->ReplaceGlobalValue(GlobalRow, VectorIndex, ScalarValue);
      }
    } else { //it's a left-hand side multivector
      double Scaling = 1.0;
      Ordering_->GlobalScaling(GlobalRow, Scaling);
      if (Scaling != 1.0) {
	  return Real_->ReplaceGlobalValue(GlobalRow, VectorIndex, ScalarValue/Scaling);
      } else {
	  return Real_->ReplaceGlobalValue(GlobalRow, VectorIndex, ScalarValue);
      }
    }
  } else { //it's two objects
    if (RHS_) {
      int Index = 0;
      int error = Ordering_->GlobalIndex(GlobalRow, Index);
      if (error != 0) { //keep, for reasons above =)
	  return error;
      }
      if (Index == -1) {
	  return Imag_->ReplaceGlobalValue(GlobalRow/2, VectorIndex, ScalarValue);
      } else {
	  return Real_->ReplaceGlobalValue(GlobalRow/2, VectorIndex, ScalarValue);
      }
    } else {  //it's a left-hand side multivector
      double Scaling = 1.0;
      Ordering_->GlobalScaling(GlobalRow, Scaling);
      if (GlobalRow % 2 == 0) {
        return Real_->ReplaceGlobalValue(GlobalRow/2, VectorIndex, ScalarValue/Scaling);
      } else {
	  return Imag_->ReplaceGlobalValue(GlobalRow/2, VectorIndex, ScalarValue/Scaling);
      }
    }
  }
}

//==========================================================================
int Komplex_MultiVector::SumIntoGlobalValue(int GlobalRow, 
							  int VectorIndex,
							  double ScalarValue) {
  if (IsOneObject_) {
    if (RHS_) { //it's a right-hand side multivector
      int Index = 0;
      int error = Ordering_->GlobalIndex(GlobalRow, Index);
      if (error != 0) { //probably want to keep this error checking in because user could have typed
                        //an index that is too big or negative or not belonging to "me" or something
        return error;
      }
      if (Index == -1) {
	  if (GlobalRow % 2 == 0) {
          return Real_->SumIntoGlobalValue(GlobalRow+1, VectorIndex, ScalarValue);
        } else {
          return Real_->SumIntoGlobalValue(GlobalRow-1, VectorIndex, ScalarValue);
        }
      } else {
        return Real_->SumIntoGlobalValue(GlobalRow, VectorIndex, ScalarValue);
      } 
    } else { //it's a left-hand side multivector
      double Scaling = 1.0;
      Ordering_->GlobalScaling(GlobalRow, Scaling);
      if (Scaling != 1.0) {
        return Real_->SumIntoGlobalValue(GlobalRow, VectorIndex, ScalarValue/Scaling);
      } else {
        return Real_->SumIntoGlobalValue(GlobalRow, VectorIndex, ScalarValue);
      }
    }
  } else { //it's two objects
    if (RHS_) {
      int Index = 0;
      int error = Ordering_->GlobalIndex(GlobalRow/2, Index);
      if (error != 0) { //keep, for reasons above =)
        return error;
      }
      if (Index == -1) {
	  return Imag_->SumIntoGlobalValue(GlobalRow/2, VectorIndex, ScalarValue);
      } else {
	  return Real_->SumIntoGlobalValue(GlobalRow/2, VectorIndex, ScalarValue);
      }
    } else {  //it's a left-hand side multivector
      double Scaling = 1.0;
      Ordering_->GlobalScaling(GlobalRow, Scaling);
      if (GlobalRow % 2 == 0) {
        return Real_->SumIntoGlobalValue(GlobalRow/2, VectorIndex, ScalarValue/Scaling);
      } else {
        return Imag_->SumIntoGlobalValue(GlobalRow/2, VectorIndex, ScalarValue/Scaling);
      }
    }
  }
}

//==========================================================================
int Komplex_MultiVector::ReplaceMyValue(int MyRow,
						    int VectorIndex,
						    double ScalarValue) {
  if (IsOneObject_) {
    if (RHS_) { //it's a right-hand side multivector
      int Index = 0;
      int error = Ordering_->MyIndex(MyRow, Index);
      if (error != 0) { //probably want to keep this error checking in because user could have typed
                        //an index that is too big or negative or not belonging to "me" or something
        return error;
      }
      if (Index == -1) {
        if (MyRow % 2 == 0) {
          return Real_->ReplaceMyValue(MyRow+1, VectorIndex, ScalarValue);
        } else {
          return Real_->ReplaceMyValue(MyRow-1, VectorIndex, ScalarValue);
        }
      } else {
        return Real_->ReplaceMyValue(MyRow, VectorIndex, ScalarValue);
      }
    } else { //it's a left-hand side multivector
      double Scaling = 1.0;
      Ordering_->MyScaling(MyRow, Scaling);
      if (Scaling != 1.0) {
        return Real_->ReplaceMyValue(MyRow, VectorIndex, ScalarValue/Scaling);
      } else {
        return Real_->ReplaceMyValue(MyRow, VectorIndex, ScalarValue);
      }
    }
  } else { //it's two objects
    if (RHS_) {
      int Index = 0;
      int error = Ordering_->MyIndex(MyRow, Index);
      if (error != 0) { //keep, for reasons above =)
        return error;
      }
      if (Index == -1) {
	  return Imag_->ReplaceMyValue(MyRow/2, VectorIndex, ScalarValue);
      } else {
        return Real_->ReplaceMyValue(MyRow/2, VectorIndex, ScalarValue);
      } 
    } else {  //it's a left-hand side multivector
      double Scaling = 1.0;
      Ordering_->MyScaling(MyRow, Scaling);
      if (MyRow % 2 == 0) {
        return Real_->ReplaceMyValue(MyRow/2, VectorIndex, ScalarValue/Scaling);
      } else {
        return Imag_->ReplaceMyValue(MyRow/2, VectorIndex, ScalarValue/Scaling);
      }
    } 
  }
}

//==========================================================================
int Komplex_MultiVector::SumIntoMyValue(int MyRow,
						    int VectorIndex,
						    double ScalarValue) {
  if (IsOneObject_) {
    if (RHS_) { //it's a right-hand side multivector
      int Index = 0;
      int error = Ordering_->MyIndex(MyRow, Index);
      if (error != 0) { //probably want to keep this error checking in because user could have typed
                        //an index that is too big or negative or not belonging to "me" or something
        return error;
      }
      if (Index == -1) {
        if (MyRow % 2 == 0) {
	    return Real_->SumIntoMyValue(MyRow+1, VectorIndex, ScalarValue);
        } else {
	    return Real_->SumIntoMyValue(MyRow-1, VectorIndex, ScalarValue);
        }
      } else {
        return Real_->SumIntoMyValue(MyRow, VectorIndex, ScalarValue);
      }
    } else { //it's a left-hand side multivector
      double Scaling = 1.0;
      Ordering_->MyScaling(MyRow, Scaling);
      if (Scaling != 1.0) {
        return Real_->SumIntoMyValue(MyRow, VectorIndex, ScalarValue/Scaling);
      } else {
        return Real_->SumIntoMyValue(MyRow, VectorIndex, ScalarValue);
      }
    }
  } else { //it's two objects
    if (RHS_) {
      int Index = 0;
      int error = Ordering_->MyIndex(MyRow, Index);
      if (error != 0) { //keep, for reasons above =)
        return error;
      }
      if (Index == -1) {
        return Imag_->SumIntoMyValue(MyRow/2, VectorIndex, ScalarValue);
      } else {
        return Real_->SumIntoMyValue(MyRow/2, VectorIndex, ScalarValue);
      } 
    } else {  //it's a left-hand side multivector
      double Scaling = 1.0;
      Ordering_->MyScaling(MyRow, Scaling);
      if (MyRow % 2 == 0) {
        return Real_->SumIntoMyValue(MyRow/2, VectorIndex, ScalarValue/Scaling);
      } else {
	return Imag_->SumIntoMyValue(MyRow/2, VectorIndex, ScalarValue/Scaling);
      }
    }
  }
}

//==========================================================================
int Komplex_MultiVector::PutScalar(double ScalarConstant) {
  if (IsOneObject_) {
    if (RHS_) {
      return Real_->PutScalar(ScalarConstant);
    } else { //could have negative scalings
      double* Scaling = new double[Real_->GlobalLength()];
      int errorScaling = Ordering_->ScalingVector(Scaling);
      for (int i = 0; i < Real_->GlobalLength(); i++) {
        if (Scaling[i] == -1.0) {
          for (int j = 0; j < Real_->NumVectors(); j++) {
            Real_->ReplaceGlobalValue(i, j, -1.0 * ScalarConstant);
          }
        } else {
          for (int j = 0; j < Real_->NumVectors(); j++) {
            Real_->ReplaceGlobalValue(i, j, ScalarConstant);
          }
        }        
      }
      return errorScaling;
    }
  } else { //it's two objects
    if (RHS_) {
      int ErrorReal = Real_->PutScalar(ScalarConstant);
      int ErrorImag = Imag_->PutScalar(ScalarConstant);
      if (ErrorReal != 0) {
        return ErrorReal;
      } else if (ErrorImag != 0) {
        return ErrorImag;
      } else {
        return(0);
      }
    } else { //it's left-hand side
      double* Scaling = new double[Real_->GlobalLength()];
      int errorScaling = Ordering_->ScalingVector(Scaling);
      for (int i = 0; i < Real_->GlobalLength(); i++) {
        if (Scaling[i] == -1.0) {
          for (int j = 0; j < Real_->NumVectors(); j++) {
            Real_->ReplaceGlobalValue(i, j, -1.0 * ScalarConstant);
            Imag_->ReplaceGlobalValue(i, j, -1.0 * ScalarConstant);
          }
        } else {
          for (int j = 0; j < Real_->NumVectors(); j++) {
            Real_->ReplaceGlobalValue(i, j, ScalarConstant);
            Imag_->ReplaceGlobalValue(i, j, ScalarConstant);
          }
        }        
      }
      return errorScaling;
    }
  }
}

//==========================================================================
int Komplex_MultiVector::Random() {
  if (IsOneObject_) {
    return Real_->Random();
  } else {
    int ErrorReal = Real_->Random();
    int ErrorImag = Imag_->Random();
    if (ErrorReal != 0) {
      return ErrorReal;
    } else if (ErrorImag != 0) {
      return ErrorImag;
    } else {
      return(0);
    }
  }
}

//==========================================================================
void Komplex_MultiVector::CreateOtherMap() {
  if (IsOneObject_) {
    CreateHalfMap();
  } else {
    CreateDoubleMap();
  }
}

//==========================================================================
int Komplex_MultiVector::Dot(const Komplex_MultiVector& A,
				     double* Result) const {
  return EpetraMultiVector()->Dot(*(A.EpetraMultiVector()), Result);
}

//==========================================================================
int Komplex_MultiVector::Abs(const Komplex_MultiVector& A) {
  if (IsOneObject_) {
    return Real_->Abs(*(A.EpetraMultiVector()));
  } else {
    int ErrorReal = Real_->Abs(*(A.RealMultiVector()));
    int ErrorImag = Imag_->Abs(*(A.ImagMultiVector()));
    if (ErrorReal != 0) {
      return ErrorReal;
    } else if (ErrorImag != 0) {
      return ErrorImag;
    } else {
      return(0);
    }
  }
}

//==========================================================================
int Komplex_MultiVector::Reciprocal(const Komplex_MultiVector& A) {
  if (IsOneObject_) {
    return Real_->Reciprocal(*(A.EpetraMultiVector()));
  } else {
    int ErrorReal = Real_->Reciprocal(*(A.RealMultiVector()));
    int ErrorImag = Imag_->Reciprocal(*(A.ImagMultiVector()));
    if (ErrorReal != 0) {
      return ErrorReal;
    } else if (ErrorImag != 0) {
      return ErrorImag;
    } else {
      return(0);
    }
  }
}

//==========================================================================
int Komplex_MultiVector::Scale(double ScalarValue) {
  if (IsOneObject_) {
    return Real_->Scale(ScalarValue);
  } else {
    int ErrorReal = Real_->Scale(ScalarValue);
    int ErrorImag = Imag_->Scale(ScalarValue);
    if (ErrorReal != 0) {
      return ErrorReal;
    } else if (ErrorImag != 0) {
      return ErrorImag;
    } else {
      return(0);
    }
  }
}

//==========================================================================
int Komplex_MultiVector::Scale(double ScalarA,
					 const Komplex_MultiVector& A) {
  if (IsOneObject_) {
    return Real_->Scale(ScalarA, *(A.EpetraMultiVector()));
  } else {
    int ErrorReal = Real_->Scale(ScalarA, *(A.RealMultiVector()));
    int ErrorImag = Imag_->Scale(ScalarA, *(A.ImagMultiVector()));
    if (ErrorReal != 0) {
      return ErrorReal;
    } else if (ErrorImag != 0) {
      return ErrorImag;
    } else {
      return(0);
    }
  }
}

//==========================================================================
int Komplex_MultiVector::Update(double ScalarA,
					  const Komplex_MultiVector& A,
					  double ScalarThis) {
  if (IsOneObject_) {
    return Real_->Update(ScalarA, *(A.EpetraMultiVector()), ScalarThis);
  } else {
    int ErrorReal = Real_->Update(ScalarA, *(A.RealMultiVector()), ScalarThis);
    int ErrorImag = Imag_->Update(ScalarA, *(A.ImagMultiVector()), ScalarThis);
    if (ErrorReal != 0) {
      return ErrorReal;
    } else if (ErrorImag != 0) {
      return ErrorImag;
    } else {
      return(0);
    }
  }
}

//==========================================================================
int Komplex_MultiVector::Update(double ScalarA,
					  const Komplex_MultiVector& A, 
					  double ScalarB,
					  const Komplex_MultiVector& B,
					  double ScalarThis) {
  if (IsOneObject_) {
    return Real_->Update(ScalarA, *(A.EpetraMultiVector()), ScalarB, *(B.EpetraMultiVector()), ScalarThis);
  } else {
    int ErrorReal = Real_->Update(ScalarA, *(A.RealMultiVector()), ScalarB, *(B.RealMultiVector()), ScalarThis);
    int ErrorImag = Imag_->Update(ScalarA, *(A.ImagMultiVector()), ScalarB, *(B.ImagMultiVector()), ScalarThis);
    if (ErrorReal != 0) {
      return ErrorReal;
    } else if (ErrorImag != 0) {
      return ErrorImag;
    } else {
      return(0);
    }
  }
}

//==========================================================================
int Komplex_MultiVector::Norm1(double* Result) const { 
  if (IsOneObject_) {
    for (int i = 0; i < Real_->NumVectors(); i++) {
      double* Values = (*Real_)[i];
      double VecNorm = 0.0;
      for (int j = 0; j < Real_->GlobalLength(); j++) {
        VecNorm += abs(Values[j]);
      }
      Result[i] = VecNorm;
    }
    return(0);
  } else { //it's two objects
    for (int i = 0; i < Real_->NumVectors(); i++) {
      double* Reals = (*Real_)[i];
      double* Imags = (*Imag_)[i];
      double VecNorm = 0.0;
      for (int j = 0; j < Real_->GlobalLength(); j++) {
        VecNorm += abs(Reals[j]) + abs(Imags[j]);
      }
      Result[i] = VecNorm;
    }
  } 
  return(0);
}

//==========================================================================
int Komplex_MultiVector::ComplexNorm1(double* Result) const {
  if (IsOneObject_) {
    for (int i = 0; i < Real_->NumVectors(); i++) {
      double* Values = (*Real_)[i];
      double VecNorm = 0.0;
      for (int j = 0; j < Real_->GlobalLength(); j+=2) {
        double tempFirst = Values[j];
        double tempSecond = Values[j+1];
        double temp = sqrt(tempFirst * tempFirst+tempSecond * tempSecond);
        VecNorm += temp;
      }
      Result[i] = VecNorm;
    }
    return(0);
  } else { //it's two objects
    for (int i = 0; i < Real_->NumVectors(); i++) {
      double* Reals = (*Real_)[i];
      double* Imags = (*Imag_)[i];
      double VecNorm = 0.0;
      for (int j = 0; j < Real_->GlobalLength(); j++) {
        double tempFirst = Reals[j];
        double tempSecond = Imags[j];
        double temp = sqrt(tempFirst * tempFirst+tempSecond * tempSecond);
        VecNorm += temp;
      }
      Result[i] = VecNorm;
    }
    return(0);
  }
}

//==========================================================================
int Komplex_MultiVector::Norm2(double* Result) const {
  if (IsOneObject_) {
    for (int i = 0; i < Real_->NumVectors(); i++) {
      double* Values = (*Real_)[i];
      double VecNorm = 0.0;
      for (int j = 0; j < Real_->GlobalLength(); j++) {
        double temp = Values[j] * Values[j];
        VecNorm += temp;
      }
      Result[i] = sqrt(VecNorm);
    }
    return(0);
  } else {
    for (int i = 0; i < Real_->NumVectors(); i++) {
      double* Reals = (*Real_)[i];
      double* Imags = (*Imag_)[i];
      double VecNorm = 0.0;
      for (int j = 0; j < Real_->GlobalLength(); j++) {
        double tempFirst = Reals[j] * Reals[j];
        double tempSecond = Imags[j] * Imags[j];
        VecNorm += (tempFirst + tempSecond);
      }
      Result[i] = sqrt(VecNorm);
    }
    return(0);
  }
}

//==========================================================================
int Komplex_MultiVector::ComplexNorm2(double* Result) const {
  return Norm2(Result);
}

//==========================================================================
int Komplex_MultiVector::NormInf(double* Result) const {
  if (IsOneObject_) {
    double* Max = new double[Real_->NumVectors()];
    double* Min = new double[Real_->NumVectors()];
    int maxError = Real_->MaxValue(Max);
    if (maxError != 0) {
      return maxError;
    }
    int minError = Real_->MinValue(Min);
    if (minError != 0) {
      return minError;
    }
    for (int i = 0; i < Real_->NumVectors(); i++) {
      Result[i] = (Max[i] >= -1.0 * Min[i]) ? Max[i] : -1.0 * Min[i];
    }
    return(0);
  } else { //two objects
    double* MaxReal = new double[Real_->NumVectors()];
    double* MaxImag = new double[Imag_->NumVectors()];
    double* MinReal = new double[Real_->NumVectors()];
    double* MinImag = new double[Imag_->NumVectors()];
    int error = Real_->MaxValue(MaxReal);
    if (error != 0) {
      return error;
    }
    error = Imag_->MaxValue(MaxImag);
    if (error != 0) {
      return error;
    }
    error = Real_->MinValue(MinReal);
    if (error != 0) {
      return error;
    }
    error = Imag_->MinValue(MinImag);
    if (error != 0) {
      return error;
    }
    for (int i = 0; i < Real_->NumVectors(); i++) {
      double max = (MaxReal[i] >= MaxImag[i]) ? MaxReal[i] : MaxImag[i];
      double min = (MinReal[i] <= MinImag[i]) ? MinReal[i] : MinImag[i];
      Result[i] = (max >= -1.0 * min) ? max : -1.0 * min;
    }
    return(0);
  }
} 

//==========================================================================
int Komplex_MultiVector::ComplexNormInf(double* Result) const {
  if (IsOneObject_) {
    for (int i = 0; i < Real_->NumVectors(); i++) {
      double* Values = (*Real_)[i];
      double CurMax = 0.0;
      for (int j = 0; j < Real_->GlobalLength(); j+=2) {
        double temp = Values[j] * Values[j] + Values[j+1] * Values[j+1];
        if (temp > CurMax) {
          CurMax = temp;
        }
      }
      Result[i] = sqrt(CurMax);
    }
    return(0);
  } else {
    for (int i = 0; i < Real_->NumVectors(); i++) {
      double* Reals = (*Real_)[i];
      double* Imags = (*Imag_)[i];
      double CurMax = 0.0;
      for (int j = 0; j < Real_->GlobalLength(); j++) {
        double temp = Reals[j] * Reals[j] + Imags[j] * Imags[j];
        if (temp > CurMax) {
          CurMax = temp;
        }
      }
      Result[i] = sqrt(CurMax);
    }
  }
  return(0);
}

//==========================================================================
int Komplex_MultiVector::NormWeighted(const Epetra_MultiVector& Weights,
						  double* Result) const {
  if (IsOneObject_) {
    if (RHS_) {
      bool OneW = false;
      if (Weights.NumVectors() == 1) {
        OneW = true; 
      } else if (Real_->NumVectors() != Weights.NumVectors()) {
        return -1;
      } else if (Real_->MyLength() != Weights.MyLength()) {
        return -2;
      }
      double* W = Weights.Values();
      double** W_Pointers = Weights.Pointers();
      int* Perms = new int[Real_->GlobalLength()];
      Ordering_->PermutationVector(Perms);
      for (int i = 0; i < Real_->NumVectors(); i++) {
        double* Values = (*Real_)[i];
        if (!OneW) {
          W = W_Pointers[i];
        }
        double CurNorm = 0.0;
        for (int j = 0; j < Real_->GlobalLength(); j+=2) {
          double temp1 = 0.0;
          double temp2 = 0.0;
          if (Perms[j] == -1) {
            temp1 = Values[j+1] / W[j];
            temp2 = Values[j] / W[j+1];
          } else {
            temp1 = Values[j] / W[j];
            temp2 = Values[j+1] / W[j+1];
          }
          CurNorm += temp1 * temp1;
          CurNorm += temp2 * temp2;
        }
        CurNorm /= Real_->NumVectors();
        Result[i] = sqrt(CurNorm);
      }
      return(0);
    } else {
      return Real_->NormWeighted(Weights, Result); 
    }
  } else {
    bool OneW = false;
    if (Weights.NumVectors() == 1) {
      OneW = true;
    } else if (Real_->NumVectors() != Weights.NumVectors()) {
      return -1;
    } else if (Real_->MyLength() != Weights.MyLength() || Imag_->MyLength() != Weights.MyLength()) {
      return -2;
    }
    double* W = Weights.Values();
    double** W_Pointers = Weights.Pointers();
    if (RHS_) {
      int* Perms = new int[Real_->GlobalLength()];
      Ordering_->PermutationVector(Perms);
      for (int i = 0; i < Real_->NumVectors(); i++) {
        double* Reals = (*Real_)[i];
        double* Imags = (*Imag_)[i];
        if (!OneW) {
          W = W_Pointers[i];
        }
        double CurNorm = 0.0;
        for (int j = 0; j < Real_->GlobalLength(); j++) {
          double temp1 = 0.0;
          double temp2 = 0.0;
          if (Perms[j] == -1) {
            temp1 = Imags[j] / W[2*j];
            temp2 = Reals[j+1] / W[2*j+1];
          } else {
            temp1 = Reals[j] / W[2*j];
            temp2 = Imags[j+1] / W[2*j+1];
          }
          CurNorm += temp1 * temp1;
          CurNorm += temp2 * temp2;
        }
        CurNorm /= Real_->NumVectors();
        Result[i] = sqrt(CurNorm);
      }
      return(0);
    } else {
      for (int i = 0; i < Real_->NumVectors(); i++) {
        double* Reals = (*Real_)[i];
        double* Imags = (*Imag_)[i];
        if (!OneW) {
          W = W_Pointers[i];
        } 
        double CurNorm = 0.0;
        for (int j = 0; j < Real_->GlobalLength(); j++) {
          double temp = Reals[j] / W[2*j];
          CurNorm += temp * temp;
          temp = Imags[j] / W[2*j+1];
          CurNorm += temp * temp;
        }
        CurNorm /= Real_->NumVectors();
        Result[i] = sqrt(CurNorm);
      }
      return(0);
    }
  }
}

//==========================================================================
int Komplex_MultiVector::MinValue(double* Result) const {
  if (IsOneObject_) {
    if (RHS_) {
      return Real_->MinValue(Result);
    } else {
      double* Scales = new double[Real_->GlobalLength()];
      int errorScale = Ordering_->ScalingVector(Scales);
      if (errorScale != 0) {
        return errorScale;
      }
      for (int i = 0; i < Real_->NumVectors(); i++) {
        double* Values = (*Real_)[i];
        double CurMin = Scales[0] * Values[0];
        for (int j = 0; j < Real_->GlobalLength(); j++) {
          double temp = Values[i] * Scales[i];
          if (temp < CurMin) {
            CurMin = temp;
          }
        }
        Result[i] = CurMin;
      }
      delete Scales;
      return(0);
    }
  } else {
    if (RHS_) {
      double* MinReal = new double[Real_->NumVectors()];
      double* MinImag = new double[Imag_->NumVectors()];
      int error = Real_->MinValue(MinReal);
      if (error != 0) {
        return error;
      }
      error = Imag_->MinValue(MinImag);
      if (error != 0) {
        return error;
      }
      for (int i = 0; i < Real_->NumVectors(); i++) {
        Result[i] = (MinReal[i] <= MinImag[i]) ? MinReal[i] : MinImag[i];
      }
      return(0);
    } else {
      double* Scales = new double[Real_->GlobalLength()];
      int errorScale = Ordering_->ScalingVector(Scales);
      if (errorScale != 0) {
        return errorScale;
      }
      for (int i = 0; i < Real_->NumVectors(); i++) {
        double* Reals = (*Real_)[i];
        double* Imags = (*Imag_)[i];
        double CurMin = Scales[0] * Reals[0];
        for (int j = 0; j < Real_->GlobalLength(); j++) {
          double tempReal = Reals[i] * Scales[i];
          double tempImag = Imags[i] * Scales[i];
          double temp = (tempReal <= tempImag) ? tempReal : tempImag;
          if (temp < CurMin) {
            CurMin = temp;
          }
        }
        Result[i] = CurMin;
      }
      delete [] Scales;
      return(0);
    }
  }
}

//==========================================================================
int Komplex_MultiVector::MaxValue(double* Result) const {
  if (IsOneObject_) {
    if (RHS_) {
      return Real_->MaxValue(Result);
    } else {
      double* Scales = new double[Real_->GlobalLength()];
      int errorScale = Ordering_->ScalingVector(Scales);
      if (errorScale != 0) {
        return errorScale;
      }
      for (int i = 0; i < Real_->NumVectors(); i++) {
        double* Values = (*Real_)[i];
        double CurMax = Values[0] * Scales[0];
        for (int j = 1; j < Real_->GlobalLength(); j++) {
          double temp = Scales[j] * Values[j];
          if (temp > CurMax) {
            CurMax = temp;
          }
        }
        Result[i] = CurMax;
      }
      delete [] Scales;
      return(0);
    }
  } else {
    if (RHS_) {
      double* MaxReal = new double[Real_->NumVectors()];
      double* MaxImag = new double[Imag_->NumVectors()];
      int error = Real_->MaxValue(MaxReal);
      if (error != 0) {
        return error;
      }
      error = Imag_->MaxValue(MaxImag);
      if (error != 0) {
        return error;
      }
      for (int i = 0; i < Real_->NumVectors(); i++) {
        Result[i] = (MaxReal[i] >= MaxImag[i]) ? MaxReal[i] : MaxImag[i];
      }
      delete [] MaxReal;
      delete [] MaxImag;
      return(0);
    } else {
      double* Scales = new double[Real_->GlobalLength()];
      int errorScale = Ordering_->ScalingVector(Scales);
      if (errorScale != 0) {
        return errorScale;
      }
      for (int i = 0; i < Real_->NumVectors(); i++) {
        double* Reals = (*Real_)[i];
        double* Imags = (*Imag_)[i];
        double CurMax = Reals[0] * Scales[0];
        for (int j = 0; j < Real_->GlobalLength(); j++) {
          double tempReal = Scales[i] * Reals[i];
          double tempImag = Scales[i] * Imags[i];
          double temp = (tempReal >= tempImag) ? tempReal : tempImag;
          if (temp > CurMax) {
            CurMax = temp;
          }
        }
        Result[i] = CurMax;
      }
      delete [] Scales;
      return(0);
    }  
  }
}

//==========================================================================
int Komplex_MultiVector::MeanValue(double* Result) const {
  if (IsOneObject_) {
    if (RHS_) {
      return Real_->MeanValue(Result);
    } else { //left-hand side
      double* Scales = new double[Real_->GlobalLength()];
      int errorScale = Ordering_->ScalingVector(Scales);
      if (errorScale != 0) {
        return errorScale;
      }
      for (int i = 0; i < Real_->NumVectors(); i++) {
        double* Values = (*Real_)[i];
        double Total = 0.0;
        for (int j = 0; j < Real_->GlobalLength(); j++) {
          Total += Values[j] * Scales[j];
        }
        Result[i] = Total / (Real_->GlobalLength());
      }
      delete [] Scales;
      return(0);
    }
  } else { //two multivectors
    if (RHS_) {
      double* MeanReal = new double[Real_->NumVectors()];
      double* MeanImag = new double[Imag_->NumVectors()];
      int error = Real_->MeanValue(MeanReal);
      if (error != 0) {
        return error;
      }
      error = Imag_->MeanValue(MeanImag);
      if (error != 0) {
        return error;
      }
      for (int i = 0; i < Real_->NumVectors(); i++) {
        Result[i] = (MeanReal[i] + MeanImag[i]) / 2.0;
      }
      return(0);
    } else { //left-hand side
      double* Scales = new double[Real_->GlobalLength()];
      int errorScale = Ordering_->ScalingVector(Scales);
      if (errorScale != 0) {
        return errorScale;
      }
      for (int i = 0; i < Real_->NumVectors(); i++) {
        double* Reals = (*Real_)[i];
        double* Imags = (*Imag_)[i];
        double Total = 0.0;
        for (int j = 0; j < Real_->GlobalLength(); j++) {
          Total += Scales[j] * (Reals[j] + Imags[j]);
        }
        Result[i] = Total / (Real_->GlobalLength());
      }
      delete [] Scales;
      return(0);
    }
  }
}

//==========================================================================
int Komplex_MultiVector::SetSeed(unsigned int Seed) {
  if (IsOneObject_) {
    return Real_->SetSeed(Seed);
  } else {
    int ErrorReal = Real_->SetSeed(Seed);
    int ErrorImag = Imag_->SetSeed(Seed);
    if (ErrorReal != 0) {
      return ErrorReal;
    } else if (ErrorImag != 0) {
      return ErrorImag;
    } else {
      return(0);
    }
  }
}

//==========================================================================
unsigned int Komplex_MultiVector::Seed()const {
  return Real_->Seed(); //because always set both Real_ and Imag_
}

//==========================================================================
Komplex_MultiVector & Komplex_MultiVector::operator=(const Komplex_MultiVector& Source) { 
  //not checking to see if the Maps are the same or not!!!!!!!!!!!!...........
  if (this != &Source) {
    if (NumVectors() != Source.NumVectors()) {
      std::cerr <<"Number of vectors incompatible.  The this MultiVector has NumVectors = " << NumVectors()
		<< ".  The other MultiVector has NumVectors = " << Source.NumVectors() << endl;
    }
    if (MyLength() != Source.MyLength()) {
      std::cerr << "Length of MultiVectors incompatible.  The this MultiVector has MyLength = " << MyLength()
		<< ".  The other MultiVector has MyLength = " + Source.MyLength() << endl;
    }
    Real_ = Source.Real_;
    Ordering_ = Source.Ordering_;
    IsOneObject_ = Source.IsOneObject_;
    if (!IsOneObject_) {
      Imag_ = Source.Imag_;
    } else { //it is one object; check to see if Imag_ had been initialized in this
      if (Imag_ != 0) {
        Imag_ = 0;
      }
      OtherMap_ = Source.OtherMap_;
      RHS_ = Source.RHS_;
      return (*this);
    }
  }
}

//==========================================================================
double * & Komplex_MultiVector::operator[](int i) {
  if (IsOneObject_) {
    if (RHS_) {
      TempDouble_ = (*Real_)[i];
      int* Perms = new int[Real_->GlobalLength()];
      Ordering_->PermutationVector(Perms);
      for (int j = 0; j < Real_->GlobalLength(); j+=2) {
        if (Perms[j] == -1) {
          double temp = TempDouble_[j];
          TempDouble_[j] = TempDouble_[j+1];
          TempDouble_[j+1] = temp;
        }
      }
      delete [] Perms;
      return TempDouble_;
    } else { //left-hand side
      TempDouble_ = (*Real_)[i];
      double* Scaling = new double[Real_->GlobalLength()];
      Ordering_->ScalingVector(Scaling);
      for (int j = 0; j < Real_->GlobalLength(); j++) {
        TempDouble_[j] *= Scaling[j];
      }
      delete [] Scaling;
      return TempDouble_;
    } 
  } else { //two objects
    if (RHS_) {
      TempDouble_ = new double[Real_->GlobalLength()*2];
      double* Reals = (*Real_)[i];
      double* Imags = (*Imag_)[i];
      int* Perms = new int[Real_->GlobalLength()];
      Ordering_->PermutationVector(Perms);
      for (int j = 0; j < Real_->GlobalLength(); j++) {
        if (Perms[j] == -1) {
	    TempDouble_[2*j] = Imags[j];
	    TempDouble_[2*j+1] = Reals[j];
	  } else {
          TempDouble_[2*j] = Reals[j];
	    TempDouble_[2*j+1] = Imags[j];
        }
      }
      delete [] Perms;
      return TempDouble_;
    } else { //left-hand side
      TempDouble_ = new double[Real_->GlobalLength()*2];
      double* Reals = (*Real_)[i];
      double* Imags = (*Imag_)[i];
      double* Scaling = new double[Real_->GlobalLength()];
      Ordering_->ScalingVector(Scaling);
      for (int j = 0; j < Real_->GlobalLength(); j++) {
        TempDouble_[j*2] = Reals[j] * Scaling[j];
        TempDouble_[j*2+1] = Imags[j] * Scaling[j];
      }
      delete [] Scaling;
      return TempDouble_;
    } 
  }
}

//==========================================================================
double * const & Komplex_MultiVector::operator[](int i) const {
  if (IsOneObject_) {
    if (RHS_) {
      TempDouble_ = (*Real_)[i];
      int* Perms = new int[Real_->GlobalLength()];
      Ordering_->PermutationVector(Perms);
      for (int j = 0; j < Real_->GlobalLength(); j+=2) {
        if (Perms[j] == -1) {
	    double temp = TempDouble_[j];
	    TempDouble_[j] = TempDouble_[j+1];
	    TempDouble_[j+1] = temp;
        }
      }
      delete [] Perms;
      return TempDouble_;
    } else { //left-hand side
      TempDouble_ = (*Real_)[i];
      double* Scaling = new double[Real_->GlobalLength()];
      Ordering_->ScalingVector(Scaling);
      for (int j = 0; j < Real_->GlobalLength(); j++) {
        TempDouble_[j] *= Scaling[j];
      }
      delete [] Scaling;
      return TempDouble_;
    } 
  } else { //two objects
    if (RHS_) {
      TempDouble_ = new double[Real_->GlobalLength()*2];
      double* Reals = (*Real_)[i];
      double* Imags = (*Imag_)[i];
      int *Perms = new int[Real_->GlobalLength()];
      Ordering_->PermutationVector(Perms);
      for (int j = 0; j < Real_->GlobalLength(); j++) {
        if (Perms[j] == -1) {
	    TempDouble_[2*j] = Imags[j];
	    TempDouble_[2*j+1] = Reals[j];
        } else {
          TempDouble_[2*j] = Reals[j];
	    TempDouble_[2*j+1] = Imags[j];
        }
      }
      delete [] Perms;
      return TempDouble_;
    } else { //left-hand side
      TempDouble_ = new double[Real_->GlobalLength() * 2];
      double* Reals = (*Real_)[i];
      double* Imags = (*Imag_)[i];
      double* Scaling = new double[Real_->GlobalLength()];
      Ordering_->ScalingVector(Scaling);
      for (int j = 0; j < Real_->GlobalLength(); j++) {
        TempDouble_[j*2] = Reals[j] * Scaling[j];
        TempDouble_[j*2+1] = Imags[j] * Scaling[j];
      }
      delete [] Scaling;
      return TempDouble_;
    } 
  }
}

//==========================================================================
/*Komplex_Vector * & Komplex_MultiVector::operator()(int i) {
  if (IsOneObject_) {
    something
  } else {
    something different
  }
  return null;
  }*/ //could I just extract the values and shoot them into a K_V using previous method?

//==========================================================================
/*const Komplex_Vector * & Komplex_MultiVector::operator()(int i) const {
  if (IsOneObject_) {
    something
  } else {
    something different
  }
  return null;
  }*/

//==========================================================================
Epetra_MultiVector * Komplex_MultiVector::EpetraMultiVector() const {
  double** ArrayOfPointers = new double*[Real_->NumVectors()];
  for (int i = 0; i < Real_->NumVectors(); i++) {
    ArrayOfPointers[i] = (*this)[i];
  }
  Epetra_MultiVector* Ret;
  if (IsOneObject_) {
    Ret = new Epetra_MultiVector(Copy, Real_->Map(), ArrayOfPointers, Real_->NumVectors());
  } else {
    Ret = new Epetra_MultiVector(Copy, *OtherMap_, ArrayOfPointers, Real_->NumVectors());
  }
  return Ret;
}

//==========================================================================
Epetra_MultiVector * Komplex_MultiVector::RealMultiVector() const {
  double** AOPReal = new double*[Real_->NumVectors()];
  for (int i = 0; i < Real_->NumVectors(); i++) {
    AOPReal[i] = RealValues(i);
  }
  Epetra_MultiVector* Ret;
  if (IsOneObject_) {
    Ret = new Epetra_MultiVector(Copy, *OtherMap_, AOPReal, Real_->NumVectors());
  } else {
    Ret = new Epetra_MultiVector(Copy, Real_->Map(), AOPReal, Real_->NumVectors());
  }
  return Ret;
}

//==========================================================================
Epetra_MultiVector * Komplex_MultiVector::ImagMultiVector() const {
  double** AOPImag = new double*[Real_->NumVectors()];
  for (int i = 0; i < Real_->NumVectors(); i++) {
    AOPImag[i] = ImagValues(i);
  }
  Epetra_MultiVector* Ret;
  if (IsOneObject_) {
    Ret = new Epetra_MultiVector(Copy, *OtherMap_, AOPImag, Real_->NumVectors());
  } else {
    Ret = new Epetra_MultiVector(Copy, Imag_->Map(), AOPImag, Imag_->NumVectors());
  }
  return Ret;
}

//==========================================================================
Epetra_Vector * Komplex_MultiVector::EpetraVector(int index) const {
  if (IsOneObject_) {
    double* Values = (*this)[index];
    Epetra_Vector* Ret = new Epetra_Vector(Copy, Real_->Map(), Values);
    return Ret;
  } else {
    double* Values = (*this)[index];
    Epetra_Vector* Ret = new Epetra_Vector(Copy, *OtherMap_, Values);
    return Ret;
  }
}

//==========================================================================
Epetra_Vector * Komplex_MultiVector::RealVector(int index) const {
  if (IsOneObject_) {
    double* Reals = RealValues(index);
    Epetra_Vector* Ret = new Epetra_Vector(Copy, *OtherMap_, Reals);
    return Ret;
  } else {
    double* Reals = RealValues(index);
    Epetra_Vector* Ret = new Epetra_Vector(Copy, Real_->Map(), Reals);
    return Ret;
  }
}

//==========================================================================
Epetra_Vector * Komplex_MultiVector::ImagVector(int index) const {
  if (IsOneObject_) {
    double* Imags = ImagValues(index);
    Epetra_Vector* Ret = new Epetra_Vector(Copy, *OtherMap_, Imags);
    return Ret;
  } else {
    double* Imags = ImagValues(index);
    Epetra_Vector* Ret = new Epetra_Vector(Copy, Imag_->Map(), Imags);
    return Ret;
  }
}

//==========================================================================
double * & Komplex_MultiVector::RealValues(int i) const {
  if (IsOneObject_) {
    if (RHS_) {
      double* Reals = new double[Real_->GlobalLength()/2];
      double* Values = (*Real_)[i];
      int* Perms = new int[Real_->GlobalLength()];
      Ordering_->PermutationVector(Perms);
      for (int j = 0; j < Real_->GlobalLength(); j+=2) {
        if (Perms[j] == -1) {
          Reals[j/2] = Values[j+1];
        } else {
          Reals[j/2] = Values[j];
        }
      }
      delete [] Perms;
      return Reals;
    } else { //left-hand side multivector
      double* Reals = new double[Real_->GlobalLength()/2];
      double* Values = (*Real_)[i];
      double* Scaling = new double[Real_->GlobalLength()];
      Ordering_->ScalingVector(Scaling);
      for (int j = 0; j < Real_->GlobalLength(); j+=2) {
        Reals[j] = Scaling[j] * Values[j];
      }
      delete [] Scaling;
      return Reals;
    }
  } else { //two objects
    if (RHS_) {
      return (*Real_)[i];
    } else { //left-hand side multivector
      double* Reals = (*Real_)[i];
      double* Scaling = new double[Real_->GlobalLength()];
      Ordering_->ScalingVector(Scaling);
      for (int j = 0; j < Real_->GlobalLength(); j++) {
        Reals[j] *= Scaling[j];
      }
      delete [] Scaling;
      return Reals;
    }
  }
}

//==========================================================================
double * & Komplex_MultiVector::ImagValues(int i) const {
  if (IsOneObject_) {
    if (RHS_) {
      double* Imags = new double[Real_->GlobalLength() / 2];
      double* Values = (*Real_)[i];
      int* Perms = new int[Real_->GlobalLength()];
      Ordering_->PermutationVector(Perms);
      for (int j = 1; j < Real_->GlobalLength(); j+=2) {
        if (Perms[j] == -1) {
          Imags[j/2] = Values[j-1];
        } else {
          Imags[j/2] = Values[j];
        }
      }
      delete [] Perms;
      return Imags;
    } else { //left-hand side multivector
      double* Imags = new double[Real_->GlobalLength()/2];
      double* Values = (*Real_)[i];
      double* Scaling = new double[Real_->GlobalLength()];
      Ordering_->ScalingVector(Scaling);
      for (int j = 0; j < Real_->GlobalLength(); j+=2) {
        Imags[j] = Scaling[j] * Values[j];
      }
      delete [] Scaling;
      return Imags;
    }
  } else { //two objects
    if (RHS_) {
      return (*Imag_)[i];
    } else { //left-hand side multivector
      double* Imags = (*Imag_)[i];
      double* Scaling = new double[Imag_->GlobalLength()];
      Ordering_->ScalingVector(Scaling);
      for (int j = 0; j < Imag_->GlobalLength(); j++) {
        Imags[j] *= Scaling[j];
      }
      delete [] Scaling;
      return Imags;
    }
  }
}

//==========================================================================
int Komplex_MultiVector::NumVectors() const {
  return Real_->NumVectors();
}

//==========================================================================
int Komplex_MultiVector::MyLength() const {
  if (IsOneObject_) {
    return Real_->MyLength();
  } else {
    return Real_->MyLength() + Imag_->MyLength(); 
  }
}

//==========================================================================
int Komplex_MultiVector::GlobalLength() const { 
  if (IsOneObject_) {
    return Real_->GlobalLength();
  } else {
    return Real_->GlobalLength() + Imag_->GlobalLength(); 
  }
}

//==========================================================================
Komplex_KForms Komplex_MultiVector::KForm() const {
  return Ordering_->KForm();
}

//==========================================================================
int Komplex_MultiVector::SwitchKForm(Komplex_KForms NewKForm) {
  return Ordering_->SwitchKForm(NewKForm);
}

//==========================================================================
int Komplex_MultiVector::ReplaceMap(const Epetra_BlockMap& map) {
  if (IsOneObject_) {
    return Real_->ReplaceMap(map);
  } else {
    int ErrorReal = Real_->ReplaceMap(map);
    int ErrorImag = Imag_->ReplaceMap(map);
    if (ErrorReal != 0) {
      return ErrorReal;
    } else if (ErrorImag != 0) {
      return ErrorImag;
    } else {
      return(0);
    }
  }
}

//==========================================================================
void Komplex_MultiVector::Print(ostream& os) const {
  int MyPID = Real_->Map().Comm().MyPID();
  int NumProc = Real_->Map().Comm().NumProc();
  double** Values = new double*[Real_->NumVectors()];
  for (int i = 0; i < Real_->NumVectors(); i++) {
    Values[i] = (*this)[i];
  }  
  int* Perms = new int[Real_->Map().NumGlobalElements()];
  Ordering_->PermutationVector(Perms);
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      os << "Currently in ";
      switch (KForm()) {
        case K1: os << "K1";
	           break;
        case K2: os << "K2";
	 	     break;
        case K3: os << "K3";
	  	     break;
        case K4: os << "K4";
	  	     break;
      } //end switch KForm()
      os << " form" << endl;
      int NumVectors1 = NumVectors();
      int NumMyElements1 = 0;
      int NumMyElementsR = 0;
      int NumMyElementsI = 0;
      int* MyGlobalElements1;
      int* MyGlobalElementsR;
      int* MyGlobalElementsI;
      if (IsOneObject_) {
	  NumMyElements1 = Real_->Map().NumMyElements();
        MyGlobalElements1 = Real_->Map().MyGlobalElements();
      } else {
	  NumMyElementsR = Real_->Map().NumMyElements();
	  NumMyElementsI = Imag_->Map().NumMyElements();
	  MyGlobalElementsR = Real_->Map().MyGlobalElements();
	  MyGlobalElementsI = Imag_->Map().MyGlobalElements();
      }
      int MaxElementSize1 = Real_->Map().MaxElementSize();
      int* FirstPointInElementList1;
      if (MaxElementSize1!=1) {
	  FirstPointInElementList1 = Real_->Map().FirstPointInElementList();
      }      
      if (MyPID==0) {
	  os.width(8);
	  os <<  "     MyPID"; os << "    ";
	  os.width(12);
	  if (MaxElementSize1==1) {
	    os <<  "GID  ";
	  } else {
	    os <<  "     GID/Point";
	  }
	  for (int j = 0; j < NumVectors1 ; j++) {   
	    os.width(20);
	    os <<  "Value  ";
	  }
	  os << endl;
      } //end MyPID==0
      if (IsOneObject_) {
        for (int i=0; i < NumMyElements1; i++) {
	    for (int ii=0; ii < Real_->Map().ElementSize(i); ii++) {
	      int iii;
	      os.width(10);
	      os <<  MyPID; os << "    ";
	      os.width(10);
	      if (MaxElementSize1==1) {
	        if (RHS_) {
		    if (Perms[i] == -1) {
		      if (i % 2 == 0) {
		        os << MyGlobalElements1[i+1] << "    ";
		      } else {	
		        os << MyGlobalElements1[i-1] << "    ";
		      }
		    } else { //Perms[i] != -1
		    os << MyGlobalElements1[i] << "    ";
		  }
		  iii = i;
	      } else { // !RHS_
		  os << MyGlobalElements1[i] << "    ";
		  iii = i;
	      } 
	    } else { //MaxElementSize1 != 1
	      os <<  MyGlobalElements1[i]<< "/" << ii << "    ";
	      iii = FirstPointInElementList1[i]+ii;
	    }
	    for (int j = 0; j < NumVectors1 ; j++) {   
	      os.width(20);
	      os << Values[j][iii];
	    }
	    os << endl;
	  } //for ii
      } //for i
    } else { //is two objects
        for (int i = 0; i < NumMyElementsR; i++) {
	    for (int ii = 0; ii < Real_->Map().ElementSize(i); ii++) {
	      int iii;
	      os.width(10);
	      os << MyPID;
	      os << "    ";
	      os.width(10);
	      if (MaxElementSize1==1) {
	        if (RHS_) {
		    if (Perms[i] == -1) {
		      os << 2 * MyGlobalElementsI[i] << "    ";
		      for (int j = 0; j < NumVectors1; j++) {
		        os.width(20);
		        os << Values[j][2*i];
		      }
		      os << endl;
		      os.width(10);
		      os << MyPID << "    ";
		      os.width(10);
		      os << 2 * MyGlobalElementsR[i]+1 << "    ";
		      for (int j = 0; j < NumVectors1; j++) {
		        os.width(20);
		        os << Values[j][2*i+1];
		      }
		      os << endl;
		    } else { //Perms[i] != -1
		      os << 2 * MyGlobalElementsR[i] << "    ";
		      for (int j = 0; j < NumVectors1; j++) {
		        os.width(20);
		        os << Values[j][2*i];
		      }
		      os << endl;
		      os.width(10);
		      os << MyPID << "    ";
		      os.width(10);
		      os << 2 * MyGlobalElementsI[i]+1 << "    ";
		      for (int j = 0; j < NumVectors1; j++) {
		        os.width(20);
		        os << Values[j][2*i+1];
		      }
		      os << endl;
		    }
	        } else { //not RHS
		    os << 2 * MyGlobalElementsR[i] << "    ";
		    for (int j = 0; j < NumVectors1; j++) {
		      os.width(20);
		      os << Values[j][i];
		    }
		    os << endl;
		    os << 2 * MyGlobalElementsI[i]+1 << "    ";
		    for (int j = 0; j < NumVectors1; j++) {
		      os.width(20);
		      os << Values[j][i];
		    }
		    os << endl;
	        }
	      } else { //MaxElementSize1 != 1, probably need to fill this in! =)
	      }
	    }
	  }	
	  os << flush; 
      }
      
      // Do a few global ops to give I/O a chance to complete
      Real_->Map().Comm().Barrier();
      Real_->Map().Comm().Barrier();
      Real_->Map().Comm().Barrier();
    }
    return; 
  }
}

//==========================================================================
void Komplex_MultiVector::CreateHalfMap() { 
  int NumMyElements = Real_->Map().NumMyElements();
  int* MyGlobalElements = new int[NumMyElements];
  int error = Real_->Map().MyGlobalElements(MyGlobalElements);
  int* NewMyGlobals = new int[(int) NumMyElements/2];
  for (int i = 0; i < NumMyElements / 2; i++) {
    NewMyGlobals[i] = MyGlobalElements[i];
  }
  OtherMap_ = new Epetra_BlockMap((Real_->Map().NumGlobalElements()) / 2, NumMyElements / 2, NewMyGlobals,
                                   Real_->Map().ElementSize(), Real_->Map().IndexBase(), Real_->Map().Comm());
}

//==========================================================================
void Komplex_MultiVector::CreateDoubleMap() {
  int NumMyElements = Real_->Map().NumMyElements();
  int* MyGlobalElements = new int[NumMyElements];
  int error = Real_->Map().MyGlobalElements(MyGlobalElements);
  int* NewMyGlobals = new int[2 * NumMyElements];
  for (int i = 0; i < NumMyElements; i++) {
    NewMyGlobals[i] = MyGlobalElements[i];
    NewMyGlobals[i+NumMyElements] = NumMyElements + MyGlobalElements[i];
  }
  OtherMap_ = new Epetra_BlockMap(2 * Real_->Map().NumGlobalElements(), 2 * NumMyElements, 
                                  NewMyGlobals, Real_->Map().ElementSize(), Real_->Map().IndexBase(), Real_->Map().Comm());
}

void Komplex_MultiVector::PrintReal(ostream& os) {
  os << (*Real_) << endl;
}

void Komplex_MultiVector::PrintImag(ostream& os) {
  os << (*Imag_) << endl;
}

