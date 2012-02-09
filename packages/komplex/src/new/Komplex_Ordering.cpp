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

#include "Komplex_KForms.hpp"
#include "Komplex_Ordering.hpp"
#include "Komplex_MultiVector.hpp"
#include "Epetra_BlockMap.h"

//==========================================================================
Komplex_Ordering::Komplex_Ordering(const Epetra_BlockMap& Map,
					     Komplex_KForms KForm, 
					     bool IsOneObject)
  : P_(Map),
    D_(Map),
    KForm_(KForm),
    StoredKForm_(K1),
    IsOneObject_(IsOneObject)
{
  if (KForm_ == K1 || KForm_ == K2) {
    P_.PutScalar(1.0);
  } else {
    P_.PutScalar(-1.0);
  }
  if (KForm_ == K1 || KForm_ == K3) {
    D_.PutScalar(1.0);
  } else {
    for (int i = 0; i < Map.NumGlobalElements(); i+=2) { 
      D_.ReplaceGlobalValue(i, 0, 1.0);
      D_.ReplaceGlobalValue(i+1, 0, -1.0); 
    }
  }
} 

//==========================================================================
Komplex_Ordering::Komplex_Ordering(Komplex_Ordering& Source)
  : P_(Source.P_),
    D_(Source.D_),
    KForm_(Source.KForm_),
    StoredKForm_(Source.StoredKForm_),
    IsOneObject_(Source.IsOneObject_)
{
}

//==========================================================================
Komplex_Ordering::~Komplex_Ordering(void) {
}

//==========================================================================
Komplex_KForms Komplex_Ordering::KForm(void){
  return KForm_;
}

//==========================================================================
int Komplex_Ordering::SwitchKForm(Komplex_KForms NewKForm) {
  if (KForm_ == NewKForm) {
    return(0);
  }
  int error = 0;  
  if (!((StoredKForm_==K1 && NewKForm==K3) || (StoredKForm_==K2 && NewKForm==K4) || 
        (StoredKForm_==K3 && NewKForm==K1) || (StoredKForm_==K4 && NewKForm==K2))) { //both even or both odd
    for (int i = 1; i < P_.GlobalLength(); i+=2) {
      D_.ReplaceGlobalValue(i, 0, -1.0);
    }
  }
  if ((StoredKForm_==K1 && NewKForm==K3) || (StoredKForm_==K1 && NewKForm==K4) || 
      (StoredKForm_==K2 && NewKForm==K3) || (StoredKForm_==K2 && NewKForm==K4)) {
    P_.PutScalar(-1.0);
  } else if ((StoredKForm_==K3 && NewKForm==K1) || (StoredKForm_==K3 && NewKForm==K2) ||
             (StoredKForm_==K4 && NewKForm==K1) || (StoredKForm_==K4 && NewKForm==K2)) {
    P_.PutScalar(1.0);
  }
  KForm_ = NewKForm;
  return error; 
}

//==========================================================================
int Komplex_Ordering::PermutationVector(int* Perms) { 
  double* PermDoubles = new double[P_.Map().NumGlobalElements()];
  int error = P_.ExtractCopy(PermDoubles);
  for (int i = 0; i < P_.Map().NumGlobalElements(); i++) {
    Perms[i] = (int) PermDoubles[i];
  }
  return error;
}

//==========================================================================
int Komplex_Ordering::ScalingVector(double* Scales) {
  return D_.ExtractCopy(Scales);
}

//==========================================================================
int Komplex_Ordering::GlobalIndex(int GlobalRow, int& Index) { 
  int error = 0;
  if (!P_.Map().MyGID(GlobalRow)) {
    error = 1;
  }
  if (GlobalRow >= P_.Map().NumGlobalElements()) {
    error = -1;
  }
  Index = (int) P_[GlobalRow];
  return error;
}

//==========================================================================
int Komplex_Ordering::GlobalScaling(int GlobalRow, double& Scalar) {
  int error = 0;
  if (!D_.Map().MyGID(GlobalRow)) {
    error = 1;
  }
  if (GlobalRow >= D_.Map().NumGlobalElements()) {
    error = -1;
  }
  Scalar = D_[GlobalRow];
  return error;
}

//==========================================================================
int Komplex_Ordering::MyIndex(int MyRow, int& Index) {
  int GID = P_.Map().GID(MyRow);
  return GlobalIndex(GID, Index);
}

//==========================================================================
int Komplex_Ordering::MyScaling(int MyRow, double& Scalar) {
  int GID = D_.Map().GID(MyRow);
  return GlobalScaling(GID, Scalar);
}

//==========================================================================
void Komplex_Ordering::Reset(Komplex_KForms NewKForm) { 
  P_.PutScalar(1.0);
  D_.PutScalar(1.0);
  KForm_ = NewKForm;
  StoredKForm_ = NewKForm;
}
