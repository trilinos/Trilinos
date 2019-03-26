/*
//@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
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
*/

#include "Ifpack_IHSS.h"
#include "Ifpack.h"
#include "Ifpack_Utils.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"



using Teuchos::RefCountPtr;
using Teuchos::rcp;


#ifdef HAVE_IFPACK_EPETRAEXT
#include "EpetraExt_MatrixMatrix.h"


Ifpack_IHSS::Ifpack_IHSS(Epetra_RowMatrix* A):
  IsInitialized_(false),
  IsComputed_(false),
  Label_(),
  EigMaxIters_(10),
  EigRatio_(30.0),
  LambdaMax_(-1.0),
  Alpha_(-1.0),
  NumSweeps_(1),

  Time_(A->Comm())
{
  Epetra_CrsMatrix *Acrs=dynamic_cast<Epetra_CrsMatrix*>(A);
  TEUCHOS_TEST_FOR_EXCEPT(!Acrs)
  A_=rcp(Acrs,false);
}

void Ifpack_IHSS::Destroy(){
}



int Ifpack_IHSS::Initialize(){
  EigMaxIters_          = List_.get("ihss: eigenvalue max iterations",EigMaxIters_);
  EigRatio_             = List_.get("ihss: ratio eigenvalue", EigRatio_);
  NumSweeps_            = List_.get("ihss: sweeps",NumSweeps_);

  // Counters
  IsInitialized_=true;
  NumInitialize_++;
  return 0;
}

int Ifpack_IHSS::SetParameters(Teuchos::ParameterList& parameterlist){
  List_=parameterlist;
  return 0;
}


int Ifpack_IHSS::Compute(){
  if(!IsInitialized_) Initialize();

  int rv;
  Ifpack Factory;
  Epetra_CrsMatrix *Askew=0,*Aherm=0;
  Ifpack_Preconditioner *Pskew=0, *Pherm=0;
  Time_.ResetStartTime();

  // Create Aherm (w/o diagonal)
  rv=EpetraExt::MatrixMatrix::Add(*A_,false,.5,*A_,true,.5,Aherm);
  Aherm->FillComplete();
  if(rv) IFPACK_CHK_ERR(-1);

  // Grab Aherm's diagonal
  Epetra_Vector avec(Aherm->RowMap());
  IFPACK_CHK_ERR(Aherm->ExtractDiagonalCopy(avec));


  // Compute alpha using the Bai, Golub & Ng 2003 formula, not the more multigrid-appropriate Hamilton, Benzi and Haber 2007.
  //  PowerMethod(Aherm, EigMaxIters_,LambdaMax_);
  //  Alpha_=LambdaMax_ / sqrt(EigRatio_);

  // Try something more Hamilton inspired, using the maximum diagonal value of Aherm.
  avec.MaxValue(&Alpha_);

  // Add alpha to the diagonal of Aherm
  for(int i=0;i<Aherm->NumMyRows();i++) avec[i]+=Alpha_;
  IFPACK_CHK_ERR(Aherm->ReplaceDiagonalValues(avec));
  Aherm_=rcp(Aherm);

  // Compute Askew (and add diagonal)
  Askew=new Epetra_CrsMatrix(Copy,A_->RowMap(),0);
  rv=EpetraExt::MatrixMatrix::Add(*A_,false,.5,*A_,true,-.5,Askew);
  if(rv) IFPACK_CHK_ERR(-2);

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Askew->RowMap().GlobalIndicesInt()) {
    for(int i=0;i<Askew->NumMyRows();i++) {
      int gid=Askew->GRID(i);
      Askew->InsertGlobalValues(gid,1,&Alpha_,&gid);
    }
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Askew->RowMap().GlobalIndicesLongLong()) {
    for(int i=0;i<Askew->NumMyRows();i++) {
      long long gid=Askew->GRID64(i);
      Askew->InsertGlobalValues(gid,1,&Alpha_,&gid);
    }
  }
  else
#endif
  throw "Ifpack_IHSS::Compute: Unable to deduce GlobalIndices type";

  Askew->FillComplete();
  Askew_=rcp(Askew);

  // Compute preconditioner for Aherm
  Teuchos::ParameterList PLh=List_.sublist("ihss: hermetian list");
  std::string htype=List_.get("ihss: hermetian type","ILU");
  Pherm= Factory.Create(htype, Aherm);
  Pherm->SetParameters(PLh);
  IFPACK_CHK_ERR(Pherm->Compute());
  Pherm_=rcp(Pherm);

  // Compute preconditoner for Askew
  Teuchos::ParameterList PLs=List_.sublist("ihss: skew hermetian list");
  std::string stype=List_.get("ihss: skew hermetian type","ILU");
  Pskew= Factory.Create(stype, Askew);
  Pskew->SetParameters(PLs);
  IFPACK_CHK_ERR(Pskew->Compute());
  Pskew_=rcp(Pskew);

  // Label
  sprintf(Label_, "IFPACK IHSS (H,S)=(%s/%s)",htype.c_str(),stype.c_str());

  // Counters
  IsComputed_=true;
  NumCompute_++;
  ComputeTime_ += Time_.ElapsedTime();
  return 0;
}


int Ifpack_IHSS::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{
  if(!IsComputed_) return -1;
  Time_.ResetStartTime();
  bool initial_guess_is_zero=false;

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RefCountPtr<const Epetra_MultiVector> Xcopy;
  Epetra_MultiVector Temp(X);
  if (X.Pointers()[0] == Y.Pointers()[0]){
    Xcopy = Teuchos::rcp( new Epetra_MultiVector(X) );
    // Since the user didn't give us anything better, our initial guess is zero.
    Y.Scale(0.0);
    initial_guess_is_zero=true;
  }
  else
    Xcopy = Teuchos::rcp( &X, false );

  Epetra_MultiVector T1(Y),T2(*Xcopy);

  // Note: Since Aherm and Askew are actually (aI+H) and (aI+S) accordingly (to save memory),
  // the application thereof needs to be a little different.
  // The published algorithm is:
  // temp = (aI+H)^{-1} [ (aI-S) y + x ]
  // y = (aI+S)^{-1} [ (aI-H) temp + x ]
  //
  // But we're doing:
  // temp = (aI+H)^{-1} [ 2a y - Shat y + x ]
  // y = (aI+S)^{-1} [ 2 a temp - Hhat temp + x ]

  for(int i=0;i<NumSweeps_;i++){
    // temp = (aI+H)^{-1} [ 2a y - Shat y + x ]
    if(!initial_guess_is_zero || i >0 ){
      Askew_->Apply(Y,T1);
      T2.Update(2*Alpha_,Y,-1,T1,1);
    }
    Pherm_->ApplyInverse(T2,Y);

    // y = (aI+S)^{-1} [ 2 a temp - Hhat temp + x ]
    Aherm_->Apply(Y,T1);
    T2.Scale(1.0,*Xcopy);
    T2.Update(2*Alpha_,Y,-1,T1,1.0);
    Pskew_->ApplyInverse(T2,Y);
  }

  // Counter update
  NumApplyInverse_++;
  ApplyInverseTime_ += Time_.ElapsedTime();
  return 0;
}


std::ostream& Ifpack_IHSS::Print(std::ostream& os) const{
  using std::endl;

  os<<"Ifpack_IHSS"<<endl;
  os<<"-Hermetian preconditioner"<<endl;
  os<<"-Skew Hermetian preconditioner"<<endl;
  os<<endl;
  return os;
}


double Ifpack_IHSS::Condest(const Ifpack_CondestType /* CT */,
                             const int /* MaxIters */,
                             const double /* Tol */,
                             Epetra_RowMatrix* /* Matrix_in */){
  return -1.0;
}


int Ifpack_IHSS::PowerMethod(Epetra_Operator * Op,const int MaximumIterations,  double& lambda_max)
{

  // this is a simple power method
  lambda_max = 0.0;
  double RQ_top, RQ_bottom, norm;
  Epetra_Vector x(Op->OperatorDomainMap());
  Epetra_Vector y(Op->OperatorRangeMap());
  x.Random();
  x.Norm2(&norm);
  if (norm == 0.0) IFPACK_CHK_ERR(-1);

  x.Scale(1.0 / norm);

  for (int iter = 0; iter < MaximumIterations; ++iter)
  {
    Op->Apply(x, y);
    IFPACK_CHK_ERR(y.Dot(x, &RQ_top));
    IFPACK_CHK_ERR(x.Dot(x, &RQ_bottom));
    lambda_max = RQ_top / RQ_bottom;
    IFPACK_CHK_ERR(y.Norm2(&norm));
    if (norm == 0.0) IFPACK_CHK_ERR(-1);
    IFPACK_CHK_ERR(x.Update(1.0 / norm, y, 0.0));
  }

  return(0);
}

#endif
