// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RBGen_ISVDUDV.h"

namespace RBGen {

  ISVDUDV::ISVDUDV() {}

  void ISVDUDV::expand(const int lup) 
  {
    //
    // build B and perform gram-schmidt expansion
    //
    //      k l
    // B = [S C] k
    //     [  Z] l
    //
    Teuchos::RCP<const Epetra_MultiVector> U1;
    Teuchos::RCP<Epetra_MultiVector> U2;
    U2 = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*U_,curRank_,lup) );
    int info = (*B_).Scale(0.0);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,
        "RBGen::ISVDUDV::expand(): error calling Epetra_SerialDenseMatrix::scale()");
    for (int i=0; i<curRank_; ++i) {
      (*B_)(i,i) = sigma_[i];
    }
    // get pointer for C,B inside of B, as Teuchos::SerialDenseMatrix objects
    Teuchos::RCP<Epetra_SerialDenseMatrix> C, Z;
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,double> > Cteuchos, Zteuchos;
    if (curRank_ > 0) {
      U1 = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*U_,0,curRank_) );
      C = Teuchos::rcp( new Epetra_SerialDenseMatrix(Epetra_DataAccess::View, &((*B_)(0,curRank_)), B_->LDA(), curRank_, lup) );
      Cteuchos = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,double>(Teuchos::View,C->A(),C->LDA(),curRank_,lup) );
    }
    Z = Teuchos::rcp( new Epetra_SerialDenseMatrix(Epetra_DataAccess::View, &((*B_)(curRank_,curRank_)), B_->LDA(), lup, lup) );
    Zteuchos = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,double>(Teuchos::View, Z->A(), Z->LDA(), lup, lup ) );
    // perform Grams-Schmidt expansion
    int newRank;
    if (curRank_ > 0) {
      newRank = ortho_->projectAndNormalize(*U2,Teuchos::tuple(U1),Teuchos::tuple(Cteuchos),Zteuchos);
    }
    else {
      newRank = ortho_->normalize(*U2,Zteuchos);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(newRank != lup,std::logic_error,
                       "RBGen::ISVDUDV::incStep(): Couldn't recover full rank basis.");
    Cteuchos = Teuchos::null;
    Zteuchos = Teuchos::null;
    C = Teuchos::null;
    Z = Teuchos::null;
    U2 = Teuchos::null;
    U1 = Teuchos::null;
    // augment V with I:
    // Vnew = [V 0]
    //        [0 I]
    //
    // Set V(0:numProc+lup-1,curRank:curRank+lup-1) = 0
    //
    // Vnew = [* 0]
    //        [* 0]
    for (int j=curRank_; j<curRank_+lup; j++) {
      for (int i=0; i<numProc_+lup; i++) {
        (*V_)[j][i] = 0.0;
        // V_->ReplaceGlobalValue(i,j,0.0);
      }
    }
    //
    // Set V(numProc:numProc+lup-1,0:curRank-1) = 0
    // 
    // Vnew = [* *]
    //        [0 *]
    for (int j=0; j<curRank_; j++) {
      for (int i=numProc_; i<numProc_+lup; i++) {
        (*V_)[j][i] = 0.0;
        // V_->ReplaceGlobalValue(i,j,0.0);
      }
    }
    //
    // Set V(numProc:numProc+lup-1,curRank:curRank+lup-1) = I
    //
    // Vnew = [*    *   ]
    //        [* diag(I)]
    for (int j=0; j<lup; j++) {
      (*V_)[curRank_+j][numProc_+j] = 1.0;
      // V_->ReplaceGlobalValue(numProc_+j,curRank_+j,1.0);
    }
    
    numProc_ += lup;
    curRank_ += lup;
  }


  void ISVDUDV::shrink(const int down, std::vector<double> &S, Epetra_SerialDenseMatrix &U, Epetra_SerialDenseMatrix &V) 
  {
    //
    // put RU1 into an Epetra MultiVector
    Epetra_LocalMap LocalMap(curRank_, 0, U_->Map().Comm());
    Epetra_MultiVector Uh1(Epetra_DataAccess::View, LocalMap, U.A(), U.LDA(), curRank_-down);
    Epetra_MultiVector Vh1(Epetra_DataAccess::View, LocalMap, V.A(), V.LDA(), curRank_-down);

    //
    // update bases
    Teuchos::RCP<Epetra_MultiVector> newwU, fullU, newU, newwV, fullV, newV;
    fullU = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*U_,0,curRank_) );
    newwU = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*workU_,0,curRank_-down) );
    // multiply by Uh1
    int info = newwU->Multiply('N','N',1.0,*fullU,Uh1,0.0);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,"ISVDUDV::shrink(): Error calling EMV::Multiply(U).");
    fullU = Teuchos::null;
    newU = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*U_,0,curRank_-down) );
    *newU = *newwU;
    newU = Teuchos::null;
    newwU = Teuchos::null;

    // multiply by Vh1
    // get multivector views of V(1:numProc,1:curRank) and workV(1:numProc,1:curRank-down)
    double *V_A, *workV_A;
    int V_LDA, workV_LDA;
    info = V_->ExtractView(&V_A,&V_LDA);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
        "RBGen::ISVDUDV::shrink(): Error calling Epetra_MultiVector::ExtractView() on V_.");
    info = workV_->ExtractView(&workV_A,&workV_LDA);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
        "RBGen::ISVDUDV::shrink(): Error calling Epetra_MultiVector::ExtractView() on workV_.");
    Epetra_LocalMap lclmap(numProc_,0,A_->Comm());
    fullV = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,lclmap,    V_A,    V_LDA,curRank_     ) );
    newwV = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,lclmap,workV_A,workV_LDA,curRank_-down) );
    // multiply workV = fullV * Vh1
    info = newwV->Multiply('N','N',1.0,*fullV,Vh1,0.0);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,"ISVDUDV::shrink(): Error calling EMV::Multiply(V).");
    fullV = Teuchos::null;
    // now set newV = workV
    newV = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,lclmap, V_A, V_LDA, curRank_-down) );
    *newV = *newwV;
    newV = Teuchos::null;
    newwV = Teuchos::null;

    // save new singular values
    for (int i=0; i<curRank_-down; i++) {
      sigma_[i] = S[i];
    }

    curRank_ = curRank_-down;
  }


  void ISVDUDV::Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params,
                            const Teuchos::RCP< const Epetra_MultiVector >& ss,
                            const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio) 
  {
    workU_ = Teuchos::rcp( new Epetra_MultiVector(ss->Map(),maxBasisSize_,false) );
    Epetra_LocalMap lclmap(ss->NumVectors(),0,ss->Comm());
    workV_ = Teuchos::rcp( new Epetra_MultiVector(lclmap,maxBasisSize_,false) );
  }

} // end of RBGen namespace

