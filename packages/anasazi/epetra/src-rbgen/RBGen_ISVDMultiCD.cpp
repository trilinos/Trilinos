// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RBGen_ISVDMultiCD.h"
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_LAPACK.h"
#include "Epetra_BLAS.h"
#include "Epetra_Comm.h"

namespace RBGen {

  ISVDMultiCD::ISVDMultiCD() {}

  void ISVDMultiCD::updateBasis(const Teuchos::RCP< Epetra_MultiVector >& update_ss ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
        "RBGen::ISVDMultiCD::updateBasis(): this routine not supported.");
  }

  void ISVDMultiCD::makePass() {
    Epetra_LAPACK lapack;
    Epetra_BLAS   blas;

    bool firstPass = (curRank_ == 0);
    const int numCols = A_->NumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION( !firstPass && (numProc_ != numCols), std::logic_error,
        "RBGen::ISVDMultiCD::makePass(): after first pass, numProc should be numCols");

    // compute W = I - Z T Z^T from current V_
    Teuchos::RCP<Epetra_MultiVector> lclAZT, lclZ;
    double *Z_A, *AZT_A;
    int Z_LDA, AZT_LDA;
    int oldRank = 0;
    double Rerr = 0.0;
    if (!firstPass) {
      // copy V_ into workZ_
      lclAZT = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*workAZT_,0,curRank_) );
      lclZ   = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*workZ_,0,curRank_) );
      {
        Epetra_MultiVector lclV(Epetra_DataAccess::View,*V_,0,curRank_);
        *lclZ = lclV;
      }
      // compute the Householder QR factorization of the current right basis
      // Vhat = W*R
      int info, lwork = curRank_;
      std::vector<double> tau(curRank_), work(lwork);
      info = lclZ->ExtractView(&Z_A,&Z_LDA);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): error calling ExtractView on Epetra_MultiVector Z.");
      lapack.GEQRF(numCols,curRank_,Z_A,Z_LDA,&tau[0],&work[0],lwork,&info);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): error calling GEQRF on current right basis while constructing next pass coefficients.");
      if (debug_) {
        // we just took the QR factorization of a set of orthonormal vectors
        // they should have an R factor which is diagonal, with unit elements (\pm 1)
        // check it
        Rerr = 0.0;
        for (int j=0; j<curRank_; j++) {
          for (int i=0; i<j; i++) {
            Rerr += abs(Z_A[j*Z_LDA+i]);
          }
          Rerr += abs(abs(Z_A[j*Z_LDA+j]) - 1.0);
        }
      }
      // compute the block representation
      // W = I - Z T Z^T
      lapack.LARFT('F','C',numCols,curRank_,Z_A,Z_LDA,&tau[0],workT_->A(),workT_->LDA());
      // LARFT left upper tri block of Z unchanged
      // note: it should currently contain R factor of V_, which is very close to
      //   diag(\pm 1, ..., \pm 1)
      //
      // we need to set it to:
      //   [1 0 0 ... 0]
      //   [  1 0 ... 0]
      //   [   ....    ]
      //   [          1]
      //
      // see documentation for LARFT
      //
      for (int j=0; j<curRank_; j++) {
        Z_A[j*Z_LDA+j] = 1.0;
        for (int i=0; i<j; i++) {
          Z_A[j*Z_LDA+i] = 0.0;
        }
      }
      // compute part of A W:  A Z T
      // put this in workAZT_
      // first, A Z
      info = lclAZT->Multiply('N','N',1.0,*A_,*lclZ,0.0);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply() for A*Z");
      // second, (A Z) T (in situ, as T is upper triangular)
      info = lclAZT->ExtractView(&AZT_A,&AZT_LDA);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): error calling ExtractView on Epetra_MultiVector AZ.");
      blas.TRMM('R','U','N','N',numCols,curRank_,1.0,workT_->A(),workT_->LDA(),AZT_A,AZT_LDA);
      // save oldRank: it tells us the width of Z
      oldRank  = curRank_;

      curRank_ = 0;
      numProc_ = 0;
    }
    else { // firstPass == true
      curRank_ = 0;
      numProc_ = 0;
    }

    while (numProc_ < numCols) {
      //
      // determine lup
      //
      // want lup >= lmin
      //      lup <= lmax
      // need lup <= numCols - numProc
      //      lup <= maxBasisSize - curRank
      //
      int lup;
      if (curRank_ == 0) {
        // first step uses startRank_
        // this is not affected by lmin,lmax
        lup = startRank_;
      }
      else {
        // this value minimizes overall complexity, assuming fixed rank
        lup = (int)(curRank_ / Teuchos::ScalarTraits<double>::squareroot(2.0));
        // contrain to [lmin,lmax]
        lup = (lup < lmin_ ? lmin_ : lup);
        lup = (lup > lmax_ ? lmax_ : lup);
      }
      //
      // now cap lup via maxBasisSize and the available data
      // these caps apply to all lup, as a result of memory and data constraints
      //
      // available data
      lup = (lup > numCols - numProc_ ? numCols - numProc_ : lup);
      // available memory
      lup = (lup > maxBasisSize_ - curRank_ ? maxBasisSize_ - curRank_ : lup);

      // get view of new vectors
      {
        const Epetra_MultiVector Aplus(Epetra_DataAccess::View,*A_,numProc_,lup);
        Epetra_MultiVector        Unew(Epetra_DataAccess::View,*U_,curRank_,lup);
        // put them in U
        if (firstPass) {
          // new vectors are just Aplus
          Unew = Aplus;
        }
        else {
          // new vectors are Aplus - (A Z T) Z_i^T
          // specifically, Aplus - (A Z T) Z(numProc:numProc+lup-1,1:oldRank)^T
          Epetra_LocalMap lclmap(lup,0,A_->Comm());
          Epetra_MultiVector Zi(Epetra_DataAccess::View,lclmap,&Z_A[numProc_],Z_LDA,oldRank);
          Unew = Aplus;
          int info = Unew.Multiply('N','T',-1.0,*lclAZT,Zi,1.0);
          TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,
              "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply() for A*Wi");
        }
      }

      // perform the incremental step
      incStep(lup);
    }

    // compute W V = V - Z T Z^T V
    // Z^T V is oldRank x curRank
    // T Z^T V is oldRank x curRank
    // we need T Z^T V in a local Epetra_MultiVector
    if (!firstPass) {
      Teuchos::RCP<Epetra_MultiVector> lclV;
      double *TZTV_A;
      int TZTV_LDA;
      int info;
      Epetra_LocalMap lclmap(oldRank,0,A_->Comm());
      // get pointer to current V
      lclV = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*V_,0,curRank_) );
      // create space for T Z^T V
      Epetra_MultiVector TZTV(lclmap,curRank_,false);
      // multiply Z^T V
      info = TZTV.Multiply('T','N',1.0,*lclZ,*lclV,0.0);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply() for Z^T V.");
      // get pointer to data in Z^T V
      info = TZTV.ExtractView(&TZTV_A,&TZTV_LDA);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): error calling ExtractView on Epetra_MultiVector TZTV.");
      // multiply T (Z^T V)
      blas.TRMM('L','U','N','N',oldRank,curRank_,1.0,workT_->A(),workT_->LDA(),TZTV_A,TZTV_LDA);
      // multiply V - Z (T Z^T V)
      info = lclV->Multiply('N','N',-1.0,*lclZ,TZTV,1.0);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply() for W V.");
    }

    //
    // compute the new residuals
    // we know that A V = U S
    // if, in addition, A^T U = V S, then have singular subspaces
    // check residuals A^T U - V S, scaling the i-th column by sigma[i]
    //
    {
      // make these static, because makePass() will be likely be called again
      static Epetra_LocalMap lclmap(A_->NumVectors(),0,A_->Comm());
      static Epetra_MultiVector ATU(lclmap,maxBasisSize_,false);

      // we know that A V = U S
      // if, in addition, A^T U = V S, then have singular subspaces
      // check residuals A^T U - V S, scaling the i-th column by sigma[i]
      Epetra_MultiVector ATUlcl(Epetra_DataAccess::View,ATU,0,curRank_);
      Epetra_MultiVector Ulcl(Epetra_DataAccess::View,*U_,0,curRank_);
      Epetra_MultiVector Vlcl(Epetra_DataAccess::View,*V_,0,curRank_);
      // compute A^T U
      int info = ATUlcl.Multiply('T','N',1.0,*A_,Ulcl,0.0);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply for A^T U.");
      Epetra_LocalMap rankmap(curRank_,0,A_->Comm());
      Epetra_MultiVector S(rankmap,curRank_,true);
      for (int i=0; i<curRank_; i++) {
        S[i][i] = sigma_[i];
      }
      // subtract V S from A^T U
      info = ATUlcl.Multiply('N','N',-1.0,Vlcl,S,1.0);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiCD::computeBasis(): Error calling Epetra_MultiVector::Multiply for V S.");
      resNorms_.resize(curRank_);
      ATUlcl.Norm2(&resNorms_[0]);
      // scale by sigmas
      for (int i=0; i<curRank_; i++) {
        if (sigma_[i] != 0.0) {
          resNorms_[i] /= sigma_[i];
        }
      }
    }

    // debugging checks
    std::vector<double> errnorms(curRank_);
    if (debug_) {
      int info;
      // Check that A V = U Sigma
      // get pointers to current U and V, create workspace for A V - U Sigma
      Epetra_MultiVector work(U_->Map(),curRank_,false), 
                         curU(Epetra_DataAccess::View,*U_,0,curRank_),
                         curV(Epetra_DataAccess::View,*V_,0,curRank_);
      // create local MV for sigmas
      Epetra_LocalMap lclmap(curRank_,0,A_->Comm());
      Epetra_MultiVector curS(lclmap,curRank_,true);
      for (int i=0; i<curRank_; i++) {
        curS[i][i] = sigma_[i];
      }
      info = work.Multiply('N','N',1.0,curU,curS,0.0);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply() for debugging U S.");
      info = work.Multiply('N','N',-1.0,*A_,curV,1.0);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::ISVDMultiCD::makePass(): Error calling Epetra_MultiVector::Multiply() for debugging U S - A V.");
      work.Norm2(&errnorms[0]);
      for (int i=0; i<curRank_; i++) {
        if (sigma_[i] != 0.0) {
          errnorms[i] /= sigma_[i];
        }
      }
    }

    // update pass counter
    curNumPasses_++;

    // print out some info
    const Epetra_Comm *comm = &A_->Comm();
    if (comm->MyPID() == 0 && verbLevel_ >= 1) {
      std::cout 
        << "------------- ISVDMultiCD::makePass() -----------" << std::endl
        << "| Number of passes: " << curNumPasses_ << std::endl
        << "|     Current rank: " << curRank_ << std::endl
        << "|   Current sigmas: " << std::endl;
      for (int i=0; i<curRank_; i++) {
        std::cout << "|             " << sigma_[i] << std::endl;
      }
      if (debug_) {
        std::cout << "|DBG   US-AV norms: " << std::endl;
        for (int i=0; i<curRank_; i++) {
          std::cout << "|DBG          " << errnorms[i] << std::endl;
        }
        if (!firstPass) {
          std::cout << "|DBG      R-I norm: " << Rerr << std::endl;
        }
      }
    }

    return;
  }

  void ISVDMultiCD::Initialize( 
      const Teuchos::RCP< Teuchos::ParameterList >& params,
      const Teuchos::RCP< const Epetra_MultiVector >& ss,
      const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio
      ) 
  {
    workAZT_ = Teuchos::rcp( new Epetra_MultiVector(ss->Map(),maxBasisSize_,false) );

    Epetra_LocalMap lclmap(ss->NumVectors(),0,ss->Comm());
    workZ_ = Teuchos::rcp( new Epetra_MultiVector(lclmap,maxBasisSize_,false) );

    workT_  = Teuchos::rcp( new Epetra_SerialDenseMatrix(maxBasisSize_,maxBasisSize_) );
  }

} // end of RBGen namespace
