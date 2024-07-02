// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RBGen_ISVDMultiSDB.h"
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_LAPACK.h"
#include "Epetra_BLAS.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"

namespace RBGen {

  ISVDMultiSDB::ISVDMultiSDB() {}

  void ISVDMultiSDB::updateBasis(const Teuchos::RCP< Epetra_MultiVector >& update_ss ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
        "RBGen::ISVDMultiSDB::updateBasis(): this routine not supported.");
  }
  void ISVDMultiSDB::makePass() {
    Epetra_LAPACK lapack;
    Epetra_BLAS   blas;

    bool firstPass = (curRank_ == 0);
    const int numCols = A_->NumVectors();

    // compute W = I - Z T Z^T from current V_
    int oldRank = 0;
    double Rerr;
    if (!firstPass) {
      // we want W = [W1 W2] = [V G]
      // second block == G is already achieved from previous pass
      // we need to put V in the first block
      {
        Epetra_MultiVector W1(Epetra_DataAccess::View,*workW_,0,curRank_);
        Epetra_MultiVector lclV(Epetra_DataAccess::View,*V_,0,curRank_);
        W1 = lclV;
      }
      // get a pointer to [W1 W2]
      Epetra_MultiVector W12(Epetra_DataAccess::View,*workW_,0,2*curRank_);

      //
      // compute the Householder QR factorization of the current right basis
      // [V G] = W*[R,0]'
      //
      int info, lwork = 2*curRank_;
      std::vector<double> tau(2*curRank_), work(lwork);
      lapack.GEQRF(numCols,2*curRank_,W12[0],numCols,&tau[0],&work[0],lwork,&info);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiSDB::makePass(): error calling GEQRF on current right basis while constructing next pass coefficients.");
      if (debug_) {
        // we just took the QR factorization of [V G]
        // V is has orthonormal columns, so that the leading part of R should 
        // be diagonal, with unit elements (\pm 1)
        // check it
        Rerr = 0;
        for (int j=0; j<curRank_; j++) {
          for (int i=0; i<j; i++) {
            Rerr += abs(W12[j][i]);
          }
          Rerr += abs(abs(W12[j][j]) - 1.0);
        }
      }
      // compute the orthonormal basis ([W_1,W_2] = [\pm V,G]) 
      // will overwrite R
      lapack.ORGQR(numCols,2*curRank_,2*curRank_,W12[0],numCols,&tau[0],&work[0],lwork,&info);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::ISVDMultiSDB::makePass(): error calling ORGQR to construct next pass coefficients.");
      // compute A [W_1 W_2]
      {
        Epetra_MultiVector lclAW(Epetra_DataAccess::View,*workAW_,0,2*curRank_);
        info = lclAW.Multiply('N','N',1.0,*A_,W12,0.0);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,
            "RBGen::ISVDMultiSDB::makePass(): Error calling Epetra_MultiVector::Multiply() for A*W");
      }
      // save oldRank, which indicates the width of V
      oldRank  = curRank_;
      curRank_ = 0;
    }

    numProc_ = 0;
    while ( ( firstPass && numProc_ <   numCols) || 
            (!firstPass && numProc_ < 2*oldRank) ) {
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
      // available data: 
      //   if firstPass, we are consuming A (numCols vectors)
      //   if !firstPass, we are consuming A*W (2*oldRank vectors)
      if (firstPass) {
        lup = (lup > numCols - numProc_ ? numCols - numProc_ : lup);
      }
      else {
        lup = (lup > 2*oldRank - numProc_ ? 2*oldRank - numProc_ : lup);
      }
      // available memory
      lup = (lup > maxBasisSize_ - curRank_ ? maxBasisSize_ - curRank_ : lup);

      // put new vectors in U
      {
        Epetra_MultiVector Unew(Epetra_DataAccess::View,*U_,curRank_,lup);
        if (firstPass) {
          // new vectors are just Aplus
          const Epetra_MultiVector Aplus(Epetra_DataAccess::View,*A_,numProc_,lup);
          Unew = Aplus;
        }
        else {
          // new vectors are AW
          const Epetra_MultiVector AWplus(Epetra_DataAccess::View,*workAW_,numProc_,lup);
          Unew = AWplus;
        }
      }

      // perform the incremental step
      incStep(lup);
    }

    // compute W V = V - Z T Z^T V
    // Z^T V is 2*oldRank x curRank
    // T Z^T V is 2*oldRank x curRank
    // we need T Z^T V in a local Epetra_MultiVector
    if (!firstPass) {
      Epetra_MultiVector lclWV(Epetra_DataAccess::View,*workWV_,0,curRank_);
      {
        // create (local) map for V(1:numProc,1:curRank)
        Epetra_LocalMap lclmap(numProc_,0,A_->Comm());
        const Epetra_MultiVector lclV(Epetra_DataAccess::View,lclmap,(*V_)[0],numCols,curRank_);
        const Epetra_MultiVector W12(Epetra_DataAccess::View,*workW_,0,2*oldRank);
        int info = lclWV.Multiply('N','N',1.0,W12,lclV,0.0);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
            "RBGen::ISVDMultiSDB::makePass(): error Multiply()ing Epetra_MultiVector for W*V.");
      }
      Epetra_MultiVector lclV(Epetra_DataAccess::View,*V_,0,curRank_);
      lclV = lclWV;
    }

    //
    // compute the new residuals
    // we know that A V = U S
    // if, in addition, A^T U = V S, then have singular subspaces
    // check residuals A^T U - V S, scaling the i-th column by sigma[i]
    //
    // store A^T U - V S into the second block of workW_
    // we will need it on the next pass 
    //
    {
      Epetra_MultiVector W2(Epetra_DataAccess::View,*workW_,curRank_,curRank_);
      Epetra_MultiVector Ulcl(Epetra_DataAccess::View,*U_,0,curRank_);
      Epetra_MultiVector Vlcl(Epetra_DataAccess::View,*V_,0,curRank_);

      //_DataAccess compute A^T U
      int info = W2.Multiply('T','N',1.0,*A_,Ulcl,0.0);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::IncSVD::computeBasis(): Error calling Epetra_MultiVector::Multiply for A^T U.");
      Epetra_LocalMap Smap(curRank_,0,A_->Comm());
      Epetra_MultiVector S(Smap,curRank_,true); // "true" inits to zero
      for (int i=0; i<curRank_; i++) {
        S[i][i] = sigma_[i];
      }
      // subtract V S from A^T U
      info = W2.Multiply('N','N',-1.0,Vlcl,S,1.0);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
          "RBGen::IncSVD::computeBasis(): Error calling Epetra_MultiVector::Multiply for V S.");

      //
      // compute residual norms
      resNorms_.resize(curRank_);
      W2.Norm2(&resNorms_[0]);
      // scale by sigmas
      for (int i=0; i<curRank_; i++) {
        if (sigma_[i] != 0.0) {
          resNorms_[i] /= sigma_[i];
        }
      }
    }

    //
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
          "RBGen::ISVDMultiSDB::makePass(): Error calling Epetra_MultiVector::Multiply() for debugging U S.");
      info = work.Multiply('N','N',-1.0,*A_,curV,1.0);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,
          "RBGen::ISVDMultiSDB::makePass(): Error calling Epetra_MultiVector::Multiply() for debugging U S - A V.");
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
        << "------------- ISVDMultiSDB::makePass() -----------" << std::endl
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

  void ISVDMultiSDB::Initialize( 
      const Teuchos::RCP< Teuchos::ParameterList >& params,
      const Teuchos::RCP< const Epetra_MultiVector >& ss,
      const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio
      ) 
  {
    // workAW has room for A * [V G], where each has maxBasisSize vectors
    // ergo, workAW has 2*maxBasisSize vectors
    workAW_ = Teuchos::rcp( new Epetra_MultiVector(ss->Map(),2*maxBasisSize_,false) );

    // workZ = [V G],
    // where V is the current right basis and G is the gradient A^T U - V S
    // ergo, workZ has 2*maxBasisSize_ vectors
    Epetra_LocalMap lclmap(ss->NumVectors(),0,ss->Comm());
    workW_ = Teuchos::rcp( new Epetra_MultiVector(lclmap,2*maxBasisSize_,false) );
    workWV_ = Teuchos::rcp( new Epetra_MultiVector(lclmap,maxBasisSize_,false) );
  }

} // end of RBGen namespace
