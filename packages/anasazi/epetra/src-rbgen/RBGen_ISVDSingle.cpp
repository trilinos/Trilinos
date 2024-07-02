// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RBGen_ISVDSingle.h"
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_Comm.h"

namespace RBGen {

  ISVDSingle::ISVDSingle() {}

  void ISVDSingle::updateBasis(const Teuchos::RCP< Epetra_MultiVector >& update_ss ) {
    // free pointer to original data set
    // switch it to new data set, momentarily
    A_ = update_ss;

    // we may be calling update basis from the beginning
    // we hope that V_ (if it exists) is tall enough, we will check this

    // starting with the current factorization, make a single pass
    const int oldNumCols = numProc_;
    const int newNumCols = A_->NumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(oldNumCols+newNumCols > maxNumCols_, std::invalid_argument,
                       "RBGen::ISVSingle::updateBasis(): number of snapshots to process has exceeded specified maximum humber of columns.");
    while (numProc_ < oldNumCols+newNumCols) {
      // determine lup
      int lup;
      if (curRank_ == 0) {
        // first step
        lup = startRank_;
      }
      else {
        // this value minimizes overall complexity for a UDV factorization, assuming fixed rank
        lup = (int)(curRank_ / Teuchos::ScalarTraits<double>::squareroot(2.0));
      }
      // now cap lup via lmin,lmax,maxBasisSize
      // want lup >= lmin
      // need lup <= numCols - numProc
      //      lup <= lmax
      //      lup <= maxBasisSize - curRank
      lup = (lup < lmin_ ? lmin_ : lup);
      lup = (lup > oldNumCols+newNumCols - numProc_ ? oldNumCols+newNumCols - numProc_ : lup);
      lup = (lup > lmax_ ? lmax_ : lup);
      lup = (lup > maxBasisSize_ - curRank_ ? maxBasisSize_ - curRank_ : lup);

      // get view of new vectors
      Teuchos::RCP<const Epetra_MultiVector> Aplus;
      Teuchos::RCP<Epetra_MultiVector> Unew;
      Aplus = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*A_,numProc_-oldNumCols,lup));
      Unew = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*U_,curRank_,lup));
      // put them in U
      *Unew = *Aplus;
      // clear the views
      Unew = Teuchos::null;
      Aplus = Teuchos::null;

      // perform the incremental step
      incStep(lup);
    }

    //
    // we can't compute residuals, because we don't have old and new snapshots
    for (int i=0; i<curRank_; i++) {
      resNorms_[i] = 0;
    }

    // print out some info
    const Epetra_Comm *comm = &A_->Comm();
    if (comm->MyPID() == 0 && verbLevel_ >= 1) {
      std::cout 
        << "------------- ISVDSingle::updateBasis() -----------" << std::endl
        << "|     Current rank: " << curRank_ << std::endl
        << "|   Current sigmas: " << std::endl;
      for (int i=0; i<curRank_; i++) {
        std::cout << "|             " << sigma_[i] << std::endl;
      }
    }

    return;
  }

  void ISVDSingle::makePass() {
    // ISVDSingle only makes a single pass
    TEUCHOS_TEST_FOR_EXCEPTION(maxNumPasses_ != 1,std::logic_error,
        "RBGen::ISVDSingle::makePass(): Max Num Passes should be 1, but is not.");
    // did we already make our one pass? we can't afford to run this again
    if (curNumPasses_ > 0) return;
    const int numCols = A_->NumVectors();
    while (numProc_ < numCols) {
      // determine lup
      int lup;
      if (curRank_ == 0) {
        // first step
        lup = startRank_;
      }
      else {
        // this value minimizes overall complexity, assuming fixed rank
        lup = (int)(curRank_ / Teuchos::ScalarTraits<double>::squareroot(2.0));
      }
      // now cap lup via lmin,lmax,maxBasisSize
      // want lup >= lmin
      // need lup <= numCols - numProc
      //      lup <= lmax
      //      lup <= maxBasisSize - curRank
      lup = (lup < lmin_ ? lmin_ : lup);
      lup = (lup > numCols - numProc_ ? numCols - numProc_ : lup);
      lup = (lup > lmax_ ? lmax_ : lup);
      lup = (lup > maxBasisSize_ - curRank_ ? maxBasisSize_ - curRank_ : lup);

      // get view of new vectors
      Teuchos::RCP<const Epetra_MultiVector> Aplus;
      Teuchos::RCP<Epetra_MultiVector> Unew;
      Aplus = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*A_,numProc_,lup));
      Unew = Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*U_,curRank_,lup));
      // put them in U
      *Unew = *Aplus;
      Unew = Teuchos::null;
      Aplus = Teuchos::null;

      // perform the incremental step
      incStep(lup);
    }

    //
    // compute the new residuals
    // we know that A V = U S
    // if, in addition, A^T U = V S, then have singular subspaces
    // check residuals A^T U - V S, scaling the i-th column by sigma[i]
    //
    {
      Epetra_LocalMap lclmap(A_->NumVectors(),0,A_->Comm());
      Epetra_MultiVector ATU(lclmap,maxBasisSize_,false);

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
        << "------------- ISVDSingle::makePass() -----------" << std::endl
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
      }
    }

    return;
  }
    
  void ISVDSingle::Initialize( 
      const Teuchos::RCP< Teuchos::ParameterList >& params,
      const Teuchos::RCP< const Epetra_MultiVector >& ss,
      const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio
      ) 
  {
    maxNumPasses_ = 1;
  }

} // end of RBGen namespace


