// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RBGen_IncSVDPOD.h"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_LAPACK.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Comm.h"

namespace RBGen {

  IncSVDPOD::IncSVDPOD() : 
    isInitialized_(false),
    filter_(Teuchos::null),
    maxBasisSize_(0),
    curRank_(0),
    sigma_(0),
    numProc_(0),
    maxNumPasses_(-1),
    curNumPasses_(0),
    tol_(1e-12),
    lmin_(0),
    lmax_(0),
    startRank_(0),
    timerComp_("Total Elapsed Time"),
    debug_(false),
    verbLevel_(0),
    resNorms_(0),
    Vlocal_(true)
  {}

  Teuchos::RCP<const Epetra_MultiVector> IncSVDPOD::getBasis() const {
    if (curRank_ == 0 || isInitialized_ == false) {
      return Teuchos::null;
    }
    return Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*U_,0,curRank_) );
  }

  Teuchos::RCP<const Epetra_MultiVector> IncSVDPOD::getRightBasis() const {
    if (curRank_ == 0 || isInitialized_ == false) {
      return Teuchos::null;
    }
    return Teuchos::rcp( new Epetra_MultiVector(Epetra_DataAccess::View,*V_,0,curRank_) );
  }

  std::vector<double> IncSVDPOD::getSingularValues() const { 
    std::vector<double> ret(sigma_.begin(),sigma_.begin()+curRank_);
    return ret;
  }

  void IncSVDPOD::Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params,
                              const Teuchos::RCP< const Epetra_MultiVector >& ss,
                              const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio ) {

    using Teuchos::rcp;

    // Get the "Reduced Basis Method" sublist.
    Teuchos::ParameterList rbmethod_params = params->sublist( "Reduced Basis Method" );

    // Get the maximum basis size
    maxBasisSize_ = rbmethod_params.get<int>("Max Basis Size");
    TEUCHOS_TEST_FOR_EXCEPTION(maxBasisSize_ < 2,std::invalid_argument,"\"Max Basis Size\" must be at least 2.");

    // Get a filter
    filter_ = rbmethod_params.get<Teuchos::RCP<Filter<double> > >("Filter",Teuchos::null);
    if (filter_ == Teuchos::null) {
      int k = rbmethod_params.get("Rank",(int)5);
      filter_ = rcp( new RangeFilter<double>(LARGEST,k,k) );
    }

    // Get convergence tolerance
    tol_ = rbmethod_params.get<double>("Convergence Tolerance",tol_);

    // Get debugging flag
    debug_ = rbmethod_params.get<bool>("IncSVD Debug",debug_);

    // Get verbosity level
    verbLevel_ = rbmethod_params.get<int>("IncSVD Verbosity Level",verbLevel_);

    // Get an Anasazi orthomanager
    if (rbmethod_params.isType<
          Teuchos::RCP< Anasazi::OrthoManager<double,Epetra_MultiVector> > 
        >("Ortho Manager")
       ) 
    {
      ortho_ = rbmethod_params.get< 
                Teuchos::RCP<Anasazi::OrthoManager<double,Epetra_MultiVector> >
               >("Ortho Manager");
      TEUCHOS_TEST_FOR_EXCEPTION(ortho_ == Teuchos::null,std::invalid_argument,"User specified null ortho manager.");
    }
    else {
      std::string omstr = rbmethod_params.get("Ortho Manager","DGKS");
      if (omstr == "SVQB") {
        ortho_ = rcp( new Anasazi::SVQBOrthoManager<double,Epetra_MultiVector,Epetra_Operator>() );
      }
      else { // if omstr == "DGKS"
        ortho_ = rcp( new Anasazi::BasicOrthoManager<double,Epetra_MultiVector,Epetra_Operator>() );
      }
    }

    // Lmin,Lmax,Kstart
    lmin_ = rbmethod_params.get("Min Update Size",1);
    TEUCHOS_TEST_FOR_EXCEPTION(lmin_ < 1 || lmin_ >= maxBasisSize_,std::invalid_argument,
                       "Method requires 1 <= min update size < max basis size.");
    lmax_ = rbmethod_params.get("Max Update Size",maxBasisSize_);
    TEUCHOS_TEST_FOR_EXCEPTION(lmin_ > lmax_,std::invalid_argument,"Max update size must be >= min update size.");

    startRank_ = rbmethod_params.get("Start Rank",lmin_);
    TEUCHOS_TEST_FOR_EXCEPTION(startRank_ < 1 || startRank_ > maxBasisSize_,std::invalid_argument,
                       "Starting rank must be in [1,maxBasisSize_)");
    // MaxNumPasses
    maxNumPasses_ = rbmethod_params.get("Maximum Number Passes",maxNumPasses_);
    TEUCHOS_TEST_FOR_EXCEPTION(maxNumPasses_ != -1 && maxNumPasses_ <= 0, std::invalid_argument,
                       "Maximum number passes must be -1 or > 0.");
    // Save the pointer to the snapshot matrix
    TEUCHOS_TEST_FOR_EXCEPTION(ss == Teuchos::null,std::invalid_argument,"Input snapshot matrix cannot be null.");
    A_ = ss;

    // MaxNumCols
    maxNumCols_ = A_->NumVectors();
    maxNumCols_ = rbmethod_params.get("Maximum Number Columns",maxNumCols_);
    TEUCHOS_TEST_FOR_EXCEPTION(maxNumCols_ < A_->NumVectors(), std::invalid_argument,
                       "Maximum number of columns must be at least as many as in the initializing data set.");

    // V locally replicated or globally distributed
    // this must be true for now
    // Vlocal_ = rbmethod_params.get("V Local",Vlocal_);

    // Allocate space for the factorization
    sigma_.reserve( maxBasisSize_ );
    U_ = Teuchos::null;
    V_ = Teuchos::null;
    U_ = rcp( new Epetra_MultiVector(ss->Map(),maxBasisSize_,false) );
    if (Vlocal_) {
      Epetra_LocalMap lclmap(maxNumCols_,0,ss->Comm());
      V_ = rcp( new Epetra_MultiVector(lclmap,maxBasisSize_,false) );
    }
    else {
      Epetra_Map gblmap(maxNumCols_,0,ss->Comm());
      V_ = rcp( new Epetra_MultiVector(gblmap,maxBasisSize_,false) );
    }
    B_ = rcp( new Epetra_SerialDenseMatrix(maxBasisSize_,maxBasisSize_) );
    resNorms_.reserve(maxBasisSize_);

    // clear counters
    numProc_ = 0;
    curNumPasses_ = 0;

    // we are now initialized, albeit with null rank
    isInitialized_ = true;
  }

  void IncSVDPOD::Reset( const Teuchos::RCP<Epetra_MultiVector>& new_ss ) {
    // Reset the pointer for the snapshot matrix
    // Note: We will not assume that it is non-null; user could be resetting our
    // pointer in order to delete the original snapshot set
    A_ = new_ss;
    isInitialized_ = false;
  }

  void IncSVDPOD::computeBasis() {

    //
    // perform enough incremental updates to consume the entire snapshot set
    Teuchos::TimeMonitor lcltimer(timerComp_);

    // check that we have a valid snapshot set: user may have cleared 
    // it using Reset()
    TEUCHOS_TEST_FOR_EXCEPTION(A_ == Teuchos::null,std::logic_error,
                       "computeBasis() requires non-null snapshot set.");

    // check that we are initialized, i.e., data structures match the data set
    TEUCHOS_TEST_FOR_EXCEPTION(isInitialized_==false,std::logic_error,
        "RBGen::IncSVDPOD::computeBasis(): Solver must be initialized.");

    // reset state
    curRank_ = 0;    
    numProc_ = 0;
    curNumPasses_ = 0;

    // print out some info
    const Epetra_Comm *comm = &A_->Comm();
    while (curNumPasses_ < maxNumPasses_ || maxNumPasses_ == -1) {

      // make pass
      makePass();

      // get residuals
      std::vector<double> resnorms = this->getResNorms();

      // check residuals for convergence
      int numConverged = 0;
      for (int i=0; i<curRank_; i++) {
        if (resnorms[i] <= tol_) {
          numConverged++;
        }
      }
      if (comm->MyPID() == 0 && verbLevel_ >= 1) {
        std::cout << "|  Num converged: " << numConverged << std::endl
             << "|    Resid norms: " << std::endl;
        for (int i=0; i<curRank_; i++) {
          std::cout << "|                   " << resnorms[i] << std::endl;
        }
      }
      if (numConverged == curRank_) break;
    }
  }

  
  void IncSVDPOD::incStep(int lup) {

    Epetra_LAPACK lapack;

    // perform gram-schmidt expansion
    // expand() will update bases U_ and V_, factor B_, as well as counters curRank_ and numProc_
    this->expand(lup);

    // compute the SVD of B
    const int lwork = 5*curRank_;
    int info;
    Epetra_SerialDenseMatrix Uhat(Epetra_DataAccess::Copy,B_->A(),B_->LDA(),curRank_,curRank_), Vhat(curRank_,curRank_);
    std::vector<double> Shat(curRank_), work(lwork);

    // Note: this actually stores Vhat^T (we remedy this below)
    lapack.GESVD('O','A',curRank_,curRank_,Uhat.A(),Uhat.LDA(),&Shat[0],Uhat.A(),Uhat.LDA(),Vhat.A(),Vhat.LDA(),&work[0],&lwork,&info);
    TEUCHOS_TEST_FOR_EXCEPTION(info!=0,std::logic_error,"RBGen::IncSVDPOD::incStep(): GESVD return info != 0");

    // use filter to determine new rank
    std::vector<int> keepind = filter_->filter(Shat);
    std::vector<int> truncind(curRank_-keepind.size());
    {
      std::vector<int> allind(curRank_);
      for (int i=0; i<curRank_; i++) {
        allind[i] = i;
      }
      set_difference(allind.begin(),allind.end(),keepind.begin(),keepind.end(),truncind.begin());
      
      // Vhat actually contains Vhat^T; correct this here
      Epetra_SerialDenseMatrix Ucopy(Uhat), Vcopy(curRank_,curRank_); 
      std::vector<double> Scopy(Shat);
      for (int j=0; j<curRank_; j++) {
        for (int i=0; i<curRank_; i++) {
          Vcopy(i,j) = Vhat(j,i);
        }
      }
      // put the desired sigmas at the front of Uhat, Vhat
      for (unsigned int j=0; j<keepind.size(); j++) {
        std::copy(&Ucopy(0,keepind[j]),&Ucopy(curRank_,keepind[j]),&Uhat(0,j));
        std::copy(&Vcopy(0,keepind[j]),&Vcopy(curRank_,keepind[j]),&Vhat(0,j));
        Shat[j] = Scopy[keepind[j]];
      }
      for (unsigned int j=0; j<truncind.size(); j++) {
        std::copy(&Ucopy(0,truncind[j]),&Ucopy(curRank_,truncind[j]),&Uhat(0,keepind.size()+j));
        std::copy(&Vcopy(0,truncind[j]),&Vcopy(curRank_,truncind[j]),&Vhat(0,keepind.size()+j));
        Shat[keepind.size()+j] = Scopy[truncind[j]];
      }
    }

    // shrink back down again
    // shrink() will update bases U_ and V_, as well as singular values sigma_ and curRank_
    this->shrink(truncind.size(),Shat,Uhat,Vhat);

    // print out some info
    const Epetra_Comm *comm = &A_->Comm();
    if (comm->MyPID() == 0 && verbLevel_ >= 2) {
      std::cout 
        << "------------- IncSVDPOD::incStep() --------------" << std::endl
        << "| Cols Processed: " << numProc_ << std::endl
        << "|    Current lup: " << lup << std::endl
        << "|  Current ldown: " << truncind.size() << std::endl
        << "|   Current rank: " << curRank_ << std::endl
        << "| Current sigmas: " << std::endl;
      for (int i=0; i<curRank_; i++) {
        std::cout << "|                   " << sigma_[i] << std::endl;
      }
    }

  }

  const std::vector<double> & IncSVDPOD::getResNorms() {
    return resNorms_;
  }
    
    
} // end of RBGen namespace


