#ifndef MUELU_SAPFACTORY_DEF_HPP
#define MUELU_SAPFACTORY_DEF_HPP

#include <Xpetra_Operator.hpp>

#include "MueLu_SaPFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Memory.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SaPFactory(RCP<const FactoryBase> InitialPFact, RCP<const FactoryBase> AFact)
    : initialPFact_(InitialPFact), AFact_(AFact),
      dampingFactor_(4./3), diagonalView_("current") {
    if(initialPFact_ == Teuchos::null)
      {
        // use tentative P factory as default
        initialPFact_ = rcp(new TentativePFactory());
      }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~SaPFactory() {}
  
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetDampingFactor(Scalar dampingFactor) {
    dampingFactor_ = dampingFactor;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetDiagonalView(std::string const& diagView) {
    diagonalView_ = diagView;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Scalar SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetDampingFactor() {
    return dampingFactor_;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetDiagonalView() {
    return diagonalView_;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    fineLevel.DeclareInput("A",AFact_.get(),this);
    coarseLevel.DeclareInput("P",initialPFact_.get(),this);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& fineLevel, Level &coarseLevel) const {
    return BuildP(fineLevel,coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildP(Level &fineLevel, Level &coarseLevel) const {
    std::ostringstream buf; buf << coarseLevel.GetLevelID();
    RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("SaPFactory::BuildP_"+buf.str()));
    timer->start(true);

    // Level Get
    RCP<Operator> A     = fineLevel.  Get< RCP<Operator> >("A", AFact_.get());
    RCP<Operator> Ptent = coarseLevel.Get< RCP<Operator> >("P", initialPFact_.get());

    if(restrictionMode_)
      A = Utils2::Transpose(A,true); // build transpose of A explicitely

    //Build final prolongator
    RCP<Operator> finalP; // output

    //FIXME Xpetra::Operator should calculate/stash max eigenvalue
    //FIXME SC lambdaMax = A->GetDinvALambda();

    if (dampingFactor_ != 0) {
      Monitor m(*this, "Prolongator smoothing");

      RCP<Teuchos::Time> sapTimer;
      //sapTimer = rcp(new Teuchos::Time("SaPFactory:I * Ptent"));
      //sapTimer->start(true);
      //Teuchos::ParameterList matrixList;
      //RCP<Operator> I = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map,CrsOperator>("Identity",fineLevel.Get< RCP<Operator> >("A")->getRowMap(),matrixList);
      //RCP<Operator> newPtent = Utils::TwoMatrixMultiply(I,false,Ptent,false);
      //Ptent = newPtent; //I tried a checkout of the original Ptent, and it seems to be gone now (which is good)
      //sapTimer->stop();
      //MemUtils::ReportTimeAndMemory(*sapTimer, *(A->getRowMap()->getComm()));

      sapTimer = rcp(new Teuchos::Time("SaPFactory:APtent_"+buf.str()));
      sapTimer->start(true);

      //JJH -- If I switch doFillComplete to false, the resulting matrix seems weird when printed with describe.
      //JJH -- The final prolongator is wrong, to boot.  So right now, I fillComplete AP, but avoid fillComplete
      //JJH -- in the scaling.  Long story short, we're doing 2 fillCompletes, where ideally we'd do just one.
      bool doFillComplete=true;
      bool optimizeStorage=false;
      RCP<Operator> AP = Utils::TwoMatrixMultiply(A,false,Ptent,false,doFillComplete,optimizeStorage);
      sapTimer->stop();
      MemUtils::ReportTimeAndMemory(*sapTimer, *(A->getRowMap()->getComm()));

      sapTimer = rcp(new Teuchos::Time("SaPFactory:Dinv_APtent_"+buf.str()));
      sapTimer->start(true);
      doFillComplete=false;
      optimizeStorage=false;
      Teuchos::ArrayRCP<SC> diag = Utils::GetMatrixDiagonal(A);
      Utils::MyOldScaleMatrix(AP,diag,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag
      sapTimer->stop();
      MemUtils::ReportTimeAndMemory(*sapTimer, *(A->getRowMap()->getComm()));

      sapTimer = rcp(new Teuchos::Time("SaPFactory:eigen_estimate_"+buf.str()));
      sapTimer->start(true);
      Scalar lambdaMax = Utils::PowerMethod(*A, true, (LO) 10,(Scalar)1e-4);
      //Scalar lambdaMax = Utils::PowerMethod(*A, true, (LO) 50,(Scalar)1e-7, true);
      sapTimer->stop();
      MemUtils::ReportTimeAndMemory(*sapTimer, *(A->getRowMap()->getComm()));
      RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();

      GetOStream(Statistics1, 0) << "Damping factor = " << dampingFactor_/lambdaMax << " (" << dampingFactor_ << " / " << lambdaMax << ")" << std::endl;

      sapTimer = rcp(new Teuchos::Time("SaPFactory:Pt_plus_DinvAPtent_"+buf.str()));
      sapTimer->start(true);

      bool doTranspose=false; 
      if (AP->isFillComplete())
        Utils::TwoMatrixAdd(Ptent,doTranspose,1.0,AP,doTranspose,-dampingFactor_/lambdaMax,finalP);
      else {
        Utils::TwoMatrixAdd(Ptent,doTranspose,1.0,AP,-dampingFactor_/lambdaMax);
        finalP = AP;
      }
      sapTimer->stop();
      MemUtils::ReportTimeAndMemory(*sapTimer, *(A->getRowMap()->getComm()));

      sapTimer = rcp(new Teuchos::Time("SaPFactory:finalP_fillComplete_"+buf.str()));
      sapTimer->start(true);
      finalP->fillComplete( Ptent->getDomainMap(), Ptent->getRangeMap() );
      sapTimer->stop();
      MemUtils::ReportTimeAndMemory(*sapTimer, *(A->getRowMap()->getComm()));
    }
    else {
      finalP = Ptent;
    }

    timer->stop();
    MemUtils::ReportTimeAndMemory(*timer, *(finalP->getRowMap()->getComm()));

    // Level Set
    if(!restrictionMode_)
      {
        // prolongation factory is in prolongation mode
        coarseLevel.Set("P", finalP, this);
      }
    else
      {
        // prolongation factory is in restriction mode
        RCP<Operator> R = Utils2::Transpose(finalP,true); // use Utils2 -> specialization for double
        coarseLevel.Set("R", R, this);
      }

  } //Build()

} //namespace MueLu

#endif // MUELU_SAPFACTORY_DEF_HPP
