#ifndef MUELU_RAPFACTORY_DEF_HPP
#define MUELU_RAPFACTORY_DEF_HPP

#include "MueLu_RAPFactory_decl.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Memory.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RAPFactory(RCP<FactoryBase> PFact, RCP<FactoryBase> RFact, RCP<FactoryBase> AFact)
    : PFact_(PFact), RFact_(RFact), AFact_(AFact), implicitTranspose_(false) {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~RAPFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    fineLevel.DeclareInput("A", AFact_.get(),this);  // AFact per default Teuchos::null -> default factory for this
    coarseLevel.DeclareInput("P",PFact_.get(),this); // transfer operators (from PRFactory, not from PFactory and RFactory!)
    coarseLevel.DeclareInput("R",RFact_.get(),this); //TODO: must be request according to (implicitTranspose flag!!!!!

    // call DeclareInput of all user-given transfer factories
    std::vector<RCP<FactoryBase> >::const_iterator it;
    for(it=TransferFacts_.begin(); it!=TransferFacts_.end(); it++) {
      (*it)->CallDeclareInput(coarseLevel);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {  //FIXME make fineLevel const!!
    std::ostringstream buf; buf << coarseLevel.GetLevelID();
    RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("RAP::Build_"+buf.str()));
    timer->start(true);

    Teuchos::OSTab tab(this->getOStream());
    RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", PFact_.get());
    RCP<Operator> A = fineLevel.Get< RCP<Operator> >("A",AFact_.get());

    RCP<Teuchos::Time> apTimer = rcp(new Teuchos::Time("RAP::A_times_P_"+buf.str()));
    apTimer->start(true);
    RCP<Operator> AP = Utils::TwoMatrixMultiply(A,false,P,false);
    apTimer->stop();
    MemUtils::ReportTimeAndMemory(*apTimer, *(P->getRowMap()->getComm()));
    //std::string filename="AP.dat";
    //Utils::Write(filename,AP);

    Monitor m(*this, "Computing Ac = RAP");

    RCP<Operator> RAP;
    if (implicitTranspose_) {
      //RCP<Operator> RA = Utils::TwoMatrixMultiply(P,true,A,false);
      //filename = "PtA.dat";
      //Utils::Write(filename,AP);

      // avoid implicitTranspose for Epetra, since EpetraExt matrix-matrix multiplication
      // with implicit transpose flags has bugs. This will hopefully be fixed, soon. (see bug #5363)
      //if(RAP->getRangeMap()->lib() == Xpetra::UseEpetra)
      GetOStream(Warnings0, 0) << "The implicitTranspose_ flag within RAPFactory for Epetra in parallel produces wrong results" << std::endl;

      RAP = Utils::TwoMatrixMultiply(P,true,AP,false);
    } else {
      RCP<Operator> R = coarseLevel.Get< RCP<Operator> >("R", RFact_.get());
      RCP<Teuchos::Time> rapTimer = rcp(new Teuchos::Time("RAP::R_times_AP_"+buf.str()));
      rapTimer->start(true);
      RAP = Utils::TwoMatrixMultiply(R,false,AP,false);
      rapTimer->stop();
      MemUtils::ReportTimeAndMemory(*rapTimer, *(P->getRowMap()->getComm()));

    }
      
    coarseLevel.Set("A", RAP, this);

    timer->stop();
    MemUtils::ReportTimeAndMemory(*timer, *(P->getRowMap()->getComm()));

    GetOStream(Statistics0, 0) << "Ac: # global rows = " << RAP->getGlobalNumRows() << ", estim. global nnz = " << RAP->getGlobalNumEntries() << std::endl;

    // call Build of all user-given transfer factories
    std::vector<RCP<FactoryBase> >::const_iterator it;
    for(it=TransferFacts_.begin(); it!=TransferFacts_.end(); it++) {
      GetOStream(Runtime0, 0) << "Ac: call transfer factory " << (*it).get() << ": " << (*it)->description() << std::endl;
      (*it)->CallBuild(coarseLevel);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetImplicitTranspose(bool const &implicit) {
    implicitTranspose_ = implicit;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddTransferFactory(const RCP<FactoryBase>& factory) {
    // check if it's a TwoLevelFactoryBase based transfer factory
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<TwoLevelFactoryBase>(factory) == Teuchos::null,Exceptions::BadCast, "Transfer factory is not derived from TwoLevelFactoryBase. This is very strange. (Note: you can remove this exception if there's a good reason for)");
    TransferFacts_.push_back(factory);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::NumTransferFactories() const {
    return TransferFacts_.size();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const FactoryBase> RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetPFactory() const {
    return PFact_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const FactoryBase> RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetRFactory() const {
    return RFact_;
  }


} //namespace MueLu

#define MUELU_RAPFACTORY_SHORT
#endif // MUELU_RAPFACTORY_DEF_HPP
