#ifndef MUELU_RAPFACTORY_DEF_HPP
#define MUELU_RAPFACTORY_DEF_HPP

#include <Xpetra_Operator.hpp>
#include <Xpetra_BlockedCrsOperator.hpp>

#include "MueLu_RAPFactory_decl.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RAPFactory(RCP<const FactoryBase> PFact, RCP<const FactoryBase> RFact, RCP<const FactoryBase> AFact)
    : PFact_(PFact), RFact_(RFact), AFact_(AFact), implicitTranspose_(false) {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~RAPFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    fineLevel.DeclareInput("A", AFact_.get(), this);  // AFact per default Teuchos::null -> default factory for this
    coarseLevel.DeclareInput("P", PFact_.get(), this); // transfer operators (from PRFactory, not from PFactory and RFactory!)
    coarseLevel.DeclareInput("R", RFact_.get(), this); //TODO: must be request according to (implicitTranspose flag!!!!!

    // call DeclareInput of all user-given transfer factories
    std::vector<RCP<FactoryBase> >::const_iterator it;
    for(it=TransferFacts_.begin(); it!=TransferFacts_.end(); it++) {
      (*it)->CallDeclareInput(coarseLevel);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {  //FIXME make fineLevel const!!
    typedef Xpetra::BlockedCrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> BlockedCrsOperatorClass; // TODO move me

    {
      int levelID = coarseLevel.GetLevelID();
      FactoryMonitor m(*this, "Computing Ac = RAP", coarseLevel);

      //
      // Inputs: A, P
      //

      RCP<Operator> A = fineLevel.Get< RCP<Operator> >("A", AFact_.get());
      RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", PFact_.get());

      //
      // Build Ac = RAP
      //
      RCP<Operator> Ac;

      const RCP<BlockedCrsOperatorClass> bA = Teuchos::rcp_dynamic_cast<BlockedCrsOperatorClass>(A);
      if( bA != Teuchos::null) {
        RCP<Operator> R = coarseLevel.Get< RCP<Operator> >("R", RFact_.get());
        Ac = BuildRAPBlock(R, A, P, levelID);
      } else {
        if (implicitTranspose_) {
          Ac = BuildRAPImplicit(A, P, levelID);
        } else {
          // explicit version
          RCP<Operator> R = coarseLevel.Get< RCP<Operator> >("R", RFact_.get());
          Ac = BuildRAPExplicit(R, A, P, levelID);
        }
      }

      TEUCHOS_TEST_FOR_EXCEPT(Ac == Teuchos::null);
      
      coarseLevel.Set("A", Ac, this);

    }

    //
    //
    //

    if (TransferFacts_.begin() != TransferFacts_.end()) {
      FactoryMonitor m(*this, "Projections", coarseLevel);

      // call Build of all user-given transfer factories
      std::vector<RCP<FactoryBase> >::const_iterator it;
      for(it = TransferFacts_.begin(); it != TransferFacts_.end(); it++) {
        GetOStream(Runtime0, 0) << "Ac: call transfer factory " << (*it).get() << ": " << (*it)->description() << std::endl;
        (*it)->CallBuild(coarseLevel);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildRAPExplicit(const RCP<Operator>& R, const RCP<Operator>& A, const RCP<Operator>& P, int levelID) const {
    SubFactoryMonitor m(*this, "Build RAP explicitly", levelID);

    RCP<Operator> AP;
    {
      SubFactoryMonitor m2(*this, "MxM: A x P", levelID);
      AP = Utils::TwoMatrixMultiply(A, false, P, false);
      //std::string filename="AP.dat";
      //Utils::Write(filename, AP);
    }

    RCP<Operator> RAP;
    {
      SubFactoryMonitor m2(*this, "MxM: R x (AP)", levelID);
      RAP = Utils::TwoMatrixMultiply(R, false, AP, false);
    }

    GetOStream(Statistics0, 0) << "Ac (explicit): # global rows = " << RAP->getGlobalNumRows() << ", estim. global nnz = " << RAP->getGlobalNumEntries() << std::endl;
    return RAP;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildRAPImplicit(const RCP<Operator>& A, const RCP<Operator>& P, int levelID) const {
    SubFactoryMonitor m(*this, "Build RAP implicitly", levelID);

    GetOStream(Warnings0, 0) << "The implicitTranspose_ flag within RAPFactory for Epetra in parallel produces wrong results" << std::endl;
    RCP<Operator> AP  = Utils::TwoMatrixMultiply(A, false, P, false);
    RCP<Operator> RAP = Utils::TwoMatrixMultiply(P, true, AP, false);

    GetOStream(Statistics0, 0) << "Ac (implicit): # global rows = " << RAP->getGlobalNumRows() << ", estim. global nnz = " << RAP->getGlobalNumEntries() << std::endl;
    return RAP;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildRAPBlock(const RCP<Operator>& R, const RCP<Operator>& A, const RCP<Operator>& P, int levelID) const {
    SubFactoryMonitor m(*this, "Build RAP block", levelID);
    
    typedef Xpetra::BlockedCrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> BlockedCrsOperatorClass;
    const RCP<BlockedCrsOperatorClass> bR = Teuchos::rcp_dynamic_cast<BlockedCrsOperatorClass>(R);
    const RCP<BlockedCrsOperatorClass> bA = Teuchos::rcp_dynamic_cast<BlockedCrsOperatorClass>(A);
    const RCP<BlockedCrsOperatorClass> bP = Teuchos::rcp_dynamic_cast<BlockedCrsOperatorClass>(P);
    TEUCHOS_TEST_FOR_EXCEPTION(bA==Teuchos::null, Exceptions::BadCast, "MueLu::RAPFactory::BuildRAPblock: input matrix A is not of type BlockedCrsOperator! error.");
    TEUCHOS_TEST_FOR_EXCEPTION(bP==Teuchos::null, Exceptions::BadCast, "MueLu::RAPFactory::BuildRAPblock: input matrix P is not of type BlockedCrsOperator! error.");
    TEUCHOS_TEST_FOR_EXCEPTION(bR==Teuchos::null, Exceptions::BadCast, "MueLu::RAPFactory::BuildRAPblock: input matrix R is not of type BlockedCrsOperator! error.");
    if(implicitTranspose_) GetOStream(Warnings0, 0) << "No support for implicitTranspose_ flag within RAPFactory for blocked matrices" << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(bA->Cols()!=bP->Rows(), Exceptions::BadCast, "MueLu::RAPFactory::BuildRAPblock: block matrix dimensions do not match. error.");
    TEUCHOS_TEST_FOR_EXCEPTION(bA->Rows()!=bR->Cols(), Exceptions::BadCast, "MueLu::RAPFactory::BuildRAPblock: block matrix dimensions do not match. error.");

    RCP<BlockedCrsOperatorClass> bAP  = Utils::TwoMatrixMultiplyBlock(bA, false, bP, false, true, true);
    RCP<BlockedCrsOperatorClass> bRAP = Utils::TwoMatrixMultiplyBlock(bR, false, bAP, false, true, true);

    GetOStream(Statistics0, 0) << "Ac (blocked): # global rows = " << bRAP->getGlobalNumRows() << ", estim. global nnz = " << bRAP->getGlobalNumEntries() << std::endl;
    return bRAP;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetImplicitTranspose(bool const &implicit) {
    implicitTranspose_ = implicit;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddTransferFactory(const RCP<FactoryBase>& factory) {
    // check if it's a TwoLevelFactoryBase based transfer factory
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast, "Transfer factory is not derived from TwoLevelFactoryBase. This is very strange. (Note: you can remove this exception if there's a good reason for)");
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
