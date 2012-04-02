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
    coarseLevel.DeclareInput("A", PFact_.get(), this); //FIXME hack
    coarseLevel.DeclareInput("Permutation", RFact_.get(), this); //FIXME hack, could result in redundant work...
    coarseLevel.DeclareInput("PermutationTrans", PFact_.get(), this); //FIXME hack, could result in redundant work...

    // call DeclareInput of all user-given transfer factories
    for(std::vector<RCP<const FactoryBase> >::const_iterator it = TransferFacts_.begin(); it!=TransferFacts_.end(); ++it) {
      (*it)->CallDeclareInput(coarseLevel);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {  //FIXME make fineLevel const!!
    typedef Xpetra::BlockedCrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> BlockedCrsOperatorClass; // TODO move me

    {
      FactoryMonitor m(*this, "Computing Ac = RAP", coarseLevel);

      int levelID = coarseLevel.GetLevelID();

      GetOStream(Runtime0, 0) << "Ac: building coarse operator" << std::endl;

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef PRINT_THE_LEVELS_FOR_DEBUGGING
          //JJH FIXME DEBUGGING PRINTS
    RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    fos->setOutputToRootOnly(0);
    *fos << "================== fine level in RAPFactory =================" << std::endl;
    for (int i=0; i<comm->getSize(); ++i) {
      if (comm->getRank() == i) {
        fos->setOutputToRootOnly(comm->getRank());
        fineLevel.print(*fos,Teuchos::VERB_EXTREME);
        fos->setShowProcRank(false);
        if (i<comm->getSize()-1)
          *fos << "        ------------------------------------------------------------------------------------------------------------------------"
               << std::endl;
      }
      comm->barrier();
    }
    fos->setOutputToRootOnly(0);
    *fos << "================= end of fine level description =========================" << std::endl;
    *fos << "================== coarse level in RAPFactory =================" << std::endl;
    for (int i=0; i<comm->getSize(); ++i) {
      if (comm->getRank() == i) {
        fos->setOutputToRootOnly(comm->getRank());
        fineLevel.print(*fos,Teuchos::VERB_EXTREME);
        fos->setShowProcRank(false);
        if (i<comm->getSize()-1)
          *fos << "        ------------------------------------------------------------------------------------------------------------------------"
               << std::endl;
      }
      comm->barrier();
    }
    fos->setOutputToRootOnly(0);
    *fos << "================= end of coarse level description =========================" << std::endl;
#endif
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

          // explicit version
          if (coarseLevel.IsAvailable("Permutation",RFact_.get())) GetOStream(Runtime0, 0) << "Permutation generated by " << RFact_.get() << " is available!!!!" << std::endl;
          if (coarseLevel.IsAvailable("PermutationTrans",PFact_.get())) GetOStream(Runtime0, 0) << "trans(Permutation) generated by " << PFact_.get() << " is available!!!!" << std::endl;
          if (coarseLevel.IsAvailable("A",PFact_.get())) GetOStream(Runtime0, 0) << "A generated by " << PFact_.get() << " is available!!!!" << std::endl;
          if ( coarseLevel.IsAvailable("A",PFact_.get()) ) {
            GetOStream(Runtime0, 0) << "retrieving Ac generated by " << PFact_.get() << std::endl;
            Ac = coarseLevel.Get< RCP<Operator> >("A", PFact_.get());
            if ( coarseLevel.IsAvailable("Permutation",RFact_.get()) &&
                 coarseLevel.IsAvailable("PermutationTrans",PFact_.get()) ) {
              GetOStream(Runtime0, 0) << "Found both Permutation and PermutationTrans!!!!" << std::endl;
              RCP<Operator> Perm = coarseLevel.Get< RCP<Operator> >("Permutation", RFact_.get());
              RCP<Operator> transPerm = coarseLevel.Get< RCP<Operator> >("PermutationTrans", PFact_.get());
              Ac = BuildRAPExplicit(Perm, Ac, transPerm, levelID);
            }
          } else {
            RCP<Operator> R = coarseLevel.Get< RCP<Operator> >("R", RFact_.get());
            Ac = BuildRAPExplicit(R, A, P, levelID);
          }

/*
          //JJH FIXME DEBUGGING PRINTS
          RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
          fos->setOutputToRootOnly(-1);
          RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
          if (comm->getRank() == 0)
            *fos << "=======================\n In RAP:  printing Acoarse\n===========================" << std::endl;
          comm->barrier();
          Ac->describe(*fos,Teuchos::VERB_EXTREME);
          comm->barrier();
          //JJH FIXME END OF DEBUGGING PRINTS
*/

   /*
   static int counter=0; counter++;

      std::string filename="RAP-Level-" + Teuchos::toString(coarseLevel.GetLevelID()) + "-Afine-id" + Teuchos::toString(counter) +
".dat";
      Utils::Write(filename, *A);

      std::string filename2="RAP-Level-" + Teuchos::toString(coarseLevel.GetLevelID()) + "-P-id" + Teuchos::toString(counter) +
".dat";
      Utils::Write(filename2, *P);

      std::string filename3="RAP-Level-" + Teuchos::toString(coarseLevel.GetLevelID()) + "-Acoarse-id" +
Teuchos::toString(counter) + ".dat";
      Utils::Write(filename3, *Ac); 
      */
        }
      }

      TEUCHOS_TEST_FOR_EXCEPT(Ac == Teuchos::null);
      
      coarseLevel.Set("A", Ac, this);

    }

    //
    //
    //

    if (TransferFacts_.begin() != TransferFacts_.end()) {
      SubFactoryMonitor m(*this, "Projections", coarseLevel);

      // call Build of all user-given transfer factories
      for(std::vector<RCP<const FactoryBase> >::const_iterator it = TransferFacts_.begin(); it != TransferFacts_.end(); ++it) {
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

    //nonzero imbalance
    size_t numMyNnz  = RAP->getNodeNumEntries();
    GO maxNnz, minNnz;
    RCP<const Teuchos::Comm<int> > comm = RAP->getRowMap()->getComm();
    maxAll(comm,(GO)numMyNnz,maxNnz);
    //min nnz over all proc (disallow any processors with 0 nnz)
    minAll(comm, (GO)((numMyNnz > 0) ? numMyNnz : maxNnz), minNnz);
    Scalar imbalance = ((SC) maxNnz) / minNnz;

    //min and max # rows per proc
    size_t numMyRows = RAP->getNodeNumRows();
    GO minNumRows, maxNumRows;
    maxAll(comm, (GO)numMyRows, maxNumRows);
    minAll(comm, (GO)((numMyRows > 0) ? numMyRows : maxNumRows), minNumRows);

    //Check whether RAP is spread over more than one process.
    GO numActiveProcesses=0;
    sumAll(comm, (GO)((RAP->getNodeNumRows() > 0) ? 1 : 0), numActiveProcesses);

    GetOStream(Statistics1,0) << "Ac (explicit) # processes with rows = " << numActiveProcesses << std::endl;
    GetOStream(Statistics1,0) << "Ac (explicit) min # rows per proc = " << minNumRows << ", max # rows per proc = " << maxNumRows << std::endl;
    GetOStream(Statistics1,0) << "Ac (explicit) nonzero imbalance = " << imbalance << std::endl;

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
  void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddTransferFactory(const RCP<const FactoryBase>& factory) {
    // check if it's a TwoLevelFactoryBase based transfer factory
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast, "Transfer factory is not derived from TwoLevelFactoryBase. This is very strange. (Note: you can remove this exception if there's a good reason for)");
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
