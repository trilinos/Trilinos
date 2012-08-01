#ifndef MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP
#define MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP

#include "Xpetra_Vector.hpp"
#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include <Xpetra_Operator.hpp>
#include <Xpetra_OperatorFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>

#include "MueLu_PermutedTransferFactory_decl.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PermutedTransferFactory(
    RCP<const FactoryBase> repartitionFact,
    RCP<const FactoryBase> initialAFact,
    RCP<const FactoryBase> initialTransferFact,
    TransferType     PorR,
    RCP<const FactoryBase> nullspaceFact,
    RCP<const FactoryBase> coordinateFact )
    : repartitionFact_(repartitionFact),
      initialAFact_(initialAFact),
      initialTransferFact_(initialTransferFact),
      PorR_(PorR),
      nullspaceFact_(nullspaceFact),
      coordinateFact_(coordinateFact)
  { }

  //-----------------------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~PermutedTransferFactory() {}

  //-----------------------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    coarseLevel.DeclareInput("A", initialAFact_.get(),this);
    if (PorR_ == MueLu::INTERPOLATION) {
      coarseLevel.DeclareInput("P",initialTransferFact_.get(),this);
    } else {
      coarseLevel.DeclareInput("R",initialTransferFact_.get(),this);
      coarseLevel.DeclareInput("Nullspace",nullspaceFact_.get(),this);
    }

    coarseLevel.DeclareInput("Importer",repartitionFact_.get(),this);
  }

  //-----------------------------------------------------------------------------------------------------

  // TODO 1) Need to pass in to ctor generating factories for nullspace (TentativePFactory) and coordinates
  // TODO (MultiVectorTransferFactory).  DeclareInput should also call the DeclareInputs of these two guys.
  // TODO 2) Also add interface (ala RAPFactory) to register additional data/generating factories that must
  // TODO be permuted.  Call those DeclareInputs, as well?

  // partially done

  //-----------------------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    static RCP<const Teuchos::Comm<int> > comm;

    if (PorR_ == MueLu::INTERPOLATION) {
      GetOStream(Warnings0, 0) <<  "Jamming A into Level " << coarseLevel.GetLevelID() << " w/ generating factory "
                               << this << std::endl;
      RCP<Operator> A = coarseLevel.Get< RCP<Operator> >("A",initialAFact_.get());
      coarseLevel.Set< RCP<Operator> >("A",A,this);
      comm = A->getRowMap()->getComm();
    }

    RCP<const Import> permImporter;
    try {
      permImporter = coarseLevel.Get< RCP<const Import> >("Importer",repartitionFact_.get());
    }
    catch(MueLu::Exceptions::HaltRepartitioning e) {
      std::string gridTransferType;
      GetOStream(Warnings0, 0) <<  "Skipping permuting of "
                               << ((PorR_ == MueLu::INTERPOLATION) ? "prolongator" : "restriction")
                               << ".  No permutation is available for the following reason:"
      << std::endl << e.what() << std::endl;
    }

    switch (PorR_) {

      case MueLu::INTERPOLATION:
        { //case scoping
          RCP<Operator> originalP = coarseLevel.Get< RCP<Operator> >("P",initialTransferFact_.get());
          if (permImporter != Teuchos::null) {
            SubFactoryMonitor m1(*this, "Rebalancing prolongator", coarseLevel.GetLevelID());
            //TODO see note below about domain/range map of R and its implications for P.

            // Now copy P so that we can give it a domain map that matches the range map of R.
            //TODO is there a better preallocation strategy?
            ArrayRCP<LocalOrdinal> nnzPerRow(originalP->getNodeNumRows(),0);
            for (size_t i=0; i<originalP->getNodeNumRows(); ++i)
              nnzPerRow[i] = originalP->getNumEntriesInLocalRow(i);
            RCP<Operator> permutedP = OperatorFactory::Build(originalP->getRowMap(), nnzPerRow, Xpetra::StaticProfile);
            //P needs to be fillCompleted again.  To achieve this, I just copy P using the trivial importer.
            RCP<Import> trivialImporter = ImportFactory::Build( originalP->getRowMap(), originalP->getRowMap());
            RCP<CrsOperator> crsOp = rcp_dynamic_cast<CrsOperator>(permutedP);
            RCP<CrsMatrix> crsMtx = crsOp->getCrsMatrix();
            RCP<CrsOperator> origOp = rcp_dynamic_cast<CrsOperator>(originalP);
            RCP<CrsMatrix> origMtx = origOp->getCrsMatrix();
            RCP<SubFactoryMonitor> m2 = rcp( new SubFactoryMonitor(*this, "Rebalancing prolongator -- import only",
                                                                   coarseLevel.GetLevelID()) );
            crsMtx->doImport(*origMtx,*trivialImporter,Xpetra::INSERT);
            crsMtx = Teuchos::null;
            //new domain (coarse) map, same range (fine) map
            m2 = rcp( new SubFactoryMonitor(*this, "Rebalancing prolongator -- fillComplete", coarseLevel.GetLevelID()) );
            permutedP->fillComplete(permImporter->getTargetMap(), originalP->getRangeMap() );
            m2 = Teuchos::null;
            originalP = Teuchos::null;
            coarseLevel.Set< RCP<Operator> >("P",permutedP,this);

            ///////////////////////// EXPERIMENTAL
            // TODO FIXME somehow we have to transfer the striding information of the permuted domain/range maps.
            // That is probably something for an external permutation factory
            //if(originalP->IsView("stridedMaps")) permutedP->CreateView("stridedMaps", originalP);
            ///////////////////////// EXPERIMENTAL

          } else {
            coarseLevel.Set< RCP<Operator> >("P",originalP,this);
          } //if (permImporter != Teuchos::null) {...} else {...}
        } //case scoping
        break;

      case MueLu::RESTRICTION:
        { //case scoping
          //TODO how do we handle implicitly transposed restriction operators?
          RCP<Operator> originalR = coarseLevel.Get< RCP<Operator> >("R",initialTransferFact_.get());
          if (permImporter != Teuchos::null) {
            SubFactoryMonitor m2(*this, "Rebalancing restriction", coarseLevel.GetLevelID());
            RCP<SubFactoryMonitor> m1 = rcp( new SubFactoryMonitor(*this, "Rebalancing restriction -- allocate new R", coarseLevel.GetLevelID()) );

            //TODO is there a better preallocation strategy?
            /*
              Answer:  YES! :
               Count nnz per row of R
               Put this in a MV
               Use the permImporter to communicate this
               Allocate permutedR according to the nnz info
               Proceed as usual
            */
            //Note, this is to avoid that Epetra only supports vectors of type int or double, not size_t, which is unsigned.
            RCP<Xpetra::Vector<double,LO,GO> > originalNnzPerRowVec = Xpetra::VectorFactory<double,LO,GO>::Build(permImporter->getSourceMap());
            ArrayRCP<double> nnzPerRow = originalNnzPerRowVec->getDataNonConst(0);
            for (size_t i=0; i<originalR->getNodeNumRows(); ++i)
              nnzPerRow[i] = originalR->getNumEntriesInLocalRow(i);
            nnzPerRow = Teuchos::null;
            RCP<Xpetra::Vector<double,LO,GO> > permutedNnzPerRowVec = Xpetra::VectorFactory<double,LO,GO>::Build(permImporter->getTargetMap());
            permutedNnzPerRowVec->doImport(*originalNnzPerRowVec,*permImporter,Xpetra::INSERT);
            ArrayRCP<const double> tmpData = permutedNnzPerRowVec->getData(0);
            ArrayRCP<LocalOrdinal> permutedNnzPerRow(permutedNnzPerRowVec->getLocalLength());
            for (size_t i=0; i<permutedNnzPerRowVec->getLocalLength(); ++i) 
              permutedNnzPerRow[i] = Teuchos::as<LocalOrdinal,double>(tmpData[i]);
            RCP<Operator> permutedR = OperatorFactory::Build(permImporter->getTargetMap(), permutedNnzPerRow, Xpetra::StaticProfile);
            permutedNnzPerRow = Teuchos::null;
            RCP<CrsOperator> crsOp = rcp_dynamic_cast<CrsOperator>(permutedR);
            RCP<CrsMatrix> crsMtx = crsOp->getCrsMatrix();
            RCP<CrsOperator> origOp = rcp_dynamic_cast<CrsOperator>(originalR);
            RCP<CrsMatrix> origMtx = origOp->getCrsMatrix();
            m1 = Teuchos::null;
            m1 = rcp( new SubFactoryMonitor(*this, "Rebalancing restriction -- import", coarseLevel.GetLevelID()) );
            crsMtx->doImport(*origMtx,*permImporter,Xpetra::INSERT);
            crsMtx = Teuchos::null;
            //TODO is the following range map correct?
            //TODO RangeMap controls where a coarse grid vector's data resides after applying R
            //TODO if the targetMap of the importer is used, then we must either resumeFill or copy/fillComplete P so
            //TODO its domain map matches the range map of R.
            m1 = Teuchos::null;
            m1 = rcp( new SubFactoryMonitor(*this, "Rebalancing restriction -- fillComplete", coarseLevel.GetLevelID()) );
            permutedR->fillComplete( originalR->getDomainMap() , permImporter->getTargetMap() );
            //coarseLevel.Set< RCP<Operator> >("newR",permutedR,this);
            m1 = Teuchos::null;
            coarseLevel.Set< RCP<Operator> >("R",permutedR,this);

            ///////////////////////// EXPERIMENTAL
            // TODO FIXME somehow we have to transfer the striding information of the permuted domain/range maps.
            // That is probably something for an external permutation factory
            //if(originalR->IsView("stridedMaps")) permutedR->CreateView("stridedMaps", originalR);
            ///////////////////////// EXPERIMENTAL
         
            coarseLevel.Set< RCP<const Import> >("Importer",permImporter,this);

            if (coarseLevel.IsAvailable("Coordinates")) //FIXME JJH
            {
              m1 = rcp( new SubFactoryMonitor(*this, "Rebalancing coordinates", coarseLevel.GetLevelID()) );
              //RCP<MultiVector> coords  = coarseLevel.Get< RCP<MultiVector> >("Coordinates",coordinateFact_.get()); //FIXME JJH
              RCP<MultiVector> coords  = coarseLevel.Get< RCP<MultiVector> >("Coordinates"); //FIXME JJH
              RCP<MultiVector> permutedCoords  = MultiVectorFactory::Build(permImporter->getTargetMap(),coords->getNumVectors());
              permutedCoords->doImport(*coords,*permImporter,Xpetra::INSERT);
              coarseLevel.Set< RCP<MultiVector> >("Coordinates",permutedCoords); //FIXME JJH no generating factory specified
              m1 = Teuchos::null;
            }
            if (coarseLevel.IsAvailable("Nullspace",nullspaceFact_.get())) {
              m1 = rcp( new SubFactoryMonitor(*this, "Rebalancing nullspace", coarseLevel.GetLevelID()) );
              RCP<MultiVector> nullspace  = coarseLevel.Get< RCP<MultiVector> >("Nullspace",nullspaceFact_.get());
              RCP<MultiVector> permutedNullspace  = MultiVectorFactory::Build(permImporter->getTargetMap(),nullspace->getNumVectors());
              permutedNullspace->doImport(*nullspace,*permImporter,Xpetra::INSERT);
              coarseLevel.Set< RCP<MultiVector> >("Nullspace",permutedNullspace); //FIXME no generating factory specified
              m1 = Teuchos::null;
            }
          } else {
            coarseLevel.Set< RCP<Operator> >("R",originalR,this);
            if (coarseLevel.IsAvailable("Nullspace",nullspaceFact_.get())) {
              RCP<MultiVector> nullspace  = coarseLevel.Get< RCP<MultiVector> >("Nullspace",nullspaceFact_.get());
              coarseLevel.Set< RCP<MultiVector> >("Nullspace",nullspace, this);
            }
          } //if (permImporter != Teuchos::null) {...} else {...}
        } //case scoping
        break;
    } //switch

  } //Build

} // namespace MueLu

#define MUELU_PERMUTEDTRANSFER_FACTORY_SHORT
#endif // MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP
