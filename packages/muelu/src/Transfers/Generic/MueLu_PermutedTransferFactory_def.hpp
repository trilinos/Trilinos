#ifndef MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP
#define MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP

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

#ifdef HAVE_MPI
    static double t0=0,t1;

    t0 = MPI_Wtime();
#endif

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
          GetOStream(Runtime0, 0) <<  "Prolongator case" << std::endl;
          RCP<Operator> originalP = coarseLevel.Get< RCP<Operator> >("P",initialTransferFact_.get());
          if (permImporter != Teuchos::null) {
            SubFactoryMonitor m1(*this, "Permuting prolongator", coarseLevel.GetLevelID());
            double tt1 = MPI_Wtime();
            //TODO see note below about domain/range map of R and its implications for P.

            // Now copy P so that we can give it a domain map that matches the range map of R.
            RCP<MultiVector> onesVec  = MultiVectorFactory::Build(originalP->getDomainMap(),1);
            RCP<MultiVector> nnzPerRow  = MultiVectorFactory::Build(originalP->getRangeMap(),1);
            originalP->apply(*onesVec,*nnzPerRow,Teuchos::NO_TRANS,1,0);
            int maxNnzPerRow=0;
            ArrayRCP<const Scalar> nnzCounts = nnzPerRow->getData(0);
            for (int i=0; i<nnzCounts.size(); ++i)
              if (nnzCounts[i] > maxNnzPerRow)
                maxNnzPerRow = (int) nnzCounts[i];
            RCP<Operator> permutedP = OperatorFactory::Build(originalP->getRowMap(), maxNnzPerRow);
            //P needs to be fillCompleted again.  To achieve this, I just copy P using the trivial importer.
            RCP<Import> trivialImporter = ImportFactory::Build( originalP->getRowMap(), originalP->getRowMap());
            RCP<CrsOperator> crsOp = rcp_dynamic_cast<CrsOperator>(permutedP);
            RCP<CrsMatrix> crsMtx = crsOp->getCrsMatrix();
            RCP<CrsOperator> origOp = rcp_dynamic_cast<CrsOperator>(originalP);
            RCP<CrsMatrix> origMtx = origOp->getCrsMatrix();
            GetOStream(Runtime0, 0) <<  "permuting P via importer" << std::endl;
            tt1 = MPI_Wtime();
            crsMtx->doImport(*origMtx,*trivialImporter,Xpetra::INSERT);
            crsMtx = Teuchos::null;
            //new domain (coarse) map, same range (fine) map
            permutedP->fillComplete(permImporter->getTargetMap(), originalP->getRangeMap() );
            tt1 = MPI_Wtime() - tt1;
            if (comm->getRank() == 0)
              std::cout << "Time to permute P with importer = " << tt1 << std::endl;
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
          GetOStream(Runtime0, 0) <<  "Restriction case" << std::endl;
          RCP<Operator> originalR = coarseLevel.Get< RCP<Operator> >("R",initialTransferFact_.get());
          if (permImporter != Teuchos::null) {
            SubFactoryMonitor m1(*this, "Permuting restriction", coarseLevel.GetLevelID());

            // In order to preallocate space for R, we do the following:
            // 1) Allocate ones vector, ovec.
            // 2) Multiply R*ovec, which gives nnz per row.
            // 3) Use permImporter to communicate nnz to new proc owners.
            // 4) Allocate new R, using nnz info. 

            double tt1 = MPI_Wtime();
            RCP<MultiVector> onesVec  = MultiVectorFactory::Build(originalR->getDomainMap(),1);
            RCP<MultiVector> nnzPerRow  = MultiVectorFactory::Build(originalR->getRangeMap(),1);
            originalR->apply(*onesVec,*nnzPerRow,Teuchos::NO_TRANS,1,0);
            int maxNnzPerRow=0;
            ArrayRCP<const Scalar> nnzCounts = nnzPerRow->getData(0);
            for (int i=0; i<nnzCounts.size(); ++i)
              if (nnzCounts[i] > maxNnzPerRow)
                maxNnzPerRow = (int) nnzCounts[i];

            RCP<Operator> permutedR = OperatorFactory::Build(permImporter->getTargetMap(), maxNnzPerRow);
            RCP<CrsOperator> crsOp = rcp_dynamic_cast<CrsOperator>(permutedR);
            RCP<CrsMatrix> crsMtx = crsOp->getCrsMatrix();
            RCP<CrsOperator> origOp = rcp_dynamic_cast<CrsOperator>(originalR);
            RCP<CrsMatrix> origMtx = origOp->getCrsMatrix();
            GetOStream(Runtime0, 0) <<  "permuting R via importer" << std::endl;
            tt1 = MPI_Wtime();
            crsMtx->doImport(*origMtx,*permImporter,Xpetra::INSERT);
            crsMtx = Teuchos::null;
            //TODO is the following range map correct?
            //TODO RangeMap controls where a coarse grid vector's data resides after applying R
            //TODO if the targetMap of the importer is used, then we must either resumeFill or copy/fillComplete P so
            //TODO its domain map matches the range map of R.
            permutedR->fillComplete( originalR->getDomainMap() , permImporter->getTargetMap() );
            tt1 = MPI_Wtime() - tt1;
            if (comm->getRank() == 0)
              std::cout << "Time to permute R with importer = " << tt1 << std::endl;
            //coarseLevel.Set< RCP<Operator> >("newR",permutedR,this);
            coarseLevel.Set< RCP<Operator> >("R",permutedR,this);

            ///////////////////////// EXPERIMENTAL
            // TODO FIXME somehow we have to transfer the striding information of the permuted domain/range maps.
            // That is probably something for an external permutation factory
            //if(originalR->IsView("stridedMaps")) permutedR->CreateView("stridedMaps", originalR);
            ///////////////////////// EXPERIMENTAL
         
            GetOStream(Runtime0, 0) <<  "Stashing 'Importer' into Level " << coarseLevel.GetLevelID() << " w/ generating factory "
                                     << this << std::endl;
            coarseLevel.Set< RCP<const Import> >("Importer",permImporter,this);

            if (coarseLevel.IsAvailable("Coordinates")) //FIXME JJH
            {
              SubFactoryMonitor m2(*this, "Permuting coordinates", coarseLevel.GetLevelID());
              //RCP<MultiVector> coords  = coarseLevel.Get< RCP<MultiVector> >("Coordinates",coordinateFact_.get()); //FIXME JJH
              RCP<MultiVector> coords  = coarseLevel.Get< RCP<MultiVector> >("Coordinates"); //FIXME JJH
              GetOStream(Runtime0, 0) <<  "importing coordinates" << std::endl;
              RCP<MultiVector> permutedCoords  = MultiVectorFactory::Build(permImporter->getTargetMap(),coords->getNumVectors());
              permutedCoords->doImport(*coords,*permImporter,Xpetra::INSERT);
              coarseLevel.Set< RCP<MultiVector> >("Coordinates",permutedCoords); //FIXME JJH no generating factory specified
            }
            if (coarseLevel.IsAvailable("Nullspace",nullspaceFact_.get())) {
              SubFactoryMonitor m2(*this, "Permuting nullspace", coarseLevel.GetLevelID());
              RCP<MultiVector> nullspace  = coarseLevel.Get< RCP<MultiVector> >("Nullspace",nullspaceFact_.get());
              GetOStream(Runtime0, 0) <<  "importing nullspace" << std::endl;
              RCP<MultiVector> permutedNullspace  = MultiVectorFactory::Build(permImporter->getTargetMap(),nullspace->getNumVectors());
              permutedNullspace->doImport(*nullspace,*permImporter,Xpetra::INSERT);
              coarseLevel.Set< RCP<MultiVector> >("Nullspace",permutedNullspace); //FIXME no generating factory specified
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

#ifdef HAVE_MPI
    t1 += MPI_Wtime() - t0;
    if (comm->getRank() == 0)
      std::cout << "cumulative PermutedTransferFactory (excluding get(\"A\")) = " << t1 << std::endl;
#endif

  } //Build

} // namespace MueLu

#define MUELU_PERMUTEDTRANSFER_FACTORY_SHORT
#endif // MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP
