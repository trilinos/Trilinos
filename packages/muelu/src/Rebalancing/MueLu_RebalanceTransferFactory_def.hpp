// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_REBALANCETRANSFERFACTORY_DEF_HPP
#define MUELU_REBALANCETRANSFERFACTORY_DEF_HPP

#include "Xpetra_Vector.hpp"
#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>

#include "MueLu_RebalanceTransferFactory_decl.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(coarseLevel, "A");
    if (PorR_ == MueLu::INTERPOLATION) {
      Input(coarseLevel, "P");
    } else {
      Input(coarseLevel, "R");
      Input(coarseLevel, "Nullspace");
    }

    Input(coarseLevel, "Importer");
  }

  //-----------------------------------------------------------------------------------------------------

  // TODO 1) Need to pass in to ctor generating factories for nullspace (TentativePFactory) and coordinates
  // TODO (MultiVectorTransferFactory).  DeclareInput should also call the DeclareInputs of these two guys.
  // TODO 2) Also add interface (ala RAPFactory) to register additional data/generating factories that must
  // TODO be permuted.  Call those DeclareInputs, as well?

  // partially done

  //-----------------------------------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RebalanceTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    RCP<const Import> rebalanceImporter;
    try {
      rebalanceImporter = Get< RCP<const Import> >(coarseLevel, "Importer");
    }
    catch(MueLu::Exceptions::HaltRepartitioning e) { //TODO: rename exceptions->HaltRebalancing
      std::string gridTransferType;
      GetOStream(Warnings0, 0) <<  "Skipping permuting of "
                               << ((PorR_ == MueLu::INTERPOLATION) ? "prolongator" : "restriction")
                               << ".  No permutation is available for the following reason:"
                               << std::endl << e.what() << std::endl;
    }

    switch (PorR_) {

    case MueLu::INTERPOLATION:
      { //case scoping
        RCP<Matrix> originalP = Get< RCP<Matrix> >(coarseLevel, "P");

        if (rebalanceImporter != Teuchos::null) {
          SubFactoryMonitor m1(*this, "Rebalancing prolongator", coarseLevel);

          // P is the tranfer operator from the coarse grid to the fine grid.
          // P must transfer the data from the newly reordered coarse A to the (unchanged) fine A.
          // This means that the domain map (coarse) of P must be changed according to the new partition. The range map (fine) is kept unchanged.
          //
          // The domain map of P must match the range map of R.
          // See also note below about domain/range map of R and its implications for P.
          //
          // To change the domain map of P, P needs to be fillCompleted again with the new domain map.
          // To achieve this, P is copied into a new matrix that is not fill-completed.
          // The doImport() operation is just used here to make a copy of P: the importer is trivial and there is no data movement involved.
          // The reordering actually happens during the fillComplete() with domainMap == rebalanceImporter->getTargetMap().

          // TODO is there a better preallocation strategy?
          ArrayRCP<size_t> nnzPerRow(originalP->getNodeNumRows(), 0);
          for (size_t i=0; i<originalP->getNodeNumRows(); ++i)
            nnzPerRow[i] = originalP->getNumEntriesInLocalRow(i);

          RCP<Matrix> rebalancedP = MatrixFactory::Build(originalP->getRowMap(), nnzPerRow, Xpetra::StaticProfile);

          // Copy of P
          {
            RCP<Import> trivialImporter = ImportFactory::Build(originalP->getRowMap(), originalP->getRowMap());
            SubFactoryMonitor m2(*this, "Rebalancing prolongator -- import only", coarseLevel);
            rebalancedP->doImport(*originalP, *trivialImporter, Xpetra::INSERT);
          }

          // New domain (coarse) map, same range (fine) map
          {
            SubFactoryMonitor m2(*this, "Rebalancing prolongator -- fillComplete", coarseLevel);
            rebalancedP->fillComplete(rebalanceImporter->getTargetMap(), originalP->getRangeMap() );
          }

          Set(coarseLevel, "P", rebalancedP);

          ///////////////////////// EXPERIMENTAL
          // TODO FIXME somehow we have to transfer the striding information of the permuted domain/range maps.
          // That is probably something for an external permutation factory
          //if(originalP->IsView("stridedMaps")) rebalancedP->CreateView("stridedMaps", originalP);
          ///////////////////////// EXPERIMENTAL

        } else {
          Set(coarseLevel, "P", originalP);
        } //if (rebalanceImporter != Teuchos::null) {...} else {...}
      } //case scoping
      break;

    case MueLu::RESTRICTION:
      { //case scoping
        //TODO how do we handle implicitly transposed restriction operators?
        RCP<Matrix> originalR = Get< RCP<Matrix> >(coarseLevel, "R");
        if (rebalanceImporter != Teuchos::null) {
          SubFactoryMonitor m2(*this, "Rebalancing restriction", coarseLevel);
          RCP<SubFactoryMonitor> m1 = rcp( new SubFactoryMonitor(*this, "Rebalancing restriction -- allocate new R", coarseLevel) );

          //TODO is there a better preallocation strategy?
          /*
            Answer:  YES! :
            Count nnz per row of R
            Put this in a MV
            Use the rebalanceImporter to communicate this
            Allocate rebalancedR according to the nnz info
            Proceed as usual
          */
          //Note is to avoid that Epetra only supports vectors of type int or double, not size_t, which is unsigned.
          RCP<Xpetra::Vector<double, LO, GO> > originalNnzPerRowVec = Xpetra::VectorFactory<double, LO, GO>::Build(rebalanceImporter->getSourceMap());
          ArrayRCP<double> nnzPerRow = originalNnzPerRowVec->getDataNonConst(0);
          for (size_t i=0; i<originalR->getNodeNumRows(); ++i)
            nnzPerRow[i] = originalR->getNumEntriesInLocalRow(i);
          nnzPerRow = Teuchos::null;
          RCP<Xpetra::Vector<double, LO, GO> > permutedNnzPerRowVec = Xpetra::VectorFactory<double, LO, GO>::Build(rebalanceImporter->getTargetMap());
          permutedNnzPerRowVec->doImport(*originalNnzPerRowVec, *rebalanceImporter, Xpetra::INSERT);

          ArrayRCP<const double> tmpData = permutedNnzPerRowVec->getData(0);
          ArrayRCP<size_t> permutedNnzPerRow(permutedNnzPerRowVec->getLocalLength());
          for (size_t i=0; i<permutedNnzPerRowVec->getLocalLength(); ++i)
            permutedNnzPerRow[i] = Teuchos::as<size_t, double>(tmpData[i]);

          RCP<Matrix> rebalancedR = MatrixFactory::Build(rebalanceImporter->getTargetMap(), permutedNnzPerRow, Xpetra::StaticProfile);
          permutedNnzPerRow = Teuchos::null;
          m1 = Teuchos::null;
          m1 = rcp( new SubFactoryMonitor(*this, "Rebalancing restriction -- import", coarseLevel) );
          rebalancedR->doImport(*originalR, *rebalanceImporter, Xpetra::INSERT);

          //TODO is the following range map correct?
          //TODO RangeMap controls where a coarse grid vector's data resides after applying R
          //TODO if the targetMap of the importer is used, then we must either resumeFill or copy/fillComplete P so
          //TODO its domain map matches the range map of R.
          m1 = Teuchos::null;
          m1 = rcp( new SubFactoryMonitor(*this, "Rebalancing restriction -- fillComplete", coarseLevel) );
          rebalancedR->fillComplete( originalR->getDomainMap(), rebalanceImporter->getTargetMap() );
          //Set(coarseLevel, "newR", rebalancedR);
          m1 = Teuchos::null;
          Set(coarseLevel, "R", rebalancedR);

          ///////////////////////// EXPERIMENTAL
          // TODO FIXME somehow we have to transfer the striding information of the permuted domain/range maps.
          // That is probably something for an external permutation factory
          //if(originalR->IsView("stridedMaps")) rebalancedR->CreateView("stridedMaps", originalR);
          ///////////////////////// EXPERIMENTAL

          if (coarseLevel.IsAvailable("Coordinates")) //FIXME JJH
            {
              SubFactoryMonitor subM(*this, "Rebalancing coordinates", coarseLevel);
              //RCP<MultiVector> coords  = Get< RCP<MultiVector> >(coarseLevel, "Coordinates"); //FIXME JJH
              RCP<MultiVector> coords  = coarseLevel.Get< RCP<MultiVector> >("Coordinates"); //FIXME JJH
              RCP<MultiVector> permutedCoords  = MultiVectorFactory::Build(rebalanceImporter->getTargetMap(), coords->getNumVectors());
              permutedCoords->doImport(*coords, *rebalanceImporter, Xpetra::INSERT);
              coarseLevel.Set("Coordinates", permutedCoords); //FIXME JJH no generating factory specified
            }
          if (IsAvailable(coarseLevel, "Nullspace")) {
            SubFactoryMonitor subM(*this, "Rebalancing nullspace", coarseLevel );
            RCP<MultiVector> nullspace  = Get< RCP<MultiVector> >(coarseLevel, "Nullspace");
            RCP<MultiVector> permutedNullspace  = MultiVectorFactory::Build(rebalanceImporter->getTargetMap(), nullspace->getNumVectors());
            permutedNullspace->doImport(*nullspace, *rebalanceImporter, Xpetra::INSERT);
            coarseLevel.Set("Nullspace", permutedNullspace); //FIXME no generating factory specified
          }
        } else {
          Set(coarseLevel, "R", originalR);
          if (IsAvailable(coarseLevel, "Nullspace")) {
            RCP<MultiVector> nullspace  = Get< RCP<MultiVector> >(coarseLevel, "Nullspace");
            Set(coarseLevel, "Nullspace", nullspace);
          }
        } //if (rebalanceImporter != Teuchos::null) {...} else {...}
      } //case scoping
      break;
    } //switch

  } //Build

} // namespace MueLu

#endif // MUELU_REBALANCETRANSFERFACTORY_DEF_HPP
