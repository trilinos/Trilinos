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
#ifndef MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP
#define MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP

#include "Xpetra_Vector.hpp"
#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>

#include "MueLu_PermutedTransferFactory_decl.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
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
  void PermutedTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    if (PorR_ == MueLu::INTERPOLATION) {
      GetOStream(Warnings0, 0) <<  "Jamming A into Level " << coarseLevel.GetLevelID() << " w/ generating factory "
                               << this << std::endl;
      RCP<Matrix> A = Get< RCP<Matrix> >(coarseLevel, "A");
      Set(coarseLevel, "A", A);
    }

    RCP<const Import> permImporter;
    try {
      permImporter = Get< RCP<const Import> >(coarseLevel, "Importer");
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
          RCP<Matrix> originalP = Get< RCP<Matrix> >(coarseLevel, "P");
          if (permImporter != Teuchos::null) {
            SubFactoryMonitor m1(*this, "Rebalancing prolongator", coarseLevel.GetLevelID());
            //TODO see note below about domain/range map of R and its implications for P.

            // Now copy P so that we can give it a domain map that matches the range map of R.
            //TODO is there a better preallocation strategy?
            ArrayRCP<size_t> nnzPerRow(originalP->getNodeNumRows(), 0);
            for (size_t i=0; i<originalP->getNodeNumRows(); ++i)
              nnzPerRow[i] = originalP->getNumEntriesInLocalRow(i);
            RCP<Matrix> permutedP = MatrixFactory::Build(originalP->getRowMap(), nnzPerRow, Xpetra::StaticProfile);
            //P needs to be fillCompleted again.  To achieve this, I just copy P using the trivial importer.
            RCP<Import> trivialImporter = ImportFactory::Build( originalP->getRowMap(), originalP->getRowMap());
            RCP<CrsMatrixWrap> crsOp = rcp_dynamic_cast<CrsMatrixWrap>(permutedP);
            RCP<CrsMatrix> crsMtx = crsOp->getCrsMatrix();
            RCP<CrsMatrixWrap> origOp = rcp_dynamic_cast<CrsMatrixWrap>(originalP);
            RCP<CrsMatrix> origMtx = origOp->getCrsMatrix();
            RCP<SubFactoryMonitor> m2 = rcp( new SubFactoryMonitor(*this, "Rebalancing prolongator -- import only",
                                                                   coarseLevel.GetLevelID()) );
            crsMtx->doImport(*origMtx, *trivialImporter, Xpetra::INSERT);
            crsMtx = Teuchos::null;
            //new domain (coarse) map, same range (fine) map
            m2 = rcp( new SubFactoryMonitor(*this, "Rebalancing prolongator -- fillComplete", coarseLevel.GetLevelID()) );
            permutedP->fillComplete(permImporter->getTargetMap(), originalP->getRangeMap() );
            m2 = Teuchos::null;
            originalP = Teuchos::null;
            Set(coarseLevel, "P", permutedP);

            ///////////////////////// EXPERIMENTAL
            // TODO FIXME somehow we have to transfer the striding information of the permuted domain/range maps.
            // That is probably something for an external permutation factory
            //if(originalP->IsView("stridedMaps")) permutedP->CreateView("stridedMaps", originalP);
            ///////////////////////// EXPERIMENTAL

          } else {
            Set(coarseLevel, "P", originalP);
          } //if (permImporter != Teuchos::null) {...} else {...}
        } //case scoping
        break;

      case MueLu::RESTRICTION:
        { //case scoping
          //TODO how do we handle implicitly transposed restriction operators?
          RCP<Matrix> originalR = Get< RCP<Matrix> >(coarseLevel, "R");
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
            //Note is to avoid that Epetra only supports vectors of type int or double, not size_t, which is unsigned.
            RCP<Xpetra::Vector<double, LO, GO> > originalNnzPerRowVec = Xpetra::VectorFactory<double, LO, GO>::Build(permImporter->getSourceMap());
            ArrayRCP<double> nnzPerRow = originalNnzPerRowVec->getDataNonConst(0);
            for (size_t i=0; i<originalR->getNodeNumRows(); ++i)
              nnzPerRow[i] = originalR->getNumEntriesInLocalRow(i);
            nnzPerRow = Teuchos::null;
            RCP<Xpetra::Vector<double, LO, GO> > permutedNnzPerRowVec = Xpetra::VectorFactory<double, LO, GO>::Build(permImporter->getTargetMap());
            permutedNnzPerRowVec->doImport(*originalNnzPerRowVec, *permImporter, Xpetra::INSERT);

            ArrayRCP<const double> tmpData = permutedNnzPerRowVec->getData(0);
            ArrayRCP<size_t> permutedNnzPerRow(permutedNnzPerRowVec->getLocalLength());
            for (size_t i=0; i<permutedNnzPerRowVec->getLocalLength(); ++i)
              permutedNnzPerRow[i] = Teuchos::as<size_t, double>(tmpData[i]);

            RCP<Matrix> permutedR = MatrixFactory::Build(permImporter->getTargetMap(), permutedNnzPerRow, Xpetra::StaticProfile);
            permutedNnzPerRow = Teuchos::null;
            RCP<CrsMatrixWrap> crsOp = rcp_dynamic_cast<CrsMatrixWrap>(permutedR);
            RCP<CrsMatrix> crsMtx = crsOp->getCrsMatrix();
            RCP<CrsMatrixWrap> origOp = rcp_dynamic_cast<CrsMatrixWrap>(originalR);
            RCP<CrsMatrix> origMtx = origOp->getCrsMatrix();
            m1 = Teuchos::null;
            m1 = rcp( new SubFactoryMonitor(*this, "Rebalancing restriction -- import", coarseLevel.GetLevelID()) );
            crsMtx->doImport(*origMtx, *permImporter, Xpetra::INSERT);
            crsMtx = Teuchos::null;
            //TODO is the following range map correct?
            //TODO RangeMap controls where a coarse grid vector's data resides after applying R
            //TODO if the targetMap of the importer is used, then we must either resumeFill or copy/fillComplete P so
            //TODO its domain map matches the range map of R.
            m1 = Teuchos::null;
            m1 = rcp( new SubFactoryMonitor(*this, "Rebalancing restriction -- fillComplete", coarseLevel.GetLevelID()) );
            permutedR->fillComplete( originalR->getDomainMap() , permImporter->getTargetMap() );
            //Set(coarseLevel, "newR", permutedR);
            m1 = Teuchos::null;
            Set(coarseLevel, "R", permutedR);

            ///////////////////////// EXPERIMENTAL
            // TODO FIXME somehow we have to transfer the striding information of the permuted domain/range maps.
            // That is probably something for an external permutation factory
            //if(originalR->IsView("stridedMaps")) permutedR->CreateView("stridedMaps", originalR);
            ///////////////////////// EXPERIMENTAL

            if (coarseLevel.IsAvailable("Coordinates")) //FIXME JJH
            {
              SubFactoryMonitor subM(*this, "Rebalancing coordinates", coarseLevel.GetLevelID());
              //RCP<MultiVector> coords  = Get< RCP<MultiVector> >(coarseLevel, "Coordinates"); //FIXME JJH
              RCP<MultiVector> coords  = coarseLevel.Get< RCP<MultiVector> >("Coordinates"); //FIXME JJH
              RCP<MultiVector> permutedCoords  = MultiVectorFactory::Build(permImporter->getTargetMap(), coords->getNumVectors());
              permutedCoords->doImport(*coords, *permImporter, Xpetra::INSERT);
              coarseLevel.Set("Coordinates", permutedCoords); //FIXME JJH no generating factory specified
            }
            if (IsAvailable(coarseLevel, "Nullspace")) {
              SubFactoryMonitor subM(*this, "Rebalancing nullspace", coarseLevel.GetLevelID());
              RCP<MultiVector> nullspace  = Get< RCP<MultiVector> >(coarseLevel, "Nullspace");
              RCP<MultiVector> permutedNullspace  = MultiVectorFactory::Build(permImporter->getTargetMap(), nullspace->getNumVectors());
              permutedNullspace->doImport(*nullspace, *permImporter, Xpetra::INSERT);
              coarseLevel.Set("Nullspace", permutedNullspace); //FIXME no generating factory specified
            }
          } else {
            Set(coarseLevel, "R", originalR);
            if (IsAvailable(coarseLevel, "Nullspace")) {
              RCP<MultiVector> nullspace  = Get< RCP<MultiVector> >(coarseLevel, "Nullspace");
              Set(coarseLevel, "Nullspace", nullspace);
            }
          } //if (permImporter != Teuchos::null) {...} else {...}
        } //case scoping
        break;
    } //switch

  } //Build

} // namespace MueLu

#define MUELU_PERMUTEDTRANSFER_FACTORY_SHORT
#endif // MUELU_PERMUTEDTRANSFER_FACTORY_DEF_HPP
