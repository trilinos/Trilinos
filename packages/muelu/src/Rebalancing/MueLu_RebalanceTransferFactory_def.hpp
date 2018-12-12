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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_REBALANCETRANSFERFACTORY_DEF_HPP
#define MUELU_REBALANCETRANSFERFACTORY_DEF_HPP

#include <Teuchos_Tuple.hpp>

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_Vector.hpp"
#include "Xpetra_VectorFactory.hpp"
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_RebalanceTransferFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"

namespace MueLu {

 template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
 RCP<const ParameterList> RebalanceTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("repartition: rebalance P and R");
    SET_VALID_ENTRY("repartition: rebalance Nullspace");
    SET_VALID_ENTRY("transpose: use implicit");
    SET_VALID_ENTRY("repartition: use subcommunicators");
#undef  SET_VALID_ENTRY

    {
      typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;
      RCP<validatorType> typeValidator = rcp (new validatorType(Teuchos::tuple<std::string>("Interpolation", "Restriction"), "type"));
      validParamList->set("type", "Interpolation", "Type of the transfer operator that need to be rebalanced (Interpolation or Restriction)", typeValidator);
    }

    validParamList->set< RCP<const FactoryBase> >("P",                   null, "Factory of the prolongation operator that need to be rebalanced (only used if type=Interpolation)");
    validParamList->set< RCP<const FactoryBase> >("R",                   null, "Factory of the restriction operator that need to be rebalanced (only used if type=Restriction)");
    validParamList->set< RCP<const FactoryBase> >("Nullspace",           null, "Factory of the nullspace that need to be rebalanced (only used if type=Interpolation)");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",         null, "Factory of the coordinates that need to be rebalanced (only used if type=Interpolation)");
    validParamList->set< RCP<const FactoryBase> >("Importer",            null, "Factory of the importer object used for the rebalancing");
    validParamList->set< int >                   ("write start",           -1, "First level at which coordinates should be written to file");
    validParamList->set< int >                   ("write end",             -1, "Last level at which coordinates should be written to file");

    // TODO validation: "P" parameter valid only for type="Interpolation" and "R" valid only for type="Restriction". Like so:
    // if (paramList.isEntry("type") && paramList.get("type) == "Interpolation) {
    //     validParamList->set< RCP<const FactoryBase> >("P",              Teuchos::null, "Factory of the prolongation operator that need to be rebalanced (only used if type=Interpolation)");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RebalanceTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    const ParameterList& pL = GetParameterList();

    if (pL.get<std::string>("type") == "Interpolation") {
      Input(coarseLevel, "P");
      if (pL.get<bool>("repartition: rebalance Nullspace"))
        Input(coarseLevel, "Nullspace");
      if (pL.get< RCP<const FactoryBase> >("Coordinates") != Teuchos::null)
        Input(coarseLevel, "Coordinates");

    } else {
      if (pL.get<bool>("transpose: use implicit") == false)
        Input(coarseLevel, "R");
    }

    Input(coarseLevel, "Importer");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RebalanceTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);
    typedef Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> xdMV;

    const ParameterList& pL = GetParameterList();

    int implicit   = !pL.get<bool>("repartition: rebalance P and R");
    int writeStart = pL.get<int> ("write start");
    int writeEnd   = pL.get<int> ("write end");

    if (writeStart == 0 && fineLevel.GetLevelID() == 0 && writeStart <= writeEnd && IsAvailable(fineLevel, "Coordinates")) {
      std::string fileName = "coordinates_level_0.m";
      RCP<xdMV> fineCoords = fineLevel.Get< RCP<xdMV> >("Coordinates");
      if (fineCoords != Teuchos::null)
        Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LO,GO,NO>::Write(fileName, *fineCoords);
    }

    RCP<const Import> importer = Get<RCP<const Import> >(coarseLevel, "Importer");
    if (implicit) {
      // Save the importer, we'll need it for solve
      coarseLevel.Set("Importer", importer, NoFactory::get());
    }

    RCP<ParameterList> params = rcp(new ParameterList());
    if (IsPrint(Statistics2)) {
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo",          true);
    }

    std::string transferType = pL.get<std::string>("type");
    if (transferType == "Interpolation") {
      RCP<Matrix> originalP = Get< RCP<Matrix> >(coarseLevel, "P");

      {
        // This line must be after the Get call
        SubFactoryMonitor m1(*this, "Rebalancing prolongator", coarseLevel);

        if (implicit || importer.is_null()) {
          GetOStream(Runtime0) << "Using original prolongator" << std::endl;
          Set(coarseLevel, "P", originalP);

        } else {
          // P is the transfer operator from the coarse grid to the fine grid.
          // P must transfer the data from the newly reordered coarse A to the
          // (unchanged) fine A.  This means that the domain map (coarse) of P
          // must be changed according to the new partition. The range map
          // (fine) is kept unchanged.
          //
          // The domain map of P must match the range map of R.  See also note
          // below about domain/range map of R and its implications for P.
          //
          // To change the domain map of P, P needs to be fillCompleted again
          // with the new domain map.  To achieve this, P is copied into a new
          // matrix that is not fill-completed.  The doImport() operation is
          // just used here to make a copy of P: the importer is trivial and
          // there is no data movement involved.  The reordering actually
          // happens during the fillComplete() with domainMap == importer->getTargetMap().
          RCP<Matrix> rebalancedP = originalP;
          RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(originalP);
          TEUCHOS_TEST_FOR_EXCEPTION(crsOp == Teuchos::null, Exceptions::BadCast, "Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

          RCP<CrsMatrix> rebalancedP2 = crsOp->getCrsMatrix();
          TEUCHOS_TEST_FOR_EXCEPTION(rebalancedP2 == Teuchos::null, std::runtime_error, "Xpetra::CrsMatrixWrap doesn't have a CrsMatrix");

          {
            SubFactoryMonitor subM(*this, "Rebalancing prolongator -- fast map replacement", coarseLevel);

            RCP<const Import> newImporter;
            {
              SubFactoryMonitor(*this, "Import construction", coarseLevel);
              newImporter = ImportFactory::Build(importer->getTargetMap(), rebalancedP->getColMap());
            }
            rebalancedP2->replaceDomainMapAndImporter(importer->getTargetMap(), newImporter);
          }

          ///////////////////////// EXPERIMENTAL
          // TODO FIXME somehow we have to transfer the striding information of the permuted domain/range maps.
          // That is probably something for an external permutation factory
          //   if (originalP->IsView("stridedMaps"))
          //     rebalancedP->CreateView("stridedMaps", originalP);
          ///////////////////////// EXPERIMENTAL

          Set(coarseLevel, "P", rebalancedP);

          if (IsPrint(Statistics2))
            GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*rebalancedP, "P (rebalanced)", params);
        }
      }

      if (importer.is_null()) {
        if (IsAvailable(coarseLevel, "Nullspace"))
          Set(coarseLevel, "Nullspace", Get<RCP<MultiVector> >(coarseLevel, "Nullspace"));

        if (pL.isParameter("Coordinates") && pL.get< RCP<const FactoryBase> >("Coordinates") != Teuchos::null)
          if (IsAvailable(coarseLevel, "Coordinates"))
            Set(coarseLevel, "Coordinates", Get< RCP<xdMV> >(coarseLevel, "Coordinates"));

        return;
      }

      if (pL.isParameter("Coordinates") &&
          pL.get< RCP<const FactoryBase> >("Coordinates") != Teuchos::null &&
          IsAvailable(coarseLevel, "Coordinates")) {
        RCP<xdMV> coords = Get<RCP<xdMV> >(coarseLevel, "Coordinates");

        // This line must be after the Get call
        SubFactoryMonitor subM(*this, "Rebalancing coordinates", coarseLevel);

        LO nodeNumElts = coords->getMap()->getNodeNumElements();

        // If a process has no matrix rows, then we can't calculate blocksize using the formula below.
        LO myBlkSize = 0, blkSize = 0;
        if (nodeNumElts > 0)
          myBlkSize = importer->getSourceMap()->getNodeNumElements() / nodeNumElts;
        MueLu_maxAll(coords->getMap()->getComm(), myBlkSize, blkSize);

        RCP<const Import> coordImporter;
        if (blkSize == 1) {
          coordImporter = importer;

        } else {
          // NOTE: there is an implicit assumption here: we assume that dof any node are enumerated consequently
          // Proper fix would require using decomposition similar to how we construct importer in the
          // RepartitionFactory
          RCP<const Map> origMap   = coords->getMap();
          GO             indexBase = origMap->getIndexBase();

          ArrayView<const GO> OEntries   = importer->getTargetMap()->getNodeElementList();
          LO                  numEntries = OEntries.size()/blkSize;
          ArrayRCP<GO> Entries(numEntries);
          for (LO i = 0; i < numEntries; i++)
            Entries[i] = (OEntries[i*blkSize]-indexBase)/blkSize + indexBase;

          RCP<const Map> targetMap = MapFactory::Build(origMap->lib(), origMap->getGlobalNumElements(), Entries(), indexBase, origMap->getComm());
          coordImporter = ImportFactory::Build(origMap, targetMap);
        }

        RCP<xdMV> permutedCoords  = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LO,GO,NO>::Build(coordImporter->getTargetMap(), coords->getNumVectors());
        permutedCoords->doImport(*coords, *coordImporter, Xpetra::INSERT);

        if (pL.isParameter("repartition: use subcommunicators") == true && pL.get<bool>("repartition: use subcommunicators") == true)
          permutedCoords->replaceMap(permutedCoords->getMap()->removeEmptyProcesses());

        if (permutedCoords->getMap() == Teuchos::null)
          permutedCoords = Teuchos::null;

        Set(coarseLevel, "Coordinates", permutedCoords);

        std::string fileName = "rebalanced_coordinates_level_" + toString(coarseLevel.GetLevelID()) + ".m";
        if (writeStart <= coarseLevel.GetLevelID() && coarseLevel.GetLevelID() <= writeEnd && permutedCoords->getMap() != Teuchos::null)
          Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LO,GO,NO>::Write(fileName, *permutedCoords);
      }

      if (IsAvailable(coarseLevel, "Nullspace")) {
        RCP<MultiVector> nullspace = Get< RCP<MultiVector> >(coarseLevel, "Nullspace");

        // This line must be after the Get call
        SubFactoryMonitor subM(*this, "Rebalancing nullspace", coarseLevel);

        RCP<MultiVector> permutedNullspace = MultiVectorFactory::Build(importer->getTargetMap(), nullspace->getNumVectors());
        permutedNullspace->doImport(*nullspace, *importer, Xpetra::INSERT);

        if (pL.get<bool>("repartition: use subcommunicators") == true)
          permutedNullspace->replaceMap(permutedNullspace->getMap()->removeEmptyProcesses());

        if (permutedNullspace->getMap() == Teuchos::null)
          permutedNullspace = Teuchos::null;

        Set(coarseLevel, "Nullspace", permutedNullspace);
      }

    } else {
      if (pL.get<bool>("transpose: use implicit") == false) {
        RCP<Matrix> originalR = Get< RCP<Matrix> >(coarseLevel, "R");

        SubFactoryMonitor m2(*this, "Rebalancing restrictor", coarseLevel);

        if (implicit || importer.is_null()) {
          GetOStream(Runtime0) << "Using original restrictor" << std::endl;
          Set(coarseLevel, "R", originalR);

        } else {
          RCP<Matrix> rebalancedR;
          {
            SubFactoryMonitor subM(*this, "Rebalancing restriction -- fusedImport", coarseLevel);

            RCP<Map> dummy;         // meaning: use originalR's domain map.
            Teuchos::ParameterList listLabel;
            listLabel.set("Timer Label","MueLu::RebalanceR-" + Teuchos::toString(coarseLevel.GetLevelID()));
            rebalancedR = MatrixFactory::Build(originalR, *importer, dummy, importer->getTargetMap(),Teuchos::rcp(&listLabel,false));
          }
          Set(coarseLevel, "R", rebalancedR);

          ///////////////////////// EXPERIMENTAL
          // TODO FIXME somehow we have to transfer the striding information of the permuted domain/range maps.
          // That is probably something for an external permutation factory
          // if (originalR->IsView("stridedMaps"))
          //   rebalancedR->CreateView("stridedMaps", originalR);
          ///////////////////////// EXPERIMENTAL

          if (IsPrint(Statistics2))
            GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*rebalancedR, "R (rebalanced)", params);
        }
      }
    }
  }

} // namespace MueLu

#endif // MUELU_REBALANCETRANSFERFACTORY_DEF_HPP
