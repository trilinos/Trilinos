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
#ifndef MUELU_MULTIVECTORAGGREGATIONTRANSFER_FACTORY_DEF_HPP
#define MUELU_MULTIVECTORAGGREGATIONTRANSFER_FACTORY_DEF_HPP

#include "Xpetra_ImportFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_IO.hpp"

#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_MultiVectorAggregationTransferFactory_decl.hpp"
//#include "MueLu_Utilities.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> MultiVectorAggregationTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< std::string >           ("Vector name",    "undefined", "Name of the vector that will be transferred on the coarse grid (level key)"); // TODO: how to set a validator without default value?
    validParamList->set< RCP<const FactoryBase> >("Vector factory", Teuchos::null, "Factory of the vector");
    validParamList->set< RCP<const FactoryBase> >("Aggregates",     Teuchos::null, "Factory of the aggregates");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVectorAggregationTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    const ParameterList & pL = GetParameterList();
    std::string vectorName   = pL.get<std::string>("Vector name");

    fineLevel.DeclareInput(vectorName, GetFactory("Vector factory").get(), this);
    Input(fineLevel, "Aggregates");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVectorAggregationTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level & fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    const ParameterList & pL = GetParameterList();
    std::string vectorName   = pL.get<std::string>("Vector name");

    RCP<MultiVector> fineVector = fineLevel.Get< RCP<MultiVector> >(vectorName, GetFactory("Vector factory").get());
    RCP<Matrix>      transferOp = Get<RCP<Matrix> >(coarseLevel, "R");

    RCP<MultiVector> coarseVector = MultiVectorFactory::Build(transferOp->getRangeMap(), fineVector->getNumVectors());
    GetOStream(Runtime0) << "Transferring multivector \"" << vectorName << "\"" << std::endl;

    RCP<MultiVector> onesVector = MultiVectorFactory::Build(transferOp->getDomainMap(), 1);
    onesVector->putScalar(Teuchos::ScalarTraits<Scalar>::one());
    RCP<MultiVector> rowSumVector = MultiVectorFactory::Build(transferOp->getRangeMap(), 1);
    transferOp->apply(*onesVector, *rowSumVector);
    transferOp->apply(*fineVector, *coarseVector);

    if (vectorName == "Coordinates")
      TEUCHOS_TEST_FOR_EXCEPTION(true,Exceptions::RuntimeError,"Use CoordinatesTransferFactory to transfer coordinates instead of MultiVectorAggregationTransferFactory.");

    Set<RCP<MultiVector> >(coarseLevel, vectorName, coarseVector);

  } // Build

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayRCP<Scalar> MultiVectorAggregationTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::expandCoordinates(ArrayRCP<SC> coordinates, LocalOrdinal blksize) {
    if (blksize == 1)
      return coordinates;

    ArrayRCP<SC> expandCoord(coordinates.size()*blksize); //TODO: how to avoid automatic initialization of the vector? using arcp()?

    for(int i=0; i<coordinates.size(); i++) {
      for(int j=0; j< blksize; j++) {
        expandCoord[i*blksize + j] = coordinates[i];
      }
    }
    return expandCoord;

  } // expandCoordinates

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVectorAggregationTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level & fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    using xdMV = Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LO,GO,NO>;

    GetOStream(Runtime0) << "Transferring coordinates" << std::endl;

    int numDimensions;
    RCP<xdMV> coarseCoords;
    RCP<xdMV> fineCoords;
    Array<GO> gCoarseNodesPerDir;
    Array<LO> lCoarseNodesPerDir;

    const ParameterList& pL = GetParameterList();

    if(pL.get<bool>("hybrid aggregation") == true) {
      std::string regionType = Get<std::string>(fineLevel,"aggregationRegionTypeCoarse");
      numDimensions          = Get<int>(fineLevel, "numDimensions");
      lCoarseNodesPerDir     = Get<Array<LO> >(fineLevel, "lCoarseNodesPerDim");
      Set<std::string>(coarseLevel, "aggregationRegionType", regionType);
      Set<int>        (coarseLevel, "numDimensions",         numDimensions);
      Set<Array<LO> > (coarseLevel, "lNodesPerDim",          lCoarseNodesPerDir);

      if((pL.get<bool>("interface aggregation") == true) && (regionType == "uncoupled")) {
        Array<LO> coarseInterfacesDimensions = Get<Array<LO> >(fineLevel, "coarseInterfacesDimensions");
        Array<LO> nodeOnCoarseInterface      = Get<Array<LO> >(fineLevel, "nodeOnCoarseInterface");
        Set<Array<LO> >(coarseLevel, "interfacesDimensions", coarseInterfacesDimensions);
        Set<Array<LO> >(coarseLevel, "nodeOnInterface",      nodeOnCoarseInterface);
      }

    } else if(pL.get<bool>("structured aggregation") == true) {
      if(pL.get<bool>("aggregation coupled") == true) {
        gCoarseNodesPerDir = Get<Array<GO> >(fineLevel, "gCoarseNodesPerDim");
        Set<Array<GO> >(coarseLevel, "gNodesPerDim", gCoarseNodesPerDir);
      }
      lCoarseNodesPerDir = Get<Array<LO> >(fineLevel, "lCoarseNodesPerDim");
      Set<Array<LO> >(coarseLevel, "lNodesPerDim", lCoarseNodesPerDir);
      numDimensions = Get<int>(fineLevel, "numDimensions");
      Set<int>(coarseLevel, "numDimensions", numDimensions);

    } else if(pL.get<bool>("Geometric") == true) {
      coarseCoords       = Get<RCP<xdMV> >(coarseLevel, "coarseCoordinates");
      gCoarseNodesPerDir = Get<Array<GO> >(coarseLevel, "gCoarseNodesPerDim");
      lCoarseNodesPerDir = Get<Array<LO> >(coarseLevel, "lCoarseNodesPerDim");
      Set<Array<GO> >(coarseLevel, "gNodesPerDim", gCoarseNodesPerDir);
      Set<Array<LO> >(coarseLevel, "lNodesPerDim", lCoarseNodesPerDir);

      Set<RCP<xdMV> >(coarseLevel, "Coordinates", coarseCoords);

    } else {
      if (coarseLevel.IsAvailable("Coordinates", this)) {
        GetOStream(Runtime0) << "Reusing coordinates" << std::endl;
        return;
      }

      RCP<Aggregates>     aggregates = Get< RCP<Aggregates> > (fineLevel, "Aggregates");
      fineCoords                     = Get< RCP<xdMV> >(fineLevel, "Coordinates");
      RCP<const Map>      coarseMap  = Get< RCP<const Map> >  (fineLevel, "CoarseMap");

      // coarseMap is being used to set up the domain map of tentative P, and therefore, the row map of Ac
      // Therefore, if we amalgamate coarseMap, logical nodes in the coordinates vector would correspond to
      // logical blocks in the matrix

      ArrayView<const GO> elementAList = coarseMap->getLocalElementList();

      LO                  blkSize      = 1;
      if (rcp_dynamic_cast<const StridedMap>(coarseMap) != Teuchos::null)
        blkSize = rcp_dynamic_cast<const StridedMap>(coarseMap)->getFixedBlockSize();

      GO                  indexBase    = coarseMap->getIndexBase();
      size_t              numElements  = elementAList.size() / blkSize;
      Array<GO>           elementList(numElements);

      // Amalgamate the map
      for (LO i = 0; i < Teuchos::as<LO>(numElements); i++)
        elementList[i] = (elementAList[i*blkSize]-indexBase)/blkSize + indexBase;

      RCP<const Map>   uniqueMap      = fineCoords->getMap();
      RCP<const Map>   coarseCoordMap = MapFactory        ::Build(coarseMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), elementList, indexBase, coarseMap->getComm());
      coarseCoords   = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LO,GO,NO>::Build(coarseCoordMap, fineCoords->getNumVectors());

      // Create overlapped fine coordinates to reduce global communication
      RCP<xdMV> ghostedCoords = fineCoords;
      if (aggregates->AggregatesCrossProcessors()) {
        RCP<const Map>    nonUniqueMap = aggregates->GetMap();
        RCP<const Import> importer     = ImportFactory::Build(uniqueMap, nonUniqueMap);

        ghostedCoords = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LO,GO,NO>::Build(nonUniqueMap, fineCoords->getNumVectors());
        ghostedCoords->doImport(*fineCoords, *importer, Xpetra::INSERT);
      }

      // Get some info about aggregates
      int                         myPID        = uniqueMap->getComm()->getRank();
      LO                          numAggs      = aggregates->GetNumAggregates();
      ArrayRCP<LO>                aggSizes     = aggregates->ComputeAggregateSizes();
      const ArrayRCP<const LO>    vertex2AggID = aggregates->GetVertex2AggId()->getData(0);
      const ArrayRCP<const LO>    procWinner   = aggregates->GetProcWinner()->getData(0);

      // Fill in coarse coordinates
      for (size_t j = 0; j < fineCoords->getNumVectors(); j++) {
        ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> fineCoordsData = ghostedCoords->getData(j);
        ArrayRCP<typename Teuchos::ScalarTraits<Scalar>::coordinateType>     coarseCoordsData = coarseCoords->getDataNonConst(j);

        for (LO lnode = 0; lnode < vertex2AggID.size(); lnode++) {
          if (procWinner[lnode] == myPID &&
              lnode < vertex2AggID.size() &&
              lnode < fineCoordsData.size() && // TAW do not access off-processor coordinates
              vertex2AggID[lnode] < coarseCoordsData.size() &&
              Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::coordinateType>::isnaninf(fineCoordsData[lnode]) == false) {
            coarseCoordsData[vertex2AggID[lnode]] += fineCoordsData[lnode];
          }
        }
        for (LO agg = 0; agg < numAggs; agg++) {
          coarseCoordsData[agg] /= aggSizes[agg];
        }
      }

      Set<RCP<xdMV> >(coarseLevel, "Coordinates", coarseCoords);
    } // if pL.get<bool>("Geometric") == true

    int writeStart = pL.get<int>("write start"), writeEnd = pL.get<int>("write end");
    if (writeStart == 0 && fineLevel.GetLevelID() == 0 && writeStart <= writeEnd) {
      std::ostringstream buf;
      buf << fineLevel.GetLevelID();
      std::string fileName = "coordinates_before_rebalance_level_" + buf.str() + ".m";
      Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LO,GO,NO>::Write(fileName,*fineCoords);
    }
    if (writeStart <= coarseLevel.GetLevelID() && coarseLevel.GetLevelID() <= writeEnd) {
      std::ostringstream buf;
      buf << coarseLevel.GetLevelID();
      std::string fileName = "coordinates_before_rebalance_level_" + buf.str() + ".m";
      Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LO,GO,NO>::Write(fileName,*coarseCoords);
    }
  }

} // namespace MueLu

#endif // MUELU_MULTIVECTORAGGREGATIONTRANSFER_FACTORY_DEF_HPP
