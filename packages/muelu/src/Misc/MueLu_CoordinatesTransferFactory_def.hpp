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
#ifndef MUELU_COORDINATESTRANSFER_FACTORY_DEF_HPP
#define MUELU_COORDINATESTRANSFER_FACTORY_DEF_HPP

#include "Xpetra_ImportFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_MapFactory.hpp"

#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_CoordinatesTransferFactory_decl.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> CoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("Coordinates",    Teuchos::null, "Factory for coordinates generation");
    validParamList->set< RCP<const FactoryBase> >("Aggregates",     Teuchos::null, "Factory for coordinates generation");
    validParamList->set< RCP<const FactoryBase> >("CoarseMap",      Teuchos::null, "Generating factory of the coarse map");
    validParamList->set< int >                   ("write start",    -1, "first level at which coordinates should be written to file");
    validParamList->set< int >                   ("write end",      -1, "last level at which coordinates should be written to file");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(fineLevel, "Coordinates");
    Input(fineLevel, "Aggregates");
    Input(fineLevel, "CoarseMap");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level & fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    GetOStream(Runtime0) << "Transferring coordinates" << std::endl;

    RCP<Aggregates>     aggregates = Get< RCP<Aggregates> > (fineLevel, "Aggregates");
    RCP<MultiVector>    fineCoords = Get< RCP<MultiVector> >(fineLevel, "Coordinates");
    RCP<const Map>      coarseMap  = Get< RCP<const Map> >  (fineLevel, "CoarseMap");

    // coarseMap is being used to set up the domain map of tentative P, and therefore, the row map of Ac
    // Therefore, if we amalgamate coarseMap, logical nodes in the coordinates vector would correspond to
    // logical blocks in the matrix

    ArrayView<const GO> elementAList = coarseMap->getNodeElementList();
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
    RCP<MultiVector> coarseCoords   = MultiVectorFactory::Build(coarseCoordMap, fineCoords->getNumVectors());

    // Create overlapped fine coordinates to reduce global communication
    RCP<MultiVector> ghostedCoords = fineCoords;
    if (aggregates->AggregatesCrossProcessors()) {
      RCP<const Map>    nonUniqueMap = aggregates->GetMap();
      RCP<const Import> importer     = ImportFactory::Build(uniqueMap, nonUniqueMap);

      ghostedCoords = MultiVectorFactory::Build(nonUniqueMap, fineCoords->getNumVectors());
      ghostedCoords->doImport(*fineCoords, *importer, Xpetra::INSERT);
    }

    // Get some info about aggregates
    int                         myPID        = uniqueMap->getComm()->getRank();
    LO                          numAggs      = aggregates->GetNumAggregates();
    ArrayRCP<LO>                aggSizes     = aggregates->ComputeAggregateSizes(true,true);
    const ArrayRCP<const LO>    vertex2AggID = aggregates->GetVertex2AggId()->getData(0);
    const ArrayRCP<const LO>    procWinner   = aggregates->GetProcWinner()->getData(0);

    // Fill in coarse coordinates
    for (size_t j = 0; j < fineCoords->getNumVectors(); j++) {
      ArrayRCP<const Scalar> fineCoordsData = ghostedCoords->getData(j);
      ArrayRCP<Scalar>     coarseCoordsData = coarseCoords->getDataNonConst(j);

      for (LO lnode = 0; lnode < vertex2AggID.size(); lnode++)
        if (procWinner[lnode] == myPID)
          coarseCoordsData[vertex2AggID[lnode]] += fineCoordsData[lnode];

      for (LO agg = 0; agg < numAggs; agg++)
        coarseCoordsData[agg] /= aggSizes[agg];
    }

    Set<RCP<MultiVector> >(coarseLevel, "Coordinates", coarseCoords);

    const ParameterList& pL = GetParameterList();
    int writeStart = pL.get<int>("write start"), writeEnd = pL.get<int>("write end");
    if (writeStart == 0 && fineLevel.GetLevelID() == 0 && writeStart <= writeEnd) {
      std::ostringstream buf;
      buf << fineLevel.GetLevelID();
      std::string fileName = "coordinates_before_rebalance_level_" + buf.str() + ".m";
      Utils::Write(fileName,*fineCoords);
    }
    if (writeStart <= coarseLevel.GetLevelID() && coarseLevel.GetLevelID() <= writeEnd) {
      std::ostringstream buf;
      buf << coarseLevel.GetLevelID();
      std::string fileName = "coordinates_before_rebalance_level_" + buf.str() + ".m";
      Utils::Write(fileName,*coarseCoords);
    }
  }

} // namespace MueLu

#endif // MUELU_COORDINATESTRANSFER_FACTORY_DEF_HPP
