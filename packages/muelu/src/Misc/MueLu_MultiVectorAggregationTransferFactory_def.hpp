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
    validParamList->set< RCP<const FactoryBase> >("CoarseMap",     Teuchos::null, "Factory of the coarse map");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVectorAggregationTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    const ParameterList & pL = GetParameterList();
    std::string vectorName   = pL.get<std::string>("Vector name");

    fineLevel.DeclareInput(vectorName, GetFactory("Vector factory").get(), this);
    Input(fineLevel, "Aggregates");
    Input(fineLevel, "CoarseMap");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVectorAggregationTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level & fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    const ParameterList & pL = GetParameterList();
    std::string vectorName   = pL.get<std::string>("Vector name");

    if (vectorName == "Coordinates")
      TEUCHOS_TEST_FOR_EXCEPTION(true,Exceptions::RuntimeError,"Use CoordinatesTransferFactory to transfer coordinates instead of MultiVectorAggregationTransferFactory.");

    int numDimensions;
    RCP<MultiVector> coarseVector;
    RCP<MultiVector> fineVector;

    if (coarseLevel.IsAvailable(vectorName, this)) {
      GetOStream(Runtime0) << "Reusing " << vectorName << std::endl;
      return;
    }

    RCP<Aggregates> aggregates = Get< RCP<Aggregates> > (fineLevel, "Aggregates");
    fineVector = fineLevel.Get<RCP<MultiVector> >(vectorName);

    if (aggregates->AggregatesCrossProcessors())
      TEUCHOS_TEST_FOR_EXCEPTION(true,Exceptions::RuntimeError,"MultiVectorAggregationTransferFactory does not support aggregates that cross processors!");

    // if the coarse map is available, use it. otherwise, construct the map based on aggregates
    RCP<const Map> coarseMap = Get< RCP<const Map> > (fineLevel, "CoarseMap");
    // if (fineLevel.IsAvailable("CoarseMap")) {
    //   coarseMap = Get< RCP<const Map> >(fineLevel, "CoarseMap");
    // } else {
    //   std::cout << "Number of local aggregates=" << aggregates->GetNumAggregates() << std::endl;
    //   coarseMap = aggregates->GetMap();
    //   std::cout << "Number of local map entries=" << coarseMap->getLocalNumElements() << std::endl;
    //   RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    //   coarseMap->describe(*fancy,Teuchos::VERB_EXTREME);
    // }

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

    RCP<const Map>   uniqueMap      = fineVector->getMap();
    RCP<const Map>   coarseMVMap = MapFactory        ::Build(coarseMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), elementList, indexBase, coarseMap->getComm());
    coarseVector   = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(coarseMVMap, fineVector->getNumVectors());

    // Create overlapped fine vector to reduce global communication
    RCP<MultiVector> ghostedVector = fineVector;
    if (aggregates->AggregatesCrossProcessors()) {
      RCP<const Map>    nonUniqueMap = aggregates->GetMap();
      RCP<const Import> importer     = ImportFactory::Build(uniqueMap, nonUniqueMap);

      ghostedVector = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(nonUniqueMap, fineVector->getNumVectors());
      ghostedVector->doImport(*fineVector, *importer, Xpetra::INSERT);
    }

    // Get some info about aggregates
    int                         myPID        = uniqueMap->getComm()->getRank();
    LO                          numAggs      = aggregates->GetNumAggregates();
    ArrayRCP<LO>                aggSizes     = aggregates->ComputeAggregateSizes();
    const ArrayRCP<const LO>    vertex2AggID = aggregates->GetVertex2AggId()->getData(0);
    const ArrayRCP<const LO>    procWinner   = aggregates->GetProcWinner()->getData(0);

    // Fill in coarse vector
    for (size_t j = 0; j < fineVector->getNumVectors(); j++) {
      ArrayRCP<const SC> fineVectorData = ghostedVector->getData(j);
      ArrayRCP<SC>     coarseVectorData = coarseVector->getDataNonConst(j);

      for (LO lnode = 0; lnode < vertex2AggID.size(); lnode++) {
        if (procWinner[lnode] == myPID &&
            lnode < vertex2AggID.size() &&
            lnode < fineVectorData.size() && // do not access off-processor entries
            vertex2AggID[lnode] < coarseVectorData.size() &&
            Teuchos::ScalarTraits<SC>::isnaninf(fineVectorData[lnode]) == false) {
          coarseVectorData[vertex2AggID[lnode]] += fineVectorData[lnode];
        }
      }
      for (LO agg = 0; agg < numAggs; agg++) {
        coarseVectorData[agg] /= aggSizes[agg];
      }
    }

    Set<RCP<MultiVector> >(coarseLevel, vectorName, coarseVector);
  }

} // namespace MueLu

#endif // MUELU_MULTIVECTORAGGREGATIONTRANSFER_FACTORY_DEF_HPP
