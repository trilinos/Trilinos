// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP
#define MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP

#include "MueLu_MultiVectorTransferFactory_decl.hpp"
#include "Tpetra_Access.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<std::string>("Vector name", "undefined", "Name of the vector that will be transferred on the coarse grid (level key)");  // TODO: how to set a validator without default value?
  validParamList->set<std::string>("Transfer name", "R", "Name of the operator that will be used to transfer data");
  validParamList->set<bool>("Normalize", false, "If a row sum normalization should be applied to preserve the mean value of the vector.");
  validParamList->set<RCP<const FactoryBase>>("Vector factory", Teuchos::null, "Factory of the vector");
  validParamList->set<RCP<const FactoryBase>>("Transfer factory", Teuchos::null, "Factory of the transfer operator");
  validParamList->set<RCP<const FactoryBase>>("CoarseMap", Teuchos::null, "Generating factory of the coarse map");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  const ParameterList &pL  = GetParameterList();
  std::string vectorName   = pL.get<std::string>("Vector name");
  std::string transferName = pL.get<std::string>("Transfer name");

  fineLevel.DeclareInput(vectorName, GetFactory("Vector factory").get(), this);
  auto transferFact             = GetFactory("Transfer factory");
  const bool isUncoupledAggFact = !Teuchos::rcp_dynamic_cast<const UncoupledAggregationFactory>(transferFact).is_null();
  if (isUncoupledAggFact) {
    fineLevel.DeclareInput(transferName, transferFact.get(), this);
    Input(fineLevel, "CoarseMap");
  } else
    coarseLevel.DeclareInput(transferName, transferFact.get(), this);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  const ParameterList &pL        = GetParameterList();
  const std::string vectorName   = pL.get<std::string>("Vector name");
  const std::string transferName = pL.get<std::string>("Transfer name");
  const bool normalize           = pL.get<bool>("Normalize");

  auto transferFact             = GetFactory("Transfer factory");
  const bool isUncoupledAggFact = !Teuchos::rcp_dynamic_cast<const UncoupledAggregationFactory>(transferFact).is_null();

  GetOStream(Runtime0) << "Transferring multivector \"" << vectorName << "\" using operator \"" << transferName << "\"" << std::endl;
  if (vectorName == "Coordinates")
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Use CoordinatesTransferFactory to transfer coordinates instead of MultiVectorTransferFactory.");

  RCP<MultiVector> fineVector = fineLevel.Get<RCP<MultiVector>>(vectorName, GetFactory("Vector factory").get());
  RCP<MultiVector> coarseVector;

  if (!isUncoupledAggFact) {
    RCP<Matrix> transferOp = coarseLevel.Get<RCP<Matrix>>(transferName, GetFactory("Transfer factory").get());
    Teuchos::ETransp transp;

    if (transferOp->getGlobalNumRows() <= transferOp->getGlobalNumCols()) {
      // R mode
      coarseVector = MultiVectorFactory::Build(transferOp->getRangeMap(), fineVector->getNumVectors());
      transp       = Teuchos::NO_TRANS;
    } else {
      // P mode
      coarseVector = MultiVectorFactory::Build(transferOp->getDomainMap(), fineVector->getNumVectors());
      transp       = Teuchos::TRANS;
    }

    transferOp->apply(*fineVector, *coarseVector, transp);

    if (normalize) {
      // Do constant row sum normalization
      RCP<MultiVector> onesVector = MultiVectorFactory::Build(fineVector->getMap(), 1);
      onesVector->putScalar(Teuchos::ScalarTraits<Scalar>::one());
      RCP<MultiVector> rowSumVector = MultiVectorFactory::Build(coarseVector->getMap(), 1);
      transferOp->apply(*onesVector, *rowSumVector, transp);

      RCP<Vector> rowSumReciprocalVector = VectorFactory::Build(coarseVector->getMap(), 1);
      rowSumReciprocalVector->reciprocal(*rowSumVector);

      RCP<MultiVector> coarseVectorNormalized = MultiVectorFactory::Build(coarseVector->getMap(), fineVector->getNumVectors());
      coarseVectorNormalized->elementWiseMultiply(1.0, *rowSumReciprocalVector, *coarseVector, 0.0);

      Set<RCP<MultiVector>>(coarseLevel, vectorName, coarseVectorNormalized);
    } else {
      Set<RCP<MultiVector>>(coarseLevel, vectorName, coarseVector);
    }
  } else {
    using execution_space = typename Node::execution_space;
#if KOKKOS_VERSION >= 40799
    using ATS = KokkosKernels::ArithTraits<Scalar>;
#else
    using ATS = Kokkos::ArithTraits<Scalar>;
#endif
    using impl_scalar_type = typename ATS::val_type;

    auto aggregates = fineLevel.Get<RCP<Aggregates>>(transferName, GetFactory("Transfer factory").get());
    TEUCHOS_ASSERT(!aggregates->AggregatesCrossProcessors());
    RCP<const Map> coarseMap = Get<RCP<const Map>>(fineLevel, "CoarseMap");

    auto aggGraph = aggregates->GetGraph();
    auto numAggs  = aggGraph.numRows();

    RCP<const Map> coarseVectorMap;

    LO blkSize = 1;
    if (rcp_dynamic_cast<const StridedMap>(coarseMap) != Teuchos::null)
      blkSize = rcp_dynamic_cast<const StridedMap>(coarseMap)->getFixedBlockSize();

    if (blkSize == 1) {
      // Scalar system
      // No amalgamation required, we can use the coarseMap
      coarseVectorMap = coarseMap;
    } else {
      // Vector system
      AmalgamationFactory<SC, LO, GO, NO>::AmalgamateMap(rcp_dynamic_cast<const StridedMap>(coarseMap), coarseVectorMap);
    }

    coarseVector = MultiVectorFactory::Build(coarseVectorMap, fineVector->getNumVectors());

    auto lcl_fineVector   = fineVector->getLocalViewDevice(Tpetra::Access::ReadOnly);
    auto lcl_coarseVector = coarseVector->getLocalViewDevice(Tpetra::Access::OverwriteAll);

    Kokkos::parallel_for(
        "MueLu:MultiVectorTransferFactory",
        Kokkos::RangePolicy<LocalOrdinal, execution_space>(0, numAggs),
        KOKKOS_LAMBDA(const LO i) {
          auto aggregate = aggGraph.rowConst(i);
          for (size_t j = 0; j < lcl_coarseVector.extent(1); j++) {
            impl_scalar_type sum = 0.0;
            for (size_t colID = 0; colID < static_cast<size_t>(aggregate.length); colID++)
              sum += lcl_fineVector(aggregate(colID), j);
            if (normalize)
              lcl_coarseVector(i, j) = sum / aggregate.length;
            else
              lcl_coarseVector(i, j) = sum;
          }
        });
    Set<RCP<MultiVector>>(coarseLevel, vectorName, coarseVector);
  }
}  // Build

}  // namespace MueLu

#endif  // MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP
