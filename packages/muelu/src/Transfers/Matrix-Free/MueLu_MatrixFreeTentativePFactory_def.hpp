// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MATRIXFREETENTATIVEPFACTORY_DEF_HPP
#define MUELU_MATRIXFREETENTATIVEPFACTORY_DEF_HPP

#include "Kokkos_UnorderedMap.hpp"

#include "MueLu_MatrixFreeTentativePFactory_decl.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_MatrixFreeTentativeP.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> MatrixFreeTentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase>>("Aggregates", Teuchos::null, "Generating factory of the aggregates");
  validParamList->set<RCP<const FactoryBase>>("Nullspace", Teuchos::null, "Generating factory of the nullspace");
  validParamList->set<RCP<const FactoryBase>>("Scaled Nullspace", Teuchos::null, "Generating factory of the scaled nullspace");
  validParamList->set<RCP<const FactoryBase>>("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
  validParamList->set<RCP<const FactoryBase>>("CoarseMap", Teuchos::null, "Generating factory of the coarse map");
  validParamList->set<RCP<const FactoryBase>>("Coordinates", Teuchos::null, "Generating factory of the coordinates");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixFreeTentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& /* coarseLevel */) const {
  const ParameterList& pL = GetParameterList();
  // NOTE: This guy can only either be 'Nullspace' or 'Scaled Nullspace' or else the validator above will cause issues
  std::string nspName = "Nullspace";
  if (pL.isParameter("Nullspace name")) nspName = pL.get<std::string>("Nullspace name");

  Input(fineLevel, "Aggregates");
  Input(fineLevel, nspName);
  Input(fineLevel, "UnAmalgamationInfo");
  Input(fineLevel, "CoarseMap");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixFreeTentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixFreeTentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  const ParameterList& pL = GetParameterList();
  std::string nspName     = "Nullspace";
  if (pL.isParameter("Nullspace name")) nspName = pL.get<std::string>("Nullspace name");

  auto aggregates                 = Get<RCP<Aggregates>>(fineLevel, "Aggregates");
  auto amalgInfo                  = Get<RCP<AmalgamationInfo>>(fineLevel, "UnAmalgamationInfo");
  auto fineNullspace              = Get<RCP<MultiVector>>(fineLevel, nspName);
  auto coarseMap                  = Get<RCP<const Map>>(fineLevel, "CoarseMap");
  Teuchos::RCP<const Map> fineMap = fineNullspace->getMap();

  // Matrix-free should never run with aggregates that cross processors
  if (aggregates->AggregatesCrossProcessors())
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MatrixFreeTentativePFactory does not support aggregates that cross processors!");

  size_t NSDim                     = fineNullspace->getNumVectors();
  RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(coarseMap, NSDim);

  Teuchos::RCP<Operator> P = Teuchos::rcp(new MatrixFreeTentativeP<Scalar, LocalOrdinal, GlobalOrdinal, Node>(coarseMap, fineMap, aggregates));
  P->apply(*fineNullspace, *coarseNullspace, Teuchos::TRANS, 1.0, 0.0);  // coarse = alpha*R*fine + beta*coarse

  Set(coarseLevel, "Nullspace", coarseNullspace);
  Set(coarseLevel, "P", P);
}

}  // namespace MueLu

#define MUELU_MATRIXFREETENTATIVEPFACTORY_SHORT
#endif  // MUELU_MATRIXFREETENTATIVEPFACTORY_DEF_HPP
