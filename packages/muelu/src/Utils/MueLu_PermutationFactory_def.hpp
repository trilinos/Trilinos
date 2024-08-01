// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * MueLu_PermutationFactory_def.hpp
 *
 *  Created on: Nov 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_PERMUTATIONFACTORY_DEF_HPP_
#define MUELU_PERMUTATIONFACTORY_DEF_HPP_

#include <vector>
#include <queue>

#include "MueLu_PermutationFactory_decl.hpp"

#include <Xpetra_Map.hpp>
#include <Xpetra_StridedMap.hpp>  // for nDofsPerNode...
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_AlgebraicPermutationStrategy.hpp"
#include "MueLu_LocalPermutationStrategy.hpp"

#undef DEBUG_OUTPUT

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PermutationFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~PermutationFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A to be permuted.");

  validParamList->set<std::string>("PermutationRowMapName", "", "Name of input row map for which rows the permutation shall be done. (default='')");
  validParamList->set<RCP<const FactoryBase> >("PermutationRowMapFactory", Teuchos::null, "Generating factory of the input row map for the permutation.");

  validParamList->set<std::string>("PermutationStrategy", "Algebraic", "Permutation strategy (default = 'Algebraic', 'Local'");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  Input(currentLevel, "A");

  const ParameterList &pL                    = GetParameterList();
  std::string mapName                        = pL.get<std::string>("PermutationRowMapName");
  Teuchos::RCP<const FactoryBase> mapFactory = GetFactory("PermutationRowMapFactory");

  if (mapName.length() > 0) {
    currentLevel.DeclareInput(mapName, mapFactory.get(), this);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
  FactoryMonitor m(*this, "Permutation Factory ", currentLevel);

  Teuchos::RCP<Matrix> A = Get<Teuchos::RCP<Matrix> >(currentLevel, "A");

  const ParameterList &pL                    = GetParameterList();
  std::string mapName                        = pL.get<std::string>("PermutationRowMapName");
  Teuchos::RCP<const FactoryBase> mapFactory = GetFactory("PermutationRowMapFactory");

  Teuchos::RCP<const Map> permRowMap = Teuchos::null;
  if (mapName.length() > 0) {
    permRowMap = currentLevel.Get<RCP<const Map> >(mapName, mapFactory.get());
  } else {
    permRowMap = A->getRowMap();  // use full row map of A
  }

  std::string strStrategy = pL.get<std::string>("PermutationStrategy");
  if (strStrategy == "Algebraic") {
    Teuchos::RCP<AlgebraicPermutationStrategy> permStrat = Teuchos::rcp(new AlgebraicPermutationStrategy());
    permStrat->BuildPermutation(A, permRowMap, currentLevel, this);
  } else if (strStrategy == "Local") {
    Teuchos::RCP<LocalPermutationStrategy> permStrat = Teuchos::rcp(new LocalPermutationStrategy());
    permStrat->BuildPermutation(A, permRowMap, currentLevel, this);
  } else
    TEUCHOS_TEST_FOR_EXCEPTION(true,
                               std::logic_error,
                               "`PermutationStrategy' has incorrect value (" << strStrategy << ") in input to PermutationFactory."
                                                                             << "Check the documentation for a list of valid choices");

  GetOStream(Runtime0) << "Using " << strStrategy << " permutation strategy." << std::endl;
}

}  // namespace MueLu

#endif /* MUELU_PERMUTATIONFACTORY_DEF_HPP_ */
