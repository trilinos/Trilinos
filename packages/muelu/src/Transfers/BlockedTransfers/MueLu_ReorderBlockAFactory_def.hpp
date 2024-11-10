// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REORDERBLOCKAFACTORY_DEF_HPP_
#define MUELU_REORDERBLOCKAFACTORY_DEF_HPP_

#include "MueLu_ReorderBlockAFactory_decl.hpp"

#include <Xpetra_BlockReorderManager.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_ReorderedBlockedCrsMatrix.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixUtils.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> ReorderBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", MueLu::NoFactory::getRCP(), "Generating factory for A.");

  validParamList->set<std::string>("Reorder Type", "", "String describing the reordering of blocks");

  // TODO not very elegant.
  validParamList->set<RCP<const FactoryBase> >("Map1", Teuchos::null, "Generating factory of the fine level map associated with the (1,1) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase> >("Map2", Teuchos::null, "Generating factory of the fine level map associated with the (2,2) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase> >("Map3", Teuchos::null, "Generating factory of the fine level map associated with the (3,3) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase> >("Map4", Teuchos::null, "Generating factory of the fine level map associated with the (4,4) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase> >("Map5", Teuchos::null, "Generating factory of the fine level map associated with the (5,5) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase> >("Map6", Teuchos::null, "Generating factory of the fine level map associated with the (6,6) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase> >("Map7", Teuchos::null, "Generating factory of the fine level map associated with the (7,7) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase> >("Map8", Teuchos::null, "Generating factory of the fine level map associated with the (8,8) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase> >("Map9", Teuchos::null, "Generating factory of the fine level map associated with the (9,9) block in your n x n block matrix.");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReorderBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReorderBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "ReorderBlockA factory", currentLevel);

  const ParameterList& pL = GetParameterList();
  std::string reorderStr  = pL.get<std::string>("Reorder Type");

  RCP<Matrix> Ain = Get<RCP<Matrix> >(currentLevel, "A");

  RCP<BlockedCrsMatrix> A = rcp_dynamic_cast<BlockedCrsMatrix>(Ain);

  // special case: we get a single block CrsMatrix object on the finest level and
  // split it into a nxn blocked operator
  if (A == Teuchos::null && currentLevel.GetLevelID() == 0) {
    GetOStream(Warnings0) << "Split input matrix (Warning: this is a rather expensive operation)" << std::endl;

    std::vector<Teuchos::RCP<const Map> > xmaps;

    for (int it = 1; it < 10; it++) {
      std::stringstream ss;
      ss << "Map" << it;
      if (currentLevel.IsAvailable(ss.str(), NoFactory::get())) {
        RCP<const Map> submap = currentLevel.Get<RCP<const Map> >(ss.str(), NoFactory::get());
        GetOStream(Runtime1) << "Use user-given submap #" << it << ": length dimension=" << submap->getGlobalNumElements() << std::endl;
        xmaps.push_back(submap);
      }
    }

    bool bThyraMode                       = false;  // no support for Thyra mode (yet)
    RCP<const MapExtractor> map_extractor = Xpetra::MapExtractorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Ain->getRowMap(), xmaps, bThyraMode);

    // split null space vectors
    // TODO: if he matrix blocks have different striding, this could be quite complicated
    // RCP<MultiVector> nullspace1 = map_extractor->ExtractVector(nullspace,0);
    // RCP<MultiVector> nullspace2 = map_extractor->ExtractVector(nullspace,1);

    Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > bOp =
        Xpetra::MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SplitMatrix(*Ain, map_extractor, map_extractor, Teuchos::null, bThyraMode);

    TEUCHOS_TEST_FOR_EXCEPTION(Ain->getGlobalNumRows() != bOp->getGlobalNumRows(), Exceptions::RuntimeError, "Split operator not consistent with input operator (different number of rows).");
    TEUCHOS_TEST_FOR_EXCEPTION(Ain->getLocalNumRows() != bOp->getLocalNumRows(), Exceptions::RuntimeError, "Split operator not consistent with input operator (different number of node rows).");
    TEUCHOS_TEST_FOR_EXCEPTION(Ain->getLocalNumEntries() != bOp->getLocalNumEntries(), Exceptions::RuntimeError, "Split operator not consistent with input operator (different number of local entries).");
    TEUCHOS_TEST_FOR_EXCEPTION(Ain->getGlobalNumEntries() != bOp->getGlobalNumEntries(), Exceptions::RuntimeError, "Split operator not consistent with input operator (different number of global entries).");

    A = bOp;
  }

  // we have a blocked operator as input
  TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(), Exceptions::BadCast, "Input matrix A is not a BlockedCrsMatrix.");
  GetOStream(Statistics1) << "Got a " << A->Rows() << "x" << A->Cols() << " blocked operator as input" << std::endl;

  // if we have a blocked operator and a reordering string, create a nested blocked operator, if not skip the process
  if (reorderStr.empty()) {
    GetOStream(Statistics1) << "No reordering information provided. Skipping reordering of A." << std::endl;
  } else {
    Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString(reorderStr);
    GetOStream(Debug) << "Reordering A using " << brm->toString() << std::endl;

    Teuchos::RCP<const ReorderedBlockedCrsMatrix> brop =
        Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(
            Xpetra::buildReorderedBlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(brm, A));

    TEUCHOS_TEST_FOR_EXCEPTION(brop.is_null(), Exceptions::RuntimeError,
                               "Block reordering of " << A->Rows() << "x" << A->Cols()
                                                      << " blocked operator failed.");

    GetOStream(Statistics1) << "Reordering A using " << brm->toString() << " block gives a " << brop->Rows() << "x"
                            << brop->Cols() << " blocked operator" << std::endl;
    GetOStream(Debug) << "Reordered operator has " << brop->getRangeMap()->getGlobalNumElements() << " rows and "
                      << brop->getDomainMap()->getGlobalNumElements() << " columns" << std::endl;
    GetOStream(Debug) << "Reordered operator: Use of Thyra style gids = "
                      << brop->getRangeMapExtractor()->getThyraMode() << std::endl;

    // get rid of const (we expect non-const operators stored in Level)
    Teuchos::RCP<ReorderedBlockedCrsMatrix> bret =
        Teuchos::rcp_const_cast<ReorderedBlockedCrsMatrix>(brop);

    A = bret;
  }

  currentLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(A), this);
}

}  // namespace MueLu

#endif /* MUELU_REORDERBLOCKAFACTORY_DEF_HPP_ */
