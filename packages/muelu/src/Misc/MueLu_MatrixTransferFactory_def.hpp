// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MATRIXTRANSFER_FACTORY_DEF_HPP
#define MUELU_MATRIXTRANSFER_FACTORY_DEF_HPP

#include "MueLu_MatrixTransferFactory_decl.hpp"

#include <Xpetra_Matrix.hpp>

#include "MueLu_MasterList.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> MatrixTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("transpose: use implicit");
  SET_VALID_ENTRY("rap: triple product");
  SET_VALID_ENTRY("rap: fix zero diagonals");
  SET_VALID_ENTRY("rap: fix zero diagonals threshold");
  SET_VALID_ENTRY("rap: fix zero diagonals replacement");
  SET_VALID_ENTRY("rap: relative diagonal floor");
#undef SET_VALID_ENTRY

  validParamList->set<std::string>("Matrix name", "A", "Name of the matrix that will be transferred on the coarse grid (level key)");
  validParamList->set<std::string>("Restrictor name", "R", "Name of the operator that will be used to transfer data");
  validParamList->set<std::string>("Prolongator name", "P", "Name of the operator that will be used to transfer data");

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");
  validParamList->set<RCP<const FactoryBase>>("P", Teuchos::null, "Prolongator factory");
  validParamList->set<RCP<const FactoryBase>>("R", Teuchos::null, "Restrictor factory");

  validParamList->set<bool>("CheckMainDiagonal", false, "Check main diagonal for zeros");
  validParamList->set<bool>("RepairMainDiagonal", false, "Repair zeros on main diagonal");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  const ParameterList &pL = GetParameterList();
  std::string matrixName  = pL.get<std::string>("Matrix name");
  fineLevel.DeclareInput(matrixName, GetFactory("A").get(), this);

  std::string prolongatorName = pL.get<std::string>("Prolongator name");
  coarseLevel.DeclareInput(prolongatorName, GetFactory("P").get(), this);

  if (!pL.get<bool>("transpose: use implicit")) {
    std::string restrictorName = pL.get<std::string>("Restrictor name");
    coarseLevel.DeclareInput(restrictorName, GetFactory("R").get(), this);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  const ParameterList &pL           = GetParameterList();
  const std::string matrixName      = pL.get<std::string>("Matrix name");
  const std::string prolongatorName = pL.get<std::string>("Prolongator name");

  RCP<Matrix> A = fineLevel.Get<RCP<Matrix>>(matrixName, GetFactory("A").get());
  RCP<Matrix> P = coarseLevel.Get<RCP<Matrix>>(prolongatorName, GetFactory("P").get());
  RCP<Matrix> R;
  if (!pL.get<bool>("transpose: use implicit")) {
    const std::string restrictorName = pL.get<std::string>("Restrictor name");
    R                                = coarseLevel.Get<RCP<Matrix>>(restrictorName, GetFactory("R").get());
  }

  RCP<Matrix> Ac;

  Teuchos::RCP<Teuchos::ParameterList> APparams, RAPparams;
  Utilities::TripleMatrixProduct(R, A, P, Ac, pL, *this, APparams, RAPparams, &coarseLevel);

  Set<RCP<Matrix>>(coarseLevel, matrixName, Ac);

}  // Build

}  // namespace MueLu

#endif
