// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TRANSPFACTORY_DEF_HPP
#define MUELU_TRANSPFACTORY_DEF_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Time.hpp>

#include <Xpetra_Matrix.hpp>

#include "MueLu_TransPFactory_decl.hpp"

#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->set<RCP<const FactoryBase> >("P", Teuchos::null, "Generating factory of the matrix P");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& /* fineLevel */, Level& coarseLevel) const {
  Input(coarseLevel, "P");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& /* fineLevel */, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Transpose P", coarseLevel);
  std::string label = "MueLu::TransP-" + Teuchos::toString(coarseLevel.GetLevelID());

  RCP<Matrix> P = Get<RCP<Matrix> >(coarseLevel, "P");
  // If we failed to create a valid P (e.g., # of global aggregates is zero), then we just bail here
  //  This level will ultimately be removed in MueLu_Hierarchy_defs.h via a resize()
  if (P == Teuchos::null) return;

  const Teuchos::ParameterList& pL = GetParameterList();

  // Reuse pattern if available (multiple solve)
  RCP<ParameterList> Tparams;
  if (pL.isSublist("matrixmatrix: kernel params"))
    Tparams = rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));
  else
    Tparams = rcp(new ParameterList);

  // By default, we don't need global constants for transpose
  Tparams->set("compute global constants: temporaries", Tparams->get("compute global constants: temporaries", false));
  Tparams->set("compute global constants", Tparams->get("compute global constants", false));

  RCP<Matrix> R = Utilities::Transpose(*P, true, label, Tparams);

  if (IsPrint(Statistics2)) {
    RCP<ParameterList> params = rcp(new ParameterList());
    params->set("printLoadBalancingInfo", true);
    params->set("printCommInfo", true);
    GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*R, "R", params);
  }

  Set(coarseLevel, "R", R);

  ///////////////////////// EXPERIMENTAL
  if (P->IsView("stridedMaps"))
    R->CreateView("stridedMaps", P, true);
  ///////////////////////// EXPERIMENTAL
}

}  // namespace MueLu

#endif  // MUELU_TRANSPFACTORY_DEF_HPP
