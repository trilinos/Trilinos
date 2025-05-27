// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_CONSTRAINTFACTORY_DEF_HPP
#define MUELU_CONSTRAINTFACTORY_DEF_HPP

#include "MueLu_ConstraintFactory_decl.hpp"

#include "MueLu_Constraint.hpp"
#include "MueLu_DenseConstraint.hpp"
#include "MueLu_SparseConstraint.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> ConstraintFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase>>("FineNullspace", Teuchos::null, "Generating factory for the nullspace");
  validParamList->set<RCP<const FactoryBase>>("CoarseNullspace", Teuchos::null, "Generating factory for the nullspace");
  validParamList->set<RCP<const FactoryBase>>("Ppattern", Teuchos::null, "Generating factory for the nonzero pattern");
  validParamList->set<RCP<const FactoryBase>>("P_nodal", Teuchos::null, "Generating factory for nodal P");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ConstraintFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(fineLevel, "Nullspace", "FineNullspace");
  Input(coarseLevel, "Nullspace", "CoarseNullspace");
  Input(coarseLevel, "Ppattern");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ConstraintFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Constraint", coarseLevel);

  auto Ppattern = Get<RCP<const CrsGraph>>(coarseLevel, "Ppattern");
  RCP<Constraint> constraint;
  if (fineLevel.IsType<RCP<MultiVector>>("Nullspace", GetFactory("FineNullspace").get()) &&
      coarseLevel.IsType<RCP<MultiVector>>("Nullspace", GetFactory("CoarseNullspace").get())) {
    RCP<MultiVector> fineNullspace   = Get<RCP<MultiVector>>(fineLevel, "Nullspace", "FineNullspace");
    RCP<MultiVector> coarseNullspace = Get<RCP<MultiVector>>(coarseLevel, "Nullspace", "CoarseNullspace");
    constraint                       = rcp(new DenseConstraint(fineNullspace, coarseNullspace, Ppattern));
  } else if (fineLevel.IsType<RCP<Matrix>>("Nullspace", GetFactory("FineNullspace").get()) &&
             coarseLevel.IsType<RCP<Matrix>>("Nullspace", GetFactory("CoarseNullspace").get())) {
    RCP<Matrix> fineNullspace   = Get<RCP<Matrix>>(fineLevel, "Nullspace", "FineNullspace");
    RCP<Matrix> coarseNullspace = Get<RCP<Matrix>>(coarseLevel, "Nullspace", "CoarseNullspace");
    RCP<Matrix> P_nodal         = Get<RCP<Matrix>>(coarseLevel, "P_nodal", "P_nodal");
    constraint                  = rcp(new SparseConstraint(P_nodal, fineNullspace, coarseNullspace, Ppattern));
  } else {
    TEUCHOS_ASSERT(false);
  }

  Set(coarseLevel, "Constraint", constraint);
}

}  // namespace MueLu

#endif  // MUELU_CONSTRAINTFACTORY_DEF_HPP
