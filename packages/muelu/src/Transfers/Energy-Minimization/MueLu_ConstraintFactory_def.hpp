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
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> ConstraintFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("Nullspace", Teuchos::null, "Generating factory for the nullspace");
  validParamList->set<RCP<const FactoryBase> >("Ppattern", Teuchos::null, "Generating factory for the nonzero pattern");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ConstraintFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(coarseLevel, "Nullspace");
  Input(coarseLevel, "Ppattern");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ConstraintFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Constraint", coarseLevel);

  RCP<MultiVector> nullspace = Get<RCP<MultiVector> >(coarseLevel, "Nullspace");

  RCP<Constraint> constraint(new Constraint);
  constraint->Setup(*nullspace, Get<RCP<const CrsGraph> >(coarseLevel, "Ppattern"));

  Set(coarseLevel, "Constraint", constraint);
}

}  // namespace MueLu

#endif  // MUELU_CONSTRAINTFACTORY_DEF_HPP
