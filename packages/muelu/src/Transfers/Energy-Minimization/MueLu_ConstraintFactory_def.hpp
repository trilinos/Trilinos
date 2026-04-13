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
#include "MueLu_Monitor.hpp"
#include "MueLu_MasterList.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> ConstraintFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("emin: constraint type");
  validParamList->getEntry("emin: constraint type").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("nullspace", "maxwell"))));

  SET_VALID_ENTRY("emin: least squares solver type");
  validParamList->getEntry("emin: least squares solver type").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("Belos", "direct"))));
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase>>("FineNullspace", Teuchos::null, "Generating factory for the nullspace");
  validParamList->set<RCP<const FactoryBase>>("CoarseNullspace", Teuchos::null, "Generating factory for the nullspace");
  validParamList->set<RCP<const FactoryBase>>("Ppattern", Teuchos::null, "Generating factory for the nonzero pattern");

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
  const ParameterList& pL = GetParameterList();

  std::string solverType = pL.get<std::string>("emin: least squares solver type");

  if (pL.get<std::string>("emin: constraint type") == "nullspace") {
    RCP<MultiVector> fineNullspace   = Get<RCP<MultiVector>>(fineLevel, "Nullspace", "FineNullspace");
    RCP<MultiVector> coarseNullspace = Get<RCP<MultiVector>>(coarseLevel, "Nullspace", "CoarseNullspace");
    {
      SubFactoryMonitor m2(*this, "Dense Constraint", coarseLevel);
      constraint = rcp(new DenseConstraint(fineNullspace, coarseNullspace, Ppattern, solverType));
    }
  } else {
    TEUCHOS_ASSERT(false);
  }

  Set(coarseLevel, "Constraint", constraint);
}

}  // namespace MueLu

#endif  // MUELU_CONSTRAINTFACTORY_DEF_HPP
