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
  validParamList->set<RCP<const FactoryBase>>("FineD0", Teuchos::null, "Generating factory for the fine discrete gradient");
  validParamList->set<RCP<const FactoryBase>>("CoarseD0", Teuchos::null, "Generating factory for the coarse discrete gradient");
  validParamList->set<RCP<const FactoryBase>>("Ppattern", Teuchos::null, "Generating factory for the nonzero pattern");
  validParamList->set<RCP<const FactoryBase>>("PnodalEmin", Teuchos::null, "Generating factory for nodal P");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ConstraintFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  const ParameterList& pL = GetParameterList();
  if (pL.get<std::string>("emin: constraint type") == "nullspace") {
    Input(fineLevel, "Nullspace", "FineNullspace");
    Input(coarseLevel, "Nullspace", "CoarseNullspace");
  } else {
    Input(fineLevel, "D0", "FineD0");
    Input(coarseLevel, "D0", "CoarseD0");
    Input(coarseLevel, "PnodalEmin", "PnodalEmin");
  }
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
    RCP<Matrix> fineD0   = Get<RCP<Matrix>>(fineLevel, "D0", "FineD0");
    RCP<Matrix> coarseD0 = Get<RCP<Matrix>>(coarseLevel, "D0", "CoarseD0");
    RCP<Matrix> Pnodal   = Get<RCP<Matrix>>(coarseLevel, "PnodalEmin", "PnodalEmin");
    RCP<SparseConstraint> sparse_constraint;
    {
      SubFactoryMonitor m2(*this, "Sparse Constraint", coarseLevel);
      sparse_constraint = rcp(new SparseConstraint(Pnodal, fineD0, coarseD0, Ppattern, solverType));
    }

    // Construct an initial guess.
    // This is different from the nullspace constraint where we use the tentative prolongator as initial guess.
    // Instead, we leverage the least-squares solve to construct a prolongator that satisfies the constraints.

    auto X             = sparse_constraint->GetConstraintMatrix();
    auto vecPProjected = MultiVectorFactory::Build(sparse_constraint->getDomainMap(), 1);
    auto temp          = MultiVectorFactory::Build(X->getRangeMap(), 1);
    auto rhs           = MultiVectorFactory::Build(X->getRangeMap(), 1);

    // We want
    //  P * coarseD0 = fineD0 * Pnodal
    // which is, after vectorization
    //  X * vec(P) = vec(fineD0 * Pnodal)
    // A solution of this is given by
    //  vec(P) = X^T * (X * X^T)^dagger * vec(fineD0 * Pnodal)

    // rhs = fineD0 * Pnodal
    {
      RCP<Matrix> fineD0_Pnodal;
      fineD0_Pnodal = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*fineD0, false, *Pnodal, false, fineD0_Pnodal, GetOStream(Runtime0), true, true);
      sparse_constraint->AssignMatrixEntriesToConstraintVector(*fineD0_Pnodal, *rhs);
    }

    // vecPProjected = X^T (X * X^T)^dagger * rhs
    sparse_constraint->LeastSquaresSolve(*rhs, *temp);
    X->apply(*temp, *vecPProjected, Teuchos::TRANS);

    auto P0 = sparse_constraint->GetMatrixWithEntriesFromVector(*vecPProjected);

    Set(coarseLevel, "P", P0);
    constraint = sparse_constraint;
  }

  Set(coarseLevel, "Constraint", constraint);
}

}  // namespace MueLu

#endif  // MUELU_CONSTRAINTFACTORY_DEF_HPP
