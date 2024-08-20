// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_EMINPFACTORY_DEF_HPP
#define MUELU_EMINPFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_StridedMapFactory.hpp>

#include "MueLu_EminPFactory_decl.hpp"

#include "MueLu_CGSolver.hpp"
#include "MueLu_Constraint.hpp"
#include "MueLu_GMRESSolver.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_SolverBase.hpp"
#include "MueLu_SteepestDescentSolver.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("emin: num iterations");
  SET_VALID_ENTRY("emin: num reuse iterations");
  SET_VALID_ENTRY("emin: iterative method");
  {
    validParamList->getEntry("emin: iterative method").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("cg", "sd", "gmres"))));
  }
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory for the matrix A used during internal iterations");
  validParamList->set<RCP<const FactoryBase> >("P", Teuchos::null, "Generating factory for the initial guess");
  validParamList->set<RCP<const FactoryBase> >("Constraint", Teuchos::null, "Generating factory for constraints");

  validParamList->set<RCP<Matrix> >("P0", Teuchos::null, "Initial guess at P");
  validParamList->set<bool>("Keep P0", false, "Keep an initial P0 (for reuse)");

  validParamList->set<RCP<Constraint> >("Constraint0", Teuchos::null, "Initial Constraint");
  validParamList->set<bool>("Keep Constraint0", false, "Keep an initial Constraint (for reuse)");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(fineLevel, "A");

  static bool isAvailableP0          = false;
  static bool isAvailableConstraint0 = false;

  // Here is a tricky little piece of code
  // We don't want to request (aka call Input) when we reuse and P0 is available
  // However, we cannot run something simple like this:
  //   if (!coarseLevel.IsAvailable("P0", this))
  //     Input(coarseLevel, "P");
  // The reason is that it works fine during the request stage, but fails
  // in the release stage as we _construct_ P0 during Build process. Therefore,
  // we need to understand whether we are in Request or Release mode
  // NOTE: This is a very unique situation, please try not to propagate the
  // mode check any further

  if (coarseLevel.GetRequestMode() == Level::REQUEST) {
    isAvailableP0          = coarseLevel.IsAvailable("P0", this);
    isAvailableConstraint0 = coarseLevel.IsAvailable("Constraint0", this);
  }

  if (isAvailableP0 == false)
    Input(coarseLevel, "P");

  if (isAvailableConstraint0 == false)
    Input(coarseLevel, "Constraint");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Prolongator minimization", coarseLevel);

  const ParameterList& pL = GetParameterList();

  // Get the matrix
  RCP<Matrix> A = Get<RCP<Matrix> >(fineLevel, "A");

  if (restrictionMode_) {
    SubFactoryMonitor m2(*this, "Transpose A", coarseLevel);

    A = Utilities::Transpose(*A, true);
  }

  // Get/make initial guess
  RCP<Matrix> P0;
  int numIts;
  if (coarseLevel.IsAvailable("P0", this)) {
    // Reuse data
    P0     = coarseLevel.Get<RCP<Matrix> >("P0", this);
    numIts = pL.get<int>("emin: num reuse iterations");
    GetOStream(Runtime0) << "Reusing P0" << std::endl;

  } else {
    // Construct data
    P0     = Get<RCP<Matrix> >(coarseLevel, "P");
    numIts = pL.get<int>("emin: num iterations");
  }
  // NOTE: the main assumption here that P0 satisfies both constraints:
  //   - nonzero pattern
  //   - nullspace preservation

  // Get/make constraint operator
  RCP<Constraint> X;
  if (coarseLevel.IsAvailable("Constraint0", this)) {
    // Reuse data
    X = coarseLevel.Get<RCP<Constraint> >("Constraint0", this);
    GetOStream(Runtime0) << "Reusing Constraint0" << std::endl;

  } else {
    // Construct data
    X = Get<RCP<Constraint> >(coarseLevel, "Constraint");
  }
  GetOStream(Runtime0) << "Number of emin iterations = " << numIts << std::endl;

  std::string solverType = pL.get<std::string>("emin: iterative method");
  RCP<SolverBase> solver;
  if (solverType == "cg")
    solver = rcp(new CGSolver(numIts));
  else if (solverType == "sd")
    solver = rcp(new SteepestDescentSolver(numIts));
  else if (solverType == "gmres")
    solver = rcp(new GMRESSolver(numIts));

  RCP<Matrix> P;
  solver->Iterate(*A, *X, *P0, P);

  // NOTE: EXPERIMENTAL and FRAGILE
  if (!P->IsView("stridedMaps")) {
    if (A->IsView("stridedMaps") == true) {
      GetOStream(Runtime1) << "Using A to fillComplete P" << std::endl;

      // FIXME: X->GetPattern() actually returns a CrsGraph.
      // CrsGraph has no knowledge of Xpetra's sup/Matrix views. As such,
      // it has no idea about strided maps. We create one, which is
      // most likely incorrect for many use cases.
      std::vector<size_t> stridingInfo(1, 1);
      RCP<const StridedMap> dMap = StridedMapFactory::Build(X->GetPattern()->getDomainMap(), stridingInfo);

      P->CreateView("stridedMaps", A->getRowMap("stridedMaps"), dMap);

    } else {
      P->CreateView("stridedMaps", P->getRangeMap(), P->getDomainMap());
    }
  }

  // Level Set
  if (!restrictionMode_) {
    // The factory is in prolongation mode
    Set(coarseLevel, "P", P);

    if (pL.get<bool>("Keep P0")) {
      // NOTE: we must do Keep _before_ set as the Needs class only sets if
      //  a) data has been requested (which is not the case here), or
      //  b) data has some keep flag
      coarseLevel.Keep("P0", this);
      Set(coarseLevel, "P0", P);
    }
    if (pL.get<bool>("Keep Constraint0")) {
      // NOTE: we must do Keep _before_ set as the Needs class only sets if
      //  a) data has been requested (which is not the case here), or
      //  b) data has some keep flag
      coarseLevel.Keep("Constraint0", this);
      Set(coarseLevel, "Constraint0", X);
    }

    if (IsPrint(Statistics2)) {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo", true);
      GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*P, "P", params);
    }

  } else {
    // The factory is in restriction mode
    RCP<Matrix> R;
    {
      SubFactoryMonitor m2(*this, "Transpose P", coarseLevel);

      R = Utilities::Transpose(*P, true);
    }

    Set(coarseLevel, "R", R);

    if (IsPrint(Statistics2)) {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo", true);
      GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*R, "R", params);
    }
  }
}

}  // namespace MueLu

#endif  // MUELU_EMINPFACTORY_DEF_HPP
