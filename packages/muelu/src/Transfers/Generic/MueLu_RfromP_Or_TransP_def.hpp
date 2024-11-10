// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_RFROMP_OR_TRANSP_DEF_HPP
#define MUELU_RFROMP_OR_TRANSP_DEF_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Time.hpp>

#include <Xpetra_Matrix.hpp>

#include "MueLu_RfromP_Or_TransP_decl.hpp"

#include "MueLu_DisableMultipleCallCheck.hpp"

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_TogglePFactory.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RfromP_Or_TransP<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->set<RCP<const FactoryBase>>("P", Teuchos::null, "Generating factory of the matrix P");
  validParamList->set<RCP<const FactoryBase>>("RfromPfactory", Teuchos::null, "Generating factory of the matrix R");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RfromP_Or_TransP<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& /* fineLevel */, Level& coarseLevel) const {
  Input(coarseLevel, "RfromPfactory");

  // Using a PgPFactory in conjunction with a TogglePFactory is a bit problematic. Normally, PgPFactory is supposed to be
  // invoked twice (once in standard mode and a second time in RestrictionMode). Unfortunately, TogglePFactory
  // is not designed to produce R. A second issue is that TogglePFactory stores prolongators in an array,
  // and there are some challenges in determining which array entry (i.e., factory) is needed when producing R.
  // The way this is addressed is a bit clumsy. RfromP_Or_TransP invokes the prolongator factory to produce R.
  // To do this, it must first check that this is needed (as opposed to just transposing P or using an already computed
  // R that might be produced by SemiCoarsenPFactory). This check is needed in both DeclareInput() and in Build().
  // The DeclareInput() check verifies that TogglePFactory was requested and that one of the prolongator factories
  // within TogglePFactory is a PgPFactory.  RfromP_Or_TransP's DeclareInput then invokes DeclareDependencies, and
  // DeclareInput for the PgPFactory. The check within Build(), looks at "RfromPFactory" to see if it is an integer. This
  // integer is used to find the prolongator factory that is invoked in RestrictionMode to produce R. Otherwise,
  // "RfromPFactory" is used to get the pre-computed restriction matrix. If "RfromPFactory" is not present, then RfromP_Or_TransP
  // just transposes P to get R.

  RCP<const FactoryBase> PFact = coarseLevel.GetFactoryManager()->GetFactory("P");
  if (PFact == Teuchos::null) {
    PFact = GetFactory("P");
  }
  coarseLevel.DeclareInput("P", PFact.get(), this);
  RCP<const MueLu::TogglePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>> myToggleFact = Teuchos::rcp_const_cast<const MueLu::TogglePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(rcp_dynamic_cast<const MueLu::TogglePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(PFact));
  if (myToggleFact != Teuchos::null) {
    for (size_t ii = 0; ii < myToggleFact->NumProlongatorFactories(); ii++) {
      RCP<const MueLu::PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>> actualPFact = Teuchos::rcp_const_cast<const MueLu::PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(rcp_dynamic_cast<const MueLu::PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(myToggleFact->getProlongatorFactory(ii)));
      if (actualPFact != Teuchos::null) {
        RCP<PFactory> subFactory = Teuchos::rcp_const_cast<PFactory>(rcp_dynamic_cast<const PFactory>(myToggleFact->getProlongatorFactory(ii)));
        ;
        bool rmode = subFactory->isRestrictionModeSet();
        subFactory->setRestrictionMode(true);
        // Force request call for actualPFact
        coarseLevel.DeclareDependencies(actualPFact.get());
        coarseLevel.DeclareInput("R", actualPFact.get(), this);
        subFactory->setRestrictionMode(rmode);
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RfromP_Or_TransP<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Transpose P", coarseLevel);
  std::string label = "MueLu::TransP-" + Teuchos::toString(coarseLevel.GetLevelID());

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

  RCP<Matrix> R;
  RCP<const FactoryBase> PFact = coarseLevel.GetFactoryManager()->GetFactory("P");
  if (PFact == Teuchos::null) {
    PFact = GetFactory("P");
  }

  RCP<Matrix> P = coarseLevel.Get<RCP<Matrix>>("P", PFact.get());

  if (coarseLevel.IsAvailable("RfromPfactory", PFact.get())) {
    std::string strType = coarseLevel.GetTypeName("RfromPfactory", PFact.get());
    // Address case where a toggle factory is used in conjunction with a PgP factory. Here,
    // we need to invoke the PgP factory a 2nd time to produce restriction. In this
    // situation the PgP factory puts an int RfromPfactory on the level.
    //
    // See comments aboove in DeclareInput() method  for more detailsd
    if (strType == "int") {
      RCP<const MueLu::TogglePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>> myToggleFact = Teuchos::rcp_const_cast<const MueLu::TogglePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(rcp_dynamic_cast<const MueLu::TogglePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(PFact));
      ;
      if (myToggleFact != Teuchos::null) {
        MueLu::DisableMultipleCallCheck check(myToggleFact);
        RCP<PFactory> actualPFact = Teuchos::rcp_const_cast<PFactory>(rcp_dynamic_cast<const PFactory>(myToggleFact->getProlongatorFactory((size_t)coarseLevel.Get<int>("RfromPfactory", PFact.get()))));
        // toggle factory sets RfromPfactory to correct index into prolongatorFactory array
        MueLu::DisableMultipleCallCheck check2(actualPFact);
        bool rmode = actualPFact->isRestrictionModeSet();
        actualPFact->setRestrictionMode(true);
        R = coarseLevel.Get<RCP<Matrix>>("R", actualPFact.get());
        actualPFact->setRestrictionMode(rmode);
      } else
        R = Utilities::Transpose(*P, true, label, Tparams);
    } else
      R = coarseLevel.Get<RCP<Matrix>>("RfromPfactory", PFact.get());
  } else
    R = Utilities::Transpose(*P, true, label, Tparams);

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

#endif  // MUELU_RFROMP_OR_TRANSP_DEF_HPP
