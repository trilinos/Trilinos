// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_RAPFACTORY_DEF_HPP
#define MUELU_RAPFACTORY_DEF_HPP

#include <sstream>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <stdexcept>

#include "MueLu_RAPFactory_decl.hpp"

#include "MueLu_Utilities.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Behavior.hpp"
#include "Teuchos_TestForException.hpp"
#include "MueLu_Behavior.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RAPFactory()
  : hasDeclaredInput_(false) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~RAPFactory() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("transpose: use implicit");
  SET_VALID_ENTRY("rap: triple product");
  SET_VALID_ENTRY("rap: fix zero diagonals");
  SET_VALID_ENTRY("rap: fix zero diagonals threshold");
  SET_VALID_ENTRY("rap: fix zero diagonals replacement");
  SET_VALID_ENTRY("rap: relative diagonal floor");
#undef SET_VALID_ENTRY
  validParamList->set<RCP<const FactoryBase> >("A", null, "Generating factory of the matrix A used during the prolongator smoothing process");
  validParamList->set<RCP<const FactoryBase> >("P", null, "Prolongator factory");
  validParamList->set<RCP<const FactoryBase> >("R", null, "Restrictor factory");

  validParamList->set<bool>("CheckMainDiagonal", false, "Check main diagonal for zeros");
  validParamList->set<bool>("RepairMainDiagonal", false, "Repair zeros on main diagonal");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  const Teuchos::ParameterList& pL = GetParameterList();
  if (!pL.get<bool>("transpose: use implicit"))
    Input(coarseLevel, "R");

  Input(fineLevel, "A");
  Input(coarseLevel, "P");

  // call DeclareInput of all user-given transfer factories
  for (std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it)
    (*it)->CallDeclareInput(coarseLevel);

  hasDeclaredInput_ = true;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  RCP<Matrix> Ac;

  {
    FactoryMonitor m(*this, "Computing Ac", coarseLevel);

    TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_ == false, Exceptions::RuntimeError,
                               "MueLu::RAPFactory::Build(): CallDeclareInput has not been called before Build!");

    RCP<Matrix> A = Get<RCP<Matrix> >(fineLevel, "A");
    RCP<Matrix> P = Get<RCP<Matrix> >(coarseLevel, "P"), AP, R;
    // We don't have a valid P (e.g., # global aggregates = 0) so we bail.
    // This level will ultimately be removed in MueLu_Hierarchy_defs.h via a resize()
    if (P.is_null()) {
      Ac = Teuchos::null;
      Set(coarseLevel, "A", Ac);
      return;
    }

    const Teuchos::ParameterList& pL = GetParameterList();
    const bool useImplicit           = pL.get<bool>("transpose: use implicit");
    bool isGPU                       = Node::is_gpu;

    Teuchos::RCP<Teuchos::ParameterList> APparams;
    Teuchos::RCP<Teuchos::ParameterList> RAPparams;
    if (coarseLevel.IsAvailable("AP reuse data", this)) {
      GetOStream(static_cast<MsgType>(Runtime0 | Test)) << "Reusing previous AP data" << std::endl;
      APparams = coarseLevel.Get<RCP<ParameterList> >("AP reuse data", this);
    } else {
      APparams = Teuchos::rcp(new Teuchos::ParameterList());
    }
    if (coarseLevel.IsAvailable("RAP reuse data", this)) {
      GetOStream(static_cast<MsgType>(Runtime0 | Test)) << "Reusing previous RAP data" << std::endl;
      RAPparams = coarseLevel.Get<RCP<ParameterList> >("RAP reuse data", this);
    } else {
      RAPparams = Teuchos::rcp(new Teuchos::ParameterList());
    }
    if (!useImplicit)
      R = Get<RCP<Matrix> >(coarseLevel, "R");
    Utilities::TripleMatrixProduct(R, A, P, Ac, pL, *this, APparams, RAPparams, &coarseLevel);

    if (!Ac.is_null()) {
      std::ostringstream oss;
      oss << "A_" << coarseLevel.GetLevelID();
      Ac->setObjectLabel(oss.str());
    }
    Set(coarseLevel, "A", Ac);
    if (!isGPU) {
      if (!pL.get<bool>("rap: triple product")) {
        TEUCHOS_TEST_FOR_EXCEPTION(!APparams->isParameter("graph"), std::runtime_error, "\"AP reuse data\" does not contain the expected reuse data.");
        Set(coarseLevel, "AP reuse data", APparams);
      }
      {
        TEUCHOS_TEST_FOR_EXCEPTION(!RAPparams->isParameter("graph"), std::runtime_error, "\"RAP reuse data\" does not contain the expected reuse data.");
        Set(coarseLevel, "RAP reuse data", RAPparams);
      }
    }
  }

  if (Behavior::debug())
    MatrixUtils::checkLocalRowMapMatchesColMap(*Ac);

  if (transferFacts_.begin() != transferFacts_.end()) {
    SubFactoryMonitor m(*this, "Projections", coarseLevel);

    // call Build of all user-given transfer factories
    for (std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it) {
      RCP<const FactoryBase> fac = *it;
      GetOStream(Runtime0) << "RAPFactory: call transfer factory: " << fac->description() << std::endl;
      fac->CallBuild(coarseLevel);
      // Coordinates transfer is marginally different from all other operations
      // because it is *optional*, and not required. For instance, we may need
      // coordinates only on level 4 if we start repartitioning from that level,
      // but we don't need them on level 1,2,3. As our current Hierarchy setup
      // assumes propagation of dependencies only through three levels, this
      // means that we need to rely on other methods to propagate optional data.
      //
      // The method currently used is through RAP transfer factories, which are
      // simply factories which are called at the end of RAP with a single goal:
      // transfer some fine data to coarser level. Because these factories are
      // kind of outside of the mainline factories, they behave different. In
      // particular, we call their Build method explicitly, rather than through
      // Get calls. This difference is significant, as the Get call is smart
      // enough to know when to release all factory dependencies, and Build is
      // dumb. This led to the following CoordinatesTransferFactory sequence:
      // 1. Request level 0
      // 2. Request level 1
      // 3. Request level 0
      // 4. Release level 0
      // 5. Release level 1
      //
      // The problem is missing "6. Release level 0". Because it was missing,
      // we had outstanding request on "Coordinates", "Aggregates" and
      // "CoarseMap" on level 0.
      //
      // This was fixed by explicitly calling Release on transfer factories in
      // RAPFactory. I am still unsure how exactly it works, but now we have
      // clear data requests for all levels.
      coarseLevel.Release(*fac);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddTransferFactory(const RCP<const FactoryBase>& factory) {
  // check if it's a TwoLevelFactoryBase based transfer factory
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast,
                             "MueLu::RAPFactory::AddTransferFactory: Transfer factory is not derived from TwoLevelFactoryBase. "
                             "This is very strange. (Note: you can remove this exception if there's a good reason for)");
  TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_, Exceptions::RuntimeError, "MueLu::RAPFactory::AddTransferFactory: Factory is being added after we have already declared input");
  transferFacts_.push_back(factory);
}

}  // namespace MueLu

#define MUELU_RAPFACTORY_SHORT
#endif  // MUELU_RAPFACTORY_DEF_HPP
