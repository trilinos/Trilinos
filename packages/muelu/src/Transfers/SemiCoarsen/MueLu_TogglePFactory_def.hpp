// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TOGGLEPFACTORY_DEF_HPP
#define MUELU_TOGGLEPFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>

#include "MueLu_TogglePFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> TogglePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("toggle: mode");
  SET_VALID_ENTRY("semicoarsen: number of levels");
#undef SET_VALID_ENTRY

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TogglePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  // request/release "P" and coarse level "Nullspace"
  for (std::vector<RCP<const FactoryBase> >::const_iterator it = prolongatorFacts_.begin(); it != prolongatorFacts_.end(); ++it) {
    coarseLevel.DeclareInput("P", (*it).get(), this);  // request/release "P" (dependencies are not affected)
    coarseLevel.DeclareInput("RfromPfactory", (*it).get(), this);
    (*it)->CallDeclareInput(coarseLevel);  // request dependencies
  }
  for (std::vector<RCP<const FactoryBase> >::const_iterator it = ptentFacts_.begin(); it != ptentFacts_.end(); ++it) {
    coarseLevel.DeclareInput("P", (*it).get(), this);  // request/release "Ptent" (dependencies are not affected)
    (*it)->CallDeclareInput(coarseLevel);              // request dependencies
  }
  for (std::vector<RCP<const FactoryBase> >::const_iterator it = nspFacts_.begin(); it != nspFacts_.end(); ++it) {
    coarseLevel.DeclareInput("Nullspace", (*it).get(), this);  // request/release coarse "Nullspace" (dependencies are not affected)
    (*it)->CallDeclareInput(coarseLevel);                      // request dependencies
  }

  // The factory needs the information about the number of z-layers. While this information is
  // provided by the user for the finest level, the factory itself is responsible to provide the
  // corresponding information on the coarser levels. Since a factory cannot be dependent on itself
  // we use the NoFactory class as generator class, but remove the UserData keep flag, such that
  // "NumZLayers" is part of the request/release mechanism.
  // Please note, that this prevents us from having several (independent) CoarsePFactory instances!
  // TODO: allow factory to dependent on self-generated data for TwoLevelFactories -> introduce ExpertRequest/Release in Level
  fineLevel.DeclareInput("NumZLayers", NoFactory::get(), this);
  fineLevel.RemoveKeepFlag("NumZLayers", NoFactory::get(), MueLu::UserData);

  hasDeclaredInput_ = true;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TogglePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Prolongator toggle", coarseLevel);
  std::ostringstream levelstr;
  levelstr << coarseLevel.GetLevelID();

  TEUCHOS_TEST_FOR_EXCEPTION(nspFacts_.size() != prolongatorFacts_.size(), Exceptions::RuntimeError, "MueLu::TogglePFactory::Build: The number of provided prolongator factories and coarse nullspace factories must be identical.");
  TEUCHOS_TEST_FOR_EXCEPTION(nspFacts_.size() != 2, Exceptions::RuntimeError, "MueLu::TogglePFactory::Build: TogglePFactory needs two different transfer operator strategies for toggling.");  // TODO adapt this/weaken this as soon as other toggling strategies are introduced.

  // decision routine which prolongator factory to be used
  int nProlongatorFactory = 0;  // default behavior: use first prolongator in list

  // extract user parameters
  const Teuchos::ParameterList& pL = GetParameterList();
  std::string mode                 = Teuchos::as<std::string>(pL.get<std::string>("toggle: mode"));
  int semicoarsen_levels           = Teuchos::as<int>(pL.get<int>("semicoarsen: number of levels"));

  TEUCHOS_TEST_FOR_EXCEPTION(mode != "semicoarsen", Exceptions::RuntimeError, "MueLu::TogglePFactory::Build: The 'toggle: mode' parameter must be set to 'semicoarsen'. No other mode supported, yet.");

  LO NumZDir = -1;
  if (fineLevel.IsAvailable("NumZLayers", NoFactory::get())) {
    NumZDir = fineLevel.Get<LO>("NumZLayers", NoFactory::get());  // obtain info
    GetOStream(Runtime1) << "Number of layers for semicoarsening: " << NumZDir << std::endl;
  }

  // Make a decision which prolongator to be used.
  if (fineLevel.GetLevelID() >= semicoarsen_levels || NumZDir == 1) {
    nProlongatorFactory = 1;
  } else {
    nProlongatorFactory = 0;
  }

  RCP<Matrix> P                    = Teuchos::null;
  RCP<Matrix> Ptent                = Teuchos::null;
  RCP<MultiVector> coarseNullspace = Teuchos::null;

  // call Build for selected transfer operator
  GetOStream(Runtime0) << "TogglePFactory: call transfer factory: " << (prolongatorFacts_[nProlongatorFactory])->description() << std::endl;
  prolongatorFacts_[nProlongatorFactory]->CallBuild(coarseLevel);
  P                = coarseLevel.Get<RCP<Matrix> >("P", (prolongatorFacts_[nProlongatorFactory]).get());
  RCP<Matrix> R    = Teuchos::null;
  int Rplaceholder = -1;  // Used to indicate that an R matrix has not been produced by a prolongator factory, but
                          // that it is capable of producing one and should be invoked a 2nd time in restrictor mode
                          // (e.g. with PgPFactory).  prolongatorFacts_[Rplaceholder] is factory that can produce R
                          // matrix, which might be later invoked by  MueLu_RfromP_Or_TransP
  if (coarseLevel.IsAvailable("RfromPfactory", (prolongatorFacts_[nProlongatorFactory]).get())) {
    std::string strType = coarseLevel.GetTypeName("RfromPfactory", (prolongatorFacts_[nProlongatorFactory]).get());
    if (strType == "int")
      Rplaceholder = nProlongatorFactory;
    else
      R = coarseLevel.Get<RCP<Matrix> >("RfromPfactory", (prolongatorFacts_[nProlongatorFactory]).get());
    // Need to get R (and set it below) so that TogglePFactory is given credit for creating R
  }
  // do not call "Build" for "Ptent" factory since it should automatically be called recursively
  // through the "Build" call for "P"
  Ptent           = coarseLevel.Get<RCP<Matrix> >("P", (ptentFacts_[nProlongatorFactory]).get());
  coarseNullspace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", (nspFacts_[nProlongatorFactory]).get());

  // Release dependencies of all prolongator and coarse level null spaces
  for (size_t t = 0; t < nspFacts_.size(); ++t) {
    coarseLevel.Release(*(prolongatorFacts_[t]));
    coarseLevel.Release(*(ptentFacts_[t]));
    coarseLevel.Release(*(nspFacts_[t]));
  }

  // store prolongator with this factory identification.
  Set(coarseLevel, "P", P);
  //  Three cases:
  //   1) R already computed and TogglePFactory takes credit for constructing it
  //   2) R not computed but prolongatorFacts_[Rplaceholder] can produce it
  //   3) R not computed  and prolongator can not produce it
  if (R != Teuchos::null)
    Set(coarseLevel, "RfromPfactory", R);
  else if (Rplaceholder != -1)
    Set(coarseLevel, "RfromPfactory", Teuchos::as<int>(Rplaceholder));
  Set(coarseLevel, "Nullspace", coarseNullspace);
  Set(coarseLevel, "Ptent", Ptent);
  Set(coarseLevel, "Chosen P", nProlongatorFactory);
}  // Build()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TogglePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddProlongatorFactory(const RCP<const FactoryBase>& factory) {
  // check if it's a TwoLevelFactoryBase based transfer factory
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast,
                             "MueLu::TogglePFactory::AddProlongatorFactory: Transfer factory is not derived from TwoLevelFactoryBase. "
                             "This is very strange. (Note: you can remove this exception if there's a good reason for)");
  TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_, Exceptions::RuntimeError, "MueLu::TogglePFactory::AddProlongatorFactory: Factory is being added after we have already declared input");
  prolongatorFacts_.push_back(factory);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TogglePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddPtentFactory(const RCP<const FactoryBase>& factory) {
  // check if it's a TwoLevelFactoryBase based transfer factory
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast,
                             "MueLu::TogglePFactory::AddPtentFactory: Transfer factory is not derived from TwoLevelFactoryBase. "
                             "This is very strange. (Note: you can remove this exception if there's a good reason for)");
  TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_, Exceptions::RuntimeError, "MueLu::TogglePFactory::AddPtentFactory: Factory is being added after we have already declared input");
  ptentFacts_.push_back(factory);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TogglePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddCoarseNullspaceFactory(const RCP<const FactoryBase>& factory) {
  // check if it's a TwoLevelFactoryBase based transfer factory
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast,
                             "MueLu::TogglePFactory::AddCoarseNullspaceFactory: Transfer factory is not derived from TwoLevelFactoryBase. Make sure you provide the factory which generates the coarse level nullspace information. Usually this is a prolongator factory."
                             "This is very strange. (Note: you can remove this exception if there's a good reason for)");
  TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_, Exceptions::RuntimeError, "MueLu::TogglePFactory::AddCoarseNullspaceFactory: Factory is being added after we have already declared input");
  nspFacts_.push_back(factory);
}

}  // namespace MueLu

#endif  // MUELU_TOGGLEPFACTORY_DEF_HPP
