// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_TimeMonitor.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

SingleLevelFactoryBase::SingleLevelFactoryBase() = default;

SingleLevelFactoryBase::~SingleLevelFactoryBase() = default;

void SingleLevelFactoryBase::CallBuild(Level& requestedLevel) const {
#ifdef HAVE_MUELU_DEBUG
  // We cannot call Build method twice for the same level, but we can call it multiple times for different levels
  int levelID = requestedLevel.GetLevelID();
  TEUCHOS_TEST_FOR_EXCEPTION((multipleCallCheck_ == ENABLED) && (multipleCallCheckGlobal_ == ENABLED) && (lastLevelID_ == levelID),
                             Exceptions::RuntimeError,
                             this->ShortClassName() << "::Build() called twice for the same level (levelID=" << levelID
                                                    << "). This is likely due to a configuration error, or calling hierarchy setup multiple times "
                                                    << "without resetting debug info through FactoryManager::ResetDebugData().");
  if (multipleCallCheck_ == FIRSTCALL)
    multipleCallCheck_ = ENABLED;

  lastLevelID_ = levelID;
#endif

  RCP<const Teuchos::Comm<int> > comm = requestedLevel.GetComm();
  if (comm.is_null()) {
    // Some factories are called before we constructed Ac, and therefore,
    // before we set the level communicator. For such factories we can get
    // the comm from the previous level, as all processes go there
    RCP<Level>& prevLevel = requestedLevel.GetPreviousLevel();
    if (!prevLevel.is_null())
      comm = prevLevel->GetComm();
  }

  int oldRank = -1;
  if (!comm.is_null())
    oldRank = SetProcRankVerbose(comm->getRank());

  // Synchronization timer
  std::string syncTimer = this->ShortClassName() + ": Build sync (level=" + toString(requestedLevel.GetLevelID()) + ")";
  if (this->timerSync_ && !comm.is_null()) {
    TimeMonitor timer(*this, syncTimer);
    comm->barrier();
  }

  Build(requestedLevel);

  if (this->timerSync_ && !comm.is_null()) {
    TimeMonitor timer(*this, syncTimer);
    comm->barrier();
  }

  if (IsPrint(Test1))
    GetOStream(Test1) << *RemoveFactoriesFromList(GetParameterList()) << std::endl;
  else
    RemoveFactoriesFromList(GetParameterList())->print(GetOStream(Test0), Teuchos::ParameterList::PrintOptions().indent(0).showFlags(true).showDefault(false));

  if (oldRank != -1)
    SetProcRankVerbose(oldRank);
}

void SingleLevelFactoryBase::CallDeclareInput(Level& requestedLevel) const {
  DeclareInput(requestedLevel);
}

}  // namespace MueLu

// TODO: code factorization between SingleLevelFactoryBase and TwoLevelFactoryBase
