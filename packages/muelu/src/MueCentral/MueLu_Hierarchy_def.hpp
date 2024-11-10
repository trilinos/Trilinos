// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_HIERARCHY_DEF_HPP
#define MUELU_HIERARCHY_DEF_HPP

#include <time.h>

#include <algorithm>
#include <sstream>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_Hierarchy_decl.hpp"

#include "MueLu_FactoryManager.hpp"
#include "MueLu_HierarchyUtils.hpp"
#include "MueLu_TopRAPFactory.hpp"
#include "MueLu_TopSmootherFactory.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_SmootherBase.hpp"

#include "Teuchos_TimeMonitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Hierarchy()
  : maxCoarseSize_(GetDefaultMaxCoarseSize())
  , implicitTranspose_(GetDefaultImplicitTranspose())
  , fuseProlongationAndUpdate_(GetDefaultFuseProlongationAndUpdate())
  , doPRrebalance_(GetDefaultPRrebalance())
  , doPRViaCopyrebalance_(false)
  , isPreconditioner_(true)
  , Cycle_(GetDefaultCycle())
  , WCycleStartLevel_(0)
  , scalingFactor_(Teuchos::ScalarTraits<double>::one())
  , lib_(Xpetra::UseTpetra)
  , isDumpingEnabled_(false)
  , dumpLevel_(-2)
  , rate_(-1)
  , sizeOfAllocatedLevelMultiVectors_(0) {
  AddLevel(rcp(new Level));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Hierarchy(const std::string& label)
  : Hierarchy() {
  setObjectLabel(label);
  Levels_[0]->setObjectLabel(label);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Hierarchy(const RCP<Matrix>& A)
  : maxCoarseSize_(GetDefaultMaxCoarseSize())
  , implicitTranspose_(GetDefaultImplicitTranspose())
  , fuseProlongationAndUpdate_(GetDefaultFuseProlongationAndUpdate())
  , doPRrebalance_(GetDefaultPRrebalance())
  , doPRViaCopyrebalance_(false)
  , isPreconditioner_(true)
  , Cycle_(GetDefaultCycle())
  , WCycleStartLevel_(0)
  , scalingFactor_(Teuchos::ScalarTraits<double>::one())
  , isDumpingEnabled_(false)
  , dumpLevel_(-2)
  , rate_(-1)
  , sizeOfAllocatedLevelMultiVectors_(0) {
  lib_ = A->getDomainMap()->lib();

  RCP<Level> Finest = rcp(new Level);
  AddLevel(Finest);

  Finest->Set("A", A);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Hierarchy(const RCP<Matrix>& A, const std::string& label)
  : Hierarchy(A) {
  setObjectLabel(label);
  Levels_[0]->setObjectLabel(label);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~Hierarchy() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddLevel(const RCP<Level>& level) {
  int levelID = LastLevelID() + 1;  // ID of the inserted level

  if (level->GetLevelID() != -1 && (level->GetLevelID() != levelID))
    GetOStream(Warnings1) << "Hierarchy::AddLevel(): Level with ID=" << level->GetLevelID() << " have been added at the end of the hierarchy\n but its ID have been redefined"
                          << " because last level ID of the hierarchy was " << LastLevelID() << "." << std::endl;

  Levels_.push_back(level);
  level->SetLevelID(levelID);
  level->setlib(lib_);

  level->SetPreviousLevel((levelID == 0) ? Teuchos::null : Levels_[LastLevelID() - 1]);
  level->setObjectLabel(this->getObjectLabel());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddNewLevel() {
  RCP<Level> newLevel = Levels_[LastLevelID()]->Build();  // new coarse level, using copy constructor
  newLevel->setlib(lib_);
  this->AddLevel(newLevel);  // add to hierarchy
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Level>& Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetLevel(const int levelID) {
  TEUCHOS_TEST_FOR_EXCEPTION(levelID < 0 || levelID > LastLevelID(), Exceptions::RuntimeError,
                             "MueLu::Hierarchy::GetLevel(): invalid input parameter value: LevelID = " << levelID);
  return Levels_[levelID];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetNumLevels() const {
  return Levels_.size();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetGlobalNumLevels() const {
  RCP<Operator> A                     = Levels_[0]->template Get<RCP<Operator> >("A");
  RCP<const Teuchos::Comm<int> > comm = A->getDomainMap()->getComm();

  int numLevels = GetNumLevels();
  int numGlobalLevels;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, numLevels, Teuchos::ptr(&numGlobalLevels));

  return numGlobalLevels;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetOperatorComplexity() const {
  double totalNnz = 0, lev0Nnz = 1;
  for (int i = 0; i < GetNumLevels(); ++i) {
    TEUCHOS_TEST_FOR_EXCEPTION(!(Levels_[i]->IsAvailable("A")), Exceptions::RuntimeError,
                               "Operator complexity cannot be calculated because A is unavailable on level " << i);
    RCP<Operator> A = Levels_[i]->template Get<RCP<Operator> >("A");
    if (A.is_null())
      break;

    RCP<Matrix> Am = rcp_dynamic_cast<Matrix>(A);
    if (Am.is_null()) {
      GetOStream(Warnings0) << "Some level operators are not matrices, operator complexity calculation aborted" << std::endl;
      return 0.0;
    }

    totalNnz += as<double>(Am->getGlobalNumEntries());
    if (i == 0)
      lev0Nnz = totalNnz;
  }
  return totalNnz / lev0Nnz;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetSmootherComplexity() const {
  double node_sc = 0, global_sc = 0;
  double a0_nnz        = 0;
  const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();
  // Get cost of fine matvec
  if (GetNumLevels() <= 0) return -1.0;
  if (!Levels_[0]->IsAvailable("A")) return -1.0;

  RCP<Operator> A = Levels_[0]->template Get<RCP<Operator> >("A");
  if (A.is_null()) return -1.0;
  RCP<Matrix> Am = rcp_dynamic_cast<Matrix>(A);
  if (Am.is_null()) return -1.0;
  a0_nnz = as<double>(Am->getGlobalNumEntries());

  // Get smoother complexity at each level
  for (int i = 0; i < GetNumLevels(); ++i) {
    size_t level_sc = 0;
    if (!Levels_[i]->IsAvailable("PreSmoother")) continue;
    RCP<SmootherBase> S = Levels_[i]->template Get<RCP<SmootherBase> >("PreSmoother");
    if (S.is_null()) continue;
    level_sc = S->getNodeSmootherComplexity();
    if (level_sc == INVALID) {
      global_sc = -1.0;
      break;
    }

    node_sc += as<double>(level_sc);
  }

  double min_sc                       = 0.0;
  RCP<const Teuchos::Comm<int> > comm = A->getDomainMap()->getComm();
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, node_sc, Teuchos::ptr(&global_sc));
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, node_sc, Teuchos::ptr(&min_sc));

  if (min_sc < 0.0)
    return -1.0;
  else
    return global_sc / a0_nnz;
}

// Coherence checks todo in Setup() (using an helper function):
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CheckLevel(Level& level, int levelID) {
  TEUCHOS_TEST_FOR_EXCEPTION(level.lib() != lib_, Exceptions::RuntimeError,
                             "MueLu::Hierarchy::CheckLevel(): wrong underlying linear algebra library.");
  TEUCHOS_TEST_FOR_EXCEPTION(level.GetLevelID() != levelID, Exceptions::RuntimeError,
                             "MueLu::Hierarchy::CheckLevel(): wrong level ID");
  TEUCHOS_TEST_FOR_EXCEPTION(levelID != 0 && level.GetPreviousLevel() != Levels_[levelID - 1], Exceptions::RuntimeError,
                             "MueLu::Hierarchy::Setup(): wrong level parent");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetMatvecParams(RCP<ParameterList> matvecParams) {
  for (int i = 0; i < GetNumLevels(); ++i) {
    RCP<Level> level = Levels_[i];
    if (level->IsAvailable("A")) {
      RCP<Operator> Aop = level->Get<RCP<Operator> >("A");
      RCP<Matrix> A     = rcp_dynamic_cast<Matrix>(Aop);
      if (!A.is_null()) {
        RCP<const Import> xpImporter = A->getCrsGraph()->getImporter();
        if (!xpImporter.is_null())
          xpImporter->setDistributorParameters(matvecParams);
        RCP<const Export> xpExporter = A->getCrsGraph()->getExporter();
        if (!xpExporter.is_null())
          xpExporter->setDistributorParameters(matvecParams);
      }
    }
    if (level->IsAvailable("P")) {
      RCP<Matrix> P                = level->Get<RCP<Matrix> >("P");
      RCP<const Import> xpImporter = P->getCrsGraph()->getImporter();
      if (!xpImporter.is_null())
        xpImporter->setDistributorParameters(matvecParams);
      RCP<const Export> xpExporter = P->getCrsGraph()->getExporter();
      if (!xpExporter.is_null())
        xpExporter->setDistributorParameters(matvecParams);
    }
    if (level->IsAvailable("R")) {
      RCP<Matrix> R                = level->Get<RCP<Matrix> >("R");
      RCP<const Import> xpImporter = R->getCrsGraph()->getImporter();
      if (!xpImporter.is_null())
        xpImporter->setDistributorParameters(matvecParams);
      RCP<const Export> xpExporter = R->getCrsGraph()->getExporter();
      if (!xpExporter.is_null())
        xpExporter->setDistributorParameters(matvecParams);
    }
    if (level->IsAvailable("Importer")) {
      RCP<const Import> xpImporter = level->Get<RCP<const Import> >("Importer");
      if (!xpImporter.is_null())
        xpImporter->setDistributorParameters(matvecParams);
    }
  }
}

// The function uses three managers: fine, coarse and next coarse
// We construct the data for the coarse level, and do requests for the next coarse
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(int coarseLevelID,
                                                                 const RCP<const FactoryManagerBase> fineLevelManager,
                                                                 const RCP<const FactoryManagerBase> coarseLevelManager,
                                                                 const RCP<const FactoryManagerBase> nextLevelManager) {
  // Use PrintMonitor/TimerMonitor instead of just a FactoryMonitor to print "Level 0" instead of Hierarchy(0)
  // Print is done after the requests for next coarse level

  TEUCHOS_TEST_FOR_EXCEPTION(LastLevelID() < coarseLevelID, Exceptions::RuntimeError,
                             "MueLu::Hierarchy:Setup(): level " << coarseLevelID << " (specified by coarseLevelID argument) "
                                                                                    "must be built before calling this function.");

  Level& level = *Levels_[coarseLevelID];

  std::string label = FormattingHelper::getColonLabel(level.getObjectLabel());
  TimeMonitor m1(*this, label + this->ShortClassName() + ": " + "Setup (total)");
  TimeMonitor m2(*this, label + this->ShortClassName() + ": " + "Setup" + " (total, level=" + Teuchos::toString(coarseLevelID) + ")");

  // TODO: pass coarseLevelManager by reference
  TEUCHOS_TEST_FOR_EXCEPTION(coarseLevelManager == Teuchos::null, Exceptions::RuntimeError,
                             "MueLu::Hierarchy::Setup(): argument coarseLevelManager cannot be null");

  typedef MueLu::TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> TopRAPFactory;
  typedef MueLu::TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> TopSmootherFactory;

  if (levelManagers_.size() < coarseLevelID + 1)
    levelManagers_.resize(coarseLevelID + 1);
  levelManagers_[coarseLevelID] = coarseLevelManager;

  bool isFinestLevel = (fineLevelManager.is_null());
  bool isLastLevel   = (nextLevelManager.is_null());

  int oldRank = -1;
  if (isFinestLevel) {
    RCP<Operator> A                     = level.Get<RCP<Operator> >("A");
    RCP<const Map> domainMap            = A->getDomainMap();
    RCP<const Teuchos::Comm<int> > comm = domainMap->getComm();

    // Initialize random seed for reproducibility
    Utilities::SetRandomSeed(*comm);

    // Record the communicator on the level (used for timers sync)
    level.SetComm(comm);
    oldRank = SetProcRankVerbose(comm->getRank());

    // Set the Hierarchy library to match that of the finest level matrix,
    // even if it was already set
    lib_ = domainMap->lib();
    level.setlib(lib_);

  } else {
    // Permeate library to a coarser level
    level.setlib(lib_);

    Level& prevLevel = *Levels_[coarseLevelID - 1];
    oldRank          = SetProcRankVerbose(prevLevel.GetComm()->getRank());
  }

  CheckLevel(level, coarseLevelID);

  // Attach FactoryManager to the fine level
  RCP<SetFactoryManager> SFMFine;
  if (!isFinestLevel)
    SFMFine = rcp(new SetFactoryManager(Levels_[coarseLevelID - 1], fineLevelManager));

  if (isFinestLevel && Levels_[coarseLevelID]->IsAvailable("Coordinates"))
    ReplaceCoordinateMap(*Levels_[coarseLevelID]);

  // Attach FactoryManager to the coarse level
  SetFactoryManager SFMCoarse(Levels_[coarseLevelID], coarseLevelManager);

  if (isDumpingEnabled_ && (dumpLevel_ == 0 || dumpLevel_ == -1) && coarseLevelID == 1)
    DumpCurrentGraph(0);

  RCP<TopSmootherFactory> coarseFact;
  RCP<TopSmootherFactory> smootherFact = rcp(new TopSmootherFactory(coarseLevelManager, "Smoother"));

  int nextLevelID = coarseLevelID + 1;

  RCP<SetFactoryManager> SFMNext;
  if (isLastLevel == false) {
    // We are not at the coarsest level, so there is going to be another level ("next coarse") after this one ("coarse")
    if (nextLevelID > LastLevelID())
      AddNewLevel();
    CheckLevel(*Levels_[nextLevelID], nextLevelID);

    // Attach FactoryManager to the next level (level after coarse)
    SFMNext = rcp(new SetFactoryManager(Levels_[nextLevelID], nextLevelManager));
    Levels_[nextLevelID]->Request(TopRAPFactory(coarseLevelManager, nextLevelManager));

    // Do smoother requests here. We don't know whether this is going to be
    // the coarsest level or not, but we need to DeclareInput before we call
    // coarseRAPFactory.Build(), otherwise some stuff may be erased after
    // level releases
    level.Request(*smootherFact);

  } else {
    // Similar to smoother above, do the coarse solver request here. We don't
    // know whether this is going to be the coarsest level or not, but we
    // need to DeclareInput before we call coarseRAPFactory.Build(),
    // otherwise some stuff may be erased after level releases. This is
    // actually evident on ProjectorSmoother. It requires both "A" and
    // "Nullspace". However, "Nullspace" is erased after all releases, so if
    // we call the coarse factory request after RAP build we would not have
    // any data, and cannot get it as we don't have previous managers. The
    // typical trace looks like this:
    //
    // MueLu::Level(0)::GetFactory(Aggregates, 0): No FactoryManager
    //   during request for data "     Aggregates" on level 0 by factory TentativePFactory
    //   during request for data "              P" on level 1 by factory EminPFactory
    //   during request for data "              P" on level 1 by factory TransPFactory
    //   during request for data "              R" on level 1 by factory RAPFactory
    //   during request for data "              A" on level 1 by factory TentativePFactory
    //   during request for data "      Nullspace" on level 2 by factory NullspaceFactory
    //   during request for data "      Nullspace" on level 2 by factory NullspacePresmoothFactory
    //   during request for data "      Nullspace" on level 2 by factory ProjectorSmoother
    //   during request for data "    PreSmoother" on level 2 by factory NoFactory
    if (coarseFact.is_null())
      coarseFact = rcp(new TopSmootherFactory(coarseLevelManager, "CoarseSolver"));
    level.Request(*coarseFact);
  }

  GetOStream(Runtime0) << std::endl;
  PrintMonitor m0(*this, "Level " + Teuchos::toString(coarseLevelID), static_cast<MsgType>(Runtime0 | Test));

  // Build coarse level hierarchy
  RCP<Operator> Ac = Teuchos::null;
  TopRAPFactory coarseRAPFactory(fineLevelManager, coarseLevelManager);

  if (level.IsAvailable("A")) {
    Ac = level.Get<RCP<Operator> >("A");
  } else if (!isFinestLevel) {
    // We only build here, the release is done later
    coarseRAPFactory.Build(*level.GetPreviousLevel(), level);
  }

  bool setLastLevelviaMaxCoarseSize = false;
  if (level.IsAvailable("A"))
    Ac = level.Get<RCP<Operator> >("A");
  RCP<Matrix> Acm = rcp_dynamic_cast<Matrix>(Ac);

  // Record the communicator on the level
  if (!Ac.is_null())
    level.SetComm(Ac->getDomainMap()->getComm());

  // Test if we reach the end of the hierarchy
  bool isOrigLastLevel = isLastLevel;
  if (isLastLevel) {
    // Last level as we have achieved the max limit
    isLastLevel = true;

  } else if (Ac.is_null()) {
    // Last level for this processor, as it does not belong to the next
    // subcommunicator. Other processors may continue working on the
    // hierarchy
    isLastLevel = true;

  } else {
    if (!Acm.is_null() && Acm->getGlobalNumRows() <= maxCoarseSize_) {
      // Last level as the size of the coarse matrix became too small
      GetOStream(Runtime0) << "Max coarse size (<= " << maxCoarseSize_ << ") achieved" << std::endl;
      isLastLevel = true;
      if (Acm->getGlobalNumRows() != 0) setLastLevelviaMaxCoarseSize = true;
    }
  }

  if (!Ac.is_null() && !isFinestLevel) {
    RCP<Operator> A = Levels_[coarseLevelID - 1]->template Get<RCP<Operator> >("A");
    RCP<Matrix> Am  = rcp_dynamic_cast<Matrix>(A);

    const double maxCoarse2FineRatio = 0.8;
    if (!Acm.is_null() && !Am.is_null() && Acm->getGlobalNumRows() > maxCoarse2FineRatio * Am->getGlobalNumRows()) {
      // We could abort here, but for now we simply notify user.
      // Couple of additional points:
      //   - if repartitioning is delayed until level K, but the aggregation
      //     procedure stagnates between levels K-1 and K.  In this case,
      //     repartitioning could enable faster coarsening once again, but the
      //     hierarchy construction will abort due to the stagnation check.
      //   - if the matrix is small enough, we could move it to one processor.
      GetOStream(Warnings0) << "Aggregation stagnated. Please check your matrix and/or adjust your configuration file."
                            << "Possible fixes:\n"
                            << "  - reduce the maximum number of levels\n"
                            << "  - enable repartitioning\n"
                            << "  - increase the minimum coarse size." << std::endl;
    }
  }

  if (isLastLevel) {
    if (!isOrigLastLevel) {
      // We did not expect to finish this early so we did request a smoother.
      // We need a coarse solver instead. Do the magic.
      level.Release(*smootherFact);
      if (coarseFact.is_null())
        coarseFact = rcp(new TopSmootherFactory(coarseLevelManager, "CoarseSolver"));
      level.Request(*coarseFact);
    }

    // Do the actual build, if we have any data.
    // NOTE: this is not a great check, we may want to call Build() regardless.
    if (!Ac.is_null())
      coarseFact->Build(level);

    // Once the dirty deed is done, release stuff. The smoother has already
    // been released.
    level.Release(*coarseFact);

  } else {
    // isLastLevel = false => isOrigLastLevel = false, meaning that we have
    // requested the smoother. Now we need to build it and to release it.
    // We don't need to worry about the coarse solver, as we didn't request it.
    if (!Ac.is_null())
      smootherFact->Build(level);

    level.Release(*smootherFact);
  }

  if (isLastLevel == true) {
    int actualNumLevels = nextLevelID;
    if (isOrigLastLevel == false) {
      // Earlier in the function, we constructed the next coarse level, and requested data for the that level,
      // assuming that we are not at the coarsest level. Now, we changed our mind, so we have to release those.
      Levels_[nextLevelID]->Release(TopRAPFactory(coarseLevelManager, nextLevelManager));

      // We truncate/resize the hierarchy and possibly remove the last created level if there is
      // something wrong with it as indicated by its P not being valid. This might happen
      // if the global number of aggregates turns out to be zero

      if (!setLastLevelviaMaxCoarseSize) {
        if (Levels_[nextLevelID - 1]->IsAvailable("P")) {
          if (Levels_[nextLevelID - 1]->template Get<RCP<Matrix> >("P") == Teuchos::null) actualNumLevels = nextLevelID - 1;
        } else
          actualNumLevels = nextLevelID - 1;
      }
    }
    if (actualNumLevels == nextLevelID - 1) {
      // Didn't expect to finish early so we requested smoother but need coarse solver instead.
      Levels_[nextLevelID - 2]->Release(*smootherFact);

      if (Levels_[nextLevelID - 2]->IsAvailable("PreSmoother")) Levels_[nextLevelID - 2]->RemoveKeepFlag("PreSmoother", NoFactory::get());
      if (Levels_[nextLevelID - 2]->IsAvailable("PostSmoother")) Levels_[nextLevelID - 2]->RemoveKeepFlag("PostSmoother", NoFactory::get());
      if (coarseFact.is_null())
        coarseFact = rcp(new TopSmootherFactory(coarseLevelManager, "CoarseSolver"));
      Levels_[nextLevelID - 2]->Request(*coarseFact);
      if (!(Levels_[nextLevelID - 2]->template Get<RCP<Matrix> >("A").is_null()))
        coarseFact->Build(*(Levels_[nextLevelID - 2]));
      Levels_[nextLevelID - 2]->Release(*coarseFact);
    }
    Levels_.resize(actualNumLevels);
  }

  // I think this is the proper place for graph so that it shows every dependence
  if (isDumpingEnabled_ && ((dumpLevel_ > 0 && coarseLevelID == dumpLevel_) || dumpLevel_ == -1))
    DumpCurrentGraph(coarseLevelID);

  if (!isFinestLevel) {
    // Release the hierarchy data
    // We release so late to help blocked solvers, as the smoothers for them need A blocks
    // which we construct in RAPFactory
    level.Release(coarseRAPFactory);
  }

  if (oldRank != -1)
    SetProcRankVerbose(oldRank);

  return isLastLevel;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupRe() {
  int numLevels = Levels_.size();
  TEUCHOS_TEST_FOR_EXCEPTION(levelManagers_.size() != numLevels, Exceptions::RuntimeError,
                             "Hierarchy::SetupRe: " << Levels_.size() << " levels, but  " << levelManagers_.size() << " level factory managers");

  const int startLevel = 0;
  Clear(startLevel);

#ifdef HAVE_MUELU_DEBUG
  // Reset factories' data used for debugging
  for (int i = 0; i < numLevels; i++)
    levelManagers_[i]->ResetDebugData();

#endif

  int levelID;
  for (levelID = startLevel; levelID < numLevels;) {
    bool r = Setup(levelID,
                   (levelID != 0 ? levelManagers_[levelID - 1] : Teuchos::null),
                   levelManagers_[levelID],
                   (levelID + 1 != numLevels ? levelManagers_[levelID + 1] : Teuchos::null));
    levelID++;
    if (r) break;
  }
  // We may construct fewer levels for some reason, make sure we continue
  // doing that in the future
  Levels_.resize(levelID);
  levelManagers_.resize(levelID);

  int sizeOfVecs = sizeOfAllocatedLevelMultiVectors_;

  AllocateLevelMultiVectors(sizeOfVecs, true);

  // since the # of levels, etc. may have changed, force re-determination of description during next call to description()
  ResetDescription();

  describe(GetOStream(Statistics0), GetVerbLevel());

  CheckForEmptySmoothersAndCoarseSolve();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(const FactoryManagerBase& manager, int startLevel, int numDesiredLevels) {
  // Use MueLu::BaseClass::description() to avoid printing "{numLevels = 1}" (numLevels is increasing...)
  PrintMonitor m0(*this, "Setup (" + this->MueLu::BaseClass::description() + ")", Runtime0);

  Clear(startLevel);

  // Check Levels_[startLevel] exists.
  TEUCHOS_TEST_FOR_EXCEPTION(Levels_.size() <= startLevel, Exceptions::RuntimeError,
                             "MueLu::Hierarchy::Setup(): fine level (" << startLevel << ") does not exist");

  TEUCHOS_TEST_FOR_EXCEPTION(numDesiredLevels <= 0, Exceptions::RuntimeError,
                             "Constructing non-positive (" << numDesiredLevels << ") number of levels does not make sense.");

  // Check for fine level matrix A
  TEUCHOS_TEST_FOR_EXCEPTION(!Levels_[startLevel]->IsAvailable("A"), Exceptions::RuntimeError,
                             "MueLu::Hierarchy::Setup(): fine level (" << startLevel << ") has no matrix A! "
                                                                                        "Set fine level matrix A using Level.Set()");

  RCP<Operator> A = Levels_[startLevel]->template Get<RCP<Operator> >("A");
  lib_            = A->getDomainMap()->lib();

  if (IsPrint(Statistics2)) {
    RCP<Matrix> Amat = rcp_dynamic_cast<Matrix>(A);

    if (!Amat.is_null()) {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo", true);

      GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*Amat, "A0", params);
    } else {
      GetOStream(Warnings1) << "Fine level operator is not a matrix, statistics are not available" << std::endl;
    }
  }

  RCP<const FactoryManagerBase> rcpmanager = rcpFromRef(manager);

  const int lastLevel = startLevel + numDesiredLevels - 1;
  GetOStream(Runtime0) << "Setup loop: startLevel = " << startLevel << ", lastLevel = " << lastLevel
                       << " (stop if numLevels = " << numDesiredLevels << " or Ac.size() < " << maxCoarseSize_ << ")" << std::endl;

  // Setup multigrid levels
  int iLevel = 0;
  if (numDesiredLevels == 1) {
    iLevel = 0;
    Setup(startLevel, Teuchos::null, rcpmanager, Teuchos::null);  // setup finest==coarsest level (first and last managers are Teuchos::null)

  } else {
    bool bIsLastLevel = Setup(startLevel, Teuchos::null, rcpmanager, rcpmanager);  // setup finest level (level 0) (first manager is Teuchos::null)
    if (bIsLastLevel == false) {
      for (iLevel = startLevel + 1; iLevel < lastLevel; iLevel++) {
        bIsLastLevel = Setup(iLevel, rcpmanager, rcpmanager, rcpmanager);  // setup intermediate levels
        if (bIsLastLevel == true)
          break;
      }
      if (bIsLastLevel == false)
        Setup(lastLevel, rcpmanager, rcpmanager, Teuchos::null);  // setup coarsest level (last manager is Teuchos::null)
    }
  }

  // TODO: some check like this should be done at the beginning of the routine
  TEUCHOS_TEST_FOR_EXCEPTION(iLevel != Levels_.size() - 1, Exceptions::RuntimeError,
                             "MueLu::Hierarchy::Setup(): number of level");

  // TODO: this is not exception safe: manager will still hold default
  // factories if you exit this function with an exception
  manager.Clean();

  describe(GetOStream(Statistics0), GetVerbLevel());

  CheckForEmptySmoothersAndCoarseSolve();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CheckForEmptySmoothersAndCoarseSolve() {
  for (LO levelNo = 0; levelNo < as<LO>(Levels_.size()); ++levelNo) {
    auto level = Levels_[levelNo];
    if ((!level->IsAvailable("PreSmoother")) && (!level->IsAvailable("PostSmoother")))
      GetOStream(Warnings1) << "No " << (levelNo == as<LO>(Levels_.size()) - 1 ? "coarse grid solver" : "smoother") << " on level " << level->GetLevelID() << std::endl;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Clear(int startLevel) {
  if (startLevel < GetNumLevels())
    GetOStream(Runtime0) << "Clearing old data (if any)" << std::endl;

  for (int iLevel = startLevel; iLevel < GetNumLevels(); iLevel++)
    Levels_[iLevel]->Clear();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ExpertClear() {
  GetOStream(Runtime0) << "Clearing old data (expert)" << std::endl;
  for (int iLevel = 0; iLevel < GetNumLevels(); iLevel++)
    Levels_[iLevel]->ExpertClear();
}

#if defined(HAVE_MUELU_EXPERIMENTAL) && defined(HAVE_MUELU_ADDITIVE_VARIANT)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ConvergenceStatus Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Iterate(const MultiVector& B, MultiVector& X, ConvData conv,
                                                                                bool InitialGuessIsZero, LO startLevel) {
  LO nIts           = conv.maxIts_;
  MagnitudeType tol = conv.tol_;

  std::string prefix       = this->ShortClassName() + ": ";
  std::string levelSuffix  = " (level=" + toString(startLevel) + ")";
  std::string levelSuffix1 = " (level=" + toString(startLevel + 1) + ")";

  using namespace Teuchos;
  RCP<Time> CompTime              = Teuchos::TimeMonitor::getNewCounter(prefix + "Computational Time (total)");
  RCP<Time> Concurrent            = Teuchos::TimeMonitor::getNewCounter(prefix + "Concurrent portion");
  RCP<Time> ApplyR                = Teuchos::TimeMonitor::getNewCounter(prefix + "R: Computational Time");
  RCP<Time> ApplyPbar             = Teuchos::TimeMonitor::getNewCounter(prefix + "Pbar: Computational Time");
  RCP<Time> CompFine              = Teuchos::TimeMonitor::getNewCounter(prefix + "Fine: Computational Time");
  RCP<Time> CompCoarse            = Teuchos::TimeMonitor::getNewCounter(prefix + "Coarse: Computational Time");
  RCP<Time> ApplySum              = Teuchos::TimeMonitor::getNewCounter(prefix + "Sum: Computational Time");
  RCP<Time> Synchronize_beginning = Teuchos::TimeMonitor::getNewCounter(prefix + "Synchronize_beginning");
  RCP<Time> Synchronize_center    = Teuchos::TimeMonitor::getNewCounter(prefix + "Synchronize_center");
  RCP<Time> Synchronize_end       = Teuchos::TimeMonitor::getNewCounter(prefix + "Synchronize_end");

  RCP<Level> Fine = Levels_[0];
  RCP<Level> Coarse;

  RCP<Operator> A                                      = Fine->Get<RCP<Operator> >("A");
  Teuchos::RCP<const Teuchos::Comm<int> > communicator = A->getDomainMap()->getComm();

  // Synchronize_beginning->start();
  // communicator->barrier();
  // Synchronize_beginning->stop();

  CompTime->start();

  SC one = STS::one(), zero = STS::zero();

  bool zeroGuess = InitialGuessIsZero;

  // =======   UPFRONT DEFINITION OF COARSE VARIABLES ===========

  // RCP<const Map> origMap;
  RCP<Operator> P;
  RCP<Operator> Pbar;
  RCP<Operator> R;
  RCP<MultiVector> coarseRhs, coarseX;
  RCP<Operator> Ac;
  RCP<SmootherBase> preSmoo_coarse, postSmoo_coarse;
  bool emptyCoarseSolve              = true;
  RCP<MultiVector> coarseX_prolonged = MultiVectorFactory::Build(X.getMap(), X.getNumVectors(), true);

  RCP<const Import> importer;

  if (Levels_.size() > 1) {
    Coarse = Levels_[1];
    if (Coarse->IsAvailable("Importer"))
      importer = Coarse->Get<RCP<const Import> >("Importer");

    R = Coarse->Get<RCP<Operator> >("R");
    P = Coarse->Get<RCP<Operator> >("P");

    // if(Coarse->IsAvailable("Pbar"))
    Pbar = Coarse->Get<RCP<Operator> >("Pbar");

    coarseRhs = MultiVectorFactory::Build(R->getRangeMap(), B.getNumVectors(), true);

    Ac = Coarse->Get<RCP<Operator> >("A");

    ApplyR->start();
    R->apply(B, *coarseRhs, Teuchos::NO_TRANS, one, zero);
    // P->apply(B, *coarseRhs, Teuchos::TRANS, one, zero);
    ApplyR->stop();

    if (doPRrebalance_ || importer.is_null()) {
      coarseX = MultiVectorFactory::Build(coarseRhs->getMap(), X.getNumVectors(), true);

    } else {
      RCP<TimeMonitor> ITime      = rcp(new TimeMonitor(*this, prefix + "Solve : import (total)", Timings0));
      RCP<TimeMonitor> ILevelTime = rcp(new TimeMonitor(*this, prefix + "Solve : import" + levelSuffix1, Timings0));

      // Import: range map of R --> domain map of rebalanced Ac (before subcomm replacement)
      RCP<MultiVector> coarseTmp = MultiVectorFactory::Build(importer->getTargetMap(), coarseRhs->getNumVectors());
      coarseTmp->doImport(*coarseRhs, *importer, Xpetra::INSERT);
      coarseRhs.swap(coarseTmp);

      coarseX = MultiVectorFactory::Build(importer->getTargetMap(), X.getNumVectors(), true);
    }

    if (Coarse->IsAvailable("PreSmoother"))
      preSmoo_coarse = Coarse->Get<RCP<SmootherBase> >("PreSmoother");
    if (Coarse->IsAvailable("PostSmoother"))
      postSmoo_coarse = Coarse->Get<RCP<SmootherBase> >("PostSmoother");
  }

  // ==========================================================

  MagnitudeType prevNorm = STS::magnitude(STS::one()), curNorm = STS::magnitude(STS::one());
  rate_ = 1.0;

  for (LO i = 1; i <= nIts; i++) {
#ifdef HAVE_MUELU_DEBUG
    if (A->getDomainMap()->isCompatible(*(X.getMap())) == false) {
      std::ostringstream ss;
      ss << "Level " << startLevel << ": level A's domain map is not compatible with X";
      throw Exceptions::Incompatible(ss.str());
    }

    if (A->getRangeMap()->isCompatible(*(B.getMap())) == false) {
      std::ostringstream ss;
      ss << "Level " << startLevel << ": level A's range map is not compatible with B";
      throw Exceptions::Incompatible(ss.str());
    }
#endif
  }

  bool emptyFineSolve = true;

  RCP<MultiVector> fineX;
  fineX = MultiVectorFactory::Build(X.getMap(), X.getNumVectors(), true);

  // Synchronize_center->start();
  // communicator->barrier();
  // Synchronize_center->stop();

  Concurrent->start();

  // NOTE: we need to check using IsAvailable before Get here to avoid building default smoother
  if (Fine->IsAvailable("PreSmoother")) {
    RCP<SmootherBase> preSmoo = Fine->Get<RCP<SmootherBase> >("PreSmoother");
    CompFine->start();
    preSmoo->Apply(*fineX, B, zeroGuess);
    CompFine->stop();
    emptyFineSolve = false;
  }
  if (Fine->IsAvailable("PostSmoother")) {
    RCP<SmootherBase> postSmoo = Fine->Get<RCP<SmootherBase> >("PostSmoother");
    CompFine->start();
    postSmoo->Apply(*fineX, B, zeroGuess);
    CompFine->stop();

    emptyFineSolve = false;
  }
  if (emptyFineSolve == true) {
    // Fine grid smoother is identity
    fineX->update(one, B, zero);
  }

  if (Levels_.size() > 1) {
    // NOTE: we need to check using IsAvailable before Get here to avoid building default smoother
    if (Coarse->IsAvailable("PreSmoother")) {
      CompCoarse->start();
      preSmoo_coarse->Apply(*coarseX, *coarseRhs, zeroGuess);
      CompCoarse->stop();
      emptyCoarseSolve = false;
    }
    if (Coarse->IsAvailable("PostSmoother")) {
      CompCoarse->start();
      postSmoo_coarse->Apply(*coarseX, *coarseRhs, zeroGuess);
      CompCoarse->stop();
      emptyCoarseSolve = false;
    }
    if (emptyCoarseSolve == true) {
      // Coarse operator is identity
      coarseX->update(one, *coarseRhs, zero);
    }
    Concurrent->stop();
    // Synchronize_end->start();
    // communicator->barrier();
    // Synchronize_end->stop();

    if (!doPRrebalance_ && !importer.is_null()) {
      RCP<TimeMonitor> ITime      = rcp(new TimeMonitor(*this, prefix + "Solve : export (total)", Timings0));
      RCP<TimeMonitor> ILevelTime = rcp(new TimeMonitor(*this, prefix + "Solve : export" + levelSuffix1, Timings0));

      // Import: range map of rebalanced Ac (before subcomm replacement) --> domain map of P
      RCP<MultiVector> coarseTmp = MultiVectorFactory::Build(importer->getSourceMap(), coarseX->getNumVectors());
      coarseTmp->doExport(*coarseX, *importer, Xpetra::INSERT);
      coarseX.swap(coarseTmp);
    }

    ApplyPbar->start();
    Pbar->apply(*coarseX, *coarseX_prolonged, Teuchos::NO_TRANS, one, zero);
    ApplyPbar->stop();
  }

  ApplySum->start();
  X.update(1.0, *fineX, 1.0, *coarseX_prolonged, 0.0);
  ApplySum->stop();

  CompTime->stop();

  // communicator->barrier();

  return (tol > 0 ? ConvergenceStatus::Unconverged : ConvergenceStatus::Undefined);
}
#else
// ---------------------------------------- Iterate -------------------------------------------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ConvergenceStatus Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Iterate(const MultiVector& B, MultiVector& X, ConvData conv,
                                                                                bool InitialGuessIsZero, LO startLevel) {
  LO nIts           = conv.maxIts_;
  MagnitudeType tol = conv.tol_;

  // These timers work as follows. "iterateTime" records total time spent in
  // iterate. "levelTime" records time on a per level basis. The label is
  // crafted to mimic the per-level messages used in Monitors. Note that a
  // given level is timed with a TimeMonitor instead of a Monitor or
  // SubMonitor. This is mainly because I want to time each level
  // separately, and Monitors/SubMonitors print "(total) xx yy zz" ,
  // "(sub,total) xx yy zz", respectively, which is subject to
  // misinterpretation.  The per-level TimeMonitors are stopped/started
  // manually before/after a recursive call to Iterate. A side artifact to
  // this approach is that the counts for intermediate level timers are twice
  // the counts for the finest and coarsest levels.

  RCP<Level> Fine = Levels_[startLevel];

  std::string label        = FormattingHelper::getColonLabel(Fine->getObjectLabel());
  std::string prefix       = label + this->ShortClassName() + ": ";
  std::string levelSuffix  = " (level=" + toString(startLevel) + ")";
  std::string levelSuffix1 = " (level=" + toString(startLevel + 1) + ")";

  bool useStackedTimer = !Teuchos::TimeMonitor::getStackedTimer().is_null();

  RCP<Monitor> iterateTime;
  RCP<TimeMonitor> iterateTime1;
  if (startLevel == 0)
    iterateTime = rcp(new Monitor(*this, "Solve", label, (nIts == 1) ? None : Runtime0, Timings0));
  else if (!useStackedTimer)
    iterateTime1 = rcp(new TimeMonitor(*this, prefix + "Solve (total, level=" + toString(startLevel) + ")", Timings0));

  std::string iterateLevelTimeLabel = prefix + "Solve" + levelSuffix;
  RCP<TimeMonitor> iterateLevelTime = rcp(new TimeMonitor(*this, iterateLevelTimeLabel, Timings0));

  bool zeroGuess = InitialGuessIsZero;

  RCP<Operator> A = Fine->Get<RCP<Operator> >("A");
  using namespace Teuchos;
  RCP<Time> CompCoarse = Teuchos::TimeMonitor::getNewCounter(prefix + "Coarse: Computational Time");

  if (A.is_null()) {
    // This processor does not have any data for this process on coarser
    // levels. This can only happen when there are multiple processors and
    // we use repartitioning.
    return ConvergenceStatus::Undefined;
  }

  // If we switched the number of vectors, we'd need to reallocate here.
  // If the number of vectors is unchanged, this is a noop.
  // NOTE: We need to check against B because the tests in AllocateLevelMultiVectors
  // will fail on Stokhos Scalar types (due to the so-called 'hidden dimension')
  const BlockedMultiVector* Bblocked = dynamic_cast<const BlockedMultiVector*>(&B);
  if (residual_.size() > startLevel &&
      ((Bblocked && !Bblocked->isSameSize(*residual_[startLevel])) ||
       (!Bblocked && !residual_[startLevel]->isSameSize(B))))
    DeleteLevelMultiVectors();
  AllocateLevelMultiVectors(X.getNumVectors());

  // Print residual information before iterating
  typedef Teuchos::ScalarTraits<typename STS::magnitudeType> STM;
  MagnitudeType prevNorm = STM::one();
  rate_                  = 1.0;
  if (IsCalculationOfResidualRequired(startLevel, conv)) {
    ConvergenceStatus convergenceStatus = ComputeResidualAndPrintHistory(*A, X, B, Teuchos::ScalarTraits<LO>::zero(), startLevel, conv, prevNorm);
    if (convergenceStatus == MueLu::ConvergenceStatus::Converged)
      return convergenceStatus;
  }

  SC one = STS::one(), zero = STS::zero();
  for (LO iteration = 1; iteration <= nIts; iteration++) {
#ifdef HAVE_MUELU_DEBUG
#if 0  // TODO fix me
      if (A->getDomainMap()->isCompatible(*(X.getMap())) == false) {
        std::ostringstream ss;
        ss << "Level " << startLevel << ": level A's domain map is not compatible with X";
        throw Exceptions::Incompatible(ss.str());
      }

      if (A->getRangeMap()->isCompatible(*(B.getMap())) == false) {
        std::ostringstream ss;
        ss << "Level " << startLevel << ": level A's range map is not compatible with B";
        throw Exceptions::Incompatible(ss.str());
      }
#endif
#endif

    if (startLevel == as<LO>(Levels_.size()) - 1) {
      // On the coarsest level, we do either smoothing (if defined) or a direct solve.
      RCP<TimeMonitor> CLevelTime = rcp(new TimeMonitor(*this, prefix + "Solve : coarse" + levelSuffix, Timings0));

      bool emptySolve = true;

      // NOTE: we need to check using IsAvailable before Get here to avoid building default smoother
      if (Fine->IsAvailable("PreSmoother")) {
        RCP<SmootherBase> preSmoo = Fine->Get<RCP<SmootherBase> >("PreSmoother");
        CompCoarse->start();
        preSmoo->Apply(X, B, zeroGuess);
        CompCoarse->stop();
        zeroGuess  = false;
        emptySolve = false;
      }
      if (Fine->IsAvailable("PostSmoother")) {
        RCP<SmootherBase> postSmoo = Fine->Get<RCP<SmootherBase> >("PostSmoother");
        CompCoarse->start();
        postSmoo->Apply(X, B, zeroGuess);
        CompCoarse->stop();
        emptySolve = false;
        zeroGuess  = false;
      }
      if (emptySolve == true) {
        // Coarse operator is identity
        X.update(one, B, zero);
      }

    } else {
      // On intermediate levels, we do cycles
      RCP<Level> Coarse = Levels_[startLevel + 1];
      {
        // ============== PRESMOOTHING ==============
        RCP<TimeMonitor> STime;
        if (!useStackedTimer)
          STime = rcp(new TimeMonitor(*this, prefix + "Solve : smoothing (total)", Timings0));
        RCP<TimeMonitor> SLevelTime = rcp(new TimeMonitor(*this, prefix + "Solve : smoothing" + levelSuffix, Timings0));

        if (Fine->IsAvailable("PreSmoother")) {
          RCP<SmootherBase> preSmoo = Fine->Get<RCP<SmootherBase> >("PreSmoother");
          preSmoo->Apply(X, B, zeroGuess);
          zeroGuess = false;
        }
      }

      RCP<MultiVector> residual;
      {
        RCP<TimeMonitor> ATime;
        if (!useStackedTimer)
          ATime = rcp(new TimeMonitor(*this, prefix + "Solve : residual calculation (total)", Timings0));
        RCP<TimeMonitor> ALevelTime = rcp(new TimeMonitor(*this, prefix + "Solve : residual calculation" + levelSuffix, Timings0));
        if (zeroGuess) {
          // If there's a pre-smoother, then zeroGuess is false.  If there isn't and people still have zeroGuess set,
          // then then X still has who knows what, so we need to zero it before we go to the coarse grid.
          X.putScalar(zero);
        }

        Utilities::Residual(*A, X, B, *residual_[startLevel]);
        residual = residual_[startLevel];
      }

      RCP<Operator> P = Coarse->Get<RCP<Operator> >("P");
      if (Coarse->IsAvailable("Pbar"))
        P = Coarse->Get<RCP<Operator> >("Pbar");

      RCP<MultiVector> coarseRhs, coarseX;
      //        const bool initializeWithZeros = true;
      {
        // ============== RESTRICTION ==============
        RCP<TimeMonitor> RTime;
        if (!useStackedTimer)
          RTime = rcp(new TimeMonitor(*this, prefix + "Solve : restriction (total)", Timings0));
        RCP<TimeMonitor> RLevelTime = rcp(new TimeMonitor(*this, prefix + "Solve : restriction" + levelSuffix, Timings0));
        coarseRhs                   = coarseRhs_[startLevel];

        if (implicitTranspose_) {
          P->apply(*residual, *coarseRhs, Teuchos::TRANS, one, zero);

        } else {
          RCP<Operator> R = Coarse->Get<RCP<Operator> >("R");
          R->apply(*residual, *coarseRhs, Teuchos::NO_TRANS, one, zero);
        }
      }

      RCP<const Import> importer;
      if (Coarse->IsAvailable("Importer"))
        importer = Coarse->Get<RCP<const Import> >("Importer");

      coarseX = coarseX_[startLevel];
      if (!doPRrebalance_ && !importer.is_null()) {
        RCP<TimeMonitor> ITime;
        if (!useStackedTimer)
          ITime = rcp(new TimeMonitor(*this, prefix + "Solve : import (total)", Timings0));
        RCP<TimeMonitor> ILevelTime = rcp(new TimeMonitor(*this, prefix + "Solve : import" + levelSuffix1, Timings0));

        // Import: range map of R --> domain map of rebalanced Ac (before subcomm replacement)
        RCP<MultiVector> coarseTmp = coarseImport_[startLevel];
        coarseTmp->doImport(*coarseRhs, *importer, Xpetra::INSERT);
        coarseRhs.swap(coarseTmp);
      }

      RCP<Operator> Ac = Coarse->Get<RCP<Operator> >("A");
      if (!Ac.is_null()) {
        RCP<const Map> origXMap   = coarseX->getMap();
        RCP<const Map> origRhsMap = coarseRhs->getMap();

        // Replace maps with maps with a subcommunicator
        coarseRhs->replaceMap(Ac->getRangeMap());
        coarseX->replaceMap(Ac->getDomainMap());

        {
          iterateLevelTime = Teuchos::null;  // stop timing this level

          Iterate(*coarseRhs, *coarseX, 1, true, startLevel + 1);
          // ^^ zero initial guess
          if (Cycle_ == WCYCLE && WCycleStartLevel_ >= startLevel)
            Iterate(*coarseRhs, *coarseX, 1, false, startLevel + 1);
          // ^^ nonzero initial guess

          iterateLevelTime = rcp(new TimeMonitor(*this, iterateLevelTimeLabel));  // restart timing this level
        }
        coarseX->replaceMap(origXMap);
        coarseRhs->replaceMap(origRhsMap);
      }

      if (!doPRrebalance_ && !importer.is_null()) {
        RCP<TimeMonitor> ITime;
        if (!useStackedTimer)
          ITime = rcp(new TimeMonitor(*this, prefix + "Solve : export (total)", Timings0));
        RCP<TimeMonitor> ILevelTime = rcp(new TimeMonitor(*this, prefix + "Solve : export" + levelSuffix1, Timings0));

        // Import: range map of rebalanced Ac (before subcomm replacement) --> domain map of P
        RCP<MultiVector> coarseTmp = coarseExport_[startLevel];
        coarseTmp->doExport(*coarseX, *importer, Xpetra::INSERT);
        coarseX.swap(coarseTmp);
      }

      {
        // ============== PROLONGATION ==============
        RCP<TimeMonitor> PTime;
        if (!useStackedTimer)
          PTime = rcp(new TimeMonitor(*this, prefix + "Solve : prolongation (total)", Timings0));
        RCP<TimeMonitor> PLevelTime = rcp(new TimeMonitor(*this, prefix + "Solve : prolongation" + levelSuffix, Timings0));
        // Update X += P * coarseX
        // Note that due to what may be round-off error accumulation, use of the fused kernel
        //    P->apply(*coarseX, X, Teuchos::NO_TRANS, one, one);
        // can in some cases result in slightly higher iteration counts.
        if (fuseProlongationAndUpdate_) {
          P->apply(*coarseX, X, Teuchos::NO_TRANS, scalingFactor_, one);
        } else {
          RCP<MultiVector> correction = correction_[startLevel];
          P->apply(*coarseX, *correction, Teuchos::NO_TRANS, one, zero);
          X.update(scalingFactor_, *correction, one);
        }
      }

      {
        // ============== POSTSMOOTHING ==============
        RCP<TimeMonitor> STime;
        if (!useStackedTimer)
          STime = rcp(new TimeMonitor(*this, prefix + "Solve : smoothing (total)", Timings0));
        RCP<TimeMonitor> SLevelTime = rcp(new TimeMonitor(*this, prefix + "Solve : smoothing" + levelSuffix, Timings0));

        if (Fine->IsAvailable("PostSmoother")) {
          RCP<SmootherBase> postSmoo = Fine->Get<RCP<SmootherBase> >("PostSmoother");
          postSmoo->Apply(X, B, false);
        }
      }
    }
    zeroGuess = false;

    if (IsCalculationOfResidualRequired(startLevel, conv)) {
      ConvergenceStatus convergenceStatus = ComputeResidualAndPrintHistory(*A, X, B, iteration, startLevel, conv, prevNorm);
      if (convergenceStatus == MueLu::ConvergenceStatus::Converged)
        return convergenceStatus;
    }
  }
  return (tol > 0 ? ConvergenceStatus::Unconverged : ConvergenceStatus::Undefined);
}
#endif

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write(const LO& start, const LO& end, const std::string& suffix) {
  LO startLevel = (start != -1 ? start : 0);
  LO endLevel   = (end != -1 ? end : Levels_.size() - 1);

  TEUCHOS_TEST_FOR_EXCEPTION(startLevel > endLevel, Exceptions::RuntimeError,
                             "MueLu::Hierarchy::Write : startLevel must be <= endLevel");

  TEUCHOS_TEST_FOR_EXCEPTION(startLevel < 0 || endLevel >= Levels_.size(), Exceptions::RuntimeError,
                             "MueLu::Hierarchy::Write bad start or end level");

  for (LO i = startLevel; i < endLevel + 1; i++) {
    RCP<Matrix> A = rcp_dynamic_cast<Matrix>(Levels_[i]->template Get<RCP<Operator> >("A")), P, R;
    if (i > 0) {
      P = rcp_dynamic_cast<Matrix>(Levels_[i]->template Get<RCP<Operator> >("P"));
      if (!implicitTranspose_)
        R = rcp_dynamic_cast<Matrix>(Levels_[i]->template Get<RCP<Operator> >("R"));
    }

    if (!A.is_null()) Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("A_" + toString(i) + suffix + ".m", *A);
    if (!P.is_null()) {
      Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("P_" + toString(i) + suffix + ".m", *P);
    }
    if (!R.is_null()) {
      Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("R_" + toString(i) + suffix + ".m", *R);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Keep(const std::string& ename, const FactoryBase* factory) {
  for (Array<RCP<Level> >::iterator it = Levels_.begin(); it != Levels_.end(); ++it)
    (*it)->Keep(ename, factory);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Delete(const std::string& ename, const FactoryBase* factory) {
  for (Array<RCP<Level> >::iterator it = Levels_.begin(); it != Levels_.end(); ++it)
    (*it)->Delete(ename, factory);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddKeepFlag(const std::string& ename, const FactoryBase* factory, KeepType keep) {
  for (Array<RCP<Level> >::iterator it = Levels_.begin(); it != Levels_.end(); ++it)
    (*it)->AddKeepFlag(ename, factory, keep);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RemoveKeepFlag(const std::string& ename, const FactoryBase* factory, KeepType keep) {
  for (Array<RCP<Level> >::iterator it = Levels_.begin(); it != Levels_.end(); ++it)
    (*it)->RemoveKeepFlag(ename, factory, keep);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  if (description_.empty()) {
    std::ostringstream out;
    out << BaseClass::description();
    out << "{#levels = " << GetGlobalNumLevels() << ", complexity = " << GetOperatorComplexity() << "}";
    description_ = out.str();
  }
  return description_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel tVerbLevel) const {
  describe(out, toMueLuVerbLevel(tVerbLevel));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::describe(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
  RCP<Operator> A0                    = Levels_[0]->template Get<RCP<Operator> >("A");
  RCP<const Teuchos::Comm<int> > comm = A0->getDomainMap()->getComm();

  int numLevels    = GetNumLevels();
  RCP<Operator> Ac = Levels_[numLevels - 1]->template Get<RCP<Operator> >("A");
  if (Ac.is_null()) {
    // It may happen that we do repartition on the last level, but the matrix
    // is small enough to satisfy "max coarse size" requirement. Then, even
    // though we have the level, the matrix would be null on all but one processors
    numLevels--;
  }
  int root = comm->getRank();

#ifdef HAVE_MPI
  int smartData = numLevels * comm->getSize() + comm->getRank(), maxSmartData;
  reduceAll(*comm, Teuchos::REDUCE_MAX, smartData, Teuchos::ptr(&maxSmartData));
  root = maxSmartData % comm->getSize();
#endif

  // Compute smoother complexity, if needed
  double smoother_comp = -1.0;
  if (verbLevel & (Statistics0 | Test))
    smoother_comp = GetSmootherComplexity();

  std::string outstr;
  if (comm->getRank() == root && verbLevel & (Statistics0 | Test)) {
    std::vector<Xpetra::global_size_t> nnzPerLevel;
    std::vector<Xpetra::global_size_t> rowsPerLevel;
    std::vector<int> numProcsPerLevel;
    bool someOpsNotMatrices             = false;
    const Xpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    for (int i = 0; i < numLevels; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION(!(Levels_[i]->IsAvailable("A")), Exceptions::RuntimeError,
                                 "Operator A is not available on level " << i);

      RCP<Operator> A = Levels_[i]->template Get<RCP<Operator> >("A");
      TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(), Exceptions::RuntimeError,
                                 "Operator A on level " << i << " is null.");

      RCP<Matrix> Am = rcp_dynamic_cast<Matrix>(A);
      if (Am.is_null()) {
        someOpsNotMatrices = true;
        nnzPerLevel.push_back(INVALID);
        rowsPerLevel.push_back(A->getDomainMap()->getGlobalNumElements());
        numProcsPerLevel.push_back(A->getDomainMap()->getComm()->getSize());
      } else {
        LO storageblocksize       = Am->GetStorageBlockSize();
        Xpetra::global_size_t nnz = Am->getGlobalNumEntries() * storageblocksize * storageblocksize;
        nnzPerLevel.push_back(nnz);
        rowsPerLevel.push_back(Am->getGlobalNumRows() * storageblocksize);
        numProcsPerLevel.push_back(Am->getRowMap()->getComm()->getSize());
      }
    }
    if (someOpsNotMatrices)
      GetOStream(Warnings0) << "Some level operators are not matrices, statistics calculation are incomplete" << std::endl;

    {
      std::string label = Levels_[0]->getObjectLabel();
      std::ostringstream oss;
      oss << std::setfill(' ');
      oss << "\n--------------------------------------------------------------------------------\n";
      oss << "---                            Multigrid Summary " << std::setw(28) << std::left << label << "---\n";
      oss << "--------------------------------------------------------------------------------" << std::endl;
      if (verbLevel & Parameters1)
        oss << "Scalar              = " << Teuchos::ScalarTraits<Scalar>::name() << std::endl;
      oss << "Number of levels    = " << numLevels << std::endl;
      oss << "Operator complexity = " << std::setprecision(2) << std::setiosflags(std::ios::fixed);
      if (!someOpsNotMatrices)
        oss << GetOperatorComplexity() << std::endl;
      else
        oss << "not available (Some operators in hierarchy are not matrices.)" << std::endl;

      if (smoother_comp != -1.0) {
        oss << "Smoother complexity = " << std::setprecision(2) << std::setiosflags(std::ios::fixed)
            << smoother_comp << std::endl;
      }

      switch (Cycle_) {
        case VCYCLE:
          oss << "Cycle type          = V" << std::endl;
          break;
        case WCYCLE:
          oss << "Cycle type          = W" << std::endl;
          if (WCycleStartLevel_ > 0)
            oss << "Cycle start level   = " << WCycleStartLevel_ << std::endl;
          break;
        default:
          break;
      };
      oss << std::endl;

      Xpetra::global_size_t tt = rowsPerLevel[0];
      int rowspacer            = 2;
      while (tt != 0) {
        tt /= 10;
        rowspacer++;
      }
      for (size_t i = 0; i < nnzPerLevel.size(); ++i) {
        tt = nnzPerLevel[i];
        if (tt != INVALID)
          break;
        tt = 100;  // This will get used if all levels are operators.
      }
      int nnzspacer = 2;
      while (tt != 0) {
        tt /= 10;
        nnzspacer++;
      }
      tt           = numProcsPerLevel[0];
      int npspacer = 2;
      while (tt != 0) {
        tt /= 10;
        npspacer++;
      }
      oss << "level " << std::setw(rowspacer) << " rows " << std::setw(nnzspacer) << " nnz "
          << " nnz/row" << std::setw(npspacer) << "  c ratio"
          << "  procs" << std::endl;
      for (size_t i = 0; i < nnzPerLevel.size(); ++i) {
        oss << "  " << i << "  ";
        oss << std::setw(rowspacer) << rowsPerLevel[i];
        if (nnzPerLevel[i] != INVALID) {
          oss << std::setw(nnzspacer) << nnzPerLevel[i];
          oss << std::setprecision(2) << std::setiosflags(std::ios::fixed);
          oss << std::setw(9) << as<double>(nnzPerLevel[i]) / rowsPerLevel[i];
        } else {
          oss << std::setw(nnzspacer) << "Operator";
          oss << std::setprecision(2) << std::setiosflags(std::ios::fixed);
          oss << std::setw(9) << "     ";
        }
        if (i)
          oss << std::setw(9) << as<double>(rowsPerLevel[i - 1]) / rowsPerLevel[i];
        else
          oss << std::setw(9) << "     ";
        oss << "    " << std::setw(npspacer) << numProcsPerLevel[i] << std::endl;
      }
      oss << std::endl;
      for (int i = 0; i < GetNumLevels(); ++i) {
        RCP<SmootherBase> preSmoo, postSmoo;
        if (Levels_[i]->IsAvailable("PreSmoother"))
          preSmoo = Levels_[i]->template Get<RCP<SmootherBase> >("PreSmoother");
        if (Levels_[i]->IsAvailable("PostSmoother"))
          postSmoo = Levels_[i]->template Get<RCP<SmootherBase> >("PostSmoother");

        if (preSmoo != null && preSmoo == postSmoo)
          oss << "Smoother (level " << i << ") both : " << preSmoo->description() << std::endl;
        else {
          oss << "Smoother (level " << i << ") pre  : "
              << (preSmoo != null ? preSmoo->description() : "no smoother") << std::endl;
          oss << "Smoother (level " << i << ") post : "
              << (postSmoo != null ? postSmoo->description() : "no smoother") << std::endl;
        }

        oss << std::endl;
      }

      outstr = oss.str();
    }
  }

#ifdef HAVE_MPI
  RCP<const Teuchos::MpiComm<int> > mpiComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
  MPI_Comm rawComm                          = (*mpiComm->getRawMpiComm())();

  int strLength = outstr.size();
  MPI_Bcast(&strLength, 1, MPI_INT, root, rawComm);
  if (comm->getRank() != root)
    outstr.resize(strLength);
  MPI_Bcast(&outstr[0], strLength, MPI_CHAR, root, rawComm);
#endif

  out << outstr;
}

// NOTE: at some point this should be replaced by a friend operator <<
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(std::ostream& out, const VerbLevel verbLevel) const {
  Teuchos::OSTab tab2(out);
  for (int i = 0; i < GetNumLevels(); ++i)
    Levels_[i]->print(out, verbLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::IsPreconditioner(const bool flag) {
  isPreconditioner_ = flag;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DumpCurrentGraph(int currLevel) const {
  if (GetProcRankVerbose() != 0)
    return;
#if defined(HAVE_MUELU_BOOST) && defined(HAVE_MUELU_BOOST_FOR_REAL) && defined(BOOST_VERSION) && (BOOST_VERSION >= 104400)

  BoostGraph graph;

  BoostProperties dp;
  dp.property("label", boost::get(boost::vertex_name, graph));
  dp.property("id", boost::get(boost::vertex_index, graph));
  dp.property("label", boost::get(boost::edge_name, graph));
  dp.property("color", boost::get(boost::edge_color, graph));

  // create local maps
  std::map<const FactoryBase*, BoostVertex> vindices;
  typedef std::map<std::pair<BoostVertex, BoostVertex>, std::string> emap;
  emap edges;

  static int call_id = 0;

  RCP<Operator> A = Levels_[0]->template Get<RCP<Operator> >("A");
  int rank        = A->getDomainMap()->getComm()->getRank();

  //    printf("[%d] CMS: ----------------------\n",rank);
  for (int i = currLevel; i <= currLevel + 1 && i < GetNumLevels(); i++) {
    edges.clear();
    Levels_[i]->UpdateGraph(vindices, edges, dp, graph);

    for (emap::const_iterator eit = edges.begin(); eit != edges.end(); eit++) {
      std::pair<BoostEdge, bool> boost_edge = boost::add_edge(eit->first.first, eit->first.second, graph);
      // printf("[%d] CMS:   Hierarchy, adding edge (%d->%d) %d\n",rank,(int)eit->first.first,(int)eit->first.second,(int)boost_edge.second);
      // Because xdot.py views 'Graph' as a keyword
      if (eit->second == std::string("Graph"))
        boost::put("label", dp, boost_edge.first, std::string("Graph_"));
      else
        boost::put("label", dp, boost_edge.first, eit->second);
      if (i == currLevel)
        boost::put("color", dp, boost_edge.first, std::string("red"));
      else
        boost::put("color", dp, boost_edge.first, std::string("blue"));
    }
  }

  std::ofstream out(dumpFile_.c_str() + std::string("_") + std::to_string(currLevel) + std::string("_") + std::to_string(call_id) + std::string("_") + std::to_string(rank) + std::string(".dot"));
  boost::write_graphviz_dp(out, graph, dp, std::string("id"));
  out.close();
  call_id++;
#else
  GetOStream(Errors) << "Dependency graph output requires boost and MueLu_ENABLE_Boost_for_real" << std::endl;
#endif
}

// Enforce that coordinate vector's map is consistent with that of A
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReplaceCoordinateMap(Level& level) {
  RCP<Operator> Ao = level.Get<RCP<Operator> >("A");
  RCP<Matrix> A    = rcp_dynamic_cast<Matrix>(Ao);
  if (A.is_null()) {
    GetOStream(Runtime1) << "Hierarchy::ReplaceCoordinateMap: operator is not a matrix, skipping..." << std::endl;
    return;
  }
  if (Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A) != Teuchos::null) {
    GetOStream(Runtime1) << "Hierarchy::ReplaceCoordinateMap: operator is a BlockedCrsMatrix, skipping..." << std::endl;
    return;
  }

  typedef Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO, NO> xdMV;

  RCP<xdMV> coords = level.Get<RCP<xdMV> >("Coordinates");

  if (A->getRowMap()->isSameAs(*(coords->getMap()))) {
    GetOStream(Runtime1) << "Hierarchy::ReplaceCoordinateMap: matrix and coordinates maps are same, skipping..." << std::endl;
    return;
  }

  if (A->IsView("stridedMaps") && rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
    RCP<const StridedMap> stridedRowMap = rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps"));

    // It is better to through an exceptions if maps may be inconsistent, than to ignore it and experience unfathomable breakdowns
    TEUCHOS_TEST_FOR_EXCEPTION(stridedRowMap->getStridedBlockId() != -1 || stridedRowMap->getOffset() != 0,
                               Exceptions::RuntimeError, "Hierarchy::ReplaceCoordinateMap: nontrivial maps (block id = " << stridedRowMap->getStridedBlockId() << ", offset = " << stridedRowMap->getOffset() << ")");
  }

  GetOStream(Runtime1) << "Replacing coordinate map" << std::endl;
  TEUCHOS_TEST_FOR_EXCEPTION(A->GetFixedBlockSize() % A->GetStorageBlockSize() != 0, Exceptions::RuntimeError, "Hierarchy::ReplaceCoordinateMap: Storage block size does not evenly divide fixed block size");

  size_t blkSize = A->GetFixedBlockSize() / A->GetStorageBlockSize();

  RCP<const Map> nodeMap = A->getRowMap();
  if (blkSize > 1) {
    // Create a nodal map, as coordinates have not been expanded to a DOF map yet.
    RCP<const Map> dofMap = A->getRowMap();
    GO indexBase          = dofMap->getIndexBase();
    size_t numLocalDOFs   = dofMap->getLocalNumElements();
    TEUCHOS_TEST_FOR_EXCEPTION(numLocalDOFs % blkSize, Exceptions::RuntimeError,
                               "Hierarchy::ReplaceCoordinateMap: block size (" << blkSize << ") is incompatible with the number of local dofs in a row map (" << numLocalDOFs);
    ArrayView<const GO> GIDs = dofMap->getLocalElementList();

    Array<GO> nodeGIDs(numLocalDOFs / blkSize);
    for (size_t i = 0; i < numLocalDOFs; i += blkSize)
      nodeGIDs[i / blkSize] = (GIDs[i] - indexBase) / blkSize + indexBase;

    Xpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    nodeMap                       = MapFactory::Build(dofMap->lib(), INVALID, nodeGIDs(), indexBase, dofMap->getComm());
  } else {
    // blkSize == 1
    // Check whether the length of vectors fits to the size of A
    // If yes, make sure that the maps are matching
    // If no, throw a warning but do not touch the Coordinates
    if (coords->getLocalLength() != A->getRowMap()->getLocalNumElements()) {
      GetOStream(Warnings) << "Coordinate vector does not match row map of matrix A!" << std::endl;
      return;
    }
  }

  Array<ArrayView<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> > coordDataView;
  std::vector<ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> > coordData;
  for (size_t i = 0; i < coords->getNumVectors(); i++) {
    coordData.push_back(coords->getData(i));
    coordDataView.push_back(coordData[i]());
  }

  RCP<xdMV> newCoords = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO, NO>::Build(nodeMap, coordDataView(), coords->getNumVectors());
  level.Set("Coordinates", newCoords);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AllocateLevelMultiVectors(int numvecs, bool forceMapCheck) {
  int N = Levels_.size();
  if (((sizeOfAllocatedLevelMultiVectors_ == numvecs && residual_.size() == N) || numvecs <= 0) && !forceMapCheck) return;

  // If, somehow, we changed the number of levels, delete everything first
  if (residual_.size() != N) {
    DeleteLevelMultiVectors();

    residual_.resize(N);
    coarseRhs_.resize(N);
    coarseX_.resize(N);
    coarseImport_.resize(N);
    coarseExport_.resize(N);
    correction_.resize(N);
  }

  for (int i = 0; i < N; i++) {
    RCP<Operator> A = Levels_[i]->template Get<RCP<Operator> >("A");
    if (!A.is_null()) {
      // This dance is because we allow A to have a BlockedMap and X/B to have (compatible) non-blocked map
      RCP<const BlockedCrsMatrix> A_as_blocked = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(A);
      RCP<const Map> Arm                       = A->getRangeMap();
      RCP<const Map> Adm                       = A->getDomainMap();
      if (!A_as_blocked.is_null()) {
        Adm = A_as_blocked->getFullDomainMap();
      }

      if (residual_[i].is_null() || !residual_[i]->getMap()->isSameAs(*Arm))
        // This is zero'd by default since it is filled via an operator apply
        residual_[i] = MultiVectorFactory::Build(Arm, numvecs, true);
      if (correction_[i].is_null() || !correction_[i]->getMap()->isSameAs(*Adm))
        correction_[i] = MultiVectorFactory::Build(Adm, numvecs, false);
    }

    if (i + 1 < N) {
      // This is zero'd by default since it is filled via an operator apply
      if (implicitTranspose_) {
        RCP<Operator> P = Levels_[i + 1]->template Get<RCP<Operator> >("P");
        if (!P.is_null()) {
          RCP<const Map> map = P->getDomainMap();
          if (coarseRhs_[i].is_null() || !coarseRhs_[i]->getMap()->isSameAs(*map))
            coarseRhs_[i] = MultiVectorFactory::Build(map, numvecs, true);
        }
      } else {
        RCP<Operator> R = Levels_[i + 1]->template Get<RCP<Operator> >("R");
        if (!R.is_null()) {
          RCP<const Map> map = R->getRangeMap();
          if (coarseRhs_[i].is_null() || !coarseRhs_[i]->getMap()->isSameAs(*map))
            coarseRhs_[i] = MultiVectorFactory::Build(map, numvecs, true);
        }
      }

      RCP<const Import> importer;
      if (Levels_[i + 1]->IsAvailable("Importer"))
        importer = Levels_[i + 1]->template Get<RCP<const Import> >("Importer");
      if (doPRrebalance_ || importer.is_null()) {
        RCP<const Map> map = coarseRhs_[i]->getMap();
        if (coarseX_[i].is_null() || !coarseX_[i]->getMap()->isSameAs(*map))
          coarseX_[i] = MultiVectorFactory::Build(map, numvecs, true);
      } else {
        RCP<const Map> map;
        map = importer->getTargetMap();
        if (coarseImport_[i].is_null() || !coarseImport_[i]->getMap()->isSameAs(*map)) {
          coarseImport_[i] = MultiVectorFactory::Build(map, numvecs, false);
          coarseX_[i]      = MultiVectorFactory::Build(map, numvecs, false);
        }
        map = importer->getSourceMap();
        if (coarseExport_[i].is_null() || !coarseExport_[i]->getMap()->isSameAs(*map))
          coarseExport_[i] = MultiVectorFactory::Build(map, numvecs, false);
      }
    }
  }
  sizeOfAllocatedLevelMultiVectors_ = numvecs;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeleteLevelMultiVectors() {
  if (sizeOfAllocatedLevelMultiVectors_ == 0) return;
  residual_.resize(0);
  coarseRhs_.resize(0);
  coarseX_.resize(0);
  coarseImport_.resize(0);
  coarseExport_.resize(0);
  correction_.resize(0);
  sizeOfAllocatedLevelMultiVectors_ = 0;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::IsCalculationOfResidualRequired(
    const LO startLevel, const ConvData& conv) const {
  return (startLevel == 0 && !isPreconditioner_ && (IsPrint(Statistics1) || conv.tol_ > 0));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ConvergenceStatus Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::IsConverged(
    const Teuchos::Array<MagnitudeType>& residualNorm, const MagnitudeType convergenceTolerance) const {
  ConvergenceStatus convergenceStatus = ConvergenceStatus::Undefined;

  if (convergenceTolerance > Teuchos::ScalarTraits<MagnitudeType>::zero()) {
    bool passed = true;
    for (LO k = 0; k < residualNorm.size(); k++)
      if (residualNorm[k] >= convergenceTolerance)
        passed = false;

    if (passed)
      convergenceStatus = ConvergenceStatus::Converged;
    else
      convergenceStatus = ConvergenceStatus::Unconverged;
  }

  return convergenceStatus;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PrintResidualHistory(
    const LO iteration, const Teuchos::Array<MagnitudeType>& residualNorm) const {
  GetOStream(Statistics1) << "iter:    "
                          << std::setiosflags(std::ios::left)
                          << std::setprecision(3) << std::setw(4) << iteration
                          << "           residual = "
                          << std::setprecision(10) << residualNorm
                          << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ConvergenceStatus Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeResidualAndPrintHistory(
    const Operator& A, const MultiVector& X, const MultiVector& B, const LO iteration,
    const LO startLevel, const ConvData& conv, MagnitudeType& previousResidualNorm) {
  Teuchos::Array<MagnitudeType> residualNorm;
  residualNorm = Utilities::ResidualNorm(A, X, B, *residual_[startLevel]);

  const MagnitudeType currentResidualNorm = residualNorm[0];
  rate_                                   = currentResidualNorm / previousResidualNorm;
  previousResidualNorm                    = currentResidualNorm;

  if (IsPrint(Statistics1))
    PrintResidualHistory(iteration, residualNorm);

  return IsConverged(residualNorm, conv.tol_);
}

}  // namespace MueLu

#endif  // MUELU_HIERARCHY_DEF_HPP
