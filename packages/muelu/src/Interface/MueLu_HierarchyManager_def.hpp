// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_HIERARCHYMANAGER_DEF_HPP
#define MUELU_HIERARCHYMANAGER_DEF_HPP

#include <string>
#include <map>

#include <Teuchos_Array.hpp>

#include <Xpetra_Operator.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_HierarchyManager_decl.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_HierarchyFactory.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_PerfUtils.hpp"

#ifdef HAVE_MUELU_INTREPID2
#include "Kokkos_DynRankView.hpp"
#endif

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::HierarchyManager(int numDesiredLevel)
  : numDesiredLevel_(numDesiredLevel)
  , maxCoarseSize_(MasterList::getDefault<int>("coarse: max size"))
  , verbosity_(Medium)
  , doPRrebalance_(MasterList::getDefault<bool>("repartition: rebalance P and R"))
  , doPRViaCopyrebalance_(MasterList::getDefault<bool>("repartition: explicit via new copy rebalance P and R"))
  , implicitTranspose_(MasterList::getDefault<bool>("transpose: use implicit"))
  , fuseProlongationAndUpdate_(MasterList::getDefault<bool>("fuse prolongation and update"))
  , suppressNullspaceDimensionCheck_(MasterList::getDefault<bool>("nullspace: suppress dimension check"))
  , sizeOfMultiVectors_(MasterList::getDefault<int>("number of vectors"))
  , graphOutputLevel_(-2) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddFactoryManager(int startLevel, int numDesiredLevel, RCP<FactoryManagerBase> manager) {
  const int lastLevel = startLevel + numDesiredLevel - 1;
  if (levelManagers_.size() < lastLevel + 1)
    levelManagers_.resize(lastLevel + 1);

  for (int iLevel = startLevel; iLevel <= lastLevel; iLevel++)
    levelManagers_[iLevel] = manager;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<FactoryManagerBase> MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetFactoryManager(int levelID) const {
  // NOTE: last levelManager is used for all the remaining levels
  return (levelID >= levelManagers_.size() ? levelManagers_[levelManagers_.size() - 1] : levelManagers_[levelID]);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNumFactoryManagers() const {
  return levelManagers_.size();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CheckConfig() {
  for (int i = 0; i < levelManagers_.size(); i++)
    TEUCHOS_TEST_FOR_EXCEPTION(levelManagers_[i] == Teuchos::null, Exceptions::RuntimeError, "MueLu:HierarchyConfig::CheckConfig(): Undefined configuration for level:");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateHierarchy() const {
  return rcp(new Hierarchy());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateHierarchy(const std::string& label) const {
  return rcp(new Hierarchy(label));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupHierarchy(Hierarchy& H) const {
  TEUCHOS_TEST_FOR_EXCEPTION(!H.GetLevel(0)->IsAvailable("A"), Exceptions::RuntimeError, "No fine level operator");

  RCP<Level> l0    = H.GetLevel(0);
  RCP<Operator> Op = l0->Get<RCP<Operator>>("A");

  // Compare nullspace dimension to NumPDEs and throw/warn based on user input
  if (l0->IsAvailable("Nullspace")) {
    RCP<Matrix> A = Teuchos::rcp_dynamic_cast<Matrix>(Op);
    if (A != Teuchos::null) {
      RCP<MultiVector> nullspace = l0->Get<RCP<MultiVector>>("Nullspace");

      if (static_cast<size_t>(A->GetFixedBlockSize()) > nullspace->getNumVectors()) {
        std::stringstream msg;
        msg << "User-provided nullspace has fewer vectors ("
            << nullspace->getNumVectors() << ") than number of PDE equations ("
            << A->GetFixedBlockSize() << "). ";

        if (suppressNullspaceDimensionCheck_) {
          msg << "It depends on the PDE, if this is a problem or not.";
          this->GetOStream(Warnings0) << msg.str() << std::endl;
        } else {
          msg << "Add the missing nullspace vectors! (You can suppress this check. See the MueLu user guide for details.)";
          TEUCHOS_TEST_FOR_EXCEPTION(static_cast<size_t>(A->GetFixedBlockSize()) > nullspace->getNumVectors(), Exceptions::RuntimeError, msg.str());
        }
      }
    } else {
      this->GetOStream(Warnings0) << "Skipping dimension check of user-supplied nullspace because user-supplied operator is not a matrix" << std::endl;
    }
  }

#ifdef HAVE_MUELU_DEBUG
  // Reset factories' data used for debugging
  for (int i = 0; i < levelManagers_.size(); i++)
    levelManagers_[i]->ResetDebugData();

#endif

  // Setup Matrix
  // TODO: I should certainly undo this somewhere...

  Xpetra::UnderlyingLib lib = Op->getDomainMap()->lib();
  H.setlib(lib);

  SetupOperator(*Op);
  SetupExtra(H);

  // Setup Hierarchy
  H.SetMaxCoarseSize(maxCoarseSize_);
  VerboseObject::SetDefaultVerbLevel(verbosity_);
  if (graphOutputLevel_ >= 0 || graphOutputLevel_ == -1)
    H.EnableGraphDumping("dep_graph", graphOutputLevel_);

  if (VerboseObject::IsPrint(Statistics2)) {
    RCP<Matrix> Amat = rcp_dynamic_cast<Matrix>(Op);

    if (!Amat.is_null()) {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo", true);

      VerboseObject::GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*Amat, "A0", params);
    } else {
      VerboseObject::GetOStream(Warnings1) << "Fine level operator is not a matrix, statistics are not available" << std::endl;
    }
  }

  H.SetPRrebalance(doPRrebalance_);
  H.SetPRViaCopyrebalance(doPRViaCopyrebalance_);
  H.SetImplicitTranspose(implicitTranspose_);
  H.SetFuseProlongationAndUpdate(fuseProlongationAndUpdate_);

  H.Clear();

  // There are few issues with using Keep in the interpreter:
  //   1. Hierarchy::Keep interface takes a name and a factory. If
  //      factories are different on different levels, the AddNewLevel() call
  //      in Hierarchy does not work properly, as it assume that factories are
  //      the same.
  //   2. FactoryManager does not have a Keep option, only Hierarchy and
  //      Level have it
  //   3. Interpreter constructs factory managers, but not levels. So we
  //      cannot set up Keep flags there.
  //
  // The solution implemented here does the following:
  //   1. Construct hierarchy with dummy levels. This avoids
  //      Hierarchy::AddNewLevel() calls which will propagate wrong
  //      inheritance.
  //   2. Interpreter constructs keep_ array with names and factories for
  //      that level
  //   3. For each level, we call Keep(name, factory) for each keep_
  for (int i = 0; i < numDesiredLevel_; i++) {
    std::map<int, std::vector<keep_pair>>::const_iterator it = keep_.find(i);
    if (it != keep_.end()) {
      RCP<Level> l                        = H.GetLevel(i);
      const std::vector<keep_pair>& keeps = it->second;
      for (size_t j = 0; j < keeps.size(); j++)
        l->Keep(keeps[j].first, keeps[j].second);
    }
    if (i < numDesiredLevel_ - 1) {
      RCP<Level> newLevel = rcp(new Level());
      H.AddLevel(newLevel);
    }
  }

  // Matrices to print
  for (auto iter = matricesToPrint_.begin(); iter != matricesToPrint_.end(); iter++)
    ExportDataSetKeepFlags(H, iter->second, iter->first);

  // Vectors, aggregates and other things that need special case handling
  ExportDataSetKeepFlags(H, nullspaceToPrint_, "Nullspace");
  ExportDataSetKeepFlags(H, coordinatesToPrint_, "Coordinates");
  // NOTE: Aggregates use the next level's Factory
  ExportDataSetKeepFlagsNextLevel(H, aggregatesToPrint_, "Aggregates");
#ifdef HAVE_MUELU_INTREPID2
  ExportDataSetKeepFlags(H, elementToNodeMapsToPrint_, "pcoarsen: element to node map");
#endif

  // Data to save only (these do not have a level, so we do all levels)
  for (int i = 0; i < dataToSave_.size(); i++)
    ExportDataSetKeepFlagsAll(H, dataToSave_[i]);

  int levelID      = 0;
  int lastLevelID  = numDesiredLevel_ - 1;
  bool isLastLevel = false;

  while (!isLastLevel) {
    bool r = H.Setup(levelID,
                     LvlMngr(levelID - 1, lastLevelID),
                     LvlMngr(levelID, lastLevelID),
                     LvlMngr(levelID + 1, lastLevelID));
    if (levelID < H.GetNumLevels())
      H.GetLevel(levelID)->print(H.GetOStream(Developer), verbosity_);

    isLastLevel = r || (levelID == lastLevelID);
    levelID++;
  }
  if (!matvecParams_.is_null())
    H.SetMatvecParams(matvecParams_);
  H.AllocateLevelMultiVectors(sizeOfMultiVectors_);
  // Set hierarchy description.
  // This is cached, but involves and MPI_Allreduce.
  H.description();
  H.describe(H.GetOStream(Runtime0), verbosity_);
  H.CheckForEmptySmoothersAndCoarseSolve();

  // When we reuse hierarchy, it is necessary that we don't
  // change the number of levels. We also cannot make requests
  // for coarser levels, because we don't construct all the
  // data on previous levels. For instance, let's say our first
  // run constructed three levels. If we try to do requests during
  // next setup for the fourth level, it would need Aggregates
  // which we didn't construct for level 3 because we reused P.
  // To fix this situation, we change the number of desired levels
  // here.
  numDesiredLevel_ = levelID;

  // Matrix prints
  for (auto iter = matricesToPrint_.begin(); iter != matricesToPrint_.end(); iter++) {
    WriteData<Matrix>(H, iter->second, iter->first);
  }

  // Vectors, aggregates and all things we need to print manually
  WriteData<MultiVector>(H, nullspaceToPrint_, "Nullspace");
  WriteData<MultiVector>(H, coordinatesToPrint_, "Coordinates");
  WriteDataAggregates(H, aggregatesToPrint_, "Aggregates");

#ifdef HAVE_MUELU_INTREPID2
  typedef Kokkos::DynRankView<LocalOrdinal, typename Node::device_type> FCi;
  WriteDataFC<FCi>(H, elementToNodeMapsToPrint_, "pcoarsen: element to node map", "el2node");
#endif

}  // SetupHierarchy

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<FactoryManagerBase> MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::LvlMngr(int levelID, int lastLevelID) const {
  // NOTE: the order of 'if' statements is important
  if (levelID == -1)  // levelID = -1 corresponds to the finest level
    return Teuchos::null;

  if (levelID == lastLevelID + 1)  // levelID = 'lastLevelID+1' corresponds to the last level (i.e., no nextLevel)
    return Teuchos::null;

  if (levelManagers_.size() == 0) {  // default factory manager.
    // The default manager is shared across levels, initialized only if needed and deleted with the HierarchyManager
    static RCP<FactoryManagerBase> defaultMngr = rcp(new FactoryManager());
    return defaultMngr;
  }

  return GetFactoryManager(levelID);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ExportDataSetKeepFlags(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name) const {
  for (int i = 0; i < data.size(); ++i) {
    if (data[i] < H.GetNumLevels()) {
      RCP<Level> L = H.GetLevel(data[i]);
      if (!L.is_null() && data[i] < levelManagers_.size())
        L->AddKeepFlag(name, &*levelManagers_[data[i]]->GetFactory(name));
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ExportDataSetKeepFlagsNextLevel(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name) const {
  for (int i = 0; i < data.size(); ++i) {
    if (data[i] < H.GetNumLevels()) {
      RCP<Level> L = H.GetLevel(data[i]);
      if (!L.is_null() && data[i] + 1 < levelManagers_.size())
        L->AddKeepFlag(name, &*levelManagers_[data[i] + 1]->GetFactory(name));
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ExportDataSetKeepFlagsAll(Hierarchy& H, const std::string& name) const {
  for (int i = 0; i < H.GetNumLevels(); i++) {
    RCP<Level> L = H.GetLevel(i);
    if (!L.is_null() && i < levelManagers_.size())
      L->AddKeepFlag(name, &*levelManagers_[i]->GetFactory(name));
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template <class T>
void MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::WriteData(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name) const {
  for (int i = 0; i < data.size(); ++i) {
    std::string fileName;
    if (H.getObjectLabel() != "")
      fileName = H.getObjectLabel() + "_" + name + "_" + Teuchos::toString(data[i]) + ".m";
    else
      fileName = name + "_" + Teuchos::toString(data[i]) + ".m";

    if (data[i] < H.GetNumLevels()) {
      RCP<Level> L = H.GetLevel(data[i]);
      if (data[i] < levelManagers_.size() && L->IsAvailable(name, &*levelManagers_[data[i]]->GetFactory(name))) {
        // Try generating factory
        RCP<T> M = L->template Get<RCP<T>>(name, &*levelManagers_[data[i]]->GetFactory(name));
        if (!M.is_null()) {
          Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write(fileName, *M);
        }
      } else if (L->IsAvailable(name)) {
        // Try nofactory
        RCP<T> M = L->template Get<RCP<T>>(name);
        if (!M.is_null()) {
          Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write(fileName, *M);
        }
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::WriteDataAggregates(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name) const {
  for (int i = 0; i < data.size(); ++i) {
    const std::string fileName = name + "_" + Teuchos::toString(data[i]) + ".m";

    if (data[i] < H.GetNumLevels()) {
      RCP<Level> L = H.GetLevel(data[i]);

      // NOTE: Aggregates use the next level's factory
      RCP<Aggregates> agg;
      if (data[i] + 1 < H.GetNumLevels() && L->IsAvailable(name, &*levelManagers_[data[i] + 1]->GetFactory(name))) {
        // Try generating factory
        agg = L->template Get<RCP<Aggregates>>(name, &*levelManagers_[data[i] + 1]->GetFactory(name));
      } else if (L->IsAvailable(name)) {
        agg = L->template Get<RCP<Aggregates>>("Aggregates");
      }
      if (!agg.is_null()) {
        std::ofstream ofs(fileName);
        Teuchos::FancyOStream fofs(rcp(&ofs, false));
        agg->print(fofs, Teuchos::VERB_EXTREME);
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template <class T>
void MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::WriteDataFC(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name, const std::string& ofname) const {
  for (int i = 0; i < data.size(); ++i) {
    const std::string fileName = ofname + "_" + Teuchos::toString(data[i]) + ".m";

    if (data[i] < H.GetNumLevels()) {
      RCP<Level> L = H.GetLevel(data[i]);

      if (L->IsAvailable(name)) {
        RCP<T> M = L->template Get<RCP<T>>(name);
        if (!M.is_null()) {
          RCP<Matrix> A          = L->template Get<RCP<Matrix>>("A");
          RCP<const CrsGraph> AG = A->getCrsGraph();
          WriteFieldContainer<T>(fileName, *M, *AG->getColMap());
        }
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template <class T>
void MueLu::HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::WriteFieldContainer(const std::string& fileName, T& fcont, const Map& colMap) const {
  size_t num_els  = (size_t)fcont.extent(0);
  size_t num_vecs = (size_t)fcont.extent(1);

  // Generate rowMap
  Teuchos::RCP<const Map> rowMap = Xpetra::MapFactory<LO, GO, NO>::Build(colMap.lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), fcont.extent(0), colMap.getIndexBase(), colMap.getComm());

  // Fill multivector to use *petra dump routines
  RCP<GOMultiVector> vec = Xpetra::MultiVectorFactory<GO, LO, GO, NO>::Build(rowMap, num_vecs);

  for (size_t j = 0; j < num_vecs; j++) {
    Teuchos::ArrayRCP<GO> v = vec->getDataNonConst(j);
    for (size_t i = 0; i < num_els; i++)
      v[i] = colMap.getGlobalElement(fcont(i, j));
  }

  Xpetra::IO<SC, LO, GO, NO>::WriteGOMV(fileName, *vec);
}

}  // namespace MueLu

#endif
