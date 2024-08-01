// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_HIERARCHYMANAGER_DECL_HPP
#define MUELU_HIERARCHYMANAGER_DECL_HPP

#include <string>
#include <map>

#include <Teuchos_Array.hpp>

#include <Xpetra_Operator.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_ConfigDefs.hpp"

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

// This class stores the configuration of a Hierarchy.
// The class also provides an algorithm to build a Hierarchy from the configuration.
//
// See also: FactoryManager
//
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class HierarchyManager : public HierarchyFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_HIERARCHYMANAGER_SHORT
#include "MueLu_UseShortNames.hpp"
  typedef std::pair<std::string, const FactoryBase*> keep_pair;

 public:
  //! Constructor
  HierarchyManager(int numDesiredLevel = MasterList::getDefault<int>("max levels"))
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

  //! Destructor
  virtual ~HierarchyManager() = default;

  //!
  void AddFactoryManager(int startLevel, int numDesiredLevel, RCP<FactoryManagerBase> manager) {
    const int lastLevel = startLevel + numDesiredLevel - 1;
    if (levelManagers_.size() < lastLevel + 1)
      levelManagers_.resize(lastLevel + 1);

    for (int iLevel = startLevel; iLevel <= lastLevel; iLevel++)
      levelManagers_[iLevel] = manager;
  }

  //!
  RCP<FactoryManagerBase> GetFactoryManager(int levelID) const {
    // NOTE: last levelManager is used for all the remaining levels
    return (levelID >= levelManagers_.size() ? levelManagers_[levelManagers_.size() - 1] : levelManagers_[levelID]);
  }

  //! returns number of factory managers stored in levelManagers_ vector.
  size_t getNumFactoryManagers() const {
    return levelManagers_.size();
  }

  //!
  void CheckConfig() {
    for (int i = 0; i < levelManagers_.size(); i++)
      TEUCHOS_TEST_FOR_EXCEPTION(levelManagers_[i] == Teuchos::null, Exceptions::RuntimeError, "MueLu:HierarchyConfig::CheckConfig(): Undefined configuration for level:");
  }

  //@{

  virtual RCP<Hierarchy> CreateHierarchy() const {
    return rcp(new Hierarchy());
  }

  virtual RCP<Hierarchy> CreateHierarchy(const std::string& label) const {
    return rcp(new Hierarchy(label));
  }

  //! Setup Hierarchy object
  virtual void SetupHierarchy(Hierarchy& H) const {
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

    // Data to keep only (these do not have a level, so we do all levels)
    for (int i = 0; i < dataToKeep_.size(); i++)
      ExportDataSetKeepFlagsAll(H, dataToKeep_[i]);

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

  //! Set the number of desired levels.
  void SetNumDesiredLevel(int numDesiredLevel) { numDesiredLevel_ = numDesiredLevel; }

  //! Get the number of desired levels.
  int GetNumDesiredLevel() { return numDesiredLevel_; }

  //@}

  typedef std::map<std::string, RCP<const FactoryBase>> FactoryMap;

 protected:  // TODO: access function
  //! Setup Matrix object
  virtual void SetupOperator(Operator& /* Op */) const {}

  //! Setup extra data
  // TODO: merge with SetupMatrix ?
  virtual void SetupExtra(Hierarchy& /* H */) const {}

  // TODO this was private
  // Used in SetupHierarchy() to access levelManagers_
  // Inputs i=-1 and i=size() are allowed to simplify calls to hierarchy->Setup()
  Teuchos::RCP<FactoryManagerBase> LvlMngr(int levelID, int lastLevelID) const {
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

  //! @group Hierarchy parameters
  //! @{
  mutable int numDesiredLevel_;
  Xpetra::global_size_t maxCoarseSize_;
  MsgType verbosity_;

  bool doPRrebalance_;
  bool doPRViaCopyrebalance_;
  bool implicitTranspose_;
  bool fuseProlongationAndUpdate_;

  /*! @brief Flag to indicate whether the check of the nullspace dimension is suppressed

  By default, we do not suppress such a check, as it acts as a safety mechanism.
  Yet, certain scenarios deliberately use nullspaces with less nullspace vectors than NumPDEs.
  Therefore, the user can suppress this check. Then, the error message is converted to a warning.
  */
  bool suppressNullspaceDimensionCheck_;

  int sizeOfMultiVectors_;

  //! -2 = no output, -1 = all levels
  int graphOutputLevel_;

  //! Lists of entities to be exported (or saved)
  // Items here get handled manually
  Teuchos::Array<int> nullspaceToPrint_;
  Teuchos::Array<int> coordinatesToPrint_;
  Teuchos::Array<int> aggregatesToPrint_;
  Teuchos::Array<int> elementToNodeMapsToPrint_;

  // Data we'll need to keep, either to dump to disk or to use post-setup
  Teuchos::Array<std::string> dataToKeep_;

  // Matrices we'll need to print
  std::map<std::string, Teuchos::Array<int>> matricesToPrint_;

  Teuchos::RCP<Teuchos::ParameterList> matvecParams_;

  std::map<int, std::vector<keep_pair>> keep_;
  //! @}

 private:
  // Set the keep flags for Export Data
  void ExportDataSetKeepFlags(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name) const {
    for (int i = 0; i < data.size(); ++i) {
      if (data[i] < H.GetNumLevels()) {
        RCP<Level> L = H.GetLevel(data[i]);
        if (!L.is_null() && data[i] < levelManagers_.size())
          L->AddKeepFlag(name, &*levelManagers_[data[i]]->GetFactory(name));
      }
    }
  }

  void ExportDataSetKeepFlagsNextLevel(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name) const {
    for (int i = 0; i < data.size(); ++i) {
      if (data[i] < H.GetNumLevels()) {
        RCP<Level> L = H.GetLevel(data[i]);
        if (!L.is_null() && data[i] + 1 < levelManagers_.size())
          L->AddKeepFlag(name, &*levelManagers_[data[i] + 1]->GetFactory(name));
      }
    }
  }

  // Set the keep flags for Export Data
  void ExportDataSetKeepFlagsAll(Hierarchy& H, const std::string& name) const {
    for (int i = 0; i < H.GetNumLevels(); i++) {
      RCP<Level> L = H.GetLevel(i);
      if (!L.is_null() && i < levelManagers_.size())
        L->AddKeepFlag(name, &*levelManagers_[i]->GetFactory(name));
    }
  }

  template <class T>
  void WriteData(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name) const {
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

  void WriteDataAggregates(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name) const {
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

  template <class T>
  void WriteDataFC(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name, const std::string& ofname) const {
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

  // For dumping an IntrepidPCoarsening element-to-node map to disk
  template <class T>
  void WriteFieldContainer(const std::string& fileName, T& fcont, const Map& colMap) const {
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

  // Levels
  Array<RCP<FactoryManagerBase>> levelManagers_;  // one FactoryManager per level (the last levelManager is used for all the remaining levels)

};  // class HierarchyManager

}  // namespace MueLu

#define MUELU_HIERARCHYMANAGER_SHORT
#endif  // MUELU_HIERARCHYMANAGER_HPP

// TODO: split into _decl/_def
//  TODO: default value for first param (FactoryManager()) should not be duplicated (code maintainability)
