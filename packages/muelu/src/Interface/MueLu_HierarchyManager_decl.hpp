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
  HierarchyManager(int numDesiredLevel = MasterList::getDefault<int>("max levels"));

  //! Destructor
  virtual ~HierarchyManager() = default;

  //!
  void AddFactoryManager(int startLevel, int numDesiredLevel, RCP<FactoryManagerBase> manager);

  //!
  RCP<FactoryManagerBase> GetFactoryManager(int levelID) const;

  //! returns number of factory managers stored in levelManagers_ vector.
  size_t getNumFactoryManagers() const;

  //!
  void CheckConfig();

  //@{

  virtual RCP<Hierarchy> CreateHierarchy() const;

  virtual RCP<Hierarchy> CreateHierarchy(const std::string& label) const;

  //! Setup Hierarchy object
  virtual void SetupHierarchy(Hierarchy& H) const;

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
  Teuchos::RCP<FactoryManagerBase> LvlMngr(int levelID, int lastLevelID) const;

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

  // Data we'll need to save, not necessarily print
  Teuchos::Array<std::string> dataToSave_;

  // Matrices we'll need to print
  std::map<std::string, Teuchos::Array<int>> matricesToPrint_;

  Teuchos::RCP<Teuchos::ParameterList> matvecParams_;

  std::map<int, std::vector<keep_pair>> keep_;
  //! @}

 private:
  // Set the keep flags for Export Data
  void ExportDataSetKeepFlags(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name) const;

  void ExportDataSetKeepFlagsNextLevel(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name) const;

  // Set the keep flags for Export Data
  void ExportDataSetKeepFlagsAll(Hierarchy& H, const std::string& name) const;

  template <class T>
  void WriteData(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name) const;

  void WriteDataAggregates(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name) const;

  template <class T>
  void WriteDataFC(Hierarchy& H, const Teuchos::Array<int>& data, const std::string& name, const std::string& ofname) const;

  // For dumping an IntrepidPCoarsening element-to-node map to disk
  template <class T>
  void WriteFieldContainer(const std::string& fileName, T& fcont, const Map& colMap) const;

  // Levels
  Array<RCP<FactoryManagerBase>> levelManagers_;  // one FactoryManager per level (the last levelManager is used for all the remaining levels)

};  // class HierarchyManager

}  // namespace MueLu

#define MUELU_HIERARCHYMANAGER_SHORT
#endif  // MUELU_HIERARCHYMANAGER_HPP

// TODO: split into _decl/_def
//  TODO: default value for first param (FactoryManager()) should not be duplicated (code maintainability)
