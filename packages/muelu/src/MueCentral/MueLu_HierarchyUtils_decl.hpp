// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_HIERARCHYUTILS_DECL_HPP
#define MUELU_HIERARCHYUTILS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_FactoryManager_fwd.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_HierarchyUtils_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Hierarchy_fwd.hpp"
#include "MueLu_HierarchyManager_fwd.hpp"

// Warning: on TopRAPFactory and TopSmootherFactory constructors, Teuchos::null doesn't mean "default factory" but "no build"

namespace MueLu {

//! An exception safe way to call the method 'Level::SetFactoryManager()'
class SetFactoryManager {
 public:
  //@{

  /*!
    @brief Constructor

    Set a given factory manager on a specific level
  */
  SetFactoryManager(const RCP<Level>& level, const RCP<const FactoryManagerBase>& factoryManager)
    : level_(level)
    , prevFactoryManager_(level->GetFactoryManager()) {
    // set new factory manager
    level->SetFactoryManager(factoryManager);
  }

  //! Destructor.
  virtual ~SetFactoryManager() {
    // restore previous factory manager
    level_->SetFactoryManager(prevFactoryManager_);
  }

  //@}

 private:
  //! needed to save & restore previous factoryManager
  const RCP<Level> level_;
  const RCP<const FactoryManagerBase> prevFactoryManager_;
};

template <class Scalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class HierarchyUtils {
#undef MUELU_HIERARCHYUTILS_SHORT
#include "MueLu_UseShortNames.hpp"
 public:
  /*!
  \brief Add non-serializable data to Hierarchy

  Add non-serializable data given level-specific sublist \c nonSerialList to the Hierarchy \c H.
  Calling \c AddLevel() along the way, if necessary.

  Non-serializable data to be added:
  - Operator "A"
  - Prolongator "P"
  - Restrictor "R"
  - "M"
  - "Mdiag"
  - "K"
  - Nullspace information "Nullspace"
  - Coordinate information "Coordinates"
  - "Node Comm"
  - Primal-to-dual node mapping "DualNodeID2PrimalNodeID"
  - "Primal interface DOF map"
  - "pcoarsen: element to node map

  This routine is used by the CreateXpetraPreconditioner() routine.

  @param HM Hierarhcy manager
  @param H Hierarchy, where non-serializable data needs to be added
  @param nonSerialList Parameter list containing non-serializable data
  */
  static void AddNonSerializableDataToHierarchy(HierarchyManager& HM, Hierarchy& H, const ParameterList& nonSerialList);
  static void CopyBetweenHierarchies(Hierarchy& fromHierarchy, Hierarchy& toHierarchy, const std::string fromLabel, const std::string toLabel, const std::string dataType);
};

}  // namespace MueLu

#define MUELU_HIERARCHYUTILS_SHORT
#endif  // MUELU_HIERARCHYUTILS_DECL_HPP
