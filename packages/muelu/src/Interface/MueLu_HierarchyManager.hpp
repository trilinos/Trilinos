// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_HIERARCHYMANAGER_DECL_HPP
#define MUELU_HIERARCHYMANAGER_DECL_HPP

#include <string>
#include <map>

#include <Teuchos_Array.hpp>

#include <Xpetra_Matrix.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_HierarchyFactory.hpp"

#include "MueLu_Hierarchy.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_Level.hpp"

namespace MueLu {

  // This class stores the configuration of a Hierarchy.
  // The class also provides an algorithm to build a Hierarchy from the configuration.
  //
  // See also: FactoryManager
  //
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void, LocalOrdinal, Node>::SparseOps>
  class HierarchyManager : public HierarchyFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> {
#undef MUELU_HIERARCHYMANAGER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //!
    HierarchyManager() : numDesiredLevel_(10), maxCoarseSize_(50), verbosity_(Medium)   // TODO: default values should be query from Hierarchy class to avoid duplication
    { }

    //!
    virtual ~HierarchyManager() { }

    // Unused
    // void AddFactoryManager(int levelID, RCP<FactoryManagerBase>& manager) {
    //   if (levelManagers_.size() < levelID + 1) levelManagers_.resize(levelID + 1);
    //   levelManagers_[levelID] = manager;
    // }

    //!
    void AddFactoryManager(int startLevel, int numDesiredLevel, RCP<FactoryManagerBase> manager) {
      const int lastLevel = startLevel + numDesiredLevel - 1;
      if (levelManagers_.size() < lastLevel + 1) levelManagers_.resize(lastLevel + 1);

      for(int iLevel = startLevel; iLevel <= lastLevel; iLevel++) {
        levelManagers_[iLevel] = manager;
      }
    }

    //!
    Teuchos::Ptr<FactoryManagerBase> GetFactoryManager(int levelID) const {
      if (levelID >= levelManagers_.size()) return levelManagers_[levelManagers_.size()-1](); // last levelManager is used for all the remaining levels.
      return levelManagers_[levelID](); // throw exception if out of bound.
    }

    //! returns number of factory managers stored in levelManagers_ vector.
    size_t getNumFactoryManagers() const {
      return levelManagers_.size();
    }

    //!
    void SetFactoryManagerCoarsestLevel(RCP<FactoryManagerBase>& manager) {
      coarsestLevelManager_ = manager;
    }

    //!
    void CheckConfig() {
      for(int i=0; i<levelManagers_.size(); i++) {
        TEUCHOS_TEST_FOR_EXCEPTION(levelManagers_[i] == Teuchos::null, Exceptions::RuntimeError, "MueLu:HierarchyConfig::CheckConfig(): Undefined configuration for level:");
      }
    }

    //@{

    virtual RCP<Hierarchy> CreateHierarchy() const {
      return rcp(new Hierarchy());
    }

    //! Setup Hierarchy object
    virtual void SetupHierarchy(Hierarchy & H) const {
      TEUCHOS_TEST_FOR_EXCEPTION(!H.GetLevel(0)->IsAvailable("A"), Exceptions::RuntimeError, "No fine level operator");

      // Setup Matrix
      // TODO: I should certainly undo this somewhere...
      RCP<Level> l = H.GetLevel(0);
      RCP<Matrix> Op = l->Get<RCP<Matrix> >("A");
      SetupMatrix(*Op);
      SetupExtra(H);

      // Setup Hierarchy
      H.SetMaxCoarseSize(maxCoarseSize_);
      H.SetDefaultVerbLevel(verbosity_);

      // TODO: coarsestLevelManager

      int  levelID     = 0;
      int  lastLevelID = numDesiredLevel_ - 1;
      bool isLastLevel = false;

      while(!isLastLevel) {
        bool r = H.Setup(levelID,
                         LvlMngr(levelID-1, lastLevelID),
                         LvlMngr(levelID,   lastLevelID),
                         LvlMngr(levelID+1, lastLevelID));

        isLastLevel = r || (levelID == lastLevelID);
        levelID++;
      }
    }

    //@}

    typedef std::map<std::string, RCP<const FactoryBase> > FactoryMap;

  protected: //TODO: access function

    //! Setup Matrix object
    virtual void SetupMatrix(Matrix & Op) const { }

    //! Setup extra data
    // TODO: merge with SetupMatrix ?
    virtual void SetupExtra(Hierarchy & H) const { }

    // TODO this was private
    // Used in SetupHierarchy() to access levelManagers_
    // Inputs i=-1 and i=size() are allowed to simplify calls to hierarchy->Setup()
    Teuchos::Ptr<FactoryManagerBase> LvlMngr(int levelID, int lastLevelID) const {

      // Please not that the order of the 'if' statements is important.

      if (levelID == -1)                    return Teuchos::null; // when this routine is called with levelID == '-1', it means that we are processing the finest Level (there is no finer level)
      if (levelID == lastLevelID+1)         return Teuchos::null; // when this routine is called with levelID == 'lastLevelID+1', it means that we are processing the last level (ie: there is no nextLevel...)

      if (0       == levelManagers_.size()) {                     // default factory manager.
        // the default manager is shared across levels, initialized only if needed and deleted with the HierarchyManager.
        static RCP<FactoryManagerBase> defaultMngr = rcp(new FactoryManager());
        return defaultMngr();
      }
      if (levelID >= levelManagers_.size()) return levelManagers_[levelManagers_.size()-1](); // last levelManager is used for all the remaining levels.

      return levelManagers_[levelID](); // throw exception if out of bound.
    }

    // Hierarchy parameters
    int                   numDesiredLevel_;
    Xpetra::global_size_t maxCoarseSize_;
    MsgType               verbosity_;

  private:
    // Levels
    Array<RCP<FactoryManagerBase> > levelManagers_;        // one FactoryManager per level. The last levelManager is used for all the remaining levels.
    RCP<FactoryManagerBase>         coarsestLevelManager_; // coarsest level manager

  }; // class HierarchyManager

} // namespace MueLu

#define MUELU_HIERARCHYMANAGER_SHORT
#endif // MUELU_HIERARCHYMANAGER_HPP

//TODO: split into _decl/_def
// TODO: default value for first param (FactoryManager()) should not be duplicated (code maintainability)
