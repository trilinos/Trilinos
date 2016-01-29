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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_HIERARCHYUTILS_DECL_HPP
#define MUELU_HIERARCHYUTILS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_FactoryManager_fwd.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_HierarchyUtils_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_SmootherBase_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_SmootherPrototype_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Hierarchy_fwd.hpp"
#include "MueLu_HierarchyManager_fwd.hpp"

// Warning: on TopRAPFactory and TopSmootherFactory constructors, Teuchos::null doesn't mean "default factory" but "no build"

namespace MueLu {

  //! An exception safe way to call the method 'Level::SetFactoryManager()'
  class SetFactoryManager {

  public:

    //@{

    //!
    SetFactoryManager(const RCP<Level> & level, const RCP<const FactoryManagerBase> & factoryManager)
      : level_(level), prevFactoryManager_(level->GetFactoryManager())
    {
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
     // needed to save & restore previous factoryManager
    const RCP<Level> level_;
    const RCP<const FactoryManagerBase> prevFactoryManager_;
  };




  template <class Scalar,
            class LocalOrdinal  = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node          = KokkosClassic::DefaultNode::DefaultNodeType>
  class HierarchyUtils {
#undef MUELU_HIERARCHYUTILS_SHORT
#include "MueLu_UseShortNames.hpp"
  public:
    /*!
     * This routine is used by the CreateE/TpetraPreconditioner routines.
     * Adds the following non-serializable data (A,P,R,Nullspace,Coordinates) from level-specific sublist nonSerialList,
     * calling AddNewLevel as appropriate.
     */
    static void AddNonSerializableDataToHierarchy(HierarchyManager& HM, Hierarchy& H, const ParameterList& nonSerialList);
  };




} // namespace MueLu

#define MUELU_HIERARCHYUTILS_SHORT
#endif // MUELU_HIERARCHYUTILS_DECL_HPP
