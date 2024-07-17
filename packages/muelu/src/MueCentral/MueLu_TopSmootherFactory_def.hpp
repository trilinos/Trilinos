// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * MueLu_TopSmootherFactory_def.hpp
 *
 *  Created on: Jan 25, 2016
 *      Author: tawiesn
 */

#ifndef PACKAGES_MUELU_SRC_MUECENTRAL_MUELU_TOPSMOOTHERFACTORY_DEF_HPP_
#define PACKAGES_MUELU_SRC_MUECENTRAL_MUELU_TOPSMOOTHERFACTORY_DEF_HPP_

#include "MueLu_ConfigDefs.hpp"

//#include "MueLu_FactoryManager_fwd.hpp"
#include "MueLu_FactoryManagerBase.hpp"
//#include "MueLu_HierarchyHelpers_fwd.hpp"
#include "MueLu_TopSmootherFactory_decl.hpp"
//#include "MueLu_SmootherBase_fwd.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_SmootherPrototype.hpp"
//#include "MueLu_Hierarchy_fwd.hpp"
//#include "MueLu_HierarchyManager_fwd.hpp"

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TopSmootherFactory(RCP<const FactoryManagerBase> parentFactoryManager, const std::string& varName) {
  TEUCHOS_TEST_FOR_EXCEPTION(varName != "CoarseSolver" && varName != "Smoother", Exceptions::RuntimeError, "varName should be either \"CoarseSolver\" or \"Smoother\"");

  if (varName == "CoarseSolver") {
    // For coarsest level, we only need one smoother/solver
    // If a user wants to do something weird there (like, solve coarsest system by using 2 forward
    // GS and 1 backward GS), one can use MergedSmoother
    RCP<const FactoryBase> coarseSolverFactory       = parentFactoryManager->GetFactory("CoarseSolver");
    RCP<const SmootherFactory> coarseSmootherFactory = Teuchos::rcp_dynamic_cast<const SmootherFactory>(coarseSolverFactory);
    if (coarseSmootherFactory != Teuchos::null) {
      RCP<SmootherPrototype> preProto;
      RCP<SmootherPrototype> postProto;
      coarseSmootherFactory->GetSmootherPrototypes(preProto, postProto);

      if (preProto == postProto)
        preSmootherFact_ = parentFactoryManager->GetFactory("CoarseSolver");
      else {
        // check whether pre- and/or post-smoothing is desired on coarsest level
        if (preProto != Teuchos::null)
          preSmootherFact_ = parentFactoryManager->GetFactory("CoarseSolver");
        if (postProto != Teuchos::null)
          postSmootherFact_ = parentFactoryManager->GetFactory("CoarseSolver");
      }
    } else  // default handling: get default direct solver as presmoother on coarsest level
      preSmootherFact_ = parentFactoryManager->GetFactory("CoarseSolver");

  } else {
    preSmootherFact_  = parentFactoryManager->GetFactory("PreSmoother");
    postSmootherFact_ = parentFactoryManager->GetFactory("PostSmoother");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~TopSmootherFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& level) const {
  if (preSmootherFact_ != Teuchos::null)
    level.DeclareInput("PreSmoother", preSmootherFact_.get());
  if (postSmootherFact_ != Teuchos::null)
    level.DeclareInput("PostSmoother", postSmootherFact_.get());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& level) const {
  if (preSmootherFact_.is_null() && postSmootherFact_.is_null())
    return;

  // NOTE 1: We need to set at least some keep flag for the smoothers, otherwise it is going to be removed as soon as all requests are released.
  // We choose to set the Final flag for the data. In addition, we allow this data to be retrieved by only using the name by the means
  // of using NoFactory. However, any data set with NoFactory gets UserData flag by default. We don't really want that flag, so we remove it.

  // NOTE 2: some smoother factories are tricky (see comments in MueLu::SmootherFactory
  // Sometimes, we don't know whether the factory is able to generate "PreSmoother" or "PostSmoother"
  // For the SmootherFactory, however, we are able to check that.

  if (!preSmootherFact_.is_null()) {
    // Checking for null is not sufficient, as SmootherFactory(null, something) does not generate "PreSmoother"
    bool isAble                  = true;
    RCP<const SmootherFactory> s = rcp_dynamic_cast<const SmootherFactory>(preSmootherFact_);
    if (!s.is_null()) {
      RCP<SmootherPrototype> pre, post;
      s->GetSmootherPrototypes(pre, post);
      if (pre.is_null())
        isAble = false;
    } else {
      // We assume that  if presmoother factory is not SmootherFactory, it *is* able to generate "PreSmoother"
    }

    if (isAble) {
      RCP<SmootherBase> Pre = level.Get<RCP<SmootherBase> >("PreSmoother", preSmootherFact_.get());

      level.Set("PreSmoother", Pre, NoFactory::get());

      level.AddKeepFlag("PreSmoother", NoFactory::get(), MueLu::Final);
      level.RemoveKeepFlag("PreSmoother", NoFactory::get(), MueLu::UserData);
    }
  }

  if (!postSmootherFact_.is_null()) {
    // Checking for null is not sufficient, as SmootherFactory(something, null) does not generate "PostSmoother"
    bool isAble                  = true;
    RCP<const SmootherFactory> s = rcp_dynamic_cast<const SmootherFactory>(postSmootherFact_);
    if (!s.is_null()) {
      RCP<SmootherPrototype> pre, post;
      s->GetSmootherPrototypes(pre, post);
      if (post.is_null())
        isAble = false;
    } else {
      // We assume that  if presmoother factory is not SmootherFactory, it *is* able to generate "PreSmoother"
    }

    if (isAble) {
      RCP<SmootherBase> Post = level.Get<RCP<SmootherBase> >("PostSmoother", postSmootherFact_.get());

      level.Set("PostSmoother", Post, NoFactory::get());

      level.AddKeepFlag("PostSmoother", NoFactory::get(), MueLu::Final);
      level.RemoveKeepFlag("PostSmoother", NoFactory::get(), MueLu::UserData);
    }
  }
}
}  // namespace MueLu

#endif /* PACKAGES_MUELU_SRC_MUECENTRAL_MUELU_TOPSMOOTHERFACTORY_DEF_HPP_ */
