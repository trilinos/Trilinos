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
#ifndef MUELU_HIERARCHYHELPERS_DEF_HPP
#define MUELU_HIERARCHYHELPERS_DEF_HPP

#include <Xpetra_Matrix.hpp>

#include "MueLu_HierarchyHelpers_decl.hpp"

#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_DirectSolver_decl.hpp"

//TODO/FIXME: DeclareInput(, **this**) cannot be used here

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManager) :
    PFact_ (parentFactoryManager->GetFactory("P")),
    RFact_ (parentFactoryManager->GetFactory("R")),
    AcFact_(parentFactoryManager->GetFactory("A"))
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManagerFine, RCP<const FactoryManagerBase> parentFactoryManagerCoarse) :
    PFact_ (parentFactoryManagerCoarse->GetFactory("P")),
    RFact_ (parentFactoryManagerCoarse->GetFactory("R")),
    AcFact_(parentFactoryManagerCoarse->GetFactory("A"))
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~TopRAPFactory() { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level & fineLevel, Level & coarseLevel) const {
    if (PFact_  != Teuchos::null)                                       coarseLevel.DeclareInput("P", PFact_.get());
    if (RFact_  != Teuchos::null)                                       coarseLevel.DeclareInput("R", RFact_.get());
    if ((AcFact_ != Teuchos::null) && (AcFact_ != NoFactory::getRCP())) coarseLevel.DeclareInput("A", AcFact_.get());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level & fineLevel, Level & coarseLevel) const {
    if (PFact_ != Teuchos::null) {
      RCP<Matrix> P = coarseLevel.Get<RCP<Matrix> >("P", PFact_.get());
      coarseLevel.Set           ("P", P, NoFactory::get());
      coarseLevel.AddKeepFlag   ("P", NoFactory::get(), MueLu::Final);    // FIXME2: Order of Remove/Add matter (data removed otherwise). Should do something about this
      coarseLevel.RemoveKeepFlag("P", NoFactory::get(), MueLu::UserData); // FIXME: This is a hack, I should change behavior of Level::Set() instead. FIXME3: Should not be removed if flag was there already
    }

    if (RFact_ != Teuchos::null) {
      RCP<Matrix> R = coarseLevel.Get<RCP<Matrix> >("R", RFact_.get());
      coarseLevel.Set           ("R", R, NoFactory::get());
      coarseLevel.AddKeepFlag   ("R", NoFactory::get(), MueLu::Final);
      coarseLevel.RemoveKeepFlag("R", NoFactory::get(), MueLu::UserData); // FIXME: This is a hack
    }

    if ((AcFact_ != Teuchos::null) && (AcFact_ != NoFactory::getRCP())) {
      RCP<Matrix> Ac = coarseLevel.Get<RCP<Matrix> >("A", AcFact_.get());
      coarseLevel.Set           ("A", Ac, NoFactory::get());
      coarseLevel.AddKeepFlag   ("A", NoFactory::get(), MueLu::Final);
      coarseLevel.RemoveKeepFlag("A", NoFactory::get(), MueLu::UserData); // FIXME: This is a hack
    }
  }

  //
  //
  //

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TopSmootherFactory(RCP<const FactoryManagerBase> parentFactoryManager, const std::string& varName) {
    TEUCHOS_TEST_FOR_EXCEPTION(varName != "CoarseSolver" && varName != "Smoother", Exceptions::RuntimeError, "varName should be either \"CoarseSolver\" or \"Smoother\"");

    if (varName == "CoarseSolver") {
      // For coarsest level, we only need one smoother (so that we don't call direct solver twice)
      // If a user wants to do something weird there (like, solve coarsest system by using 2 forward
      // GS and 1 backward GS), one can use MergedSmoother
      // However, the user may have special requirements. So, check the CoarseSolver object to be of type
      // SmootherFactory and analyze the pre- and post smoother prototype.
      // make sure that a direct solver is only applied once (as "presmoother")
      // a direct solver as pre- or post-smoother in combination with another type of post- or pre-smoother
      // would be allowed (but probably makes no sense, of course...)
      RCP<const FactoryBase> coarseSolverFactory = parentFactoryManager->GetFactory("CoarseSolver");
      RCP<const SmootherFactory> coarseSmootherFactory = Teuchos::rcp_dynamic_cast<const SmootherFactory>(coarseSolverFactory);
      if (coarseSmootherFactory != Teuchos::null) {
        RCP<SmootherPrototype> preProto;
        RCP<SmootherPrototype> postProto;
        coarseSmootherFactory->GetSmootherPrototypes(preProto, postProto);


        if (preProto == postProto)
          preSmootherFact_  = parentFactoryManager->GetFactory("CoarseSolver");
        else {
          // check whether pre- and/or post-smoothing is desired on coarsest level
          if(preProto != Teuchos::null)
            preSmootherFact_  = parentFactoryManager->GetFactory("CoarseSolver");
          if(postProto != Teuchos::null)
            postSmootherFact_  = parentFactoryManager->GetFactory("CoarseSolver");
        }

        // TODO remove DirectSolver header...
        /*bool bPreDirectSolver  = Teuchos::rcp_dynamic_cast<DirectSolver>(preProto)  != Teuchos::null ? true : false;
        bool bPostDirectSolver = Teuchos::rcp_dynamic_cast<DirectSolver>(postProto) != Teuchos::null ? true : false;

        if (bPreDirectSolver == true && bPostDirectSolver == true) {
          // make sure that a direct solver is only applied once!
          preSmootherFact_  = parentFactoryManager->GetFactory("CoarseSolver");
        } else {

        }*/
      }
      else  // default handling: get default direct solver as presmoother on coarsest level
        preSmootherFact_  = parentFactoryManager->GetFactory("CoarseSolver");

    } else {
      preSmootherFact_  = parentFactoryManager->GetFactory("PreSmoother");
      postSmootherFact_ = parentFactoryManager->GetFactory("PostSmoother");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~TopSmootherFactory() { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level & level) const {
    if (preSmootherFact_  != Teuchos::null)
      level.DeclareInput("PreSmoother",  preSmootherFact_.get());
    if (postSmootherFact_ != Teuchos::null)
      level.DeclareInput("PostSmoother", postSmootherFact_.get());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level & level) const {
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
      bool isAble = true;
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
        RCP<SmootherBase> Pre  = level.Get<RCP<SmootherBase> >("PreSmoother", preSmootherFact_.get());

        level.Set           ("PreSmoother", Pre, NoFactory::get());

        level.AddKeepFlag   ("PreSmoother", NoFactory::get(), MueLu::Final);
        level.RemoveKeepFlag("PreSmoother", NoFactory::get(), MueLu::UserData);
      }
    }

    if (!postSmootherFact_.is_null()) {
      // Checking for null is not sufficient, as SmootherFactory(something, null) does not generate "PostSmoother"
      bool isAble = true;
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

        level.Set           ("PostSmoother", Post, NoFactory::get());

        level.AddKeepFlag   ("PostSmoother", NoFactory::get(), MueLu::Final);
        level.RemoveKeepFlag("PostSmoother", NoFactory::get(), MueLu::UserData);
      }
    }
  }

} // namespace MueLu

#define MUELU_HIERARCHY_HELPERS_SHORT
#endif // MUELU_HIERARCHYHELPERS_DEF_HPP
