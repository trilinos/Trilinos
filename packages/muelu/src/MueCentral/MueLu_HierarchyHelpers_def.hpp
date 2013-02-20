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
#ifndef MUELU_HIERARCHYHELPERS_DEF_HPP
#define MUELU_HIERARCHYHELPERS_DEF_HPP

#include <Xpetra_Matrix.hpp>

#include "MueLu_HierarchyHelpers_decl.hpp"

#include "MueLu_SmootherBase.hpp"

//TODO/FIXME: DeclareInput(, **this**) cannot be used here

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManager)
    : PFact_(parentFactoryManager->GetFactory("P")), RFact_(parentFactoryManager->GetFactory("R")), AcFact_(parentFactoryManager->GetFactory("A"))
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManagerFine, RCP<const FactoryManagerBase> parentFactoryManagerCoarse)
    : PFact_(parentFactoryManagerCoarse->GetFactory("P")), RFact_(parentFactoryManagerCoarse->GetFactory("R")), AcFact_(parentFactoryManagerCoarse->GetFactory("A"))
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~TopRAPFactory() { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & fineLevel, Level & coarseLevel) const {
    if (PFact_  != Teuchos::null)                                       coarseLevel.DeclareInput("P", PFact_.get());
    if (RFact_  != Teuchos::null)                                       coarseLevel.DeclareInput("R", RFact_.get());
    if ((AcFact_ != Teuchos::null) && (AcFact_ != NoFactory::getRCP())) coarseLevel.DeclareInput("A", AcFact_.get());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TopRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level & coarseLevel) const {
    if (PFact_ != Teuchos::null) {
      RCP<Matrix> P = coarseLevel.Get<RCP<Matrix> >("P", PFact_.get());
      coarseLevel.Set("P", P, NoFactory::get());
      coarseLevel.AddKeepFlag("P", NoFactory::get(), MueLu::Final);       // FIXME2: Order of Remove/Add matter (data removed otherwise). Should do something about this
      coarseLevel.RemoveKeepFlag("P", NoFactory::get(), MueLu::UserData); // FIXME: This is a hack, I should change behavior of Level::Set() instead. FIXME3: Should not be removed if flag was there already
    }

    if (RFact_ != Teuchos::null) {
      RCP<Matrix> R = coarseLevel.Get<RCP<Matrix> >("R", RFact_.get());
      coarseLevel.Set("R", R, NoFactory::get());
      coarseLevel.AddKeepFlag("R", NoFactory::get(), MueLu::Final);
      coarseLevel.RemoveKeepFlag("R", NoFactory::get(), MueLu::UserData); // FIXME: This is a hack
    }

    if ((AcFact_ != Teuchos::null) && (AcFact_ != NoFactory::getRCP())) {
      RCP<Matrix> Ac = coarseLevel.Get<RCP<Matrix> >("A", AcFact_.get());
      coarseLevel.Set("A", Ac, NoFactory::get());
      coarseLevel.AddKeepFlag("A", NoFactory::get(), MueLu::Final);
      coarseLevel.RemoveKeepFlag("A", NoFactory::get(), MueLu::UserData); // FIXME: This is a hack
    }
  }

  //
  //
  //

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TopSmootherFactory(RCP<const FactoryManagerBase> parentFactoryManager, const std::string& varName)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(varName == "CoarseSolver" || varName == "Smoother", Exceptions::RuntimeError, "varName should be either \"CoarseSolver\" or \"Smoother\"");

    if (varName == "CoarseSolver") {
      presmootherFact_  = parentFactoryManager->GetFactory("CoarseSolver");
      postsmootherFact_ = parentFactoryManager->GetFactory("CoarseSolver");

    } else {
      presmootherFact_  = parentFactoryManager->GetFactory("PreSmoother");
      postsmootherFact_ = parentFactoryManager->GetFactory("PostSmoother");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~TopSmootherFactory() { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & level) const {

    if (presmootherFact_  != Teuchos::null)
      level.DeclareInput("PreSmoother",  presmootherFact_.get());
    if (postsmootherFact_ != Teuchos::null)
      level.DeclareInput("PostSmoother", postsmootherFact_.get());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & level) const {
    typedef MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> SmootherBase2; //TODO

    if (presmootherFact_.get() == postsmootherFact_.get()) {
      // Only call factory if at least one smoother is missing (mimic behavior of level.Get<> but level.Get<> cannot be used here as we don't know if the factory will produce both Pre and Post smoother)
      if (!level.IsAvailable("PreSmoother", presmootherFact_.get()) || !level.IsAvailable("PostSmoother", postsmootherFact_.get())) {
        presmootherFact_->CallBuild(level);
      }
    } else {
      if (!level.IsAvailable("PreSmoother", presmootherFact_.get()))
        presmootherFact_->CallBuild(level);
      if (!level.IsAvailable("PostSmoother", postsmootherFact_.get()))
        postsmootherFact_->CallBuild(level);
    }

      // Factory might or might not have created a pre smoother
      if (level.IsAvailable("PreSmoother", presmootherFact_.get())) {
        RCP<SmootherBase2> Pre  = level.Get<RCP<SmootherBase2> >("PreSmoother", presmootherFact_.get());
        level.Set           ("PreSmoother", Pre, NoFactory::get());
        level.AddKeepFlag   ("PreSmoother", NoFactory::get(), MueLu::Final);
        level.RemoveKeepFlag("PreSmoother", NoFactory::get(), MueLu::UserData); // FIXME: This is a hack
      }

      // Factory might or might not have created a post smoother
      if (level.IsAvailable("PostSmoother", postsmootherFact_.get())) {
        RCP<SmootherBase2> Post = level.Get<RCP<SmootherBase2> >("PostSmoother", postsmootherFact_.get());
        level.Set           ("PostSmoother", Post, NoFactory::get());
        level.AddKeepFlag   ("PostSmoother", NoFactory::get(), MueLu::Final);
        level.RemoveKeepFlag("PostSmoother", NoFactory::get(), MueLu::UserData); // FIXME: This is a hack
      }

    }
} // namespace MueLu

#define MUELU_HIERARCHY_HELPERS_SHORT
#endif // MUELU_HIERARCHYHELPERS_DEF_HPP
