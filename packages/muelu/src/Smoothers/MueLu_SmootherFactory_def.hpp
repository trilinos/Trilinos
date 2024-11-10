// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SMOOTHERFACTORY_DEF_HPP
#define MUELU_SMOOTHERFACTORY_DEF_HPP

#include "MueLu_SmootherFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_SmootherPrototype.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~SmootherFactory() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SmootherFactory(RCP<SmootherPrototype> preAndPostSmootherPrototype) {
  SetSmootherPrototypes(preAndPostSmootherPrototype);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SmootherFactory(RCP<SmootherPrototype> preSmootherPrototype,
                                                                            RCP<SmootherPrototype> postSmootherPrototype) {
  SetSmootherPrototypes(preSmootherPrototype, postSmootherPrototype);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetSmootherPrototypes(RCP<SmootherPrototype> preAndPostSmootherPrototype) {
  preSmootherPrototype_ = postSmootherPrototype_ = preAndPostSmootherPrototype;
  CheckPrototypes();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetSmootherPrototypes(RCP<SmootherPrototype> preSmootherPrototype,
                                                                                       RCP<SmootherPrototype> postSmootherPrototype) {
  preSmootherPrototype_  = preSmootherPrototype;
  postSmootherPrototype_ = postSmootherPrototype;
  CheckPrototypes();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<bool>("keep smoother data", false, "Keep constructed smoothers for later reuse");

  validParamList->set<RCP<SmootherPrototype> >("PreSmoother data", null, "Pre-smoother data for reuse");
  validParamList->set<RCP<SmootherPrototype> >("PostSmoother data", null, "Post-smoother data for reuse");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CheckPrototypes() const {
  TEUCHOS_TEST_FOR_EXCEPTION(preSmootherPrototype_ != Teuchos::null && preSmootherPrototype_->IsSetup() == true,
                             Exceptions::RuntimeError, "preSmoother prototype is not a smoother prototype (IsSetup() == true)");
  TEUCHOS_TEST_FOR_EXCEPTION(postSmootherPrototype_ != Teuchos::null && postSmootherPrototype_->IsSetup() == true,
                             Exceptions::RuntimeError, "postSmoother prototype is not a smoother prototype (IsSetup() == true)");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetSmootherPrototypes(RCP<SmootherPrototype>& preSmootherPrototype,
                                                                                       RCP<SmootherPrototype>& postSmootherPrototype) const {
  preSmootherPrototype  = preSmootherPrototype_;
  postSmootherPrototype = postSmootherPrototype_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  if (preSmootherPrototype_ != Teuchos::null)
    preSmootherPrototype_->DeclareInput(currentLevel);

  if ((postSmootherPrototype_ != Teuchos::null) && (preSmootherPrototype_ != postSmootherPrototype_))
    postSmootherPrototype_->DeclareInput(currentLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  return BuildSmoother(currentLevel, BOTH);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildSmoother(Level& currentLevel, PreOrPost const preOrPost) const {
  // SmootherFactory is quite tricky because of the fact that one of the smoother prototypes may be zero.
  // The challenge is that we have no way of knowing how user uses this factory. For instance, lets say
  // user wants to use s1 prototype as a presmoother, and s2 as a postsmoother. They could do:
  //   (a) create SmootherFactory(s1, s2), or
  //   (b) create SmootherFactory(s1, null) and SmootherFactory(null, s2)
  // It may also happen that somewhere somebody set presmoother factory = postsmoother factory = (a)
  // How do you do DeclareInput in this case? It could easily introduce a bug if a user does not check
  // whether presmoother = postsmoother. A buggy code could look like that:
  //   RCP<SmootherFactory> s = rcp(new SmootherFactory(s1,s2));
  //   level.Request("PreSmoother",  s.get());
  //   level.Request("PostSmoother", s.get());
  //   Get<RCP<SmootherBase> > pre  = Get<RCP<SmootherBase> >("PreSmoother",  s.get());
  //   Get<RCP<SmootherBase> > post = Get<RCP<SmootherBase> >("PostSmoother", s.get());
  // This code would call DeclareInput in request mode twice, but as the Build method generates both Pre and Post
  // smoothers, it would call DelcareInput in release mode only once, leaving requests.
  // This code has another problem if s2 = Teuchos::null. In that case, despite the request for PostSmoother, the factory
  // would not generate one, and second Get would throw. The real issue here is that given a Factory pointer
  // there is no way to be sure that this factory would generate any of "PreSmoother" or "PostSmoother", unless you are
  // able to cast it to SmootherFactory, do GetPrototypes and to check whether any of those is Teuchos::null.

  const Teuchos::ParameterList& pL = GetParameterList();

  RCP<SmootherPrototype> preSmoother, postSmoother;
  ParameterList preSmootherParams, postSmootherParams;

  if ((preOrPost & PRE) && !preSmootherPrototype_.is_null()) {
    if (currentLevel.IsAvailable("PreSmoother data", this))
      preSmoother = currentLevel.Get<RCP<SmootherPrototype> >("PreSmoother data", this);
    else
      preSmoother = preSmootherPrototype_->Copy();

    int oldRank = -1;
    if (!currentLevel.GetComm().is_null())
      oldRank = preSmoother->SetProcRankVerbose(currentLevel.GetComm()->getRank());

    preSmoother->Setup(currentLevel);
    preSmootherParams = preSmoother->GetParameterList();

    if (oldRank != -1)
      preSmoother->SetProcRankVerbose(oldRank);

    currentLevel.Set<RCP<SmootherBase> >("PreSmoother", preSmoother, this);

    if (pL.get<bool>("keep smoother data"))
      Set(currentLevel, "PreSmoother data", preSmoother);
  }

  if ((preOrPost & POST) && !postSmootherPrototype_.is_null()) {
    if (preOrPost == BOTH && preSmootherPrototype_ == postSmootherPrototype_) {
      // Simple reuse
      // Same prototypes for pre- and post-smoothers mean that we only need to call Setup only once
      postSmoother = preSmoother;

      //               else if (preOrPost == BOTH &&
      //                        preSmootherPrototype_ != Teuchos::null &&
      //                        preSmootherPrototype_->GetType() == postSmootherPrototype_->GetType()) {

      //               // More complex reuse case: need implementation of CopyParameters() and a smoothers smart enough to know when parameters affect the setup phase.

      //               // YES: post-smoother == pre-smoother
      //               // => copy the pre-smoother to avoid the setup phase of the post-smoother.
      //               postSmoother = preSmoother->Copy();
      //               // If the post-smoother parameters are different from
      //               // pre-smoother, the parameters stored in the post-smoother
      //               // prototype are copied in the new post-smoother object.
      //               postSmoother->CopyParameters(postSmootherPrototype_);
      //               // If parameters don't influence the Setup phase (it is the case
      //               // for Jacobi, Chebyshev...), PostSmoother is already setup. Nothing
      //               // more to do. In the case of ILU, parameters of the smoother
      //               // are in fact the parameters of the Setup phase. The call to
      //               // CopyParameters resets the smoother (only if parameters are
      //               // different) and we must call Setup() again.
      //               postSmoother->Setup(currentLevel);
      //               }

      //               // TODO: if CopyParameters do not exist, do setup twice.

    } else {
      if (currentLevel.IsAvailable("PostSmoother data", this)) {
        postSmoother = currentLevel.Get<RCP<SmootherPrototype> >("PostSmoother data", this);
      } else {
        // No reuse:
        //  - either we only do postsmoothing without any presmoothing
        //  - or our postsmoother is different from presmoother
        postSmoother = postSmootherPrototype_->Copy();
      }

      int oldRank = -1;
      if (!currentLevel.GetComm().is_null())
        oldRank = postSmoother->SetProcRankVerbose(GetProcRankVerbose());

      postSmoother->Setup(currentLevel);
      postSmootherParams = postSmoother->GetParameterList();

      if (oldRank != -1)
        postSmoother->SetProcRankVerbose(oldRank);
    }

    currentLevel.Set<RCP<SmootherBase> >("PostSmoother", postSmoother, this);

    if (pL.get<bool>("keep smoother data"))
      Set(currentLevel, "PostSmoother data", postSmoother);
  }

  ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());
  if (postSmoother == preSmoother && !preSmoother.is_null()) {
    paramList.sublist("smoother", false) = preSmoother->GetParameterList();

  } else {
    if (!preSmoother.is_null())
      paramList.sublist("presmoother", false) = preSmootherParams;

    if (!postSmoother.is_null())
      paramList.sublist("postsmoother", false) = postSmootherParams;
  }

}  // Build()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  out << SmootherFactoryBase::description();
  std::string preStr  = (preSmootherPrototype_ == Teuchos::null) ? "null" : preSmootherPrototype_->description();
  std::string postStr = (preSmootherPrototype_ == postSmootherPrototype_) ? "pre" : ((postSmootherPrototype_ == Teuchos::null) ? "null" : postSmootherPrototype_->description());
  out << "{pre = " << preStr << ", post = " << postStr << "}";
  return out.str();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::describe(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;

  if (verbLevel & Parameters0) {
    out0 << "PreSmoother : ";
    if (preSmootherPrototype_.is_null()) {
      out0 << "null" << std::endl;
    } else {
      Teuchos::OSTab tab2(out);
      preSmootherPrototype_->describe(out, verbLevel);
    }

    out0 << "PostSmoother: ";
    if (postSmootherPrototype_ == preSmootherPrototype_) {
      out0 << "same as PreSmoother" << std::endl;
    } else if (postSmootherPrototype_ == Teuchos::null) {
      out0 << "null" << std::endl;
    } else {
      Teuchos::OSTab tab2(out);
      postSmootherPrototype_->describe(out, verbLevel);
      out0 << "PostSmoother is different than PreSmoother (not the same object)" << std::endl;
    }
  }

  if (verbLevel & Debug) {
    if (preSmootherPrototype_ != Teuchos::null || postSmootherPrototype_ != Teuchos::null) {
      out0 << "-" << std::endl;
    }
    if (preSmootherPrototype_ != Teuchos::null) {
      out0 << "RCP<preSmootherPrototype_> : " << preSmootherPrototype_ << std::endl;
    }
    if (postSmootherPrototype_ != Teuchos::null) {
      out0 << "RCP<postSmootherPrototype_>: " << postSmootherPrototype_ << std::endl;
    }
  }
}

}  // namespace MueLu

// TODO: doc: setup done twice if PostSmoother object != PreSmoother object and no adv. reused capability

// TODO ReUse:   If only one smoother is missing, SmootherFactory can be smart and build only the missing smoother.
// TODO (optim): we can also reuse if preOrPost = post and preSmoother available in Level
//               we can also reuse if preOrPost = pre  and postSmoother available in Level

#endif  // MUELU_SMOOTHERFACTORY_DEF_HPP
