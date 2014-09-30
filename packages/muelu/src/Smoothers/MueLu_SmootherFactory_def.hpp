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
#ifndef MUELU_SMOOTHERFACTORY_DEF_HPP
#define MUELU_SMOOTHERFACTORY_DEF_HPP

#include "MueLu_SmootherFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Ifpack2Smoother.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SmootherFactory(RCP<SmootherPrototype> preAndPostSmootherPrototype) {
    SetSmootherPrototypes(preAndPostSmootherPrototype);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SmootherFactory(RCP<SmootherPrototype> preSmootherPrototype,
                                                                              RCP<SmootherPrototype> postSmootherPrototype) {
    SetSmootherPrototypes(preSmootherPrototype, postSmootherPrototype);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetSmootherPrototypes(RCP<SmootherPrototype> preAndPostSmootherPrototype) {
    preSmootherPrototype_  = postSmootherPrototype_ = preAndPostSmootherPrototype;
    CheckPrototypes();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetSmootherPrototypes(RCP<SmootherPrototype> preSmootherPrototype,
                                                                                         RCP<SmootherPrototype> postSmootherPrototype) {
    preSmootherPrototype_  = preSmootherPrototype;
    postSmootherPrototype_ = postSmootherPrototype;
    CheckPrototypes();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CheckPrototypes() const {
    TEUCHOS_TEST_FOR_EXCEPTION(preSmootherPrototype_  != Teuchos::null && preSmootherPrototype_->IsSetup()  == true,
                               Exceptions::RuntimeError, "preSmoother prototype is not a smoother prototype (IsSetup() == true)");
    TEUCHOS_TEST_FOR_EXCEPTION(postSmootherPrototype_ != Teuchos::null && postSmootherPrototype_->IsSetup() == true,
                               Exceptions::RuntimeError, "postSmoother prototype is not a smoother prototype (IsSetup() == true)");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetSmootherPrototypes(RCP<SmootherPrototype>& preSmootherPrototype,
                                                                                         RCP<SmootherPrototype>& postSmootherPrototype) const {
    preSmootherPrototype  = preSmootherPrototype_;
    postSmootherPrototype = postSmootherPrototype_;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    if (preSmootherPrototype_ != Teuchos::null)
      preSmootherPrototype_->DeclareInput(currentLevel);

    if ((postSmootherPrototype_ != Teuchos::null) && (preSmootherPrototype_ != postSmootherPrototype_))
      postSmootherPrototype_->DeclareInput(currentLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    return BuildSmoother(currentLevel, BOTH);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildSmoother(Level& currentLevel, PreOrPost const preOrPost) const {
    // SmootherFactory is quite tricky because of the fact that one of the smoother prototypes may be zero.
    // The challenge is that we have no way of knowing how user uses this factory. For instance, lets say
    // user wants to use s1 prototype as a presmoother, and s2 as a postsmoother. He could do:
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

    RCP<SmootherPrototype> preSmoother, postSmoother;
    ParameterList preSmootherParams, postSmootherParams;

    if ((preOrPost & PRE) && !preSmootherPrototype_.is_null()) {
      preSmoother = preSmootherPrototype_->Copy();
      preSmoother->Setup(currentLevel);
      preSmootherParams = preSmoother->GetParameterList();

      currentLevel.Set<RCP<SmootherBase> >("PreSmoother", preSmoother, this);
    }

    if ((preOrPost & POST) && !postSmootherPrototype_.is_null()) {
      if (preOrPost == BOTH && preSmootherPrototype_ == postSmootherPrototype_) {
        // Simple reuse
        // Same prototypes for pre- and post-smoothers mean that we only need to call Setup only once
        postSmoother = preSmoother;

        //            }  else if (preOrPost == BOTH &&
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

        //               // TODO: if CopyParameters do not exist, do setup twice.

      } else {
        // No reuse:
        //  - either we only do postsmoothing without any presmoothing
        //  - or our postsmoother is different from presmoother
        postSmoother = postSmootherPrototype_->Copy();
        postSmoother->Setup(currentLevel);
      }
      postSmootherParams = postSmoother->GetParameterList();

      currentLevel.Set<RCP<SmootherBase> >("PostSmoother", postSmoother, this);
    }

    ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());
    if (postSmoother == preSmoother && !preSmoother.is_null()) {
      paramList = preSmoother->GetParameterList();

    } else {
      if (!preSmoother.is_null()) {
        ParameterList& preList = paramList.sublist("presmoother", false);
        preList = preSmootherParams;
      }

      if (!postSmoother.is_null()) {
        ParameterList& postList = paramList.sublist("postsmoother", false);
        postList = postSmootherParams;
      }
    }

  } // Build()

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
    std::ostringstream out;
    out << SmootherFactoryBase::description();
    std::string preStr  = (preSmootherPrototype_ == Teuchos::null)          ? "null" : preSmootherPrototype_->description();
    std::string postStr = (preSmootherPrototype_ == postSmootherPrototype_) ? "pre"  : ( (postSmootherPrototype_ == Teuchos::null) ? "null" : postSmootherPrototype_->description() );
    out << "{pre = "  << preStr << ", post = "<< postStr << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
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
      if      (postSmootherPrototype_ == preSmootherPrototype_) { out0 << "same as PreSmoother" << std::endl; }
      else if (postSmootherPrototype_ == Teuchos::null)         { out0 << "null" << std::endl; }
      else {
        Teuchos::OSTab tab2(out);
        postSmootherPrototype_->describe(out, verbLevel);
        out0 << "PostSmoother is different than PreSmoother (not the same object)" << std::endl;
      }
    }

    if (verbLevel & Debug) {
      if (preSmootherPrototype_  != Teuchos::null || postSmootherPrototype_ != Teuchos::null) { out0 << "-" << std::endl; }
      if (preSmootherPrototype_  != Teuchos::null) { out0 << "RCP<preSmootherPrototype_> : " << preSmootherPrototype_  << std::endl; }
      if (postSmootherPrototype_ != Teuchos::null) { out0 << "RCP<postSmootherPrototype_>: " << postSmootherPrototype_ << std::endl; }
    }
  }


} // namespace MueLu

//TODO: doc: setup done twice if PostSmoother object != PreSmoother object and no adv. reused capability

// TODO ReUse:   If only one smoother is missing, SmootherFactory can be smart and build only the missing smoother.
// TODO (optim): we can also reuse if preOrPost = post and preSmoother available in Level
//               we can also reuse if preOrPost = pre  and postSmoother available in Level

#endif // MUELU_SMOOTHERFACTORY_DEF_HPP
