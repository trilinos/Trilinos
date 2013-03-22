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
#ifndef MUELU_SMOOTHERFACTORY_DEF_HPP
#define MUELU_SMOOTHERFACTORY_DEF_HPP

#include "MueLu_SmootherFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_SmootherPrototype.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SmootherFactory(RCP<SmootherPrototype> preAndPostSmootherPrototype)
    : preSmootherPrototype_(preAndPostSmootherPrototype), postSmootherPrototype_(preAndPostSmootherPrototype)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(preAndPostSmootherPrototype != Teuchos::null && preAndPostSmootherPrototype->IsSetup() == true, Exceptions::RuntimeError, "preAndPostSmootherPrototype is not a smoother prototype (IsSetup() == true)");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SmootherFactory(RCP<SmootherPrototype> preSmootherPrototype, RCP<SmootherPrototype> postSmootherPrototype)
    : preSmootherPrototype_(preSmootherPrototype), postSmootherPrototype_(postSmootherPrototype)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(preSmootherPrototype  != Teuchos::null && preSmootherPrototype->IsSetup()  == true, Exceptions::RuntimeError, "preSmootherPrototype is not a smoother prototype (IsSetup() == true)");
    TEUCHOS_TEST_FOR_EXCEPTION(postSmootherPrototype != Teuchos::null && postSmootherPrototype->IsSetup() == true, Exceptions::RuntimeError, "postSmootherPrototype is not a smoother prototype (IsSetup() == true)");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~SmootherFactory() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetSmootherPrototypes(RCP<SmootherPrototype> preAndPostSmootherPrototype) {
    TEUCHOS_TEST_FOR_EXCEPTION(preAndPostSmootherPrototype != Teuchos::null && preAndPostSmootherPrototype->IsSetup() == true, Exceptions::RuntimeError, "preAndPostSmootherPrototype is not a smoother prototype (IsSetup() == true)");

    preSmootherPrototype_ = preAndPostSmootherPrototype;
    postSmootherPrototype_ = preAndPostSmootherPrototype;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetSmootherPrototypes(RCP<SmootherPrototype> preSmootherPrototype, RCP<SmootherPrototype> postSmootherPrototype) {
    TEUCHOS_TEST_FOR_EXCEPTION(preSmootherPrototype  != Teuchos::null && preSmootherPrototype->IsSetup()  == true, Exceptions::RuntimeError, "preSmootherPrototype is not a smoother prototype (IsSetup() == true)");
    TEUCHOS_TEST_FOR_EXCEPTION(postSmootherPrototype != Teuchos::null && postSmootherPrototype->IsSetup() == true, Exceptions::RuntimeError, "postSmootherPrototype is not a smoother prototype (IsSetup() == true)");
    preSmootherPrototype_ = preSmootherPrototype;
    postSmootherPrototype_ = postSmootherPrototype;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetSmootherPrototypes(RCP<SmootherPrototype> & preSmootherPrototype, RCP<SmootherPrototype> & postSmootherPrototype) const {
    preSmootherPrototype = preSmootherPrototype_;
    postSmootherPrototype = postSmootherPrototype_;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    if (preSmootherPrototype_ != Teuchos::null) {
      preSmootherPrototype_->DeclareInput(currentLevel);
    }
    if ((postSmootherPrototype_ != Teuchos::null) && (preSmootherPrototype_ != postSmootherPrototype_)) {
      postSmootherPrototype_->DeclareInput(currentLevel);
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& currentLevel) const {
    return BuildSmoother(currentLevel, BOTH);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildSmoother(Level & currentLevel, PreOrPost const preOrPost) const {
    RCP<SmootherPrototype> preSmoother;
    RCP<SmootherPrototype> postSmoother;

    if ((preOrPost == BOTH || preOrPost == PRE) && (preSmootherPrototype_ != Teuchos::null)) {
      preSmoother = preSmootherPrototype_->Copy();
      //preSmoother = rcp( new SmootherPrototype(preSmootherPrototype_) );
      //TODO if outputlevel high enough
      //TODO preSmoother.Print();
      preSmoother->Setup(currentLevel);
      currentLevel.Release(*preSmoother);

      // Level Set
      currentLevel.Set<RCP<SmootherBase> >("PreSmoother", preSmoother, this);
    }

    if ((preOrPost == BOTH || preOrPost == POST) && (postSmootherPrototype_ != Teuchos::null))
      {
        if (preOrPost == BOTH && preSmootherPrototype_ == postSmootherPrototype_) {

          // Very simple reuse. TODO: should be done in MueMat too
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

          // NO reuse: preOrPost==POST or post-smoother != pre-smoother
          // Copy the prototype and run the setup phase.
          postSmoother = postSmootherPrototype_->Copy();
          postSmoother->Setup(currentLevel);
          currentLevel.Release(*postSmoother);
        }

        // Level Set
        currentLevel.Set<RCP<SmootherBase> >("PostSmoother", postSmoother, this);
      }

  } // Build()

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherFactoryBase::description();
    std::string preStr  = (preSmootherPrototype_ == Teuchos::null) ? "null" : preSmootherPrototype_->description();
    std::string postStr = (preSmootherPrototype_ == postSmootherPrototype_) ? "pre" : ( (postSmootherPrototype_ == Teuchos::null) ? "null" : preSmootherPrototype_->description() );
    out << "{pre = "  << preStr << ", post = "<< postStr << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "PreSmoother : "; if (preSmootherPrototype_  == Teuchos::null) { out0 << "null" << std::endl; } else { Teuchos::OSTab tab2(out); preSmootherPrototype_->describe(out, verbLevel); }

      out0 << "PostSmoother: ";
      if      (postSmootherPrototype_ == preSmootherPrototype_) { out0 << "same as PreSmoother" << std::endl; }
      else if (postSmootherPrototype_ == Teuchos::null)         { out0 << "null" << std::endl; }
      else {
        { Teuchos::OSTab tab2(out); postSmootherPrototype_->describe(out, verbLevel); }
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
