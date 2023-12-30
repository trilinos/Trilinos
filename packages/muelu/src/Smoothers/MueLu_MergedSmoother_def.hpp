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
#ifndef MUELU_MERGEDSMOOTHER_DEF_HPP
#define MUELU_MERGEDSMOOTHER_DEF_HPP

#include "MueLu_MergedSmoother_decl.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MergedSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MergedSmoother(ArrayRCP<RCP<SmootherPrototype> >& smootherList, bool verbose)
  : smootherList_(smootherList)
  , reverseOrder_(false)
  , verbose_(verbose) {
  // TODO: check that on each method TEUCHOS_TEST_FOR_EXCEPTION(smootherList == Teuchos::null, MueLu::Exceptions::RuntimeError, "");

  SmootherPrototype::IsSetup(false);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MergedSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MergedSmoother(const MergedSmoother& src)
  : reverseOrder_(src.reverseOrder_)
  , verbose_(src.verbose_) {
  // Deep copy of src.smootherList_
  smootherList_ = SmootherListDeepCopy(src.GetSmootherList());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MergedSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetFactory(const std::string& varName, const RCP<const FactoryBase>& factory) {
  // We need to propagate SetFactory to proper place
  for (typename ArrayView<RCP<SmootherPrototype> >::iterator it = smootherList_.begin(); it != smootherList_.end(); it++)
    (*it)->SetFactory(varName, factory);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MergedSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  for (typename ArrayView<RCP<SmootherPrototype> >::iterator it = smootherList_.begin(); it != smootherList_.end(); ++it)
    (*it)->DeclareInput(currentLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MergedSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& level) {
  if (SmootherPrototype::IsSetup() == true)
    this->GetOStream(Warnings0) << "MueLu::MergedSmoother::Setup(): Setup() has already been called";

  for (typename ArrayView<RCP<SmootherPrototype> >::iterator it = smootherList_.begin(); it != smootherList_.end(); ++it) {
    try {
      (*it)->Setup(level);

    } catch (MueLu::Exceptions::RuntimeError& e) {
      std::string msg = "MueLu::MergedSmoother<>::Setup(): Runtime Error.\n One of the underlying smoother throwed the following exception: \n";
      msg += e.what();
      throw MueLu::Exceptions::RuntimeError(msg);
    }
  }

  SmootherPrototype::IsSetup(true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MergedSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, MueLu::Exceptions::RuntimeError, "MueLu::MergedSmoother<>:Apply(): Setup() has not been called");

  typedef typename ArrayRCP<RCP<SmootherPrototype> >::size_type sz_t;
  sz_t n = smootherList_.size(), c = (reverseOrder_ ? n - 1 : 0);
  char d = (reverseOrder_ ? -1 : 1);

  for (sz_t i = 0; i < n; i++)  // loop unifying both forward and reverse order
    try {
      // Be careful with nonnegative numbers
      smootherList_[c + d * Teuchos::as<char>(i)]->Apply(X, B, InitialGuessIsZero);

      // For second and later iterations, initial guess = previous result
      InitialGuessIsZero = false;

    } catch (MueLu::Exceptions::RuntimeError& e) {
      std::string msg = "MueLu::MergedSmoother<>::Apply(): Runtime Error. One of the underlying smoothers throws the following exception: \n";
      msg += e.what();
      throw MueLu::Exceptions::RuntimeError(msg);
    }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MergedSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& /* out */, const VerbLevel /* verbLevel */) const {
  throw Exceptions::NotImplemented("MueLu::MergedSmoother<>::Print() is not implemented");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MergedSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CopyParameters(RCP<SmootherPrototype> src) {  // TODO: wrong prototype. We do not need an RCP here.
  RCP<MergedSmoother> srcMergedSmoother = rcp_dynamic_cast<MergedSmoother>(src);                              // TODO: check if dynamic cast fails

  reverseOrder_ = srcMergedSmoother->GetReverseOrder();

  {
    const ArrayRCP<const RCP<SmootherPrototype> >& srcSmootherList  = srcMergedSmoother->GetSmootherList();
    const ArrayRCP<const RCP<SmootherPrototype> >& thisSmootherList = smootherList_;
    TEUCHOS_TEST_FOR_EXCEPTION(srcSmootherList == Teuchos::null, MueLu::Exceptions::RuntimeError,  // might be allowed later if useful
                               "MueLu::MergedSmoother<>:CopyParameters(): thisSmootherList == Teuchos::null");

    // If the smootherList of 'this' and 'src' contains the same type of smoothers,
    // we can transfert parameters from src to 'this' in order to tentatively reuse
    // the current setup information of each smoothers. Note that the reuse of the
    // setup phase of the MergedSmoother 'src' can be implemented for a larger set
    // of cases (and more complicated cases), but it does not seems useful for now.

    bool reuse = true;  // true == can we transfert parameters of smoothers one by one or do we have to copy the whole list of src?

    // test 1: same list size
    reuse = reuse && (thisSmootherList.size() == srcSmootherList.size());

    if (reuse) {
      //  test 2: one-by-one comparison of smoother types
      for (typename ArrayRCP<RCP<SmootherPrototype> >::size_type i = 0; i < srcSmootherList.size(); i++) {
        // The following test should never throw in our use cases because 'src' is a prototype and
        // 'this' is a real smoother so they don't share any data. We may allow such case later if useful.
        TEUCHOS_TEST_FOR_EXCEPTION((thisSmootherList[i] == srcSmootherList[i]) && (thisSmootherList[i] != Teuchos::null), MueLu::Exceptions::RuntimeError,
                                   "MueLu::MergedSmoother<>:CopyParameters(): internal logic error");

        // TODO
        //  reuse = reuse && ((thisSmootherList[i] == Teuchos::null && srcSmootherList[i] == Teuchos::null) ||
        //                    thisSmootherList[i]->GetType() == srcSmootherList[i]->GetType());
      }
    }

    reuse = false;  // TODO: temporary disactivated.

    if (reuse) {
      bool isSetup = true;

      // Call CopyParameters for each smoothers and update IsSetup status of the MergedSmoother
      for (typename ArrayRCP<RCP<SmootherPrototype> >::size_type i = 0; i < srcSmootherList.size(); i++) {
        if (srcSmootherList[i] != Teuchos::null) {
          TEUCHOS_TEST_FOR_EXCEPTION(srcSmootherList[i] == Teuchos::null, MueLu::Exceptions::RuntimeError, "MueLu::MergedSmoother<>:CopyParameters(): internal logic error");

          // TODO              thisSmootherList[i]->CopyParameters(srcSmootherList[i]);
          isSetup = isSetup && thisSmootherList[i]->IsSetup();
        }
        SmootherPrototype::IsSetup(isSetup);
      }

    } else {
      // No reuse: copy srcSmootherList.
      smootherList_ = Teuchos::null;
      smootherList_ = SmootherListDeepCopy(srcSmootherList);
      SmootherPrototype::IsSetup(false);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
MergedSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Copy() const {
  return rcp(new MergedSmoother(*this));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >
MergedSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    SmootherListDeepCopy(const ArrayRCP<const RCP<SmootherPrototype> >& srcSmootherList) {
  ArrayRCP<RCP<SmootherPrototype> > newSmootherList(srcSmootherList.size());

  for (typename ArrayRCP<RCP<SmootherPrototype> >::size_type i = 0; i < srcSmootherList.size(); i++)
    newSmootherList[i] = srcSmootherList[i]->Copy();

  return newSmootherList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t MergedSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getNodeSmootherComplexity() const {
  // FIXME: This is a placeholder
  return Teuchos::OrdinalTraits<size_t>::invalid();
}

}  // namespace MueLu

#endif  // MUELU_MERGEDSMOOTHER_DEF_HPP
