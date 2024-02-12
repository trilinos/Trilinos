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
#ifndef MUELU_FAKESMOOTHERPROTOTYPE_HPP
#define MUELU_FAKESMOOTHERPROTOTYPE_HPP

// This is a minimal implementation of the SmootherPrototype interface. Used by unit tests. Do not use elsewhere.

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherPrototype.hpp"

namespace MueLu {

class Level;

template <class Scalar        = SmootherPrototype<>::scalar_type,
          class LocalOrdinal  = typename SmootherPrototype<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class FakeSmootherPrototype : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#include "MueLu_UseShortNames.hpp"

 public:
  FakeSmootherPrototype(int param = 0)
    : param_(param)
    , numOfSetup_(0)
    , numOfSetupCall_(0) {}

  virtual ~FakeSmootherPrototype() {}

  virtual RCP<SmootherPrototype> Copy() const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == true, Exceptions::RuntimeError, "Not a prototype. Do not copy");  // test not mandatory, but it is the only use case that we need.
    return rcp(new FakeSmootherPrototype(*this));
  }

  void DeclareInput(Level& currentLevel) const {}

  void Setup(Level&) {
    numOfSetupCall_++;
    if (SmootherPrototype::IsSetup()) return;

    numOfSetup_++;

    SmootherPrototype::IsSetup(true);
  }

  void Apply(MultiVector& x, const MultiVector& rhs, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(1, Exceptions::NotImplemented, "MueLu::FakeSmootherPrototype()::Apply(): this class is for test purpose only.")
  }

  void SetParam(int param) { param_ = param; }

  int GetParam() const { return param_; }

  int GetNumOfSetup() const { return numOfSetup_; }
  int GetNumOfSetupCall() const { return numOfSetupCall_; }

  size_t getNodeSmootherComplexity() const { return Teuchos::OrdinalTraits<size_t>::invalid(); }

 private:
  int param_;

  int numOfSetup_;
  int numOfSetupCall_;
};

}  // namespace MueLu

#define MUELU_FAKESMOOTHERPROTOTYPE_SHORT

#endif  // ifndef MUELU_FAKESMOOTHERPROTOTYPE_HPP
