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
#ifndef MUELU_SMOOTHERBASE_HPP
#define MUELU_SMOOTHERBASE_HPP

#include "Xpetra_MultiVector.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

namespace MueLu {
/*!
  @class SmootherBase
  @ingroup MueLuSmootherClasses
  @brief Base class for smoothers

  This has the signature for the required Apply function and contains data that is generic across all
  smoothers.
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class SmootherBase : public virtual BaseClass {
 public:
  typedef Scalar scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

 private:
#undef MUELU_SMOOTHERBASE_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //@{ Constructors/Destructors.
  SmootherBase() {}

  virtual ~SmootherBase() {}
  //@}

  //! @name Apply methods.
  //@{

  //! Apply smoother.
  virtual void Apply(MultiVector& x, const MultiVector& rhs, bool InitialGuessIsZero = false) const = 0;

  //! Compute a rough estimate of the cost to apply this smoother on this MPI rank.  Return Teuchos::OrdinalTraits<size_t>::invalid() if such an estimate cannot be computed.
  virtual size_t getNodeSmootherComplexity() const = 0;

  void declareConstructionOutcome(bool fail, std::string msg) {
    constructionSuccessful_ = !fail;
    if (!fail)
      constructionErrorMsg_ = "";
    else
      constructionErrorMsg_ = msg;
  };
  bool constructionSuccessful() { return constructionSuccessful_; }
  std::string constructionErrorMsg() { return constructionErrorMsg_; }

  //@}

 private:
  bool constructionSuccessful_;
  std::string constructionErrorMsg_;

};  // class SmootherBase
}  // namespace MueLu

#define MUELU_SMOOTHERBASE_SHORT

#endif  // ifndef MUELU_SMOOTHERBASE_HPP

// SmootherBase = Interface used by Hierarchy.Iterate(). Minimal condition to be used as smoother.
// SmootherPrototype = minimal interface used by the generic SmootherFactory.
// Note that one can implements and use his own SmootherFactory. In this case, SmootherBase is enough.
// AdvSmootherPrototype = for more complex case of reusing setup between presmoother and postsmoother
