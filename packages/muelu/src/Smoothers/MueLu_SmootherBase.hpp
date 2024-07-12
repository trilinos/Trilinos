// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
