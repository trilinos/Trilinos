// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_DENSECONSTRAINT_DECL_HPP
#define MUELU_DENSECONSTRAINT_DECL_HPP

#include "Teuchos_ScalarTraits.hpp"

#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsGraph_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Constraint.hpp"

namespace MueLu {

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class DenseConstraint
  : public Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_DENSECONSTRAINT_SHORT
#include "MueLu_UseShortNames.hpp"
 public:
  /*!
    @class Dense constraint class.
    @brief Class which contains the constraint space details
    */

  using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

  //! @name Constructors
  //@{

  DenseConstraint() = default;

  /*! @brief Constructor
      \param B -- Fine nullspace vectors
      \param Bc -- Coarse nullspace vectors
      \param Ppattern -- Nonzero sparsity pattern for the prolongator
   */
  DenseConstraint(const RCP<MultiVector>& B, const RCP<MultiVector>& Bc, RCP<const CrsGraph> Ppattern, const std::string& solverType);

  //@}

  void Setup();

  //! Compute norm of residual B - P*Bc.
  MagnitudeType ResidualNorm(RCP<const Matrix> P) const override;

 private:
  //! Fine nullspace
  RCP<MultiVector> B_;

  //! Coarse nullspace
  RCP<MultiVector> Bc_;
};

}  // namespace MueLu

#define MUELU_DENSECONSTRAINT_SHORT
#endif
