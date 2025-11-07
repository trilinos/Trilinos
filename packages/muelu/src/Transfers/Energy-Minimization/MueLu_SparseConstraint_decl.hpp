// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SPARSECONSTRAINT_DECL_HPP
#define MUELU_SPARSECONSTRAINT_DECL_HPP

#include "Teuchos_ScalarTraits.hpp"

#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsGraph_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include <MueLu_Utilities_fwd.hpp>
#include "MueLu_Constraint.hpp"

namespace MueLu {

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class SparseConstraint
  : public Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_SPARSECONSTRAINT_SHORT
#include "MueLu_UseShortNames.hpp"
 public:
  /*!
    @class Sparse constraint class.
    @brief Class which contains the constraint space details
    */

  using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

  //! @name Constructors
  //@{

  SparseConstraint() = default;

  /*! @brief Constructor
      \param P_nodal -- Nodal prolongator matrix
      \param D -- Fine matrix
      \param Dc -- Coarse matrix
      \param Ppattern -- Nonzero sparsity pattern for the prolongator
   */
  SparseConstraint(const RCP<Matrix>& P_nodal, const RCP<Matrix>& D, const RCP<Matrix>& Dc, RCP<const CrsGraph> Ppattern, const std::string& solverType);

  //@}

  void Setup();

  MagnitudeType ResidualNorm(RCP<const Matrix> P) const override;

  void AssignMatrixEntriesToConstraintVector(const Matrix& A, MultiVector& vecC) const;

 private:
  //! Nodal prolongator
  RCP<Matrix> P_nodal_;

  //! Fine nullspace
  RCP<Matrix> D_;

  //! Coarse nullspace
  RCP<Matrix> Dc_;

  //! Pattern for RHS
  RCP<const CrsGraph> RHS_pattern_;

  //! block_is_singular(i) indicates whether matrix for ith constraint is singular
    Kokkos::View<bool*> block_is_singular_;
};

}  // namespace MueLu

#define MUELU_SPARSECONSTRAINT_SHORT
#endif
