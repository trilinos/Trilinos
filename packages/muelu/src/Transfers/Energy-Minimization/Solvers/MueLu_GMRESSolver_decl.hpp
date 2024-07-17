// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_GMRESSOLVER_DECL_HPP
#define MUELU_GMRESSOLVER_DECL_HPP

#include <Xpetra_CrsMatrixFactory_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SolverBase.hpp"
#include "MueLu_Constraint_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

/*!
  @class GMRESSolver class.
  @brief Implements conjugate gradient algorithm for energy-minimization

  This is GMRES applied to the problem
  \f[ \min_P  \frac12 \sum_i (p_i)^T A p_i \f]
  subject to
  \f[ P(i,j) = 0 \quad \mbox{if} \; (i,j) \mbox{ is not in SparsityPattern} \f]
  and
  \f[ P cnull =  fnull. \f]

  \note
      For now we assume, that matrices have real values. The case of complex values is not explored.
  */

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class GMRESSolver : public SolverBase<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_GMRESSOLVER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{

  /*! Constructor
    \param Its -- Number of performed iterations
    */
  GMRESSolver(size_t Its);

  //@}

  //! @name Iterate methods.
  //@{

  //! Iterate
  void Iterate(const Matrix& A, const Constraint& C, const Matrix& P0, RCP<Matrix>& P) const;

  //@}

 private:
  void givapp(SC* c, SC* s, SC* v, int k) const;

 private:
  size_t nIts_;  //!< Number of performed iterations
};

}  // namespace MueLu

#define MUELU_GMRESSOLVER_SHORT
#endif  // MUELU_GMRESSOLVER_DECL_HPP
