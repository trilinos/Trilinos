// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_SOLVERMAP_LINEARPROBLEM_DECL_HPP
#define TPETRA_SOLVERMAP_LINEARPROBLEM_DECL_HPP

/// \file Tpetra_SolverMap_LinearProblem_decl.hpp
/// \brief Declaration of the Tpetra::SolverMap_LinearProblem class

#include <Tpetra_Transform.hpp>
#include <Tpetra_LinearProblem.hpp>
#include <Tpetra_SolverMap_CrsMatrix.hpp>

namespace Tpetra {

///
/** Constructs a LinearProblem with a "fixed" Column Map for the CrsMatrix.
 *  Almost entirely a view except for the "fixed" Tpetra_CrsGraph.
 */
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class SolverMap_LinearProblem : public StructuralSameTypeTransform<LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> > {
 public:
  using NewType      = typename StructuralSameTypeTransform<LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::NewType;
  using OriginalType = typename StructuralSameTypeTransform<LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::OriginalType;

  ///
  /** Constructor
   */
  SolverMap_LinearProblem();

  ///
  /** Destructor
   */
  ~SolverMap_LinearProblem();

  ///
  /** Constructs "fixed" Tpetra::LinearProblem
   */
  NewType operator()(OriginalType const& origProblem);

 private:
  SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> solverMapCrsMatrixTrans_;
};

}  // namespace Tpetra

#endif  // TPETRA_SOLVERMAP_LINEARPROBLEM_DECL_HPP
