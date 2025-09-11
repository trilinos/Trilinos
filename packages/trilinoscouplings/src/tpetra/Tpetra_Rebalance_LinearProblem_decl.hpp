// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_REBALANCE_LINEARPROBLEM_DECL_HPP
#define TPETRA_REBALANCE_LINEARPROBLEM_DECL_HPP

/// \file Tpetra_Rebalance_LinearProblem_decl.hpp
/// \brief Declaration of the Tpetra::Rebalance_LinearProblem class

#include <Tpetra_Transform.hpp>
#include <Tpetra_LinearProblem.hpp>

namespace Tpetra {

///
/** Given and input Tpetra LinearProblem, a "rebalanced" version will be returned.
 *  The data in the new T_LP is a "rebalanced" view of the original.
 */
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class Rebalance_LinearProblem : public ViewTransform< LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
{
public:
using NewType      = typename ViewTransform< LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::NewType;
using OriginalType = typename ViewTransform< LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::OriginalType;

  ///
  /** Constructor
   */
  Rebalance_LinearProblem();

  ///
  /** Destructor
   */
  ~Rebalance_LinearProblem();

  ///
  /** Constructs a new rebalanced view the original LP.
   */
  NewType operator()( OriginalType const & origProblem );
};

} // namespace Tpetra

#endif // TPETRA_REBALANCE_LINEARPROBLEM_DECL_HPP
