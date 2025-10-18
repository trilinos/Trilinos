// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TRILINOSCOUPLINGS_REBALANCE_LINEARPROBLEM_DECL_HPP
#define TRILINOSCOUPLINGS_REBALANCE_LINEARPROBLEM_DECL_HPP

/// \file TrilinosCouplings_Rebalance_LinearProblem_decl.hpp
/// \brief Declaration of the TrilinosCouplings::Rebalance_LinearProblem class

#include <Tpetra_Transform.hpp>
#include <Tpetra_LinearProblem.hpp>

namespace TrilinosCouplings {

///
/** Given and input Tpetra LinearProblem, a "rebalanced" version will be returned.
 *  The data in the new T_LP is a "rebalanced" view of the original.
 */
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class Rebalance_LinearProblem : public Tpetra::ViewTransform< Tpetra::LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
{
public:
using NewType      = typename Tpetra::ViewTransform< Tpetra::LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::NewType;
using OriginalType = typename Tpetra::ViewTransform< Tpetra::LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::OriginalType;

  ///
  /** Constructor
   */
  Rebalance_LinearProblem( Teuchos::RCP< Teuchos::ParameterList > paramListForZoltan2PartitioningProblem );

  ///
  /** Destructor
   */
  ~Rebalance_LinearProblem();

  ///
  /** Constructs a new rebalanced view the original LP.
   */
  NewType operator()( OriginalType const & origProblem );

private:
  Teuchos::RCP< Teuchos::ParameterList > paramListForZoltan2PartitioningProblem_;
};

} // namespace TrilinosCouplings

#endif // TRILINOSCOUPLINGS_REBALANCE_LINEARPROBLEM_DECL_HPP
