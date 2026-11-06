// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CREATECRSGRAPHFROMSTENCIL_DECL_HPP
#define TPETRA_CREATECRSGRAPHFROMSTENCIL_DECL_HPP

#include "Tpetra_CrsGraph_decl.hpp"
#include "Tpetra_Map_fwd.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_RCP.hpp"

namespace Tpetra {
namespace Experimental {

/// \brief A relative Cartesian grid position in at most three dimensions.
///
/// Components belonging to dimensions not present in the grid must be zero.
struct StencilOffset {
  int x = 0;
  int y = 0;
  int z = 0;
};

enum class StencilNumbering {
  globalLexicographic,
  rankBlocked
};

/// \brief Create a fill-complete graph directly in device-resident CSR form.
/// This implementation supports radius-two stencils, at most three
/// dimensions, and at most three degrees of freedom per grid point.
template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>
createCrsGraphFromStencil(
    const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
    const Teuchos::ArrayView<const LocalOrdinal>& localGridDimensions,
    const Teuchos::ArrayView<const StencilOffset>& stencilOffsets,
    const LocalOrdinal dofsPerNode, const StencilNumbering numbering,
    const bool debug);

}  // namespace Experimental
}  // namespace Tpetra

#endif  // TPETRA_CREATECRSGRAPHFROMSTENCIL_DECL_HPP
