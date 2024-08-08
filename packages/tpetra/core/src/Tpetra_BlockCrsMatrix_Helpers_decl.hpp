// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKCRSMATRIX_HELPERS_DECL_HPP
#define TPETRA_BLOCKCRSMATRIX_HELPERS_DECL_HPP

/// \file Tpetra_BlockCrsMatrix_Helpers_decl.hpp

#include "Tpetra_BlockCrsMatrix_fwd.hpp"
#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Tpetra_CrsGraph_fwd.hpp"
#include "Tpetra_Map_fwd.hpp"
#include "Teuchos_RCP.hpp"
#include <string>
#include <ostream>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  class ParameterList; // forward declaration
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {

  /// \brief Helper function to write a BlockCrsMatrix.  Calls the 3-argument version.
  template<class Scalar, class LO, class GO, class Node>
  void blockCrsMatrixWriter(BlockCrsMatrix<Scalar,LO,GO,Node> const &A, std::string const &fileName);

  /// \brief Helper function to write a BlockCrsMatrix.  Calls the 3-argument version.
  template<class Scalar, class LO, class GO, class Node>
  void blockCrsMatrixWriter(BlockCrsMatrix<Scalar,LO,GO,Node> const &A, std::string const &fileName, Teuchos::ParameterList const &params);

  /// \brief Helper function to write a BlockCrsMatrix.  Calls the 3-argument version.
  template<class Scalar, class LO, class GO, class Node>
  void blockCrsMatrixWriter(BlockCrsMatrix<Scalar,LO,GO,Node> const &A, std::ostream &os);

  /*! \brief Helper function to write a BlockCrsMatrix.

    Writes the block matrix to the specified ostream in point form.  The following parameter list options are available:

    - "always use parallel algorithm" : on one process, this forces the use of the parallel strip-mining algorithm (default=false)
    - "print MatrixMarket header"     : if false, don't print the MatrixMarket header (default=true)
    - "precision"                     : precision to be used in printing matrix entries (default=C++ default)
    - "zero-based indexing"           : if true, print the matrix with 0-based indexing (default=false)
  */
  template<class Scalar, class LO, class GO, class Node>
  void blockCrsMatrixWriter(BlockCrsMatrix<Scalar,LO,GO,Node> const &A, std::ostream &os, Teuchos::ParameterList const &params);

  /*! @brief Helper function called by blockCrsMatrixWriter.

  This function should not be called directly.
  */
  template<class Scalar, class LO, class GO, class Node>
  void writeMatrixStrip(BlockCrsMatrix<Scalar,LO,GO,Node> const &A, std::ostream &os, Teuchos::ParameterList const &params);

  /// \brief Non-member constructor that creates the CrsGraph of a BlockCrsMatrix
  /// from an existing point CrsMatrix and a block size.
  ///
  /// This function accepts an already constructed point version of the block matrix.
  /// Assumptions:
  ///   - All point entries in a logical block must be stored in the CrsMatrix, even
  ///     if the values are zero.
  ///   - Point rows corresponding to a particular mesh node must be stored consecutively.
  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<Tpetra::CrsGraph<LO, GO, Node>>
  getBlockCrsGraph(const Tpetra::CrsMatrix<Scalar, LO, GO, Node>& pointMatrix, const LO &blockSize, bool use_local_ID=true);

  /// \brief Non-member constructor that creates a BlockCrsMatrix from an existing point CrsMatrix.
  ///
  /// This function accepts an already constructed point version of the block matrix.
  /// Assumptions:
  ///   - All point entries in a logical block must be stored in the CrsMatrix, even
  ///     if the values are zero.
  ///   - Point rows corresponding to a particular mesh node must be stored consecutively.
  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<BlockCrsMatrix<Scalar, LO, GO, Node>>
  convertToBlockCrsMatrix(const Tpetra::CrsMatrix<Scalar, LO, GO, Node>& pointMatrix, const LO &blockSize, bool use_local_ID=true);

  /// \brief Fill all point entries in a logical block of a CrsMatrix with zeroes. This should be called
  /// before convertToBlockCrsMatrix if your Crs is not filled. Performed on host.
  ///
  /// This function accepts an already constructed point version of the block matrix.
  /// Assumptions:
  ///   - Point rows corresponding to a particular mesh node must be stored consecutively.
  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>>
  fillLogicalBlocks(const Tpetra::CrsMatrix<Scalar, LO, GO, Node>& pointMatrix, const LO &blockSize);

  /// \brief Unfill all point entries in a logical block of a CrsMatrix with zeroes. This can be called
  /// after convertToCrsMatrix if you don't want to keep the filled values. Performed on host.
  ///
  /// This function accepts an already constructed point version of the block matrix. It is expected
  /// that, given a crs A
  /// Afill = fillLogicalBlocks(A)
  /// Abcrs = convertToBlockCrsMatrix(Afill)
  /// Acrs  = convertToCrsMatrix(Abcrs);
  /// Aorig = unfillFormerBlockCrs(Acrs);
  /// Expect A and Aorig are identical
  /// Assumptions:
  ///   - Point rows corresponding to a particular mesh node must be stored consecutively.
  ///   - Assumes all zero entries are filled values that should be removed.
  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>>
  unfillFormerBlockCrs(const Tpetra::CrsMatrix<Scalar, LO, GO, Node>& pointMatrix);

  /// @brief Helper function to generate a mesh map from a point map.
  /// Important! It's assumed that point GIDs associated with a single mesh GID appear consecutively in pointMap.
  template<class LO, class GO, class Node>
  Teuchos::RCP<const Tpetra::Map<LO,GO,Node>>
  createMeshMap(LO const &blockSize, const Tpetra::Map<LO,GO,Node> &pointMap, bool use_local_ID=false);

  /// \brief Non-member constructor that creates a point CrsMatrix from an existing BlockCrsMatrix.
  ///
  /// This function accepts an already constructed block matrix.
  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<CrsMatrix<Scalar, LO, GO, Node>>
  convertToCrsMatrix(const Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node>& blockMatrix);

  /// @brief Helper function to generate a point map from a block map (with a given block size)
  /// GIDs associated with a single mesh GID appear consecutively in the output map
  template<class LO, class GO, class Node>
  Teuchos::RCP<const Tpetra::Map<LO,GO,Node>>
  createPointMap(LO const &blockSize, const Tpetra::Map<LO,GO,Node> &blockMap);


} // namespace Tpetra

#endif // TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_HELPERS_DECL_HPP
