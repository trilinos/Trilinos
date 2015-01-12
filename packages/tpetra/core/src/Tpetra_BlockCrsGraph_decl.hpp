// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKCRSGRAPH_DECL_HPP
#define TPETRA_BLOCKCRSGRAPH_DECL_HPP

/// \file Tpetra_BlockCrsGraph_decl.hpp
/// \brief Declarations for the class Tpetra::BlockCrsGraph.

#include <Tpetra_ConfigDefs.hpp>

#ifndef HAVE_TPETRA_CLASSIC_VBR
#  error "It is an error to include this file if VBR (variable-block-size) sparse matrix support is disabled in Tpetra.  If you would like to enable VBR support, please reconfigure Trilinos with the CMake option Tpetra_ENABLE_CLASSIC_VBR set to ON, and rebuild Trilinos."
#else

#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_BlockMap.hpp>

namespace Tpetra {

/// \class BlockCrsGraph
/// \brief Block-entry counterpart to Tpetra::CrsGraph.
///
/// BlockCrsGraph doesn't inherit from Tpetra::CrsGraph, but always
/// holds a Tpetra::CrsGraph as an instance variable.
///
/// The reason BlockCrsGraph exists is to create and hold the block
/// versions (Tpetra::BlockMap) of the Tpetra::Map objects that
/// CrsGraph holds.
///
/// BlockCrsGraph is used by Tpetra::VbrMatrix (variable block row matrix).
///
/// \warning This class is DEPRECATED.  There are known outstanding
///   bugs with the current implementations of variable-block-size
///   sparse matrices and related classes in Tpetra.
template <class LocalOrdinal = Details::DefaultTypes::local_ordinal_type,
          class GlobalOrdinal = Details::DefaultTypes::global_ordinal_type,
          class Node = Details::DefaultTypes::node_type>
class TPETRA_DEPRECATED BlockCrsGraph : public Teuchos::Describable {
public:
  typedef LocalOrdinal  local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node          node_type;
  typedef BlockMap<LocalOrdinal, GlobalOrdinal, Node> block_map_type;

  //! @name Constructor/Destructor Methods
  //@{

  /*! \brief BlockCrsGraph constructor specifying a block-row-map and max-num-block-entries-per-row.
   */
  BlockCrsGraph (const Teuchos::RCP<const block_map_type>& blkRowMap,
                 size_t maxNumEntriesPerRow,
                 ProfileType pftype = DynamicProfile);

  //! BlockCrsGraph destructor.
  ~BlockCrsGraph ()
  {}

  //@}
  //! @name Insertion/Extraction Methods
  //@{

  //! Submit graph indices, using global IDs.
  void insertGlobalIndices (GlobalOrdinal row, const Teuchos::ArrayView<const GlobalOrdinal>& indices);

  //! Get row-offsets. (This is the bptr array in VBR terminology.)
  /*! Returns null if optimizeStorage has not been called.
   */
  Teuchos::ArrayRCP<const size_t> getNodeRowOffsets() const;

  //! Get packed-col-indices. (This is the bindx array in VBR terminology.)
  /*! Returns null if optimizeStorage has not been called.
   */
  Teuchos::ArrayRCP<const LocalOrdinal> getNodePackedIndices() const;

  //@}
  //! @name Transformational Methods
  //@{

  //! Communicate non-local contributions to other nodes.
  void globalAssemble ();

  /*! \brief Signal that data entry is complete, specifying domain and range maps.
      Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
      If \c OptimizeStorage is true, then optimizeStorage() is called as well.
   */
  void
  fillComplete (const Teuchos::RCP<const block_map_type>& blkDomainMap,
                const Teuchos::RCP<const block_map_type>& blkRangeMap,
                OptimizeOption os = DoOptimizeStorage);

  /*! \brief Signal that data entry is complete.
      Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
      If \c OptimizeStorage is true, then optimizeStorage() is called as well.
      \note This method calls fillComplete( getRowMap(), getRowMap(), OptimizeStorage ).
   */
  void fillComplete (OptimizeOption os = DoOptimizeStorage);

  //! \brief Re-allocate the data into contiguous storage.
  void optimizeStorage ();

  //@}

  //! @name Attribute Accessor Methods
  //@{

  //! Whether fillComplete () has been called.
  bool isFillComplete () const;

  //! Whether the column indices are stored as local indices.
  bool isLocallyIndexed() const;

  /// \brief Whether the graph is upper triangular on the calling process.
  ///
  /// Upper or lower triangularity is a local (per-process) attribute.
  bool isUpperTriangular() const;

  /// \brief Whether the graph is upper triangular on the calling process.
  ///
  /// Upper or lower triangularity is a local (per-process) attribute.
  bool isLowerTriangular() const;

  //! The number of block rows owned by the calling process.
  size_t getNodeNumBlockRows() const;

  //! The global number of block rows.
  size_t getGlobalNumBlockRows() const;

  //! The number of (block) diagonal entries owned by the calling process.
  size_t getNodeNumBlockDiags() const;

  //! The number of (block) entries owned by the calling process.
  size_t getNodeNumBlockEntries() const;

  //! Returns the number of block-columns in the specified global block row.
  size_t getGlobalBlockRowLength(GlobalOrdinal row) const;

  //! Returns a read-only view of the block-column-indices for the specified global block row.
  void
  getGlobalBlockRowView (GlobalOrdinal row,
                         Teuchos::ArrayView<const GlobalOrdinal>& blockCols) const;

  //! Returns a read-only view of the block-column-indices for the specified local block row.
  void
  getLocalBlockRowView (LocalOrdinal row,
                        Teuchos::ArrayView<const LocalOrdinal>& blockCols) const;

  //! Returns the block-row map.
  Teuchos::RCP<const block_map_type> getBlockRowMap () const;

  //! Returns the block-column map.
  Teuchos::RCP<const block_map_type> getBlockColMap () const;

  //! Returns the block-domain map.
  Teuchos::RCP<const block_map_type> getBlockDomainMap () const;

  //! Returns the block-range map.
  Teuchos::RCP<const block_map_type> getBlockRangeMap () const;

  //@}

private:
  typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node> crs_graph_type;
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> point_map_type;

  Teuchos::RCP<crs_graph_type> ptGraph_;
  Teuchos::RCP<const block_map_type> blkRowMap_;
  Teuchos::RCP<const block_map_type> blkColMap_;
  Teuchos::RCP<const block_map_type> blkDomainMap_;
  Teuchos::RCP<const block_map_type> blkRangeMap_;
};
} // namespace Tpetra

#endif // ! HAVE_TPETRA_CLASSIC_VBR
#endif // ! TPETRA_BLOCKCRSGRAPH_DECL_HPP

