// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_BLOCKCRSGRAPH_DECL_HPP
#define TPETRA_BLOCKCRSGRAPH_DECL_HPP

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_BlockMap.hpp"

/** \file Tpetra_BlockCrsGraph_decl.hpp

  Declarations for the class Tpetra::BlockCrsGraph.
*/
namespace Tpetra {

/** \brief Block-entry counterpart to Tpetra::CrsGraph.

  BlockCrsGraph doesn't inherit Tpetra::CrsGraph, but always holds a
  Tpetra::CrsGraph as a class-member attribute.

  The reason BlockCrsGraph exists is to create and hold the block-versions
  (Tpetra::BlockMap) of the Tpetra::Map objects that CrsGraph holds.

  BlockCrsGraph is used by Tpetra::VbrMatrix (variable block row matrix).
*/
template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
class BlockCrsGraph : public Teuchos::Describable {
 public:
  typedef LocalOrdinal  local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node          node_type;

  //! @name Constructor/Destructor Methods
  //@{

  /*! \brief BlockCrsGraph constructor specifying a block-row-map and max-num-block-entries-per-row.
   */
  BlockCrsGraph(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blkRowMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

  //! BlockCrsGraph destructor.
  ~BlockCrsGraph(){}

  //@}

  //! @name Insertion/Extraction Methods
  //@{

  //! Submit graph indices, using global IDs.
  void insertGlobalIndices(GlobalOrdinal row, const Teuchos::ArrayView<const GlobalOrdinal> &indices);

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

  //! \brief Communicate non-local contributions to other nodes.
  void globalAssemble();

  /*! \brief Signal that data entry is complete, specifying domain and range maps. 
      Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
      If \c OptimizeStorage is true, then optimizeStorage() is called as well.
   */
  void fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &blkDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &blkRangeMap, OptimizeOption os = DoOptimizeStorage);

  /*! \brief Signal that data entry is complete. 
      Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
      If \c OptimizeStorage is true, then optimizeStorage() is called as well.
      \note This method calls fillComplete( getRowMap(), getRowMap(), OptimizeStorage ).
   */
  void fillComplete(OptimizeOption os = DoOptimizeStorage);

  //! \brief Re-allocate the data into contiguous storage.
  void optimizeStorage();

  //@}

  //! @name Attribute Accessor Methods
  //@{

  //! Returns \c true if fillComplete() has been called.
  bool isFillComplete() const;

  //! \brief If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */
  bool isLocallyIndexed() const;

  //! \brief true if graph is upper-triangular.
  bool isUpperTriangular() const;

  //! \brief true if graph is lower-triangular.
  bool isLowerTriangular() const;

  //! Returns the number of block rows owned on the calling node.
  size_t getNodeNumBlockRows() const;

  //! Returns the global number of block rows.
  size_t getGlobalNumBlockRows() const;

  //! Returns the number of diagonal entries on the calling node.
  size_t getNodeNumBlockDiags() const;

  //! Returns the local number of entries in the graph.
  size_t getNodeNumBlockEntries() const;

  //! Returns the number of block-columns in the specified global block row.
  size_t getGlobalBlockRowLength(GlobalOrdinal row) const;

  //! Returns a read-only view of the block-column-indices for the specified global block row.
  void getGlobalBlockRowView(GlobalOrdinal row,
                             Teuchos::ArrayView<const GlobalOrdinal>& blockCols) const;

  //! Returns a read-only view of the block-column-indices for the specified local block row.
  void getLocalBlockRowView(LocalOrdinal row,
                             Teuchos::ArrayView<const LocalOrdinal>& blockCols) const;

  //! Returns the block-row map.
  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockRowMap() const;

  //! Returns the block-column map.
  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockColMap() const;

  //! Returns the block-domain map.
  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockDomainMap() const;

  //! Returns the block-range map.
  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > & getBlockRangeMap() const;

  //@}

 private:
  Teuchos::RCP<CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > ptGraph_;

  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blkRowMap_;
  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blkColMap_;
  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blkDomainMap_;
  Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > blkRangeMap_;
};//class BlockCrsGraph
}//namespace Tpetra

#endif

