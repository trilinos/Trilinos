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

#ifndef TPETRA_BLOCKMAP_DECL_HPP
#define TPETRA_BLOCKMAP_DECL_HPP

#include <map>

#include "Tpetra_Map.hpp"

/** \file Tpetra_BlockMap_decl.hpp

  Declarations for the class Tpetra::BlockMap.
*/
namespace Tpetra {

/** \brief Block-entry counterpart to Tpetra::Map.

  BlockMap doesn't inherit Tpetra::Map, but always holds a Tpetra::Map as
  a class-member attribute.

  Tpetra::BlockMap essentially holds information about how the point-entries
  in Tpetra::Map are grouped together in blocks. A block-entry consists of
  1 or more point-entries.

  Example usage: If a solution-space consists of multiple degrees-of-freedom
  at each finite-element node in a mesh, such as a displacement vector, it
  might be described as having a block of size 3 (in 3D) at each mesh node.
  Thus for a mesh with N nodes, the point-entry map will have N*3 entries,
  whereas the block-map will have N blocks, each of size 3.
*/
template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
class BlockMap : public Teuchos::Describable {
 public:
  typedef LocalOrdinal  local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node          node_type;

  //! @name Constructor/Destructor Methods
  //@{

  /*! \brief BlockMap constructor specifying numGlobalBlocks and constant blockSize.
   */
  BlockMap(global_size_t numGlobalBlocks,
           LocalOrdinal blockSize,
           GlobalOrdinal indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());

  /*! \brief BlockMap constructor specifying num global and local blocks, and constant blockSize.
   */
  BlockMap(global_size_t numGlobalBlocks,
           size_t numLocalBlocks,
           LocalOrdinal blockSize,
           GlobalOrdinal indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());

  /*! \brief BlockMap constructor specifying numGlobalBlocks and lists of local blocks first-global-point-in-blocks, and blockSizes.
   */
  BlockMap(global_size_t numGlobalBlocks,
      const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalBlockIDs,
      const Teuchos::ArrayView<const GlobalOrdinal>& myFirstGlobalPointInBlocks,
      const Teuchos::ArrayView<const LocalOrdinal>& myBlockSizes,
      GlobalOrdinal indexBase,
      const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
      const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());

  /*! \brief BlockMap constructor which takes a point-entry Map.
   * The arrays myGlobalBlockIDs and myBlockSizes must be the same length, and
   * sum(myBlockSizes) must equal pointMap->getNodeNumElements().
   * If these arrays are different lengths or sum(myBlockSizes) is incorrect,
   * then std::runtime_error is thrown.
   */
  BlockMap(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& pointMap,
           const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalBlockIDs,
           const Teuchos::ArrayView<const LocalOrdinal>& myBlockSizes,
           const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());

  //! BlockMap destructor.
  ~BlockMap(){}

  //@}

  //! @name Attribute Accessor Methods
  //@{

  //! Return this block-map's point-entry map attribute.
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& getPointMap() const
    { return pointMap_; }

  //! Return global number of blocks.
  global_size_t getGlobalNumBlocks() const;

  //! Return number of blocks on the local processor.
  size_t getNodeNumBlocks() const;

  //! Return array-view of block-ids for this local processor.
  Teuchos::ArrayView<const GlobalOrdinal> getNodeBlockIDs() const;

  //! Return true if all blocks have the same size.
  bool isBlockSizeConstant() const;

  //! Return ArrayRCP of first-local-point in local blocks.
  Teuchos::ArrayRCP<const LocalOrdinal> getNodeFirstPointInBlocks() const;

  //! Return device-resident ArrayRCP of first-local-point in local blocks.
  /*! This version of this method is primarily used internally by VbrMatrix
      for passing data to the matrix-vector-product kernel.
  */
  Teuchos::ArrayRCP<const LocalOrdinal> getNodeFirstPointInBlocks_Device() const;

  //! Return the globalBlockID corresponding to the given localBlockID
  /*! If localBlockID is not present on this processor, returns Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
  */
  GlobalOrdinal getGlobalBlockID(LocalOrdinal localBlockID) const;

  //! Return the localBlockID corresponding to the given globalBlockID
  /*! If globalBlockID is not present on this processor, returns Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
  */
  LocalOrdinal getLocalBlockID(GlobalOrdinal globalBlockID) const;

  //! Return the block-size for localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  LocalOrdinal getLocalBlockSize(LocalOrdinal localBlockID) const;

  //! Return the first local point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  LocalOrdinal getFirstLocalPointInLocalBlock(LocalOrdinal localBlockID) const;

  //! Return the first global point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  GlobalOrdinal getFirstGlobalPointInLocalBlock(LocalOrdinal localBlockID) const;

  //! Return the first-global-point-in-block and block-sizes for a list of block-IDs on remote processors.
  void getRemoteBlockInfo(const Teuchos::ArrayView<const GlobalOrdinal>& GBIDs,
                          const Teuchos::ArrayView<GlobalOrdinal>& firstGlobalPointInBlocks,
                          const Teuchos::ArrayView<LocalOrdinal>& blockSizes) const;
  //@}

 private:
  void setup_noncontig_mapping();

  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > pointMap_;
  global_size_t globalNumBlocks_;
  Teuchos::Array<GlobalOrdinal> myGlobalBlockIDs_;
  Teuchos::ArrayRCP<LocalOrdinal> pbuf_firstPointInBlock_;
  Teuchos::ArrayRCP<const LocalOrdinal> view_firstPointInBlock_;
  bool blockIDsAreContiguous_;
  LocalOrdinal constantBlockSize_;
  /// \brief Global-to-local index lookup table.
  ///
  /// TODO: Use Tpetra::Details::HashTable here instead.
  std::map<GlobalOrdinal,LocalOrdinal> map_global_to_local_;
};//class BlockMap

//-----------------------------------------------------------------
template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
convertBlockMapToPointMap(const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>& blockMap)
{
  global_size_t numGlobalElems = Teuchos::OrdinalTraits<global_size_t>::invalid();
  GlobalOrdinal indexBase = blockMap.getPointMap()->getIndexBase();
  const Teuchos::RCP<const Teuchos::Comm<int> >& comm = blockMap.getPointMap()->getComm();
  const Teuchos::RCP<Node>& node = blockMap.getPointMap()->getNode();

  //Create a point-entry map where each point
  //corresponds to a block in the block-map:
  return Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElems, blockMap.getNodeBlockIDs(), indexBase, comm, node));
}

}//namespace Tpetra

#endif

