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

#ifndef TPETRA_BLOCKMAP_DECL_HPP
#define TPETRA_BLOCKMAP_DECL_HPP

#include "Tpetra_Map.hpp"

/** \file Tpetra_BlockMap_decl.hpp

  Declarations for the class Tpetra::BlockMap.
*/
namespace Tpetra {

/** \brief Block-entry counterpart to Tpetra::Map.

  BlockMap doesn't inherit Tpetra::Map, but always holds a Tpetra::Map as
  a class-member attribute.
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

  /*! \brief BlockMap constructor which takes a "regular" Map.
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

  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& getPointMap() const
    { return pointMap_; }

  global_size_t getGlobalNumBlocks() const;

  //! Return number of blocks on the local processor.
  size_t getNodeNumBlocks() const;

  Teuchos::ArrayView<const GlobalOrdinal> getNodeBlockIDs() const;

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

  //@}

 private:
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > pointMap_;
  global_size_t globalNumBlocks_;
  Teuchos::Array<GlobalOrdinal> myGlobalBlockIDs_;
  Teuchos::ArrayRCP<LocalOrdinal> pbuf_firstPointInBlock_;
  Teuchos::ArrayRCP<const LocalOrdinal> view_firstPointInBlock_;
  bool blockIDsAreContiguous_;
  LocalOrdinal constantBlockSize_;
};//class BlockMap
}//namespace Tpetra

#endif

