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

namespace Tpetra {

template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
class BlockMap : public Teuchos::Describable {
 public:
  //! @name Constructor/Destructor Methods
  //@{

  /*! \brief BlockMap constructor specifying numGlobalBlocks and constant blockSize.
   */
  BlockMap(global_size_t numGlobalBlocks, LocalOrdinal blockSize,
      GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
      const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());

  /*! \brief BlockMap constructor specifying num global and local blocks, and constant blockSize.
   */
  BlockMap(global_size_t numGlobalBlocks, size_t numLocalBlocks,
      LocalOrdinal blockSize, GlobalOrdinal indexBase,
      const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
      const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());

  /*! \brief BlockMap constructor specifying numGlobalBlocks and lists of local blocks and blockSizes.
   */
  BlockMap(global_size_t numGlobalBlocks,
      const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalBlockIDs,
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
  BlockMap(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& pointMap, const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalBlockIDs, const Teuchos::ArrayView<const LocalOrdinal>& myBlockSizes);

  //! BlockMap destructor.
  ~BlockMap(){}

  //@}

  //! @name Attribute Accessor Methods
  //@{

  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& getPointMap() const
    { return pointMap_; }

  size_t getNodeNumBlocks() const;

  Teuchos::ArrayView<const GlobalOrdinal> getBlockIDs() const;

  Teuchos::ArrayView<const LocalOrdinal> getBlockSizes() const;

  Teuchos::ArrayView<const LocalOrdinal> getFirstPointInBlocks() const;

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
  LocalOrdinal getBlockSize(LocalOrdinal localBlockID) const;

  //! Return the first local point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  LocalOrdinal getFirstLocalPointInBlock(LocalOrdinal localBlockID) const;

  //! Return the first local point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  GlobalOrdinal getFirstGlobalPointInBlock(LocalOrdinal localBlockID) const;

  //@}

 private:
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > pointMap_;
  Teuchos::Array<GlobalOrdinal> myGlobalBlockIDs_;
  Teuchos::Array<LocalOrdinal> blockSizes_;
  Teuchos::Array<LocalOrdinal> firstPointInBlock_;
  bool blockIDsAreContiguous_;
  LocalOrdinal constantBlockSize_;
};//class BlockMap
}//namespace Tpetra

#endif

