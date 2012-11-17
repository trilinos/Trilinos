// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_TPETRABLOCKMAP_HPP
#define XPETRA_TPETRABLOCKMAP_HPP

#include "Xpetra_ConfigDefs.hpp"

#ifndef HAVE_XPETRA_TPETRA
#error This file should be included only if HAVE_XPETRA_TPETRA is defined.
#endif

#include "Xpetra_Map.hpp"
#include "Xpetra_BlockMap.hpp"

#include "Tpetra_BlockMap.hpp"
#include "Tpetra_Map.hpp"

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
class TpetraBlockMap
  : public BlockMap<LocalOrdinal,GlobalOrdinal,Node>
{
 public:

  // The following typedef are used by the XPETRA_DYNAMIC_CAST() macro.
  typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass;

  //! @name Constructor/Destructor Methods
  //@{

  /*! \brief TpetraBlockMap constructor specifying numGlobalBlocks and constant blockSize.
   */
  TpetraBlockMap(global_size_t numGlobalBlocks,
                 LocalOrdinal blockSize,
                 GlobalOrdinal indexBase,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                 const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) : map_(rcp(new Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalBlocks, blockSize, indexBase, comm, node))) {  }

  /*! \brief TpetraBlockMap constructor specifying num global and local blocks, and constant blockSize.
   */
  TpetraBlockMap(global_size_t numGlobalBlocks,
                 size_t numLocalBlocks,
                 LocalOrdinal blockSize,
                 GlobalOrdinal indexBase,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                 const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) : map_(rcp(new Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalBlocks, numLocalBlocks, blockSize, indexBase, comm, node))) {  }

  /*! \brief TpetraBlockMap constructor specifying numGlobalBlocks and lists of local blocks first-global-point-in-blocks, and blockSizes.
   */
  TpetraBlockMap(global_size_t numGlobalBlocks,
                 const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalBlockIDs,
                 const Teuchos::ArrayView<const GlobalOrdinal>& myFirstGlobalPointInBlocks,
                 const Teuchos::ArrayView<const LocalOrdinal>& myBlockSizes,
                 GlobalOrdinal indexBase,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                 const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) : map_(rcp(new Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalBlocks, myGlobalBlockIDs, myFirstGlobalPointInBlocks, myBlockSizes, indexBase, comm, node))) {  }

  /*! \brief TpetraBlockMap constructor which takes a "regular" Map.
   * The arrays myGlobalBlockIDs and myBlockSizes must be the same length, and
   * sum(myBlockSizes) must equal pointMap->getNodeNumElements().
   * If these arrays are different lengths or sum(myBlockSizes) is incorrect,
   * then std::runtime_error is thrown.
   */
  TpetraBlockMap(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& pointMap,
                 const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalBlockIDs,
                 const Teuchos::ArrayView<const LocalOrdinal>& myBlockSizes,
                 const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) {

    XPETRA_RCP_DYNAMIC_CAST(const TpetraMapClass, pointMap, tPointMap, "Xpetra::TpetraBlockMap constructors only accept Xpetra::TpetraMap as input arguments.");
    map_ = rcp(new Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node>(tPointMap->getTpetra_Map(), myGlobalBlockIDs, myBlockSizes, node));
  }

  TpetraBlockMap(const Teuchos::RCP<const Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node> > &map) : map_(map) {  }

  //! TpetraBlockMap destructor.
  virtual ~TpetraBlockMap(){  }

  //@}

  //! @name Attribute Accessor Methods
  //@{

  global_size_t getGlobalNumBlocks() const {  return map_->getGlobalNumBlocks(); }

  //! Return number of blocks on the local processor.
  size_t getNodeNumBlocks() const {  return map_->getNodeNumBlocks(); }

  Teuchos::ArrayView<const GlobalOrdinal> getNodeBlockIDs() const {  return map_->getNodeBlockIDs(); }

  bool isBlockSizeConstant() const {  return map_->isBlockSizeConstant(); }

  //! Return ArrayRCP of first-local-point in local blocks.
  Teuchos::ArrayRCP<const LocalOrdinal> getNodeFirstPointInBlocks() const {  return map_->getNodeFirstPointInBlocks(); }

  //! Return device-resident ArrayRCP of first-local-point in local blocks.
  /*! This version of this method is primarily used internally by VbrMatrix
      for passing data to the matrix-vector-product kernel.
  */
  Teuchos::ArrayRCP<const LocalOrdinal> getNodeFirstPointInBlocks_Device() const {  return map_->getNodeFirstPointInBlocks_Device(); }

  //! Return the globalBlockID corresponding to the given localBlockID
  /*! If localBlockID is not present on this processor, returns Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
  */
  GlobalOrdinal getGlobalBlockID(LocalOrdinal localBlockID) const {  return map_->getGlobalBlockID(localBlockID); }

  //! Return the localBlockID corresponding to the given globalBlockID
  /*! If globalBlockID is not present on this processor, returns Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
  */
  LocalOrdinal getLocalBlockID(GlobalOrdinal globalBlockID) const {  return map_->getLocalBlockID(globalBlockID); }

  //! Return the block-size for localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  LocalOrdinal getLocalBlockSize(LocalOrdinal localBlockID) const {  return map_->getLocalBlockSize(localBlockID); }

  //! Return the first local point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  LocalOrdinal getFirstLocalPointInLocalBlock(LocalOrdinal localBlockID) const {  return map_->getFirstLocalPointInLocalBlock(localBlockID); }

  //! Return the first global point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  GlobalOrdinal getFirstGlobalPointInLocalBlock(LocalOrdinal localBlockID) const {  return map_->getFirstGlobalPointInLocalBlock(localBlockID); }

  //@}

  RCP< const Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node> > getTpetra_BlockMap() const {  return map_; }

private:
  RCP< const Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node> > map_;


};//class TpetraBlockMap

}//namespace Xpetra

#define XPETRA_TPETRABLOCKMAP_SHORT
#endif
