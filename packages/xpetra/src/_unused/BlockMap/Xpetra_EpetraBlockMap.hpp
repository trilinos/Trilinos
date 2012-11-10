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
#ifndef XPETRA_EPETRABLOCKMAP_HPP
#define XPETRA_EPETRABLOCKMAP_HPP

#include "Xpetra_ConfigDefs.hpp"

#ifndef HAVE_XPETRA_EPETRA
#error This file should be included only if HAVE_XPETRA_EPETRA is defined.
#endif

#include "Xpetra_Map.hpp"
#include "Xpetra_BlockMap.hpp"

#include "Epetra_BlockMap.h"
#include "Xpetra_EpetraUtils.hpp"
#include "Xpetra_EpetraExceptions.hpp"

namespace Xpetra {

/*!
  @class EpetraBlockMap
  @brief Block-entry counterpart to Xpetra::Map.

  BlockMap doesn't inherit Xpetra::Map
*/
class EpetraBlockMap : public Xpetra::BlockMap<int,int> {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;

 public:

  //! @name Constructor/Destructor Methods
  //@{

  /*! \brief EpetraBlockMap constructor specifying numGlobalBlocks and constant blockSize.
   */
  EpetraBlockMap(global_size_t numGlobalBlocks,
                 int blockSize,
                 int indexBase,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                 const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) : map_(rcp(new Epetra_BlockMap(numGlobalBlocks, blockSize, indexBase, *Teuchos2Epetra_Comm(comm)))) {  }

  //TODO
  //CATCH_EPETRA_EXCEPTION_AND_THROW_INVALID_ARG((map_ = (rcp(new Epetra_Map(numGlobalElements, indexBase, *Teuchos2Epetra_Comm(comm))))););

  /*! \brief EpetraBlockMap constructor specifying num global and local blocks, and constant blockSize.
   */
  EpetraBlockMap(global_size_t numGlobalBlocks,
                 size_t numLocalBlocks,
                 int blockSize,
                 int indexBase,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                 const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) : map_(rcp(new Epetra_BlockMap(numGlobalBlocks, numLocalBlocks, blockSize, indexBase, *Teuchos2Epetra_Comm(comm)))) {  }

  /*! \brief EpetraBlockMap constructor specifying numGlobalBlocks and lists of local blocks first-global-point-in-blocks, and blockSizes.
   */
 // EpetraBlockMap(global_size_t numGlobalBlocks,
//                  const Teuchos::ArrayView<const int>& myGlobalBlockIDs,
//                  const Teuchos::ArrayView<const int>& myFirstGlobalPointInBlocks,
//                  const Teuchos::ArrayView<const int>& myBlockSizes,
//                  int indexBase,
//                  const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
//                  const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) {
//
//     map_ = rcp(new Epetra_BlockMap(numGlobalBlocks, myGlobalBlockIDs, myFirstGlobalPointInBlocks, myBlockSizes, indexBase, *Teuchos2Epetra_Comm(comm)));
//   }

  /*! \brief EpetraBlockMap constructor which takes a "regular" Map.
   * The arrays myGlobalBlockIDs and myBlockSizes must be the same length, and
   * sum(myBlockSizes) must equal pointMap->getNodeNumElements().
   * If these arrays are different lengths or sum(myBlockSizes) is incorrect,
   * then std::runtime_error is thrown.
   */
//   EpetraBlockMap(const Teuchos::RCP<const Map<int,int> >& pointMap,
//                  const Teuchos::ArrayView<const int>& myGlobalBlockIDs,
//                  const Teuchos::ArrayView<const int>& myBlockSizes,
//                  const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) {
//
//     map_ = rcp(new Epetra_BlockMap(pointMap, myGlobalBlockIDs, myBlockSizes));
//   }

  EpetraBlockMap(const Teuchos::RCP<const Epetra_BlockMap > &map) : map_(map) {  }

  //! EpetraBlockMap destructor.
  virtual ~EpetraBlockMap(){  }

  //@}

  //! @name Attribute Accessor Methods
  //@{
#ifdef XPETRA_NOT_IMPLEMENTED
  const Teuchos::RCP<const Map<int,int> >& getPointMap() const {  return map_->getPointMap(); }
#endif

  global_size_t getGlobalNumBlocks() const {  return map_->NumGlobalElements(); }

  //! Return number of blocks on the local processor.
  size_t getNodeNumBlocks() const {  return map_->NumMyElements(); }

  //TODO Teuchos::ArrayView<const int> getNodeBlockIDs() const {  return map_->getNodeBlockIDs(); }

  bool isBlockSizeConstant() const {  return map_->ConstantElementSize(); }

  //! Return ArrayRCP of first-local-point in local blocks.
  //TODO Teuchos::ArrayRCP<const int> getNodeFirstPointInBlocks() const {  return map_->FirstPointInElementList(); }

#ifdef XPETRA_NOT_IMPLEMENTED_FOR_EPETRA
  //! Return device-resident ArrayRCP of first-local-point in local blocks.
  /*! This version of this method is primarily used internally by VbrMatrix
      for passing data to the matrix-vector-product kernel.
  */
  Teuchos::ArrayRCP<const int> getNodeFirstPointInBlocks_Device() const {  return map_->getNodeFirstPointInBlocks_Device(); }
#endif

  //! Return the globalBlockID corresponding to the given localBlockID
  /*! If localBlockID is not present on this processor, returns Teuchos::OrdinalTraits<int>::invalid().
  */
  int getGlobalBlockID(int localBlockID) const {  return map_->GID(localBlockID); }

  //! Return the localBlockID corresponding to the given globalBlockID
  /*! If globalBlockID is not present on this processor, returns Teuchos::OrdinalTraits<int>::invalid().
  */
  int getLocalBlockID(int globalBlockID) const {  return map_->LID(globalBlockID); }

  //! Return the block-size for localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  int getLocalBlockSize(int localBlockID) const {  return map_->ElementSize(localBlockID); }

  //! Return the first local point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  int getFirstLocalPointInLocalBlock(int localBlockID) const {  return map_->FirstPointInElement(localBlockID); }

#ifdef XPETRA_NOT_IMPLEMENTED_FOR_EPETRA
  //! Return the first global point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown. //TODO: throw !
   */
  int getFirstGlobalPointInLocalBlock(int localBlockID) const {  return; }
#endif

  //@}

  RCP< const Epetra_BlockMap > getEpetra_BlockMap() const {  return map_; }

private:
  const RCP< const Epetra_BlockMap > map_;


};//class EpetraBlockMap

// //-----------------------------------------------------------------
// template<class int,class int,class Node>
// Teuchos::RCP<const Xpetra::Map<int,int> >
// convertEpetraBlockMapToEpetraPointMap(const Teuchos::RCP<const Xpetra::EpetraBlockMap<int,int> >& blockMap) {
//   return rcp(new EpetraMap(convertEpetraBlockMapToEpetraPointMap(blockMap.getEpetra_BlockMap())));
// }

}//namespace Xpetra

#endif
