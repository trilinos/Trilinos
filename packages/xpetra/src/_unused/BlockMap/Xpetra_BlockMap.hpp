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
#ifndef XPETRA_BLOCKMAP_HPP
#define XPETRA_BLOCKMAP_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Map.hpp"

namespace Xpetra {

/*!
  @class BlockMap
  @brief Block-entry counterpart to Xpetra::Map.

  Note: BlockMap doesn't inherit from Xpetra::Map
*/
template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
class BlockMap : public Teuchos::Describable {
 public:

  //! @name Constructor/Destructor Methods
  //@{

  //! BlockMap destructor.
  virtual ~BlockMap(){  }

  //@}

#ifdef XPETRA_NOT_IMPLEMENTED_FOR_EPETRA

  //! @name Attribute Accessor Methods
  //@{

  virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& getPointMap() const =0;

  virtual global_size_t getGlobalNumBlocks() const =0;

  //! Return number of blocks on the local processor.
  virtual size_t getNodeNumBlocks() const =0;

  virtual Teuchos::ArrayView<const GlobalOrdinal> getNodeBlockIDs() const =0;

  virtual bool isBlockSizeConstant() const =0;

  //! Return ArrayRCP of first-local-point in local blocks.
  virtual Teuchos::ArrayRCP<const LocalOrdinal> getNodeFirstPointInBlocks() const =0;

#ifdef XPETRA_NOT_IMPLEMENTED_FOR_EPETRA
  //! Return device-resident ArrayRCP of first-local-point in local blocks.
  /*! This version of this method is primarily used internally by VbrMatrix
      for passing data to the matrix-vector-product kernel.
  */
  virtual Teuchos::ArrayRCP<const LocalOrdinal> getNodeFirstPointInBlocks_Device() const =0;
#endif

  //! Return the globalBlockID corresponding to the given localBlockID
  /*! If localBlockID is not present on this processor, returns Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
  */
  virtual GlobalOrdinal getGlobalBlockID(LocalOrdinal localBlockID) const =0;

  //! Return the localBlockID corresponding to the given globalBlockID
  /*! If globalBlockID is not present on this processor, returns Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
  */
  virtual LocalOrdinal getLocalBlockID(GlobalOrdinal globalBlockID) const =0;

  //! Return the block-size for localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  virtual LocalOrdinal getLocalBlockSize(LocalOrdinal localBlockID) const =0;

  //! Return the first local point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  virtual LocalOrdinal getFirstLocalPointInLocalBlock(LocalOrdinal localBlockID) const =0;

  //! Return the first global point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  virtual GlobalOrdinal getFirstGlobalPointInLocalBlock(LocalOrdinal localBlockID) const =0;

  //@}

#endif // XPETRA_NOT_IMPLEMENTED_FOR_EPETRA
};//class BlockMap

//-----------------------------------------------------------------
// template<class LocalOrdinal,class GlobalOrdinal,class Node>
// Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
// convertBlockMapToPointMap(const Teuchos::RCP<const Xpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockMap)

}//namespace Xpetra


#define XPETRA_BLOCKMAP_SHORT
#endif

