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

