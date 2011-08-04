#ifndef CTHULHU_EPETRABLOCKMAP_HPP
#define CTHULHU_EPETRABLOCKMAP_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_EPETRA
#error This file should be included only if HAVE_CTHULHU_EPETRA is defined.
#endif

#include "Cthulhu_Map.hpp"
#include "Cthulhu_BlockMap.hpp"

#include "Epetra_BlockMap.h"
#include "Cthulhu_Comm.hpp"
#include "Cthulhu_EpetraExceptions.hpp"

#include "Cthulhu_Debug.hpp"

namespace Cthulhu {

/*!
  @class EpetraBlockMap
  @brief Block-entry counterpart to Cthulhu::Map.

  BlockMap doesn't inherit Cthulhu::Map
*/
class EpetraBlockMap : public Cthulhu::BlockMap<int,int> {
 public:

  //! @name Constructor/Destructor Methods
  //@{

  /*! \brief EpetraBlockMap constructor specifying numGlobalBlocks and constant blockSize.
   */
  EpetraBlockMap(global_size_t numGlobalBlocks,
                 int blockSize,
                 int indexBase,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                 const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) : map_(rcp(new Epetra_BlockMap(numGlobalBlocks, blockSize, indexBase, *Teuchos2Epetra_Comm(comm)))) { CTHULHU_DEBUG_ME; }

  //TODO 
  //CATCH_EPETRA_EXCEPTION_AND_THROW_INVALID_ARG((map_ = (rcp(new Epetra_Map(numGlobalElements, indexBase, *Teuchos2Epetra_Comm(comm))))););

  /*! \brief EpetraBlockMap constructor specifying num global and local blocks, and constant blockSize.
   */
  EpetraBlockMap(global_size_t numGlobalBlocks,
                 size_t numLocalBlocks,
                 int blockSize,
                 int indexBase,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                 const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) : map_(rcp(new Epetra_BlockMap(numGlobalBlocks, numLocalBlocks, blockSize, indexBase, *Teuchos2Epetra_Comm(comm)))) { CTHULHU_DEBUG_ME; }

  /*! \brief EpetraBlockMap constructor specifying numGlobalBlocks and lists of local blocks first-global-point-in-blocks, and blockSizes.
   */
 // EpetraBlockMap(global_size_t numGlobalBlocks,
//                  const Teuchos::ArrayView<const int>& myGlobalBlockIDs,
//                  const Teuchos::ArrayView<const int>& myFirstGlobalPointInBlocks,
//                  const Teuchos::ArrayView<const int>& myBlockSizes,
//                  int indexBase,
//                  const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
//                  const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) {
//     CTHULHU_DEBUG_ME; 
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
//     CTHULHU_DEBUG_ME;
//     map_ = rcp(new Epetra_BlockMap(pointMap, myGlobalBlockIDs, myBlockSizes));
//   }

  EpetraBlockMap(const Teuchos::RCP<const Epetra_BlockMap > &map) : map_(map) { CTHULHU_DEBUG_ME; }

  //! EpetraBlockMap destructor.
  virtual ~EpetraBlockMap(){ CTHULHU_DEBUG_ME; }

  //@}

  //! @name Attribute Accessor Methods
  //@{
#ifdef CTHULHU_NOT_IMPLEMENTED
  inline const Teuchos::RCP<const Map<int,int> >& getPointMap() const { CTHULHU_DEBUG_ME; return map_->getPointMap(); }
#endif

  inline global_size_t getGlobalNumBlocks() const { CTHULHU_DEBUG_ME; return map_->NumGlobalElements(); }

  //! Return number of blocks on the local processor.
  inline size_t getNodeNumBlocks() const { CTHULHU_DEBUG_ME; return map_->NumMyElements(); }

  //TODO inline Teuchos::ArrayView<const int> getNodeBlockIDs() const { CTHULHU_DEBUG_ME; return map_->getNodeBlockIDs(); }

  inline bool isBlockSizeConstant() const { CTHULHU_DEBUG_ME; return map_->ConstantElementSize(); }

  //! Return ArrayRCP of first-local-point in local blocks.
  //TODO inline Teuchos::ArrayRCP<const int> getNodeFirstPointInBlocks() const { CTHULHU_DEBUG_ME; return map_->FirstPointInElementList(); }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
  //! Return device-resident ArrayRCP of first-local-point in local blocks.
  /*! This version of this method is primarily used internally by VbrMatrix
      for passing data to the matrix-vector-product kernel.
  */
  inline Teuchos::ArrayRCP<const int> getNodeFirstPointInBlocks_Device() const { CTHULHU_DEBUG_ME; return map_->getNodeFirstPointInBlocks_Device(); }
#endif

  //! Return the globalBlockID corresponding to the given localBlockID
  /*! If localBlockID is not present on this processor, returns Teuchos::OrdinalTraits<int>::invalid().
  */
  inline int getGlobalBlockID(int localBlockID) const { CTHULHU_DEBUG_ME; return map_->GID(localBlockID); }

  //! Return the localBlockID corresponding to the given globalBlockID
  /*! If globalBlockID is not present on this processor, returns Teuchos::OrdinalTraits<int>::invalid().
  */
  inline int getLocalBlockID(int globalBlockID) const { CTHULHU_DEBUG_ME; return map_->LID(globalBlockID); }

  //! Return the block-size for localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  inline int getLocalBlockSize(int localBlockID) const { CTHULHU_DEBUG_ME; return map_->ElementSize(localBlockID); }

  //! Return the first local point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  inline int getFirstLocalPointInLocalBlock(int localBlockID) const { CTHULHU_DEBUG_ME; return map_->FirstPointInElement(localBlockID); }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
  //! Return the first global point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown. //TODO: throw !
   */
  inline int getFirstGlobalPointInLocalBlock(int localBlockID) const { CTHULHU_DEBUG_ME; return; }
#endif

  //@}

  RCP< const Epetra_BlockMap > getEpetra_BlockMap() const { CTHULHU_DEBUG_ME; return map_; }

private:
  const RCP< const Epetra_BlockMap > map_;


};//class EpetraBlockMap

// //-----------------------------------------------------------------
// template<class int,class int,class Node>
// Teuchos::RCP<const Cthulhu::Map<int,int> >
// convertEpetraBlockMapToEpetraPointMap(const Teuchos::RCP<const Cthulhu::EpetraBlockMap<int,int> >& blockMap) { CTHULHU_DEBUG_ME;
//   return rcp(new EpetraMap(convertEpetraBlockMapToEpetraPointMap(blockMap.getEpetra_BlockMap())));
// }

}//namespace Cthulhu

#endif
