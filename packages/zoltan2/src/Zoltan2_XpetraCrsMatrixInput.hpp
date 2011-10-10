// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraCrsMatrixInput.hpp

    \brief An input adapter for a Xpetra::CrsMatrix.
*/

#ifndef _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_
#define _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_

#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Zoltan2_MatrixInput.hpp>

namespace Zoltan2 {

/*! Zoltan2::XpetraCrsMatrixInput
    \brief Provides access for Zoltan2 to Xpetra::CrsMatrix data.

    The template parameter is the weight type.  Xpetra local and global IDs 
    are ints.
    TODO: we assume FillComplete has been called.  We should support
                objects that are not FillCompleted.

*/

template <Z2CLASS_TEMPLATE>
class XpetraCrsMatrixInput : public MatrixInput<Z2PARAM_TEMPLATE> {
private:

  typedef Xpetra::CrsMatrix<Scalar, LNO, GNO, Node> xmatrixType;

  RCP<const xmatrixType > matrix_;
  RCP<const Xpetra::Map<LID, GID, Node> > rowMap_;
  RCP<const Xpetra::Map<LID, GID, Node> > colMap_;
  LID base_;
  ArrayRCP<LNO> offset_;
  ArrayRCP<GNO> columnIds_;

public:

  /*! Name of input adapter type
   */
  std::string inputAdapterName()const {return std::string("XpetraCrsMatrix");}

  /*! Destructor
   */
  ~XpetraCrsMatrixInput() { }

  /*! Constructor 
   */
  XpetraCrsMatrixInput(const RCP<const xmatrixType > matrix):
    matrix_(), rowMap_(), colMap_(), base_(), offset_(), columnIds_()
  {
   matrix_ = matrix;
   rowMap_ = matrix_->getRowMap();
   colMap_ = matrix_->getColMap();
   base_ = rowMap_->getIndexBase();

   size_t nrows = matrix_->getNodeNumRows();
   size_t nnz = matrix_->getNodeNumEntries();

    offset_.resize(nrows+1, LNO(0));
    columnIds_.resize(nnz);
    ArrayView<const LNO> indices;
    ArrayView<const Scalar> nzs;
    LNO next = 0;
    for (unsigned i=0; i < nrows; i++){
      LNO row = i + base_;
      LNO nnz = matrix_->getNumEntriesInLocalRow(row);
      matrix_->getLocalRowView(row, indices, nzs);
      for (LNO j=0; j < nnz; j++){
        // TODO - this will be slow
        //   Is it possible that global columns ids might be stored in order?
        columnIds_[next++] = colMap_->getGlobalElement(indices[j]);
      }
      offset_[i+1] = offset_[i] + nnz;
    }
  }

  ////////////////////////////////////////////////////
  // The MatrixInput interface.
  ////////////////////////////////////////////////////

  /*! Returns the number rows on this process.
   */
  size_t getLocalNumRows() const { 
    return matrix_->getNodeNumRows();
  }

  /*! Returns the number rows in the entire matrix.
   */
  global_size_t getGlobalNumRows() const { 
    return matrix_->getGlobalNumRows();
  }

  /*! Return whether input adapter wants to use local IDs.
   */

  bool haveLocalIds() const {return true;}

  /*! Return whether local ids are consecutive and if so the base.
   */

  bool haveConsecutiveLocalIds (size_t &base) const
  {
    base = static_cast<size_t>(base_);
    return true;
  }

  /*! Returns the number columns on this process.
   */
  size_t getLocalNumColumns() const { 
    return matrix_->getNodeNumCols();
  }

  /*! Returns the number columns on this entire matrix.
   *    what about directional columns, count twice?
   */
  global_size_t getGlobalNumColumns() const { 
    return matrix_->getGlobalNumCols();
  }

  /*! Get copy of matrix entries on local process
   */
  void getRowListCopy(std::vector<GID> &rowIds,
    std::vector<LID> &localIds, std::vector<LNO> &offsets,
    std::vector<GID> &colIds) const
  {
    size_t nrows = getLocalNumRows();
    size_t nnz = matrix_->getNodeNumEntries();
    size_t maxrow = matrix_->getNodeMaxNumRowEntries();
    size_t next = 0;

    rowIds.resize(nrows);
    offsets.resize(nrows+1);
    colIds.resize(nnz);
    localIds.clear();   // consecutive integer IDs implied

    Teuchos::Array<LNO> indices(maxrow);
    Teuchos::Array<Scalar> nzs(maxrow);

    offsets[0] = 0;

    for (unsigned i=0; i < nrows; i++){
      LNO row = i + base_;
      LNO nnz = matrix_->getNumEntriesInLocalRow(row);
      size_t n;
      matrix_->getLocalRowCopy(row, indices.view(0,nnz), nzs.view(0,nnz), n);
      for (LNO j=0; j < nnz; j++){
        colIds[next++] = colMap_->getGlobalElement(indices[j]);
      }
      rowIds[i] = rowMap_->getGlobalElement(row);
      offsets[i+1] = offsets[i] + nnz;
    }
  } 

  /*! Access to xpetra matrix
   */

  RCP<const xmatrixType> getMatrix() const
  {
    return matrix_;
  }



  /*! Return a read only view of the data.
     \param rowIds  Global row ids.  The memory for the global 
          row IDs persists until the underlying Xpetra::CrsMatrix is deleted.
     \param localIds on return is NULL, signifying that local IDs are
           contiguous integers starting at 0.
     \param offsets The columns for rowIds[i] begin at colIds[offsets[i]].  There are
           numRows+1 offsets.  The last offset is the length of the colIds array.
           The memory pointed to by offsets persists in the GraphModel is deleted.
     \param colIds The global column Ids. The memory pointed to by colIds 
          persists until the GraphModel is deleted.

     \return  The number rows in the rowIds list is returned.
   */

  size_t getRowListView(const GID *&rowIds, const LID *&localIds,
    const LNO *&offsets, const GID *& colIds) const
  {
    size_t nrows = getLocalNumRows();

    ArrayView<const GID> rowView = rowMap_->getNodeElementList();
    rowIds = rowView.getRawPtr();
   
    localIds = NULL;   // Implies consecutive integers

    offsets = offset_.getRawPtr();
    colIds = columnIds_.getRawPtr();
    return nrows;
  }
};
  
}  //namespace Zoltan2
  
#endif
