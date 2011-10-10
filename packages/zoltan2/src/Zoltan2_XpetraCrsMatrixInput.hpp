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

  RCP<const xmatrixType > _matrix;
  RCP<const Xpetra::Map<LID, GID, Node> > _rowMap;
  RCP<const Xpetra::Map<LID, GID, Node> > _colMap;
  LID _base;
  ArrayRCP<LNO> _offsets;
  ArrayRCP<GNO> _columnIds;

public:

  std::string inputAdapterName()const {return std::string("XpetraCrsMatrix");}

  ~XpetraCrsMatrixInput() { }

  /*! Constructor with an Xpetra::CrsMatrix
   */
  XpetraCrsMatrixInput(const RCP<const xmatrixType > matrix):
    _matrix(), _rowMap(), _colMap(), _base(), _offsets(), _columnIds()
  {
   _matrix = matrix;
   _rowMap = _matrix->getRowMap();
   _colMap = _matrix->getColMap();
   _base = _rowMap->getIndexBase();

   size_t nrows = _matrix->getNodeNumRows();
   size_t nnz = _matrix->getNodeNumEntries();

    _offsets.resize(nrows+1, LNO(0));
    _columnIds.resize(nnz);
    ArrayView<const LNO> indices;
    ArrayView<const Scalar> nzs;
    LNO next = 0;
    for (unsigned i=0; i < nrows; i++){
      LNO row = i + _base;
      LNO nnz = _matrix->getNumEntriesInLocalRow(row);
      _matrix->getLocalRowView(row, indices, nzs);
      for (LNO j=0; j < nnz; j++){
        // TODO - this will be slow
        //   Is it possible that global columns ids might be stored in order?
        _columnIds[next++] = _colMap->getGlobalElement(indices[j]);
      }
      _offsets[i+1] = _offsets[i] + nnz;
    }
  }

  ////////////////////////////////////////////////////
  // The MatrixInput interface.
  ////////////////////////////////////////////////////

  /*! Returns the number rows on this process.
   */
  size_t getLocalNumRows() const { 
    return _matrix->getNodeNumRows();
  }

  /*! Returns the number rows in the entire matrix.
   */
  global_size_t getGlobalNumRows() const { 
    return _matrix->getGlobalNumRows();
  }

  /*! Return whether input adapter wants to use local IDs.
   */

  bool haveLocalIds() const {return true;}

  /*! Return whether local ids are consecutive and if so the base.
   */

  bool haveConsecutiveLocalIds (size_t &base) const
  {
    base = static_cast<size_t>(_base);
    return true;
  }

  /*! Returns the number columns on this process.
   */
  size_t getLocalNumColumns() const { 
    return _matrix->getNodeNumCols();
  }

  /*! Returns the number columns on this entire matrix.
   *    what about directional columns, count twice?
   */
  global_size_t getGlobalNumColumns() const { 
    return _matrix->getGlobalNumCols();
  }

  /*! Get copy of matrix entries on local process
   */
  void getRowListCopy(std::vector<GID> &rowIds,
    std::vector<LID> &localIds, std::vector<LNO> &offsets,
    std::vector<GID> &colIds) const
  {
    size_t nrows = getLocalNumRows();
    size_t nnz = _matrix->getNodeNumEntries();
    size_t maxrow = _matrix->getNodeMaxNumRowEntries();
    size_t next = 0;

    rowIds.resize(nrows);
    offsets.resize(nrows+1);
    colIds.resize(nnz);
    localIds.clear();   // consecutive integer IDs implied

    Teuchos::Array<LNO> indices(maxrow);
    Teuchos::Array<Scalar> nzs(maxrow);

    offsets[0] = 0;

    for (unsigned i=0; i < nrows; i++){
      LNO row = i + _base;
      LNO nnz = _matrix->getNumEntriesInLocalRow(row);
      size_t n;
      _matrix->getLocalRowCopy(row, indices.view(0,nnz), nzs.view(0,nnz), n);
      for (LNO j=0; j < nnz; j++){
        colIds[next++] = _colMap->getGlobalElement(indices[j]);
      }
      rowIds[i] = _rowMap->getGlobalElement(row);
      offsets[i+1] = offsets[i] + nnz;
    }
  } 

  /*! Access to xpetra matrix
   */

  RCP<const xmatrixType> getMatrix() const
  {
    return _matrix;
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

    ArrayView<const GID> rowView = _rowMap->getNodeElementList();
    rowIds = rowView.getRawPtr();
   
    localIds = NULL;   // Implies consecutive integers

    offsets = _offsets.getRawPtr();
    colIds = _columnIds.getRawPtr();
    return nrows;
  }
};
  
}  //namespace Zoltan2
  
#endif
