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

public:

  std::string inputAdapterName()const {return std::string("XpetraCrsMatrix");}

  ~XpetraCrsMatrixInput() { }

  /*! Constructor with an Xpetra::CrsMatrix
   */
  XpetraCrsMatrixInput(const RCP<const xmatrixType > matrix):
    _matrix(), _rowMap(), _colMap(), _base()
  {
   _matrix = matrix;
   _rowMap = _matrix->getRowMap();
   _colMap = _matrix->getColMap();
   _base = _rowMap->getIndexBase();
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
    std::vector<LID> &localIds, std::vector<LNO> &rowSize,
    std::vector<GID> &colIds) const
  {
    size_t nrows = getLocalNumRows();
    size_t nnz = _matrix->getNodeNumEntries();
    size_t maxrow = _matrix->getNodeMaxNumRowEntries();
    size_t next = 0;

    rowIds.resize(nrows);
    rowSize.resize(nrows);
    colIds.resize(nnz);
    localIds.clear();   // consecutive IDs implied

    Teuchos::Array<LNO> indices(maxrow);
    Teuchos::Array<Scalar> nzs(maxrow);

    for (unsigned i=0; i < nrows; i++){
      LNO row = i + _base;
      LNO nnz = _matrix->getNumEntriesInLocalRow(row);
      size_t n;
      _matrix->getLocalRowCopy(row, indices.view(0,nnz), nzs.view(0,nnz), n);
      for (LNO j=0; j < nnz; j++){
        colIds[next++] = _colMap->getGlobalElement(indices[j]);
      }
      rowIds[i] = _rowMap->getGlobalElement(row);
      rowSize[i] = nnz;
    }
  } 

  /*! Access to xpetra matrix
   */

  RCP<const xmatrixType> getMatrix()
  {
    return _matrix;
  }



  /*! Return a read only view of the data.
      We don't have a view of global data, only of local data.
   */
  //LNO getRowListView(GID *&rowIds, LID *&localIds,
  //  LNO *&rowSize, GID *& colIds)


};
  
}  //namespace Zoltan2
  
#endif
