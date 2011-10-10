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

template <User>
class XpetraCrsMatrixInput : public MatrixInput<User> {
public:

  // Adapter must define these types for Zoltan2.
  typedef typename User::Scalar        scalar_t;
  typedef typename User::GlobalOrdinal gno_t;
  typedef typename User::LocalOrdinal  lno_t;
  typedef typename User::Node          node_t;
  typedef typename gno_t gid_t;   // GNO and GID are same in Xpetra.
  typedef typename lno_t lid_t;   // LNO and LID are same in Xpetra.

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
  void getRowListCopy(std::vector<gid_t> &rowIds,
    std::vector<lid_t> &localIds, std::vector<lno_t> &rowSize,
    std::vector<gid_t> &colIds) const
  {
    size_t nrows = getLocalNumRows();
    size_t nnz = _matrix->getNodeNumEntries();
    size_t maxrow = _matrix->getNodeMaxNumRowEntries();
    size_t next = 0;

    rowIds.resize(nrows);
    rowSize.resize(nrows);
    colIds.resize(nnz);
    localIds.clear();   // consecutive IDs implied

    Teuchos::Array<lno_t> indices(maxrow);
    Teuchos::Array<Scalar> nzs(maxrow);

    for (unsigned i=0; i < nrows; i++){
      lno_t row = i + _base;
      lno_t nnz = _matrix->getNumEntriesInLocalRow(row);
      size_t n;
      _matrix->getLocalRowCopy(row, indices.view(0,nnz), nzs.view(0,nnz), n);
      for (lno_t j=0; j < nnz; j++){
        colIds[next++] = _colMap->getGlobalElement(indices[j]);
      }
      rowIds[i] = _rowMap->getGlobalElement(row);
      rowSize[i] = nnz;
    }
  } 

  /*! Access to xpetra matrix
   */

  RCP<const xmatrixType> getMatrix() const
  {
    return _matrix;
  }



  /*! Return a read only view of the data.
      We don't have a view of global data, only of local data.
   */
  //lno_t getRowListView(gid_t *&rowIds, lid_t *&localIds,
  //  lno_t *&rowSize, gid_t *& colIds) const

private:

  typedef Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> xmatrixType;

  RCP<const xmatrixType > _matrix;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > _rowMap;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > _colMap;
  lno_t _base;


};
  
}  //namespace Zoltan2
  
#endif
