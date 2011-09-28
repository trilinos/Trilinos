// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraCrsMatrixInput.hpp

    \brief An input adapter for a Xpetra::CrsMatrix.

    \author Siva Rajamanickam
*/

// TODO:  Should LID and GID be used anywhere except the interface 
// TODO:  in this class?  KDDKDD

#ifndef _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_
#define _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_

#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_Standards.hpp>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////////////
/*! Zoltan2::XpetraCrsMatrixInput
    \brief This objects provides access for Zoltan2 to Xpetra::CrsMatrix data.

*/
//////////////////////////////////////////////////////////////////////////////

CONSISTENT_CLASS_TEMPLATE_LINE
class XpetraCrsMatrixInput : public MatrixInput<CONSISTENT_TEMPLATE_PARAMS>
{
private:
  RCP<Xpetra::CrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> > _xmatrix;
  RCP<const Xpetra::Map<LNO, GNO, Node> > _rowMap;
  RCP<const Xpetra::Map<LNO, GNO, Node> > _colMap;
  LNO _base;

public:

  ///////////////////////////////////////////////////////////////////////////
  /*! Default constructor  TODO - get rid of this?
   */
  XpetraCrsMatrixInput(): _xmatrix(), _rowMap(),_colMap(),  _base() { }

  ///////////////////////////////////////////////////////////////////////////
  /*! Constructor with an Xpetra::CrsMatrix
   */
  XpetraCrsMatrixInput(RCP<Xpetra::CrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> >
                            matrix) : _xmatrix(matrix), _rowMap(),_colMap(),  _base()
  {
    _rowMap = _xmatrix->getRowMap();
    _colMap = _xmatrix->getColMap();
    _base = _rowMap->getIndexBase();
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Constructor with an Xpetra::TpetraCrsMatrix
   */
  XpetraCrsMatrixInput(RCP<Xpetra::TpetraCrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> >
                            matrix) : _xmatrix(), _rowMap(), _colMap(), _base()
  {
    _xmatrix = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix
                             <CONSISTENT_TRILINOS_TEMPLATE_PARAMS> > (matrix);
    _rowMap = _xmatrix->getRowMap();
    _colMap = _xmatrix->getColMap();
    _base = _rowMap->getIndexBase();
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Constructor with an Xpetra::EpetraCrsMatrix
   */
  XpetraCrsMatrixInput(RCP<Xpetra::EpetraCrsMatrix> matrix) : 
    _xmatrix(), _rowMap(), _colMap(), _base()
  {
    _xmatrix = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix
                             <CONSISTENT_TRILINOS_TEMPLATE_PARAMS> > (matrix);
    _rowMap = _xmatrix->getRowMap();
    _colMap = _xmatrix->getColMap();
    _base = _rowMap->getIndexBase();
  }

  ///////////////////////////////////////////////////////////////////////////
  ~XpetraCrsMatrixInput() {}

  ///////////////////////////////////////////////////////////////////////////
  //! String identifying the input adapter; useful for output. 
  std::string inputAdapterName() const {
    return std::string("XpetraCrsMatrixInput"); 
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // The MatrixInput adapter interface
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  /*! Returns the number rows on this process.
   */
  size_t getLocalNumRows() const
  {
    return _xmatrix->getNodeNumRows();
  }

  /*! Returns the number rows in the entire matrix.
   */
  global_size_t getGlobalNumRows() const
  {
    return _xmatrix->getGlobalNumRows();
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Returns true if input adapter uses local Ids.
   */
  bool haveLocalIds() const { return true; }

  ///////////////////////////////////////////////////////////////////////////
  /*! Return true if local Ids are consecutive integral
   *   values and supply the base.  Providing this information
   *   can save memory, making local Id lists unneccesary.
   */
  bool haveConsecutiveLocalIds(size_t &base) const{
    base = static_cast<size_t>(_base); 
    return true;
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Returns the number columns used by rows on this process
   */
  size_t getLocalNumColumns() const
  {
    return _xmatrix->getNodeNumCols();
  }

  /*! Returns the number columns in the entire matrix.
   */
  global_size_t getGlobalNumColumns() const
  {
    return _xmatrix->getGlobalNumCols();
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Return the total number of non-zero entries on this process.
  */
  size_t getLocalNumNonZeros() const
  {
    return _xmatrix->getNodeNumEntries();
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Return the maximum number of non-zero entries in any local row.
  */
  size_t getLocalMaxNumNonZeros() const
  {
    return _xmatrix->getNodeMaxNumRowEntries();
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Return the local row information
      \param Ids will on return hold a list of the global Ids for
        each row on this process.
      \param localIds can, optionally, on return hold a list of locally
        relevant values that the process will use to refer to the objects
        listed in the first list.  If localIds are omitted and
        haveConsecutiveLocalIds is true, it is assumed that the
        global Ids are in local Id order.
      \param nnz will on return hold the number of non zeros in the
         cooresponding row.
  */
  void getRowListCopy(std::vector<GID> &Ids,
    std::vector<LID> &localIds, std::vector<size_t> &nnz)
  {
    size_t numRows = this->getLocalNumRows();
    Ids.resize(numRows,0);
    localIds.clear();   // don't need it
    nnz.resize(numRows,0);
    LNO end = _base + numRows;

    for (LNO i= _base, j=0; i < end; i++, j++){
      Ids[j] = _rowMap->getGlobalElement(i);
      nnz[j] = _xmatrix->getNumEntriesInLocalRow(i);
    }
  }

  /*!  This optional method is not implemented because we don't have list of nnz.
   */
  // getRowListView(GID *&Ids, LID *&localIds, size_t *nnz)

  ///////////////////////////////////////////////////////////////////////////
  /*! Return the column Ids of the non-zeros for the given row.
      \param Id  global Id for a row on this process
      \param localId  app's local Id, if any, associated with this row
      \param columnId on return will contain the list of global column Ids
   
       TODO will model ever want the non zero values?
   */

  void getRowNonZeroCopy(GID Id, LID localId, std::vector<GID> &columnId)
  {
    size_t nCols = _xmatrix->getNumEntriesInLocalRow(localId);
    columnId.resize(nCols);

    Array<LNO> columnLid(nCols);
    Array<Scalar> temp(nCols);
    _xmatrix->getLocalRowView(localId, columnLid, temp);

    for (int i=0; i < nCols; i++){
       columnId[i] = _colMap->getGlobalElement(columnLid[i]);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Not defined because underlying xpetra object has local IDs now.
      TODO explain better.
   */
  //size_t getRowNonZeroView(GID Id, LID localId, GID *&columnId) const

  ///////////////////////////////////////////////////////////////////////////
  /*! Return true of matrix is globally lower triangular.   TODO
   */
  bool isLowerTriagular() const { return false; }
//TODO Not supported by Xpetra KDDKDD  bool isLowerTriangular() const { return _xmatrix->isLowerTriangular(); }

  ///////////////////////////////////////////////////////////////////////////
  /*! Return true of matrix is globally upper triangular.   TODO
   */
  bool isUpperTriangular() const { return false; }
//TODO Not supported by Xpetra KDDKDD  bool isUpperTriangular() const { return _xmatrix->isUpperTriangular(); }

  ///////////////////////////////////////////////////////////////////////////
  /*! Return true of matrix globally has any diagonal entries.
   */
  bool hasDiagonalEntries() const { return (_xmatrix->getGlobalNumDiags() > 0);}
};
  
  
}  //namespace Zoltan2
  
#endif
