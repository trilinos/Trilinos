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

// TODO:  Should LID and GID be used at all in this class?  KDDKDD

#ifndef _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_
#define _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_

#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_TemplateMacros.hpp>

#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>

#include <Teuchos_RCP.hpp>

using Teuchos::RCP;

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
  bool _valid;
  Teuchos::RCP<Xpetra::CrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> > _xmatrix;
  Teuchos::RCP<const Xpetra::Map<LNO, GNO, Node> > _rowMap;
  LNO _base;

public:

  ///////////////////////////////////////////////////////////////////////////
  /*! Default constructor
   */
  XpetraCrsMatrixInput(): _valid(false), _xmatrix(), _rowMap(), _base() { }

  ///////////////////////////////////////////////////////////////////////////
  /*! Constructor with an Xpetra::CrsMatrix
   */
  XpetraCrsMatrixInput(RCP<Xpetra::CrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> >
                            matrix) : _valid(false), _xmatrix(matrix), _rowMap(), _base()
  {
    _rowMap = _xmatrix->getRowMap();
    _base = _rowMap->getIndexBase();
    _valid = true;
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Constructor with an Xpetra::TpetraCrsMatrix
   */
  XpetraCrsMatrixInput(RCP<Xpetra::TpetraCrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> >
                            matrix) : _valid(false), _xmatrix(), _rowMap(), _base()
  {
    _xmatrix = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix
                             <CONSISTENT_TRILINOS_TEMPLATE_PARAMS> > (matrix);
    _rowMap = _xmatrix->getRowMap();
    _base = _rowMap->getIndexBase();
    _valid = true;
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Constructor with an Xpetra::EpetraCrsMatrix
   */
  XpetraCrsMatrixInput(RCP<Xpetra::EpetraCrsMatrix> matrix) : 
    _valid(false), _xmatrix(), _rowMap(), _base()
  {
    _xmatrix = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix
                             <CONSISTENT_TRILINOS_TEMPLATE_PARAMS> > (matrix);
    _rowMap = _xmatrix->getRowMap();
    _base = _rowMap->getIndexBase();
    _valid = true;
  }

  ///////////////////////////////////////////////////////////////////////////
  ~XpetraCrsMatrixInput() {}

  ///////////////////////////////////////////////////////////////////////////
  //! String identifying the input adapter; useful for output. 
  std::string inputAdapterName() const {
    return std::string("XpetraCrsMatrixInput"); 
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Post construction setting of the underlying matrix.
   */
  void setMatrix(RCP<Xpetra::CrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> >
                            matrix)
  {
    if (_valid)
      throw std::runtime_error("Can not change the matrix in the adapter");
    _xmatrix = matrix;
    _rowMap = _xmatrix->getRowMap();
    _base = _rowMap->getIndexBase();
    _valid = true;
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Post construction setting of the underlying matrix.
   */
  void setMatrix(RCP<Xpetra::TpetraCrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> >
                            matrix)
  {
    if (_valid)
      throw std::runtime_error("Can not change the matrix in the adapter");
    _xmatrix = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix
                             <CONSISTENT_TRILINOS_TEMPLATE_PARAMS> > (matrix);
    _rowMap = _xmatrix->getRowMap();
    _base = _rowMap->getIndexBase();
    _valid = true;
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Post construction setting of the underlying matrix.
   */
  void setMatrix(RCP<Xpetra::EpetraCrsMatrix> matrix)
  {
    if (_valid)
      throw std::runtime_error("Can not change the matrix in the adapter");
    _xmatrix = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix
                             <CONSISTENT_TRILINOS_TEMPLATE_PARAMS> > (matrix);
    _rowMap = _xmatrix->getRowMap();
    _base = _rowMap->getIndexBase();
    _valid = true;
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // The MatrixInput adapter interface
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  /*! Returns the number rows on this process.
   */
  LNO getLocalNumRows() const
  {
    return _xmatrix->getNodeNumRows();
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
  bool haveConsecutiveLocalIds(LID &base) const{
    base = _base; 
    return true;
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Returns the number columns used by rows on this process
   */
  LNO getLocalNumColumns() const
  {
    return _xmatrix->getNodeNumCols();
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Return the total number of non-zero entries on this process.
  */
  LNO getLocalNumNonZeros() const
  {
    return _xmatrix->getNodeNumEntries();
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Return the maximum number of non-zero entries in any local row.
  */
  LNO getLocalMaxRowNumNonZeros() const
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
    std::vector<LID> &localIds, std::vector<LNO> &nnz)
  {
    LNO numRows = this->getLocalNumRows();
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
  // LNO getRowListView(GID *&Ids, LID *&localIds, LNO *nnz)

  ///////////////////////////////////////////////////////////////////////////
  /*! Return the column Ids of the non-zeros for the given row.
      \param Id  global Id for a row on this process
      \param localId  app's local Id, if any, associated with this row
      \param columnId on return will contain the list of global column Ids
   */

  void getRowNonZeroCopy(GID Id, LID localId, std::vector<GID> &columnId)
  {
    size_t nCols = _xmatrix->getNumEntriesInLocalRow(localId);  // TODO:  Is localID the correct argument here??  Previously, argument was "i", which didn't compile. KDDKDD   
    columnId.resize(nCols);

    // TODO:  This line will not compile for me; compiler complains about the template argument to Teuchos::Array.  KDDKDD  Teuchos::Array<std::vector<GID>::iterator> cols(columnId.begin(), columnId.end());
    Teuchos::Array<Scalar> values(nCols);

    // TODO:  I don't see this method in the Xpetra CrsMatrix class.  KDDKDD _xmatrix->getGlobalRowCopy(Id, cols, values);
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Obtain a read-only view, if possible, of the column Ids of the
      input row.
      \param Id  global Id for a row on this process
      \param localId  if input adapter supplied local Ids, this
         is that localId
      \param columnId on return will point a list of global column global Ids.
      \return The number of ids in the columnId list.
   */
  LNO getRowNonZeroView(GID Id, LID localId, GID *&columnId) const
  {
    Teuchos::ArrayView<GID> cols;
    Teuchos::ArrayView<GID> values;

    _xmatrix->getGlobalRowView(Id, cols, values);
    columnId = cols->getRawPtr();
    return 0;
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Return maximum number of nonzeros across all local rows.
   */
  LNO getLocalMaxNumNonZeros() const { 
   return -1; //TODO Placeholder to allow compilation
  }

  ///////////////////////////////////////////////////////////////////////////
  /*! Return true of matrix is globally lower triangular.
   */
//TODO Not supported by Xpetra KDDKDD  bool isLowerTriangular() const { return _xmatrix->isLowerTriangular(); }

  ///////////////////////////////////////////////////////////////////////////
  /*! Return true of matrix is globally upper triangular.
   */
//TODO Not supported by Xpetra KDDKDD  bool isUpperTriangular() const { return _xmatrix->isUpperTriangular(); }

  ///////////////////////////////////////////////////////////////////////////
  /*! Return true of matrix globally has any diagonal entries.
   */
  bool hasDiagonalEntries() const { return (_xmatrix->getGlobalNumDiags() > 0);}
};
  
  
}  //namespace Zoltan2
  
#endif
