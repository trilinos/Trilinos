// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraCrsMatrixInput.hpp
    \brief Defines the XpetraCrsMatrixInput adapter class.
*/

#ifndef _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_
#define _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_

#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_Util.hpp>

#include <Xpetra_CrsMatrix.hpp>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////////////
/*!  \brief Provides access for Zoltan2 to Xpetra::CrsMatrix data.

    \todo we assume FillComplete has been called.  We should support
                objects that are not FillCompleted.
    \todo add RowMatrix

    The template parameter is the user's input object:
     \li Tpetra::CrsGraph
     \li Xpetra::CrsGraph
     \li Epetra_CrsGraph

    The \c scalar_t type, representing use data such as matrix values, is
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.

*/

template <typename User>
  class XpetraCrsMatrixInput : public MatrixInput<User> {
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> xmatrix_t;
  typedef MatrixInput<User>       base_adapter_t;
  typedef User user_t;
#endif

  /*! \brief Destructor
   */
  ~XpetraCrsMatrixInput() { }

  /*! \brief Constructor   
   *    \param inmatrix The users Epetra, Tpetra, or Xpetra CrsMatrix object 
   *    \param coordDim Some algorithms can use row or column geometric
   *            information if it is available.  If coordinates will be
   *            supplied in setRowCoordinates() 
   *            then provide the dimension of the coordinates here.
   */
  XpetraCrsMatrixInput(const RCP<const User> &inmatrix, int coordDim=0);

  /*! \brief Specify geometric coordinates for matrix rows.
   *    \param dim  A value between zero and one less that the \c coordDim
   *                  argument to the constructor.
   *    \param coordVal  A pointer to the coordinates.
   *    \stride          A stride to be used in reading the values.  The
   *        dimension \c dim coordinate for row \k should be found at
   *        <tt>coordVal[k*stride]</tt>.
   *
   * The order of coordinates should correspond to the order of rows
   * returned by
   *   \code
   *       theMatrix->getRowMap()->getNodeElementList();
   *   \endcode
   */
  void setRowCoordinates(int dim, const scalar_t *coordVal, int stride);

  /*! \brief Access to Xpetra-wrapped user matrix. 
   */

  const RCP<const xmatrix_t> &getMatrix() const
  {
    return matrix_;
  }

  ////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////

  string inputAdapterName()const { return string("XpetraCrsMatrix");}

  size_t getLocalNumberOfObjects() const { return getLocalNumRows();}

  int getNumberOfWeightsPerObject() const { return 0;}

  ////////////////////////////////////////////////////
  // The MatrixInput interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumRows() const { 
    return matrix_->getNodeNumRows();
  }

  global_size_t getGlobalNumRows() const { 
    return matrix_->getGlobalNumRows();
  }

  size_t getLocalNumColumns() const { 
    return matrix_->getNodeNumCols();
  }

  global_size_t getGlobalNumColumns() const { 
    return matrix_->getGlobalNumCols();
  }

  bool diagonalEntriesMayBePresent() const {
    return ((matrix_->getGlobalNumCols() == matrix_->getGlobalNumRows()) && 
            (matrix_->getGlobalNumDiags() > 0));
  }

  size_t getRowListView(const gid_t *&rowIds,
    const lno_t *&offsets, const gid_t *& colIds) const
  {
    size_t nrows = getLocalNumRows();

    ArrayView<const gid_t> rowView = rowMap_->getNodeElementList();
    rowIds = rowView.getRawPtr();
    offsets = offset_.getRawPtr();
    colIds = columnIds_.getRawPtr();
    return nrows;
  }

  int getCoordinateDimension() const {return coordinateDim_;}

  size_t getRowCoordinates(int dim,
    const scalar_t *&coords, int &stride) const
  {
    env_->localInputAssertion(__FILE__, __LINE__,
      "invalid coordinate dimension",
      dim >= 0 && dim < coordinateDim_, BASIC_ASSERTION);

    size_t length;
    rowCoords_[dim].getStridedList(length, coords, stride);
    return length;
  }

  ////////////////////////////////////////////////////
  // End of MatrixInput interface.
  ////////////////////////////////////////////////////

  template <typename Adapter>
    size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<Adapter> &solution) const;

private:

  RCP<Environment> env_;    // for error messages, etc.

  RCP<const User> inmatrix_;
  RCP<const xmatrix_t> matrix_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > rowMap_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > colMap_;
  lno_t base_;
  ArrayRCP<lno_t> offset_;
  ArrayRCP<gno_t> columnIds_;

  int coordinateDim_;
  ArrayRCP<StridedData<lno_t, scalar_t> > rowCoords_;

};

/////////////////////////////////////////////////////////////////
// Definitions
/////////////////////////////////////////////////////////////////

template <typename User>
  XpetraCrsMatrixInput<User>::XpetraCrsMatrixInput(
    const RCP<const User> &inmatrix, int coordDim):
      env_(rcp(new Environment)),
      inmatrix_(inmatrix), matrix_(), rowMap_(), colMap_(), base_(),
      offset_(), columnIds_(),
      coordinateDim_(coordDim), rowCoords_()
{
  typedef StridedData<lno_t,scalar_t> input_t;
  matrix_ = XpetraTraits<User>::convertToXpetra(inmatrix);
  rowMap_ = matrix_->getRowMap();
  colMap_ = matrix_->getColMap();
  base_ = rowMap_->getIndexBase();

  size_t nrows = matrix_->getNodeNumRows();
  size_t nnz = matrix_->getNodeNumEntries();
 
  offset_.resize(nrows+1, 0);
  columnIds_.resize(nnz);
  ArrayView<const lno_t> indices;
  ArrayView<const scalar_t> nzs;
  lno_t next = 0;
  for (size_t i=0; i < nrows; i++){
    lno_t row = i + base_;
    lno_t nnz = matrix_->getNumEntriesInLocalRow(row);
    matrix_->getLocalRowView(row, indices, nzs);
    for (lno_t j=0; j < nnz; j++){
      // TODO - this will be slow
      //   Is it possible that global columns ids might be stored in order?
      columnIds_[next++] = colMap_->getGlobalElement(indices[j]);
    }
    offset_[i+1] = offset_[i] + nnz;
  } 

  if (coordinateDim_ > 0)
    rowCoords_ = arcp(new input_t [coordinateDim_], 0, coordinateDim_, true);
}

// TODO (from 3/21/12 mtg):  Consider changing interface to take an XpetraMultivector
template <typename User>
  void XpetraCrsMatrixInput<User>::setRowCoordinates(int dim,
    const scalar_t *coordVal, int stride)
{
  typedef StridedData<lno_t,scalar_t> input_t;

  env_->localInputAssertion(__FILE__, __LINE__, 
    "invalid row coordinate dimension",
    dim >= 0 && dim < coordinateDim_, BASIC_ASSERTION);

  size_t nvtx = getLocalNumRows();

  ArrayRCP<const scalar_t> coordV(coordVal, 0, nvtx, false);
  rowCoords_[dim] = input_t(coordV, stride);
}

template <typename User>
  template <typename Adapter>
    size_t XpetraCrsMatrixInput<User>::applyPartitioningSolution(
      const User &in, User *&out, 
      const PartitioningSolution<Adapter> &solution) const
{ 
  // Get an import list

  size_t len = solution.getLocalNumberOfIds();
  const gid_t *gids = solution.getIdList();
  const partId_t *parts = solution.getPartList();
  ArrayRCP<gid_t> gidList = arcp(const_cast<gid_t *>(gids), 0, len, false); 
  ArrayRCP<partId_t> partList = arcp(const_cast<partId_t *>(parts), 0, len, 
    false); 
  ArrayRCP<lno_t> dummyIn;
  ArrayRCP<gid_t> importList;
  ArrayRCP<lno_t> dummyOut;
  size_t numNewRows;
  const RCP<const Comm<int> > comm = matrix_->getRowMap()->getComm();

  try{
    numNewRows = convertSolutionToImportList<Adapter, lno_t>(
      solution, dummyIn, importList, dummyOut);
  }
  Z2_FORWARD_EXCEPTIONS;

for (int i=0; i < numNewRows; i++){
std::cout << importList[i] << " ";
}
std::cout << std::endl;

  gno_t lsum = numNewRows;
  gno_t gsum = 0;
  reduceAll<int, gno_t>(*comm, Teuchos::REDUCE_SUM, 1, &lsum, &gsum);

  RCP<const User> inPtr = rcp(&in, false);

  RCP<const User> outPtr = XpetraTraits<User>::doMigration(
   inPtr, lsum, importList.getRawPtr());

  out = const_cast<User *>(outPtr.get());
  outPtr.release();
  return numNewRows;
}

}  //namespace Zoltan2
  
#endif
