// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_XpetraCrsMatrixAdapter.hpp
    \brief Defines the XpetraCrsMatrixAdapter class.
*/

#ifndef _ZOLTAN2_XPETRACRSMATRIXADAPTER_HPP_
#define _ZOLTAN2_XPETRACRSMATRIXADAPTER_HPP_

#include <Zoltan2_MatrixAdapter.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_PartitioningHelpers.hpp>

#include <Xpetra_CrsMatrix.hpp>

#include <iostream>
#include <cassert>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////////////
/*!  \brief Provides access for Zoltan2 to Xpetra::CrsMatrix data.

    \todo we assume FillComplete has been called.  We should support
                objects that are not FillCompleted.
    \todo add RowMatrix

    The template parameter is the user's input object:
     \li Tpetra::CrsMatrix
     \li Xpetra::CrsMatrix
     \li Epetra_CrsMatrix

    The \c scalar_t type, representing use data such as matrix values, is
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.

*/

template <typename User, typename UserCoord=User>
  class XpetraCrsMatrixAdapter : public MatrixAdapter<User, UserCoord> {
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  using scalar_t = typename InputTraits<User>::scalar_t;
  using lno_t = typename InputTraits<User>::lno_t;
  using gno_t = typename InputTraits<User>::gno_t;
  using part_t = typename InputTraits<User>::part_t;
  using node_t = typename InputTraits<User>::node_t;
  using offset_t = typename InputTraits<User>::offset_t;
  using xmatrix_t = Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t>;

  using userCoord_t = UserCoord;
  using user_t = User;
#endif

  /*! \brief Constructor
   *    \param inmatrix The user's Epetra, Tpetra, or Xpetra CrsMatrix object
   *    \param nWeightsPerRow If row weights will be provided in setRowWeights(),
   *        then set \c nWeightsPerRow to the number of weights per row.
   */
  XpetraCrsMatrixAdapter(const RCP<const User> &inmatrix,
                         int nWeightsPerRow=0);

  /*! \brief Specify a weight for each entity of the primaryEntityType.
   *    \param weightVal A pointer to the weights for this index.
   *    \stride          A stride to be used in reading the values.  The
   *        index \c idx weight for entity \k should be found at
   *        <tt>weightVal[k*stride]</tt>.
   *    \param idx  A value between zero and one less that the \c nWeightsPerRow
   *                  argument to the constructor.
   *
   * The order of weights should correspond to the order of the primary
   * entity type; see, e.g.,  setRowWeights below.
   */

  void setWeights(const scalar_t *weightVal, int stride, int idx = 0);

  /*! \brief Specify a weight for each row.
   *    \param weightVal A pointer to the weights for this index.
   *    \stride          A stride to be used in reading the values.  The
   *        index \c idx weight for row \k should be found at
   *        <tt>weightVal[k*stride]</tt>.
   *    \param idx  A value between zero and one less that the \c nWeightsPerRow
   *                  argument to the constructor.
   *
   * The order of weights should correspond to the order of rows
   * returned by
   *   \code
   *       theMatrix->getRowMap()->getLocalElementList();
   *   \endcode
   */

  void setRowWeights(const scalar_t *weightVal, int stride, int idx = 0);

  /*! \brief Specify an index for which the weight should be
              the degree of the entity
   *    \param idx Zoltan2 will use the entity's
   *         degree as the entity weight for index \c idx.
   */
  void setWeightIsDegree(int idx);

  /*! \brief Specify an index for which the row weight should be
              the global number of nonzeros in the row
   *    \param idx Zoltan2 will use the global number of nonzeros in a row
   *         as the row weight for index \c idx.
   */
  void setRowWeightIsNumberOfNonZeros(int idx);

  ////////////////////////////////////////////////////
  // The MatrixAdapter interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumRows() const {
    return matrix_->getLocalNumRows();
  }

  size_t getLocalNumColumns() const {
    return matrix_->getLocalNumCols();
  }

  size_t getLocalNumEntries() const {
    return matrix_->getLocalNumEntries();
  }

  bool CRSViewAvailable() const { return true; }

  void getRowIDsView(const gno_t *&rowIds) const
  {
    ArrayView<const gno_t> rowView = rowMap_->getLocalElementList();
    rowIds = rowView.getRawPtr();
  }

  void getCRSView(ArrayRCP<const offset_t> &offsets,
                  ArrayRCP<const gno_t> &colIds) const
  {
    ArrayRCP< const lno_t > localColumnIds;
    ArrayRCP<const scalar_t> values;
    matrix_->getAllValues(offsets,localColumnIds,values);
    colIds = columnIds_;
  }

  void getCRSView(ArrayRCP<const offset_t> &offsets,
                  ArrayRCP<const gno_t> &colIds,
                  ArrayRCP<const scalar_t> &values) const {
    ArrayRCP< const lno_t > localColumnIds;
    matrix_->getAllValues(offsets,localColumnIds,values);
    colIds = columnIds_;
  }


  int getNumWeightsPerRow() const { return nWeightsPerRow_; }

  void getRowWeightsView(const scalar_t *&weights, int &stride,
                           int idx = 0) const
  {
    if(idx<0 || idx >= nWeightsPerRow_)
    {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid row weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
    }

    size_t length;
    rowWeights_[idx].getStridedList(length, weights, stride);
  }

  bool useNumNonzerosAsRowWeight(int idx) const { return numNzWeight_[idx];}

  template <typename Adapter>
    void applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<Adapter> &solution) const;

  template <typename Adapter>
    void applyPartitioningSolution(const User &in, RCP<User> &out,
         const PartitioningSolution<Adapter> &solution) const;

private:

  RCP<const User> inmatrix_;
  RCP<const xmatrix_t> matrix_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > rowMap_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > colMap_;
  lno_t base_;
  ArrayRCP<gno_t> columnIds_;  // TODO:  Refactor adapter to localColumnIds_

  int nWeightsPerRow_;
  ArrayRCP<StridedData<lno_t, scalar_t> > rowWeights_;
  ArrayRCP<bool> numNzWeight_;

  bool mayHaveDiagonalEntries;
};

/////////////////////////////////////////////////////////////////
// Definitions
/////////////////////////////////////////////////////////////////

template <typename User, typename UserCoord>
  XpetraCrsMatrixAdapter<User,UserCoord>::XpetraCrsMatrixAdapter(
    const RCP<const User> &inmatrix, int nWeightsPerRow):
      inmatrix_(inmatrix), matrix_(), rowMap_(), colMap_(),
      columnIds_(),
      nWeightsPerRow_(nWeightsPerRow), rowWeights_(), numNzWeight_(),
      mayHaveDiagonalEntries(true)
{
  typedef StridedData<lno_t,scalar_t> input_t;
  try {
    matrix_ = rcp_const_cast<const xmatrix_t>(
           XpetraTraits<User>::convertToXpetra(rcp_const_cast<User>(inmatrix)));
  }
  Z2_FORWARD_EXCEPTIONS

  rowMap_ = matrix_->getRowMap();
  colMap_ = matrix_->getColMap();

  size_t nrows = matrix_->getLocalNumRows();
  size_t nnz = matrix_->getLocalNumEntries();

  // Get ArrayRCP pointers to the structures in the underlying matrix
  ArrayRCP< const offset_t > offset;
  ArrayRCP< const lno_t > localColumnIds;
  ArrayRCP< const scalar_t > values;
  matrix_->getAllValues(offset,localColumnIds,values);
  columnIds_.resize(nnz, 0);

  for(offset_t i = 0; i < offset[nrows]; i++){
    columnIds_[i] = colMap_->getGlobalElement(localColumnIds[i]);
  }

  if (nWeightsPerRow_ > 0){
    rowWeights_ = arcp(new input_t [nWeightsPerRow_], 0, nWeightsPerRow_, true);
    numNzWeight_ = arcp(new bool [nWeightsPerRow_], 0, nWeightsPerRow_, true);
    for (int i=0; i < nWeightsPerRow_; i++)
      numNzWeight_[i] = false;
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  void XpetraCrsMatrixAdapter<User,UserCoord>::setWeights(
    const scalar_t *weightVal, int stride, int idx)
{
  if (this->getPrimaryEntityType() == MATRIX_ROW)
    setRowWeights(weightVal, stride, idx);
  else {
    // TODO:  Need to allow weights for columns and/or nonzeros
    std::ostringstream emsg;
    emsg << __FILE__ << "," << __LINE__
         << " error:  setWeights not yet supported for"
         << " columns or nonzeros."
         << std::endl;
    throw std::runtime_error(emsg.str());
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  void XpetraCrsMatrixAdapter<User,UserCoord>::setRowWeights(
    const scalar_t *weightVal, int stride, int idx)
{
  typedef StridedData<lno_t,scalar_t> input_t;
  if(idx<0 || idx >= nWeightsPerRow_)
  {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid row weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
  }

  size_t nvtx = getLocalNumRows();
  ArrayRCP<const scalar_t> weightV(weightVal, 0, nvtx*stride, false);
  rowWeights_[idx] = input_t(weightV, stride);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  void XpetraCrsMatrixAdapter<User,UserCoord>::setWeightIsDegree(
    int idx)
{
  if (this->getPrimaryEntityType() == MATRIX_ROW)
    setRowWeightIsNumberOfNonZeros(idx);
  else {
    // TODO:  Need to allow weights for columns and/or nonzeros
    std::ostringstream emsg;
    emsg << __FILE__ << "," << __LINE__
         << " error:  setWeightIsNumberOfNonZeros not yet supported for"
         << " columns" << std::endl;
    throw std::runtime_error(emsg.str());
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  void XpetraCrsMatrixAdapter<User,UserCoord>::setRowWeightIsNumberOfNonZeros(
    int idx)
{
  if(idx<0 || idx >= nWeightsPerRow_)
  {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid row weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
  }


  numNzWeight_[idx] = true;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  template <typename Adapter>
    void XpetraCrsMatrixAdapter<User,UserCoord>::applyPartitioningSolution(
      const User &in, User *&out,
      const PartitioningSolution<Adapter> &solution) const
{
  // Get an import list (rows to be received)
  size_t numNewRows;
  ArrayRCP<gno_t> importList;
  try{
    numNewRows = Zoltan2::getImportList<Adapter,
                                        XpetraCrsMatrixAdapter<User,UserCoord> >
                                       (solution, this, importList);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Move the rows, creating a new matrix.
  RCP<User> outPtr = XpetraTraits<User>::doMigration(in, numNewRows,
                                                     importList.getRawPtr());
  out = const_cast<User *>(outPtr.get());
  outPtr.release();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  template <typename Adapter>
    void XpetraCrsMatrixAdapter<User,UserCoord>::applyPartitioningSolution(
      const User &in, RCP<User> &out,
      const PartitioningSolution<Adapter> &solution) const
{
  // Get an import list (rows to be received)
  size_t numNewRows;
  ArrayRCP<gno_t> importList;
  try{
    numNewRows = Zoltan2::getImportList<Adapter,
                                        XpetraCrsMatrixAdapter<User,UserCoord> >
                                       (solution, this, importList);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Move the rows, creating a new matrix.
  out = XpetraTraits<User>::doMigration(in, numNewRows,
                                        importList.getRawPtr());
}

}  //namespace Zoltan2

#endif
