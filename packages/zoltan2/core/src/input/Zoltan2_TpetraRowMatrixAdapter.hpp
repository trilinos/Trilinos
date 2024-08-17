// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_TpetraRowMatrixAdapter.hpp
    \brief Defines the TpetraRowMatrixAdapter class.
*/

#ifndef _ZOLTAN2_TPETRAROWMATRIXADAPTER_HPP_
#define _ZOLTAN2_TPETRAROWMATRIXADAPTER_HPP_

#include <Zoltan2_MatrixAdapter.hpp>
#include <Zoltan2_PartitioningHelpers.hpp>
#include <Zoltan2_StridedData.hpp>

#include <Tpetra_RowMatrix.hpp>

#include <vector>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////////////
/*!  \brief Provides access for Zoltan2 to Tpetra::RowMatrix data.

    The \c scalar_t type, representing user data such as matrix values, is
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::RowMatrix) have an inherent scalar type,
    and some
    (like Tpetra::RowGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.

*/

template <typename User, typename UserCoord = User>
class TpetraRowMatrixAdapter : public MatrixAdapter<User, UserCoord> {
public:


#ifndef DOXYGEN_SHOULD_SKIP_THIS
  using scalar_t = typename InputTraits<User>::scalar_t;
  using offset_t = typename InputTraits<User>::offset_t;
  using lno_t = typename InputTraits<User>::lno_t;
  using gno_t = typename InputTraits<User>::gno_t;
  using part_t = typename InputTraits<User>::part_t;
  using node_t = typename InputTraits<User>::node_t;
  using device_t = typename node_t::device_type;
  using host_t = typename Kokkos::HostSpace::memory_space;
  using user_t = User;
  using userCoord_t = UserCoord;

  using Base = MatrixAdapter<User, UserCoord>;
#endif

  /*! \brief Constructor
   *    \param inmatrix The user's Tpetra RowMatrix object
   *    \param nWeightsPerRow If row weights will be provided in setRowWeights(),
   *        then set \c nWeightsPerRow to the number of weights per row.
   */
  TpetraRowMatrixAdapter(const RCP<const User> &inmatrix,
                         int nWeightsPerRow = 0);

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

  /*! \brief Provide a device view of weights for the primary entity type.
   *    \param val A view to the weights for index \c idx.
   *    \param idx A number from 0 to one less than
   *          weight idx specified in the constructor.
   *
   *  The order of the weights should match the order that
   *  entities appear in the input data structure.
   */

  void setWeightsDevice(typename Base::ConstWeightsDeviceView1D val, int idx);

  /*! \brief Provide a host view of weights for the primary entity type.
   *    \param val A view to the weights for index \c idx.
   *    \param idx A number from 0 to one less than
   *          weight idx specified in the constructor.
   *
   *  The order of the weights should match the order that
   *  entities appear in the input data structure.
   */

  void setWeightsHost(typename Base::ConstWeightsHostView1D val, int idx);

  /*! \brief Specify a weight for each row.
   *    \param weightVal A pointer to the weights for this index.
   *    \param stride    A stride to be used in reading the values.
   *                     The index \c idx weight for row \k should be found at
   *                     <tt>weightVal[k*stride]</tt>.
   *    \param idx  A value between zero and one less that the \c nWeightsPerRow
   *                  argument to the constructor.
   *
   * The order of weights should correspond to the order of rows
   * returned by
   *   \code
   *       TheMatrix->getRowMap()->getLocalElementList();
   *   \endcode
   */

  void setRowWeights(const scalar_t *weightVal, int stride, int idx = 0);

  /*! \brief Provide a device view to row weights.
   *    \param val A pointer to the weights for index \c idx.
   *    \param idx A number from 0 to one less than
   *          number of row weights specified in the constructor.
   *
   *  The order of the row weights should match the order that
   *  rows appear in the input data structure.
   *     \code
   *       TheMatrix->getRowMap()->getLocalElementList()
   *     \endcode
   */

  void setRowWeightsDevice(typename Base::ConstWeightsDeviceView1D val,
                           int idx);

  /*! \brief Provide a host view to row weights.
   *    \param val A pointer to the weights for index \c idx.
   *    \param idx A number from 0 to one less than
   *          number of row weights specified in the constructor.
   *
   *  The order of the row weights should match the order that
   *  rows appear in the input data structure.
   *     \code
   *       TheMatrix->getRowMap()->getLocalElementList()
   *     \endcode
   */

  void setRowWeightsHost(typename Base::ConstWeightsHostView1D val,
                         int idx);


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

/////////////////////////////////////////////////////////////////
// The MatrixAdapter Interface
/////////////////////////////////////////////////////////////////

  size_t getLocalNumRows() const;

  size_t getLocalNumColumns() const;

  size_t getLocalNumEntries() const;

  bool CRSViewAvailable() const;

  void getRowIDsView(const gno_t *&rowIds) const override;

  void getRowIDsHostView(
      typename Base::ConstIdsHostView &rowIds) const override;

  void getRowIDsDeviceView(
      typename Base::ConstIdsDeviceView &rowIds) const override;

  void getCRSView(ArrayRCP<const offset_t> &offsets,
                  ArrayRCP<const gno_t> &colIds) const;

  void getCRSHostView(
      typename Base::ConstOffsetsHostView &offsets,
      typename Base::ConstIdsHostView &colIds) const override;

  void getCRSDeviceView(
      typename Base::ConstOffsetsDeviceView &offsets,
      typename Base::ConstIdsDeviceView &colIds) const override;

  void getCRSView(ArrayRCP<const offset_t> &offsets,
                  ArrayRCP<const gno_t> &colIds,
                  ArrayRCP<const scalar_t> &values) const;

  void getCRSHostView(
      typename Base::ConstOffsetsHostView &offsets,
      typename Base::ConstIdsHostView &colIds,
      typename Base::ConstScalarsHostView &values) const override;

  void getCRSDeviceView(
      typename Base::ConstOffsetsDeviceView &offsets,
      typename Base::ConstIdsDeviceView &colIds,
      typename Base::ConstScalarsDeviceView &values) const override;

  int getNumWeightsPerRow() const;

  void getRowWeightsView(const scalar_t *&weights, int &stride,
                         int idx = 0) const;

  void getRowWeightsDeviceView(typename Base::WeightsDeviceView1D &weights,
                                  int idx = 0) const;

  void getRowWeightsDeviceView(
      typename Base::WeightsDeviceView &weights) const override;

  void getRowWeightsHostView(typename Base::WeightsHostView1D &weights,
                                int idx = 0) const;

  void getRowWeightsHostView(
      typename Base::WeightsHostView &weights) const override;

  bool useNumNonzerosAsRowWeight(int idx) const;

  template <typename Adapter>
  void applyPartitioningSolution(
      const User &in, User *&out,
      const PartitioningSolution<Adapter> &solution) const;

  template <typename Adapter>
  void applyPartitioningSolution(
      const User &in, RCP<User> &out,
      const PartitioningSolution<Adapter> &solution) const;

protected:
  // Used by TpetraCrsMatrixAdapter
  TpetraRowMatrixAdapter(int nWeightsPerRow,
                         const RCP<const User> &inmatrix)
        : matrix_(inmatrix), nWeightsPerRow_(nWeightsPerRow) {}

  RCP<const User> matrix_;

  ArrayRCP<offset_t> offset_;
  ArrayRCP<gno_t> columnIds_;
  ArrayRCP<scalar_t> values_;

  typename Base::ConstOffsetsHostView offsHost_;
  typename Base::ConstIdsHostView colIdsHost_;
  typename Base::ScalarsHostView valuesHost_;

  typename Base::ConstOffsetsDeviceView offsDevice_;
  typename Base::ConstIdsDeviceView colIdsDevice_;
  typename Base::ScalarsDeviceView valuesDevice_;

  int nWeightsPerRow_;
  ArrayRCP<StridedData<lno_t, scalar_t>> rowWeights_;
  typename Base::WeightsDeviceView rowWeightsDevice_;
  Kokkos::View<bool *, host_t> numNzWeight_;

  bool mayHaveDiagonalEntries;

  virtual RCP<User> doMigration(const User &from, size_t numLocalRows,
                        const gno_t *myNewRows) const;
};

/////////////////////////////////////////////////////////////////
// Definitions
/////////////////////////////////////////////////////////////////

template <typename User, typename UserCoord>
TpetraRowMatrixAdapter<User, UserCoord>::TpetraRowMatrixAdapter(
    const RCP<const User> &inmatrix, int nWeightsPerRow):
      matrix_(inmatrix), offset_(), columnIds_(),
      nWeightsPerRow_(nWeightsPerRow), rowWeights_(),
      mayHaveDiagonalEntries(true) {
  using strided_t = StridedData<lno_t, scalar_t>;
  using localInds_t = typename User::nonconst_local_inds_host_view_type;
  using localVals_t = typename User::nonconst_values_host_view_type;

  const auto nrows = matrix_->getLocalNumRows();
  const auto nnz = matrix_->getLocalNumEntries();
  auto maxNumEntries = matrix_->getLocalMaxNumRowEntries();

  // Unfortunately we have to copy the offsets, column Ids, and vals
  // because column Ids are not usually stored in row id order.

  colIdsHost_ = typename Base::ConstIdsHostView("colIdsHost_", nnz);
  offsHost_ = typename Base::ConstOffsetsHostView("offsHost_", nrows + 1);
  valuesHost_ = typename Base::ScalarsHostView("valuesHost_", nnz);

  localInds_t localColInds("localColInds", maxNumEntries);
  localVals_t localVals("localVals", maxNumEntries);

  for (size_t r = 0; r < nrows; r++) {
    size_t numEntries = 0;
    matrix_->getLocalRowCopy(r, localColInds, localVals, numEntries); // Diff from CrsGraph

    offsHost_(r + 1) = offsHost_(r) + numEntries;
    for (offset_t e = offsHost_(r), i = 0; e < offsHost_(r + 1); e++) {
      colIdsHost_(e) = matrix_->getColMap()->getGlobalElement(localColInds(i++));
    }
    for (size_t j = 0; j < nnz; j++) {
      valuesHost_(r) = localVals[j];
    }
  }
  offsDevice_ = Kokkos::create_mirror_view_and_copy(
                        typename Base::device_t(), offsHost_);
  colIdsDevice_ = Kokkos::create_mirror_view_and_copy(
                        typename Base::device_t(), colIdsHost_);
  valuesDevice_ = Kokkos::create_mirror_view_and_copy(
                        typename Base::device_t(), valuesHost_);

  if (nWeightsPerRow_ > 0) {
    rowWeights_ =
        arcp(new strided_t[nWeightsPerRow_], 0, nWeightsPerRow_, true);

    rowWeightsDevice_ = typename Base::WeightsDeviceView(
        "rowWeightsDevice_", nrows, nWeightsPerRow_);

    numNzWeight_ =  Kokkos::View<bool *, host_t>(
        "numNzWeight_", nWeightsPerRow_);

    for (int i = 0; i < nWeightsPerRow_; ++i) {
      numNzWeight_(i) = false;
    }
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::setWeights(
    const scalar_t *weightVal, int stride, int idx) {
  if (this->getPrimaryEntityType() == MATRIX_ROW)
    setRowWeights(weightVal, stride, idx);
  else {
    // TODO:  Need to allow weights for columns and/or nonzeros
    std::ostringstream emsg;
    emsg << __FILE__ << "," << __LINE__
         << " error:  setWeights not yet supported for"
         << " columns or nonzeros." << std::endl;
    throw std::runtime_error(emsg.str());
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::setWeightsDevice(
    typename Base::ConstWeightsDeviceView1D val, int idx) {
  if (this->getPrimaryEntityType() == MATRIX_ROW)
    setRowWeightsDevice(val, idx);
  else {
    // TODO:  Need to allow weights for columns and/or nonzeros
    std::ostringstream emsg;
    emsg << __FILE__ << "," << __LINE__
         << " error:  setWeights not yet supported for"
         << " columns or nonzeros." << std::endl;
    throw std::runtime_error(emsg.str());
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::setWeightsHost(
    typename Base::ConstWeightsHostView1D val, int idx) {
  if (this->getPrimaryEntityType() == MATRIX_ROW)
    setRowWeightsHost(val, idx);
  else {
    // TODO:  Need to allow weights for columns and/or nonzeros
    std::ostringstream emsg;
    emsg << __FILE__ << "," << __LINE__
         << " error:  setWeights not yet supported for"
         << " columns or nonzeros." << std::endl;
    throw std::runtime_error(emsg.str());
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::setRowWeights(
    const scalar_t *weightVal, int stride, int idx) {
  typedef StridedData<lno_t, scalar_t> input_t;
  AssertCondition((idx >= 0) and (idx < nWeightsPerRow_),
                  "Invalid row weight index: " + std::to_string(idx));

  size_t nrows = getLocalNumRows();
  ArrayRCP<const scalar_t> weightV(weightVal, 0, nrows * stride, false);
  rowWeights_[idx] = input_t(weightV, stride);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::setRowWeightsDevice(
    typename Base::ConstWeightsDeviceView1D weights, int idx) {

  AssertCondition((idx >= 0) and (idx < nWeightsPerRow_),
                  "Invalid row weight index: " + std::to_string(idx));

    Kokkos::parallel_for(
      rowWeightsDevice_.extent(0), KOKKOS_CLASS_LAMBDA(const int rowID) {
        rowWeightsDevice_(rowID, idx) = weights(rowID);
      });

  Kokkos::fence();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::setRowWeightsHost(
    typename Base::ConstWeightsHostView1D weightsHost, int idx) {
  AssertCondition((idx >= 0) and (idx < nWeightsPerRow_),
                  "Invalid row weight index: " + std::to_string(idx));

  auto weightsDevice = Kokkos::create_mirror_view_and_copy(
      typename Base::device_t(), weightsHost);

  setRowWeightsDevice(weightsDevice, idx);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::setWeightIsDegree(int idx) {
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
void TpetraRowMatrixAdapter<User, UserCoord>::setRowWeightIsNumberOfNonZeros(
    int idx) {
  if (idx < 0 || idx >= nWeightsPerRow_) {
    std::ostringstream emsg;
    emsg << __FILE__ << ":" << __LINE__ << "  Invalid row weight index " << idx
         << std::endl;
    throw std::runtime_error(emsg.str());
  }

  numNzWeight_(idx) = true;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
size_t TpetraRowMatrixAdapter<User, UserCoord>::getLocalNumRows() const {
  return matrix_->getLocalNumRows();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
size_t TpetraRowMatrixAdapter<User, UserCoord>::getLocalNumColumns() const {
  return matrix_->getLocalNumCols();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
size_t TpetraRowMatrixAdapter<User, UserCoord>::getLocalNumEntries() const {
  return matrix_->getLocalNumEntries();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
bool TpetraRowMatrixAdapter<User, UserCoord>::CRSViewAvailable() const { return true; }

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getRowIDsView(const gno_t *&rowIds) const {
  ArrayView<const gno_t> rowView = matrix_->getRowMap()->getLocalElementList();
  rowIds = rowView.getRawPtr();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getRowIDsHostView(
    typename Base::ConstIdsHostView &rowIds) const {
  auto idsDevice = matrix_->getRowMap()->getMyGlobalIndices();
  auto tmpIds = typename Base::IdsHostView("", idsDevice.extent(0));

  Kokkos::deep_copy(tmpIds, idsDevice);

  rowIds = tmpIds;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getRowIDsDeviceView(
    typename Base::ConstIdsDeviceView &rowIds) const {

  auto idsDevice = matrix_->getRowMap()->getMyGlobalIndices();
  auto tmpIds = typename Base::IdsDeviceView("", idsDevice.extent(0));

  Kokkos::deep_copy(tmpIds, idsDevice);

  rowIds = tmpIds;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getCRSView(ArrayRCP<const offset_t> &offsets,
                ArrayRCP<const gno_t> &colIds) const {
  offsets = offset_;
  colIds = columnIds_;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getCRSHostView(
    typename Base::ConstOffsetsHostView &offsets,
    typename Base::ConstIdsHostView &colIds) const {
  auto hostOffsets = Kokkos::create_mirror_view(offsDevice_);
  Kokkos::deep_copy(hostOffsets, offsDevice_);
  offsets = hostOffsets;

  auto hostColIds = Kokkos::create_mirror_view(colIdsDevice_);
  Kokkos::deep_copy(hostColIds, colIdsDevice_);
  colIds = hostColIds;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getCRSDeviceView(
    typename Base::ConstOffsetsDeviceView &offsets,
    typename Base::ConstIdsDeviceView &colIds) const {
  offsets = offsDevice_;
  colIds = colIdsDevice_;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getCRSView(ArrayRCP<const offset_t> &offsets,
                ArrayRCP<const gno_t> &colIds,
                ArrayRCP<const scalar_t> &values) const {
  offsets = offset_;
  colIds = columnIds_;
  values = values_;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getCRSHostView(
    typename Base::ConstOffsetsHostView &offsets,
    typename Base::ConstIdsHostView &colIds,
    typename Base::ConstScalarsHostView &values) const {
  auto hostOffsets = Kokkos::create_mirror_view(offsDevice_);
  Kokkos::deep_copy(hostOffsets, offsDevice_);
  offsets = hostOffsets;

  auto hostColIds = Kokkos::create_mirror_view(colIdsDevice_);
  Kokkos::deep_copy(hostColIds, colIdsDevice_);
  colIds = hostColIds;

  auto hostValues = Kokkos::create_mirror_view(valuesDevice_);
  Kokkos::deep_copy(hostValues, valuesDevice_);
  values = hostValues;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getCRSDeviceView(
                  typename Base::ConstOffsetsDeviceView &offsets,
                  typename Base::ConstIdsDeviceView &colIds,
                  typename Base::ConstScalarsDeviceView &values) const {
  offsets = offsDevice_;
  colIds = colIdsDevice_;
  values = valuesDevice_;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
int TpetraRowMatrixAdapter<User, UserCoord>::getNumWeightsPerRow() const { return nWeightsPerRow_; }

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getRowWeightsView(const scalar_t *&weights, int &stride,
                        int idx) const {
  if (idx < 0 || idx >= nWeightsPerRow_) {
    std::ostringstream emsg;
    emsg << __FILE__ << ":" << __LINE__ << "  Invalid row weight index "
          << idx << std::endl;
    throw std::runtime_error(emsg.str());
  }

  size_t length;
  rowWeights_[idx].getStridedList(length, weights, stride);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getRowWeightsDeviceView(
    typename Base::WeightsDeviceView1D &weights, int idx) const {
  AssertCondition((idx >= 0) and (idx < nWeightsPerRow_),
                  "Invalid row weight index.");

  const auto size = rowWeightsDevice_.extent(0);
  weights = typename Base::WeightsDeviceView1D("weights", size);

  Kokkos::parallel_for(
      size, KOKKOS_CLASS_LAMBDA(const int id) {
        weights(id) = rowWeightsDevice_(id, idx);
      });

  Kokkos::fence();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getRowWeightsDeviceView(
    typename Base::WeightsDeviceView &weights) const {

  weights = rowWeightsDevice_;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getRowWeightsHostView(
    typename Base::WeightsHostView1D &weights, int idx) const {
  AssertCondition((idx >= 0) and (idx < nWeightsPerRow_),
                  "Invalid row weight index.");

  auto weightsDevice = typename Base::WeightsDeviceView1D(
      "weights", rowWeightsDevice_.extent(0));
  getRowWeightsDeviceView(weightsDevice, idx);

  weights = Kokkos::create_mirror_view(weightsDevice);
  Kokkos::deep_copy(weights, weightsDevice);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowMatrixAdapter<User, UserCoord>::getRowWeightsHostView(
    typename Base::WeightsHostView &weights) const {

  weights = Kokkos::create_mirror_view(rowWeightsDevice_);
  Kokkos::deep_copy(weights, rowWeightsDevice_);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
bool TpetraRowMatrixAdapter<User, UserCoord>::useNumNonzerosAsRowWeight(int idx) const { return numNzWeight_[idx]; }

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
template <typename Adapter>
void TpetraRowMatrixAdapter<User, UserCoord>::applyPartitioningSolution(
    const User &in, User *&out,
    const PartitioningSolution<Adapter> &solution) const {
  // Get an import list (rows to be received)
  size_t numNewRows;
  ArrayRCP<gno_t> importList;
  try {
    numNewRows =
        Zoltan2::getImportList<Adapter, TpetraRowMatrixAdapter<User, UserCoord>>(
            solution, this, importList);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Move the rows, creating a new matrix.
  RCP<User> outPtr = doMigration(in, numNewRows, importList.getRawPtr());
  out = outPtr.get();
  outPtr.release();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
template <typename Adapter>
void TpetraRowMatrixAdapter<User, UserCoord>::applyPartitioningSolution(
    const User &in, RCP<User> &out,
    const PartitioningSolution<Adapter> &solution) const {
  // Get an import list (rows to be received)
  size_t numNewRows;
  ArrayRCP<gno_t> importList;
  try {
    numNewRows =
        Zoltan2::getImportList<Adapter, TpetraRowMatrixAdapter<User, UserCoord>>(
            solution, this, importList);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Move the rows, creating a new matrix.
  out = doMigration(in, numNewRows, importList.getRawPtr());
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
RCP<User> TpetraRowMatrixAdapter<User, UserCoord>::doMigration(
    const User &from, size_t numLocalRows, const gno_t *myNewRows) const {
  typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;
  typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tcrsmatrix_t;

  // We cannot create a Tpetra::RowMatrix, unless the underlying type is
  // something we know (like Tpetra::CrsMatrix).
  // If the underlying type is something different, the user probably doesn't
  // want a Tpetra::CrsMatrix back, so we throw an error.

  // Try to cast "from" matrix to a TPetra::CrsMatrix
  // If that fails we throw an error.
  // We could cast as a ref which will throw std::bad_cast but with ptr
  // approach it might be clearer what's going on here
  const tcrsmatrix_t *pCrsMatrix = dynamic_cast<const tcrsmatrix_t *>(&from);

  if (!pCrsMatrix) {
    throw std::logic_error("TpetraRowMatrixAdapter cannot migrate data for "
                           "your RowMatrix; it can migrate data only for "
                           "Tpetra::CrsMatrix.  "
                           "You can inherit from TpetraRowMatrixAdapter and "
                           "implement migration for your RowMatrix.");
  }

  // source map
  const RCP<const map_t> &smap = from.getRowMap();
  gno_t numGlobalRows = smap->getGlobalNumElements();
  gno_t base = smap->getMinAllGlobalIndex();

  // target map
  ArrayView<const gno_t> rowList(myNewRows, numLocalRows);
  const RCP<const Teuchos::Comm<int>> &comm = from.getComm();
  RCP<const map_t> tmap = rcp(new map_t(numGlobalRows, rowList, base, comm));

  // importer
  Tpetra::Import<lno_t, gno_t, node_t> importer(smap, tmap);

  int oldNumElts = smap->getLocalNumElements();
  int newNumElts = numLocalRows;

  // number of non zeros in my new rows
  typedef Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> vector_t;
  vector_t numOld(smap); // TODO These vectors should have scalar=size_t,
  vector_t numNew(tmap); // but ETI does not yet support that.
  for (int lid = 0; lid < oldNumElts; lid++) {
    numOld.replaceGlobalValue(smap->getGlobalElement(lid),
                              scalar_t(from.getNumEntriesInLocalRow(lid)));
  }
  numNew.doImport(numOld, importer, Tpetra::INSERT);

  // TODO Could skip this copy if could declare vector with scalar=size_t.
  ArrayRCP<size_t> nnz(newNumElts);
  if (newNumElts > 0) {
    ArrayRCP<scalar_t> ptr = numNew.getDataNonConst(0);
    for (int lid = 0; lid < newNumElts; lid++) {
      nnz[lid] = static_cast<size_t>(ptr[lid]);
    }
  }

  RCP<tcrsmatrix_t> M = rcp(new tcrsmatrix_t(tmap, nnz()));

  M->doImport(from, importer, Tpetra::INSERT);
  M->fillComplete();

  return Teuchos::rcp_dynamic_cast<User>(M);
}

} // namespace Zoltan2

#endif // _ZOLTAN2_TPETRAROWMATRIXADAPTER_HPP_
