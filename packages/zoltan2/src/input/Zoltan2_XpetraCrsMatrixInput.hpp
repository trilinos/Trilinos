// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_XpetraCrsMatrixAdapter.hpp
    \brief Defines the XpetraCrsMatrixAdapter class.
*/

#ifndef _ZOLTAN2_XPETRACRSMATRIXADAPTER_HPP_
#define _ZOLTAN2_XPETRACRSMATRIXADAPTER_HPP_

#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

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
  class XpetraCrsMatrixAdapter : public MatrixAdapter<User> {
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> xmatrix_t;
  typedef MatrixAdapter<User> base_adapter_t;
  typedef User user_t;
#endif

  /*! \brief Destructor
   */
  ~XpetraCrsMatrixAdapter() { }

  /*! \brief Constructor   
   *    \param inmatrix The users Epetra, Tpetra, or Xpetra CrsMatrix object 
   *    \param coordDim Some algorithms can use row geometric
   *            information if it is available.  If coordinates will be
   *            supplied in setRowCoordinates() 
   *            then provide the dimension of the coordinates here.
   *    \param weightDim If row weights will be provided in setRowWeights(),
   *        the set \c weightDim to the number of weights per row.
   */
  XpetraCrsMatrixAdapter(const RCP<const User> &inmatrix, int coordDim=0, 
    int weightDim=0);

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

  /*! \brief Specify a weight for each row.
   *    \param dim  A value between zero and one less that the \c weightDim 
   *                  argument to the constructor.
   *    \param weightVal A pointer to the weights for this dimension.
   *    \stride          A stride to be used in reading the values.  The
   *        dimension \c dim weight for row \k should be found at
   *        <tt>weightVal[k*stride]</tt>.
   *
   * The order of weights should correspond to the order of rows
   * returned by
   *   \code
   *       theMatrix->getRowMap()->getNodeElementList();
   *   \endcode
   */

  void setRowWeights(int dim, const scalar_t *weightVal, int stride);

  /*! \brief Specify whether or not row weights for a dimension should be
              the count of row non zeros.
   *    \param dim If true, Zoltan2 will automatically us the number of
   *         non zeros in an row as the row's weight for dimension \c dim.
   */
  void setRowWeightIsNumberOfNonZeros(int dim);

  /*! \brief Access to Xpetra-wrapped user matrix. 
   */

  const RCP<const xmatrix_t> &getMatrix() const
  {
    return matrix_;
  }

  ////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumberOfObjects() const { return getLocalNumRows();}

  int getNumberOfWeightsPerObject() const { return 0;}

  size_t getObjectWeights(int dim, const scalar_t *&wgt, int &stride) const
  {
    return getRowWeights(dim, wgt, stride);
  }

  ////////////////////////////////////////////////////
  // The MatrixAdapter interface.
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

  int getRowWeightDimension() const
  {
    return weightDim_;
  }

  size_t getRowWeights(int dim,
     const scalar_t *&weights, int &stride) const
  {
    env_->localInputAssertion(__FILE__, __LINE__,
      "invalid weight dimension",
      dim >= 0 && dim < weightDim_, BASIC_ASSERTION);

    size_t length;
    rowWeights_[dim].getStridedList(length, weights, stride);
    return length;
  }

  bool getRowWeightIsNumberOfNonZeros(int dim) const { return numNzWeight_[dim];}

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
  // End of MatrixAdapter interface.
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

  int weightDim_;
  ArrayRCP<StridedData<lno_t, scalar_t> > rowWeights_;
  ArrayRCP<bool> numNzWeight_;

  bool mayHaveDiagonalEntries;
};

/////////////////////////////////////////////////////////////////
// Definitions
/////////////////////////////////////////////////////////////////

template <typename User>
  XpetraCrsMatrixAdapter<User>::XpetraCrsMatrixAdapter(
    const RCP<const User> &inmatrix, int coordDim, int weightDim):
      env_(rcp(new Environment)),
      inmatrix_(inmatrix), matrix_(), rowMap_(), colMap_(), base_(),
      offset_(), columnIds_(),
      coordinateDim_(coordDim), rowCoords_(),
      weightDim_(weightDim), rowWeights_(), numNzWeight_(),
      mayHaveDiagonalEntries(true)
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
    nnz = matrix_->getNumEntriesInLocalRow(row);
    matrix_->getLocalRowView(row, indices, nzs);
    for (size_t j=0; j < nnz; j++){
      // TODO - this will be slow
      //   Is it possible that global columns ids might be stored in order?
      columnIds_[next++] = colMap_->getGlobalElement(indices[j]);
    }
    offset_[i+1] = offset_[i] + nnz;
  } 

  if (coordinateDim_ > 0)
    rowCoords_ = arcp(new input_t [coordinateDim_], 0, coordinateDim_, true);

  if (weightDim_ > 0){
    rowWeights_ = arcp(new input_t [weightDim_], 0, weightDim_, true);
    numNzWeight_ = arcp(new bool [weightDim_], 0, weightDim_, true);
    for (int i=0; i < weightDim_; i++)
      numNzWeight_[i] = false;
  }
}

// TODO (from 3/21/12 mtg):  Consider changing interface to take an XpetraMultivector
template <typename User>
  void XpetraCrsMatrixAdapter<User>::setRowCoordinates(int dim,
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
  void XpetraCrsMatrixAdapter<User>::setRowWeights(int dim,
    const scalar_t *weightVal, int stride)
{
  typedef StridedData<lno_t,scalar_t> input_t;

  env_->localInputAssertion(__FILE__, __LINE__,
    "invalid row weight dimension",
    dim >= 0 && dim < weightDim_, BASIC_ASSERTION);

  size_t nvtx = getLocalNumRows();

  ArrayRCP<const scalar_t> weightV(weightVal, 0, nvtx, false);
  rowWeights_[dim] = input_t(weightV, stride);
}

template <typename User>
  void XpetraCrsMatrixAdapter<User>::setRowWeightIsNumberOfNonZeros(int dim)
{
  env_->localInputAssertion(__FILE__, __LINE__,
    "invalid row weight dimension",
    dim >= 0 && dim < weightDim_, BASIC_ASSERTION);

  numNzWeight_[dim] = true;
}

template <typename User>
  template <typename Adapter>
    size_t XpetraCrsMatrixAdapter<User>::applyPartitioningSolution(
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
    numNewRows = solution.convertSolutionToImportList(
      0, dummyIn, importList, dummyOut);
  }
  Z2_FORWARD_EXCEPTIONS;

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
