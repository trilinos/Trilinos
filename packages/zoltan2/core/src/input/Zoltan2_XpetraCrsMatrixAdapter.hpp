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

#include <Zoltan2_MatrixAdapter.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_PartitioningHelpers.hpp>

#include <Xpetra_CrsMatrix.hpp>

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
  class XpetraCrsMatrixAdapter : public MatrixAdapter<User,UserCoord> {
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::part_t   part_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef typename InputTraits<User>::offset_t offset_t;
  typedef Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> xmatrix_t;
  typedef User user_t;
  typedef UserCoord userCoord_t;
#endif

  /*! \brief Destructor
   */
  ~XpetraCrsMatrixAdapter() { }

  /*! \brief Constructor   
   *    \param inmatrix The users Epetra, Tpetra, or Xpetra CrsMatrix object 
   *    \param nWeightsPerRow If row weights will be provided in setRowWeights(),
   *        the set \c nWeightsPerRow to the number of weights per row.
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
   *       theMatrix->getRowMap()->getNodeElementList();
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
    return matrix_->getNodeNumRows();
  }

  size_t getLocalNumColumns() const { 
    return matrix_->getNodeNumCols();
  }

  size_t getLocalNumEntries() const {
    return matrix_->getNodeNumEntries();
  }

  bool CRSViewAvailable() const { return true; }

  void getRowIDsView(const gno_t *&rowIds) const 
  {
    ArrayView<const gno_t> rowView = rowMap_->getNodeElementList();
    rowIds = rowView.getRawPtr();
  }

  void getCRSView(const offset_t *&offsets, const gno_t *&colIds) const
  {
    offsets = offset_.getRawPtr();
    colIds = columnIds_.getRawPtr();
  }

  void getCRSView(const offset_t *&offsets, const gno_t *&colIds,
                    const scalar_t *&values) const
  {
    offsets = offset_.getRawPtr();
    colIds = columnIds_.getRawPtr();
    values = values_.getRawPtr();
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
  ArrayRCP< const offset_t > offset_;
  ArrayRCP< const lno_t > localColumnIds_;
  ArrayRCP< const scalar_t > values_;
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
      offset_(), columnIds_(),
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

  size_t nrows = matrix_->getNodeNumRows();
  size_t nnz = matrix_->getNodeNumEntries();

  // Get ArrayRCP pointers to the structures in the underlying matrix
  matrix_->getAllValues(offset_,localColumnIds_,values_);
  columnIds_.resize(nnz, 0);

  for(offset_t i = 0; i < offset_[nrows]; i++) {
    columnIds_[i] = colMap_->getGlobalElement(localColumnIds_[i]);
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
