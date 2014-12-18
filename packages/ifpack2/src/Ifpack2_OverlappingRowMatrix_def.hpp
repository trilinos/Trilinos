/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_OVERLAPPINGROWMATRIX_DEF_HPP
#define IFPACK2_OVERLAPPINGROWMATRIX_DEF_HPP

#include <Ifpack2_OverlappingRowMatrix_decl.hpp>
#include <Ifpack2_Details_OverlappingRowGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>

namespace Ifpack2 {

template<class MatrixType>
OverlappingRowMatrix<MatrixType>::
OverlappingRowMatrix (const Teuchos::RCP<const row_matrix_type>& A,
                      const int overlapLevel) :
  A_ (A),
  OverlapLevel_ (overlapLevel),
  UseSubComm_ (false)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Array;
  using Teuchos::outArg;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  typedef Tpetra::global_size_t GST;
  typedef Tpetra::CrsGraph<local_ordinal_type,
                           global_ordinal_type, node_type> crs_graph_type;
  TEUCHOS_TEST_FOR_EXCEPTION(
    OverlapLevel_ <= 0, std::runtime_error,
    "Ifpack2::OverlappingRowMatrix: OverlapLevel must be > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_->getComm()->getSize() == 1, std::runtime_error,
    "Ifpack2::OverlappingRowMatrix: Matrix must be "
    "distributed over more than one MPI process.");

  RCP<const crs_matrix_type> ACRS =
    rcp_dynamic_cast<const crs_matrix_type, const row_matrix_type> (A_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    ACRS.is_null (), std::runtime_error,
    "Ifpack2::OverlappingRowMatrix: The input matrix must be a Tpetra::"
    "CrsMatrix with matching template parameters.  This class currently "
    "requires that CrsMatrix's fifth template parameter be the default.");
  RCP<const crs_graph_type> A_crsGraph = ACRS->getCrsGraph ();

  const size_t numMyRowsA = A_->getNodeNumRows ();
  const global_ordinal_type global_invalid =
    Teuchos::OrdinalTraits<global_ordinal_type>::invalid ();

  // Temp arrays
  Array<global_ordinal_type> ExtElements;
  RCP<map_type>        TmpMap;
  RCP<crs_graph_type>  TmpGraph;
  RCP<import_type>     TmpImporter;
  RCP<const map_type>  RowMap, ColMap;

  // The big import loop
  for (int overlap = 0 ; overlap < OverlapLevel_ ; ++overlap) {
    // Get the current maps
    if (overlap == 0) {
      RowMap = A_->getRowMap ();
      ColMap = A_->getColMap ();
    }
    else {
      RowMap = TmpGraph->getRowMap ();
      ColMap = TmpGraph->getColMap ();
    }

    const size_t size = ColMap->getNodeNumElements () - RowMap->getNodeNumElements ();
    Array<global_ordinal_type> mylist (size);
    size_t count = 0;

    // define the set of rows that are in ColMap but not in RowMap
    for (local_ordinal_type i = 0 ; (size_t) i < ColMap->getNodeNumElements() ; ++i) {
      const global_ordinal_type GID = ColMap->getGlobalElement (i);
      if (A_->getRowMap ()->getLocalElement (GID) == global_invalid) {
        typedef typename Array<global_ordinal_type>::iterator iter_type;
        const iter_type end = ExtElements.end ();
        const iter_type pos = std::find (ExtElements.begin (), end, GID);
        if (pos == end) {
          ExtElements.push_back (GID);
          mylist[count] = GID;
          ++count;
        }
      }
    }

    // mfh 24 Nov 2013: We don't need TmpMap, TmpGraph, or
    // TmpImporter after this loop, so we don't have to construct them
    // on the last round.
    if (overlap + 1 < OverlapLevel_) {
      // Allocate & import new matrices, maps, etc.
      //
      // FIXME (mfh 24 Nov 2013) Do we always want to use index base
      // zero?  It doesn't really matter, since the actual index base
      // (in the current implementation of Map) will always be the
      // globally least GID.
      TmpMap = rcp (new map_type (global_invalid, mylist (0, count),
                                  Teuchos::OrdinalTraits<global_ordinal_type>::zero (),
                                  A_->getComm (), A_->getNode ()));
      TmpGraph = rcp (new crs_graph_type (TmpMap, 0));
      TmpImporter = rcp (new import_type (A_->getRowMap (), TmpMap));

      TmpGraph->doImport (*A_crsGraph, *TmpImporter, Tpetra::INSERT);
      TmpGraph->fillComplete (A_->getDomainMap (), TmpMap);
    }
  }

  // build the map containing all the nodes (original
  // matrix + extended matrix)
  Array<global_ordinal_type> mylist (numMyRowsA + ExtElements.size ());
  for (local_ordinal_type i = 0; (size_t)i < numMyRowsA; ++i) {
    mylist[i] = A_->getRowMap ()->getGlobalElement (i);
  }
  for (local_ordinal_type i = 0; i < ExtElements.size (); ++i) {
    mylist[i + numMyRowsA] = ExtElements[i];
  }

  RowMap_ = rcp (new map_type (global_invalid, mylist (),
                               Teuchos::OrdinalTraits<global_ordinal_type>::zero (),
                               A_->getComm (), A_->getNode ()));
  ColMap_ = RowMap_;

  // now build the map corresponding to all the external nodes
  // (with respect to A().RowMatrixRowMap().
  ExtMap_ = rcp (new map_type (global_invalid, ExtElements (),
                               Teuchos::OrdinalTraits<global_ordinal_type>::zero (),
                               A_->getComm (), A_->getNode ()));
  ExtMatrix_ = rcp (new crs_matrix_type (ExtMap_, ColMap_, 0));
  ExtImporter_ = rcp (new import_type (A_->getRowMap (), ExtMap_));

  RCP<crs_matrix_type> ExtMatrixCRS =
    rcp_dynamic_cast<crs_matrix_type, row_matrix_type> (ExtMatrix_);
  ExtMatrixCRS->doImport (*ACRS, *ExtImporter_, Tpetra::INSERT);
  ExtMatrixCRS->fillComplete (A_->getDomainMap (), RowMap_);

  Importer_ = rcp (new import_type (A_->getRowMap (), RowMap_));

  // fix indices for overlapping matrix
  const size_t numMyRowsB = ExtMatrix_->getNodeNumRows ();

  GST NumMyNonzeros_tmp = A_->getNodeNumEntries () + ExtMatrix_->getNodeNumEntries ();
  GST NumMyRows_tmp = numMyRowsA + numMyRowsB;
  {
    GST inArray[2], outArray[2];
    inArray[0] = NumMyNonzeros_tmp;
    inArray[1] = NumMyRows_tmp;
    outArray[0] = 0;
    outArray[1] = 0;
    reduceAll<int, GST> (* (A_->getComm ()), REDUCE_SUM, 2, inArray, outArray);
    NumGlobalNonzeros_ = outArray[0];
    NumGlobalRows_ = outArray[1];
  }
  // reduceAll<int, GST> (* (A_->getComm ()), REDUCE_SUM, NumMyNonzeros_tmp,
  //                      outArg (NumGlobalNonzeros_));
  // reduceAll<int, GST> (* (A_->getComm ()), REDUCE_SUM, NumMyRows_tmp,
  //                      outArg (NumGlobalRows_));

  MaxNumEntries_ = A_->getNodeMaxNumRowEntries ();
  if (MaxNumEntries_ < ExtMatrix_->getNodeMaxNumRowEntries ()) {
    MaxNumEntries_ = ExtMatrix_->getNodeMaxNumRowEntries ();
  }

  // Create the graph (returned by getGraph()).
  typedef Details::OverlappingRowGraph<row_graph_type> row_graph_impl_type;
  RCP<row_graph_impl_type> graph =
    rcp (new row_graph_impl_type (A_->getGraph (),
                                  ExtMatrix_->getGraph (),
                                  RowMap_,
                                  ColMap_,
                                  NumGlobalRows_,
                                  NumGlobalRows_, // # global cols == # global rows
                                  NumGlobalNonzeros_,
                                  MaxNumEntries_,
                                  rcp_const_cast<const import_type> (Importer_),
                                  rcp_const_cast<const import_type> (ExtImporter_)));
  graph_ = rcp_const_cast<const row_graph_type> (rcp_implicit_cast<row_graph_type> (graph));
  // Resize temp arrays
  Indices_.resize (MaxNumEntries_);
  Values_.resize (MaxNumEntries_);
}


template<class MatrixType>
OverlappingRowMatrix<MatrixType>::
OverlappingRowMatrix (const Teuchos::RCP<const row_matrix_type>& A,
                      const int overlapLevel,
                      const int subdomainID)
{
  //FIXME
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error,
    "Ifpack2::OverlappingRowMatrix: Subdomain code not implemented yet.");
}


template<class MatrixType>
OverlappingRowMatrix<MatrixType>::~OverlappingRowMatrix() {}


template<class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
OverlappingRowMatrix<MatrixType>::getComm () const
{
  return A_->getComm ();
}


template<class MatrixType>
Teuchos::RCP<typename MatrixType::node_type>
OverlappingRowMatrix<MatrixType>::getNode () const
{
  return A_->getNode();
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> >
OverlappingRowMatrix<MatrixType>::getRowMap () const
{
  // FIXME (mfh 12 July 2013) Is this really the right Map to return?
  return RowMap_;
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> >
OverlappingRowMatrix<MatrixType>::getColMap () const
{
  // FIXME (mfh 12 July 2013) Is this really the right Map to return?
  return ColMap_;
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> >
OverlappingRowMatrix<MatrixType>::getDomainMap () const
{
  return A_->getDomainMap();
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> >
OverlappingRowMatrix<MatrixType>::getRangeMap () const
{
  return A_->getRangeMap ();
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::RowGraph<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> >
OverlappingRowMatrix<MatrixType>::getGraph() const
{
  return graph_;
}


template<class MatrixType>
global_size_t OverlappingRowMatrix<MatrixType>::getGlobalNumRows() const
{
  return NumGlobalRows_;
}


template<class MatrixType>
global_size_t OverlappingRowMatrix<MatrixType>::getGlobalNumCols() const
{
  return NumGlobalRows_;
}


template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getNodeNumRows() const
{
  return A_->getNodeNumRows () + ExtMatrix_->getNodeNumRows ();
}


template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getNodeNumCols() const
{
  return this->getNodeNumRows ();
}


template<class MatrixType>
typename MatrixType::global_ordinal_type
OverlappingRowMatrix<MatrixType>::getIndexBase () const
{
  return A_->getIndexBase();
}


template<class MatrixType>
Tpetra::global_size_t OverlappingRowMatrix<MatrixType>::getGlobalNumEntries() const
{
  return NumGlobalNonzeros_;
}


template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getNodeNumEntries() const
{
  return A_->getNodeNumEntries () + ExtMatrix_->getNodeNumEntries ();
}


template<class MatrixType>
size_t
OverlappingRowMatrix<MatrixType>::
getNumEntriesInGlobalRow (global_ordinal_type globalRow) const
{
  const local_ordinal_type localRow = RowMap_->getLocalElement (globalRow);
  if (localRow == Teuchos::OrdinalTraits<local_ordinal_type>::invalid ()) {
    return Teuchos::OrdinalTraits<size_t>::invalid();
  } else {
    return getNumEntriesInLocalRow (localRow);
  }
}


template<class MatrixType>
size_t
OverlappingRowMatrix<MatrixType>::
getNumEntriesInLocalRow (local_ordinal_type localRow) const
{
  using Teuchos::as;
  const size_t numMyRowsA = A_->getNodeNumRows ();
  if (as<size_t> (localRow) < numMyRowsA) {
    return A_->getNumEntriesInLocalRow (localRow);
  } else {
    return ExtMatrix_->getNumEntriesInLocalRow (as<local_ordinal_type> (localRow - numMyRowsA));
  }
}


template<class MatrixType>
global_size_t OverlappingRowMatrix<MatrixType>::getGlobalNumDiags() const
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix::getGlobalNumDiags() not supported.");
}


template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getNodeNumDiags() const
{
  return A_->getNodeNumDiags();
}


template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getGlobalMaxNumRowEntries() const
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix::getGlobalMaxNumRowEntries() not supported.");
}


template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getNodeMaxNumRowEntries() const
{
  return MaxNumEntries_;
}


template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::hasColMap() const
{
  return true;
}


template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::isLowerTriangular() const
{
  return A_->isLowerTriangular();
}


template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::isUpperTriangular() const
{
  return A_->isUpperTriangular();
}


template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::isLocallyIndexed() const
{
  return true;
}


template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::isGloballyIndexed() const
{
  return false;
}


template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::isFillComplete() const
{
  return true;
}


template<class MatrixType>
void
OverlappingRowMatrix<MatrixType>::
getGlobalRowCopy (global_ordinal_type GlobalRow,
                  const Teuchos::ArrayView<global_ordinal_type> &Indices,
                  const Teuchos::ArrayView<scalar_type>& Values,
                  size_t& NumEntries) const
{
  const local_ordinal_type LocalRow = RowMap_->getLocalElement (GlobalRow);
  if (LocalRow == Teuchos::OrdinalTraits<local_ordinal_type>::invalid ()) {
    NumEntries = Teuchos::OrdinalTraits<size_t>::invalid ();
  } else {
    if (Teuchos::as<size_t> (LocalRow) < A_->getNodeNumRows ()) {
      A_->getGlobalRowCopy (GlobalRow, Indices, Values, NumEntries);
    } else {
      ExtMatrix_->getGlobalRowCopy (GlobalRow, Indices, Values, NumEntries);
    }
  }
}


template<class MatrixType>
void
OverlappingRowMatrix<MatrixType>::
getLocalRowCopy (local_ordinal_type LocalRow,
                 const Teuchos::ArrayView<local_ordinal_type> &Indices,
                 const Teuchos::ArrayView<scalar_type> &Values,
                 size_t &NumEntries) const
{
  using Teuchos::as;
  const size_t numMyRowsA = A_->getNodeNumRows ();
  if (as<size_t> (LocalRow) < numMyRowsA) {
    A_->getLocalRowCopy (LocalRow, Indices, Values, NumEntries);
  } else {
    ExtMatrix_->getLocalRowCopy (LocalRow - as<local_ordinal_type> (numMyRowsA),
                                 Indices, Values, NumEntries);
  }
}


template<class MatrixType>
void
OverlappingRowMatrix<MatrixType>::
getGlobalRowView (global_ordinal_type GlobalRow,
                  Teuchos::ArrayView<const global_ordinal_type>& indices,
                  Teuchos::ArrayView<const scalar_type>& values) const
{
  const local_ordinal_type LocalRow = RowMap_->getLocalElement (GlobalRow);
  if (LocalRow == Teuchos::OrdinalTraits<local_ordinal_type>::invalid())  {
    indices = Teuchos::null;
    values = Teuchos::null;
  } else {
    if (Teuchos::as<size_t> (LocalRow) < A_->getNodeNumRows ()) {
      A_->getGlobalRowView (GlobalRow, indices, values);
    } else {
      ExtMatrix_->getGlobalRowView (GlobalRow, indices, values);
    }
  }
}


template<class MatrixType>
void
OverlappingRowMatrix<MatrixType>::
getLocalRowView (local_ordinal_type LocalRow,
                 Teuchos::ArrayView<const local_ordinal_type>& indices,
                 Teuchos::ArrayView<const scalar_type>& values) const
{
  using Teuchos::as;
  const size_t numMyRowsA = A_->getNodeNumRows ();
  if (as<size_t> (LocalRow) < numMyRowsA) {
    A_->getLocalRowView (LocalRow, indices, values);
  } else {
    ExtMatrix_->getLocalRowView (LocalRow - as<local_ordinal_type> (numMyRowsA),
                                 indices, values);
  }
}


template<class MatrixType>
void
OverlappingRowMatrix<MatrixType>::
getLocalDiagCopy (Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& diag) const
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix::getLocalDiagCopy not supported.");
}


template<class MatrixType>
void
OverlappingRowMatrix<MatrixType>::
leftScale (const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& x)
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix does not support leftScale.");
}


template<class MatrixType>
void
OverlappingRowMatrix<MatrixType>::
rightScale (const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& x)
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix does not support leftScale.");
}


template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
OverlappingRowMatrix<MatrixType>::getFrobeniusNorm () const
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix does not support getFrobeniusNorm.");
}


template<class MatrixType>
void
OverlappingRowMatrix<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::ArrayRCP;
  using Teuchos::as;
  typedef scalar_type RangeScalar;
  typedef scalar_type DomainScalar;

  typedef Teuchos::ScalarTraits<RangeScalar> STRS;
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
    "Ifpack2::OverlappingRowMatrix::apply: The input X and the output Y must "
    "have the same number of columns.  X.getNumVectors() = "
    << X.getNumVectors() << " != Y.getNumVectors() = " << Y.getNumVectors()
    << ".");

  // FIXME (mfh 13 July 2013) This would be a good candidate for a
  // local parallel operator implementation.  That would obviate the
  // need for getting views of the data and make the code below a lot
  // simpler.

  const RangeScalar zero = STRS::zero ();
  ArrayRCP<ArrayRCP<const DomainScalar> > x_ptr = X.get2dView();
  ArrayRCP<ArrayRCP<RangeScalar> >        y_ptr = Y.get2dViewNonConst();
  Y.putScalar(zero);
  size_t NumVectors = Y.getNumVectors();

  const size_t numMyRowsA = A_->getNodeNumRows ();
  for (size_t i = 0; i < numMyRowsA; ++i) {
    size_t Nnz;
    // Use this class's getrow to make the below code simpler
    A_->getLocalRowCopy (i, Indices_ (),Values_ (), Nnz);
    if (mode == Teuchos::NO_TRANS) {
      for (size_t j = 0; j < Nnz; ++j)
        for (size_t k = 0; k < NumVectors; ++k)
          y_ptr[k][i] += as<RangeScalar> (Values_[j]) *
            as<RangeScalar> (x_ptr[k][Indices_[j]]);
    }
    else if (mode == Teuchos::TRANS){
      for (size_t j = 0; j < Nnz; ++j)
        for (size_t k = 0; k < NumVectors; ++k)
          y_ptr[k][Indices_[j]] += as<RangeScalar> (Values_[j]) *
            as<RangeScalar> (x_ptr[k][i]);
    }
    else { // mode == Teuchos::CONJ_TRANS
      for (size_t j = 0; j < Nnz; ++j)
        for (size_t k = 0; k < NumVectors; ++k)
          y_ptr[k][Indices_[j]] +=
            STRS::conjugate (as<RangeScalar> (Values_[j])) *
            as<RangeScalar> (x_ptr[k][i]);
    }
  }

  const size_t numMyRowsB = ExtMatrix_->getNodeNumRows ();
  for (size_t i = 0 ; i < numMyRowsB ; ++i) {
    size_t Nnz;
    // Use this class's getrow to make the below code simpler
    ExtMatrix_->getLocalRowCopy (i, Indices_ (), Values_ (), Nnz);
    if (mode == Teuchos::NO_TRANS) {
      for (size_t j = 0; j < Nnz; ++j)
        for (size_t k = 0; k < NumVectors; ++k)
          y_ptr[k][numMyRowsA+i] += as<RangeScalar> (Values_[j]) *
            as<RangeScalar> (x_ptr[k][Indices_[j]]);
    }
    else if (mode == Teuchos::TRANS) {
      for (size_t j = 0; j < Nnz; ++j)
        for (size_t k = 0; k < NumVectors; ++k)
          y_ptr[k][numMyRowsA+Indices_[j]] += as<RangeScalar> (Values_[j]) *
            as<RangeScalar> (x_ptr[k][i]);
    }
    else { // mode == Teuchos::CONJ_TRANS
      for (size_t j = 0; j < Nnz; ++j)
        for (size_t k = 0; k < NumVectors; ++k)
          y_ptr[k][numMyRowsA+Indices_[j]] +=
            STRS::conjugate (as<RangeScalar> (Values_[j])) *
            as<RangeScalar> (x_ptr[k][i]);
    }
  }
}


template<class MatrixType>
void
OverlappingRowMatrix<MatrixType>::
importMultiVector (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
                   Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &OvX,
                   Tpetra::CombineMode CM)
{
  OvX.doImport (X, *Importer_, CM);
}


template<class MatrixType>
void
OverlappingRowMatrix<MatrixType>::
exportMultiVector (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &OvX,
                   Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
                   Tpetra::CombineMode CM)
{
  X.doExport (OvX, *Importer_, CM);
}


template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::hasTransposeApply () const
{
  return true;
}


template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::supportsRowViews () const
{
  return false;
}

template<class MatrixType>
std::string OverlappingRowMatrix<MatrixType>::description() const
{
  std::ostringstream oss;
  if (isFillComplete()) {
    oss << "{ isFillComplete: true"
        << ", global rows: " << getGlobalNumRows()
        << ", global columns: " << getGlobalNumCols()
        << ", global entries: " << getGlobalNumEntries()
        << " }";
  }
  else {
    oss << "{ isFillComplete: false"
        << ", global rows: " << getGlobalNumRows()
        << " }";
  }
  return oss.str();
}

template<class MatrixType>
void OverlappingRowMatrix<MatrixType>::describe(Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
{
    using std::endl;
    using std::setw;
    using Teuchos::as;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    using Teuchos::RCP;
    using Teuchos::null;
    using Teuchos::ArrayView;

    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) {
      vl = VERB_LOW;
    }
    RCP<const Teuchos::Comm<int> > comm = this->getComm();
    const int myRank = comm->getRank();
    const int numProcs = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t> (width, as<size_t> (11)) + 2;
    Teuchos::OSTab tab(out);
    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per process
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print indices and values
    //
    // for medium and higher, print constituent objects at specified verbLevel
    if (vl != VERB_NONE) {
      if (myRank == 0) {
        out << this->description() << std::endl;
      }
      // O(1) globals, minus what was already printed by description()
      //if (isFillComplete() && myRank == 0) {
      //  out << "Global number of diagonal entries: " << getGlobalNumDiags() << std::endl;
      //  out << "Global max number of entries in a row: " << getGlobalMaxNumRowEntries() << std::endl;
      //}
      // constituent objects
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        if (myRank == 0) {
          out << endl << "Row map:" << endl;
        }
        getRowMap()->describe(out,vl);
        //
        if (getColMap() != null) {
          if (getColMap() == getRowMap()) {
            if (myRank == 0) {
              out << endl << "Column map is row map.";
            }
          }
          else {
            if (myRank == 0) {
              out << endl << "Column map:" << endl;
            }
            getColMap()->describe(out,vl);
          }
        }
        if (getDomainMap() != null) {
          if (getDomainMap() == getRowMap()) {
            if (myRank == 0) {
              out << endl << "Domain map is row map.";
            }
          }
          else if (getDomainMap() == getColMap()) {
            if (myRank == 0) {
              out << endl << "Domain map is column map.";
            }
          }
          else {
            if (myRank == 0) {
              out << endl << "Domain map:" << endl;
            }
            getDomainMap()->describe(out,vl);
          }
        }
        if (getRangeMap() != null) {
          if (getRangeMap() == getDomainMap()) {
            if (myRank == 0) {
              out << endl << "Range map is domain map." << endl;
            }
          }
          else if (getRangeMap() == getRowMap()) {
            if (myRank == 0) {
              out << endl << "Range map is row map." << endl;
            }
          }
          else {
            if (myRank == 0) {
              out << endl << "Range map: " << endl;
            }
            getRangeMap()->describe(out,vl);
          }
        }
        if (myRank == 0) {
          out << endl;
        }
      }
      // O(P) data
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int curRank = 0; curRank < numProcs; ++curRank) {
          if (myRank == curRank) {
            out << "Process rank: " << curRank << std::endl;
            out << "  Number of entries: " << getNodeNumEntries() << std::endl;
            if (isFillComplete()) {
              out << "  Number of diagonal entries: " << getNodeNumDiags() << std::endl;
            }
            out << "  Max number of entries per row: " << getNodeMaxNumRowEntries() << std::endl;
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
      // O(N) and O(NNZ) data
      if (vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int curRank = 0; curRank < numProcs; ++curRank) {
          if (myRank == curRank) {
            out << std::setw(width) << "Proc Rank"
                << std::setw(width) << "Global Row"
                << std::setw(width) << "Num Entries";
            if (vl == VERB_EXTREME) {
              out << std::setw(width) << "(Index,Value)";
            }
            out << endl;
            for (size_t r = 0; r < getNodeNumRows (); ++r) {
              const size_t nE = getNumEntriesInLocalRow(r);
              typename MatrixType::global_ordinal_type gid = getRowMap()->getGlobalElement(r);
              out << std::setw(width) << myRank
                  << std::setw(width) << gid
                  << std::setw(width) << nE;
              if (vl == VERB_EXTREME) {
                if (isGloballyIndexed()) {
                  ArrayView<const typename MatrixType::global_ordinal_type> rowinds;
                  ArrayView<const typename MatrixType::scalar_type> rowvals;
                  getGlobalRowView (gid, rowinds, rowvals);
                  for (size_t j = 0; j < nE; ++j) {
                    out << " (" << rowinds[j]
                        << ", " << rowvals[j]
                        << ") ";
                  }
                }
                else if (isLocallyIndexed()) {
                  ArrayView<const typename MatrixType::local_ordinal_type> rowinds;
                  ArrayView<const typename MatrixType::scalar_type> rowvals;
                  getLocalRowView (r, rowinds, rowvals);
                  for (size_t j=0; j < nE; ++j) {
                    out << " (" << getColMap()->getGlobalElement(rowinds[j])
                        << ", " << rowvals[j]
                        << ") ";
                  }
                } // globally or locally indexed
              } // vl == VERB_EXTREME
              out << endl;
            } // for each row r on this process

          } // if (myRank == curRank)
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }

        out.setOutputToRootOnly(0);
        out << "===========\nlocal matrix\n=================" << std::endl;
        out.setOutputToRootOnly(-1);
        A_->describe(out,Teuchos::VERB_EXTREME);
        out.setOutputToRootOnly(0);
        out << "===========\nend of local matrix\n=================" << std::endl;
        comm->barrier();
        out.setOutputToRootOnly(0);
        out << "=================\nghost matrix\n=================" << std::endl;
        out.setOutputToRootOnly(-1);
        ExtMatrix_->describe(out,Teuchos::VERB_EXTREME);
        out.setOutputToRootOnly(0);
        out << "===========\nend of ghost matrix\n=================" << std::endl;
      }
    }
}


} // namespace Ifpack2

#define IFPACK2_OVERLAPPINGROWMATRIX_INSTANT(S,LO,GO,N)                 \
  template class Ifpack2::OverlappingRowMatrix< Tpetra::CrsMatrix<S, LO, GO, N> >;

#endif // IFPACK2_OVERLAPPINGROWMATRIX_DEF_HPP
