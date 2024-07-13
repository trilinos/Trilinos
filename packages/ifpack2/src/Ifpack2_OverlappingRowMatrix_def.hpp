// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_OVERLAPPINGROWMATRIX_DEF_HPP
#define IFPACK2_OVERLAPPINGROWMATRIX_DEF_HPP

#include <sstream>

#include <Ifpack2_OverlappingRowMatrix_decl.hpp>
#include <Ifpack2_Details_OverlappingRowGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Import.hpp>
#include "Tpetra_Map.hpp"
#include <Teuchos_CommHelpers.hpp>
#include <unordered_set>

namespace Ifpack2 {

template<class MatrixType>
OverlappingRowMatrix<MatrixType>::
OverlappingRowMatrix (const Teuchos::RCP<const row_matrix_type>& A,
                      const int overlapLevel) :
  A_ (Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A, true)),
  OverlapLevel_ (overlapLevel)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Array;
  using Teuchos::outArg;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  typedef Tpetra::global_size_t GST;
  typedef Tpetra::CrsGraph<local_ordinal_type,
                           global_ordinal_type, node_type> crs_graph_type;
  TEUCHOS_TEST_FOR_EXCEPTION(
    OverlapLevel_ <= 0, std::runtime_error,
    "Ifpack2::OverlappingRowMatrix: OverlapLevel must be > 0.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error,
     "Ifpack2::OverlappingRowMatrix: The input matrix must be a "
     "Tpetra::CrsMatrix with the same scalar_type, local_ordinal_type, "
     "global_ordinal_type, and device_type typedefs as MatrixType.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_->getComm()->getSize() == 1, std::runtime_error,
    "Ifpack2::OverlappingRowMatrix: Matrix must be "
    "distributed over more than one MPI process.");


  RCP<const crs_graph_type> A_crsGraph = A_->getCrsGraph ();
  const size_t numMyRowsA = A_->getLocalNumRows ();
  const global_ordinal_type global_invalid =
    Teuchos::OrdinalTraits<global_ordinal_type>::invalid ();

  // Temp arrays
  Array<global_ordinal_type> ExtElements;
  // Use an unordered_set to efficiently keep track of which GIDs have already
  // been added to ExtElements. Still need ExtElements because we also want a
  // list of the GIDs ordered LID in the ColMap.
  std::unordered_set<global_ordinal_type> ExtElementSet;
  RCP<map_type>        TmpMap;
  RCP<crs_graph_type>  TmpGraph;
  RCP<import_type>     TmpImporter;
  RCP<const map_type>  RowMap, ColMap;
  Kokkos::resize(ExtHaloStarts_, OverlapLevel_+1);
  ExtHaloStarts_h = Kokkos::create_mirror_view(ExtHaloStarts_);

  // The big import loop
  for (int overlap = 0 ; overlap < OverlapLevel_ ; ++overlap) {
    ExtHaloStarts_h(overlap) = (size_t) ExtElements.size();

    // Get the current maps
    if (overlap == 0) {
      RowMap = A_->getRowMap ();
      ColMap = A_->getColMap ();
    }
    else {
      RowMap = TmpGraph->getRowMap ();
      ColMap = TmpGraph->getColMap ();
    }

    const size_t size = ColMap->getLocalNumElements () - RowMap->getLocalNumElements ();
    Array<global_ordinal_type> mylist (size);
    size_t count = 0;

    // define the set of rows that are in ColMap but not in RowMap
    for (local_ordinal_type i = 0 ; (size_t) i < ColMap->getLocalNumElements() ; ++i) {
      const global_ordinal_type GID = ColMap->getGlobalElement (i);
      if (A_->getRowMap ()->getLocalElement (GID) == global_invalid) {
        // unordered_set insert can return a pair, where the second element is
        // true if a new element was inserted, false otherwise.
        if(ExtElementSet.insert(GID).second)
        {
          ExtElements.push_back (GID);
          mylist[count] = GID;
          ++count;
        }
      }
    }

    // On last import round, TmpMap, TmpGraph, and TmpImporter are unneeded,
    // so don't build them.
    if (overlap + 1 < OverlapLevel_) {
      //map consisting of GIDs that are in the current halo level-set
      TmpMap = rcp (new map_type (global_invalid, mylist (0, count),
                                  Teuchos::OrdinalTraits<global_ordinal_type>::zero (),
                                  A_->getComm ()));
      //graph whose rows are the current halo level-set to import
      TmpGraph = rcp (new crs_graph_type (TmpMap, 0));
      TmpImporter = rcp (new import_type (A_->getRowMap (), TmpMap));

      //import from original matrix graph to current halo level-set graph
      TmpGraph->doImport (*A_crsGraph, *TmpImporter, Tpetra::INSERT);
      TmpGraph->fillComplete (A_->getDomainMap (), TmpMap);
    }
  } // end overlap loop
  ExtHaloStarts_h[OverlapLevel_] = (size_t) ExtElements.size();
  Kokkos::deep_copy(ExtHaloStarts_,ExtHaloStarts_h);

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
                               A_->getComm ()));
  Importer_ = rcp (new import_type (A_->getRowMap (), RowMap_));
  ColMap_ = RowMap_;

  // now build the map corresponding to all the external nodes
  // (with respect to A().RowMatrixRowMap().
  ExtMap_ = rcp (new map_type (global_invalid, ExtElements (),
                               Teuchos::OrdinalTraits<global_ordinal_type>::zero (),
                               A_->getComm ()));
  ExtImporter_ = rcp (new import_type (A_->getRowMap (), ExtMap_));

  {
    auto ExtMatrixDynGraph = rcp (new crs_matrix_type (ExtMap_, ColMap_, 0));
    ExtMatrixDynGraph->doImport (*A_, *ExtImporter_, Tpetra::INSERT);
    ExtMatrixDynGraph->fillComplete (A_->getDomainMap (), RowMap_);
    auto ExtLclMatrix = ExtMatrixDynGraph->getLocalMatrixDevice();
    auto ExtMatrixStaticGraph = rcp (new crs_graph_type(ExtLclMatrix.graph,
      ExtMap_,
      ColMap_,
      ExtMatrixDynGraph->getDomainMap(),
      ExtMatrixDynGraph->getRangeMap()));
    ExtMatrix_ = rcp (new crs_matrix_type(ExtMatrixStaticGraph, ExtLclMatrix.values));
    ExtMatrix_->fillComplete ();
  }

  // fix indices for overlapping matrix
  const size_t numMyRowsB = ExtMatrix_->getLocalNumRows ();

  GST NumMyNonzeros_tmp = A_->getLocalNumEntries () + ExtMatrix_->getLocalNumEntries ();
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

  MaxNumEntries_ = A_->getLocalMaxNumRowEntries ();
  if (MaxNumEntries_ < ExtMatrix_->getLocalMaxNumRowEntries ()) {
    MaxNumEntries_ = ExtMatrix_->getLocalMaxNumRowEntries ();
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
                                  Importer_,
                                  ExtImporter_));
  graph_ = Teuchos::rcp_const_cast<const row_graph_type>
    (Teuchos::rcp_implicit_cast<row_graph_type> (graph));
  // Resize temp arrays
  Kokkos::resize(Indices_,MaxNumEntries_);
  Kokkos::resize(Values_,MaxNumEntries_);
}


template<class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
OverlappingRowMatrix<MatrixType>::getComm () const
{
  return A_->getComm ();
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
  // The original matrix's domain map is irrelevant; we want the map associated
  // with the overlap. This can then be used by LocalFilter, for example, while
  // letting LocalFilter still filter based on domain and range maps (instead of
  // column and row maps).
  // FIXME Ideally, this would be the same map but restricted to a local
  // communicator. If replaceCommWithSubset were free, that would be the way to
  // go. That would require a new Map ctor. For now, we'll stick with ColMap_'s
  // global communicator.
  return ColMap_;
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> >
OverlappingRowMatrix<MatrixType>::getRangeMap () const
{
  return RowMap_;
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
size_t OverlappingRowMatrix<MatrixType>::getLocalNumRows() const
{
  return A_->getLocalNumRows () + ExtMatrix_->getLocalNumRows ();
}


template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getLocalNumCols() const
{
  return this->getLocalNumRows ();
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
size_t OverlappingRowMatrix<MatrixType>::getLocalNumEntries() const
{
  return A_->getLocalNumEntries () + ExtMatrix_->getLocalNumEntries ();
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
  const size_t numMyRowsA = A_->getLocalNumRows ();
  if (as<size_t> (localRow) < numMyRowsA) {
    return A_->getNumEntriesInLocalRow (localRow);
  } else {
    return ExtMatrix_->getNumEntriesInLocalRow (as<local_ordinal_type> (localRow - numMyRowsA));
  }
}


template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getGlobalMaxNumRowEntries() const
{
  return std::max<size_t>(A_->getGlobalMaxNumRowEntries(), ExtMatrix_->getGlobalMaxNumRowEntries());
}


template<class MatrixType>
size_t OverlappingRowMatrix<MatrixType>::getLocalMaxNumRowEntries() const
{
  return MaxNumEntries_;
}

template<class MatrixType>
typename MatrixType::local_ordinal_type OverlappingRowMatrix<MatrixType>::getBlockSize() const
{
  return A_->getBlockSize();
}

template<class MatrixType>
bool OverlappingRowMatrix<MatrixType>::hasColMap() const
{
  return true;
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
                   nonconst_global_inds_host_view_type &Indices,
                   nonconst_values_host_view_type &Values,
                   size_t& NumEntries) const
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix::getGlobalRowCopy() not supported.");
}

template<class MatrixType>
void
OverlappingRowMatrix<MatrixType>::
  getLocalRowCopy (local_ordinal_type LocalRow,
                   nonconst_local_inds_host_view_type &Indices,
                   nonconst_values_host_view_type &Values,
                   size_t& NumEntries) const
{
  using Teuchos::as;
  const size_t numMyRowsA = A_->getLocalNumRows ();
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
                  global_inds_host_view_type &indices,
                  values_host_view_type &values) const {
  const local_ordinal_type LocalRow = RowMap_->getLocalElement (GlobalRow);
  if (LocalRow == Teuchos::OrdinalTraits<local_ordinal_type>::invalid())  {
    indices = global_inds_host_view_type();
    values = values_host_view_type();
  } else {
    if (Teuchos::as<size_t> (LocalRow) < A_->getLocalNumRows ()) {
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
                   local_inds_host_view_type & indices,
                   values_host_view_type & values) const {
  using Teuchos::as;
  const size_t numMyRowsA = A_->getLocalNumRows ();
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
  using Teuchos::Array;

  //extract diagonal of original matrix
  vector_type baseDiag(A_->getRowMap());         // diagonal of original matrix A_
  A_->getLocalDiagCopy(baseDiag);
  Array<scalar_type> baseDiagVals(baseDiag.getLocalLength());
  baseDiag.get1dCopy(baseDiagVals());
  //extra diagonal of ghost matrix
  vector_type extDiag(ExtMatrix_->getRowMap());
  ExtMatrix_->getLocalDiagCopy(extDiag);
  Array<scalar_type> extDiagVals(extDiag.getLocalLength());
  extDiag.get1dCopy(extDiagVals());

  Teuchos::ArrayRCP<scalar_type> allDiagVals = diag.getDataNonConst();
  if (allDiagVals.size() != baseDiagVals.size() + extDiagVals.size()) {
    std::ostringstream errStr;
    errStr << "Ifpack2::OverlappingRowMatrix::getLocalDiagCopy : Mismatch in diagonal lengths, "
           << allDiagVals.size() << " != " << baseDiagVals.size() << "+" << extDiagVals.size();
    throw std::runtime_error(errStr.str());
  }
  for (Teuchos::Ordinal i=0; i<baseDiagVals.size(); ++i)
    allDiagVals[i] = baseDiagVals[i];
  Teuchos_Ordinal offset=baseDiagVals.size();
  for (Teuchos::Ordinal i=0; i<extDiagVals.size(); ++i)
    allDiagVals[i+offset] = extDiagVals[i];
}


template<class MatrixType>
void
OverlappingRowMatrix<MatrixType>::
leftScale (const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& /* x */)
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix does not support leftScale.");
}


template<class MatrixType>
void
OverlappingRowMatrix<MatrixType>::
rightScale (const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& /* x */)
{
  throw std::runtime_error("Ifpack2::OverlappingRowMatrix does not support leftScale.");
}


template<class MatrixType>
typename OverlappingRowMatrix<MatrixType>::mag_type
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
  using MV = Tpetra::MultiVector<scalar_type, local_ordinal_type,
                                 global_ordinal_type, node_type>;
  TEUCHOS_TEST_FOR_EXCEPTION
    (X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Ifpack2::OverlappingRowMatrix::apply: X.getNumVectors() = "
     << X.getNumVectors() << " != Y.getNumVectors() = " << Y.getNumVectors()
     << ".");
  // If X aliases Y, we'll need to copy X.
  bool aliases = X.aliases(Y);
  if (aliases) {
    MV X_copy (X, Teuchos::Copy);
    this->apply (X_copy, Y, mode, alpha, beta);
    return;
  }

  const auto& rowMap0 = * (A_->getRowMap ());
  const auto& colMap0 = * (A_->getColMap ());
  MV X_0 (X, mode == Teuchos::NO_TRANS ? colMap0 : rowMap0, 0);
  MV Y_0 (Y, mode == Teuchos::NO_TRANS ? rowMap0 : colMap0, 0);
  A_->localApply (X_0, Y_0, mode, alpha, beta);

  const auto& rowMap1 = * (ExtMatrix_->getRowMap ());
  const auto& colMap1 = * (ExtMatrix_->getColMap ());
  MV X_1 (X, mode == Teuchos::NO_TRANS ? colMap1 : rowMap1, 0);
  MV Y_1 (Y, mode == Teuchos::NO_TRANS ? rowMap1 : colMap1, A_->getLocalNumRows ());
  ExtMatrix_->localApply (X_1, Y_1, mode, alpha, beta);
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
            out << "  Number of entries: " << getLocalNumEntries() << std::endl;
            out << "  Max number of entries per row: " << getLocalMaxNumRowEntries() << std::endl;
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
            for (size_t r = 0; r < getLocalNumRows (); ++r) {
              const size_t nE = getNumEntriesInLocalRow(r);
              typename MatrixType::global_ordinal_type gid = getRowMap()->getGlobalElement(r);
              out << std::setw(width) << myRank
                  << std::setw(width) << gid
                  << std::setw(width) << nE;
              if (vl == VERB_EXTREME) {
                if (isGloballyIndexed()) {
                  global_inds_host_view_type rowinds;
                  values_host_view_type rowvals;
                  getGlobalRowView (gid, rowinds, rowvals);
                  for (size_t j = 0; j < nE; ++j) {
                    out << " (" << rowinds[j]
                        << ", " << rowvals[j]
                        << ") ";
                  }
                }
                else if (isLocallyIndexed()) {
                  local_inds_host_view_type rowinds;
                  values_host_view_type rowvals;
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

template<class MatrixType>
Teuchos::RCP<const Tpetra::CrsMatrix<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> >
OverlappingRowMatrix<MatrixType>::getUnderlyingMatrix() const
{
  return A_;
}

template<class MatrixType>
Teuchos::RCP<const Tpetra::CrsMatrix<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> >
OverlappingRowMatrix<MatrixType>::getExtMatrix() const
{
  return ExtMatrix_;
}

template<class MatrixType>
Kokkos::View<size_t*, typename OverlappingRowMatrix<MatrixType>::device_type>
OverlappingRowMatrix<MatrixType>::getExtHaloStarts() const
{
  return ExtHaloStarts_;
}

template<class MatrixType>
typename Kokkos::View<size_t*, typename OverlappingRowMatrix<MatrixType>::device_type>::HostMirror
OverlappingRowMatrix<MatrixType>::getExtHaloStartsHost() const
{
  return ExtHaloStarts_h;
}

template<class MatrixType>
void OverlappingRowMatrix<MatrixType>::doExtImport()
{
  ExtMatrix_->resumeFill();
  ExtMatrix_->doImport (*A_, *ExtImporter_, Tpetra::REPLACE);
  ExtMatrix_->fillComplete (A_->getDomainMap (), RowMap_);
}

} // namespace Ifpack2

#define IFPACK2_OVERLAPPINGROWMATRIX_INSTANT(S,LO,GO,N)                 \
  template class Ifpack2::OverlappingRowMatrix< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif // IFPACK2_OVERLAPPINGROWMATRIX_DEF_HPP
