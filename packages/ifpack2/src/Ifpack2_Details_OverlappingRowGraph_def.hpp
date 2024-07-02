// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_OVERLAPPINGROWGRAPH_DEF_HPP
#define IFPACK2_OVERLAPPINGROWGRAPH_DEF_HPP

#include <Ifpack2_Details_OverlappingRowGraph_decl.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Export.hpp>


namespace Ifpack2 {
namespace Details {

template<class GraphType>
OverlappingRowGraph<GraphType>::
OverlappingRowGraph (const Teuchos::RCP<const row_graph_type>& nonoverlappingGraph,
                     const Teuchos::RCP<const row_graph_type>& overlappingGraph,
                     const Teuchos::RCP<const map_type>& rowMap,
                     const Teuchos::RCP<const map_type>& colMap,
                     const Tpetra::global_size_t numGlobalRows,
                     const Tpetra::global_size_t numGlobalCols,
                     const Tpetra::global_size_t numGlobalNonzeros,
                     const size_t maxNumEntries,
                     const Teuchos::RCP<const import_type>& nonoverlappingImporter,
                     const Teuchos::RCP<const import_type>& overlappingImporter) :
  nonoverlappingGraph_ (nonoverlappingGraph),
  overlappingGraph_ (overlappingGraph),
  rowMap_ (rowMap),
  colMap_ (colMap),
  numGlobalRows_ (numGlobalRows),
  numGlobalCols_ (numGlobalCols),
  numGlobalNonzeros_ (numGlobalNonzeros),
  maxNumEntries_ (maxNumEntries),
  nonoverlappingImporter_ (nonoverlappingImporter),
  overlappingImporter_ (overlappingImporter)
{}


template<class GraphType>
OverlappingRowGraph<GraphType>::~OverlappingRowGraph() {}


template<class GraphType>
Teuchos::RCP<const Teuchos::Comm<int> > 
OverlappingRowGraph<GraphType>::getComm () const
{
  return nonoverlappingGraph_->getComm ();
}
  

  

template<class GraphType>
Teuchos::RCP<const Tpetra::Map<typename GraphType::local_ordinal_type, typename GraphType::global_ordinal_type, typename GraphType::node_type> > 
OverlappingRowGraph<GraphType>::getRowMap () const
{
  return rowMap_;
}
  

template<class GraphType>
Teuchos::RCP<const Tpetra::Map<typename GraphType::local_ordinal_type, typename GraphType::global_ordinal_type, typename GraphType::node_type> > 
OverlappingRowGraph<GraphType>::getColMap () const
{
  return colMap_;
}


template<class GraphType>
Teuchos::RCP<const Tpetra::Map<typename GraphType::local_ordinal_type, typename GraphType::global_ordinal_type, typename GraphType::node_type> > 
OverlappingRowGraph<GraphType>::getDomainMap () const
{
  return nonoverlappingGraph_->getDomainMap ();
}


template<class GraphType>
Teuchos::RCP<const Tpetra::Map<typename GraphType::local_ordinal_type, typename GraphType::global_ordinal_type, typename GraphType::node_type> >
OverlappingRowGraph<GraphType>::getRangeMap () const
{
  return nonoverlappingGraph_->getRangeMap ();
}


template<class GraphType>
Teuchos::RCP<const Tpetra::Import<typename GraphType::local_ordinal_type, typename GraphType::global_ordinal_type, typename GraphType::node_type> >
OverlappingRowGraph<GraphType>::getImporter () const
{
  return nonoverlappingImporter_;
}


template<class GraphType>
Teuchos::RCP<const Tpetra::Export<typename GraphType::local_ordinal_type, typename GraphType::global_ordinal_type, typename GraphType::node_type> >
OverlappingRowGraph<GraphType>::getExporter () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}
  

template<class GraphType>
global_size_t OverlappingRowGraph<GraphType>::getGlobalNumRows () const
{
  return numGlobalRows_;
}
  

template<class GraphType>
global_size_t OverlappingRowGraph<GraphType>::getGlobalNumCols () const
{
  return numGlobalCols_;
}
  

template<class GraphType>
size_t OverlappingRowGraph<GraphType>::getLocalNumRows () const
{
  return nonoverlappingGraph_->getLocalNumRows () + 
    overlappingGraph_->getLocalNumRows ();
}
  

template<class GraphType>
size_t OverlappingRowGraph<GraphType>::getLocalNumCols () const
{
  return this->getLocalNumRows ();
}
  

template<class GraphType>
typename GraphType::global_ordinal_type 
OverlappingRowGraph<GraphType>::getIndexBase () const
{
  return nonoverlappingGraph_->getIndexBase ();
}
  

template<class GraphType>
Tpetra::global_size_t OverlappingRowGraph<GraphType>::getGlobalNumEntries () const
{
  return numGlobalNonzeros_;
}
  

template<class GraphType>
size_t OverlappingRowGraph<GraphType>::getLocalNumEntries () const
{
  return nonoverlappingGraph_->getLocalNumEntries () + 
    overlappingGraph_->getLocalNumEntries ();
}
  

template<class GraphType>
size_t
OverlappingRowGraph<GraphType>::
getNumEntriesInGlobalRow (global_ordinal_type globalRow) const
{
  const local_ordinal_type localRow = rowMap_->getLocalElement (globalRow);
  if (localRow == Teuchos::OrdinalTraits<local_ordinal_type>::invalid ()) {
    return Teuchos::OrdinalTraits<size_t>::invalid();
  } else {
    return getNumEntriesInLocalRow (localRow);
  }
}

  
template<class GraphType>
size_t
OverlappingRowGraph<GraphType>::
getNumEntriesInLocalRow (local_ordinal_type localRow) const
{
  using Teuchos::as;
  const size_t numMyRowsA = nonoverlappingGraph_->getLocalNumRows ();
  if (as<size_t> (localRow) < numMyRowsA) {
    return nonoverlappingGraph_->getNumEntriesInLocalRow (localRow);
  } else {
    return overlappingGraph_->getNumEntriesInLocalRow (as<local_ordinal_type> (localRow - numMyRowsA));
  }
}
  

template<class GraphType>
size_t OverlappingRowGraph<GraphType>::getGlobalMaxNumRowEntries () const
{
  throw std::runtime_error("Ifpack2::OverlappingRowGraph::getGlobalMaxNumRowEntries() not supported.");
}
  

template<class GraphType>
size_t OverlappingRowGraph<GraphType>::getLocalMaxNumRowEntries () const
{
  return maxNumEntries_;
}
  

template<class GraphType>
bool OverlappingRowGraph<GraphType>::hasColMap () const
{
  return true;
}
  

template<class GraphType>
bool OverlappingRowGraph<GraphType>::isLocallyIndexed () const
{
  return true;
}
   

template<class GraphType>
bool OverlappingRowGraph<GraphType>::isGloballyIndexed () const
{
  return false;
}
  

template<class GraphType>
bool OverlappingRowGraph<GraphType>::isFillComplete () const
{
  return true;
}
 
template<class GraphType>
void
OverlappingRowGraph<GraphType>::
  getGlobalRowCopy (global_ordinal_type globalRow,
                    nonconst_global_inds_host_view_type& indices,
                    size_t& numIndices) const
{
  const local_ordinal_type localRow = rowMap_->getLocalElement (globalRow);
  if (localRow == Teuchos::OrdinalTraits<local_ordinal_type>::invalid ()) {
    numIndices = Teuchos::OrdinalTraits<size_t>::invalid ();
  } else {
    if (Teuchos::as<size_t> (localRow) < nonoverlappingGraph_->getLocalNumRows ()) {
      nonoverlappingGraph_->getGlobalRowCopy (globalRow, indices, numIndices);
    } else {
      overlappingGraph_->getGlobalRowCopy (globalRow, indices, numIndices);
    }
  }
}


template<class GraphType>
void
OverlappingRowGraph<GraphType>::
getLocalRowCopy (local_ordinal_type localRow,
                 nonconst_local_inds_host_view_type& indices,
                 size_t& numIndices) const
{
  using Teuchos::as;
  const size_t numMyRowsA = nonoverlappingGraph_->getLocalNumRows ();
  if (as<size_t> (localRow) < numMyRowsA) {
    nonoverlappingGraph_->getLocalRowCopy (localRow, indices, numIndices);
  } else {
    const local_ordinal_type localRowOffset = 
      localRow - as<local_ordinal_type> (numMyRowsA);
    overlappingGraph_->getLocalRowCopy (localRowOffset, indices, numIndices);
  }
}


template<class GraphType>
void
OverlappingRowGraph<GraphType>::
getGlobalRowView (global_ordinal_type GlobalRow,
                  global_inds_host_view_type &indices) const {
  const local_ordinal_type LocalRow = rowMap_->getLocalElement (GlobalRow);
  if (LocalRow == Teuchos::OrdinalTraits<local_ordinal_type>::invalid())  {
    indices = global_inds_host_view_type();
  } else {
    if (Teuchos::as<size_t> (LocalRow) < nonoverlappingGraph_->getLocalNumRows ()) {
      nonoverlappingGraph_->getGlobalRowView (GlobalRow, indices);
    } else {
      overlappingGraph_->getGlobalRowView (GlobalRow, indices);
    }
  }
}


template<class GraphType>
void
OverlappingRowGraph<GraphType>::
  getLocalRowView (local_ordinal_type LocalRow,
                   local_inds_host_view_type & indices) const {
  using Teuchos::as;
  const size_t numMyRowsA = nonoverlappingGraph_->getLocalNumRows ();
  if (as<size_t> (LocalRow) < numMyRowsA) {
    nonoverlappingGraph_->getLocalRowView (LocalRow, indices);
  } else {
    overlappingGraph_->getLocalRowView (LocalRow - as<local_ordinal_type> (numMyRowsA),
                                 indices);
  }

}


} // namespace Details
} // namespace Ifpack2

#define IFPACK2_DETAILS_OVERLAPPINGROWGRAPH_INSTANT(LO,GO,N) \
  template class Ifpack2::Details::OverlappingRowGraph<Tpetra::CrsGraph< LO, GO, N > >; \
  template class Ifpack2::Details::OverlappingRowGraph<Tpetra::RowGraph< LO, GO, N > >;

#endif // IFPACK2_OVERLAPPINGROWGRAPH_DEF_HPP
