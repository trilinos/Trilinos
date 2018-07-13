/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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
Teuchos::RCP<typename GraphType::node_type> 
OverlappingRowGraph<GraphType>::getNode () const
{
  return nonoverlappingGraph_->getNode();
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
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
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
size_t OverlappingRowGraph<GraphType>::getNodeNumRows () const
{
  return nonoverlappingGraph_->getNodeNumRows () + 
    overlappingGraph_->getNodeNumRows ();
}
  

template<class GraphType>
size_t OverlappingRowGraph<GraphType>::getNodeNumCols () const
{
  return this->getNodeNumRows ();
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
size_t OverlappingRowGraph<GraphType>::getNodeNumEntries () const
{
  return nonoverlappingGraph_->getNodeNumEntries () + 
    overlappingGraph_->getNodeNumEntries ();
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
  const size_t numMyRowsA = nonoverlappingGraph_->getNodeNumRows ();
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
size_t OverlappingRowGraph<GraphType>::getNodeMaxNumRowEntries () const
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
                  const Teuchos::ArrayView<global_ordinal_type>& indices,
                  size_t& numIndices) const
{
  const local_ordinal_type localRow = rowMap_->getLocalElement (globalRow);
  if (localRow == Teuchos::OrdinalTraits<local_ordinal_type>::invalid ()) {
    numIndices = Teuchos::OrdinalTraits<size_t>::invalid ();
  } else {
    if (Teuchos::as<size_t> (localRow) < nonoverlappingGraph_->getNodeNumRows ()) {
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
                 const Teuchos::ArrayView<local_ordinal_type>& indices, 
                 size_t& numIndices) const
{
  using Teuchos::as;
  const size_t numMyRowsA = nonoverlappingGraph_->getNodeNumRows ();
  if (as<size_t> (localRow) < numMyRowsA) {
    nonoverlappingGraph_->getLocalRowCopy (localRow, indices, numIndices);
  } else {
    const local_ordinal_type localRowOffset = 
      localRow - as<local_ordinal_type> (numMyRowsA);
    overlappingGraph_->getLocalRowCopy (localRowOffset, indices, numIndices);
  }
}
  
} // namespace Details
} // namespace Ifpack2

#define IFPACK2_DETAILS_OVERLAPPINGROWGRAPH_INSTANT(LO,GO,N) \
  template class Ifpack2::Details::OverlappingRowGraph<Tpetra::CrsGraph< LO, GO, N > >; \
  template class Ifpack2::Details::OverlappingRowGraph<Tpetra::RowGraph< LO, GO, N > >;

#endif // IFPACK2_OVERLAPPINGROWGRAPH_DEF_HPP
