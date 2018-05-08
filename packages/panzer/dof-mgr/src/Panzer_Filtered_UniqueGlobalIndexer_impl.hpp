// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_Filtered_UniqueGlobalIndexer_impl_hpp__
#define __Panzer_Filtered_UniqueGlobalIndexer_impl_hpp__

#include <unordered_set>

#include "PanzerDofMgr_config.hpp"
#include "Panzer_NodeType.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Vector.hpp"

namespace panzer {

template <typename LocalOrdinalT,typename GlobalOrdinalT>
Filtered_UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT>::
Filtered_UniqueGlobalIndexer()
{ }

///////////////////////////////////////////////////////////////////////////////
//
//  initialize()
//
///////////////////////////////////////////////////////////////////////////////
template<typename LocalOrdinalT, typename GlobalOrdinalT>
void
Filtered_UniqueGlobalIndexer<LocalOrdinalT, GlobalOrdinalT>::
initialize(
  const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT, GlobalOrdinalT>> & ugi,
  const std::vector<GlobalOrdinalT>& filtered)
{
  typedef GlobalOrdinalT GO;
  typedef LocalOrdinalT LO;
  typedef panzer::TpetraNodeType Node;
  typedef Tpetra::Map<LO, GO, Node> Map;
  typedef Tpetra::Vector<GO,LO,GO,Node> Vector;
  typedef Tpetra::Export<LO,GO,Node> Export;

  using std::size_t;
  using std::vector;
  using HashTable = std::unordered_set<GlobalOrdinalT>;
  using Teuchos::RCP;

  owned_.clear();
  ghosted_.clear();
  base_ = ugi;

  // From the base global indexer, build the filtered owned indices.
  vector<GlobalOrdinalT> baseOwned, baseGhosted;
  base_->getOwnedIndices(baseOwned);
  base_->getGhostedIndices(baseGhosted);

  RCP<const Map> ownedMap 
      = Tpetra::createNonContigMap<LO,GO>(baseOwned,getComm());
  RCP<const Map> ghostedMap 
      = Tpetra::createNonContigMap<LO,GO>(baseGhosted,getComm());

  Vector ownedFiltered(ownedMap);
  Vector ghostedFiltered(ghostedMap);

  ownedFiltered.putScalar(0.0);
  ghostedFiltered.putScalar(0.0);

  for(GlobalOrdinalT f : filtered) {
    bool isOwned = std::find(baseOwned.begin(),baseOwned.end(),f)!=baseOwned.end();
    bool isGhosted = std::find(baseGhosted.begin(),baseGhosted.end(),f)!=baseGhosted.end();
 
    if(isOwned)
      ownedFiltered.replaceGlobalValue(f,1.0);
    else if(isGhosted)
      ghostedFiltered.replaceGlobalValue(f,1.0);
    // else no one cares...
  }

  Export exporter(ghostedMap,ownedMap);
  ownedFiltered.doExport(ghostedFiltered, exporter, Tpetra::ADD);

  Teuchos::ArrayRCP<const GlobalOrdinalT> data = ownedFiltered.getData();

  // Build a hash table for fast searching.
  HashTable filteredHash;
  for(int i(0); i < data.size(); ++i) {
    if(data[i]!=0)
      filteredHash.insert(baseOwned[i]);
  }
  // for (size_t i(0); i < filtered.size(); ++i) {
  //   filteredHash.insert(filtered[i]);
  // }

  // Search for indices in the filtered array; add to owned_ if not found, and
  // add to ghosted_ otherwise.
  for (size_t i(0); i < baseOwned.size(); ++i)
  {
    auto itr = filteredHash.find(baseOwned[i]);
    if (itr == filteredHash.end())
      owned_.push_back(baseOwned[i]);
    else
      ghosted_.push_back(baseOwned[i]);
  }
  ghosted_.insert(ghosted_.end(), baseGhosted.begin(), baseGhosted.end());

  // Now that we've change the owned_ and ghosted_ vectors, we need to rebuild
  // the local IDs.
  this->buildLocalIds();
} // end of initialize()

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void 
Filtered_UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT>::
getOwnedAndGhostedNotFilteredIndicator(std::vector<int> & indicator) const
{
  using Teuchos::RCP;

  typedef GlobalOrdinalT GO;
  typedef LocalOrdinalT LO;
  typedef panzer::TpetraNodeType Node;
  typedef Tpetra::Map<LO, GO, Node> Map;
  typedef Tpetra::Vector<GO,LO,GO,Node> Vector;
  typedef Tpetra::Import<LO,GO,Node> Import;

  std::vector<GlobalOrdinalT> ownedIndices;
  std::vector<GlobalOrdinalT> ghostedIndices;

  // build owned and ghosted maps
  getOwnedIndices(ownedIndices);
  getOwnedAndGhostedIndices(ghostedIndices);

  RCP<const Map> ownedMap 
      = Tpetra::createNonContigMap<LO,GO>(ownedIndices,getComm());
  RCP<const Map> ghostedMap 
      = Tpetra::createNonContigMap<LO,GO>(ghostedIndices,getComm());

  // allocate the owned vector, mark those GIDs as unfiltered
  // (they are by definition)
  Vector ownedActive(ownedMap);
  ownedActive.putScalar(1);

  // Initialize all indices to zero
  Vector ghostedActive(ghostedMap);
  ghostedActive.putScalar(0);

  // do communication, marking unfiltered indices as 1 (filtered
  // indices locally are marked as zero)
  Import importer(ownedMap,ghostedMap);
  ghostedActive.doImport(ownedActive,importer,Tpetra::INSERT);

  Teuchos::ArrayRCP<const GO> data = ghostedActive.getData();

  // copy communicated data (clear it out first)
  indicator.clear();
  indicator.insert(indicator.end(),data.begin(),data.end()); 
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void 
Filtered_UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT>::
getFilteredOwnedAndGhostedIndices(std::vector<GlobalOrdinalT> & indices) const
{
  using Teuchos::RCP;

  // get filtered/unfiltered indicator vector
  std::vector<int> indicators;
  getOwnedAndGhostedNotFilteredIndicator(indicators);

  // build ghosted maps
  std::vector<GlobalOrdinalT> ghostedIndices;
  getOwnedAndGhostedIndices(ghostedIndices);

  // filtered out filtered indices (isn't that a useful comment)
  for(std::size_t i=0;i<indicators.size();i++) {
    if(indicators[i]==1)
      indices.push_back(ghostedIndices[i]);
  }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void 
Filtered_UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT>::
ownedIndices(const std::vector<GlobalOrdinalT> & indices,std::vector<bool> & isOwned) const
{
  //Resizes the isOwned array.
  if(indices.size()!=isOwned.size())
    isOwned.resize(indices.size(),false);
  typename std::vector<GlobalOrdinalT>::const_iterator endOf = owned_.end();
  for (std::size_t i = 0; i < indices.size(); ++i) {
    isOwned[i] = ( std::find(owned_.begin(), owned_.end(), indices[i])!=endOf );
  }
}

}

#endif
