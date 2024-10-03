// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_Filtered_GlobalIndexer_impl_hpp__
#define __Panzer_Filtered_GlobalIndexer_impl_hpp__

#include <unordered_set>

#include "PanzerDofMgr_config.hpp"
#include "Panzer_NodeType.hpp"
#include "Panzer_Filtered_GlobalIndexer.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Vector.hpp"

namespace panzer {

Filtered_GlobalIndexer::Filtered_GlobalIndexer(){}

///////////////////////////////////////////////////////////////////////////////
void
Filtered_GlobalIndexer::
initialize(const Teuchos::RCP<const GlobalIndexer> & ugi,
           const std::vector<panzer::GlobalOrdinal>& filtered)
{
  using GO = panzer::GlobalOrdinal;
  using LO = panzer::LocalOrdinal;
  using Node = panzer::TpetraNodeType;
  using Map = Tpetra::Map<LO, GO, Node>;
  using Vector = Tpetra::Vector<GO,LO,GO,Node>;
  using Export = Tpetra::Export<LO,GO,Node>;

  using std::size_t;
  using std::vector;
  using HashTable = std::unordered_set<panzer::GlobalOrdinal>;
  using Teuchos::RCP;

  owned_.clear();
  ghosted_.clear();
  base_ = ugi;

  // From the base global indexer, build the filtered owned indices.
  vector<panzer::GlobalOrdinal> baseOwned, baseGhosted;
  base_->getOwnedIndices(baseOwned);
  base_->getGhostedIndices(baseGhosted);

  RCP<const Map> ownedMap 
      = Tpetra::createNonContigMapWithNode<LO,GO,Node>(baseOwned,getComm());
  RCP<const Map> ghostedMap 
      = Tpetra::createNonContigMapWithNode<LO,GO,Node>(baseGhosted,getComm());

  Vector ownedFiltered(ownedMap);
  Vector ghostedFiltered(ghostedMap);

  ownedFiltered.putScalar(0.0);
  ghostedFiltered.putScalar(0.0);

  {
    auto ownedFiltered_host = ownedFiltered.getLocalViewHost(Tpetra::Access::ReadWrite);
    auto ghostedFiltered_host = ghostedFiltered.getLocalViewHost(Tpetra::Access::ReadWrite);

    for(panzer::GlobalOrdinal f : filtered) {
      bool isOwned = std::find(baseOwned.begin(),baseOwned.end(),f)!=baseOwned.end();
      bool isGhosted = std::find(baseGhosted.begin(),baseGhosted.end(),f)!=baseGhosted.end();

      if(isOwned)
        ownedFiltered_host(ownedMap->getLocalElement(f),0) = 1.0;
      else if(isGhosted)
        ghostedFiltered_host(ghostedMap->getLocalElement(f),0) = 1.0;
      // else no one cares...
    }
  }

  Export exporter(ghostedMap,ownedMap);
  ownedFiltered.doExport(ghostedFiltered, exporter, Tpetra::ADD);

  Teuchos::ArrayRCP<const panzer::GlobalOrdinal> data = ownedFiltered.getData();

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

///////////////////////////////////////////////////////////////////////////////
void 
Filtered_GlobalIndexer::
getOwnedAndGhostedNotFilteredIndicator(std::vector<int> & indicator) const
{
  using Teuchos::RCP;

  using GO = panzer::GlobalOrdinal;
  using LO = panzer::LocalOrdinal;
  using Node = panzer::TpetraNodeType;
  using Map = Tpetra::Map<LO, GO, Node>;
  using Vector = Tpetra::Vector<GO,LO,GO,Node>;
  using Import = Tpetra::Import<LO,GO,Node>;

  std::vector<panzer::GlobalOrdinal> ownedIndices;
  std::vector<panzer::GlobalOrdinal> ghostedIndices;

  // build owned and ghosted maps
  getOwnedIndices(ownedIndices);
  getOwnedAndGhostedIndices(ghostedIndices);

  RCP<const Map> ownedMap 
      = Tpetra::createNonContigMapWithNode<LO,GO,Node>(ownedIndices,getComm());
  RCP<const Map> ghostedMap 
      = Tpetra::createNonContigMapWithNode<LO,GO,Node>(ghostedIndices,getComm());

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

///////////////////////////////////////////////////////////////////////////////
void 
Filtered_GlobalIndexer::
getFilteredOwnedAndGhostedIndices(std::vector<panzer::GlobalOrdinal> & indices) const
{
  using Teuchos::RCP;

  // get filtered/unfiltered indicator vector
  std::vector<int> indicators;
  getOwnedAndGhostedNotFilteredIndicator(indicators);

  // build ghosted maps
  std::vector<panzer::GlobalOrdinal> ghostedIndices;
  getOwnedAndGhostedIndices(ghostedIndices);

  // filtered out filtered indices (isn't that a useful comment)
  for(std::size_t i=0;i<indicators.size();i++) {
    if(indicators[i]==1)
      indices.push_back(ghostedIndices[i]);
  }
}

///////////////////////////////////////////////////////////////////////////////
void 
Filtered_GlobalIndexer::
ownedIndices(const std::vector<panzer::GlobalOrdinal> & indices,std::vector<bool> & isOwned) const
{
  //Resizes the isOwned array.
  if(indices.size()!=isOwned.size())
    isOwned.resize(indices.size(),false);
  typename std::vector<panzer::GlobalOrdinal>::const_iterator endOf = owned_.end();
  for (std::size_t i = 0; i < indices.size(); ++i) {
    isOwned[i] = ( std::find(owned_.begin(), owned_.end(), indices[i])!=endOf );
  }
}

}

#endif
