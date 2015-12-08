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
#include "Tpetra_Vector.hpp"

namespace panzer {

template <typename LocalOrdinalT,typename GlobalOrdinalT>
Filtered_UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT>::
Filtered_UniqueGlobalIndexer()
{ }

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void
Filtered_UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT>::
initialize(const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > & ugi,
           const std::vector<GlobalOrdinalT> & filtered)
{ 
  typedef std::unordered_set<GlobalOrdinalT> HashTable;

  base_ = ugi;

  // ensure the localIDs match with the users 
  // this is essential for a class to be a decorator
  this->shareLocalIDs(*base_);

  // from base global indexer build the filtered owned indices
  std::vector<GlobalOrdinalT> baseOwned;
  base_->getOwnedIndices(baseOwned);

  // build a hash table for fast searching
  HashTable filteredHash;
  for(std::size_t i=0;i<filtered.size();i++)
    filteredHash.insert(filtered[i]);

  // search for indices in filtered array, add to owned_ if not found
  for(std::size_t i=0;i<baseOwned.size();i++) {
    typename HashTable::const_iterator itr = filteredHash.find(baseOwned[i]);    

    if(itr==filteredHash.end())
      owned_.push_back(baseOwned[i]);
  }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void 
Filtered_UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT>::
getOwnedAndSharedNotFilteredIndicator(std::vector<int> & indicator) const
{
  using Teuchos::RCP;

  typedef GlobalOrdinalT GO;
  typedef LocalOrdinalT LO;
  typedef panzer::TpetraNodeType Node;
  typedef Tpetra::Map<LO, GO, Node> Map;
  typedef Tpetra::Vector<GO,LO,GO,Node> Vector;
  typedef Tpetra::Import<LO,GO,Node> Import;

  std::vector<GlobalOrdinalT> uniqueIndices;
  std::vector<GlobalOrdinalT> ghostedIndices;

  // build unique and ghosted maps
  getOwnedIndices(uniqueIndices);
  getOwnedAndSharedIndices(ghostedIndices);

  RCP<const Map> uniqueMap 
      = Tpetra::createNonContigMap<LO,GO>(uniqueIndices,getComm());
  RCP<const Map> ghostedMap 
      = Tpetra::createNonContigMap<LO,GO>(ghostedIndices,getComm());

  // allocate the unique vector, mark those GIDs as unfiltered
  // (they are by definition)
  Vector uniqueActive(uniqueMap);
  uniqueActive.putScalar(1);

  // Initialize all indices to zero
  Vector ghostedActive(ghostedMap);
  ghostedActive.putScalar(0);

  // do communication, marking unfiltered indices as 1 (filtered
  // indices locally are marked as zero)
  Import importer(uniqueMap,ghostedMap);
  ghostedActive.doImport(uniqueActive,importer,Tpetra::INSERT);

  Teuchos::ArrayRCP<const GO> data = ghostedActive.getData();

  // copy communicated data (clear it out first)
  indicator.clear();
  indicator.insert(indicator.end(),data.begin(),data.end()); 
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void 
Filtered_UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT>::
getFilteredOwnedAndSharedIndices(std::vector<GlobalOrdinalT> & indices) const
{
  using Teuchos::RCP;

  // get filtered/unfiltered indicator vector
  std::vector<int> indicators;
  getOwnedAndSharedNotFilteredIndicator(indicators);

  // build ghosted maps
  std::vector<GlobalOrdinalT> ghostedIndices;
  getOwnedAndSharedIndices(ghostedIndices);

  // filtered out filtered indices (isn't that a useful comment)
  for(std::size_t i=0;i<indicators.size();i++) {
    if(indicators[i]==1)
      indices.push_back(ghostedIndices[i]);
  }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void 
Filtered_UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT>::
getOwnedIndices(std::vector<GlobalOrdinalT> & indices) const
{
  indices.resize(owned_.size());
  for (size_t i = 0; i < owned_.size(); ++i) {
    indices[i]=owned_[i];
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
