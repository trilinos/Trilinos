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
#ifndef _ZOLTAN2_ALGDEFAULTMAPPING_HPP_
#define _ZOLTAN2_ALGDEFAULTMAPPING_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Util.hpp>


//////////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgDefaultMapping.hpp
//! \brief Define a default mapping of parts to processors
//  
//////////////////////////////////////////////////////////////////////////////

namespace Zoltan2 {

template <typename Adapter, typename MachineRep>
class AlgDefaultMapping : public Algorithm<Adapter>
{

private:

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::userCoord_t userCoord_t;

  const RCP<const Environment> env;
  const RCP<const Comm<int> > problemComm;
  const RCP<const typename Adapter::base_adapter_t> adapter;

public:

  /*! Constructor
   *  \param env  parameters for the problem and library configuration
   *  \param problemComm  the communicator for the problem
   *  \param adapter  the user's input adapter
   */
  AlgDefaultMapping(
    const Teuchos::RCP <const Teuchos::Comm<int> > &comm_,
    const Teuchos::RCP <const MachineRep> &machine_,
    const Teuchos::RCP <const Adapter> &adapter_,
    const Teuchos::RCP <const Zoltan2::PartitioningSolution<Adapter> > &psoln_,
    const Teuchos::RCP <const Environment> &envConst)
  : env(envConst), adapter(adapter_),
    partIsRank(false), haveContiguousParts(false)
    rankForPart(Teuchos::null)
  { 
    int nRanks = comm->getSize();

    // Get set of parts from partitioning solution or adapter, if provided
    // Otherwise, we'll assume set of parts == set of ranks
    const *partList;
    if (psoln_ != Teuchos::null) {
      partList = psoln_->getPartListView();
    }
    else {
      adapter->getPartsView(partList);
    }

    // Compute maxPart, nParts
    typedef typename Tpetra::Map<lno_t, part_t> tpetraMap_t
    Teuchos::RCP<tpetraMap_t> tmap;

    part_t minPart, maxPart;

    if (partList == NULL) {
      // Parts for IDs are the same as ranks
      nParts = nRanks;
      maxPart = nRanks-1;
      minPart = 0;
    }
    else {
      // Parts were provided by partitioning solution or input adapter

      // Find unique parts on this rank, 
      // as Tpetra::Map does not like duplicate GIDs within a rank

      size_t nLocal = adapter->getLocalNumIDs();

      std::set<part_t> unique(nLocal);
      for (size_t i; i < adapter->getLocalNumIDs(); i++)
        unique.insert(partList[i]);

      size_t nUnique = unique.size();
      Array<const part_t> uniquePartList(nUnique);
      size_t k = 0;
      for (typename std::set<part_t>::iterator it = set.begin();
           it != set.end(); it++)
        uniquePartList[k++] = *it;
        
      // Use Tpetra::Map to find the max, min, total part.

      global_size_t nGlobalElts = 
                    Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
      tmap = rcp(new tpetraMap_t(nGlobalElts, uniquePartList(), 0, comm));

      nParts = as<part_t>(tmap->getGlobalNumElements());
      minPart = tmap->getMinAllGlobalIndex();
      maxPart = tmap->getMaxAllGlobalIndex();
    }

    nParts_Div_nRanks = nParts / nRanks;
    nParts_Mod_nRanks = nParts % nRanks;

    // Determine number of unique parts, as well as min and max part
    // Can probably use a Tpetra Map.

    if (maxPart < nRanks)
      partIsRank = true;  // Most common case
      
    if ((minPart == 0) && (maxPart == nParts-1))
      // Have a contiguous set of parts; can use simplest default mapping
      haveContiguousParts = true;

    if (!partIsRank && !haveContiguousParts) {
      // Use Tpetra createOneToOne to map parts to ranks
      Teuchos::RCP<tpetraMap_t> oneToOneMap = Tpetra::createOneToOne(tmap);

      // Each rank needs names of all parts
      // Should we gather map to one rank, or copy to all?
      // I don't know how to do it.
      Teuchos::RCP<tpetraMap_t> gatheredMap = ;
     
      Teuchos::ArrayView<const part_t> allParts = 
                                       gatheredMap->getNodeElementList();

      // Look up the rank for all parts assigned by oneToOneMap
      Teuchos::Array<int> allRanks(allParts.size());
      oneToOneMap->getRemoveIndexList(allParts, allRanks());
      
      for (size_t i = 0; i < allPart.size())
        (*rankForPart)[allParts[i]] = allRanks[i];
    }
  }

  void map(Teuchos::RCP<MappingSolution<Adapter> > msoln) {
    // Mapping was computed in constructor; 
    // nothing to do here since we implement getRankForPart and
    // getMyPartsView.  (If we didn't we would need to store
    // the mapping in the soln here.)
  }

  int getRankForPart(part_t p) {
    if (p < 0 || p > maxPart)
      throw std::runtime_error("Invalid part in getRankForPart");

    // Most common case:  at most one part per rank
    if (partIsRank) 
      return as<int>(p);

    // Multiple parts per rank
    if (haveContiguousParts) { 
      int tmp = p / nParts_Div_nRanks;
      while (firstContiguousPart(tmp) > p) tmp--;
      while (firstContiguousPart(tmp+1) < p) tmp++;
      return tmp;
    }

    // None of the simple schemes apply; use the unordered_map.
    return rankForPart[p];
  }

  void getMyPartsView(part_t &numParts, part_t *&parts)
  {
    // Will need to construct the arrays if this function is called.
    // Will not work if this is a const method.
    partsForRankIdx = rcp(new part_t[nRanks+1]);
    partsForRank = rcp(new part_t[nParts]);
      
  }

  void map(const RCP<PartitioningSolution<Adapter> > &solution);

private:

  inline part_t firstContiguousPart(int rank) {
    // For contiguous part assignments, the first part assigned to rank
    return (rank * nParts_Div_nRanks + min(rank, nParts_Mod_nRanks));
  }
  bool partIsRank;
  bool haveContiguousParts;
  Teuchos::RCP<rankForPart_t> rankForPart;
};


} // namespace Zoltan2

#endif
