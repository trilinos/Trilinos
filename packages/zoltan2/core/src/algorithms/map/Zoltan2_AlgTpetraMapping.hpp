// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

template <typename Adapter>
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
    const Teuchos::RCP <const MachineRepresentation<pcoord_t,part_t> > &machine_,
    const Teuchos::RCP <const Adapter> &adapter_,
    const Teuchos::RCP <const Zoltan2::PartitioningSolution<Adapter> > &psoln_,
    const Teuchos::RCP <const Environment> &envConst):
    env(env__), comm(comm_), adapter(adapter_),
    partIsRank(false), haveContigousParts(false)
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
                                       gatheredMap->getLocalElementList();

      // Look up the rank for all parts assigned by oneToOneMap
      Teuchos::Array<int> allRanks(allParts.size());
      oneToOneMap->getRemoveIndexList(allParts, allRanks());
      
      for (size_t i = 0; i < allPart.size())
        (*rankForPart)[allParts[i]] = allRanks[i];
    }
  }

  void map(Teuchos::RCP<MappingSolution<Adapter> > msoln) {
    // Mapping was computed in constructor; 
    // nothing to do here sine we implement getRankForPart and
    // getPartsForRankView.  (If we didn't we would need to store
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

  void getPartsForRankView(int rank, part_t &numParts, part_t *&parts)
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
};


} // namespace Zoltan2

#endif
