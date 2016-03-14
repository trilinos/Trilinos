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

/*! \file Zoltan2_MappingSolution.hpp
    \brief Defines the MappingSolution class.
*/

#ifndef _ZOLTAN2_MAPPINGSOLUTION_HPP_
#define _ZOLTAN2_MAPPINGSOLUTION_HPP_
#include "Teuchos_Comm.hpp"
#include "Zoltan2_Environment.hpp"
#include "Zoltan2_MachineRepresentation.hpp"

namespace Zoltan2 {

/*! \brief PartitionMapping maps a solution or an input distribution to ranks.
*/

template <typename Adapter>
class MappingSolution : public Solution
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::part_t part_t;
#endif

/*! \brief Constructor 
 */
  MappingSolution() : maxPart(0), nRanks(1) {}

  ~MappingSolution() {}

  typedef std::unordered_map<part_t, int> rankmap_t;

  /*! \brief Get the parts belonging to a process.
   *  \param rank a process rank
   *  \param numParts on return, set to the number of parts assigned to rank.
   *  \param parts on return, pointer (view) to the parts assigned to rank
   *
   */
  // TODO:  KDDKDD Decide whether information should be avail for any process
  // TODO:  KDDKDD (requiring more storage or a directory) or only for the 
  // TODO:  KDDKDD local process.
  // TODO:  KDDKDD Could require O(nprocs) storage
  void getPartsForRankView(int rank, part_t &numParts, part_t *&parts)
  {
    if ((partsForRankIdx == Teuchos::null) && (rankForPart == Teuchos::null)) {
      numParts = 0;
      parts = NULL;
      throw std::runtime_error("No mapping solution available.");
    }

    if (rank < 0 || rank >= nRanks) {
      numParts = 0;
      parts = NULL;
      throw std::runtime_error("Invalid rank input to getPartsForRankView");
    }

    if (partsForRankIdx == Teuchos::null) {
      // Need data stored in CRS format; create the arrays.
      partsForRankIdx = arcp(new part_t[nRanks+1], 0, nRanks+1, true);
      for (int i = 0; i <= nRanks; i++) partsForRankIdx[i] = 0;

      size_t nParts = rankForPart->size();
      partsForRank = arcp(new part_t[nParts], 0, nParts, true);

      for (typename rankmap_t::iterator it = rankForPart->begin();
           it != rankForPart->end(); it++) {
        partsForRankIdx[it->second+1]++;       
      }
      for (int i = 1; i <= nRanks; i++)
        partsForRankIdx[i] += partsForRankIdx[i-1];
      for (typename rankmap_t::iterator it = rankForPart->begin();
           it != rankForPart->end(); it++) {
        partsForRank[partsForRankIdx[it->second]] = it->first;
        partsForRankIdx[it->second]++;
      }
      for (int i = nRanks; i > 0; i--)
        partsForRankIdx[i] = partsForRankIdx[i-1];
      partsForRankIdx[0] = 0;
    }

    numParts = partsForRankIdx[rank+1] - partsForRankIdx[rank];
    parts = &(partsForRank[rank]);
  }

  /*! \brief Get the rank containing a part.
   *  Simplifying assumption:  a part is wholy assigned to a rank; it is not
   *  spread across ranks.
   *  \param part Id of the part whose rank is sought
   *  \return rank to which part is assigned
   */
  // TODO:  KDDKDD Decide how to handle reduced storage case, where entire
  // TODO:  map is not stored on each processor.

  int getRankForPart(part_t part) {

    if ((partsForRankIdx == Teuchos::null) && (rankForPart == Teuchos::null)) {
      throw std::runtime_error("No mapping solution available.");
    }

    if (part < 0 || part > maxPart) {
      throw std::runtime_error("Invalid part number input to getRankForPart");
    }

    if (rankForPart == Teuchos::null) {
      // Need data stored in unordered_map; create it
      rankForPart = rcp(new rankmap_t(partsForRankIdx[nRanks]));
      for (int i = 0; i < nRanks; i++)
        for (part_t j = partsForRankIdx[i]; j < partsForRankIdx[i+1]; j++)
          (*rankForPart)[partsForRank[j]] = i;
    }

    typename rankmap_t::iterator it;
    if ((it = rankForPart->find(part)) != rankForPart->end())
      return it->second;
    else 
      throw std::runtime_error("Invalid part number input to getRankForPart");
  }

  ///////////////////////////////////////////////////////////////////////

  void setMap_PartsForRank(ArrayRCP<int> &idx, ArrayRCP<part_t> &parts) {
    nRanks = idx.size() - 1;
    partsForRankIdx = idx;
    partsForRank = parts;
    size_t nparts = parts.size();
    for (size_t i = 0; i < nparts; i++)
      if (parts[i] > maxPart) maxPart = parts[i];
  }

  void setMap_RankForPart(ArrayRCP<part_t> &parts, ArrayRCP<int> &ranks) {
    size_t nparts = parts.size();
    int maxRank = 0;
    for (size_t i = 0; i < nparts; i++) {
      (*rankForPart)[parts[i]] = ranks[i];
      if (parts[i] > maxPart) maxPart = parts[i];
      if (ranks[i] > maxRank) maxRank = ranks[i];
    }
    nRanks = maxRank+1;
  }

  void setMap_RankForPart(RCP<rankmap_t> &rankmap) {
    rankForPart = rankmap;
    int maxRank = 0;
    for (typename rankmap_t::iterator it = rankForPart->begin();
         it != rankForPart->end(); it++) {
      if (it->first > maxPart) maxPart = it->first;
      if (it->second > maxRank) maxRank = it->second;
    }
    nRanks = maxRank+1;
  }
  // TODO:  can add other methods for setting the map, particularly if decide
  // TODO:  to store only local procs and parts info rather than global info.

private:

  part_t maxPart;
  int nRanks;
  ArrayRCP<part_t> partsForRankIdx;
  ArrayRCP<part_t> partsForRank;
  RCP<rankmap_t> rankForPart;


};

}  // namespace Zoltan2

#endif
