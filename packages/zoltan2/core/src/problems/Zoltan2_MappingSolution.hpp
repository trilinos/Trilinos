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
#include "Zoltan2_PartitioningSolution.hpp"
#include <unordered_map>

namespace Zoltan2 {

/*! \brief PartitionMapping maps a solution or an input distribution to ranks.
*/

template <typename Adapter>
class MappingSolution : public PartitioningSolution<Adapter>
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::scalar_t t_scalar_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
#endif

/*! \brief Constructor 
 */
  MappingSolution(
      const RCP<const Environment> &env,
      const RCP<const Comm<int> > &comm,
      const RCP<Algorithm<Adapter> > &algorithm = Teuchos::null)
    :PartitioningSolution <Adapter> (
      env, comm, 1, algorithm),
      nParts(0), nRanks(1), myRank(comm->getRank()), maxPart(0),
      mapping_algorithm(algorithm) {}

  virtual ~MappingSolution() {}

  typedef std::unordered_map<lno_t, int> rankmap_t;

  /*! \brief Get the parts belonging to this rank
   *  \param numParts on return, set to the number of parts assigned to rank.
   *  \param parts on return, pointer (view) to the parts assigned to rank
   */
  void getMyPartsView(part_t &numParts, part_t *&parts)
  {
    bool useAlg = true;

    // Check first whether this algorithm answers getMyPartsView.
    if (mapping_algorithm != Teuchos::null) {
      try {
        mapping_algorithm->getMyPartsView(numParts, parts);
      }
      catch (NotImplemented &e) {
        // It is OK if the algorithm did not implement this method;
        // we'll get the information from the solution below.
        useAlg = false;
      }
      Z2_FORWARD_EXCEPTIONS;
    }

    if (!useAlg) {  

      // Algorithm did not implement this method.

      // Did the algorithm register a result with the solution?
      if ((partsForRank==Teuchos::null) && (rankForPart==Teuchos::null)) {
        numParts = 0;
        parts = NULL;
        throw std::runtime_error("No mapping solution available.");
      }
  
      if (partsForRank == Teuchos::null) {
        // Need to create the array since we haven't created it before.
        Teuchos::Array<part_t> tmp;

        part_t cnt = 0;
        for (typename rankmap_t::iterator it = rankForPart->begin();
             it != rankForPart->end(); it++) {
          if (it->second == myRank) {
            tmp.push_back(it->first); 
            cnt++;
          }
        }
        if (cnt)
          partsForRank = arcp(&tmp[0], 0, cnt, true);
      }

      numParts = partsForRank.size();
      if (numParts)
        parts = partsForRank.getRawPtr();
      else 
        parts = NULL;
    }
  }

  /*! \brief Get the rank containing a part.
   *  Simplifying assumption:  a part is wholy assigned to a rank; it is not
   *  spread across ranks.
   *  \param part Id of the part whose rank is sought
   *  \return rank to which part is assigned
   */
  int getRankForPart(part_t part) {

    int r = -1;

    // Check first whether this algorithm answers getRankForPart.
    // Some algorithms can compute getRankForPart implicitly, without having
    // to store the mapping explicitly.  It is more efficient for them
    // to implement getRankForPart themselves.
    if (mapping_algorithm != Teuchos::null) {
      try {
        r = mapping_algorithm->getRankForPart(part);
      }
      catch (NotImplemented &e) {
        // It is OK if the algorithm did not implement this method;
        // we'll get the information from the solution below.
      }
      Z2_FORWARD_EXCEPTIONS;
    }

    if (r == -1) {  // Algorithm did not implement this method
      if (rankForPart==Teuchos::null) {
        throw std::runtime_error("No mapping solution available.");
      }

      if (part < 0 || part > maxPart) {
        throw std::runtime_error("Invalid part number input to getRankForPart");
      }


      typename rankmap_t::iterator it;
      if ((it = rankForPart->find(part)) != rankForPart->end())
        r = it->second;
      else 
        throw std::runtime_error("Invalid part number input to getRankForPart");
    }
    return r;
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // Methods for storing mapping data in the solution.
  // Algorithms can store their data in the solution, or implement
  // getRankForPart and getMyPartsView themselves.

  void setMap_PartsForRank(ArrayRCP<int> &idx, ArrayRCP<part_t> &parts) {
    nRanks = idx.size() - 1;
    nParts = parts.size();

    // Need data stored in unordered_map; create it
    rankForPart = rcp(new rankmap_t(idx[nRanks]));

    maxPart = 0;
    for (int i = 0; i < nRanks; i++) {
      for (part_t j = idx[i]; j < idx[i+1]; j++) {
        (*rankForPart)[parts[j]] = i;
        if (parts[j] > maxPart) maxPart = parts[j];
      }
    }

    // Parts for this rank are already contiguous in parts arcp.  
    // Keep a view of them.
    partsForRank = parts.persistingView(idx[myRank],idx[myRank+1]-idx[myRank]);
  }

  /**
   * This is a mapping from local elements to ranks.
   * "parts" in the other functions should also mean the local elements,
   * since we allow the direct mapping call with local elements as well.
   * ranks[i] hold the mappend rank for local element i.
   * Function will fail if part_t != int
   */
  void setMap_RankForLocalElements(ArrayRCP<int> &ranks) {
    this->setParts(ranks);
  }


  void setMap_RankForPart(ArrayRCP<part_t> &parts, ArrayRCP<int> &ranks) {
    nParts = parts.size();
    int maxRank = 0;

    // Need data stored in unordered_map; create it
    rankForPart = rcp(new rankmap_t(parts.size()));

    for (size_t i = 0; i < nParts; i++) {
      (*rankForPart)[parts[i]] = ranks[i];
      if (parts[i] > maxPart) maxPart = parts[i];
      if (ranks[i] > maxRank) maxRank = ranks[i];
    }
    nRanks = maxRank+1;
  }

  void setMap_RankForPart(RCP<rankmap_t> &rankmap) {
    rankForPart = rankmap;
    nParts = rankForPart.size();
    int maxRank = 0;
    typename rankmap_t::iterator it;
    for (it = rankForPart->begin(); it != rankForPart->end(); it++) {
      if (it->first > maxPart) maxPart = it->first;
      if (it->second > maxRank) maxRank = it->second;
    }
    nRanks = maxRank+1;
  }
  // TODO:  can add other methods for setting the map, particularly if decide
  // TODO:  to store only local procs and parts info rather than global info.

private:

  part_t nParts;  // Global number of parts
  int nRanks;     // Global number of ranks
  int myRank;     // This ranks
  part_t maxPart; // Maximum part number

  // Ways to access the answer:  it can be stored in MappingSolution or
  // provided by the Algorithm.
  ArrayRCP<part_t> partsForRankIdx;
  ArrayRCP<part_t> partsForRank;
  RCP<rankmap_t> rankForPart;

  const RCP<Algorithm<Adapter> > mapping_algorithm;

};

}  // namespace Zoltan2

#endif
