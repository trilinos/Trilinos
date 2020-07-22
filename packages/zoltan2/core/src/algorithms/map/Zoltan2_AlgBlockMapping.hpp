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
#ifndef _ZOLTAN2_ALGBLOCKMAPPING_HPP_
#define _ZOLTAN2_ALGBLOCKMAPPING_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Util.hpp>
#include <Tpetra_Map.hpp>
#include <set>
#include <cmath>

//////////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgBlockMapping.hpp
//! \brief Define a simple mapping of parts to processors assuming parts
//         are contiguously numbered from 0 to nParts-1
//         This use case is common; because it requires no explicit storage
//         of the parts to ranks,
//         we separate it from methods that do need extra storage.
//////////////////////////////////////////////////////////////////////////////

namespace Zoltan2 {

template <typename Adapter, typename MachineRep>
class AlgBlockMapping : public Algorithm<Adapter>
{

private:

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::userCoord_t userCoord_t;

public:

  /*! Constructor that can be accessed directly by user through MappingProblem.
   *  \param env  parameters for the problem and library configuration
   *  \param comm  the communicator for the problem
   *  \param adapter  the user's input adapter
   */
  AlgBlockMapping(
    const Teuchos::RCP <const Teuchos::Comm<int> > &comm_,
    const Teuchos::RCP <const MachineRep> &machine_,
    const Teuchos::RCP <const Adapter> &adapter_,
    const Teuchos::RCP <const Zoltan2::PartitioningSolution<Adapter> > &psoln_,
    const Teuchos::RCP <const Environment> &envConst):
    nRanks(comm_->getSize()), myRank(comm_->getRank()),
    nMyParts(0), myParts(Teuchos::null)
  { 
    // Get set of parts from partitioning solution or adapter, if provided
    // Otherwise, we'll assume set of parts == set of ranks
    const part_t *partList;
    if (psoln_ != Teuchos::null) {
      partList = psoln_->getPartListView();
    }
    else {
      adapter_->getPartsView(partList);
    }

    // Compute maxPart, nParts

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

      size_t nLocal = adapter_->getLocalNumIDs();

      // Ideally, we'd use part_t as global ID in the map, but that won't
      // work if Tpetra is compiled without global ID = int.
      //    typedef part_t use_this_gno_t;
      // We'll use Tpetra's default instead
      typedef Tpetra::Map<>::global_ordinal_type use_this_gno_t;

      std::set<use_this_gno_t> unique;
      for (size_t i = 0; i < nLocal; i++)
        unique.insert(partList[i]);

      size_t nUnique = unique.size();
      Array<use_this_gno_t> uniquePartList(nUnique);
      size_t k = 0;
      for (typename std::set<use_this_gno_t>::iterator it = unique.begin();
           it != unique.end(); it++)
        uniquePartList[k++] = *it;
        
      // Use Tpetra::Map to find the max, min, total part.

      global_size_t nGlobalElts = 
                    Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
      Tpetra::Map<lno_t, use_this_gno_t> tmap(nGlobalElts, uniquePartList(),
                                              0, comm_);

      nParts = Teuchos::as<part_t>(tmap.getGlobalNumElements());
      minPart = tmap.getMinAllGlobalIndex();
      maxPart = tmap.getMaxAllGlobalIndex();
    }

    // Determine number of unique parts, as well as min and max part
    // Can probably use a Tpetra Map.

    if ((minPart != 0) || (maxPart != nParts-1)) {
      // Cannot use this mapping method
      throw std::runtime_error("Cannot use mapping_algorithm = contiguous "
                               "unless parts are numbered from 0 to nParts-1");
    }

    sharedConstructor();
  }

  //! \brief Constructor that allows this mapping method to be used as a 
  //         submethod of the default
  //         Assume that calling method has already confirmed assumptions
  //         about part numbering.
  AlgBlockMapping(
    const Teuchos::RCP <const Teuchos::Comm<int> > comm_,
    const part_t nparts) :
    nRanks(comm_->getSize()),
    nParts(nparts),  
    myRank(comm_->getRank()),
    nMyParts(0),
    myParts(Teuchos::null)
  {
    sharedConstructor();
  }

  void sharedConstructor()
  {
    nParts_Div_nRanks = nParts / nRanks;
    nParts_Mod_nRanks = nParts % nRanks;
    nMyParts = nParts_Div_nRanks + (myRank < nParts_Mod_nRanks);
  }

  void map(const Teuchos::RCP<MappingSolution<Adapter> > &msoln) {
    // Mapping is computed implicitly; 
    // nothing to do here since we implement 
    // getRankForPart and getMyPartsView.
  }

  int getRankForPart(part_t p) {
    if (p < 0 || p >= nParts)
      throw std::runtime_error("Invalid part in getRankForPart");

    int tmp = p / nParts_Div_nRanks;
    while (firstPart(tmp) > p) tmp--;
    while (firstPart(tmp+1) < p) tmp++;
    return tmp;
  }

  void getMyPartsView(part_t &numParts, part_t *&parts)
  {
    // Will need to construct the arrays if this function is called.
    // Will not work if this is a const method.
    numParts = nMyParts;
    if (nMyParts) {
      if (myParts == Teuchos::null) {
        myParts = arcp(new part_t[nMyParts], 0, nMyParts, true);
        for (part_t i = 0; i < nMyParts; i++)
          myParts[i] = firstPart(myRank) + i;
      }
      parts = myParts.getRawPtr();
    }
    else
      parts = NULL;
  }

private:

  inline part_t firstPart(int rank) {
    // For contiguous part assignments, the first part assigned to rank
    return (rank * nParts_Div_nRanks + std::min(rank, nParts_Mod_nRanks));
  }

  const int nRanks;            // Global number of ranks
  part_t nParts;               // Global number of parts
  part_t nParts_Div_nRanks;    // (nParts/nRanks) precomputed for frequent use
  part_t nParts_Mod_nRanks;    // (nParts%nRanks) precomputed for frequent use

  const int myRank;            // Local rank
  part_t nMyParts;       // Local number of parts
  ArrayRCP<part_t> myParts;  // Array of my parts; created only if 
                                   // getMyPartsView is called.
}; 


} // namespace Zoltan2

#endif
