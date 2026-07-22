// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef KDTREE_STK_INTERFACE_H_
#define KDTREE_STK_INTERFACE_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stk_search/kdtree/KDTree_BoundingBox.hpp"
#include "stk_search/kdtree/KDTree_ParallelConsistencyUtils.hpp"
#include "stk_util/environment/Env.hpp"
#include "stk_util/environment/WallTime.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"
#include <stk_search/CommonSearchUtil.hpp>
#include <stk_search/Sphere.hpp>
#include <stk_search/IdentProc.hpp>

namespace stk::search {

//
//  More general search for an arbitrary range type
//
template <typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeObjType>
inline void coarse_search_kdtree(std::vector< std::pair<DomainObjType, DomainIdentifier> > const & local_domain,
                                 std::vector< std::pair<RangeObjType,  RangeIdentifier > > const & local_range,
                                 MPI_Comm comm,
                                 std::vector<std::pair<DomainIdentifier, RangeIdentifier> >& searchResults,
                                 bool enforceSearchResultSymmetry = true,
                                 bool sortSearchResults = false)
{

  int num_procs = -1;
  int proc_id   = -1;
  MPI_Comm_rank(comm, &proc_id);
  MPI_Comm_size(comm, &num_procs);

  searchResults.clear();
  std::vector<RangeObjType> rangeObjs( local_range.size() );

  std::vector<RangeIdentifier> rangeGhostIdentifiers;

  stk::search::ComputeRangeWithGhostsForCoarseSearch(local_domain, local_range,
                                                     num_procs, rangeObjs, rangeGhostIdentifiers, comm);

#ifdef _OPENMP
  std::vector<std::vector<std::pair<DomainIdentifier, RangeIdentifier> > >
      threadLocalSearchResults( omp_get_max_threads() );
#endif

  if ((local_domain.size() > 0) && (rangeObjs.size() > 0)) {

    //
    //  Need to convert range objects to actual box type objects for proximity search
    //

    using rangeValueType = typename RangeObjType::value_type;
    using RangeBox       = stk::search::Box<rangeValueType>;

    std::vector<RangeBox> rangeBoxes;
    rangeBoxes.reserve( rangeObjs.size() );

    for(auto& p : rangeObjs) {
      auto& obj = p;
      rangeBoxes.emplace_back(RangeBox(obj.get_x_min(), obj.get_y_min(), obj.get_z_min(),
                                       obj.get_x_max(), obj.get_y_max(), obj.get_z_max()));

    }

    const stk::search::ProximitySearchTree_T<RangeBox> proxSearch(rangeBoxes);
    const unsigned numBoxDomain = local_domain.size();

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
    {
      //
      //  Set the known return vector sizes
      //

      std::vector<int> overlapList;
#ifdef _OPENMP
      std::vector<std::pair<DomainIdentifier, RangeIdentifier> >& interList   = threadLocalSearchResults[omp_get_thread_num()];
#else
      std::vector<std::pair<DomainIdentifier, RangeIdentifier> >& interList   = searchResults;
#endif
      //
      //  Create an array to store interactions returned by the recursive search routines.  There are at maximum
      //  N interactions per object when searching N objects
      //

      //
      //  Loop over all boxAs in group1 and search them against those objects in group2
      //
#ifdef _OPENMP
#pragma omp for
#endif
      for(unsigned int iboxDomain = 0; iboxDomain < numBoxDomain; ++iboxDomain) {
        proxSearch.SearchForOverlap(local_domain[iboxDomain].first, overlapList);
        for(auto&& jboxRange : overlapList) {
          if(intersects(local_domain[iboxDomain].first, rangeObjs[jboxRange])) {
            if(jboxRange < (int)local_range.size()) {
              interList.emplace_back( local_domain[iboxDomain].second, local_range[jboxRange].second );
            } else {
              interList.emplace_back( local_domain[iboxDomain].second, rangeGhostIdentifiers[jboxRange-local_range.size()] );
            }
          }
        }
      }
    }
  }
#ifdef _OPENMP
  stk::search::concatenate_thread_lists(threadLocalSearchResults, searchResults);
#endif

  if(enforceSearchResultSymmetry) {
    stk::search::communicate_vector(comm, searchResults, enforceSearchResultSymmetry);
  }

  if (sortSearchResults) {
    std::sort(searchResults.begin(), searchResults.end());
  }
}



//
//  Most optimal search specific to actual box arguments
//
template <typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RBoxNumType>
inline void coarse_search_kdtree(std::vector< std::pair<DomainObjType, DomainIdentifier> > const & local_domain,
                                 std::vector< std::pair<stk::search::Box<RBoxNumType>,  RangeIdentifier > > const & local_range,
                                 MPI_Comm comm,
                                 std::vector<std::pair<DomainIdentifier, RangeIdentifier> >& searchResults,
                                 bool enforceSearchResultSymmetry=true,
                                 bool sortSearchResults = false)
{
  int num_procs = -1;
  int proc_id   = -1;
  MPI_Comm_rank(comm, &proc_id);
  MPI_Comm_size(comm, &num_procs);

  searchResults.clear();

#ifdef _OPENMP
  std::vector<std::vector<std::pair<DomainIdentifier, RangeIdentifier> > >
      threadLocalSearchResults( omp_get_max_threads() );
#endif

  {
    std::vector<stk::search::Box<RBoxNumType> > rangeBoxes( local_range.size() );

    std::vector<RangeIdentifier> rangeGhostIdentifiers;

    stk::search::ComputeRangeWithGhostsForCoarseSearch(local_domain, local_range,
                                                       num_procs, rangeBoxes, rangeGhostIdentifiers, comm);

    if ((local_domain.size() > 0) && (rangeBoxes.size() > 0)) {

      const stk::search::ProximitySearchTree_T<stk::search::Box<RBoxNumType> > proxSearch(rangeBoxes);
      const unsigned numBoxDomain = local_domain.size();

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
      {
        //
        //  Set the known return vector sizes
        //

        std::vector<int> overlapList;
#ifdef _OPENMP
        std::vector<std::pair<DomainIdentifier, RangeIdentifier> >& interList   = threadLocalSearchResults[omp_get_thread_num()];
#else
        std::vector<std::pair<DomainIdentifier, RangeIdentifier> >& interList   = searchResults;
#endif

        //
        //  Create an array to store interactions returned by the recursive search routines.  There are at maximum
        //  N interactions per object when searching N objects
        //

        //
        //  Loop over all boxAs in group1 and search them against those objects in group2
        //
        interList.reserve((numBoxDomain*3)/2);
#ifdef _OPENMP
#pragma omp for
#endif
        for(unsigned int iboxDomain = 0; iboxDomain < numBoxDomain; ++iboxDomain) {
          proxSearch.SearchForOverlap(local_domain[iboxDomain].first, overlapList);
          for(auto&& jboxRange : overlapList) {
            if(jboxRange < (int)local_range.size()) {
              interList.emplace_back( local_domain[iboxDomain].second, local_range[jboxRange].second );
            } else {
              interList.emplace_back( local_domain[iboxDomain].second, rangeGhostIdentifiers[jboxRange-local_range.size()] );
            }
          }
        }
      }
    }
  }
#ifdef _OPENMP
  stk::search::concatenate_thread_lists(threadLocalSearchResults, searchResults);
#endif
  if(enforceSearchResultSymmetry) {
    stk::search::communicate_vector(comm, searchResults, enforceSearchResultSymmetry);
  }

  if (sortSearchResults) {
    std::sort(searchResults.begin(), searchResults.end());
  }
}

template <typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeObjType>
inline void coarse_search_kdtree_driver(std::vector< std::pair<DomainObjType, DomainIdentifier> > const & local_domain,
                                        std::vector< std::pair<RangeObjType,  RangeIdentifier > > const & local_range,
                                        MPI_Comm comm,
                                        std::vector<std::pair<DomainIdentifier, RangeIdentifier> >& searchResults,
                                        bool enforceSearchResultSymmetry = true,
                                        bool sortSearchResults = false)
{
  const size_t local_sizes[2] = {local_domain.size(), local_range.size()};
  size_t global_sizes[2];
  all_reduce_sum(comm, local_sizes, global_sizes, 2);
  const bool domain_has_more_boxes = (global_sizes[0] >= global_sizes[1]);
  if(domain_has_more_boxes)
  {
    coarse_search_kdtree(local_domain, local_range, comm, searchResults, enforceSearchResultSymmetry, sortSearchResults);
  }
  else
  {
    std::vector<std::pair<RangeIdentifier, DomainIdentifier> > tempSearchResults;
    coarse_search_kdtree(local_range, local_domain, comm, tempSearchResults, enforceSearchResultSymmetry, sortSearchResults);
    const int p_rank = stk::parallel_machine_rank(comm);
    searchResults.reserve(tempSearchResults.size());
    for(size_t i=0; i<tempSearchResults.size(); ++i)
    {
      size_t idx = tempSearchResults.size() - i - 1;
      auto&& pair = tempSearchResults[idx];
      if (enforceSearchResultSymmetry || get_proc<DomainIdentifier>()(pair.second) == p_rank)
        searchResults.emplace_back(pair.second, pair.first);
    }
  
    if (sortSearchResults) {
      std::sort(searchResults.begin(), searchResults.end());
    }
  }
}

}

#endif
