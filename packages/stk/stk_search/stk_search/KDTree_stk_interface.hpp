
#ifndef KDTREE_STK_INTERFACE_H_
#define KDTREE_STK_INTERFACE_H_

#ifdef _OPENMP
#include <omp.h>
#endif


#include <stk_search/OctTreeOps.hpp>
#include "stk_util/environment/WallTime.hpp"
#include "stk_search/KDTree_BoundingBox.hpp"
#include <stk_search/CommonSearchUtil.hpp>
#include <stk_search/Sphere.hpp>


namespace stk {
  namespace search {


   template <typename DomainIdentifier, typename RangeIdentifier, typename DomainBoxType, typename RangeBoxType>
      inline void kdtree_search(std::vector< std::pair<DomainBoxType, DomainIdentifier> > const & local_domain,
                                std::vector< std::pair<RangeBoxType,  RangeIdentifier > > const & local_range,
                                MPI_Comm comm,
                                std::vector<std::pair<DomainIdentifier, RangeIdentifier> >& searchResults,
                                bool communicateRangeBoxInfo=true)
    {
      using value_type = typename DomainBoxType::value_type;
      using Box        = stk::search::Box<value_type>;

      std::vector< std::pair<Box, DomainIdentifier> > domainBoxes;
      domainBoxes.reserve( local_domain.size() );

      for(auto& p : local_domain) {

        auto& boxAble = p.first;

        domainBoxes.push_back( std::make_pair(Box(boxAble.get_x_min(), boxAble.get_y_min(), boxAble.get_z_min(), 
                                                  boxAble.get_x_max(), boxAble.get_y_max(), boxAble.get_z_max()), p.second)
                             );
      }

      std::vector< std::pair<Box,  RangeIdentifier > > rangeBoxes;
      rangeBoxes.reserve( local_range.size() );

      for(auto& p : local_range) {
        auto& boxAble = p.first;
        rangeBoxes.push_back( std::make_pair(Box(boxAble.get_x_min(), boxAble.get_y_min(), boxAble.get_z_min(), 
                                                 boxAble.get_x_max(), boxAble.get_y_max(), boxAble.get_z_max()), p.second)
                            );
      }

      kdtree_search(domainBoxes, rangeBoxes, comm, searchResults, communicateRangeBoxInfo);
    }

    //
    //  Wrapping of the KDTree search for use in the general STK coarse_search algorithm.  Find the overlaps between
    //  boxes in domain and boxes in range
    //
    template <typename DomainIdentifier, typename RangeIdentifier, typename DomainBoxType, typename RangeBoxType>
    inline void kdtree_search(std::vector< std::pair<stk::search::Box<DomainBoxType>, DomainIdentifier> > const & local_domain,
                              std::vector< std::pair<stk::search::Box<RangeBoxType>,  RangeIdentifier > > const & local_range,
                               MPI_Comm comm,
                               std::vector<std::pair<DomainIdentifier, RangeIdentifier> >& searchResults,
                               bool communicateRangeBoxInfo=true)
    {
      int num_procs = -1;
      int proc_id   = -1;
      MPI_Comm_rank(comm, &proc_id);
      MPI_Comm_size(comm, &num_procs);

      std::vector<stk::search::Box<RangeBoxType> > rangeBoxes( local_range.size() );

      std::vector<RangeIdentifier> rangeGhostIdentifiers;

      stk::search::ComputeRangeWithGhostsForCoarseSearch(local_domain, local_range,
                                            num_procs, rangeBoxes, rangeGhostIdentifiers, comm);

#ifdef _OPENMP
      std::vector<std::vector<std::pair<DomainIdentifier, RangeIdentifier> > >
        threadLocalSearchResults( omp_get_max_threads() );
#endif

      if ((local_domain.size() > 0) && (rangeBoxes.size() > 0)) {

        const stk::search::ProximitySearchTree_T<stk::search::Box<RangeBoxType> > proxSearch(rangeBoxes);
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
            for(unsigned ilist = 0; ilist < overlapList.size(); ++ilist) {
              const int jboxRange = overlapList[ilist];
              if(jboxRange < (const int)local_range.size()) {
                interList.emplace_back( local_domain[iboxDomain].second, local_range[jboxRange].second );
              } else {
                interList.emplace_back( local_domain[iboxDomain].second, rangeGhostIdentifiers[jboxRange-local_range.size()] );
              }
            }
          }
        }
      }
#ifdef _OPENMP
      stk::search::ConcatenateThreadLists(threadLocalSearchResults, searchResults);
#endif

      if(communicateRangeBoxInfo) {
        std::vector <std::pair<DomainIdentifier,RangeIdentifier> > tmp;
        tmp.reserve(searchResults.size());
        stk::search::communicateVector(comm, searchResults, tmp, communicateRangeBoxInfo);
        searchResults=tmp;
        std::sort(searchResults.begin(), searchResults.end());
      }
    }

  }
}

#endif
