/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_search_CoarseSearch_hpp
#define stk_search_CoarseSearch_hpp

#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/CoarseSearchBoostRTree.hpp>
#include <stk_search/OctTreeOps.hpp>
#include <stk_search/SearchMethod.hpp>

#include <vector>
#include <utility>


namespace stk { namespace search {

inline
std::ostream& operator<<(std::ostream &out, SearchMethod method)
{
  switch( method )   {
  case BOOST_RTREE: out << "BOOST_RTREE"; break;
  case OCTREE:      out << "OCTREE"; break;
  }
  return out;
}


template <typename DomainBox, typename DomainIdent, typename RangeBox, typename RangeIdent>
void coarse_search( std::vector<std::pair<DomainBox,DomainIdent> > const& domain,
                    std::vector<std::pair<RangeBox,RangeIdent> >   const& range,
                    SearchMethod                                          method,
                    stk::ParallelMachine                                  comm,
                    std::vector< std::pair< DomainIdent, RangeIdent> > &  intersections
                  )
{
  switch( method )
  {
  case BOOST_RTREE:
    coarse_search_boost_rtree(domain,range,comm,intersections);
    break;
  case OCTREE:
    coarse_search_octree(domain,range,comm,intersections);
    break;
  }
}

}} // namespace stk::search

#endif // stk_search_CoarseSearch_hpp
