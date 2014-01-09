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

// EXPERIMENTAL
//
// intersections returned will be those resulting from the (local) domain boxes
// intersecting range boxes from the entire distributed set.  NEEDS TO BE CALLED
// WITH EXPLICIT TEMPLATE ARGUMENTS.
template <typename DomainBox, typename DomainIdent, typename RangeBox, typename RangeIdent>
typename boost::disable_if< boost::mpl::or_<typename get_proc<DomainIdent>::supported,
                                            typename get_proc<RangeIdent>::supported> >::type
coarse_search_nonIdentProc(
                    std::vector<std::pair<DomainBox,DomainIdent> > const& domain,
                    std::vector<std::pair<RangeBox,RangeIdent> >   const& range,
                    SearchMethod                                          method,
                    stk::ParallelMachine                                  comm,
                    std::vector< std::pair< DomainIdent, RangeIdent> > &  intersections
                  )
{
  switch( method )
  {
  case BOOST_RTREE:
    coarse_search_boost_rtree_output_locally<DomainBox, DomainIdent, RangeBox, RangeIdent>(domain,range,comm,intersections);
    break;
  case OCTREE:
    std::cerr << "coarse_search_octree(..) does not support std::search::coarse_search_nonIdentProc(..) yet" << std::endl;
    std::abort();
    // coarse_search_octree(domain,range,comm,intersections);
    break;
  }
}


// THIS MIGHT BE WHAT WE ACTUALLY WANT FOR THE INTERFACE.
template <typename DomainBox, typename DomainIdent, typename RangeBox, typename RangeIdent>
void coarse_search_nonIdentProc(
                    std::vector<std::pair<DomainBox,DomainIdent> > const& domain,
                    std::vector<std::pair<RangeBox, RangeIdent> >  const& range,
                    SearchMethod                                          method,
                    stk::ParallelMachine                                  comm,
                    std::vector< std::pair< IdentProc<DomainIdent, unsigned int>,
                                            IdentProc<RangeIdent, unsigned int> > > &  intersections
                  )
{
  std::cerr << "Future version of coarse_search called" << std::endl;
  abort();
}

// intersections will be those of distributed domain boxes associated with this
// processor rank via get_proc<DomainIdent>(.) that intersect distributed range
// boxes.  Optionally, also include intersections of distributed domain boxes
// with distributed range boxes associated with this processor rank via
// get_proc<RangeIdent>(.).
template <typename DomainBox, typename DomainIdent, typename RangeBox, typename RangeIdent>
typename boost::enable_if< boost::mpl::and_<typename get_proc<DomainIdent>::supported,
                                            typename get_proc<RangeIdent>::supported> >::type
coarse_search( std::vector<std::pair<DomainBox,DomainIdent> > const& domain,
               std::vector<std::pair<RangeBox,RangeIdent> >   const& range,
               SearchMethod                                          method,
               stk::ParallelMachine                                  comm,
               std::vector< std::pair< DomainIdent, RangeIdent> > &  intersections,
               bool communicateRangeBoxInfo=true
             )
{
  switch( method )
  {
  case BOOST_RTREE:
    coarse_search_boost_rtree(domain,range,comm,intersections,communicateRangeBoxInfo);
    break;
  case OCTREE:
    coarse_search_octree(domain,range,comm,intersections,communicateRangeBoxInfo);
    break;
  }
}


}} // namespace stk::search

#endif // stk_search_CoarseSearch_hpp
