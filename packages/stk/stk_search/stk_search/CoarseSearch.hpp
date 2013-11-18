/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_search_CoarseSearch_hpp
#define stk_search_CoarseSearch_hpp

#include <vector>
#include <utility>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/BihTree.hpp>
#include <stk_search/BihTreeParallelOps.hpp>
#include <stk_search/OctTreeOps.hpp>
#include <stk_search/CoarseSearchBoostRTree.hpp>

namespace stk {
namespace search {

template <typename DomainBox, typename RangeBox>
void coarse_search_bihtree( std::vector<DomainBox> const& domain,
                            std::vector<RangeBox>  const& range,
                            stk::ParallelMachine   comm,
                            std::vector< std::pair< typename DomainBox::Key, typename RangeBox::Key> > & intersections
                          )
{
  const unsigned p_size = parallel_machine_size( comm );
  //find the global box
  if (p_size == 1) {
    bih::BihTree<RangeBox> tree(range.begin(), range.end());
    for (unsigned i = 0, ie = domain.size(); i < ie ; ++i) {
      tree.intersect(domain[i], intersections);
    }
  }
  else {
    std::vector<float> global_box(6);
    stk::search::box_global_bounds( comm, domain.size(), &*domain.begin(), range.size(), &*range.begin(), &*global_box.begin() );

    stk::search::oct_tree_bih_tree_proximity_search(
        comm,
        &*global_box.begin(),
        domain.size(),
        &*domain.begin(),
        range.size(),
        &*range.begin(),
        intersections
        );
  }
}

template <typename DomainBox, typename RangeBox>
void coarse_search_octree( std::vector<DomainBox> const& domain,
                           std::vector<RangeBox>  const& range,
                           stk::ParallelMachine   comm,
                           std::vector< std::pair< typename DomainBox::Key, typename RangeBox::Key> > & intersections
                         )
{

  std::vector<float> global_box(6);
  stk::search::box_global_bounds( comm, domain.size(), &*domain.begin(), range.size(), &*range.begin(), &*global_box.begin() );

  stk::search::oct_tree_proximity_search(
      comm,
      &*global_box.begin(),
      domain.size(),
      &*domain.begin(),
      range.size(),
      &*range.begin(),
      intersections
      );
}


enum SearchMethod {
  BOOST_RTREE,
  OCTREE,
  BIHTREE
};

inline
std::ostream& operator<<(std::ostream &out, SearchMethod method)
{
  switch( method )   {
  case BOOST_RTREE: out << "BOOST_RTREE"; break;
  case OCTREE:      out << "OCTREE"; break;
  case BIHTREE:     out << "BIHTREE"; break;
  }
  return out;
}


namespace impl {

template <typename DomainBox, typename RangeBox, SearchMethod Method, size_t Dimension = DomainBox::DIMENSION>
struct coarse_search_impl
{
  typedef std::vector<DomainBox> DomainVector;
  typedef std::vector<DomainBox> RangeVector;
  typedef std::vector< std::pair< typename DomainBox::Key, typename RangeBox::Key> >  IntersectionVector;

  void operator()(DomainVector const&, std::vector<RangeBox> const&, stk::ParallelMachine , IntersectionVector &) const
  {
    ThrowRequireMsg(false, "Coarse search with " << Method << " not implemented for " << Dimension << "D");
  }
};

template <typename DomainBox, typename RangeBox, size_t Dimension>
struct coarse_search_impl<DomainBox,RangeBox,BOOST_RTREE,Dimension>
{
  typedef std::vector<DomainBox> DomainVector;
  typedef std::vector<DomainBox> RangeVector;
  typedef std::vector< std::pair< typename DomainBox::Key, typename RangeBox::Key> >  IntersectionVector;

  void operator()(DomainVector const& domain, std::vector<RangeBox> const& range, stk::ParallelMachine comm, IntersectionVector & intersections) const
  {
    coarse_search_boost_rtree(domain,range,comm,intersections);
  }
};

// currently Octree only supports 3D searches
template <typename DomainBox, typename RangeBox>
struct coarse_search_impl<DomainBox,RangeBox,OCTREE,3>
{
  typedef std::vector<DomainBox> DomainVector;
  typedef std::vector<DomainBox> RangeVector;
  typedef std::vector< std::pair< typename DomainBox::Key, typename RangeBox::Key> >  IntersectionVector;

  void operator()(DomainVector const& domain, std::vector<RangeBox> const& range, stk::ParallelMachine comm, IntersectionVector & intersections) const
  {
    coarse_search_octree(domain,range,comm,intersections);
  }
};

// currently Octree only supports 3D searches
template <typename DomainBox, typename RangeBox>
struct coarse_search_impl<DomainBox,RangeBox,BIHTREE,3>
{
  typedef std::vector<DomainBox> DomainVector;
  typedef std::vector<DomainBox> RangeVector;
  typedef std::vector< std::pair< typename DomainBox::Key, typename RangeBox::Key> >  IntersectionVector;

  void operator()(DomainVector const& domain, std::vector<RangeBox> const& range, stk::ParallelMachine comm, IntersectionVector & intersections) const
  {
    coarse_search_bihtree(domain,range,comm,intersections);
  }
};

} // namespace impl



template <typename DomainBox, typename RangeBox>
void coarse_search( std::vector<DomainBox> const& domain,
                    std::vector<RangeBox>  const& range,
                    SearchMethod                  method,
                    stk::ParallelMachine          comm,
                    std::vector< std::pair< typename DomainBox::Key, typename RangeBox::Key> > & intersections
                  )
{
  switch( method )
  {
  case BOOST_RTREE:
    impl::coarse_search_impl<DomainBox,RangeBox,BOOST_RTREE>()(domain,range,comm,intersections);
    break;
  case OCTREE:
    impl::coarse_search_impl<DomainBox,RangeBox,OCTREE>()(domain,range,comm,intersections);
    break;
  case BIHTREE:
    impl::coarse_search_impl<DomainBox,RangeBox,BIHTREE>()(domain,range,comm,intersections);
    break;
  }
}

} // namespace search
} // namespace stk

#endif // stk_search_CoarseSearch_hpp
