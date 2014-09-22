/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_SEARCH_COARSE_SEARCH_BOOST_RTREE_HPP
#define STK_SEARCH_COARSE_SEARCH_BOOST_RTREE_HPP

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/ReportHandler.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/timer.hpp>
#include <boost/math/special_functions/fpclassify.hpp> //for isnan
#include <boost/utility/enable_if.hpp>

#include <stk_search/CoarseSearchGeometryToolkit.hpp>

#include <vector>
#include <utility>


namespace stk { namespace search {

namespace impl {

template <typename Point, unsigned Index = 0, unsigned Dimension = boost::geometry::dimension<Point>::value>
struct fill_point_array
{
  template <typename InputIterator>
  void operator()(Point const& p, InputIterator itr) const
  {
    namespace bg = boost::geometry;

    BOOST_STATIC_ASSERT(Index < Dimension);

    *itr = bg::get<Index>(p);
    ++itr;

    fill_point_array<Point, Index + 1>()(p, itr);
  }
};

template <typename Point, unsigned Dimension>
struct fill_point_array<Point, Dimension, Dimension>
{
  template <typename InputIterator>
  void operator()(Point const&, InputIterator) const
  {}
};

template <typename Point, unsigned Index = 0, unsigned Dimension = boost::geometry::dimension<Point>::value>
struct set_point_impl
{
  template <typename Iterator>
  void operator()(Point & p, Iterator itr) const
  {
    namespace bg = boost::geometry;

    BOOST_STATIC_ASSERT(Index < Dimension);

    bg::set<Index>(p, *itr);
    ++itr;

    set_point_impl<Point, Index + 1>()(p, itr);
  }
};

template <typename Point, unsigned Dimension>
struct set_point_impl<Point, Dimension, Dimension>
{
  template <typename Iterator>
  void operator()(Point & p, Iterator itr) const
  {}
};

template <typename Point, typename InputIterator>
typename boost::enable_if< boost::is_same< typename boost::geometry::traits::tag<Point>::type, boost::geometry::point_tag>, void>::type
fill_array(Point const& point, InputIterator itr)
{
  fill_point_array<Point>()(point, itr);
}

template <typename Box, typename InputIterator>
typename boost::enable_if< boost::is_same< typename boost::geometry::traits::tag<Box>::type, boost::geometry::box_tag>, void>::type
fill_array(Box const& box, InputIterator itr)
{
  namespace bg = boost::geometry;

  fill_array(box.min_corner(), itr);
  fill_array(box.max_corner(), itr + bg::dimension<typename bg::point_type<Box>::type>::value);
}

template <typename Point, typename Iterator>
typename boost::enable_if< boost::is_same< typename boost::geometry::traits::tag<Point>::type, boost::geometry::point_tag>, void>::type
set_point(Point & point, Iterator itr)
{
  set_point_impl<Point>()(point, itr);
}

template <typename Box, typename Iterator>
typename boost::enable_if< boost::is_same< typename boost::geometry::traits::tag<Box>::type, boost::geometry::box_tag>, void>::type
set_box(Box & box, Iterator itr)
{
  namespace bg = boost::geometry;

  set_point(box.min_corner(), itr);
  set_point(box.max_corner(), itr + bg::dimension<typename bg::point_type<Box>::type>::value);
}

inline
ParallelDatatype get_mpi_type(double)
{
#ifdef STK_HAS_MPI
  return MPI_DOUBLE;
#else
  return 0;
#endif
}

inline
ParallelDatatype get_mpi_type(float)
{
#ifdef STK_HAS_MPI
  return MPI_FLOAT;
#else
  return 1;
#endif
}

template <typename RangeBox>
struct IntersectPredicate
{
  RangeBox const& range;

  IntersectPredicate(const RangeBox& x_range) : range(x_range) {}

  // For PointBoundingBox and AxisAlignedBoundingBox, normal boost::intersects is fine
  template <typename DomainBox, typename DomainIdent>
  bool operator()(std::pair<DomainBox,DomainIdent> const& domain) const { return intersects(domain.first,range); }
};

} // namespace impl

// #if defined(__APPLE__) && !defined(isnan)
//   /* for Mac OSX 10, this isnan function may need to be manually declared;
//    * however, on some g++ versions, isnan is a macro so it doesn't need
//    * to be manually declared...
//    */
//   extern "C" int isnan(double value);
// #endif
using boost::math::isnan;

template <typename CoordType, int Dimension>
inline
bool invalid_box(CoordType *raw_box)
{
  bool retval = false;
  for (int i = 0; i < Dimension; ++i)  {
    retval |= isnan(raw_box[i]);
    retval |= isnan(raw_box[i + Dimension]);
    retval |= ((raw_box[i + Dimension] - raw_box[i]) < 0);
  }
  return retval;
}

template <typename Box, typename SpatialIndex>
void // TODO: enable if box matches spatialindex
create_global_spatial_index(SpatialIndex& index, Box const& local_bounding_box, MPI_Comm comm)
{
  namespace bg = boost::geometry;

  typedef typename bg::point_type<Box>::type point_t;
  typedef typename bg::coordinate_type<point_t>::type coordinate_t;

  const unsigned dimension = bg::dimension<point_t>::value;

  const int data_per_proc = dimension * 2;
  coordinate_t local_box[data_per_proc];
  impl::fill_array(local_bounding_box, local_box);

  int size = stk::parallel_machine_size(comm);

  std::vector<coordinate_t> recv(data_per_proc * size, 0);
#ifdef STK_HAS_MPI
  MPI_Allgather(local_box, data_per_proc, impl::get_mpi_type(coordinate_t()),
                &*recv.begin(), data_per_proc, impl::get_mpi_type(coordinate_t()), comm);
#else
  for (int i=0; i < data_per_proc; i++) {
    recv[i] = local_box[i];
  }
#endif
  for (int p = 0; p < size; ++p) {
    coordinate_t *raw_box = &*(recv.begin() + p * data_per_proc);
    if (invalid_box<coordinate_t, bg::dimension<point_t>::value>(raw_box)) {
      continue;
    }
    Box temp;
    impl::set_box(temp, recv.begin() + p * data_per_proc);
    index.insert(std::make_pair(temp, p));
  }
}

// useful for nearest negihbor, not periodc bc
template <typename SpatialIndex>
void create_parallel_domain(SpatialIndex& local_domain, stk::ParallelMachine comm, double box_inflation=0.0)
{
  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  typedef typename SpatialIndex::bounds_type box_type;
  typedef typename bg::point_type<box_type>::type  point_type;
  typedef typename bg::coordinate_type<point_type>::type coordinate_type;
  typedef typename SpatialIndex::value_type value_type;
  typedef std::pair<box_type,int> box_proc_type;
  typedef bgi::rtree< box_proc_type, bgi::quadratic<16> > global_index_type;

  // inflate local box to ensure overlap
  box_type local_box = local_domain.bounds();
  bg::add_value(local_box.min_corner(), -box_inflation);
  bg::add_value(local_box.max_corner(),  box_inflation);

  // compute global domain
  global_index_type global_domain;
  create_global_spatial_index(global_domain, local_box, comm);

  // compute overlaps between global_domain and local box
  std::vector<box_proc_type> intersecting_procs;
  bgi::query(global_domain, bgi::intersects(local_box), std::back_inserter(intersecting_procs));

  // Communicate overlapping domains
  stk::CommAll comm_all( comm );
  int size = stk::parallel_machine_size(comm);
  for (int phase = 0; phase < 2; ++phase) {
    for (size_t i=0, ie=intersecting_procs.size(); i < ie; ++i) {
      box_type box;
      int proc;
      boost::tie(box, proc) = intersecting_procs[i];

      std::vector<value_type> overlaps;
      bgi::query(local_domain, bgi::intersects(box), std::back_inserter(overlaps));

      stk::CommBuffer & buff = comm_all.send_buffer(proc);
      buff.pack<value_type>(&*overlaps.begin(), overlaps.size());
    }
    if (phase == 0) {
      comm_all.allocate_buffers( size / 4, false /*not symmetric*/ );
    }
  }

  comm_all.communicate();

  // Add overlapping processes to local domain
  for ( int p = 0 ; p < size ; ++p ) {
    stk::CommBuffer & buf = comm_all.recv_buffer( p );

    const int num_recv = buf.remaining() / sizeof(value_type);
    std::vector<value_type> values(num_recv);
    buf.unpack<>( &*values.begin(), num_recv );

    ThrowRequireMsg(buf.remaining() == 0, buf.remaining());

    local_domain.insert(values.begin(), values.end());
  }
}

namespace impl {
template <typename DomainBox, typename DomainIdent, typename RangeBox, typename RangeIdent, size_t MaxVolumesPerNode>
int coarse_search_boost_rtree_gather_range( std::vector< std::pair<RangeBox,RangeIdent> > const& local_range,
                                            boost::geometry::index::rtree< std::pair<DomainBox,DomainIdent>,
                                                                           boost::geometry::index::quadratic<MaxVolumesPerNode> > &local_domain_tree,
                                            std::vector<std::pair<RangeBox,RangeIdent> > &gather_range,
                                            stk::ParallelMachine comm)
{
  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  typedef std::pair<RangeBox,RangeIdent>  RangeValue;
  typedef boost::geometry::index::rtree< std::pair<DomainBox,DomainIdent>,
                                        boost::geometry::index::quadratic<MaxVolumesPerNode> > LocalDomainTree;
  typedef typename LocalDomainTree::bounds_type GlobalBox;
  typedef std::pair<GlobalBox,int> GlobalBoxProc;
  typedef bgi::rtree< GlobalBoxProc, bgi::quadratic<MaxVolumesPerNode> > GlobalDomainTree;

  GlobalDomainTree global_domain_tree;
  create_global_spatial_index(global_domain_tree, local_domain_tree.bounds(), comm);

  // outer index is proc
  int p_size = stk::parallel_machine_size(comm);
  std::vector<std::vector<RangeValue> > range_send(p_size);

  // compute what part of range to send to each process
  {
    std::vector<GlobalBoxProc> potential_process_intersections;
    for (size_t i = 0, ie = local_range.size(); i < ie; ++i) {
      potential_process_intersections.clear();
      bgi::query(global_domain_tree, bgi::intersects(local_range[i].first), std::back_inserter(potential_process_intersections));

      for (size_t j = 0, je = potential_process_intersections.size(); j < je; ++j) {
        int proc = potential_process_intersections[j].second;

        range_send[proc].push_back(local_range[i]);
      }
    }
  }

  // gather all range that can potentially intersect my local domain
  {
    stk::CommAll comm_all( comm );
    for (int phase = 0; phase < 2; ++phase) {
      for (int p = 0; p < p_size; ++p) {
        comm_all.send_buffer(p).pack<RangeValue>(&*range_send[p].begin(), range_send[p].size());
      }
      if (phase == 0) {
        comm_all.allocate_buffers( p_size / 4, false /*not symmetric*/ );
      }
    }

    comm_all.communicate();

    for ( int p = 0 ; p < p_size ; ++p ) {
      stk::CommBuffer & buf = comm_all.recv_buffer( p );

      const int num_recv = buf.remaining() / sizeof(RangeValue);
      std::vector<RangeValue> values(num_recv);
      gather_range.resize(gather_range.size() + num_recv);
      buf.unpack<RangeValue>( &*gather_range.end() - num_recv, num_recv );

      ThrowRequireMsg(buf.remaining() == 0, buf.remaining());
    }
  }
  return p_size;
}

} //namespace impl



template <typename DomainBox, typename DomainIdent, typename RangeBox, typename RangeIdent>
void coarse_search_boost_rtree_output_locally( std::vector< std::pair<DomainBox,DomainIdent> > const& local_domain,
                                          std::vector< std::pair<RangeBox,RangeIdent> > const& local_range,
                                          stk::ParallelMachine comm,
                                          std::vector<std::pair<DomainIdent, RangeIdent> >& output
                                         )
{
  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  const unsigned MaxVolumesPerNode = 16;
  typedef std::pair<DomainBox,DomainIdent> DomainValue;
  typedef std::pair<RangeBox,RangeIdent>  RangeValue;
  typedef bgi::rtree< DomainValue, bgi::quadratic<MaxVolumesPerNode> > LocalDomainTree;
  typedef std::pair<DomainIdent, RangeIdent> Output;
  typedef std::vector<Output> OutputVector;

  LocalDomainTree local_domain_tree(local_domain.begin(), local_domain.end());
  std::vector<RangeValue> gather_range;

  impl::coarse_search_boost_rtree_gather_range(local_range, local_domain_tree,
                                               gather_range, comm);

  // Gather results into output
  {
    for (size_t r = 0, re = gather_range.size(); r < re; ++r) {
      RangeBox const& range = gather_range[r].first;
      std::vector<DomainValue> domain_intersections;
      bgi::query(local_domain_tree, bgi::intersects(range) && bgi::satisfies(impl::IntersectPredicate<RangeBox>(range)), std::back_inserter(domain_intersections));

      for (int i = 0, ie = domain_intersections.size(); i < ie; ++i) {
        DomainIdent domain_id = domain_intersections[i].second;
        RangeIdent  range_id  = gather_range[r].second;
        Output temp(domain_id, range_id);
        output.push_back(temp);
      }
    }
  }

  std::sort(output.begin(), output.end());
  typename OutputVector::iterator eitr = std::unique(output.begin(), output.end());
  output.erase(eitr, output.end());
}


template <typename DomainBox, typename DomainIdent, typename RangeBox, typename RangeIdent>
void coarse_search_boost_rtree( std::vector< std::pair<DomainBox,DomainIdent> > const& local_domain,
                                std::vector< std::pair<RangeBox,RangeIdent> > const& local_range,
                                stk::ParallelMachine comm,
                                std::vector<std::pair<DomainIdent, RangeIdent> >& output,
                                bool communicateRangeBoxInfo
                              )
{
  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  const unsigned MaxVolumesPerNode = 16;
  typedef std::pair<DomainBox,DomainIdent> DomainValue;
  typedef std::pair<RangeBox,RangeIdent>  RangeValue;
  typedef bgi::rtree< DomainValue, bgi::quadratic<MaxVolumesPerNode> > LocalDomainTree;
  typedef std::pair<DomainIdent, RangeIdent> Output;
  typedef std::vector<Output> OutputVector;

  LocalDomainTree local_domain_tree(local_domain.begin(), local_domain.end());
  std::vector<RangeValue> gather_range;

  int p_size =
      impl::coarse_search_boost_rtree_gather_range(local_range, local_domain_tree,
                                                   gather_range, comm);

  // Gather results into output
  {
    int p_rank = stk::parallel_machine_rank(comm);

    stk::CommAll comm_all( comm );
    std::vector<OutputVector> send_matches(p_size);
    for (size_t r = 0, re = gather_range.size(); r < re; ++r) {
      RangeBox const& range = gather_range[r].first;
      std::vector<DomainValue> domain_intersections;
      bgi::query(local_domain_tree, bgi::intersects(range) && bgi::satisfies(impl::IntersectPredicate<RangeBox>(range)), std::back_inserter(domain_intersections));

      for (int i = 0, ie = domain_intersections.size(); i < ie; ++i) {
        DomainIdent domain_id = domain_intersections[i].second;
        RangeIdent  range_id  = gather_range[r].second;
        Output temp(domain_id, range_id);
        output.push_back(temp);

        if ((p_size > 1) && ( get_proc<DomainIdent>()(domain_id) != p_rank || (communicateRangeBoxInfo &&  get_proc<RangeIdent>()(range_id) != p_rank )))
        {
          int other_proc = get_proc<DomainIdent>()(domain_id) == p_rank ? get_proc<RangeIdent>()(range_id) : get_proc<DomainIdent>()(domain_id);
          send_matches[other_proc].push_back(temp);
        }
      }
    }

    if (p_size > 1)
    {
      for (int phase = 0; phase < 2; ++phase) {
        for (int p = 0; p < p_size; ++p) {
          comm_all.send_buffer(p).pack<Output>(&*send_matches[p].begin(), send_matches[p].size());
        }
        if (phase == 0) {
          comm_all.allocate_buffers( p_size / 4, false /*not symmetric*/ );
        }
      }

      comm_all.communicate();

      for ( int p = 0 ; p < p_size ; ++p ) {
        stk::CommBuffer & buf = comm_all.recv_buffer( p );

        const int num_recv = buf.remaining() / sizeof(Output);
        OutputVector values(num_recv);
        buf.unpack<Output>( &*values.begin(), num_recv );

        ThrowRequireMsg(buf.remaining() == 0, buf.remaining());

        output.insert(output.end(), values.begin(), values.end());
      }
    }
  }

  std::sort(output.begin(), output.end());
  typename OutputVector::iterator eitr = std::unique(output.begin(), output.end());
  output.erase(eitr, output.end());
}

}} // namespace stk::search

#endif // STK_SEARCH_COARSE_SEARCH_BOOST_RTREE_HPP
