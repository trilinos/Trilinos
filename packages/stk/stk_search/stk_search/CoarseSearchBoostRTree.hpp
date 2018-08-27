// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
// 

#ifndef STK_SEARCH_COARSE_SEARCH_BOOST_RTREE_HPP
#define STK_SEARCH_COARSE_SEARCH_BOOST_RTREE_HPP

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/BoostRTreeInterface.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <assert.h> // for static_assert
#include <cmath>    // for isnan
#include <vector>   // for vector
#include <utility>

namespace boost {
namespace geometry {
namespace traits {

// traits for stk::search::Box<T>
template <typename T> struct tag< stk::search::Box<T> > { typedef box_tag type; };
template <typename T> struct point_type< stk::search::Box<T> > { typedef stk::search::Point<T> type; };

template <typename T, size_t Index>
struct indexed_access< stk::search::Box<T>, min_corner, Index >
{
  static_assert((Index < 3)," Index is required to be less than 3");
  static inline T const& get( stk::search::Box<T> const& s) { return s.min_corner()[Index]; }
};

template <typename T, size_t Index>
struct indexed_access< stk::search::Box<T>, max_corner, Index >
{
  static_assert((Index < 3)," Index is required to be less than 3");
  static inline T const& get( stk::search::Box<T> const& s) { return s.max_corner()[Index]; }
};

}}} // namespace boost::geometry::traits


namespace boost { namespace geometry { namespace traits {

// traits for stk::search::Point<T>
template <typename T> struct tag< stk::search::Point<T> > { typedef point_tag type; };
template <typename T> struct coordinate_type< stk::search::Point<T> > { typedef T type; };
template <typename T> struct coordinate_system< stk::search::Point<T> > { typedef cs::cartesian type; };
template <typename T> struct dimension< stk::search::Point<T> > : public boost::mpl::int_<3> {};

template <typename T, size_t Index>
struct access< stk::search::Point<T>, Index >
{
  static_assert((Index < 3)," Index is required to be less than 3");
  static inline T const& get( stk::search::Point<T> const& p) { return p[Index]; }
  static inline void set( stk::search::Point<T> const& p, T const& v) { p[Index] = v; }
};

}}} // namespace boost::geometry::traits

namespace boost { namespace geometry { namespace traits {

// traits for stk::search::Sphere<T>
template <typename T> struct tag< stk::search::Sphere<T> > { typedef box_tag type; };
template <typename T> struct point_type< stk::search::Sphere<T> > { typedef stk::search::Point<T> type; };

template <typename T, size_t Index>
struct indexed_access< stk::search::Sphere<T>, min_corner, Index >
{
  static_assert((Index < 3)," Index is required to be less than 3");
  static inline T const get( stk::search::Sphere<T> const& s) { return s.center()[Index] - s.radius(); }
};

template <typename T, size_t Index>
struct indexed_access< stk::search::Sphere<T>, max_corner, Index >
{
  static_assert((Index < 3)," Index is required to be less than 3");
  static inline T const get( stk::search::Sphere<T> const& s) { return s.center()[Index] + s.radius(); }
};

}}} // namespace boost::geometry::traits


namespace stk { namespace search {

namespace impl {

template <typename Point, unsigned Index = 0, unsigned Dimension = boost::geometry::dimension<Point>::value>
struct fill_point_array
{
  template <typename InputIterator>
  void operator()(Point const& p, InputIterator itr) const
  {
    namespace bg = boost::geometry;

    static_assert(Index < Dimension, "Index is required to be less than Dimension");

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

    static_assert(Index < Dimension, "Index is required to be less than Dimension");

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
typename boost::enable_if< std::is_same< typename boost::geometry::traits::tag<Point>::type, boost::geometry::point_tag>, void>::type
fill_array(Point const& point, InputIterator itr)
{
  fill_point_array<Point>()(point, itr);
}

template <typename Box, typename InputIterator>
typename boost::enable_if< std::is_same< typename boost::geometry::traits::tag<Box>::type, boost::geometry::box_tag>, void>::type
fill_array(Box const& box, InputIterator itr)
{
  namespace bg = boost::geometry;

  fill_array(box.min_corner(), itr);
  fill_array(box.max_corner(), itr + bg::dimension<typename bg::point_type<Box>::type>::value);
}

template <typename Point, typename Iterator>
typename boost::enable_if< std::is_same< typename boost::geometry::traits::tag<Point>::type, boost::geometry::point_tag>, void>::type
set_point(Point & point, Iterator itr)
{
  set_point_impl<Point>()(point, itr);
}

template <typename Box, typename Iterator>
typename boost::enable_if< std::is_same< typename boost::geometry::traits::tag<Box>::type, boost::geometry::box_tag>, void>::type
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

template <typename CoordType, int Dimension>
inline
bool invalid_box(CoordType *raw_box)
{
  bool retval = false;
  for (int i = 0; i < Dimension; ++i)  {
    retval |= std::isnan(raw_box[i]);
    retval |= std::isnan(raw_box[i + Dimension]);
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
  stk::CommSparse comm_sparse( comm );
  int size = stk::parallel_machine_size(comm);
  for (int phase = 0; phase < 2; ++phase) {
    for (size_t i=0, ie=intersecting_procs.size(); i < ie; ++i) {
      box_type box;
      int proc;
      boost::tie(box, proc) = intersecting_procs[i];

      std::vector<value_type> overlaps;
      bgi::query(local_domain, bgi::intersects(box), std::back_inserter(overlaps));

      stk::CommBuffer & buff = comm_sparse.send_buffer(proc);
      buff.pack<value_type>(&*overlaps.begin(), overlaps.size());
    }
    if (phase == 0) {
      comm_sparse.allocate_buffers();
    }
  }

  comm_sparse.communicate();

  // Add overlapping processes to local domain
  for ( int p = 0 ; p < size ; ++p ) {
    stk::CommBuffer & buf = comm_sparse.recv_buffer( p );

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
    stk::CommSparse comm_sparse( comm );
    for (int phase = 0; phase < 2; ++phase) {
      for (int p = 0; p < p_size; ++p) {
        comm_sparse.send_buffer(p).pack<RangeValue>(&*range_send[p].begin(), range_send[p].size());
      }
      if (phase == 0) {
        comm_sparse.allocate_buffers();
      }
    }

    comm_sparse.communicate();

    for ( int p = 0 ; p < p_size ; ++p ) {
      stk::CommBuffer & buf = comm_sparse.recv_buffer( p );

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

namespace {
template <typename DomainBox, typename DomainIdent, typename RangeBox, typename RangeIdent>
void coarse_search_boost_rtree_impl( std::vector< std::pair<DomainBox,DomainIdent> > const& local_domain,
                                     std::vector< std::pair<RangeBox,RangeIdent> > const& local_range,
                                     stk::ParallelMachine comm,
                                     std::vector<std::pair<DomainIdent, RangeIdent> >& output,
                                     bool communicateDomainBoxInfo,
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

    stk::CommSparse comm_sparse( comm );
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

        if ((p_size > 1) &&
            ((communicateDomainBoxInfo && get_proc<DomainIdent>()(domain_id) != p_rank ) ||
             (communicateRangeBoxInfo && get_proc<RangeIdent>()(range_id) != p_rank )))
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
          comm_sparse.send_buffer(p).pack<Output>(&*send_matches[p].begin(), send_matches[p].size());
        }
        if (phase == 0) {
          comm_sparse.allocate_buffers();
        }
      }

      comm_sparse.communicate();

      for ( int p = 0 ; p < p_size ; ++p ) {
        stk::CommBuffer & buf = comm_sparse.recv_buffer( p );

        const int num_recv = buf.remaining() / sizeof(Output);
        OutputVector values(num_recv);
        buf.unpack<Output>( &*values.begin(), num_recv );

        ThrowRequireMsg(buf.remaining() == 0, buf.remaining());

        output.insert(output.end(), values.begin(), values.end());
      }
    }
  }
}
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

  // The search implementation needs to communicate all domain bounding boxes that intersect the
  // bounding box for all range bounding boxes on a processor to that processor. In order to minimize
  // the number of bounding boxes that are communicated we use whichever of the domain or range has
  // fewer bounding boxes in it as the "domain" for the call to coarse_search_boost_rtree_impl().
  const size_t local_sizes[2] = {local_domain.size(), local_range.size()};
  size_t global_sizes[2];
  all_reduce_sum(comm, local_sizes, global_sizes, 2);
  const bool domain_has_more_boxes = global_sizes[0] > global_sizes[1];

  if(domain_has_more_boxes)
  {
    coarse_search_boost_rtree_impl(local_domain, local_range, comm, output, true, communicateRangeBoxInfo);
  }
  else
  {
    std::vector<std::pair<RangeIdent, DomainIdent>> temp_output;
    coarse_search_boost_rtree_impl(local_range, local_domain, comm, temp_output, communicateRangeBoxInfo, true);

    const int p_rank = stk::parallel_machine_rank(comm);
    output.reserve(temp_output.size());
    for(auto && pair : temp_output)
    {
      if(communicateRangeBoxInfo || get_proc<DomainIdent>()(pair.second) == p_rank)
        output.emplace_back(pair.second, pair.first);
    }
  }

  std::sort(output.begin(), output.end());
  auto eitr = std::unique(output.begin(), output.end());
  output.erase(eitr, output.end());
}

}} // namespace stk::search

#endif // STK_SEARCH_COARSE_SEARCH_BOOST_RTREE_HPP
