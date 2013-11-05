/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_SEARCH_COARSE_SEARCH_2_HPP
#define STK_SEARCH_COARSE_SEARCH_2_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/ReportHandler.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/timer.hpp>
#include <boost/utility/enable_if.hpp>

#include <vector>
#include <utility>

namespace stk {
namespace search {

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
MPI_Datatype get_mpi_type(double)
{
  return MPI_DOUBLE;
}

inline
MPI_Datatype get_mpi_type(float)
{
  return MPI_FLOAT;
}

template <typename Ident>
struct get_proc;

template <typename T>
struct get_proc<std::pair<T, int> >
{
  int operator()(std::pair<T, int> const& id) const
  {
    return id.second;
  }
};

}

template <typename CoordType, int Dimension>
inline
bool invalid_box(CoordType *raw_box)
{
  bool retval = false;
  for (int i = 0; i < Dimension; ++i)  {
    retval |= ((raw_box[i + Dimension] - raw_box[i]) < 0);
  }
  return retval;
};

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

  MPI_Allgather(local_box, data_per_proc, impl::get_mpi_type(coordinate_t()),
                &*recv.begin(), data_per_proc, impl::get_mpi_type(coordinate_t()), comm);

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

template <typename DomainBox, typename DomainIdent, typename RangeBox, typename RangeIdent>
void coarse_search2( std::vector<std::pair<DomainBox, DomainIdent> > const& local_domain,
                     std::vector<std::pair<RangeBox, RangeIdent> > const& local_range,
                     stk::ParallelMachine comm,
                     std::vector<std::pair<DomainIdent, RangeIdent> >& output)
{
  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  const unsigned MaxVolumesPerNode = 16;
  typedef std::pair<DomainBox, DomainIdent> DomainValue;
  typedef std::pair<RangeBox, RangeIdent> RangeValue;
  typedef bgi::rtree< DomainValue, bgi::quadratic<MaxVolumesPerNode> > LocalDomainTree;
  typedef typename LocalDomainTree::bounds_type GlobalBox;
  typedef std::pair<GlobalBox,int> GlobalBoxProc;
  typedef bgi::rtree< GlobalBoxProc, bgi::quadratic<MaxVolumesPerNode> > GlobalDomainTree;
  typedef typename bg::point_type<RangeBox>::type Point;
  typedef std::pair<DomainIdent, RangeIdent> Output;

  LocalDomainTree local_domain_tree(local_domain.begin(), local_domain.end());

  GlobalDomainTree global_domain_tree;
  create_global_spatial_index(global_domain_tree, local_domain_tree.bounds(), comm);

  // outer index is proc
  int size = stk::parallel_machine_size(comm);
  int rank = stk::parallel_machine_rank(comm);
  std::vector<std::vector<RangeValue> > range_send(size);

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
  std::vector<RangeValue> gather_range;
  {
    stk::CommAll comm_all( comm );
    for (int phase = 0; phase < 2; ++phase) {
      for (int p = 0; p < size; ++p) {
        comm_all.send_buffer(p).pack<RangeValue>(&*range_send[p].begin(), range_send[p].size());
      }
      if (phase == 0) {
        comm_all.allocate_buffers( size / 4, false /*not symmetric*/ );
      }
    }

    comm_all.communicate();

    for ( int p = 0 ; p < size ; ++p ) {
      stk::CommBuffer & buf = comm_all.recv_buffer( p );

      const int num_recv = buf.remaining() / sizeof(RangeValue);
      std::vector<RangeValue> values(num_recv);
      gather_range.resize(gather_range.size() + num_recv);
      buf.unpack<RangeValue>( &*gather_range.end() - num_recv, num_recv );

      ThrowRequireMsg(buf.remaining() == 0, buf.remaining());
    }
  }

  // Gather results into output
  {
    stk::CommAll comm_all( comm );
    std::vector<std::vector<Output> > send_matches(size);
    for (size_t r = 0, re = gather_range.size(); r < re; ++r) {
      RangeBox const& range = gather_range[r].first;
      Point center = range.min_corner();
      bg::add_point(center, range.max_corner());
      bg::divide_value(center, 2);
      std::vector<DomainValue> domain_intersections;
      bgi::query(local_domain_tree, bgi::intersects(range), std::back_inserter(domain_intersections));

      if (!domain_intersections.empty()) {
        DomainIdent domain_id = domain_intersections[0].second;
        RangeIdent  range_id  = gather_range[r].second;
        Output temp(domain_id, range_id);
        output.push_back(temp);

        if ((size > 1)
            && (impl::get_proc<DomainIdent>()(domain_id) != rank || impl::get_proc<RangeIdent>()(range_id) != rank)) {
          int other_proc = impl::get_proc<DomainIdent>()(domain_id) == rank ? impl::get_proc<RangeIdent>()(range_id) : impl::get_proc<DomainIdent>()(domain_id);
          send_matches[other_proc].push_back(temp);
        }
      }
    }

    if (size > 1)
    {
      for (int phase = 0; phase < 2; ++phase) {
        for (int p = 0; p < size; ++p) {
          comm_all.send_buffer(p).pack<Output>(&*send_matches[p].begin(), send_matches[p].size());
        }
        if (phase == 0) {
          comm_all.allocate_buffers( size / 4, false /*not symmetric*/ );
        }
      }

      comm_all.communicate();

      for ( int p = 0 ; p < size ; ++p ) {
        stk::CommBuffer & buf = comm_all.recv_buffer( p );

        const int num_recv = buf.remaining() / sizeof(Output);
        std::vector<Output> values(num_recv);
        buf.unpack<Output>( &*values.begin(), num_recv );

        ThrowRequireMsg(buf.remaining() == 0, buf.remaining());

        output.insert(output.end(), values.begin(), values.end());
      }
    }
  }
}

} // namespace search
} // namespace stk

#endif // STK_SEARCH_COARSE_SEARCH_2_HPP
