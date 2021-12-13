// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Faceted_Surface.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_Facet.hpp>
#include <Akri_Utility.hpp>

#include <stk_util/environment/EnvData.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>

namespace krino{

Faceted_Surface::Faceted_Surface(const std::string & sn)
   : SurfaceThatTakesAdvantageOfNarrowBandAndThereforeMightHaveWrongSign(),
     my_name(sn) {}

void
Faceted_Surface::parallel_distribute_facets(const size_t batch_size, const std::vector<BoundingBox> & proc_bboxes)
{
  const int num_procs = stk::EnvData::parallel_size();
  if ( num_procs == 1 ) return;

  const int me = stk::EnvData::parallel_rank();
  ThrowRequire(me != 0 || batch_size <= my_local_facets.size());

  stk::CommSparse comm_sparse(stk::EnvData::parallel_comm());

  std::vector<int> dest_procs;
  std::vector<size_t> proc_facet_counts;
  size_t start = 0;
  if (me == 0)
  {
    dest_procs.resize(batch_size, 0);
    proc_facet_counts.resize(num_procs, 0);
    start = my_local_facets.size() - batch_size;
    for ( size_t index = 0; index < batch_size; ++index )
    {
      const Facet * facet = my_local_facets[start+index].get();
      for ( int dest_proc = 1; dest_proc < num_procs; ++dest_proc )
      {
        if ( proc_bboxes[dest_proc].intersects(facet->bounding_box()) )
        {
          dest_procs[index] = dest_proc;
          ++(proc_facet_counts[dest_proc]);
          break;
        }
      }
      if (dest_procs[index] == 0)
      {
        ++(proc_facet_counts[0]);
      }
    }
  }

  // Communication involves two steps, the first one sizes the messages, the second one actually packs and sends the messages
  for ( int comm_step = 0; comm_step < 2; ++comm_step)
  {
    if (me == 0)
    {
      for ( int dest_proc = 0; dest_proc < num_procs; ++dest_proc )
      {
        if (comm_step == 0 && krinolog.shouldPrint(LOG_DEBUG))
        {
          krinolog << "P" << me << ":"
            << " Packaging " << proc_facet_counts[dest_proc]
            << " facets for proc#" << dest_proc << stk::diag::dendl;
        }

        stk::CommBuffer & b = comm_sparse.send_buffer(dest_proc);
        b.pack(proc_facet_counts[dest_proc]);

        for ( size_t index = 0; index < batch_size; ++index )
        {
          if (dest_procs[index] == dest_proc)
          {
            const Facet * facet = my_local_facets[start+index].get();
            facet->pack_into_buffer(b);
          }
        }
      }
    }
    if (comm_step == 0)
    { //allocation step
      comm_sparse.allocate_buffers();
    }
    else
    { //communication step
      comm_sparse.communicate();
    }
  }

  if (me == 0)
  {
    // delete facets that have been sent away (even if they are headed back to 0)
    my_local_facets.erase(my_local_facets.begin()+start, my_local_facets.end());
  }

  // unload, creating locally owned copy of facet
  const int recv_proc = 0;
  stk::CommBuffer & b = comm_sparse.recv_buffer(recv_proc);

  size_t proc_num_facets_recvd = 0;
  b.unpack(proc_num_facets_recvd);

  if(krinolog.shouldPrint(LOG_DEBUG))
    krinolog << "P" << stk::EnvData::parallel_rank() << ":" << " Receiving " << proc_num_facets_recvd << " facets from proc#" << recv_proc << stk::diag::dendl;

  for ( size_t n = 0; n < proc_num_facets_recvd; ++n )
  {
    std::unique_ptr<Facet> facet = Facet::unpack_from_buffer( b );
    my_local_facets.emplace_back( std::move(facet) );
  }
  ThrowAssert( 0 == b.remaining() );
}

void
Faceted_Surface::pack_into_buffer(stk::CommBuffer & b) const
{
  ThrowRuntimeError("pack_into_buffer should not be called for a Faceted_Surface.");
}

void 
Faceted_Surface::prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length)
{ /* %TRACE[ON]% */ Trace trace__("krino::Faceted_Surface::prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length)"); /* %TRACE% */
  
  build_local_facets(point_bbox);

  if (my_transformation != nullptr)
  {
    for ( auto&& facet : my_local_facets )
    {
      facet->apply_transformation(*my_transformation);
    }
  }

  my_bounding_box.clear();
  for (auto && facet : my_local_facets)
  {
    my_bounding_box.accommodate(facet->bounding_box());
  }

  my_all_facets.clear();
  for ( auto&& facet : my_local_facets )
  {
    my_all_facets.push_back(facet.get());
  }

  // Get all remote facets that might be closest to this processors query points
  gather_nonlocal_facets(point_bbox, truncation_length);

  if(krinolog.shouldPrint(LOG_DEBUG))
    krinolog << "P" << stk::EnvData::parallel_rank() << ":" << " Building facet tree for " << my_all_facets.size() << " facets." << stk::diag::dendl;

  my_facet_tree = std::make_unique<SearchTree<Facet*>>( my_all_facets, Facet::get_bounding_box );

  if ( krinolog.shouldPrint(LOG_DEBUG) )
    krinolog << "P" << stk::EnvData::parallel_rank() << ": After building search tree, storage size is " << storage_size()/(1024.*1024.) << " Mb." << stk::diag::dendl;
}

void
Faceted_Surface::gather_nonlocal_facets(const BoundingBox & point_bbox, const double truncation_length)
{ /* %TRACE[ON]% */ Trace trace__("krino::Faceted_Surface::get_nonlocal_descendants(const BoundingBox & point_bbox, const double truncation_length)"); /* %TRACE% */

  // If truncation length is specified, get all of the facets that are within
  // our padded processor's bounding box.
  // To do this,  see if any local facets lie in the padded nodal bounding box
  // of another proc. if so, send them a copy of those facets.

  my_nonlocal_facets.clear();
  
  const int num_procs = stk::EnvData::parallel_size();
  if ( num_procs == 1) return;

  // We estimate the padding based on the uppoer bound of the distance from this processor's nodal bounding box
  // to the other processors bounding box of facets.  This calculation requires that that local
  // descendants are already gathered.

  std::vector<BoundingBox> facet_bboxes;
  BoundingBox::gather_bboxes( my_bounding_box, facet_bboxes );

  double bbox_padding = std::numeric_limits<double>::max();
  for ( int p = 0; p < num_procs; ++p )
  {
    const BoundingBox & proc_facet_bbox = facet_bboxes[p];
    if (proc_facet_bbox.valid())
    {
      const double upperBnd = std::sqrt(proc_facet_bbox.SqrDistUpperBnd(point_bbox));
      bbox_padding = std::min(bbox_padding,upperBnd);
    }
  }
  if (std::numeric_limits<double>::max() == bbox_padding)
  {
    bbox_padding = truncation_length; // only should happen for no facets anywhere
  }
  if (truncation_length > 0.0)
  {
    bbox_padding = std::min(bbox_padding, truncation_length);
  }

  // gather the bounding box sizes for all procs
  BoundingBox local_bbox = point_bbox;
  std::vector<BoundingBox> bboxes;
  local_bbox.pad(bbox_padding);
  BoundingBox::gather_bboxes( local_bbox, bboxes );

  // determine which facets will be sent to each processor and formulate message sizes
  const int me = stk::EnvData::parallel_rank();
  size_t me_intersecting_facet_counts = 0;
  std::vector<Facet*> intersecting_facets;
  intersecting_facets.reserve(my_all_facets.size());

  // Perform communication in stages.  In the nth stage, processor, p,
  // sends facets to processor p+n and receives from p-n.
  // In the 0th stage, each processor count how many of its own facets are
  // are within that processor's nodal bounding box.
  for ( int comm_partner = 0; comm_partner < num_procs; ++comm_partner )
  {
    stk::CommSparse comm_sparse(stk::EnvData::parallel_comm());

    const int dest_proc = (me+comm_partner) % num_procs;

    // Communication involves two steps, the first one sizes the messages, the second one actually packs and sends the messages
    for ( int comm_step = 0; comm_step < 2; ++comm_step)
    {
      if (comm_step == 0)
      {
        const BoundingBox & proc_bbox = bboxes[dest_proc];
        intersecting_facets.clear();
        if (proc_bbox.intersects(my_bounding_box))
        {
          for ( auto&& facet : my_all_facets )
          {
            if ( proc_bbox.intersects(facet->bounding_box()) )
            {
              intersecting_facets.push_back(facet);
            }
          }
        }
        if (dest_proc == me) me_intersecting_facet_counts = intersecting_facets.size();

        if(krinolog.shouldPrint(LOG_DEBUG))
          krinolog << "P" << me << ":" << " Packaging " << intersecting_facets.size() << " facets for proc#" << dest_proc << stk::diag::dendl;
      }

      stk::CommBuffer & b = comm_sparse.send_buffer(dest_proc);
      if (dest_proc != me) // Don't talk to yourself, it's embarrassing
      {
        const size_t intersecting_facets_size = intersecting_facets.size();
        b.pack(intersecting_facets_size);
        for ( auto&& facet : intersecting_facets )
        {
          facet->pack_into_buffer(b);
        }
      }

      if (comm_step == 0)
      { //allocation step
        comm_sparse.allocate_buffers();
      }
      else
      { //communication step
        comm_sparse.communicate();
      }
    }

    // unload, creating local copies of nonlocal facets

    const int recv_proc = (me+num_procs-comm_partner) % num_procs;

    if (recv_proc != me)
    {
      stk::CommBuffer & b = comm_sparse.recv_buffer(recv_proc);

      size_t proc_num_facets_recvd = 0;
      b.unpack(proc_num_facets_recvd);

      if(krinolog.shouldPrint(LOG_DEBUG))
        krinolog << "P" << stk::EnvData::parallel_rank() << ":" << " Receiving " << proc_num_facets_recvd << " facets from proc#" << recv_proc << stk::diag::dendl;

      for ( size_t n = 0; n < proc_num_facets_recvd; ++n )
      {
        std::unique_ptr<Facet> facet = Facet::unpack_from_buffer( b );
        my_nonlocal_facets.emplace_back( std::move(facet) );
      }
      ThrowAssert( 0 == b.remaining() );
    }
  }

  // only retain intersecting local descendants
  if (my_all_facets.size() != me_intersecting_facet_counts)
  {
    FacetVec local_descendants;
    local_descendants.reserve(me_intersecting_facet_counts);
    for ( auto&& facet : my_all_facets )
    {
      if ( local_bbox.intersects(facet->bounding_box()) )
      {
        local_descendants.push_back( facet );
      }
    }
    my_all_facets.swap(local_descendants);
  }

  // copy nonlocal facets into my_descendants
  my_all_facets.reserve(my_all_facets.size()+ my_nonlocal_facets.size());
  for (auto && nonlocal_descendant : my_nonlocal_facets) { my_all_facets.push_back(nonlocal_descendant.get()); }
}

double
Faceted_Surface::point_distance(const Vector3d &x, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance) const
{ /* %TRACE% */  /* %TRACE% */
  ThrowAssertMsg(my_facet_tree, "ERROR: Empty facet tree");

  // get all facets we need to check
  FacetVec facets;
  my_facet_tree->find_closest_entities( x, facets, narrow_band_size );

  if (facets.empty())
  {
    ThrowRequire( 0.0 != narrow_band_size || my_facet_tree->empty() );
    return far_field_value;
  }

  double dist = 0.0;
  if (compute_signed_distance)
  {
    dist = compute_point_to_facets_distance_by_average_normal(x, facets);
    if (0.0 != narrow_band_size && std::abs(dist) > narrow_band_size)
    {
      dist = far_field_value;
    }
  }
  else
  {
    double min_sqr_dist = std::numeric_limits<double>::max();
    for ( auto&& facet : facets )
    {
      const double sqr_dist = facet->point_distance_squared(x);
      if (sqr_dist < min_sqr_dist)
      {
        min_sqr_dist = sqr_dist;
      }
    }
    if (0.0 != narrow_band_size && min_sqr_dist > narrow_band_size*narrow_band_size)
    {
      dist = far_field_value;
    }
    else
    {
      dist = std::sqrt(min_sqr_dist);
    }
  }

  return dist;
}

double
Faceted_Surface::compute_point_to_facets_distance_by_average_normal(const Vector3d &x, const FacetVec & facets) const
{ /* %TRACE% */  /* %TRACE% */

  // If the closest_point weights are all larger than this value, then the closest point
  // is considered to be on the face of the closest facet rather than on the edges of the facet, and
  // therefore only the closest facet is considered in the distance calculation.  Otherwise, all of the
  // facets are considered to compute an average normal in order to compute the distance.
  const double edge_tol = 1.e-6;

  std::vector<FacetDistanceQuery> facet_queries;
  facet_queries.reserve(facets.size());
  for ( auto&& facet : facets )
  {
    if (facet->degenerate()) continue; // Skip zero-sized facets
    facet_queries.emplace_back(*facet, x);
  }

  ThrowRequireMsg(!facet_queries.empty(), "All facets are degenerate in compute_point_to_facets_distance_by_average_normal.");

  unsigned nearest = 0;
  for ( unsigned index=0; index<facet_queries.size(); ++index )
  {
    if ( facet_queries[index].distance_squared() < facet_queries[nearest].distance_squared() )
    {
      nearest = index;
    }
  }

  if ( facet_queries[nearest].distance_squared() == 0. )
  {
    return 0.0;
  }

  const Vector3d closest_pt_wts = facet_queries[nearest].closest_point_weights();
  const int dim = dynamic_cast<krino::Facet3d *>(facets[nearest]) ? 3 : 2;

  bool closest_point_on_edge = false;
  for (int d=0; d<dim; ++d)
  {
    if (closest_pt_wts[d] < edge_tol)
    {
      closest_point_on_edge = true;
      break;
    }
  }

  if (!closest_point_on_edge) return facet_queries[nearest].signed_distance();

  const double min_sqr_dist = facet_queries[nearest].distance_squared();
  const Vector3d pseudo_normal = compute_pseudo_normal(dim, facet_queries, nearest);

  if (pseudo_normal.length_squared() == 0.0)
  {
    krinolog << "Warning:  Cannot determine the average facet normal for computing the level set distance at point " << x
        << ".  This can happen when faceted facet include coincident facets with opposite normals.  Arbitrarily setting the distance to be positive." << stk::diag::dendl;
    return std::sqrt(min_sqr_dist);
  }
  else
  {
    if (Dot(pseudo_normal, x-facet_queries[nearest].closest_point()) > 0)
    {
      return std::sqrt(min_sqr_dist);
    }
    else
    {
      return -std::sqrt(min_sqr_dist);
    }
  }
}

Vector3d
Faceted_Surface::compute_pseudo_normal(const unsigned dim, const std::vector<FacetDistanceQuery> & facet_queries, const unsigned nearest) const
{
  const double tol = 1.e-6;
  Vector3d pseudo_normal = Vector3d::ZERO;
  Vector3d average_normal = Vector3d::ZERO;

  const Vector3d nearest_closest_point = facet_queries[nearest].closest_point();
  const double nearest_size2 = facet_queries[nearest].facet().bounding_box().SqrSize();
  unsigned close_count = 0;
  for ( auto&& query : facet_queries )
  {
    const Vector3d closest_point = query.closest_point();
    const double dist2_from_nearest = (closest_point-nearest_closest_point).length_squared();
    if (dist2_from_nearest < tol*tol*nearest_size2)
    {
      ++close_count;
      const Facet & facet = query.facet();

      average_normal += facet.facet_normal();

      if (3 == dim)
      {
        const Vector3d closest_pt_wts = query.closest_point_weights();
        const int closest_node = (closest_pt_wts[0] > closest_pt_wts[1]) ? ((closest_pt_wts[0] > closest_pt_wts[2]) ? 0 : 2) : ((closest_pt_wts[1] > closest_pt_wts[2]) ? 1 : 2);

        const int n0 = closest_node;
        const int n1 = (n0<2) ? (n0+1) : 0;
        const int n2 = (n1<2) ? (n1+1) : 0;
        const Vector3d edge0 = facet.facet_vertex(n1) - facet.facet_vertex(n0);
        const Vector3d edge1 = facet.facet_vertex(n2) - facet.facet_vertex(n0);
        const double facet_angle = std::acos(Dot(edge0, edge1)/(edge0.length()*edge1.length()));

        pseudo_normal += facet.facet_normal()*facet_angle;
      }
    }
  }
  ThrowRequireMsg(close_count>0,"Issue with tolerance in compute_pseudo_normal.  No facet found within tolerance of closest point.");

  return (3 == dim && close_count > 2) ? pseudo_normal : average_normal;
}

size_t
Faceted_Surface::storage_size() const
{
  size_t store_size = sizeof(Faceted_Surface);
  if (my_facet_tree)
  {
    store_size += my_facet_tree->storage_size();
  }

  store_size += utility::storage_size(my_local_facets);
  store_size += utility::storage_size(my_nonlocal_facets);
  store_size += utility::storage_size(my_all_facets);

  return store_size;
}

} // namespace krino
