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
#include <Akri_FacetedSurfaceCalcs.hpp>

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
  STK_ThrowRequire(me != 0 || batch_size <= my_local_facets.size());

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
        if ( facet->does_intersect(proc_bboxes[dest_proc]) )
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
  STK_ThrowAssert( 0 == b.remaining() );
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
    facet->insert_into(my_bounding_box);
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

  my_facet_tree = std::make_unique<SearchTree<Facet*>>( my_all_facets, Facet::get_centroid, Facet::insert_into_bounding_box );

  if ( krinolog.shouldPrint(LOG_DEBUG) )
    krinolog << "P" << stk::EnvData::parallel_rank() << ": After building search tree, storage size is " << storage_size()/(1024.*1024.) << " Mb." << stk::diag::dendl;
}

void Faceted_Surface::prepare_to_compute(const BoundingBox & point_bbox, const double truncation_length)
{
  STK_ThrowAssert(nullptr == my_transformation);
  prepare_to_compute(0.0, point_bbox, truncation_length);
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

  const std::vector<BoundingBox> procPaddedQueryBboxes = fill_processor_bounding_boxes(my_bounding_box, point_bbox, truncation_length);

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
        const BoundingBox & proc_bbox = procPaddedQueryBboxes[dest_proc];
        intersecting_facets.clear();
        if (proc_bbox.intersects(my_bounding_box))
        {
          for ( auto&& facet : my_all_facets )
          {
            if ( facet->does_intersect(proc_bbox) )
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
      STK_ThrowAssert( 0 == b.remaining() );
    }
  }

  // only retain intersecting local descendants
  if (my_all_facets.size() != me_intersecting_facet_counts)
  {
    FacetVec local_descendants;
    local_descendants.reserve(me_intersecting_facet_counts);
    for ( auto&& facet : my_all_facets )
    {
      if ( facet->does_intersect(procPaddedQueryBboxes[me]) )
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
Faceted_Surface::point_distance(const stk::math::Vector3d &x, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance) const
{ /* %TRACE% */  /* %TRACE% */
  STK_ThrowAssertMsg(my_facet_tree, "ERROR: Empty facet tree");

  // get all facets we need to check
  FacetVec nearestFacets;
  my_facet_tree->find_closest_entities( x, nearestFacets, narrow_band_size );
  STK_ThrowRequire( !nearestFacets.empty() || 0.0 != narrow_band_size || my_facet_tree->empty() );

  return point_distance_given_nearest_facets(x, nearestFacets, narrow_band_size, far_field_value, compute_signed_distance);
}

size_t
Faceted_Surface::storage_size() const
{
  size_t store_size = sizeof(Faceted_Surface);
  if (my_facet_tree)
  {
    store_size += my_facet_tree->storage_size();
  }

  store_size += krino::storage_size(my_local_facets);
  store_size += krino::storage_size(my_nonlocal_facets);
  store_size += krino::storage_size(my_all_facets);

  return store_size;
}

} // namespace krino
