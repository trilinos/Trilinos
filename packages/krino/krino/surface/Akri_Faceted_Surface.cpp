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
#include <Akri_Sign.hpp>

#include <stk_util/environment/EnvData.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

namespace krino{

Faceted_Surface::Faceted_Surface(const std::string & sn)
   : SurfaceThatTakesAdvantageOfNarrowBandAndThereforeMightHaveWrongSign(),
     my_name(sn) {}

static int find_destination_proc_for_facet(const std::vector<BoundingBox> & procBboxes, const Facet * facet)
{
  int numProcs = procBboxes.size();
  for ( int destProc = 1; destProc < numProcs; ++destProc )
    if ( facet->does_intersect(procBboxes[destProc]) )
      return destProc;
  return 0;
}

static void unpack_and_append_facets_from_proc(stk::CommSparse & commSparse, const int recvProc, FacetOwningVec & facetVec)
{
  stk::CommBuffer & b = commSparse.recv_buffer(recvProc);
  if (b.remaining())
  {
    size_t numProcFacets = 0;
    b.unpack(numProcFacets);

    if(krinolog.shouldPrint(LOG_DEBUG))
      krinolog << "P" << stk::EnvData::parallel_rank() << ":" << " Receiving " << numProcFacets << " facets from proc#" << recvProc << stk::diag::dendl;

    for ( size_t n = 0; n < numProcFacets; ++n )
    {
      std::unique_ptr<Facet> facet = Facet::unpack_from_buffer( b );
      facetVec.emplace_back( std::move(facet) );
    }
    STK_ThrowAssert( 0 == b.remaining() );
  }
}

void
Faceted_Surface::parallel_distribute_facets(const size_t batch_size, const std::vector<BoundingBox> & procBboxes)
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
      dest_procs[index] = find_destination_proc_for_facet(procBboxes, my_local_facets[start+index].get());
      ++(proc_facet_counts[dest_procs[index]]);
    }
  }

  // Communication involves two steps, the first one sizes the messages, the second one actually packs and sends the messages
  for ( int comm_step = 0; comm_step < 2; ++comm_step)
  {
    if (me == 0)
    {
      for ( int dest_proc = 0; dest_proc < num_procs; ++dest_proc )
      {
        const size_t numProcFacets = proc_facet_counts[dest_proc];
        if (numProcFacets > 0)
        {
          if (comm_step == 0 && krinolog.shouldPrint(LOG_DEBUG))
            krinolog << "P" << me << ":" << " Packaging " << numProcFacets << " facets for proc#" << dest_proc << stk::diag::dendl;

          stk::CommBuffer & b = comm_sparse.send_buffer(dest_proc);
          b.pack(numProcFacets);

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
    }

    if (comm_step == 0)
      comm_sparse.allocate_buffers();
    else
      comm_sparse.communicate();
  }

  if (me == 0)
  {
    // delete facets that have been sent away (even if they are headed back to 0)
    my_local_facets.erase(my_local_facets.begin()+start, my_local_facets.end());
  }

  unpack_and_append_facets_from_proc(comm_sparse, 0, my_local_facets);
}

void 
Faceted_Surface::prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length)
{
  
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

static void append_intersecting_facets(const BoundingBox & sourceBbox, const std::vector<Facet*> & sourceFacets, const BoundingBox & targetBbox, std::vector<Facet*> & targetFacets)
{
  targetFacets.reserve(targetFacets.size() + sourceFacets.size());
  if (targetBbox.intersects(sourceBbox))
  {
    for ( auto&& facet : sourceFacets )
      if ( facet->does_intersect(targetBbox) )
        targetFacets.push_back(facet);
  }
}

static void retain_intersecting_facets(const BoundingBox & bbox, std::vector<Facet*> & searchFacets)
{
  std::vector<Facet*> targetFacets;
  append_intersecting_facets(bbox, searchFacets, bbox, targetFacets);
  searchFacets.swap(targetFacets);
}

static size_t get_global_num_facets(const stk::ParallelMachine comm, const std::vector<Facet*> & localFacets)
{
  const size_t localNumFacets = localFacets.size();
  size_t globalNumFacets = 0;
  stk::all_reduce_sum(comm, &localNumFacets, &globalNumFacets, 1);
  return globalNumFacets;
}

void
Faceted_Surface::gather_nonlocal_facets(const BoundingBox & point_bbox, const double truncation_length)
{

  // If truncation length is specified, get all of the facets that are within
  // our padded processor's bounding box.
  // To do this,  see if any local facets lie in the padded nodal bounding box
  // of another proc. if so, send them a copy of those facets.

  my_nonlocal_facets.clear();
  
  const stk::ParallelMachine comm = stk::EnvData::parallel_comm();
  const int numProcs = stk::parallel_machine_size(comm);
  if ( numProcs == 1) return;

  const std::vector<BoundingBox> procPaddedQueryBboxes = fill_processor_bounding_boxes(my_bounding_box, point_bbox, truncation_length);

  const int me = stk::parallel_machine_rank(comm);
  const size_t numFacetsPerBatch = 1 + get_global_num_facets(comm, my_all_facets)/numProcs; // min size of 1
  std::vector<Facet*> intersectingFacets;
  std::vector<std::pair<int,size_t>> procAndNumFacets;

  // Perform communication in batches, communicating with as many neighbors as possible until the batch size is reached.
  // Proc p starts by communicating with proc p+1 and incrementing until all neighbors are processed.
  size_t numBatches = 0;
  size_t maxFacetsPerBatch = 0;
  int commPartner = 1;
  while ( stk::is_true_on_any_proc(comm, commPartner < numProcs) )
  {
    intersectingFacets.clear();
    procAndNumFacets.clear();

    while (commPartner < numProcs && intersectingFacets.size() < numFacetsPerBatch)
    {
      const int destProc = (me+commPartner) % numProcs;

      const size_t numPreviousIntersectingFacets = intersectingFacets.size();
      append_intersecting_facets(my_bounding_box, my_all_facets, procPaddedQueryBboxes[destProc], intersectingFacets);
      const size_t procNumIntersectingFacets = intersectingFacets.size() - numPreviousIntersectingFacets;
      procAndNumFacets.emplace_back(destProc, procNumIntersectingFacets);

      if(krinolog.shouldPrint(LOG_DEBUG))
        krinolog << "P" << me << ":" << " Packaging " << procNumIntersectingFacets << " facets for proc#" << destProc << stk::diag::dendl;

      ++commPartner;
    }
    maxFacetsPerBatch = std::max(maxFacetsPerBatch, intersectingFacets.size());

    // Communication involves two steps, the first one sizes the messages, the second one actually packs and sends the messages
    stk::CommSparse commSparse(comm);
    for ( int commStep = 0; commStep < 2; ++commStep)
    {
      size_t iFacet = 0;
      for (const auto & [destProc, numProcFacets] : procAndNumFacets)
      {
        stk::CommBuffer & b = commSparse.send_buffer(destProc);

        b.pack(numProcFacets);
        for ( size_t iProcFacet=0; iProcFacet<numProcFacets; ++iProcFacet )
          intersectingFacets[iFacet++]->pack_into_buffer(b);
      }

      if (commStep == 0)
        commSparse.allocate_buffers();
      else
        commSparse.communicate();
    }

    for(int recvProc=0; recvProc<numProcs; ++recvProc)
    {
      unpack_and_append_facets_from_proc(commSparse, recvProc, my_nonlocal_facets);
    }

    ++numBatches;
    if(krinolog.shouldPrint(LOG_DEBUG))
      krinolog << "P" << me << ":" << " Facet communication after batch #" << numBatches << ", still need to communicate with " << numProcs-commPartner << " neighboring procs." << stk::diag::dendl;
  }

  if(krinolog.shouldPrint(LOG_DEBUG))
    krinolog << "P" << me << ":" << " Facet communication required " << numBatches << " batches of up to " << maxFacetsPerBatch << " facets per batch." << stk::diag::dendl;

  retain_intersecting_facets(procPaddedQueryBboxes[me], my_all_facets);

  // copy pointers to nonlocal facets into my_all_facets
  my_all_facets.reserve(my_all_facets.size()+ my_nonlocal_facets.size());
  for (auto && nonlocal_descendant : my_nonlocal_facets) { my_all_facets.push_back(nonlocal_descendant.get()); }
}

double
Faceted_Surface::point_distance(const stk::math::Vector3d &x, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance) const
{
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

static std::vector<Facet*> find_candidate_surface_facets_for_intersection_with_segment(SearchTree<Facet*> & facetTree, const stk::math::Vector3d & edgePt0, const stk::math::Vector3d & edgePt1)
{
  BoundingBox facetBbox;
  facetBbox.accommodate(edgePt0);
  facetBbox.accommodate(edgePt1);
  facetBbox.pad_epsilon();
  std::vector<Facet*> results;
  facetTree.get_intersecting_entities(facetBbox, results);
  return results;
}

std::pair<int, double> Faceted_Surface::compute_intersection_with_segment(const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double edgeCrossingTol) const
{
  const std::vector<Facet*> candidates = find_candidate_surface_facets_for_intersection_with_segment(*my_facet_tree, pt0, pt1);
  return compute_intersection_between_and_surface_facets_and_edge(candidates, pt0, pt1);
}

} // namespace krino
