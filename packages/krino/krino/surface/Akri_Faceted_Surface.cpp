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
#include <Akri_String_Function_Expression.hpp>
#include <Akri_SurfaceIntersectionFromSignedDistance.hpp>
#include <krino/mesh_utils/Akri_AllReduce.hpp>

#include <stk_util/environment/EnvData.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace krino{

template<class FACET>
static int find_destination_proc_for_facet(const std::vector<BoundingBox> & procBboxes, const FACET & facet)
{
  int numProcs = procBboxes.size();
  for ( int destProc = 1; destProc < numProcs; ++destProc )
    if ( facet.does_intersect(procBboxes[destProc]) )
      return destProc;
  return 0;
}

template<class FACET>
static void unpack_and_append_facets_from_proc(stk::CommSparse & commSparse, const int recvProc, std::vector<FACET> & facetVec)
{
  stk::CommBuffer & b = commSparse.recv_buffer(recvProc);
  if (b.remaining())
  {
    size_t numProcFacets = 0;
    b.unpack(numProcFacets);

    if(krinolog.shouldPrint(LOG_DEBUG))
      krinolog << "P" << stk::EnvData::parallel_rank() << ":" << " Receiving " << numProcFacets << " facets from proc#" << recvProc << stk::diag::dendl;

    for ( size_t n = 0; n < numProcFacets; ++n )
      facetVec.emplace_back(b);
    STK_ThrowAssert( 0 == b.remaining() );
  }
}

const std::vector<Facet2d> & FacetedSurfaceBase::get_facets_2d() const { return as_derived_type<Facet2d>().get_facets(); }
const std::vector<Facet3d> & FacetedSurfaceBase::get_facets_3d() const { return as_derived_type<Facet3d>().get_facets(); }
std::vector<Facet2d> & FacetedSurfaceBase::get_facets_2d() { return as_derived_type<Facet2d>().get_facets(); }
std::vector<Facet3d> & FacetedSurfaceBase::get_facets_3d() { return as_derived_type<Facet3d>().get_facets(); }

template<class FACET>
void Faceted_Surface<FACET>::swap(FacetedSurfaceBase & rhs)
{
  Faceted_Surface<FACET> * other = dynamic_cast<Faceted_Surface<FACET>*>(&rhs);
  STK_ThrowRequire(nullptr != other);
  myLocalFacets.swap(other->myLocalFacets);
}

template<class FACET>
void Faceted_Surface<FACET>::parallel_distribute_facets(const size_t batch_size, const std::vector<BoundingBox> & procBboxes)
{
  const int num_procs = stk::EnvData::parallel_size();
  if ( num_procs == 1 ) return;

  const int me = stk::EnvData::parallel_rank();
  STK_ThrowRequire(me != 0 || batch_size <= myLocalFacets.size());

  stk::CommSparse comm_sparse(stk::EnvData::parallel_comm());

  std::vector<int> dest_procs;
  std::vector<size_t> proc_facet_counts;
  size_t start = 0;
  if (me == 0)
  {
    dest_procs.resize(batch_size, 0);
    proc_facet_counts.resize(num_procs, 0);
    start = myLocalFacets.size() - batch_size;
    for ( size_t index = 0; index < batch_size; ++index )
    {
      dest_procs[index] = find_destination_proc_for_facet(procBboxes, myLocalFacets[start+index]);
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
              const FACET & facet = myLocalFacets[start+index];
              facet.pack_into_buffer(b);
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
    myLocalFacets.erase(myLocalFacets.begin()+start, myLocalFacets.end());
  }

  unpack_and_append_facets_from_proc<FACET>(comm_sparse, 0, myLocalFacets);
}

template<class FACET>
static void append_intersecting_facets(const BoundingBox & sourceBbox, const std::vector<FACET> & sourceFacets, const BoundingBox & targetBbox, std::vector<const FACET*> & targetFacets)
{
  if (targetBbox.intersects(sourceBbox))
  {
    for ( auto&& facet : sourceFacets )
      if ( facet.does_intersect(targetBbox) )
        targetFacets.push_back(&facet);
  }
}

template<class FACET>
static void store_pointers_to_local_intersecting_facets_and_nonlocal_facets(const BoundingBox & bbox, const std::vector<FACET> & localFacets, const std::vector<FACET> & nonLocalFacets, std::vector<const FACET*> & allFacetPtrs)
{
  allFacetPtrs.clear();
  allFacetPtrs.reserve(localFacets.size()+nonLocalFacets.size());
  append_intersecting_facets(bbox, localFacets, bbox, allFacetPtrs);

  for ( auto&& facet : nonLocalFacets )
    allFacetPtrs.push_back(&facet);
}

template<class FACET>
static size_t get_global_num_facets(const stk::ParallelMachine comm, const std::vector<FACET> & localFacets)
{
  const size_t localNumFacets = localFacets.size();
  size_t globalNumFacets = 0;
  stk::all_reduce_sum(comm, &localNumFacets, &globalNumFacets, 1);
  return globalNumFacets;
}

template<class FACET>
static void fill_intersecting_nonlocal_facets(const BoundingBox & localFacetBBox, const std::vector<BoundingBox> & procPaddedQueryBboxes, const std::vector<FACET> & localFacets, std::vector<FACET> & nonLocalFacetsThatIntersectProcPaddedQueryBbox)
{
  nonLocalFacetsThatIntersectProcPaddedQueryBbox.clear();

  const stk::ParallelMachine comm = stk::EnvData::parallel_comm();
  const int numProcs = stk::parallel_machine_size(comm);
  if ( numProcs == 1) return;

  // If truncation length is specified, get all of the facets that are within
  // our padded processor's bounding box.
  // To do this,  see if any local facets lie in the padded nodal bounding box
  // of another proc. if so, send them a copy of those facets.

  const int me = stk::parallel_machine_rank(comm);
  const size_t numFacetsPerBatch = 1 + get_global_num_facets(comm, localFacets)/numProcs; // min size of 1
  std::vector<const FACET*> intersectingFacets;
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
      append_intersecting_facets(localFacetBBox, localFacets, procPaddedQueryBboxes[destProc], intersectingFacets);
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
      unpack_and_append_facets_from_proc(commSparse, recvProc, nonLocalFacetsThatIntersectProcPaddedQueryBbox);
    }

    ++numBatches;
    if(krinolog.shouldPrint(LOG_DEBUG))
      krinolog << "P" << me << ":" << " Facet communication after batch #" << numBatches << ", still need to communicate with " << numProcs-commPartner << " neighboring procs." << stk::diag::dendl;
  }

  if(krinolog.shouldPrint(LOG_DEBUG))
    krinolog << "P" << me << ":" << " Facet communication required " << numBatches << " batches of up to " << maxFacetsPerBatch << " facets per batch." << stk::diag::dendl;
}

template<class FACET>
void Faceted_Surface<FACET>::prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length)
{
  build_local_facets(point_bbox);

  if (my_transformation != nullptr)
  {
    for ( auto&& facet : myLocalFacets )
    {
      facet.apply_transformation(*my_transformation);
    }
  }

  my_bounding_box.clear();
  for (auto && facet : myLocalFacets)
    facet.insert_into(my_bounding_box);

  const std::vector<BoundingBox> procPaddedQueryBboxes = fill_processor_bounding_boxes(my_bounding_box, point_bbox, truncation_length);

  // Get all remote facets that might be closest to this processors query points
  fill_intersecting_nonlocal_facets(my_bounding_box, procPaddedQueryBboxes, myLocalFacets, myNonLocalFacets);

  const int me = stk::EnvData::parallel_rank();
  store_pointers_to_local_intersecting_facets_and_nonlocal_facets(procPaddedQueryBboxes[me], myLocalFacets, myNonLocalFacets, myAllFacetPtrs);

  if(krinolog.shouldPrint(LOG_DEBUG))
    krinolog << "P" << stk::EnvData::parallel_rank() << ":" << " Building facet tree for " << myAllFacetPtrs.size() << " facets." << stk::diag::dendl;

  my_facet_tree = std::make_unique<SearchTree<const FACET*>>( myAllFacetPtrs, FACET::get_centroid, FACET::insert_into_bounding_box );

  if (krinolog.shouldPrint(LOG_DEBUG))
    krinolog << "P" << stk::EnvData::parallel_rank() << ": After building search tree, storage size is " << storage_size()/(1024.*1024.) << " Mb." << stk::diag::dendl;
}

template<class FACET>
double Faceted_Surface<FACET>::point_distance(const stk::math::Vector3d &x, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance) const
{
  STK_ThrowAssertMsg(my_facet_tree, "ERROR: Empty facet tree");

  // get all facets we need to check
  std::vector<const FACET*> nearestFacets;
  my_facet_tree->find_closest_entities( x, nearestFacets, narrow_band_size );
  STK_ThrowRequire( !nearestFacets.empty() || 0.0 != narrow_band_size || my_facet_tree->empty() );

  return point_distance_given_nearest_facets(x, nearestFacets, narrow_band_size, far_field_value, compute_signed_distance);
}

template<class FACET>
static const FACET * find_closest_facet(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets)
{
  const FACET* nearest = nullptr;
  double minSqrDist = std::numeric_limits<double>::max();
  for ( auto&& facet : nearestFacets )
  {
    const double sqrDist = facet->point_distance_squared(x);
    if (sqrDist < minSqrDist)
    {
      minSqrDist = sqrDist;
      nearest = facet;
    }
  }
  return nearest;
}

template<class FACET>
const FACET * Faceted_Surface<FACET>::get_closest_facet(const stk::math::Vector3d &x) const
{
  STK_ThrowAssertMsg(my_facet_tree, "ERROR: Empty facet tree");

  std::vector<const FACET*> nearestFacets;
  my_facet_tree->find_closest_entities( x, nearestFacets );
  STK_ThrowRequire( !nearestFacets.empty() || my_facet_tree->empty() );

  return find_closest_facet(x, nearestFacets);
}

template<class FACET>
stk::math::Vector3d Faceted_Surface<FACET>::pseudo_normal_at_closest_point(const stk::math::Vector3d &x) const
{
  if (my_facet_tree->empty())
    return stk::math::Vector3d::ZERO;

  std::vector<const FACET*> nearestFacets;
  my_facet_tree->find_closest_entities( x, nearestFacets, 0. );
  STK_ThrowRequire( !nearestFacets.empty() );

  return compute_pseudo_normal(x, nearestFacets);
}

template<class FACET>
stk::math::Vector3d Faceted_Surface<FACET>::closest_point(const stk::math::Vector3d &x) const
{
  if (my_facet_tree->empty())
    return stk::math::Vector3d::ZERO;

  std::vector<const FACET*> nearestFacets;
  my_facet_tree->find_closest_entities( x, nearestFacets, 0. );
  STK_ThrowRequire( !nearestFacets.empty() );

  return compute_closest_point(x, nearestFacets);
}

template<class FACET>
size_t Faceted_Surface<FACET>::storage_size() const
{
  size_t store_size = sizeof(Faceted_Surface<FACET>);
  if (my_facet_tree)
  {
    store_size += my_facet_tree->storage_size();
  }

  store_size += krino::storage_size(myLocalFacets);
  store_size += krino::storage_size(myNonLocalFacets);
  store_size += krino::storage_size(myAllFacetPtrs);

  return store_size;
}

template<class FACET>
static std::vector<const FACET*> find_candidate_surface_facets_for_intersection_with_segment(SearchTree<const FACET*> & facetTree, const stk::math::Vector3d & edgePt0, const stk::math::Vector3d & edgePt1)
{
  BoundingBox facetBbox;
  facetBbox.accommodate(edgePt0);
  facetBbox.accommodate(edgePt1);
  facetBbox.pad_epsilon();
  std::vector<const FACET*> results;
  facetTree.get_intersecting_entities(facetBbox, results);
  return results;
}

template<class FACET>
std::pair<int, double> Faceted_Surface<FACET>::compute_intersection_with_segment(const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double edgeCrossingTol) const
{
  const double dist0 = point_signed_distance(pt0);
  const double dist1 = point_signed_distance(pt1);
  if (sign_change(dist0, dist1))
  {
    constexpr double tightEnoughTolToCalcDistanceAtSignedDistCrossing = 0.001;
    constexpr double looseEnoughTolToHandleSmallGaps = 0.05;
    const std::vector<const FACET*> candidates = find_candidate_surface_facets_for_intersection_with_segment(*my_facet_tree, pt0, pt1);
    const double intersectionLoc = compute_intersection_between_surface_facets_and_edge(candidates, pt0, pt1);
    if (intersectionLoc >= 0.)
      return {sign(dist1), intersectionLoc};
    auto [crossingSign, linearCrossingLoc] = compute_surface_intersection_with_crossed_segment_from_signed_distance(*this, pt0, pt1, dist0, dist1, tightEnoughTolToCalcDistanceAtSignedDistCrossing);
    const stk::math::Vector3d linearCrossingCoords = (1.-linearCrossingLoc)*pt0 + linearCrossingLoc*pt1;
    const double crossingDist = point_signed_distance(linearCrossingCoords);
    if (std::abs(crossingDist) < looseEnoughTolToHandleSmallGaps*(pt1-pt0).length())
      return {sign(dist1), linearCrossingLoc};

//    krinolog << "Kind of surprised " << crossingDist << " " << looseEnoughTolToHandleSmallGaps*(pt1-pt0).length() << " " << (pt1-pt0).length() << " " << dist0 << " " << dist1 << " " << linearCrossingLoc << stk::diag::dendl;
//    krinolog << pt0 << " " << pt1 << stk::diag::dendl;
  }
  return {0, -1.};
}

template<class FACET>
std::string Faceted_Surface<FACET>::print_sizes() const
{
  const stk::ParallelMachine comm = stk::EnvData::parallel_comm();
  unsigned numFacets = size();
  all_reduce_sum(comm, numFacets);

  // iterate local facets and calc area

  double facets_totalArea = 0.0; // total surface area of interface from facets
  double facets_maxArea = -1.0; // area of largest facet on interface
  double facets_minArea = std::numeric_limits<double>::max();  // area of smallest facet on interface

  // loop over facets
  for ( auto&& facet : myLocalFacets )
  {
    double area = facet.facet_area();
    facets_totalArea += area;

    facets_minArea = std::min(area, facets_minArea);
    facets_maxArea = std::max(area, facets_maxArea);
  }

  all_reduce_min(comm, facets_minArea);
  all_reduce_max(comm, facets_maxArea);
  all_reduce_sum(comm, facets_totalArea);

  std::ostringstream out;
  out << "\t " << "Global sizes: { " << numFacets << " facets on }" << std::endl;
  out << "\t " << "Local sizes: { " << size() << " facets on }" << std::endl;
  out << "\t Global areas: { Min = " << facets_minArea << ", Max = " << facets_maxArea
  << ", Total = " << facets_totalArea << " }" << std::endl;
  return out.str();
}

// Explicit template instantiation
template class Faceted_Surface<Facet3d>;
template class Faceted_Surface<Facet2d>;
template class Faceted_Surface<FacetWithVelocity2d>;
template class Faceted_Surface<FacetWithVelocity3d>;

} // namespace krino
