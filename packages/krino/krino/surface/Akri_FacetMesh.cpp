// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_FacetMesh.hpp>

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>

namespace krino
{

FacetMesh::FacetMesh(
    const stk::ParallelMachine comm, const unsigned dim, const FacetMeshOptions & options)
    : myComm(comm), myDim(dim), myOptions(options)
{
  STK_ThrowRequireMsg(
      dim == 2 || dim == 3, "FacetMesh supports dimension 2 or 3, not " << dim << ".");
}

unsigned FacetMesh::add_vertex(const FacetVertexId id,
    const stk::math::Vector3d & coords,
    const stk::math::Vector3d & velocity)
{
  const auto [iter, inserted] = myVertexIndexOfId.try_emplace(id, myVertexIds.size());
  if (inserted)
  {
    myVertexIds.push_back(id);
    myVertexCoords.push_back(coords);
    if (myOptions.withVertexVelocities) myVertexVelocities.push_back(velocity);
  }
  return iter->second;
}

void FacetMesh::append_facet_with_own_vertices(
    const stk::math::Vector3d * coords, const stk::math::Vector3d * velocities)
{
  std::array<unsigned, 3> facet{{0, 0, 0}};
  for (unsigned n = 0; n < myDim; ++n)
  {
    facet[n] = myVertexCoords.size();
    myVertexCoords.push_back(coords[n]);
    if (myOptions.withVertexVelocities) myVertexVelocities.push_back(velocities[n]);
  }
  myFacetVertices.push_back(facet);
}

void FacetMesh::set_owned_facets(const VertexIdsAndCoords & vertexIdsAndCoords,
    const FacetVertexIds & facetVertexIds,
    const std::vector<stk::math::Vector3d> & vertexVelocities)
{
  STK_ThrowRequireMsg(
      vertexVelocities.size() == (myOptions.withVertexVelocities ? vertexIdsAndCoords.size() : 0),
      "FacetMesh: " << vertexVelocities.size() << " vertex velocities given for "
                    << vertexIdsAndCoords.size() << " vertices of a facet mesh "
                    << (myOptions.withVertexVelocities ? "with" : "without")
                    << " vertex velocities.");

  // The index in vertexIdsAndCoords (and vertexVelocities) of each vertex id.
  std::unordered_map<FacetVertexId, size_t> givenVertexIndexOfId;
  givenVertexIndexOfId.reserve(vertexIdsAndCoords.size());
  for (size_t i = 0; i < vertexIdsAndCoords.size(); ++i)
    givenVertexIndexOfId.emplace(vertexIdsAndCoords[i].first, i);

  myVertexIds.clear();
  myVertexCoords.clear();
  myVertexVelocities.clear();
  myVertexIndexOfId.clear();
  myFacetVertices.clear();
  myFacetVertices.reserve(facetVertexIds.size());
  for (auto & facetIds : facetVertexIds)
  {
    std::array<unsigned, 3> facet{{0, 0, 0}};
    std::array<stk::math::Vector3d, 3> coords;
    std::array<stk::math::Vector3d, 3> velocities;
    for (unsigned n = 0; n < myDim; ++n)
    {
      const auto indexIter = givenVertexIndexOfId.find(facetIds[n]);
      STK_ThrowRequireMsg(indexIter != givenVertexIndexOfId.end(),
          "FacetMesh: no coordinates given for vertex " << facetIds[n] << ".");
      coords[n] = vertexIdsAndCoords[indexIter->second].second;
      velocities[n] = myOptions.withVertexVelocities ? vertexVelocities[indexIter->second] :
                                                       stk::math::Vector3d::ZERO;
      if (myOptions.sharedVertices) facet[n] = add_vertex(facetIds[n], coords[n], velocities[n]);
    }
    if (myOptions.sharedVertices)
      myFacetVertices.push_back(facet);
    else
      append_facet_with_own_vertices(coords.data(), velocities.data());
  }
  myNumOwnedFacets = myFacetVertices.size();
  myNumOwnedFacetVertices = myVertexCoords.size();
  build_vertex_facets();
}

void FacetMesh::set_owned_facets_from_coordinates(
    const std::vector<stk::math::Vector3d> & facetVertexCoords,
    const std::vector<stk::math::Vector3d> & facetVertexVelocities)
{
  STK_ThrowRequireMsg(facetVertexCoords.size() % myDim == 0,
      "FacetMesh: " << facetVertexCoords.size() << " facet vertex coordinates are not " << myDim
                    << " per facet.");
  STK_ThrowRequireMsg(facetVertexVelocities.size() ==
          (myOptions.withVertexVelocities ? facetVertexCoords.size() : 0),
      "FacetMesh: " << facetVertexVelocities.size() << " facet vertex velocities given for "
                    << facetVertexCoords.size() << " facet vertices of a facet mesh "
                    << (myOptions.withVertexVelocities ? "with" : "without")
                    << " vertex velocities.");

  if (myOptions.sharedVertices)
  {
    // New ids, unique over the processors: the processor in the high bits and a local count.
    const FacetVertexId firstId =
        (static_cast<FacetVertexId>(stk::parallel_machine_rank(myComm)) << 40) + 1;
    VertexIdsAndCoords vertexIdsAndCoords;
    FacetVertexIds facetVertexIds(facetVertexCoords.size() / myDim, {{0, 0, 0}});
    for (size_t i = 0; i < facetVertexCoords.size(); ++i)
    {
      vertexIdsAndCoords.emplace_back(firstId + i, facetVertexCoords[i]);
      facetVertexIds[i / myDim][i % myDim] = firstId + i;
    }
    set_owned_facets(vertexIdsAndCoords, facetVertexIds, facetVertexVelocities);
    return;
  }

  myVertexIds.clear();
  myVertexCoords.clear();
  myVertexVelocities.clear();
  myVertexIndexOfId.clear();
  myFacetVertices.clear();
  myFacetVertices.reserve(facetVertexCoords.size() / myDim);
  myVertexCoords.reserve(facetVertexCoords.size());
  for (size_t i = 0; i < facetVertexCoords.size(); i += myDim)
    append_facet_with_own_vertices(&facetVertexCoords[i],
        myOptions.withVertexVelocities ? &facetVertexVelocities[i] : nullptr);
  myNumOwnedFacets = myFacetVertices.size();
  myNumOwnedFacetVertices = myVertexCoords.size();
}

void FacetMesh::remove_owned_facets(const std::vector<bool> & isRemoved)
{
  STK_ThrowRequire(isRemoved.size() == myNumOwnedFacets);

  if (!myOptions.sharedVertices)
  {
    // Keep the remaining owned facets, each with its own vertices, in order.
    size_t numKeptFacets = 0;
    for (unsigned facet = 0; facet < myNumOwnedFacets; ++facet)
    {
      if (isRemoved[facet]) continue;
      for (unsigned n = 0; n < myDim; ++n)
      {
        myVertexCoords[myDim * numKeptFacets + n] = myVertexCoords[myDim * facet + n];
        if (myOptions.withVertexVelocities)
          myVertexVelocities[myDim * numKeptFacets + n] = myVertexVelocities[myDim * facet + n];
        myFacetVertices[numKeptFacets][n] = myDim * numKeptFacets + n;
      }
      ++numKeptFacets;
    }
    myFacetVertices.resize(numKeptFacets);
    myVertexCoords.resize(myDim * numKeptFacets);
    if (myOptions.withVertexVelocities) myVertexVelocities.resize(myDim * numKeptFacets);
    myNumOwnedFacets = numKeptFacets;
    myNumOwnedFacetVertices = myVertexCoords.size();
    return;
  }

  // Keep the remaining owned facets, then their vertices, in order, renumbering both.
  size_t numKeptFacets = 0;
  std::vector<bool> isVertexUsed(myNumOwnedFacetVertices, false);
  for (unsigned facet = 0; facet < myNumOwnedFacets; ++facet)
  {
    if (isRemoved[facet]) continue;
    for (unsigned n = 0; n < myDim; ++n)
      isVertexUsed[myFacetVertices[facet][n]] = true;
    myFacetVertices[numKeptFacets++] = myFacetVertices[facet];
  }
  myFacetVertices.resize(numKeptFacets);

  std::vector<unsigned> newIndex(myNumOwnedFacetVertices, 0);
  size_t numKeptVertices = 0;
  myVertexIndexOfId.clear();
  for (size_t v = 0; v < myNumOwnedFacetVertices; ++v)
  {
    if (!isVertexUsed[v]) continue;
    newIndex[v] = numKeptVertices;
    myVertexIds[numKeptVertices] = myVertexIds[v];
    myVertexCoords[numKeptVertices] = myVertexCoords[v];
    if (myOptions.withVertexVelocities) myVertexVelocities[numKeptVertices] = myVertexVelocities[v];
    myVertexIndexOfId.emplace(myVertexIds[v], numKeptVertices);
    ++numKeptVertices;
  }
  myVertexIds.resize(numKeptVertices);
  myVertexCoords.resize(numKeptVertices);
  if (myOptions.withVertexVelocities) myVertexVelocities.resize(numKeptVertices);
  for (auto & facet : myFacetVertices)
    for (unsigned n = 0; n < myDim; ++n)
      facet[n] = newIndex[facet[n]];

  myNumOwnedFacets = numKeptFacets;
  myNumOwnedFacetVertices = numKeptVertices;
  build_vertex_facets();
}

size_t FacetMesh::global_num_vertices() const
{
  require_shared_vertices("global_num_vertices");
  // Each id goes to a rendezvous processor (id % numProcs), which counts it once.
  const int numProcs = stk::parallel_machine_size(myComm);
  stk::CommSparse toRendezvous(myComm);
  stk::pack_and_communicate(toRendezvous,
      [&]()
      {
        for (size_t v = 0; v < myNumOwnedFacetVertices; ++v)
          toRendezvous.send_buffer(myVertexIds[v] % numProcs).pack(myVertexIds[v]);
      });

  std::unordered_set<FacetVertexId> rendezvousIds;
  stk::unpack_communications(toRendezvous,
      [&](int proc)
      {
        stk::CommBuffer & buffer = toRendezvous.recv_buffer(proc);
        while (buffer.remaining())
        {
          FacetVertexId id;
          buffer.unpack(id);
          rendezvousIds.insert(id);
        }
      });

  const size_t localCount = rendezvousIds.size();
  size_t globalCount = 0;
  stk::all_reduce_sum(myComm, &localCount, &globalCount, 1);
  return globalCount;
}

std::vector<int> FacetMesh::compute_vertex_owners() const
{
  require_shared_vertices("compute_vertex_owners");
  // Each id goes to a rendezvous processor (id % numProcs), which returns to every processor that
  // sent it the lowest of them.
  const int numProcs = stk::parallel_machine_size(myComm);
  stk::CommSparse toRendezvous(myComm);
  stk::pack_and_communicate(toRendezvous,
      [&]()
      {
        for (size_t v = 0; v < myNumOwnedFacetVertices; ++v)
          toRendezvous.send_buffer(myVertexIds[v] % numProcs).pack(myVertexIds[v]);
      });

  std::unordered_map<FacetVertexId, int> ownerOfRendezvousId;
  std::vector<std::vector<FacetVertexId>> rendezvousIdsFromProc(numProcs);
  stk::unpack_communications(toRendezvous,
      [&](int proc)
      {
        stk::CommBuffer & buffer = toRendezvous.recv_buffer(proc);
        while (buffer.remaining())
        {
          FacetVertexId id;
          buffer.unpack(id);
          rendezvousIdsFromProc[proc].push_back(id);
          const auto [iter, inserted] = ownerOfRendezvousId.try_emplace(id, proc);
          if (!inserted) iter->second = std::min(iter->second, proc);
        }
      });

  stk::CommSparse fromRendezvous(myComm);
  stk::pack_and_communicate(fromRendezvous,
      [&]()
      {
        for (int proc = 0; proc < numProcs; ++proc)
          for (FacetVertexId id : rendezvousIdsFromProc[proc])
            fromRendezvous.send_buffer(proc).pack(id).pack(ownerOfRendezvousId.at(id));
      });

  std::vector<int> owners(myNumOwnedFacetVertices, -1);
  stk::unpack_communications(fromRendezvous,
      [&](int proc)
      {
        stk::CommBuffer & buffer = fromRendezvous.recv_buffer(proc);
        while (buffer.remaining())
        {
          FacetVertexId id;
          int owner;
          buffer.unpack(id).unpack(owner);
          owners[myVertexIndexOfId.at(id)] = owner;
        }
      });
  return owners;
}

BoundingBox FacetMesh::facet_bounding_box(const unsigned facet) const
{
  BoundingBox bbox;
  for (unsigned n = 0; n < myDim; ++n)
    bbox.accommodate(myVertexCoords[myFacetVertices[facet][n]]);
  return bbox;
}

stk::math::Vector3d FacetMesh::facet_centroid(const unsigned facet) const
{
  stk::math::Vector3d centroid = stk::math::Vector3d::ZERO;
  for (unsigned n = 0; n < myDim; ++n)
    centroid += myVertexCoords[myFacetVertices[facet][n]];
  return centroid / myDim;
}

BoundingBox FacetMesh::owned_facets_bounding_box() const
{
  BoundingBox bbox;
  for (size_t v = 0; v < myNumOwnedFacetVertices; ++v)
    bbox.accommodate(myVertexCoords[v]);
  return bbox;
}

void FacetMesh::update_ghost_facets(const std::vector<BoundingBox> & procBboxes)
{
  // Replace the previous ghosts: drop them and the vertices that only they used.
  myFacetVertices.resize(myNumOwnedFacets);
  for (size_t v = myNumOwnedFacetVertices; v < myVertexIds.size(); ++v)
    myVertexIndexOfId.erase(myVertexIds[v]);
  myVertexIds.resize(myNumOwnedFacetVertices);
  myVertexCoords.resize(myNumOwnedFacetVertices);
  if (myOptions.withVertexVelocities) myVertexVelocities.resize(myNumOwnedFacetVertices);

  const int me = stk::parallel_machine_rank(myComm);
  stk::CommSparse commSparse(myComm);
  stk::pack_and_communicate(commSparse,
      [&]()
      {
        for (unsigned facet = 0; facet < myNumOwnedFacets; ++facet)
        {
          const BoundingBox facetBbox = facet_bounding_box(facet);
          for (size_t proc = 0; proc < procBboxes.size(); ++proc)
          {
            if (int(proc) == me || !procBboxes[proc].intersects(facetBbox)) continue;
            stk::CommBuffer & buffer = commSparse.send_buffer(proc);
            for (unsigned n = 0; n < myDim; ++n)
            {
              const unsigned v = myFacetVertices[facet][n];
              if (myOptions.sharedVertices) buffer.pack(myVertexIds[v]);
              buffer.pack(myVertexCoords[v]);
              if (myOptions.withVertexVelocities) buffer.pack(myVertexVelocities[v]);
            }
          }
        }
      });

  stk::unpack_communications(commSparse,
      [&](int proc)
      {
        stk::CommBuffer & buffer = commSparse.recv_buffer(proc);
        while (buffer.remaining())
        {
          std::array<unsigned, 3> facet{{0, 0, 0}};
          std::array<stk::math::Vector3d, 3> coords;
          std::array<stk::math::Vector3d, 3> velocities;
          for (unsigned n = 0; n < myDim; ++n)
          {
            FacetVertexId id = 0;
            velocities[n] = stk::math::Vector3d::ZERO;
            if (myOptions.sharedVertices) buffer.unpack(id);
            buffer.unpack(coords[n]);
            if (myOptions.withVertexVelocities) buffer.unpack(velocities[n]);
            if (myOptions.sharedVertices) facet[n] = add_vertex(id, coords[n], velocities[n]);
          }
          if (myOptions.sharedVertices)
            myFacetVertices.push_back(facet);
          else
            append_facet_with_own_vertices(coords.data(), velocities.data());
        }
      });

  build_vertex_facets();
}

void FacetMesh::build_vertex_facets()
{
  if (!myOptions.sharedVertices) return;
  myVertexFacetOffsets.assign(myVertexIds.size() + 1, 0);
  for (auto & facet : myFacetVertices)
    for (unsigned n = 0; n < myDim; ++n)
      ++myVertexFacetOffsets[facet[n] + 1];
  for (size_t v = 0; v < myVertexIds.size(); ++v)
    myVertexFacetOffsets[v + 1] += myVertexFacetOffsets[v];
  myVertexFacets.resize(myVertexFacetOffsets.back());
  std::vector<unsigned> next(myVertexFacetOffsets.begin(), myVertexFacetOffsets.end() - 1);
  for (unsigned facet = 0; facet < myFacetVertices.size(); ++facet)
    for (unsigned n = 0; n < myDim; ++n)
      myVertexFacets[next[myFacetVertices[facet][n]]++] = facet;
}

void FacetMesh::throw_no_shared_vertices(const char * function) const
{
  STK_ThrowErrorMsg(
      "FacetMesh::" << function << " requires shared vertices, which this facet mesh has not.");
}

void FacetMesh::throw_wrong_facet_dimension(const int facetDim) const
{
  STK_ThrowErrorMsg("FacetMesh::make_facet: a " << facetDim << "D facet requested from a " << myDim
                                                << "D facet mesh.");
}

void FacetMesh::throw_no_vertex_velocities() const
{
  STK_ThrowErrorMsg(
      "FacetMesh::make_facet: a facet with velocity requested from a facet mesh without vertex "
      "velocities.");
}

size_t FacetMesh::storage_size() const
{
  // An unordered_map entry is a node (key, value and next pointer) plus a bucket pointer.
  const size_t idMapEntrySize = sizeof(FacetVertexId) + sizeof(unsigned) + 2 * sizeof(void *);
  return myVertexIds.capacity() * sizeof(FacetVertexId) +
      (myVertexCoords.capacity() + myVertexVelocities.capacity()) * sizeof(stk::math::Vector3d) +
      myVertexIndexOfId.size() * idMapEntrySize +
      myFacetVertices.capacity() * sizeof(std::array<unsigned, 3>) +
      (myVertexFacetOffsets.capacity() + myVertexFacets.capacity()) * sizeof(unsigned);
}

} // namespace krino
