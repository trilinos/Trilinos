// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_FacetFileSurface.hpp>
#include <Akri_Facet.hpp>
#include <Akri_FacetedSurfaceCalcs.hpp>
#include <Akri_FacetFileReader.hpp>

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <algorithm>
#include <unordered_set>

namespace krino
{

FacetFileSurface::FacetFileSurface(const stk::ParallelMachine comm,
    const std::string & format,
    const std::string & filename,
    const int sign,
    const stk::math::Vector3d & scale,
    const bool mergeCoincidentVertices)
    : MeshFacetedSurface(comm, 3, mergeCoincidentVertices),
      myFormat(format),
      myFilename(filename),
      mySign(sign),
      myScale(scale),
      myMergeCoincidentVertices(mergeCoincidentVertices)
{
  STK_ThrowRequireMsg(is_facet_file_format(format),
      "Connected facets cannot be read from format "
          << format << "; supported formats are FAC, PLY, EXO and STL.");
}

void FacetFileSurface::insert_into(BoundingBox & bbox) const
{
  bbox.accommodate(
      compute_facet_file_bounding_box(parallel_comm(), myFormat, myFilename, mySign, myScale));
}

// Each vertex takes the id of the first (in file order) vertex with exactly the same coordinates.
static void give_coincident_vertices_the_same_id(
    const std::vector<stk::math::Vector3d> & coords, std::vector<stk::mesh::EntityId> & ids)
{
  std::vector<size_t> order(coords.size());
  for (size_t i = 0; i < order.size(); ++i)
    order[i] = i;
  std::stable_sort(order.begin(),
      order.end(),
      [&](const size_t a, const size_t b) { return coords[a] < coords[b]; });
  for (size_t i = 1; i < order.size(); ++i)
    if (coords[order[i]] == coords[order[i - 1]]) ids[order[i]] = ids[order[i - 1]];
}

void FacetFileSurface::build_local_facets(const std::vector<BoundingBox> & procQueryBboxes,
    VertexIdsAndCoords & vertexIdsAndCoords,
    FacetVertexIds & facetVertexIds)
{
  const stk::ParallelMachine comm = parallel_comm();
  const int me = stk::parallel_machine_rank(comm);

  std::unique_ptr<FacetFileReader> reader;
  std::vector<stk::math::Vector3d> coords;
  std::vector<stk::mesh::EntityId> ids;
  size_t numFacets = 0;
  read_on_proc0_and_check_on_all(comm,
      "reading points",
      myFilename,
      [&]()
      {
        reader = create_facet_file_reader(
            myFormat, myFilename, mySign, myScale, myMergeCoincidentVertices);
        reader->read_vertices(coords, ids);
        numFacets = reader->num_facets();
        if (myMergeCoincidentVertices) give_coincident_vertices_the_same_id(coords, ids);
      });

  size_t batchSize = 0;
  size_t numBatches = 0;
  compute_facet_file_read_batch_size(comm, numFacets, batchSize, numBatches);

  // Each vertex is added to vertexIdsAndCoords once, however many of this processor's facets use
  // it.
  std::unordered_set<stk::mesh::EntityId> addedVertexIds;
  auto add_vertex = [&](const stk::mesh::EntityId id, const stk::math::Vector3d & vertexCoords)
  {
    if (addedVertexIds.insert(id).second) vertexIdsAndCoords.emplace_back(id, vertexCoords);
  };

  std::vector<std::array<size_t, 3>> batch;
  size_t numRead = 0;
  for (size_t b = 0; b < numBatches; ++b)
  {
    const size_t currentBatchSize =
        (numRead < numFacets) ? std::min(batchSize, numFacets - numRead) : 0;
    batch.clear();
    read_on_proc0_and_check_on_all(comm,
        "reading facets",
        myFilename,
        [&]()
        {
          for (size_t i = 0; i < currentBatchSize; ++i)
            batch.push_back(reader->read_facet());
        });
    numRead += currentBatchSize;

    // Facets with a repeated vertex (after merging) have no area, and are not kept.
    std::erase_if(batch,
        [&](const std::array<size_t, 3> & facet)
        {
          return ids[facet[0]] == ids[facet[1]] || ids[facet[1]] == ids[facet[2]] ||
              ids[facet[2]] == ids[facet[0]];
        });

    std::vector<int> destinations(batch.size());
    for (size_t i = 0; i < batch.size(); ++i)
    {
      const auto & facet = batch[i];
      destinations[i] = find_destination_proc_for_facet(
          procQueryBboxes, Facet3d(coords[facet[0]], coords[facet[1]], coords[facet[2]]));
      if (destinations[i] == me)
      {
        facetVertexIds.push_back({{ids[facet[0]], ids[facet[1]], ids[facet[2]]}});
        for (unsigned n = 0; n < 3; ++n)
          add_vertex(ids[facet[n]], coords[facet[n]]);
      }
    }

    // pack_and_communicate calls the packing function twice (to size, then to fill the buffers), so
    // it must only pack.
    stk::CommSparse commSparse(comm);
    stk::pack_and_communicate(commSparse,
        [&]()
        {
          for (size_t i = 0; i < batch.size(); ++i)
          {
            if (destinations[i] == me) continue;
            stk::CommBuffer & buffer = commSparse.send_buffer(destinations[i]);
            for (unsigned n = 0; n < 3; ++n)
            {
              buffer.pack(ids[batch[i][n]]);
              buffer.pack(coords[batch[i][n]]);
            }
          }
        });

    stk::unpack_communications(commSparse,
        [&](int proc)
        {
          stk::CommBuffer & buffer = commSparse.recv_buffer(proc);
          while (buffer.remaining())
          {
            std::array<stk::mesh::EntityId, 3> facet;
            for (unsigned n = 0; n < 3; ++n)
            {
              stk::math::Vector3d vertexCoords;
              buffer.unpack(facet[n]);
              buffer.unpack(vertexCoords);
              add_vertex(facet[n], vertexCoords);
            }
            facetVertexIds.push_back(facet);
          }
        });
  }

  // The file's vertices (held on processor 0) are no longer needed.
  reader.reset();
  std::vector<stk::math::Vector3d>().swap(coords);
  std::vector<stk::mesh::EntityId>().swap(ids);
}

} // namespace krino
