// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_MESH_SURFACE_AKRI_FACETFILEREADER_HPP_
#define KRINO_KRINO_MESH_SURFACE_AKRI_FACETFILEREADER_HPP_

#include <Akri_BoundingBox.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace krino
{

// Reads a facet file as its vertices and then its facets as indices into them.  FAC, PLY and EXO
// files list them that way; an STL file lists each facet's own vertices, which are given 3 per
// facet (or merged, see create_facet_file_reader).  FacetFileSurface reads through it.
// A reader is used once, in order: read_vertices, then read_facet num_facets() times.  Called on
// processor 0 only; errors are thrown, so callers make them collective (with
// read_on_proc0_and_check_on_all).
class FacetFileReader
{
public:
  virtual ~FacetFileReader() = default;

  // Opens the file and reads the vertices: their coordinates (scaled), and a global id for each.
  virtual void read_vertices(
      std::vector<stk::math::Vector3d> & coords, std::vector<stk::mesh::EntityId> & ids) = 0;
  // Valid after read_vertices.
  virtual size_t num_facets() const = 0;
  // The next facet, as indices into the vertices, with its orientation reversed when sign is -1.
  virtual std::array<size_t, 3> read_facet() = 0;
};

// Whether the format (upper case) is FAC, PLY, EXO or STL.
bool is_facet_file_format(const std::string & format);

// Collective.  Runs read on processor 0, then throws on every processor if it threw there, so that
// a bad file stops all processors rather than leaving the others waiting.
void read_on_proc0_and_check_on_all(const stk::ParallelMachine comm,
    const std::string & description,
    const std::string & filename,
    const std::function<void(void)> & read);

// Collective.  The batches in which processor 0 reads and distributes a file's facets: numFacets
// (on each processor) are summed, and each batch has at least min(2048, total) and at least
// total/numProcs facets, so that no processor receives far more than its share at a time.
void compute_facet_file_read_batch_size(const stk::ParallelMachine comm,
    const size_t numFacets,
    size_t & batchSize,
    size_t & numBatches);

// Collective.  The bounding box of the vertices of the file, read on processor 0.
BoundingBox compute_facet_file_bounding_box(const stk::ParallelMachine comm,
    const std::string & format,
    const std::string & filename,
    const int sign,
    const stk::math::Vector3d & scale);

// The reader for the format (upper case: FAC, PLY, EXO or STL; others throw), scaling the
// coordinates and reversing the facets for sign -1.  With mergeCoincidentVertices, a reader may
// give vertices with exactly the same coordinates as one vertex (the STL reader does, as it reads),
// which its facets then share.
std::unique_ptr<FacetFileReader> create_facet_file_reader(const std::string & format,
    const std::string & filename,
    const int sign,
    const stk::math::Vector3d & scale,
    const bool mergeCoincidentVertices = false);

} // namespace krino

#endif /* KRINO_KRINO_MESH_SURFACE_AKRI_FACETFILEREADER_HPP_ */
