// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_MESH_SURFACE_AKRI_FACETFILESURFACE_HPP_
#define KRINO_KRINO_MESH_SURFACE_AKRI_FACETFILESURFACE_HPP_

#include <Akri_MeshFacetedSurface.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <string>
#include <vector>

namespace krino
{

// A MeshFacetedSurface read from a FAC, PLY, EXO or STL file.  The file is read by the first
// prepare_to_compute: processor 0 reads it in batches, sending each facet to the first processor
// p >= 1 whose query box it intersects, else 0, so no processor holds the whole surface.
class FacetFileSurface : public MeshFacetedSurface
{
public:
  // mergeCoincidentVertices gives vertices with exactly the same coordinates one id, so the facets
  // meeting there are connected even when the file lists the vertex more than once (seams). Without
  // it, the surface is not taken to be connected (see MeshFacetedSurface).
  FacetFileSurface(const stk::ParallelMachine comm,
      const std::string & format,
      const std::string & filename,
      const int sign,
      const stk::math::Vector3d & scale,
      const bool mergeCoincidentVertices = true);

  // Collective.  Reads the file's vertices on processor 0.
  virtual void insert_into(BoundingBox & bbox) const override;

protected:
  virtual void build_local_facets(const std::vector<BoundingBox> & procQueryBboxes,
      VertexIdsAndCoords & vertexIdsAndCoords,
      FacetVertexIds & facetVertexIds) override;

private:
  std::string myFormat;
  std::string myFilename;
  int mySign;
  stk::math::Vector3d myScale;
  bool myMergeCoincidentVertices;
};

} // namespace krino

#endif /* KRINO_KRINO_MESH_SURFACE_AKRI_FACETFILESURFACE_HPP_ */
