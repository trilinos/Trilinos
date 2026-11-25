// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_MeshSurface_h
#define Akri_MeshSurface_h

#include <Akri_MeshFacetedSurface.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Types.hpp>

#include <memory>
#include <vector>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class MetaData; } }

namespace krino {

// The active, locally owned sides of the mesh selected by surfaceSelector as a connected
// MeshFacetedSurface whose vertices are their nodes, rebuilt from the current coordinates in every
// prepare_to_compute (so it follows mesh refinement and motion).  The sides of each part of a union
// selector are oriented by the blocks touching it (see should_side_be_flipped); sign = -1 reverses
// them.  Triangles and segments are facets; quadratic sides are split at their mid nodes, and
// quadrilaterals into 2 triangles if planar, otherwise 4 about their centroid.
//
// With skinnedBlocks, the surface is instead the skin of those blocks (see build_skin_mesh_surface).
class MeshSurface : public MeshFacetedSurface
{
public:
  MeshSurface(const stk::mesh::MetaData & meta,
      const stk::mesh::Field<double> & coordsField,
      const stk::mesh::Selector & surfaceSelector,
      const int sign,
      const stk::mesh::PartVector & skinnedBlocks = {});

protected:
  virtual void build_local_facets(const std::vector<BoundingBox> & procQueryBboxes,
      VertexIdsAndCoords & vertexIdsAndCoords,
      FacetVertexIds & facetVertexIds) override;

private:
  const stk::mesh::MetaData & myMeta;
  const stk::mesh::Field<double> & myCoordsField;
  const stk::mesh::Selector mySurfaceSelector;
  int mySign;
  stk::mesh::PartVector mySkinnedBlocks;
};

// Whether a side's nodes are reversed when making a facet from it, so that facets are oriented
// consistently with respect to the volume elements selected by volumeSelector (or keep the side's
// own orientation if volumeSelector is null).  A negative sign reverses the orientation.
bool should_side_be_flipped(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & volumeSelector,
    const stk::mesh::Entity side,
    const int sign);

std::unique_ptr<MeshSurface> build_mesh_surface(const stk::mesh::MetaData & meta,
    const stk::mesh::Field<double> & coordsField,
    const stk::mesh::Selector & surfaceSelector,
    const int sign);

// The skin of the blocks: the sides of their elements that are not between two of them, each facing
// out of the blocks (sign = -1 reverses them), so that a side between one of the blocks and another
// block faces out of the one.  The sides must already exist, as made by
// create_exposed_block_boundary_sides and create_interior_block_boundary_sides.  A side is between
// two of the blocks if it is in both of their parts, which needs no elements from other processors.
std::unique_ptr<MeshSurface> build_skin_mesh_surface(const stk::mesh::MetaData & meta,
    const stk::mesh::Field<double> & coordsField,
    const stk::mesh::PartVector & skinnedBlocks,
    const int sign);

// Appends the facets of the skin of the blocks (as in build_skin_mesh_surface), from the sides selected
// by ownedSideSelector, with the vertex ids and coordinates they use.
void append_skin_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Field<double> & coordsField,
    const stk::mesh::PartVector & skinnedBlocks,
    const stk::mesh::Selector & ownedSideSelector,
    const int sign,
    FacetMesh::VertexIdsAndCoords & vertexIdsAndCoords,
    FacetMesh::FacetVertexIds & facetVertexIds);

// Whether the sides in the bucket are in the parts of at least two of the blocks.
bool are_sides_between_two_of_blocks(const stk::mesh::Bucket & sideBucket, const stk::mesh::PartVector & blocks);

} // namespace krino

#endif // Akri_MeshSurface_h
