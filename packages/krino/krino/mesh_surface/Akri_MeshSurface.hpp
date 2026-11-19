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
class MeshSurface : public MeshFacetedSurface
{
public:
  MeshSurface(const stk::mesh::MetaData & meta,
      const stk::mesh::Field<double> & coordsField,
      const stk::mesh::Selector & surfaceSelector,
      const int sign);

protected:
  virtual void build_local_facets(const std::vector<BoundingBox> & procQueryBboxes,
      VertexIdsAndCoords & vertexIdsAndCoords,
      FacetVertexIds & facetVertexIds) override;

private:
  const stk::mesh::MetaData & myMeta;
  const stk::mesh::Field<double> & myCoordsField;
  const stk::mesh::Selector mySurfaceSelector;
  int mySign;
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

} // namespace krino

#endif // Akri_MeshSurface_h
