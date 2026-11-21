// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MeshSurface.hpp>
#include <Akri_DiagWriter.hpp>

#include <stk_util/environment/EnvData.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <array>
#include <utility>

namespace krino{

std::unique_ptr<MeshSurface> build_mesh_surface(const stk::mesh::MetaData & meta,
    const stk::mesh::Field<double> & coordsField,
    const stk::mesh::Selector & surfaceSelector,
    const int sign)
{
  return std::make_unique<MeshSurface>(meta, coordsField, surfaceSelector, sign);
}

static stk::ParallelMachine get_parallel_comm(const stk::mesh::MetaData & meta)
{
  return meta.has_mesh() ? meta.mesh_bulk_data().parallel() : stk::EnvData::parallel_comm();
}

MeshSurface::MeshSurface(const stk::mesh::MetaData & meta,
    const stk::mesh::Field<double> & coordsField,
    const stk::mesh::Selector & surfaceSelector,
    const int sign)
    : MeshFacetedSurface(
          get_parallel_comm(meta), meta.spatial_dimension(), true, FacetMeshOptions{}, true),
      myMeta(meta),
      myCoordsField(coordsField),
      mySurfaceSelector(surfaceSelector),
      mySign(sign)
{
}

static bool determine_polarity_for_side(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & volumeSelector, const stk::mesh::Entity side)
{
  const stk::topology sideTopology = mesh.bucket(side).topology();
  const unsigned numSideElems = mesh.num_elements(side);
  const stk::mesh::Entity * sideElems = mesh.begin_elements(side);
  const stk::mesh::Permutation * sideElemPermutatons = mesh.begin_permutations(side, stk::topology::ELEMENT_RANK);

  for (unsigned iElem = 0; iElem < numSideElems; ++iElem)
    if (volumeSelector(mesh.bucket(sideElems[iElem])))
      return sideTopology.is_positive_polarity(sideElemPermutatons[iElem]);

  // We can get here in parallel without aura.  Check for orientation wrt non-selected elements
  for (unsigned iElem = 0; iElem < numSideElems; ++iElem)
    if (!volumeSelector(mesh.bucket(sideElems[iElem])))
      return !sideTopology.is_positive_polarity(sideElemPermutatons[iElem]);

  STK_ThrowRequireMsg(false, "determine_polarity_for_negative_side_of_interface has no touching element.");
  return false;
}

bool should_side_be_flipped(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & volumeSelector,
    const stk::mesh::Entity side,
    const int sign)
{
  if (volumeSelector.is_null() || determine_polarity_for_side(mesh, volumeSelector, side))
    return sign < 0;
  return sign >= 0;
}

namespace
{

// Appends the facets of sides, as node ids (and, for the centroid of a non-planar quadrilateral, an
// id of its own), collecting the nodes they use.
template <class COORDS_DATA> class SideFacetAppender
{
public:
  SideFacetAppender(const stk::mesh::BulkData & mesh,
      const COORDS_DATA & coordsData,
      const int sign,
      FacetMesh::VertexIdsAndCoords & centroidIdsAndCoords,
      FacetMesh::FacetVertexIds & facetVertexIds)
      : myMesh(mesh),
        myCoordsData(coordsData),
        mySign(sign),
        myDim(mesh.mesh_meta_data().spatial_dimension()),
        myCentroidIdsAndCoords(centroidIdsAndCoords),
        myFacetVertexIds(facetVertexIds)
  {
  }

  void add_selected_side_facets(
      const stk::mesh::Selector & volumeSelector, const stk::mesh::Selector & sideSelector)
  {
    for (auto * bucket : myMesh.get_buckets(myMesh.mesh_meta_data().side_rank(), sideSelector))
    {
      const stk::topology sideTopology = bucket->topology();
      for (auto side : *bucket)
        add_side_facet(volumeSelector, sideTopology, side);
    }
  }

  const std::vector<stk::mesh::Entity> & used_nodes()
  {
    stk::util::sort_and_unique(myUsedNodes);
    return myUsedNodes;
  }

  stk::math::Vector3d node_coordinates(const stk::mesh::Entity node) const
  {
    return stk::math::Vector3d(myCoordsData.entity_values(node), myDim);
  }

private:
  void add_side_facet(const stk::mesh::Selector & volumeSelector,
      const stk::topology sideTopology,
      const stk::mesh::Entity side)
  {
    const bool doFlip = should_side_be_flipped(myMesh, volumeSelector, side, mySign);
    const stk::mesh::Entity * sideNodes = myMesh.begin_nodes(side);
    if (stk::topology::TRI_3 == sideTopology)
    {
      add_facet3d(sideNodes, doFlip, 0, 1, 2);
    }
    else if (stk::topology::QUAD_4 == sideTopology)
    {
      add_quad3d(side, sideNodes, doFlip, 0, 1, 2, 3);
    }
    else if (stk::topology::TRI_6 == sideTopology)
    {
      add_facet3d(sideNodes, doFlip, 0, 3, 5);
      add_facet3d(sideNodes, doFlip, 1, 4, 3);
      add_facet3d(sideNodes, doFlip, 2, 5, 4);
      add_facet3d(sideNodes, doFlip, 3, 4, 5);
    }
    else if (stk::topology::LINE_2 == sideTopology)
    {
      add_facet2d(sideNodes, doFlip, 0, 1);
    }
    else if (stk::topology::LINE_3 == sideTopology)
    {
      add_facet2d(sideNodes, doFlip, 0, 2);
      add_facet2d(sideNodes, doFlip, 2, 1);
    }
    else
    {
      STK_ThrowErrorMsg("Elements with side topology "
          << sideTopology.name() << " not supported for mesh surface initialization.");
    }
  }

  FacetVertexId node_vertex_id(const stk::mesh::Entity node)
  {
    myUsedNodes.push_back(node);
    return myMesh.identifier(node);
  }

  void add_facet2d(const stk::mesh::Entity * sideNodes, const bool doFlip, unsigned p0, unsigned p1)
  {
    STK_ThrowAssert(2 == myDim);
    if (doFlip) std::swap(p0, p1); // permute to flip normal direction
    myFacetVertexIds.push_back({{node_vertex_id(sideNodes[p0]), node_vertex_id(sideNodes[p1]), 0}});
  }

  void add_facet3d(
      const stk::mesh::Entity * sideNodes, const bool doFlip, unsigned p0, unsigned p1, unsigned p2)
  {
    STK_ThrowAssert(3 == myDim);
    if (doFlip) std::swap(p1, p2); // permute to flip normal direction
    myFacetVertexIds.push_back({{node_vertex_id(sideNodes[p0]),
        node_vertex_id(sideNodes[p1]),
        node_vertex_id(sideNodes[p2])}});
  }

  void add_quad3d(const stk::mesh::Entity side,
      const stk::mesh::Entity * sideNodes,
      const bool doFlip,
      unsigned p0,
      unsigned p1,
      unsigned p2,
      unsigned p3)
  {
    STK_ThrowAssert(3 == myDim);
    if (doFlip) std::swap(p1, p3); // permute to flip normal direction

    const std::array<stk::mesh::Entity, 4> nodes{
        {sideNodes[p0], sideNodes[p1], sideNodes[p2], sideNodes[p3]}};
    std::array<stk::math::Vector3d, 4> pt;
    std::array<FacetVertexId, 4> id;
    for (unsigned n = 0; n < 4; ++n)
    {
      pt[n] = node_coordinates(nodes[n]);
      id[n] = node_vertex_id(nodes[n]);
    }

    // scale is the RMS of the diagonal lengths, sqr_scale is scale*scale, or the average of the
    // lengths squared
    const double sqr_scale =
        0.5 * ((pt[0] - pt[2]).length_squared() + (pt[1] - pt[3]).length_squared());
    // tol is 1e-4, sqr_tol is 1e-8
    const double sqr_tol = 1.e-8;
    // Test is that the distance between the midpoints of the 2 diagonals is less than tol*scale
    const bool is_planar =
        (0.5 * (pt[0] + pt[2]) - 0.5 * (pt[1] + pt[3])).length_squared() < sqr_tol * sqr_scale;

    if (is_planar)
    {
      myFacetVertexIds.push_back({{id[0], id[1], id[2]}});
      myFacetVertexIds.push_back({{id[0], id[2], id[3]}});
    }
    else
    {
      // The centroid is not a node: its id is the side's id with the top bit set, so it is unique
      // and cannot be a node id.
      constexpr FacetVertexId centroidBit = FacetVertexId(1) << 63;
      const FacetVertexId sideId = myMesh.identifier(side);
      STK_ThrowRequireMsg(sideId < centroidBit,
          "Side id " << sideId << " is too large for a mesh surface quadrilateral centroid id.");
      const FacetVertexId centroidId = centroidBit | sideId;
      myCentroidIdsAndCoords.emplace_back(centroidId, 0.25 * (pt[0] + pt[1] + pt[2] + pt[3]));
      myFacetVertexIds.push_back({{id[0], id[1], centroidId}});
      myFacetVertexIds.push_back({{id[1], id[2], centroidId}});
      myFacetVertexIds.push_back({{id[2], id[3], centroidId}});
      myFacetVertexIds.push_back({{id[3], id[0], centroidId}});
    }
  }

  const stk::mesh::BulkData & myMesh;
  const COORDS_DATA & myCoordsData;
  int mySign;
  unsigned myDim;
  FacetMesh::VertexIdsAndCoords & myCentroidIdsAndCoords;
  FacetMesh::FacetVertexIds & myFacetVertexIds;
  std::vector<stk::mesh::Entity> myUsedNodes;
};

} // namespace

void MeshSurface::build_local_facets(const std::vector<BoundingBox> & /*procQueryBboxes*/,
    VertexIdsAndCoords & vertexIdsAndCoords,
    FacetVertexIds & facetVertexIds)
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  stk::mesh::Selector active_locally_owned_selector =
      (NULL != myMeta.get_part("ACTIVE_CONTEXT_BIT")) ?
      (*myMeta.get_part("ACTIVE_CONTEXT_BIT") & myMeta.locally_owned_part()) :
      stk::mesh::Selector(myMeta.locally_owned_part());

  const auto coordsData = myCoordsField.data();
  SideFacetAppender appender(mesh, coordsData, mySign, vertexIdsAndCoords, facetVertexIds);

  if (mySurfaceSelector.is_all_unions())
  {
    stk::mesh::PartVector surfaceParts;
    mySurfaceSelector.get_parts(surfaceParts);
    for (auto * surfacePart : surfaceParts)
    {
      const stk::mesh::Selector activeOwnedSurfSelector = active_locally_owned_selector & *surfacePart;
      const auto touchingVols = myMeta.get_blocks_touching_surface(surfacePart);
      const stk::mesh::Selector volSelector = stk::mesh::selectUnion(touchingVols);
      krinolog << "Building facets for " << surfacePart->name() << " oriented with respect to " << volSelector << stk::diag::dendl;
      appender.add_selected_side_facets(volSelector, activeOwnedSurfSelector);
    }
  }
  else
  {
    krinolog << "Not considering side part connectivity or polarity for mesh surface with "
                "non-union surface selector "
             << mySurfaceSelector << stk::diag::dendl;
    const stk::mesh::Selector emptySelector;
    const stk::mesh::Selector activeOwnedSurfSelector =
        active_locally_owned_selector & mySurfaceSelector;
    appender.add_selected_side_facets(emptySelector, activeOwnedSurfSelector);
  }

  // The centroids are already in vertexIdsAndCoords; add the nodes.
  for (auto node : appender.used_nodes())
    vertexIdsAndCoords.emplace_back(mesh.identifier(node), appender.node_coordinates(node));
}

} // namespace krino
