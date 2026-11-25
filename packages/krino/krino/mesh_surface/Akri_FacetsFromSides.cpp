/*
 * Akri_FacetsFromSides.cpp
 *
 *  Created on: Sep 8, 2025
 *      Author: drnoble
 */
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_Facet.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshFacetedSurface.hpp>
#include <Akri_OrientedSideNodes.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_util/util/SortAndUnique.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace krino {

// The interface velocity at the node: NP1, or with two states the average of N and NP1.
static stk::math::Vector3d interface_velocity_at_node(const stk::mesh::BulkData & mesh,
    const FieldRef interfaceVelocity,
    const unsigned numVelocityStates,
    const stk::mesh::Entity node)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  const stk::math::Vector3d velocity = get_vector_field(mesh, interfaceVelocity, node, dim);
  if (1 == numVelocityStates) return velocity;
  STK_ThrowAssert(2 == numVelocityStates && interfaceVelocity.number_of_states() > 1);
  return 0.5 *
      (velocity +
          get_vector_field(mesh, interfaceVelocity.field_state(stk::mesh::StateOld), node, dim));
}

// Calls fn with the nodes (dim of them) of each of this processor's owned interface sides, oriented
// out of negativeSideBlockSelector.
template <class FN>
static void for_each_owned_oriented_interface_side(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & interfaceSelector,
    const stk::mesh::Selector & negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FN & fn)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  const unsigned dim = meta.spatial_dimension();
  const stk::topology sideTopology = (3 == dim) ? stk::topology::TRIANGLE_3 : stk::topology::LINE_2;
  for (auto * bucket : mesh.get_buckets(
           meta.side_rank(), interfaceSelector & activePart & meta.locally_owned_part()))
  {
    STK_ThrowRequire(bucket->topology() == sideTopology);
    for (auto side : *bucket)
    {
      std::array<stk::mesh::Entity, 3> sideNodes;
      if (3 == dim)
        sideNodes = get_oriented_triangle_side_nodes(mesh, negativeSideBlockSelector, side);
      else
      {
        const std::array<stk::mesh::Entity, 2> lineNodes =
            get_oriented_line_side_nodes(mesh, negativeSideBlockSelector, side);
        sideNodes = {{lineNodes[0], lineNodes[1], stk::mesh::Entity()}};
      }
      fn(sideNodes);
    }
  }
}

// Appends the node ids of each of this processor's owned interface sides, oriented out of
// negativeSideBlockSelector, to facetVertexIds, and their nodes to usedNodes.
static void append_owned_oriented_interface_side_node_ids(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & interfaceSelector,
    const stk::mesh::Selector & negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    FacetMesh::FacetVertexIds & facetVertexIds,
    std::vector<stk::mesh::Entity> & usedNodes)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  for_each_owned_oriented_interface_side(mesh,
      interfaceSelector,
      negativeSideBlockSelector,
      activePart,
      [&](const std::array<stk::mesh::Entity, 3> & sideNodes)
      {
        std::array<FacetVertexId, 3> facet{{0, 0, 0}};
        for (unsigned n = 0; n < dim; ++n)
        {
          facet[n] = mesh.identifier(sideNodes[n]);
          usedNodes.push_back(sideNodes[n]);
        }
        facetVertexIds.push_back(facet);
      });
}

// The ids and coordinates of the nodes of this processor's owned sides of each surface, oriented
// out of its negativeSideBlockSelector, and each side's node ids.
static void collect_owned_oriented_interface_side_vertices(const stk::mesh::BulkData & mesh,
    const std::vector<std::pair<stk::mesh::Selector, stk::mesh::Selector>> &
        interfaceAndNegativeSideBlockSelectors,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    FacetMesh::VertexIdsAndCoords & vertexIdsAndCoords,
    FacetMesh::FacetVertexIds & facetVertexIds)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::mesh::Entity> usedNodes;
  for (auto & [interfaceSelector, negativeSideBlockSelector] :
      interfaceAndNegativeSideBlockSelectors)
    append_owned_oriented_interface_side_node_ids(
        mesh, interfaceSelector, negativeSideBlockSelector, activePart, facetVertexIds, usedNodes);
  stk::util::sort_and_unique(usedNodes, stk::mesh::EntityLess(mesh));

  vertexIdsAndCoords.reserve(usedNodes.size());
  for (auto node : usedNodes)
    vertexIdsAndCoords.emplace_back(
        mesh.identifier(node), get_vector_field(mesh, coordsField, node, dim));
}

FacetMesh build_facet_mesh_from_interface_sides(const stk::mesh::BulkData & mesh,
    const std::vector<std::pair<stk::mesh::Selector, stk::mesh::Selector>> &
        interfaceAndNegativeSideBlockSelectors,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField)
{
  FacetMesh::VertexIdsAndCoords vertexIdsAndCoords;
  FacetMesh::FacetVertexIds facetVertexIds;
  collect_owned_oriented_interface_side_vertices(mesh,
      interfaceAndNegativeSideBlockSelectors,
      activePart,
      coordsField,
      vertexIdsAndCoords,
      facetVertexIds);

  FacetMesh facetMesh(mesh.parallel(), mesh.mesh_meta_data().spatial_dimension());
  facetMesh.set_owned_facets(vertexIdsAndCoords, facetVertexIds);
  return facetMesh;
}

void set_interface_conforming_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector interfaceSelector,
    const stk::mesh::Selector negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    MeshFacetedSurface & facets)
{
  FacetMesh::VertexIdsAndCoords vertexIdsAndCoords;
  FacetMesh::FacetVertexIds facetVertexIds;
  collect_owned_oriented_interface_side_vertices(mesh,
      {{interfaceSelector, negativeSideBlockSelector}},
      activePart,
      coordsField,
      vertexIdsAndCoords,
      facetVertexIds);
  facets.set_facets(vertexIdsAndCoords, facetVertexIds);
}

std::unique_ptr<MeshFacetedSurface>
build_mesh_faceted_surface_with_velocity_from_interface_sides(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector interfaceSelector,
    const stk::mesh::Selector negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    const FieldRef interfaceVelocity,
    const unsigned numVelocityStates)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::mesh::Entity> usedNodes;
  FacetMesh::FacetVertexIds facetVertexIds;
  append_owned_oriented_interface_side_node_ids(
      mesh, interfaceSelector, negativeSideBlockSelector, activePart, facetVertexIds, usedNodes);
  stk::util::sort_and_unique(usedNodes, stk::mesh::EntityLess(mesh));

  FacetMesh::VertexIdsAndCoords vertexIdsAndCoords;
  std::vector<stk::math::Vector3d> vertexVelocities;
  vertexIdsAndCoords.reserve(usedNodes.size());
  vertexVelocities.reserve(usedNodes.size());
  for (auto node : usedNodes)
  {
    vertexIdsAndCoords.emplace_back(
        mesh.identifier(node), get_vector_field(mesh, coordsField, node, dim));
    vertexVelocities.push_back(
        interface_velocity_at_node(mesh, interfaceVelocity, numVelocityStates, node));
  }
  return std::make_unique<MeshFacetedSurface>(
      mesh.parallel(), dim, vertexIdsAndCoords, facetVertexIds, vertexVelocities);
}

double compute_max_interface_velocity_magnitude(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector interfaceSelector,
    const stk::mesh::Selector negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FieldRef interfaceVelocity,
    const unsigned numVelocityStates)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  double maxSqrMag = 0.;
  for_each_owned_oriented_interface_side(mesh,
      interfaceSelector,
      negativeSideBlockSelector,
      activePart,
      [&](const std::array<stk::mesh::Entity, 3> & sideNodes)
      {
        for (unsigned n = 0; n < dim; ++n)
          maxSqrMag = std::max(maxSqrMag,
              interface_velocity_at_node(mesh, interfaceVelocity, numVelocityStates, sideNodes[n])
                  .length_squared());
      });
  return std::sqrt(maxSqrMag);
}
}

