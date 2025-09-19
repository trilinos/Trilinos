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
#include <Akri_Faceted_Surface.hpp>
#include <Akri_Surface_Identifier.hpp>
#include <Akri_OrientedSideNodes.hpp>

namespace krino {

static std::array<stk::math::Vector3d,2> get_line_side_vector(const stk::mesh::BulkData & mesh, const FieldRef vecField, const std::array<stk::mesh::Entity,2> lineNodes)
{
  return {{ get_vector_field(mesh, vecField, lineNodes[0], 2), get_vector_field(mesh, vecField, lineNodes[1], 2) }};
}

static void append_facet_from_triangle_side(const stk::mesh::BulkData & mesh, const FieldRef coords, const stk::mesh::Selector & /*interfaceSelector*/, const stk::mesh::Selector & negativeSideElementSelector, const stk::mesh::Entity side, std::vector<Facet3d> & facets)
{
  const std::array<stk::mesh::Entity,3> orientedSideNodes = get_oriented_triangle_side_nodes(mesh, negativeSideElementSelector, side);
  const std::array<stk::math::Vector3d,3> sideNodeCoords = get_triangle_vector(mesh, coords, orientedSideNodes);
  facets.emplace_back( sideNodeCoords[0], sideNodeCoords[1], sideNodeCoords[2] );
}

static void append_facet_from_line_side(const stk::mesh::BulkData & mesh, const FieldRef coords, const stk::mesh::Selector & /*sideSelector*/, const stk::mesh::Selector & negativeSideElementSelector, const stk::mesh::Entity side, std::vector<Facet2d> & facets)
{
  const std::array<stk::mesh::Entity,2> orientedSideNodes = get_oriented_line_side_nodes(mesh, negativeSideElementSelector, side);
  const std::array<stk::math::Vector3d,2> sideNodeCoords = get_line_side_vector(mesh, coords, orientedSideNodes);
  facets.emplace_back(sideNodeCoords[0], sideNodeCoords[1]);
}

static void append_facet_with_velocity_from_triangle_side(const stk::mesh::BulkData & mesh, const FieldRef coords, const FieldRef interfaceVelocity, const unsigned /*numVelocityStates*/, const stk::mesh::Selector & /*interfaceSelector*/, const stk::mesh::Selector & negativeSideElementSelector, const stk::mesh::Entity side, std::vector<FacetWithVelocity3d> & facets)
{
  const std::array<stk::mesh::Entity,3> orientedSideNodes = get_oriented_triangle_side_nodes(mesh, negativeSideElementSelector, side);
  const std::array<stk::math::Vector3d,3> sideNodeCoords = get_triangle_vector(mesh, coords, orientedSideNodes);
  const std::array<stk::math::Vector3d,3> sideNodeVelocity = get_triangle_vector(mesh, interfaceVelocity, orientedSideNodes);
  facets.emplace_back( sideNodeCoords[0], sideNodeCoords[1], sideNodeCoords[2], sideNodeVelocity[0], sideNodeVelocity[1], sideNodeVelocity[2] );
}

static void append_facet_with_velocity_from_line_side(const stk::mesh::BulkData & mesh, const FieldRef coords, const FieldRef interfaceVelocity, const unsigned numVelocityStates, const stk::mesh::Selector & /*sideSelector*/, const stk::mesh::Selector & negativeSideElementSelector, const stk::mesh::Entity side, std::vector<FacetWithVelocity2d> & facets)
{
  const std::array<stk::mesh::Entity,2> orientedSideNodes = get_oriented_line_side_nodes(mesh, negativeSideElementSelector, side);
  const std::array<stk::math::Vector3d,2> sideNodeCoords = get_line_side_vector(mesh, coords, orientedSideNodes);
  const std::array<stk::math::Vector3d,2> sideNodeVelocity = get_line_side_vector(mesh, interfaceVelocity, orientedSideNodes);
  if (1 == numVelocityStates)
  {
    facets.emplace_back( sideNodeCoords[0], sideNodeCoords[1], sideNodeVelocity[0], sideNodeVelocity[1] );
  }
  else
  {
    STK_ThrowAssert(2==numVelocityStates && interfaceVelocity.number_of_states() > 1);
    const std::array<stk::math::Vector3d,2> sideNodeVelocityOld = get_line_side_vector(mesh, interfaceVelocity.field_state(stk::mesh::StateOld), orientedSideNodes);
    facets.emplace_back( sideNodeCoords[0], sideNodeCoords[1], 0.5*(sideNodeVelocity[0]+sideNodeVelocityOld[0]), 0.5*(sideNodeVelocity[1]+sideNodeVelocityOld[1]) );
  }
}

static void append_owned_facets_from_triangle_sides(const stk::mesh::BulkData & mesh,
    const FieldRef cooordsField,
    const stk::mesh::Selector & sideSelector,
    const stk::mesh::Selector & negativeSideElementSelector,
    std::vector<Facet3d> & facets)
{
  for ( auto & bucket : mesh.get_buckets( mesh.mesh_meta_data().side_rank(), sideSelector & mesh.mesh_meta_data().locally_owned_part()) )
  {
    STK_ThrowRequire(bucket->topology() == stk::topology::TRIANGLE_3);
    for (auto & side : *bucket)
      append_facet_from_triangle_side(mesh, cooordsField, sideSelector, negativeSideElementSelector, side, facets);
  }
}

static void append_owned_facets_from_line_sides(const stk::mesh::BulkData & mesh,
    const FieldRef cooordsField,
    const stk::mesh::Selector & sideSelector,
    const stk::mesh::Selector & negativeSideElementSelector,
    std::vector<Facet2d> & facets)
{
  for ( auto & bucket : mesh.get_buckets( mesh.mesh_meta_data().side_rank(), sideSelector & mesh.mesh_meta_data().locally_owned_part()) )
  {
    STK_ThrowRequire(bucket->topology() == stk::topology::LINE_2);
    for (auto & side : *bucket)
      append_facet_from_line_side(mesh, cooordsField, sideSelector, negativeSideElementSelector, side, facets);
  }
}

static void append_owned_facets_with_velocity_from_triangle_sides(const stk::mesh::BulkData & mesh,
    const FieldRef cooordsField,
    const FieldRef interfaceVelocity,
    const unsigned numVelocityStates,
    const stk::mesh::Selector & sideSelector,
    const stk::mesh::Selector & negativeSideElementSelector,
    std::vector<FacetWithVelocity3d> & facets)
{
  for ( auto & bucket : mesh.get_buckets( mesh.mesh_meta_data().side_rank(), sideSelector & mesh.mesh_meta_data().locally_owned_part()) )
  {
    STK_ThrowRequire(bucket->topology() == stk::topology::TRIANGLE_3);
    for (auto & side : *bucket)
      append_facet_with_velocity_from_triangle_side(mesh, cooordsField, interfaceVelocity, numVelocityStates, sideSelector, negativeSideElementSelector, side, facets);
  }
}

static void append_owned_facets_with_velocity_from_line_sides(const stk::mesh::BulkData & mesh,
    const FieldRef cooordsField,
    const FieldRef interfaceVelocity,
    const unsigned numVelocityStates,
    const stk::mesh::Selector & sideSelector,
    const stk::mesh::Selector & negativeSideElementSelector,
    std::vector<FacetWithVelocity2d> & facets)
{
  for ( auto & bucket : mesh.get_buckets( mesh.mesh_meta_data().side_rank(), sideSelector & mesh.mesh_meta_data().locally_owned_part()) )
  {
    STK_ThrowRequire(bucket->topology() == stk::topology::LINE_2);
    for (auto & side : *bucket)
      append_facet_with_velocity_from_line_side(mesh, cooordsField, interfaceVelocity, numVelocityStates, sideSelector, negativeSideElementSelector, side, facets);
  }
}

void build_interface_conforming_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector interfaceSelector,
    const stk::mesh::Selector negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    const Surface_Identifier lsIdentifier,
    FacetedSurfaceBase & facets)
{
  const stk::mesh::Selector sideSelector = interfaceSelector & activePart;
  const stk::mesh::Selector ownedSideSelector = sideSelector & mesh.mesh_meta_data().locally_owned_part();

  facets.clear();
  if (3 == mesh.mesh_meta_data().spatial_dimension())
    append_owned_facets_from_triangle_sides(mesh, coordsField, sideSelector, negativeSideBlockSelector, facets.as_derived_type<Facet3d>().get_facets());
  else
    append_owned_facets_from_line_sides(mesh, coordsField, sideSelector, negativeSideBlockSelector, facets.as_derived_type<Facet2d>().get_facets());
}

void build_interface_conforming_facets_with_interface_velocity(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector interfaceSelector,
    const stk::mesh::Selector negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    const FieldRef interfaceVelocity,
    const unsigned numVelocityStates,
    const Surface_Identifier lsIdentifier,
    FacetedSurfaceBase & facets)
{
  const stk::mesh::Selector sideSelector = interfaceSelector & activePart;
  const stk::mesh::Selector ownedSideSelector = sideSelector & mesh.mesh_meta_data().locally_owned_part();

  facets.clear();
  if (3 == mesh.mesh_meta_data().spatial_dimension())
    append_owned_facets_with_velocity_from_triangle_sides(mesh, coordsField, interfaceVelocity, numVelocityStates, sideSelector, negativeSideBlockSelector, facets.as_derived_type<FacetWithVelocity3d>().get_facets());
  else
    append_owned_facets_with_velocity_from_line_sides(mesh, coordsField, interfaceVelocity, numVelocityStates, sideSelector, negativeSideBlockSelector, facets.as_derived_type<FacetWithVelocity2d>().get_facets());
}

}

