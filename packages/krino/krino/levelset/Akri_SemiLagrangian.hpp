#ifndef KRINO_KRINO_LEVELSET_AKRI_SEMILAGRANGIAN_HPP_
#define KRINO_KRINO_LEVELSET_AKRI_SEMILAGRANGIAN_HPP_

#include <Akri_BoundingBox.hpp>
#include <Akri_FacetAppender.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_String_Function_Expression.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Types.hpp>

namespace krino {

class CSG_Surface;
class MeshFacetedSurface;

using ExtensionVelocityFunction = std::function<stk::math::Vector3d(const double, const stk::math::Vector3d&)>;

ExtensionVelocityFunction build_extension_velocity_at_closest_point_using_string_expressions(const std::vector<String_Function_Expression> & interfaceVelocity);
ExtensionVelocityFunction build_extension_velocity_at_closest_point_using_facets_with_velocity(
    const MeshFacetedSurface & facets);
ExtensionVelocityFunction build_extension_velocity_using_velocity_at_closest_point(
    const MeshFacetedSurface & facets, const ExtensionVelocityFunction & velAtClosestPt);

BoundingBox compute_padded_node_bounding_box_for_semilagrangian_using_string_velocity_expressions(
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const std::vector<String_Function_Expression> & interfaceVelocity,
    const MeshFacetedSurface & facets);

// The bounding box of the nodes, padded for the distance travelled in dt at maxVelocityMagnitude.
BoundingBox compute_padded_node_bounding_box_for_semilagrangian(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double dt,
    const FieldRef coordsField,
    const double maxVelocityMagnitude);

void build_nonadaptive_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const FieldRef coordsField,
    const FieldRef distField,
    const double lengthScale,
    FacetAppender & facets);

void build_initial_adaptive_facets_after_nodal_distance_is_initialized_from_initial_surfaces(
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double time,
    const FieldRef coordsField,
    const FieldRef distField,
    const double avgEdgeLength,
    const CSG_Surface & initSurfaces,
    FacetAppender & facets);

void calc_single_step_nonadaptive_semilagrangian_nodal_distance_and_build_facets(
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const ExtensionVelocityFunction & extension_velocity,
    const double narrowBandSize,
    const double avgEdgeLength,
    const MeshFacetedSurface & facetsN,
    FacetAppender & facetsNp1);

void calc_single_step_semilagrangian_nodal_distance_and_build_facets(
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const ExtensionVelocityFunction & extension_velocity,
    const double narrowBandSize,
    const double avgEdgeLength,
    const MeshFacetedSurface & facetsN,
    FacetAppender & facetsNp1);

void predict_semilagrangian_nodal_distance_and_build_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const ExtensionVelocityFunction & extension_velocity,
    const double narrowBandSize,
    const double avgEdgeLength,
    const MeshFacetedSurface & facetsN,
    FacetAppender & facetsPred);

void correct_semilagrangian_nodal_distance_and_build_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const ExtensionVelocityFunction & extension_velocity_old,
    const ExtensionVelocityFunction & extension_velocity_pred,
    const double narrowBandSize,
    const double avgEdgeLength,
    const MeshFacetedSurface & facetsN,
    FacetAppender & facetsNp1);
}

#endif /* KRINO_KRINO_LEVELSET_AKRI_SEMILAGRANGIAN_HPP_ */
