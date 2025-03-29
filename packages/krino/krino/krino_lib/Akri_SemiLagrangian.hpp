#ifndef KRINO_KRINO_KRINO_LIB_AKRI_SEMILAGRANGIAN_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_SEMILAGRANGIAN_HPP_

#include <Akri_Faceted_Surface.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_String_Function_Expression.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Types.hpp>

namespace krino {

class Composite_Surface;

using ExtensionVelocityFunction = std::function<stk::math::Vector3d(const double, const stk::math::Vector3d&)>;

ExtensionVelocityFunction build_extension_velocity_at_closest_point_using_string_expressions(const std::vector<String_Function_Expression> & interfaceVelocity);
ExtensionVelocityFunction build_extension_velocity_at_closest_point_using_facets_with_velocity(const int dim, const FacetedSurfaceBase & facets);
ExtensionVelocityFunction build_extension_velocity_using_velocity_at_closest_point(const FacetedSurfaceBase & facets, const ExtensionVelocityFunction & velAtClosestPt);

BoundingBox compute_padded_node_bounding_box_for_semilagrangian_using_string_velocity_expressions(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const std::vector<String_Function_Expression> & interfaceVelocity,
    const FacetedSurfaceBase & facets);

BoundingBox compute_padded_node_bounding_box_for_semilagrangian_using_facets_with_velocity(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double dt,
    const FieldRef coordsField,
    const FacetedSurfaceBase & facets);

void build_nonadaptive_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const FieldRef coordsField,
    const FieldRef distField,
    const double lengthScale,
    FacetedSurfaceBase & facets);

void build_initial_adaptive_facets_after_nodal_distance_is_initialized_from_initial_surfaces(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double time,
    const FieldRef coordsField,
    const FieldRef distField,
    const double avgEdgeLength,
    const Composite_Surface & initSurfaces,
    FacetedSurfaceBase & facets);

void calc_single_step_nonadaptive_semilagrangian_nodal_distance_and_build_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const ExtensionVelocityFunction & extension_velocity,
    const double narrowBandSize,
    const double avgEdgeLength,
    const FacetedSurfaceBase & facetsN,
    FacetedSurfaceBase & facetsNp1);

void calc_single_step_semilagrangian_nodal_distance_and_build_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const ExtensionVelocityFunction & extension_velocity,
    const double narrowBandSize,
    const double avgEdgeLength,
    const FacetedSurfaceBase & facetsN,
    FacetedSurfaceBase & facetsNp1);

void predict_semilagrangian_nodal_distance_and_build_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const ExtensionVelocityFunction & extension_velocity,
    const double narrowBandSize,
    const double avgEdgeLength,
    const FacetedSurfaceBase & facetsN,
    FacetedSurfaceBase & facetsPred);

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
    const FacetedSurfaceBase & facetsN,
    FacetedSurfaceBase & facetsNp1);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_SEMILAGRANGIAN_HPP_ */
