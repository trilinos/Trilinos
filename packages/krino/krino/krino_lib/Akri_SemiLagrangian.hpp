#ifndef KRINO_KRINO_KRINO_LIB_AKRI_SEMILAGRANGIAN_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_SEMILAGRANGIAN_HPP_

#include <Akri_Faceted_Surface.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_String_Function_Expression.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Types.hpp>

namespace krino {

class Composite_Surface;

BoundingBox compute_padded_node_bounding_box_for_semilagrangian(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
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
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
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
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
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
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
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
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
    const double narrowBandSize,
    const double avgEdgeLength,
    const FacetedSurfaceBase & facetsN,
    const FacetedSurfaceBase & facetsPred,
    FacetedSurfaceBase & facetsNp1);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_SEMILAGRANGIAN_HPP_ */
