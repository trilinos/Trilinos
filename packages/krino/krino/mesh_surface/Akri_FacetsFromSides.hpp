#ifndef KRINO_KRINO_MESH_SURFACE_AKRI_FACETSFROMSIDES_HPP_
#define KRINO_KRINO_MESH_SURFACE_AKRI_FACETSFROMSIDES_HPP_

#include <stk_mesh/base/Types.hpp>
#include <Akri_FieldRef.hpp>

namespace krino {

class FacetedSurfaceBase;

void append_interface_conforming_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector interfaceSelector,
    const stk::mesh::Selector negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    FacetedSurfaceBase & facets);

void build_interface_conforming_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector interfaceSelector,
    const stk::mesh::Selector negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    FacetedSurfaceBase & facets);

void build_interface_conforming_facets_with_interface_velocity(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector interfaceSelector,
    const stk::mesh::Selector negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    const FieldRef interfaceVelocity,
    const unsigned numVelocityStates,
    FacetedSurfaceBase & facets);

}

#endif /* KRINO_KRINO_MESH_SURFACE_AKRI_FACETSFROMSIDES_HPP_ */
