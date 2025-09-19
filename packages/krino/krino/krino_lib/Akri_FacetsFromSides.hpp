#ifndef KRINO_KRINO_KRINO_LIB_AKRI_FACETSFROMSIDES_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_FACETSFROMSIDES_HPP_

#include <stk_mesh/base/Types.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_Surface_Identifier.hpp>

namespace krino {

class FacetedSurfaceBase;

void build_interface_conforming_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector interfaceSelector,
    const stk::mesh::Selector negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    const Surface_Identifier lsIdentifier,
    FacetedSurfaceBase & facets);

void build_interface_conforming_facets_with_interface_velocity(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector interfaceSelector,
    const stk::mesh::Selector negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    const FieldRef interfaceVelocity,
    const unsigned numVelocityStates,
    const Surface_Identifier lsIdentifier,
    FacetedSurfaceBase & facets);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_FACETSFROMSIDES_HPP_ */
