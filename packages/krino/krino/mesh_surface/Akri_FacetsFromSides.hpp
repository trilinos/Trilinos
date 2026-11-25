#ifndef KRINO_KRINO_MESH_SURFACE_AKRI_FACETSFROMSIDES_HPP_
#define KRINO_KRINO_MESH_SURFACE_AKRI_FACETSFROMSIDES_HPP_

#include <memory>
#include <utility>
#include <vector>

#include <stk_mesh/base/Types.hpp>
#include <Akri_FacetMesh.hpp>
#include <Akri_FieldRef.hpp>

namespace krino {

class MeshFacetedSurface;

// Replaces the facets of the surface with this processor's owned interface sides (interfaceSelector
// & activePart, oriented out of negativeSideBlockSelector), connected: their vertices are their
// nodes, with the node ids.
void set_interface_conforming_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector interfaceSelector,
    const stk::mesh::Selector negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    MeshFacetedSurface & facets);

// This processor's owned sides of each surface (interfaceSelector & activePart, oriented out of its
// negativeSideBlockSelector) as a connected FacetMesh whose vertices are their nodes.  A side of
// two surfaces is two facets.
FacetMesh build_facet_mesh_from_interface_sides(const stk::mesh::BulkData & mesh,
    const std::vector<std::pair<stk::mesh::Selector, stk::mesh::Selector>> &
        interfaceAndNegativeSideBlockSelectors,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField);

// This processor's owned interface sides (interfaceSelector & activePart, oriented out of
// negativeSideBlockSelector) as a connected MeshFacetedSurface whose vertices are their nodes, each
// with the interface velocity: NP1, or with two states the average of N and NP1.
std::unique_ptr<MeshFacetedSurface>
build_mesh_faceted_surface_with_velocity_from_interface_sides(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector interfaceSelector,
    const stk::mesh::Selector negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    const FieldRef interfaceVelocity,
    const unsigned numVelocityStates);

// The largest magnitude of the interface velocity, as above, at the nodes of this processor's owned
// interface sides.
double compute_max_interface_velocity_magnitude(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector interfaceSelector,
    const stk::mesh::Selector negativeSideBlockSelector,
    const stk::mesh::Part & activePart,
    const FieldRef interfaceVelocity,
    const unsigned numVelocityStates);
}

#endif /* KRINO_KRINO_MESH_SURFACE_AKRI_FACETSFROMSIDES_HPP_ */
