#ifndef KRINO_KRINO_FIELD_CALCS_AKRI_SURFACENODALNORMAL_HPP_
#define KRINO_KRINO_FIELD_CALCS_AKRI_SURFACENODALNORMAL_HPP_
#include <stk_mesh/base/Types.hpp>
#include <string>

namespace krino
{

class FieldRef;

// Area-weighted nodal normal of a mesh surface: the sum over the node's sides of the side's length
// (2D) or area (3D) times its unit normal.  At a node of a closed surface this is 2 times (2D) or
// 3 times (3D) the derivative of the enclosed volume with respect to the node's coordinates.
// The field is allocated only on the nodes of surfaceSelector, which must be parts that induce
// membership onto their nodes.
FieldRef register_surface_nodal_normal_field(stk::mesh::MetaData & meta,
    const std::string & name,
    const stk::mesh::Selector & surfaceSelector);

// Each side in sideSelector adds its area vector, oriented out of each attached element in
// elementSelector, to its nodes.  Nodes with no such side get zero.  Not normalized.  Linear line
// or tri sides of tri or tet elements only.  Collective; needs no aura.
void update_surface_nodal_normal(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef normalField,
    const stk::mesh::Selector & sideSelector,
    const stk::mesh::Selector & elementSelector);
} // namespace krino

#endif /* KRINO_KRINO_FIELD_CALCS_AKRI_SURFACENODALNORMAL_HPP_ */
