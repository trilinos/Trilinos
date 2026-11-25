#ifndef KRINO_KRINO_LEVELSET_SHAPE_SENS_AKRI_INTERFACEINVERSEGRADIENTAREA_HPP_
#define KRINO_KRINO_LEVELSET_SHAPE_SENS_AKRI_INTERFACEINVERSEGRADIENTAREA_HPP_

#include <Akri_LS_Field.hpp>
#include <stk_mesh/base/Types.hpp>

namespace krino
{

class FieldRef;

// The interface inverse gradient area of a level set: at each node of its negative interface sides,
// sum_s A_s / |grad phi_s| over the node's sides s, with A_s the side's length (2D) or area (3D)
// and grad phi_s the level set gradient of the side's negative-phase parent element.  Each term is
// the side's normal speed for a uniform level set increase times its area, its contribution to the
// enclosed volume's derivative.  Using the negative-phase element makes this the derivative for
// increasing level set where the gradient jumps across a side on nodes with zero level set.  The
// areas use the current (for snapping, snapped) coordinates; the gradients use the original level
// set at the unsnapped parent nodes.

// Registers the field on the nodes of the level set's negative interface sides.
FieldRef register_interface_inverse_gradient_area_field(
    stk::mesh::MetaData & meta, const LS_Field & lsField);

FieldRef get_interface_inverse_gradient_area_field(
    const stk::mesh::MetaData & meta, const LS_Field & lsField);

// Collective; needs no aura.  Linear simplex parent elements only.  The sides are the level set's
// interface sides, or those of interfaceSideSelector, a subset of them.
FieldRef
update_interface_inverse_gradient_area(const stk::mesh::BulkData & mesh, const LS_Field & lsField);
FieldRef update_interface_inverse_gradient_area(const stk::mesh::BulkData & mesh,
    const LS_Field & lsField,
    const stk::mesh::Selector & interfaceSideSelector);

} // namespace krino

#endif /* KRINO_KRINO_LEVELSET_SHAPE_SENS_AKRI_INTERFACEINVERSEGRADIENTAREA_HPP_ */
