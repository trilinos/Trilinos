#ifndef KRINO_KRINO_KRINO_LIB_AKRI_NODALGRADIENT_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_NODALGRADIENT_HPP_
#include <stk_mesh/base/Types.hpp>
#include <string>

namespace krino {

class FieldRef;

std::string get_nodal_gradient_field_name(const std::string & scalarFieldName);
FieldRef get_nodal_gradient_for_scalar_field(const stk::mesh::MetaData & meta, const FieldRef scalarField);

FieldRef register_nodal_gradient_for_scalar_field(stk::mesh::MetaData & meta, FieldRef scalarField);

void update_nodal_gradient(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef scalarField);
void update_nodal_gradient(const stk::mesh::BulkData & mesh, const FieldRef coords, const FieldRef scalarField, const FieldRef gradientField);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_NODALGRADIENT_HPP_ */
