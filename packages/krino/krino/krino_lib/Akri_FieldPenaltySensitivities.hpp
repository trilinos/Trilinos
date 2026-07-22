#include "stk_mesh/base/BulkData.hpp"

namespace krino
{
double compute_l2_penalty_between_fields(const stk::mesh::BulkData & mesh, const stk::mesh::FieldBase & fldWSens, 
  const stk::mesh::FieldBase & otherFld, double normalization, const stk::mesh::Selector & s, std::map<stk::mesh::EntityId, double> & sens);

double compute_gradient_penalty_between_fields(const stk::mesh::BulkData & mesh, const stk::mesh::FieldBase & fldWSens, 
  const stk::mesh::FieldBase & otherFld, const stk::mesh::Selector & s, std::map<stk::mesh::EntityId, double> & sens);
}