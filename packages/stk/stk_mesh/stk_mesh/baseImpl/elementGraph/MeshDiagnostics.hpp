#ifndef MESH_DIAGNOSTICS_HPP
#define MESH_DIAGNOSTICS_HPP

#include <map>
#include "../../base/Types.hpp"

namespace stk { namespace mesh {

class BulkData;

std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > get_split_coincident_elements(stk::mesh::BulkData& bulkData);
void print_and_throw_if_elements_are_split(const stk::mesh::BulkData& bulkData, const std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > &splitCoincidentElements);

} }

#endif
