#ifndef MESH_DIAGNOSTICS_HPP
#define MESH_DIAGNOSTICS_HPP

#include <map>
#include "../../base/Types.hpp"

namespace stk { namespace mesh {

class BulkData;

std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > get_split_coincident_elements(stk::mesh::BulkData& bulkData);
std::string get_message_for_split_coincident_elements(const stk::mesh::BulkData& bulkData, const std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > & splitCoincidentElements);

std::vector<stk::mesh::EntityKeyProc> get_non_unique_key_procs(const stk::mesh::BulkData& bulkData);
std::string get_non_unique_key_messages(const stk::mesh::BulkData& bulkData, const std::vector<stk::mesh::EntityKeyProc> &badKeyProcs);

void throw_if_any_proc_has_false(MPI_Comm comm, bool is_all_ok_locally);

} }

#endif
