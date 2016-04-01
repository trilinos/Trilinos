#ifndef MESH_DIAGNOSTICS_HPP
#define MESH_DIAGNOSTICS_HPP

#include <map>
#include "Types.hpp"

namespace stk { namespace mesh {

typedef std::map<stk::mesh::EntityId, std::vector<std::pair<stk::mesh::EntityId, int>>> SplitCoincidentInfo;

class BulkData;

stk::mesh::SplitCoincidentInfo get_split_coincident_elements(stk::mesh::BulkData& bulkData);
std::vector<std::string> get_messages_for_split_coincident_elements(const stk::mesh::BulkData& bulkData, const stk::mesh::SplitCoincidentInfo & splitCoincidentElements);

std::vector<stk::mesh::EntityKeyProc> get_non_unique_key_procs(const stk::mesh::BulkData& bulkData);
std::vector<std::string> get_non_unique_key_messages(const stk::mesh::BulkData& bulkData, const std::vector<stk::mesh::EntityKeyProc> &badKeyProcs);

std::vector<stk::mesh::Entity> get_orphaned_owned_sides(const stk::mesh::BulkData& bulkData);
std::vector<stk::mesh::Entity> get_orphaned_sides_with_attached_element_on_different_proc(const stk::mesh::BulkData& bulkData);
std::vector<std::string> get_messages_for_orphaned_owned_sides(const stk::mesh::BulkData& bulkData, std::vector<stk::mesh::Entity>& keys);

void throw_if_any_proc_has_false(MPI_Comm comm, bool is_all_ok_locally);

} }

#endif
