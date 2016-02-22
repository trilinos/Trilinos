#include "MeshDiagnosticObserver.hpp"
#include <map>
#include "../../base/BulkData.hpp"
#include "MeshDiagnostics.hpp"

namespace stk { namespace mesh {

void MeshDiagnosticObserver::finished_modification_end_notification()
{
    int my_proc_id = mBulkData.parallel_rank();
    std::ofstream out("mesh_diagnostics_failures_" + std::to_string(my_proc_id) + ".txt", std::ios::app);

    std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > splitCoincidentElements = stk::mesh::get_split_coincident_elements(mBulkData);
    out << stk::mesh::get_message_for_split_coincident_elements(mBulkData, splitCoincidentElements);

    std::vector<stk::mesh::EntityKeyProc> badKeyProcs = stk::mesh::get_non_unique_key_procs(mBulkData);
    out << stk::mesh::get_non_unique_key_messages(mBulkData, badKeyProcs);

    std::vector<stk::mesh::EntityKey> orphanedSides = stk::mesh::get_orphaned_owned_sides(mBulkData);
    out << stk::mesh::get_messages_for_orphaned_owned_nodes(mBulkData, orphanedSides);

    out.close();

    //throw_if_any_proc_has_false(mBulkData.parallel(), badKeyProcs.empty() && splitCoincidentElements.empty() && orphanedSides.empty());
}

}}
