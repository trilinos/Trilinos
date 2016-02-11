#include "MeshDiagnosticObserver.hpp"
#include <map>
#include "../../base/BulkData.hpp"
#include "MeshDiagnostics.hpp"

namespace stk { namespace mesh {

void MeshDiagnosticObserver::finished_modification_end_notification()
{
    int my_proc_id = mBulkData.parallel_rank();
    std::ofstream out("mesh_diagnostics_failures_" + std::to_string(my_proc_id) + ".txt");

    std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > splitCoincidentElements = get_split_coincident_elements(mBulkData);
    print_and_throw_if_elements_are_split(out, mBulkData, splitCoincidentElements);

    std::vector<stk::mesh::EntityKeyProc> badKeyProcs = stk::mesh::get_non_unique_key_procs(mBulkData);
    stk::mesh::print_and_throw_if_entities_are_not_unique(out, mBulkData, badKeyProcs);

    out.close();
}

}}
