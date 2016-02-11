#include "MeshDiagnosticObserver.hpp"
#include <map>
#include "../../base/BulkData.hpp"
#include "MeshDiagnostics.hpp"

namespace stk { namespace mesh {

void MeshDiagnosticObserver::finished_modification_end_notification()
{
    std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > splitCoincidentElements = get_split_coincident_elements(mBulkData);
    print_and_throw_if_elements_are_split(mBulkData, splitCoincidentElements);
}

}}
