#ifndef KRINO_KRINO_MESH_UTILS_AKRI_SMALLCLUSTERCLEANUP_HPP_
#define KRINO_KRINO_MESH_UTILS_AKRI_SMALLCLUSTERCLEANUP_HPP_
#include <tuple>
#include <vector>

#include <stk_mesh/base/Types.hpp>

namespace krino {

class Phase_Support;
class NamedPhase;

void cleanup_small_clusters_of_phase(stk::mesh::BulkData & mesh, const Phase_Support & phaseSupport, const NamedPhase & fromPhase, const NamedPhase & toPhase, const size_t smallClusterSize);
void cleanup_small_clusters_of_phases(stk::mesh::BulkData & mesh, const Phase_Support & phaseSupport, const std::vector<std::tuple<NamedPhase,NamedPhase,unsigned>> & fromPhaseToPhaseAndClusterSize);

}

#endif /* KRINO_KRINO_MESH_UTILS_AKRI_SMALLCLUSTERCLEANUP_HPP_ */
