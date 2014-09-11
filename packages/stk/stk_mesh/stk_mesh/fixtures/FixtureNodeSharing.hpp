#ifndef Stk_Mesh_Fixtures_FixtureNodeSharing_hpp
#define Stk_Mesh_Fixtures_FixtureNodeSharing_hpp

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/Types.hpp>      // for EntityId, EntityRank
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine


namespace stk {
namespace mesh {
namespace fixtures {

typedef std::multimap<EntityId, int> NodeToProcsMMap;

void AddToNodeProcsMMap(NodeToProcsMMap &nodes_to_procs, EntityId node_id, int proc_rank);

void DoAddNodeSharings(BulkData &bulk_data, NodeToProcsMMap &nodes_to_procs, EntityId node_id, Entity node);

} // fixtures
} // mesh
} // stk
#endif
