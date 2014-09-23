/*
 * FixtureNodeSharing.cpp
 *
 *  Created on: Sep 12, 2014
 *      Author: pgxavie
 */


#include "stk_mesh/fixtures/FixtureNodeSharing.hpp"


namespace stk {
namespace mesh {
namespace fixtures {

void AddToNodeProcsMMap(NodeToProcsMMap &nodes_to_procs, EntityId node_id, int proc_rank)
{
  bool need_add = true;
  NodeToProcsMMap::iterator np_j = nodes_to_procs.find(node_id);
  NodeToProcsMMap::iterator np_e = nodes_to_procs.end();
  for (; (np_j != np_e) && (np_j->first == node_id); ++np_j ) {
    if (np_j->second == proc_rank) {
      need_add = false;
      break;
    }
  }
  if (need_add) {
    nodes_to_procs.insert(std::pair<EntityId, int>(node_id, proc_rank));
  }
}

void DoAddNodeSharings(BulkData &bulk_data, NodeToProcsMMap &nodes_to_procs, EntityId node_id, Entity node)
{
  const int p_rank = bulk_data.parallel_rank();
  NodeToProcsMMap::iterator np_j = nodes_to_procs.find(node_id);
  NodeToProcsMMap::iterator np_e = nodes_to_procs.end();
  for (; (np_j != np_e) && (np_j->first == node_id); ++np_j ) {
    int other_proc = np_j->second;
    if (p_rank != other_proc) {
      bulk_data.add_node_sharing(node, other_proc);
    }
  }
}

}}}
