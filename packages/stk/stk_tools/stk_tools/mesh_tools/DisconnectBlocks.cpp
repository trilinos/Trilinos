#include "stk_tools/mesh_tools/DisconnectBlocks.hpp"
#include "stk_tools/mesh_tools/DisconnectBlocksImpl.hpp"
#include "stk_mesh/base/BulkData.hpp"

namespace stk
{
namespace tools
{

void disconnect_all_blocks(stk::mesh::BulkData & bulk)
{
  std::vector<impl::BlockPairType> blockPairsToDisconnect = impl::get_block_pairs_to_disconnect(bulk);

  bulk.modification_begin();

  impl::NodeMapType nodeMap;
  for (size_t i = 0; i < blockPairsToDisconnect.size(); ++i) {
    impl::add_nodes_to_disconnect(bulk, blockPairsToDisconnect[i], nodeMap);
  }

  impl::create_new_duplicate_node_IDs(bulk, nodeMap);

  impl::communicate_shared_node_information(bulk, nodeMap);

  for (size_t i = 0; i < blockPairsToDisconnect.size(); ++i) {
    impl::disconnect_elements(bulk, blockPairsToDisconnect[i], nodeMap);
  }

  bulk.modification_end();

  // Trigger a re-generation of the element graph because the relation changes made here
  // do not trigger any of our notifiers to update the graph.
  if (bulk.has_face_adjacent_element_graph()) {
    bulk.delete_face_adjacent_element_graph();
    bulk.initialize_face_adjacent_element_graph();
  }
}

}
}
