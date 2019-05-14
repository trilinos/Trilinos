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
  std::vector<impl::SideSetType> sideSetsToDisconnect = impl::get_sidesets_to_disconnect(bulk, blockPairsToDisconnect);

  bulk.modification_begin();

  for (size_t i = 0; i < blockPairsToDisconnect.size(); ++i) {
//    std::cout << "Disconnecting block pair: (" << blockPairsToDisconnect[i].first->name() << ", " << blockPairsToDisconnect[i].second->name() << ")" << std::endl;
    impl::NodeMapType nodeMap = impl::get_nodes_to_disconnect(bulk, blockPairsToDisconnect[i], sideSetsToDisconnect[i]);
    impl::create_new_duplicate_nodes(bulk, nodeMap);
    impl::disconnect_face_adjacent_elements(bulk, blockPairsToDisconnect[i], nodeMap);
  }

  impl::NodeMapType nonFaceAdjacentNodeMap = impl::get_non_face_adjacent_nodes(bulk, blockPairsToDisconnect);

  for (size_t i = 0; i < blockPairsToDisconnect.size(); ++i) {
    impl::disconnect_non_face_adjacent_elements(bulk, blockPairsToDisconnect[i], nonFaceAdjacentNodeMap);
  }

  bulk.modification_end();
}

}
}
