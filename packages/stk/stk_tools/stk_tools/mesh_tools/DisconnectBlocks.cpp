#include "stk_tools/mesh_tools/DisconnectBlocks.hpp"
#include <stk_util/environment/RuntimeWarning.hpp>
#include "stk_mesh/base/BulkData.hpp"
#include "stk_tools/mesh_tools/CustomAura.hpp"
#include "stk_tools/mesh_tools/DisconnectBlocksImpl.hpp"

namespace stk
{
namespace tools
{

void disconnect_all_blocks(stk::mesh::BulkData & bulk, bool preserveOrphans)
{
  if(bulk.parallel_rank() == 0) {
    std::cout << "Constructing block pairs for disconnect" << std::endl;
  }
  std::vector<impl::BlockPairType> blockPairsToDisconnect = impl::get_block_pairs_to_disconnect(bulk);

  impl::LinkInfo info;
  info.preserveOrphans = preserveOrphans;

  disconnect_block_pairs(bulk, blockPairsToDisconnect, info);
}

void disconnect_block_pairs(stk::mesh::BulkData& bulk, const std::vector<impl::BlockPairType>& blockPairsToDisconnect,
                            impl::LinkInfo& info)
{
  if(blockPairsToDisconnect.empty()) {
    stk::RuntimeWarningP0() << "No block pairs to disconnect" << std::endl;
  }

  bulk.modification_begin();

  if((bulk.parallel_rank() == 0) && (info.debugLevel > 0)) {
    std::cout << "Adding nodes for disconnect" << std::endl;
  }

  for (size_t i = 0; i < blockPairsToDisconnect.size(); ++i) {
    info.os << "P" << bulk.parallel_rank()
         << ": First " << blockPairsToDisconnect[i].first->name() << " second " <<
         blockPairsToDisconnect[i].second->name() << std::endl;
    impl::add_nodes_to_disconnect(bulk, blockPairsToDisconnect[i], info);
  }

  if((bulk.parallel_rank() == 0) && (info.debugLevel > 0)) {
    info.os << "Creating new duplicate node IDs" << std::endl;
  }
  impl::create_new_duplicate_node_IDs(bulk, info);

  if((bulk.parallel_rank() == 0) && (info.debugLevel > 0)) {
    info.os << "Communicating shared node info" << std::endl;
  }
  impl::communicate_shared_node_information(bulk, info);

  if((bulk.parallel_rank() == 0) && (info.debugLevel > 0)) {
    info.os << "Disconnecting elements" << std::endl;
  }
  for (size_t i = 0; i < blockPairsToDisconnect.size(); ++i) {
    impl::disconnect_elements(bulk, blockPairsToDisconnect[i], info);
  }

  info.flush();

  bulk.modification_end();

  if (bulk.has_face_adjacent_element_graph()) {
    bulk.delete_face_adjacent_element_graph();
    bulk.initialize_face_adjacent_element_graph();
  }
}

void reconnect_block_pairs(stk::mesh::BulkData& bulk, const std::vector<impl::BlockPairType>& blockPairsToReconnect,
                            impl::LinkInfo& info)
{
  bulk.modification_begin();

  impl::sanitize_node_map(info.preservedNodeMap, info.os);

  impl::determine_reconnect_node_id(bulk, blockPairsToReconnect, info);

  for(const impl::BlockPairType & blockPair : blockPairsToReconnect) {
    impl::reconnect_block_pair(bulk, blockPair, info);
  }
  impl::fix_indirect_node_sharing(bulk, blockPairsToReconnect, info);

  info.flush();

  bulk.modification_end();
}

}
}
