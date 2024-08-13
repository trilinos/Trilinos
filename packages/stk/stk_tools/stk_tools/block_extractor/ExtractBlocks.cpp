#include "mpi.h"
#include <stk_util/stk_config.h>
#include <stk_tools/block_extractor/ExtractBlocks.hpp>
#include "stk_tools/transfer_utils/TransientFieldTransferById.hpp"
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_tools/mesh_clone/MeshClone.hpp>
#include "stk_mesh/base/Selector.hpp"
#include <stk_io/StkMeshIoBroker.hpp>
#include "stk_io/IossBridge.hpp"
#include "stk_util/parallel/ParallelReduceBool.hpp"

namespace stk {
namespace tools {

bool part_selects_entities_globally(stk::mesh::BulkData & bulk, stk::mesh::Part & part)
{
  stk::mesh::BucketVector const & buckets = bulk.get_buckets(part.primary_entity_rank(), part);
  bool hasEntities = !buckets.empty();
  return stk::is_true_on_any_proc(bulk.parallel(), hasEntities);
}

void remove_io_attribute_from_empty_parts(stk::mesh::BulkData & bulk)
{
  stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  for (stk::mesh::Part * part : meta.get_mesh_parts()) {
    if (stk::io::has_io_part_attribute(*part)) {
      if (!part_selects_entities_globally(bulk, *part)) {
        stk::io::remove_io_part_attribute(*part);
      }
    }
  }
}

std::vector<std::string> GetEntityNamesFromIDs(const stk::mesh::BulkData & meshBulk, stk::topology::rank_t  stkTop, const std::vector<int> & theIDs)
{

  std::vector<std::string> entity_names;
  const stk::mesh::PartVector &parts = meshBulk.mesh_meta_data().get_mesh_parts();
  for (size_t IdIndex = 0; IdIndex < theIDs.size(); IdIndex++)
  {
    for(size_t i = 0; i < parts.size(); i++)
    {

      if ( parts[i]->primary_entity_rank() == stkTop && parts[i]->id() == theIDs[IdIndex] )
      {
        entity_names.push_back(parts[i]->name());
      }
    }
  }
  return entity_names;
}

std::vector<std::string> GetBlockNamesFromIDs(const stk::mesh::BulkData & meshBulk, const std::vector<int> & block_ids)
{
  return (GetEntityNamesFromIDs( meshBulk, stk::topology::ELEMENT_RANK, block_ids));
}

std::vector<std::string> find_nodeset_names_from_id(const stk::mesh::BulkData & meshBulk, const std::vector<int> & nodeset_ids)
{
  return (GetEntityNamesFromIDs( meshBulk, stk::topology::NODE_RANK, nodeset_ids));
}

void GetPartsByName(std::vector<stk::mesh::Part*> & parts,
                    const stk::mesh::BulkData& inBulk,
                    std::vector < std::string > names)
{
  for(size_t i = 0; i < names.size(); i++)
  {
    parts.push_back(inBulk.mesh_meta_data().get_part(names[i]));
    STK_ThrowRequireMsg(parts[i] != nullptr, "Can't find " <<  names[i] << " in mesh.\n");
  }
}

stk::mesh::Selector GetBlockAndNodesetSelector(const stk::mesh::BulkData & inBulk,
                                               const std::vector<std::string>& nodesetNames,
                                               const std::vector<std::string>& blockNames)
{

  stk::mesh::PartVector parts;
  GetPartsByName(parts, inBulk, nodesetNames);
  GetPartsByName(parts, inBulk, blockNames);
  stk::mesh::Selector nodeset_and_block_selector = stk::mesh::selectUnion(parts);

  return nodeset_and_block_selector;
}


void extract_blocks_and_ns_from_file(const std::string &inFile,
                                     const std::string &outFile,
                                     const std::vector<int> &blockIDs,
                                     const std::vector<int> &nodesetIDs,
                                     MPI_Comm comm)
{
  stk::mesh::MeshBuilder builder(comm);
  builder.set_aura_option(stk::mesh::BulkData::AUTO_AURA);
  std::shared_ptr<stk::mesh::BulkData> inBulk = builder.create();
  std::shared_ptr<stk::mesh::BulkData> outBulk = builder.create();

  stk::io::StkMeshIoBroker stkInput;
  stk::io::fill_mesh_preexisting(stkInput, inFile, *inBulk);


  stk::mesh::Selector nothingSelector_byDefaultConstruction;
  stk::mesh::Selector allSelector(!nothingSelector_byDefaultConstruction);

  std::vector<std::string> blockNames = GetBlockNamesFromIDs(*inBulk, blockIDs);

  stk::mesh::Selector nodeset_and_block_selector = allSelector;
  // if user asks for nodeset, use selector to output nodesets and blocks (this ensures we actually get the nodes associated with the nodesets
  if (nodesetIDs.size() > 0)
  {
    std::vector < std::string > nodesetNames = find_nodeset_names_from_id(*inBulk, nodesetIDs);
    nodeset_and_block_selector = GetBlockAndNodesetSelector(*inBulk, nodesetNames, blockNames);
    stk::tools::copy_mesh(*inBulk, allSelector, *outBulk);
  }
  else // if user only asks for blocks, use extract_blocks, which completely deletes all other blocks and associated nodes
  {
    extract_blocks(*inBulk, *outBulk, blockNames);
    remove_io_attribute_from_empty_parts(*outBulk);
  }

  stk::io::StkMeshIoBroker stkOutput;
  stkOutput.set_bulk_data(*outBulk);
  stkOutput.set_attribute_field_ordering_stored_by_part_ordinal(stkInput.get_attribute_field_ordering_stored_by_part_ordinal());

  stk::transfer_utils::TransientFieldTransferById transfer(stkInput, stkOutput);
  transfer.transfer_and_write_transient_fields(outFile,nodeset_and_block_selector);
}



void extract_blocks(stk::mesh::BulkData &oldBulk, stk::mesh::BulkData &newBulk, const std::vector<std::string> &blockNames)
{
  stk::mesh::PartVector parts(blockNames.size());
  for(size_t i=0; i<blockNames.size(); i++)
  {
    parts[i] = oldBulk.mesh_meta_data().get_part(blockNames[i]);
    STK_ThrowRequireMsg(parts[i] != nullptr, "Can't find block " << blockNames[i] << " in mesh.\n");
  }

  stk::mesh::Selector selector = stk::mesh::selectUnion(parts);
  stk::tools::copy_mesh(oldBulk, selector, newBulk);
}

}

}
