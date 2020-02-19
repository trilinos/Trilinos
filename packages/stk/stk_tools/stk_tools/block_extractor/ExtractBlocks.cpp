#include "mpi.h"
#include <stk_util/stk_config.h>
#include <stk_tools/block_extractor/ExtractBlocks.hpp>
#include "stk_tools/transfer_utils/TransientFieldTransferById.hpp"
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
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

void extract_blocks_from_file(const std::string &inFile,
                              const std::string &outFile,
                              const std::vector<std::string> &blockNames,
                              MPI_Comm comm)
{
    stk::mesh::MetaData inMeta;
    stk::mesh::BulkData inBulk(inMeta, comm, stk::mesh::BulkData::AUTO_AURA
#ifdef SIERRA_MIGRATION
                               , false
#endif
                               , (stk::mesh::FieldDataManager*)nullptr);

    stk::io::StkMeshIoBroker stkInput;
    stk::io::fill_mesh_preexisting(stkInput, inFile, inBulk);

    stk::mesh::MetaData outMeta;
    stk::mesh::BulkData outBulk(outMeta, comm, stk::mesh::BulkData::AUTO_AURA
#ifdef SIERRA_MIGRATION
                                , false
#endif
                                , (stk::mesh::FieldDataManager*)nullptr);


    extract_blocks(inBulk, outBulk, blockNames);


    remove_io_attribute_from_empty_parts(outBulk);

    stk::io::StkMeshIoBroker stkOutput;
    stkOutput.set_bulk_data(outBulk);
    stkOutput.set_attribute_field_ordering_stored_by_part_ordinal(stkInput.get_attribute_field_ordering_stored_by_part_ordinal());

    stk::transfer_utils::TransientFieldTransferById transfer(stkInput, stkOutput);
    transfer.transfer_and_write_transient_fields(outFile);
}

void extract_blocks_and_ns_from_file(const std::string &inFile,
                              const std::string &outFile,
                              const std::vector<std::string> &blockNames,
							  const std::vector<std::string> &nodesetNames,
                              MPI_Comm comm)
{
    stk::mesh::MetaData inMeta;
    stk::mesh::BulkData inBulk(inMeta, comm, stk::mesh::BulkData::AUTO_AURA
#ifdef SIERRA_MIGRATION
                               , false
#endif
                               , (stk::mesh::FieldDataManager*)nullptr);

    stk::mesh::MetaData outMeta;
    stk::mesh::BulkData outBulk(outMeta, comm, stk::mesh::BulkData::AUTO_AURA
#ifdef SIERRA_MIGRATION
                                , false
#endif
                                , (stk::mesh::FieldDataManager*)nullptr);

    stk::io::StkMeshIoBroker stkInput;
    stk::io::fill_mesh_preexisting(stkInput, inFile, inBulk);

    stk::mesh::PartVector parts;
    for(size_t i=0; i<nodesetNames.size(); i++)
    {
        parts.push_back(inBulk.mesh_meta_data().get_part(nodesetNames[i]));
        ThrowRequireMsg(parts[i] != nullptr, "Can't find nodeset " << nodesetNames[i] << " in mesh.\n");
    }

    for(size_t i=0; i<blockNames.size(); i++)
        {
            parts.push_back(inBulk.mesh_meta_data().get_part(blockNames[i]));
            ThrowRequireMsg(parts[i] != nullptr, "Can't find block " << blockNames[i] << " in mesh.\n");
        }


    stk::mesh::Selector nodeset_and_block_selector = stk::mesh::selectUnion(parts);

    stk::mesh::Selector nothingSelector_byDefaultConstruction;
    stk::mesh::Selector allSelector(!nothingSelector_byDefaultConstruction);

    stk::tools::copy_mesh(inBulk, allSelector, outBulk);

    stk::io::StkMeshIoBroker stkOutput;
    stkOutput.set_bulk_data(outBulk);
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
        ThrowRequireMsg(parts[i] != nullptr, "Can't find block " << blockNames[i] << " in mesh.\n");
    }

    stk::mesh::Selector selector = stk::mesh::selectUnion(parts);
    stk::tools::copy_mesh(oldBulk, selector, newBulk);
}

}
}
