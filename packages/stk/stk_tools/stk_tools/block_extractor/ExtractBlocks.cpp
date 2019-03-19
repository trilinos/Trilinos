#include "mpi.h"
#include <stk_util/stk_config.h>
#include <stk_tools/block_extractor/ExtractBlocks.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_tools/mesh_clone/MeshClone.hpp>
#include "stk_mesh/base/Selector.hpp"
#include <stk_io/StkMeshIoBroker.hpp>

namespace stk {
namespace tools {

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

    int numSteps = 9;
    double maxTime = 0;
    stk::io::fill_mesh_save_step_info(inFile, inBulk, numSteps, maxTime);

    stk::mesh::MetaData outMeta;
    stk::mesh::BulkData outBulk(outMeta, comm, stk::mesh::BulkData::AUTO_AURA
#ifdef SIERRA_MIGRATION
, false
#endif
, (stk::mesh::FieldDataManager*)nullptr);
    extract_blocks(inBulk, outBulk, blockNames);

    stk::io::write_mesh_with_fields(outFile, outBulk, numSteps, maxTime);
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
