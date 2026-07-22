
// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "unittestMeshUtils.hpp"
#include <stddef.h>                    // for size_t
#include <sstream>                     // for basic_ostream::operator<<, etc
#include <stk_mesh/base/BulkData.hpp>  // for BulkData
#include <stk_mesh/base/Entity.hpp>    // for Entity
#include <stk_mesh/base/MetaData.hpp>  // for MetaData
#include <stk_topology/topology.hpp>   // for operator++, topology, etc
#include <string>                      // for allocator, operator<<, string, etc
#include <vector>                      // for vector
#include "stk_mesh/base/Types.hpp"     // for PartVector, EntityVector, etc
namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class Part; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk
{
namespace unit_test_util
{


void put_mesh_into_part(stk::mesh::BulkData& bulkData, stk::mesh::Part& part)
{
  stk::mesh::Selector all = bulkData.mesh_meta_data().universal_part();
  stk::mesh::PartVector addParts = {&part};
  stk::mesh::PartVector emptyRemoveParts;
  for(stk::topology::rank_t rank=stk::topology::BEGIN_RANK; rank < bulkData.mesh_meta_data().entity_rank_count(); rank++)
  {
    bulkData.batch_change_entity_parts(all, rank, addParts, emptyRemoveParts);
  }
}


std::string get_name_of_generated_mesh(int xdim, int ydim, int zdim, const std::string &options)
{
    std::ostringstream os;
    os << "generated:" << xdim << "x" << ydim << "x" << zdim << options;
    return os.str();
}


void move_killed_elements_out_of_parts(stk::mesh::BulkData& bulkData,
                                  const stk::mesh::EntityVector& killedElements,
                                  const stk::mesh::PartVector& removeParts)
{
    stk::mesh::PartVector emptyAddParts;

    bulkData.batch_change_entity_parts(killedElements, emptyAddParts, removeParts);
}

void put_elements_into_part(stk::mesh::BulkData& bulk, const std::vector<ElementAndPart> & entries)
{
  stk::mesh::MetaData & meta = bulk.mesh_meta_data();

  for (const ElementAndPart & entry : entries) {
    const stk::mesh::Part * part = meta.get_part(entry.partName);
    if (part == nullptr) {
      meta.declare_part(entry.partName, stk::topology::ELEM_RANK);
    }
  }

  bulk.modification_begin();
  for (const ElementAndPart & entry : entries) {
    const stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, entry.id);
    if (bulk.is_valid(elem) && bulk.bucket(elem).owned()) {
      const stk::mesh::Part * part = meta.get_part(entry.partName);
      bulk.change_entity_parts(elem, stk::mesh::ConstPartVector{part});
    }
  }
  bulk.modification_end();
}

} // namespace unit_test_util
} // namespace stk

