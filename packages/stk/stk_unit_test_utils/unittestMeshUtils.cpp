
#include "unittestMeshUtils.hpp"

#include <ext/alloc_traits.h>
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <sstream>                      // for basic_ostream::operator<<, etc
#include <stk_io/IossBridge.hpp>        // for put_io_part_attribute
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_topology/topology.hpp>    // for topology, operator++, etc
#include <vector>                       // for vector
#include "ioUtils.hpp"                  // for fill_mesh_using_stk_io, etc
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/CoordinateSystems.hpp"  // for Cartesian2d
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityVector, etc
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowRequireMsg
namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class Part; } }

namespace stk
{
namespace unit_test_util
{


void put_mesh_into_part(stk::mesh::BulkData& bulkData, stk::mesh::Part& part)
{
    stk::mesh::EntityVector entitiesToMakeActive;
    std::vector<stk::mesh::PartVector> add_parts;
    std::vector<stk::mesh::PartVector> rm_parts;
    for(stk::topology::rank_t rank=stk::topology::BEGIN_RANK; rank < bulkData.mesh_meta_data().entity_rank_count(); rank++)
    {
        const stk::mesh::BucketVector &buckets = bulkData.get_buckets(rank, bulkData.mesh_meta_data().locally_owned_part());
        for(const stk::mesh::Bucket *bucket : buckets)
        {
            for(stk::mesh::Entity entity : *bucket)
            {
                entitiesToMakeActive.push_back(entity);
                add_parts.push_back(stk::mesh::PartVector(1, &part));
                rm_parts.push_back(stk::mesh::PartVector());
            }
        }
    }
    bulkData.batch_change_entity_parts(entitiesToMakeActive, add_parts, rm_parts);
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
    std::vector<stk::mesh::PartVector> add_parts(killedElements.size());
    std::vector<stk::mesh::PartVector> rm_parts(killedElements.size());

    for (size_t j=0;j<killedElements.size();++j)
    {
        rm_parts[j] = removeParts;
    }

    bulkData.batch_change_entity_parts(killedElements, add_parts, rm_parts);
}

} // namespace unit_test_util
} // namespace stk

