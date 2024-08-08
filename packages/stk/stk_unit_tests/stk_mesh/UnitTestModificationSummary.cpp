#include <gtest/gtest.h>                // for TEST
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/ModificationSummary.hpp>
#include <string>                       // for string
#include <utility>                      // for make_pair
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Types.hpp"      // for EntityProc, PartVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Ghosting; } }

TEST(ModificationSummary, testString)
{
  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (numprocs == 1 )
  {
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    const std::string generatedMeshSpecification = "generated:1x1x3";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    ////////////////////////

    stk::ModificationSummary writer(stkMeshBulkData);

    int mod_count = 0;

    stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
    stk::mesh::Entity node2 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 2);

    writer.track_destroy_entity(node1);

    const stk::mesh::PartVector &addParts = stkMeshBulkData.bucket(node1).supersets();
    stk::mesh::OrdinalVector addPartsOrdinals, rmPartsOrdinals;
    stkMeshBulkData.bucket(node1).supersets(addPartsOrdinals);
    stkMeshBulkData.bucket(node2).supersets(rmPartsOrdinals);

    writer.track_change_entity_parts(node2, addPartsOrdinals, rmPartsOrdinals);
    writer.track_change_entity_parts(node1, addPartsOrdinals, rmPartsOrdinals);

    std::vector<stk::mesh::EntityProc> changes;
    changes.push_back(std::make_pair(node1, 2));
    changes.push_back(std::make_pair(node2, 4));

    writer.track_change_entity_owner(changes);

    stk::mesh::EntityId newId = 12;
    writer.track_change_entity_id(newId, node1);

    writer.track_declare_entity(stk::topology::NODE_RANK, newId, addParts);

    stk::mesh::Entity element1 = stkMeshBulkData.get_entity(stk::topology::ELEM_RANK, 1);
    stk::mesh::Permutation permut = stk::mesh::Permutation::INVALID_PERMUTATION;
    stk::mesh::RelationIdentifier rel = 0;

    writer.track_declare_relation(element1, node1, rel, permut);
    writer.track_declare_relation(element1, node2, rel, permut);

    writer.track_destroy_relation(element1, node1, rel);

    const stk::mesh::Ghosting &aura = *stkMeshBulkData.ghostings()[1];

    std::vector<stk::mesh::Entity> remove_receive;
    remove_receive.push_back(stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1));
    remove_receive.push_back(stkMeshBulkData.get_entity(stk::topology::ELEM_RANK, 1));

    writer.track_change_ghosting(aura, changes , remove_receive);

    writer.write_summary(mod_count);
    unlink("modification_cycle_P000_B000_C000000.txt");
  }
}

