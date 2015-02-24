#include <gtest/gtest.h>
#include <mpi.h>
#include <stk_mesh/base/ModificationSummary.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

TEST(ModificationSummary, testString)
{
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    const std::string generatedMeshSpecification = "generated:1x1x3";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    if (stkMeshBulkData.parallel_size() == 1 )
    {
        ////////////////////////

        stk::ModificationSummary writer(stkMeshBulkData);

        int mod_count = 0;

        stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
        stk::mesh::Entity node2 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 2);

        writer.store_stack_trace();

        writer.track_destroy_entity(node1);

        const stk::mesh::PartVector &addParts = stkMeshBulkData.bucket(node1).supersets();
        const stk::mesh::PartVector &rmParts = stkMeshBulkData.bucket(node2).supersets();

        writer.track_change_entity_parts(node2, addParts, rmParts);
        writer.track_change_entity_parts(node1, addParts, rmParts);

        std::vector<stk::mesh::EntityProc> changes;
        changes.push_back(std::make_pair(node1, 2));
        changes.push_back(std::make_pair(node2, 4));

        writer.track_change_entity_owner(changes);

        stk::mesh::EntityId newId = 12;
        writer.track_change_entity_id(newId, node1);

        writer.track_declare_entity(stk::topology::NODE_RANK, newId, addParts);

        stk::mesh::Entity element1 = stkMeshBulkData.get_entity(stk::topology::ELEM_RANK, 1);
        stk::mesh::Permutation permut = static_cast<stk::mesh::Permutation>(0);
        stk::mesh::RelationIdentifier rel = 0;

        writer.track_declare_relation(element1, node1, rel, permut);
        writer.track_declare_relation(element1, node2, rel, permut);

        writer.track_destroy_relation(element1, node1, rel);

        const stk::mesh::Ghosting &aura = *stkMeshBulkData.ghostings()[1];

        std::vector<stk::mesh::EntityKey> remove_receive;
        remove_receive.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 1));
        remove_receive.push_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK, 1));

        writer.track_change_ghosting(aura, changes , remove_receive);

        writer.write_summary(mod_count);
    }
}

