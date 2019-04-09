#include <gtest/gtest.h>
#include "mpi.h"
#include "Ioss_Region.h"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_io/ProcessSetsOrBlocks.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_io/InputFile.hpp"         // for InputFile
#include "stk_unit_test_utils/ioUtils.hpp"
#include <stk_balance/fixSplitCoincidentElements.hpp>

class StkMeshIoBrokerTester : public stk::io::StkMeshIoBroker
{
public:

    StkMeshIoBrokerTester(stk::ParallelMachine comm) :
        stk::io::StkMeshIoBroker(comm)
    {
        this->set_sideset_face_creation_behavior_for_testing(STK_IO_SIDE_CREATION_USING_GRAPH_TEST);
    }

    virtual ~StkMeshIoBrokerTester() {}

    virtual stk::mesh::EntityIdProcMap move_elements()
    {
        stk::mesh::EntityIdProcMap elemIdMovedToProc;
        stk::mesh::EntityProcVec elemToMove;
        stk::mesh::Entity elem2 = bulk_data().get_entity(stk::topology::ELEM_RANK, 2);
        if(bulk_data().bucket(elem2).owned())
        {
            int otherProc = 1 - bulk_data().parallel_rank();
            elemToMove.push_back({elem2, otherProc});
            elemIdMovedToProc[bulk_data().identifier(elem2)] = otherProc;
        }
        bulk_data().change_entity_owner(elemToMove);
        return elemIdMovedToProc;
    }

    virtual void populate_mesh(bool delay_field_data_allocation)
    {
        validate_input_file_index(m_activeMeshIndex);

        create_bulk_data();

        if (delay_field_data_allocation) {
          bulk_data().deactivate_field_updating();
        }

        bool i_started_modification_cycle = bulk_data().modification_begin("Mesh Read - step 1");

        Ioss::Region *region = m_inputFiles[m_activeMeshIndex]->get_input_io_region().get();
        bool ints64bit = stk::io::db_api_int_size(region) == 8;
        if (ints64bit) {
          stk::io::process_nodeblocks<int64_t>(*region,    bulk_data());
          stk::io::process_elementblocks<int64_t>(*region, bulk_data());
          stk::io::process_nodesets<int64_t>(*region,      bulk_data());
        } else {
            stk::io::process_nodeblocks<int>(*region,    bulk_data());
            stk::io::process_elementblocks<int>(*region, bulk_data());
            stk::io::process_nodesets<int>(*region,      bulk_data());
        }

        stk_mesh_resolve_node_sharing();

        bulk_data().modification_end();

        stk::mesh::EntityIdProcMap elemIdMovedToProc = move_elements();

        bulk_data().modification_begin();

        bulk_data().initialize_face_adjacent_element_graph();
        stk::io::process_sidesets(*region, bulk_data(), elemIdMovedToProc, m_sidesetFaceCreationBehavior);

        bulk_data().modification_end();

        // stk_mesh_modification_end_after_node_sharing_resolution();

        // Not sure if this is needed anymore. Don't think it'll be called with a nested modification cycle
        if(!i_started_modification_cycle)
            bulk_data().modification_begin();
    }

};

void read_mesh_using_iobroker(stk::io::StkMeshIoBroker & broker, const std::string& fileName, stk::mesh::BulkData &mesh, const std::string &decompositionMethod)
{
    broker.set_bulk_data(mesh);
    broker.property_add(Ioss::Property("DECOMPOSITION_METHOD", decompositionMethod));
    broker.add_mesh_database(fileName, stk::io::READ_MESH);
    broker.create_input_mesh();
    broker.populate_bulk_data();
}

void read_from_serial_file_and_decompose_tester(const std::string& fileName, stk::mesh::BulkData &mesh, const std::string &decompositionMethod)
{
    StkMeshIoBrokerTester broker(MPI_COMM_NULL);
    read_mesh_using_iobroker(broker, fileName, mesh, decompositionMethod);
}

class Pmr1StkMeshIoBrokerTester : public StkMeshIoBrokerTester
{
public:
    Pmr1StkMeshIoBrokerTester(stk::ParallelMachine comm) : StkMeshIoBrokerTester(comm) { }

    virtual ~Pmr1StkMeshIoBrokerTester() {}

    virtual stk::mesh::EntityIdProcMap move_elements()
    {
        return stk::balance::make_mesh_consistent_with_parallel_mesh_rule1(bulk_data());
    }
};

void read_from_serial_file_and_decompose_using_pmr_tester(const std::string& fileName, stk::mesh::BulkData &mesh, const std::string &decompositionMethod)
{
    Pmr1StkMeshIoBrokerTester broker(MPI_COMM_NULL);
    read_mesh_using_iobroker(broker, fileName, mesh, decompositionMethod);
}

TEST(MeshWithElementAndSide, elementIsMovedDuringReadUsingChangeEntityOwner_sideIsAlsoMoved)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
        read_from_serial_file_and_decompose_tester("generated:1x1x3|sideset:Y", bulk, "cyclic");
        stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, 2);
        if(bulk.is_valid(elem2))
        {
            EXPECT_EQ(1u, bulk.num_faces(elem2));
        }
    }
}

TEST(MeshWithElementAndSide, elementIsMovedDuringReadWhileEnforcingPMR1_sideIsAlsoMoved)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
        EXPECT_NO_THROW(read_from_serial_file_and_decompose_using_pmr_tester("ARefLA.e", bulk, "cyclic"));

        stk::mesh::EntityVector elements;
        stk::mesh::get_selected_entities(meta.locally_owned_part(), bulk.buckets(stk::topology::ELEM_RANK), elements);
        for(stk::mesh::Entity element : elements)
        {
            if(bulk.bucket(element).topology().is_shell())
            {
                EXPECT_EQ(2u, bulk.num_faces(element));
            }
            else
            {
                EXPECT_EQ(1u, bulk.num_faces(element));
            }
        }
    }
}

