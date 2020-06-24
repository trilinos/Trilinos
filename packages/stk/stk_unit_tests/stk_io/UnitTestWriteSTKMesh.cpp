#include <gtest/gtest.h>
#include "mpi.h"
#include <Ioss_IOFactory.h>
#include <Ioss_Region.h>
#include <Ioss_DBUsage.h>
#include <Ioss_PropertyManager.h>
#include <string>
#include <init/Ionit_Initializer.h>

#include <Ioss_ElementBlock.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_SideSet.h>
#include <Ioss_SideBlock.h>
#include <Ioss_NodeSet.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/memory_util.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>

#include <iostream>
#include <unistd.h>                     // for unlink

#include "stk_util/environment/Env.hpp"
namespace
{

Ioss::DatabaseIO* create_output_db_io(const std::string &filename)
{
    Ioss::Init::Initializer init_db;

    Ioss::DatabaseUsage db_usage = Ioss::WRITE_RESULTS;
    MPI_Comm communicator = MPI_COMM_WORLD;
    Ioss::PropertyManager properties;

    properties.add(Ioss::Property("INTEGER_SIZE_DB",  8));
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));

    Ioss::DatabaseIO *db_io = Ioss::IOFactory::create("exodus", filename, db_usage,
						      communicator, properties);
    return db_io;
}

//BeginDocTest1
TEST(StkIo, write_stk_mesh_to_file)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    std::string file_written = "out.exo";

    if(stk::parallel_machine_size(comm) == 1)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulkData(meta, comm);
        stk::io::fill_mesh("generated:2x2x2|sideset:xX|nodeset:x", bulkData);

        const stk::mesh::PartVector & all_parts = meta.get_parts();

        Ioss::DatabaseIO* db_io = create_output_db_io(file_written);
        Ioss::Region output_region(db_io);
        EXPECT_TRUE(db_io->ok());

        ////////////////////////////////////////////////////////////

        output_region.begin_mode(Ioss::STATE_DEFINE_MODEL);

        std::string NodeBlockName = "nodeblock_1";

        std::vector<size_t> entity_counts;
        stk::mesh::comm_mesh_counts(bulkData, entity_counts);

        int64_t num_nodes = entity_counts[stk::topology::NODE_RANK];
        int spatial_dim = meta.spatial_dimension();

        Ioss::NodeBlock *output_node_block = new Ioss::NodeBlock(db_io, NodeBlockName, num_nodes, spatial_dim);
        output_region.add(output_node_block);

        for(stk::mesh::PartVector::const_iterator i = all_parts.begin(); i != all_parts.end(); ++i)
        {
            stk::mesh::Part * const part = *i;

            if(stk::io::is_part_io_part(*part)) // this means it is an io_part
            {
                if(part->primary_entity_rank() == stk::topology::NODE_RANK)
                {
                }
                else if(part->primary_entity_rank() == stk::topology::ELEMENT_RANK)
                {
                    //
                    stk::mesh::EntityVector entities;
                    const stk::mesh::BucketVector &input_buckets = bulkData.buckets(stk::topology::ELEMENT_RANK);
                    stk::mesh::get_selected_entities(*part, input_buckets, entities);
                    Ioss::ElementBlock *output_element_block = new Ioss::ElementBlock(db_io, part->name(), part->topology().name(), entities.size());

                    output_element_block->property_add(Ioss::Property("original_topology_type", part->topology().name()));
                    output_element_block->property_add(Ioss::Property("id", part->id()));
                    output_region.add(output_element_block);

                    // how about attributes?

                }
                else if(part->primary_entity_rank() == stk::topology::FACE_RANK)
                {
                    //
                }
                //else if(part->primary_entity_rank() == stk::topology::EDGE_RANK)
                //{
                //    //
                //}
            }
        }

        output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

        ////////////////////////////////////////////////////////////

        output_region.begin_mode(Ioss::STATE_MODEL);

        Ioss::NodeBlock *node_block = output_region.get_node_blocks()[0];

        stk::mesh::Field<double,stk::mesh::Cartesian> * coordField =
          meta.get_field<stk::mesh::Field<double,stk::mesh::Cartesian> >(stk::topology::NODE_RANK, "coordinates");

        ASSERT_TRUE(coordField != NULL);

        std::vector<double> coordinates(spatial_dim*num_nodes);
        std::vector<int64_t> node_ids(num_nodes);

        stk::mesh::Selector local_nodes = meta.locally_owned_part();

        const stk::mesh::BucketVector &input_buckets = bulkData.get_buckets(stk::topology::NODE_RANK, local_nodes);

        int node_counter = 0;
        for(size_t i=0;i<input_buckets.size();++i)
        {
            const stk::mesh::Bucket &bucket = *input_buckets[i];
            for(size_t j=0;j<bucket.size();++j)
            {
                stk::mesh::Entity node = bucket[j];
                int node_id = bulkData.identifier(node);
                node_ids[node_counter] = node_id;

                double* coords = stk::mesh::field_data(*coordField, node);
                for(int k=0;k<spatial_dim;++k)
                {
                    coordinates[spatial_dim*node_counter+k] = coords[k];
                }
                node_counter++;
            }
        }

        node_block->put_field_data("mesh_model_coordinates", coordinates);
        node_block->put_field_data("ids", node_ids);

        for(stk::mesh::PartVector::const_iterator i = all_parts.begin(); i != all_parts.end(); ++i)
        {
            stk::mesh::Part * const part = *i;

            if(stk::io::is_part_io_part(*part)) // this means it is an io_part
            {
                if(part->primary_entity_rank() == stk::topology::ELEMENT_RANK)
                {
                    //
                    stk::mesh::EntityVector entities;
                    const stk::mesh::BucketVector &input_bucketsA = bulkData.buckets(stk::topology::ELEMENT_RANK);
                    stk::mesh::get_selected_entities(*part, input_bucketsA, entities);
                    Ioss::ElementBlock *output_element_block = output_region.get_element_block(part->id());

                    std::vector<int64_t> elem_ids(entities.size());
                    unsigned connectivity_size = entities.size()*part->topology().num_nodes();
                    std::vector<int64_t> connectivity(connectivity_size);
                    unsigned conn_counter = 0;

                    for(size_t j=0;j<entities.size();++j)
                    {
                        elem_ids[j] = bulkData.identifier(entities[j]);
                        unsigned num_nodes_per = bulkData.num_nodes(entities[j]);
                        const stk::mesh::Entity *nodes = bulkData.begin_nodes(entities[j]);
                        for(unsigned k=0;k<num_nodes_per;++k)
                        {
                            connectivity[conn_counter] = bulkData.identifier(nodes[k]);
                            conn_counter++;
                        }
                    }

                    output_element_block->put_field_data("connectivity_raw", connectivity);
                    output_element_block->put_field_data("ids", elem_ids);
                }
            }
        }

        output_region.end_mode(Ioss::STATE_MODEL);
        ////////////////////////////////////////////////////////////
    }

    if(stk::parallel_machine_size(comm) == 1)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulkData(meta, comm);
        stk::io::fill_mesh(file_written, bulkData);

        std::vector<size_t> entity_counts;
        stk::mesh::comm_mesh_counts(bulkData, entity_counts);
        EXPECT_EQ(27u, entity_counts[stk::topology::NODE_RANK]);
    }

    unlink(file_written.c_str());

}
//EndDocTest1

void print_memory(size_t &current, size_t &hwm)
{
    static int counter = 1;
    std::cerr << "Current(" << counter << ") : " << current << "\thwm: " << hwm << std::endl;
    ++counter;
}

TEST(StkIo, check_memory)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    size_t current_usage = 0, hwm_usage = 0;
    current_usage = stk::get_memory_usage_now();
    print_memory(current_usage, hwm_usage);

    if(stk::parallel_machine_size(comm) == 1)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulkData(meta, comm);
        std::string filename = stk::unit_test_util::get_mesh_spec("--dim");

        size_t current_usage2 = 0, hwm_usage2 = 0;
        {
            size_t current_usage1 = 0, hwm_usage1 = 0;
            current_usage1 = stk::get_memory_usage_now();
            print_memory(current_usage1, hwm_usage1);

            stk::io::StkMeshIoBroker exodusFileReader(comm);
            exodusFileReader.set_bulk_data(bulkData);
            exodusFileReader.add_mesh_database(filename, stk::io::READ_MESH);
            exodusFileReader.create_input_mesh();
            exodusFileReader.populate_bulk_data();

            current_usage2 = stk::get_memory_usage_now();
            print_memory(current_usage2, hwm_usage2);
        }

        size_t current_usage3 = 0, hwm_usage3 = 0;
        current_usage3 = stk::get_memory_usage_now();
        print_memory(current_usage3, hwm_usage3);
        size_t padDueToMemNoise = 8192;
        EXPECT_LE(current_usage2, (current_usage3+padDueToMemNoise));
    }
}

class StkIoResultsOutput : public stk::unit_test_util::MeshFixture
{
protected:
    void setup_mesh(const std::string & meshSpec, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        setup_empty_mesh(auraOption);

        stk::mesh::Field<int> & field = get_meta().declare_field<stk::mesh::Field<int>>(stk::topology::NODE_RANK, "nodal_field");
        const int initValue = 0;
        stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &initValue);

        stk::io::fill_mesh(meshSpec, get_bulk());
    }

    const stk::mesh::Part& setup_mesh_with_part(const std::string & meshSpec, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        setup_empty_mesh(auraOption);

        stk::mesh::Part& surface_part = get_meta().declare_part_with_topology("surface_1", stk::topology::QUADRILATERAL_4);
        stk::io::put_io_part_attribute(surface_part);

        stk::mesh::Field<int> & field = get_meta().declare_field<stk::mesh::Field<int>>(stk::topology::NODE_RANK, "nodal_field");
        const int initValue = 0;
        stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &initValue);

        stk::io::fill_mesh(meshSpec, get_bulk());

        return surface_part;
    }
};

TEST_F(StkIoResultsOutput, write_nodal_face_variable_multiple_procs)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

    std::string meshSpec = stk::unit_test_util::get_option("--mesh-spec", "generated:1x1x2|sideset:z");

    setup_mesh(meshSpec, stk::mesh::BulkData::NO_AUTO_AURA);

    const std::string fileName = "nodal_field_as_face_variable.e";
    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(get_bulk());
    size_t outputFileIndex = stkIo.create_output_mesh(fileName, stk::io::WRITE_RESULTS);
    stkIo.use_nodeset_for_sideset_nodes_fields(outputFileIndex, true);
    stkIo.check_field_existence_when_creating_nodesets(outputFileIndex, false);

    stk::mesh::FieldBase * nodalField = get_meta().get_field(stk::topology::NODE_RANK, "nodal_field");
    ASSERT_TRUE(nodalField != nullptr);
    stkIo.add_field(outputFileIndex, *nodalField, stk::topology::FACE_RANK, "nodal_field");
    stkIo.write_output_mesh(outputFileIndex);
    stkIo.begin_output_step(outputFileIndex, 0.0);
    stkIo.write_defined_output_fields(outputFileIndex);
    stkIo.end_output_step(outputFileIndex);

    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
    EXPECT_NO_THROW(stk::io::fill_mesh(fileName, bulk));
    unlink(fileName.c_str());
}

TEST_F(StkIoResultsOutput, no_reconstruct_on_input)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

    std::string meshSpec = stk::unit_test_util::get_option("--mesh-spec", "generated:1x1x2|sideset:z");

    setup_mesh(meshSpec, stk::mesh::BulkData::NO_AUTO_AURA);

    const stk::mesh::Part* surface_1 = get_meta().get_part("surface_1");
    EXPECT_TRUE(surface_1 != nullptr);

    const stk::mesh::BulkData& bulk = get_bulk();
    EXPECT_TRUE(bulk.does_sideset_exist(*surface_1));

    EXPECT_FALSE( stk::io::should_reconstruct_sideset(bulk, *surface_1) );
}

TEST_F(StkIoResultsOutput, reconstruct_on_creating_sideset)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

    std::string meshSpec = stk::unit_test_util::get_option("--mesh-spec", "generated:1x1x1");

    const stk::mesh::Part& surface_part = setup_mesh_with_part(meshSpec, stk::mesh::BulkData::NO_AUTO_AURA);
    const stk::mesh::Part& block_1 = *get_meta().get_part("block_1");
    get_meta().set_surface_to_block_mapping(&surface_part, {&block_1});

    stk::mesh::BulkData& bulk = get_bulk();
    EXPECT_FALSE(bulk.does_sideset_exist(surface_part));

    EXPECT_TRUE( stk::io::should_reconstruct_sideset(bulk, surface_part) );

    bulk.create_sideset(surface_part);

    EXPECT_FALSE( stk::io::should_reconstruct_sideset(bulk, surface_part) );

    stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 1);
    EXPECT_TRUE(bulk.is_valid(elem));

    bulk.modification_begin();
    bulk.declare_element_side(elem, 1, stk::mesh::ConstPartVector{&surface_part});
    bulk.modification_end();

    EXPECT_TRUE(bulk.does_sideset_exist(surface_part));
    stk::mesh::SideSet& ss = bulk.get_sideset(surface_part);
    EXPECT_EQ(1u, ss.size());
    stk::mesh::SideSetEntry entry = ss[0];
    EXPECT_EQ(elem, entry.element);
    EXPECT_EQ(1, entry.side);

    EXPECT_FALSE( stk::io::should_reconstruct_sideset(bulk, surface_part) );
}

TEST(TestStkIo, readWrite)
{
    std::string meshSpec = stk::unit_test_util::get_option("--mesh", "none specified");
    if (meshSpec == "none specified") {
        if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
            std::cout<<"No mesh specified, exiting."<<std::endl;
        }
        return;
    }

    std::string autoDecomp = stk::unit_test_util::get_option("--auto-decomp", "false");

    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

    if (autoDecomp == "false") {
        stk::io::fill_mesh(meshSpec, bulk);
    }
    else {
        stk::io::fill_mesh_with_auto_decomp(meshSpec, bulk);
    }

    stk::io::write_mesh("readWriteTest.exo", bulk);
}

}
