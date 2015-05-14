#include <gtest/gtest.h>
#include <exo_fpp/Iofx_DatabaseIO.h>
#include <mpi.h>
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
#include "stk_mesh/base/Field.hpp"

#include <stk_util/parallel/Parallel.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>

#include <iostream>
#include <unistd.h>                     // for unlink

namespace
{

Iofx::DatabaseIO* create_output_db_io(const std::string &filename)
{
    Ioss::Init::Initializer init_db;

    Ioss::DatabaseUsage db_usage = Ioss::WRITE_RESULTS;
    MPI_Comm communicator = MPI_COMM_WORLD;
    Ioss::PropertyManager properties;

    properties.add(Ioss::Property("INTEGER_SIZE_DB",  8));
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));

    Iofx::DatabaseIO *db_io = new Iofx::DatabaseIO(NULL, filename, db_usage, communicator, properties);
    return db_io;
}

TEST(StkIo, write_stk_mesh_to_file)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    std::string file_written = "out.exo";

    if(stk::parallel_machine_size(comm) == 1)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulkData(meta, comm);
        stk::unit_test_util::fill_mesh_using_stk_io("generated:2x2x2|sideset:xX|nodeset:x", bulkData, comm);

        const stk::mesh::PartVector & all_parts = meta.get_parts();

        Iofx::DatabaseIO* db_io = create_output_db_io(file_written);
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

            if(NULL != part->attribute<Ioss::GroupingEntity>()) // this means it is an io_part
            {
                if(part->primary_entity_rank() == stk::topology::NODE_RANK)
                {
                    //
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
                else if(part->primary_entity_rank() == stk::topology::EDGE_RANK)
                {
                    //
                }
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

            if(NULL != part->attribute<Ioss::GroupingEntity>()) // this means it is an io_part
            {
                if(part->primary_entity_rank() == stk::topology::ELEMENT_RANK)
                {
                    //
                    stk::mesh::EntityVector entities;
                    const stk::mesh::BucketVector &input_buckets = bulkData.buckets(stk::topology::ELEMENT_RANK);
                    stk::mesh::get_selected_entities(*part, input_buckets, entities);
                    Ioss::ElementBlock *output_element_block = output_region.get_element_block(part->id());

                    std::vector<int64_t> elem_ids(entities.size());
                    unsigned connectivity_size = entities.size()*part->topology().num_nodes();
                    std::vector<int64_t> connectivity(connectivity_size);
                    unsigned conn_counter = 0;

                    for(size_t j=0;j<entities.size();++j)
                    {
                        elem_ids[j] = bulkData.identifier(entities[j]);
                        unsigned num_nodes = bulkData.num_nodes(entities[j]);
                        const stk::mesh::Entity *nodes = bulkData.begin_nodes(entities[j]);
                        for(unsigned k=0;k<num_nodes;++k)
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

        if(stk::parallel_machine_size(comm) == 1)
        {
            stk::mesh::MetaData meta;
            stk::mesh::BulkData bulkData(meta, comm);
            stk::unit_test_util::fill_mesh_using_stk_io(file_written, bulkData, comm);

            std::vector<size_t> entity_counts;
            stk::mesh::comm_mesh_counts(bulkData, entity_counts);
            EXPECT_EQ(27u, entity_counts[stk::topology::NODE_RANK]);
        }

        unlink(file_written.c_str());

        ////////////////////////////////////////////////////////////
    }
}


}
