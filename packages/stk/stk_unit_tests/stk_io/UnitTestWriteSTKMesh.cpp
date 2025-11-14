#include <gtest/gtest.h>
#include <Ioss_IOFactory.h>
#include <Ioss_Region.h>
#include <Ioss_DBUsage.h>
#include <Ioss_PropertyManager.h>
#include <Ionit_Initializer.h>

#include <Ioss_ElementBlock.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_SideSet.h>
#include <Ioss_NodeSet.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_mesh/base/SideSetUtil.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/memory_util.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>

#include <string>
#include <iostream>
#include <unistd.h>                     // for unlink

#include <stk_unit_test_utils/BuildMesh.hpp>

namespace {

using stk::unit_test_util::build_mesh;
using stk::unit_test_util::build_mesh_no_simple_fields;

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

void verify_num_nodes_in_file(MPI_Comm comm,
                              const std::string& meshFileName,
                              unsigned expectedNumNodes)
{
    std::shared_ptr<stk::mesh::BulkData> bulkData = build_mesh(comm);
    stk::io::fill_mesh(meshFileName, *bulkData);

    std::vector<size_t> entity_counts;
    stk::mesh::comm_mesh_counts(*bulkData, entity_counts);
    EXPECT_EQ(expectedNumNodes, entity_counts[stk::topology::NODE_RANK]);
}

void fill_node_ids_and_coords(const stk::mesh::BulkData& bulk,
                              std::vector<int64_t>& node_ids,
                              std::vector<double>& coordinates)
{
    const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    stk::mesh::Field<double> * coordField = meta.get_field<double>(stk::topology::NODE_RANK, "coordinates");
    int spatial_dim = meta.spatial_dimension();

    STK_ThrowAssert(coordField != NULL);
    auto coordFieldData = coordField->data();

    stk::mesh::Selector locallyOwned = meta.locally_owned_part();

    int node_counter = 0;
    stk::mesh::for_each_entity_run_no_threads(bulk, stk::topology::NODE_RANK, locallyOwned,
      [&](const stk::mesh::BulkData& mesh, stk::mesh::Entity node)
      {
          int node_id = mesh.identifier(node);
          node_ids[node_counter] = node_id;

          auto coords = coordFieldData.entity_values(node);
          for(stk::mesh::ComponentIdx k=0_comp;k<spatial_dim;++k) {
              coordinates[spatial_dim*node_counter+k] = coords(k);
          }
          node_counter++;
      });
}

void fill_elem_ids_and_connectivity(const stk::mesh::BulkData& bulkData,
                                    const stk::mesh::Part* elemBlock,
                                    std::vector<int64_t>& elem_ids,
                                    std::vector<int64_t>& connectivity)
{
    stk::mesh::EntityVector elems;
    stk::mesh::get_entities(bulkData, stk::topology::ELEM_RANK, *elemBlock, elems);

    elem_ids.resize(elems.size());
    const unsigned connectivity_size = elems.size()*elemBlock->topology().num_nodes();
    connectivity.resize(connectivity_size);

    unsigned conn_counter = 0;
    for(size_t j=0;j<elems.size();++j) {
        elem_ids[j] = bulkData.identifier(elems[j]);
        const stk::mesh::ConnectedEntities nodes = bulkData.get_connected_entities(elems[j], stk::topology::NODE_RANK);
        for(unsigned k=0;k<nodes.size();++k) {
            connectivity[conn_counter++] = bulkData.identifier(nodes[k]);
        }
    }
}

TEST(StkIo, write_stk_mesh_using_ioss_instead_of_stk_io)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    std::string file_written = "out.exo";

    if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

    {
        std::shared_ptr<stk::mesh::BulkData> bulkData = build_mesh(comm);
        stk::mesh::MetaData& meta = bulkData->mesh_meta_data();

        stk::io::fill_mesh("generated:2x2x2|sideset:xX|nodeset:x", *bulkData);

        Ioss::DatabaseIO* db_io = create_output_db_io(file_written);
        EXPECT_TRUE(db_io->ok());
        Ioss::Region output_region(db_io);

        ////////////////////////////////////////////////////////////

        output_region.begin_mode(Ioss::STATE_DEFINE_MODEL);

        std::string NodeBlockName = "nodeblock_1";

        std::vector<size_t> entity_counts;
        stk::mesh::comm_mesh_counts(*bulkData, entity_counts);

        int64_t num_nodes = entity_counts[stk::topology::NODE_RANK];
        int spatial_dim = meta.spatial_dimension();

        Ioss::NodeBlock *output_node_block = new Ioss::NodeBlock(db_io, NodeBlockName, num_nodes, spatial_dim);
        output_region.add(output_node_block);

        stk::mesh::PartVector elemBlockParts;
        stk::mesh::fill_element_block_parts(meta, stk::topology::HEX_8, elemBlockParts);

        for(const stk::mesh::Part* elemBlock : elemBlockParts) {
            STK_ThrowRequireMsg(stk::io::is_part_io_part(*elemBlock),"element-block-part "<<elemBlock->name()<<" is not an IO part.");
            unsigned numElems = stk::mesh::count_entities(*bulkData, stk::topology::ELEM_RANK, *elemBlock);
            Ioss::ElementBlock *output_element_block = new Ioss::ElementBlock(db_io, elemBlock->name(), elemBlock->topology().name(), numElems);

            output_element_block->property_add(Ioss::Property("original_topology_type", elemBlock->topology().name()));
            output_element_block->property_add(Ioss::Property("id", elemBlock->id()));
            output_region.add(output_element_block);

            // how about attributes?
        }

        output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

        ////////////////////////////////////////////////////////////

        output_region.begin_mode(Ioss::STATE_MODEL);

        Ioss::NodeBlock *node_block = output_region.get_node_blocks()[0];

        std::vector<double> coordinates(spatial_dim*num_nodes);
        std::vector<int64_t> node_ids(num_nodes);

        fill_node_ids_and_coords(*bulkData, node_ids, coordinates);

        node_block->put_field_data("mesh_model_coordinates", coordinates);
        node_block->put_field_data("ids", node_ids);

        for(const stk::mesh::Part* elemBlock : elemBlockParts) {
            std::vector<int64_t> elem_ids;
            std::vector<int64_t> connectivity;

            fill_elem_ids_and_connectivity(*bulkData, elemBlock, elem_ids, connectivity);

            Ioss::ElementBlock *output_element_block = output_region.get_element_block(elemBlock->id());
            output_element_block->put_field_data("connectivity_raw", connectivity);
            output_element_block->put_field_data("ids", elem_ids);
        }

        output_region.end_mode(Ioss::STATE_MODEL);
        ////////////////////////////////////////////////////////////
    }

    verify_num_nodes_in_file(comm, file_written, 27);
    unlink(file_written.c_str());
}

class StkIoResultsOutput : public stk::unit_test_util::MeshFixture
{
protected:
    void setup_mesh(const std::string & meshSpec,
                    stk::mesh::BulkData::AutomaticAuraOption auraOption,
                    unsigned initialBucketCapacity = stk::mesh::get_default_initial_bucket_capacity(),
                    unsigned maximumBucketCapacity = stk::mesh::get_default_maximum_bucket_capacity()) override
    {
        setup_empty_mesh(auraOption, initialBucketCapacity, maximumBucketCapacity);

        stk::mesh::Field<int> & field = get_meta().declare_field<int>(stk::topology::NODE_RANK, "nodal_field");
        stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);

        stk::io::fill_mesh(meshSpec, get_bulk());
    }

    const stk::mesh::Part& setup_mesh_with_part(const std::string & meshSpec, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        setup_empty_mesh(auraOption);

        stk::mesh::Part& surface_part = get_meta().declare_part_with_topology("surface_1", stk::topology::QUADRILATERAL_4);
        stk::io::put_io_part_attribute(surface_part);

        stk::mesh::Field<int> & field = get_meta().declare_field<int>(stk::topology::NODE_RANK, "nodal_field");
        stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);

        stk::io::fill_mesh(meshSpec, get_bulk());

        return surface_part;
    }
};


TEST_F(StkIoResultsOutput, close_output_mesh_makes_it_invalid) {
    if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

    std::string meshSpec = stk::unit_test_util::get_option("--mesh-spec", "generated:1x1x2|sideset:z");

    setup_mesh(meshSpec, stk::mesh::BulkData::NO_AUTO_AURA);

    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(get_bulk());

    std::string fileName1 = "output1.e";

    stk::mesh::FieldBase * nodalField = get_meta().get_field(stk::topology::NODE_RANK, "nodal_field");

    size_t outputFileIndex = stkIo.create_output_mesh(fileName1, stk::io::WRITE_RESULTS);
    stkIo.add_field(outputFileIndex, *nodalField, stk::topology::NODE_RANK, "nodal_field");

    stkIo.close_output_mesh(outputFileIndex);
    unlink((fileName1).c_str());

    EXPECT_THROW(stkIo.add_field(outputFileIndex, *nodalField, stk::topology::NODE_RANK, "nodal_field2"), std::runtime_error);


}

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

    const unsigned spatialDim = 3;
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

    EXPECT_NO_THROW(stk::io::fill_mesh(fileName, *bulk));
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

    EXPECT_FALSE( stk::mesh::should_reconstruct_sideset(bulk, *surface_1) );
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

    EXPECT_TRUE( stk::mesh::should_reconstruct_sideset(bulk, surface_part) );

    bulk.create_sideset(surface_part);

    EXPECT_FALSE( stk::mesh::should_reconstruct_sideset(bulk, surface_part) );

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

    EXPECT_FALSE( stk::mesh::should_reconstruct_sideset(bulk, surface_part) );
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

    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

    if (autoDecomp == "false") {
        stk::io::fill_mesh(meshSpec, *bulk);
    }
    else {
        stk::io::fill_mesh_with_auto_decomp(meshSpec, *bulk);
    }

    stk::io::write_mesh("readWriteTest.exo", *bulk);
}

}
