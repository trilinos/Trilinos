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

namespace {

void test_topology(const Ioss::ElementTopology* topology, const std::string &gold_top, const int parameteric_dim,
        const int num_vertices, const int num_nodes, const int num_edges, const int num_faces, const int num_boundaries);

TEST(Iofx, test_constructor)
{
    Ioss::Init::Initializer init_db;

    const std::string filename = "ADeDA.e";
    Ioss::DatabaseUsage db_usage = Ioss::READ_MODEL;
    MPI_Comm communicator = MPI_COMM_WORLD;
    Ioss::PropertyManager properties;

    properties.add(Ioss::Property("INTEGER_SIZE_DB",  8));
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));

    Iofx::DatabaseIO *db_io = new Iofx::DatabaseIO(NULL, filename, db_usage, communicator, properties);

    // Creation of region reads meta data using db_io
    Ioss::Region region(db_io);

    EXPECT_TRUE(db_io->ok());

    int64_t id = 1;
    std::vector<std::string> block_membership;
    db_io->compute_block_membership(id, block_membership);

    ASSERT_EQ(2u, block_membership.size());
    EXPECT_EQ("block_1", block_membership[0]);
    EXPECT_EQ("block_2", block_membership[1]);

    // title comparison did not work
    //    std::string title()               const     {return databaseTitle;}

    int spatial_dim = db_io->spatial_dimension();
    EXPECT_EQ(3, spatial_dim);

    int64_t num_nodes = db_io->node_count();
    EXPECT_EQ(12, num_nodes);

    int64_t num_elems = db_io->element_count();
    EXPECT_EQ(3, num_elems);

    int num_node_blocks = db_io->node_block_count();
    EXPECT_EQ(1, num_node_blocks);

    int num_elem_blocks = db_io->element_block_count();
    EXPECT_EQ(2, num_elem_blocks);

    int num_sidesets = db_io->sideset_count();
    EXPECT_EQ(2, num_sidesets);

    int num_nodesets = db_io->nodeset_count();
    EXPECT_EQ(0, num_nodesets);

    const std::vector<Ioss::ElementBlock*> &element_blocks = region.get_element_blocks();
    EXPECT_EQ(2u, element_blocks.size());

    std::vector<std::string> gold_strings{"block_2", "block_1"};
    std::vector<std::string> gold_top_names{"hex8", "shell4"};
    std::vector<int> parametric_dim{3, 2};
    std::vector<int> num_vertices{8, 4};
    std::vector<int> number_nodes{8, 4};
    std::vector<int> num_edges{12, 4};
    std::vector<int> num_faces{6, 2};
    std::vector<int> num_boundaries{6, 6};
    std::vector<int> gold_conn_size{16, 4};
    std::vector<std::vector<int> > gold_connectivity{{1,2,3,4,5,6,7,8,5,6,7,8,9,10,11,12}, {5,6,7,8}};

    std::vector<size_t> gold_num_elements_per_block{2, 1};
    std::vector<std::vector<int64_t> > gold_ids{{1,2},{3}};

    std::vector<bool> attributes_exist{false, true};
    std::vector<int> num_attributes{0, 1};

    // element block testing

    for (size_t i=0;i<element_blocks.size();++i)
    {
        std::vector<std::string> empty;
        db_io->get_block_adjacencies(element_blocks[i], empty);
        ASSERT_TRUE(!empty.empty());
        EXPECT_EQ(gold_strings[i], empty[0]);

        const Ioss::ElementTopology *topology = element_blocks[i]->topology();
        test_topology(topology, gold_top_names[i], parametric_dim[i], num_vertices[i], number_nodes[i],
                num_edges[i], num_faces[i], num_boundaries[i]);

        std::vector<int64_t> connectivity;
        element_blocks[i]->get_field_data("connectivity_raw", connectivity);

        EXPECT_EQ(gold_conn_size[i], connectivity.size());
        for (int j=0;j<gold_conn_size[i];++j)
        {
            EXPECT_EQ(gold_connectivity[i][j], connectivity[j]);
        }

        std::vector<int64_t>  ids;
        element_blocks[i]->get_field_data("ids", ids);

        EXPECT_EQ(gold_num_elements_per_block[i], ids.size());
        for (size_t j=0;j<ids.size();++j)
        {
            EXPECT_EQ(gold_ids[i][j], ids[j]);
        }

        std::vector<double> attributeValues;
        EXPECT_EQ(attributes_exist[i], element_blocks[i]->field_exists("attribute"));
        if (attributes_exist[i])
        {
            element_blocks[i]->get_field_data("attribute", attributeValues);
            EXPECT_EQ(num_attributes[i], attributeValues.size()) << gold_top_names[i];
        }
    }

    // node block testing

    std::vector<double> gold_coordinates{-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 1.5,
        -0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5,
         0.5, -0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5 };

    Ioss::NodeBlock *nb = region.get_node_blocks()[0];
    std::vector<double> coordinates;
    nb->get_field_data("mesh_model_coordinates", coordinates);
    size_t num_coordinates = num_nodes*spatial_dim;

    ASSERT_TRUE(coordinates.size() == num_coordinates);
    for (int i=0;i<num_nodes;++i)
    {
        for (int j=0;j<spatial_dim;++j)
        {
            int index_coordinate = spatial_dim*i+j;
            int index_gold_coordinate = num_nodes*j +i;
            // gold_coordinates = { all_x, all_y, all_z };
            // coordinates = { x1, y1, z1, x2, y2, z2, ..., };
            EXPECT_EQ(gold_coordinates[index_gold_coordinate], coordinates[index_coordinate]);
        }
    }

    // sidesets

    const std::vector<Ioss::SideSet*> &sidesets = region.get_sidesets();
    EXPECT_EQ(2u, sidesets.size());

    // std::vector<bool> gold_df_exist{true, true};

    // sidesets might be sorted by local element id

    std::vector<std::vector<int64_t> > gold_element_ids = { {1, 3}, {2, 3} };
    std::vector<std::vector<int64_t> > gold_side_ids = { {6, 2}, {5, 1} };

    for (size_t i=0;i<sidesets.size();++i)
    {
        std::vector<int64_t> element_ids;
        std::vector<int64_t> side_ids;

        for (size_t j=0;j<sidesets[i]->block_count();++j)
        {
            Ioss::SideBlock *block = sidesets[i]->get_block(j);
            std::vector<int64_t> face_ids_per_block;
            std::vector<int64_t> side_ids_per_block;

            ASSERT_TRUE(block->field_exists("ids"));
            block->get_field_data("ids", face_ids_per_block);

            ASSERT_TRUE(block->field_exists("element_side"));
            block->get_field_data("element_side", side_ids_per_block);
            for (size_t k=0;k<side_ids_per_block.size();k+=2)
            {
                element_ids.push_back(side_ids_per_block[k]);
                side_ids.push_back(side_ids_per_block[k+1]);
            }
        }

        ASSERT_EQ(element_ids.size(), side_ids.size());
        for (size_t j=0;j<element_ids.size();++j)
        {
            EXPECT_EQ(gold_element_ids[i][j], element_ids[j]) << sidesets[i]->name();
            EXPECT_EQ(gold_side_ids[i][j], side_ids[j]) << sidesets[i]->name();
        }

        /*
        Ioss::NameList names;
        sidesets[i]->field_describe( &names );
        for (size_t j=0;j<names.size();++j)
        {
            std::cerr << "for " << i << " name is " << names[j] << std::endl;
        }

        std::string df_name = "distribution_factors_" + sidesets[i]->name();
        EXPECT_EQ(gold_df_exist[i], sidesets[i]->field_exists(df_name));
        if (gold_df_exist[i])
        {
            std::vector<double> df;
            sidesets[i]->get_field_data(df_name, df);
            ASSERT_EQ(8, df.size());
            for (size_t j=0;j<df.size();++j)
            {
                EXPECT_EQ(1, df.size());
            }
        }
        */
    }
}

void test_topology(const Ioss::ElementTopology* topology, const std::string &gold_top, const int parameteric_dim,
        const int num_vertices, const int num_nodes, const int num_edges, const int num_faces, const int num_boundaries)
{
    const std::string &name = topology->name();
    EXPECT_EQ(gold_top, name);

    EXPECT_TRUE(topology->is_element());
    EXPECT_EQ(3, topology->spatial_dimension()) << gold_top;
    EXPECT_EQ(parameteric_dim, topology->parametric_dimension()) << gold_top;
    EXPECT_EQ(1, topology->order()) << gold_top;

    EXPECT_TRUE(topology->edges_similar());
    EXPECT_TRUE(topology->faces_similar());

    EXPECT_EQ(num_vertices, topology->number_corner_nodes());
    EXPECT_EQ(num_nodes, topology->number_nodes());
    EXPECT_EQ(num_edges, topology->number_edges());
    EXPECT_EQ(num_faces, topology->number_faces());
    EXPECT_EQ(num_boundaries, topology->number_boundaries()) << gold_top;

    for (int i=0;i<num_edges;++i)
    {
        EXPECT_EQ(2, topology->number_nodes_edge(i));
    }

    for (int i=0;i<num_faces;++i)
    {
        EXPECT_EQ(4, topology->number_nodes_face(i));
        EXPECT_EQ(4, topology->number_edges_face(i));
    }

    std::vector<int> element_connectivity = topology->element_connectivity();
    EXPECT_EQ(num_nodes, element_connectivity.size()) << gold_top;
}


}
