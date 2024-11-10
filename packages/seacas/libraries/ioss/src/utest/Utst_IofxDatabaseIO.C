// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ionit_Initializer.h"
#include "Ioss_DBUsage.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_Hex8.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_Region.h"
#include "Ioss_Shell4.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "exodus/Ioex_DatabaseIO.h"

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <fmt/core.h>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <vector>

#include "Ioss_CodeTypes.h"
#include "Ioss_Field.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Property.h"
#include "Ioss_ScopeGuard.h"
#include "Ioss_State.h"
#include "Ioss_SurfaceSplit.h"

namespace {
  std::string input_filename = "ADeDA.e";

  void test_topology(const Ioss::ElementTopology *topology, const std::string &gold_top,
                     const int parameteric_dim, const int num_vertices, const int num_nodes,
                     const int num_edges, const int num_faces, const int num_boundaries);

  Ioex::DatabaseIO *create_input_db_io(const std::string &filename)
  {
    Ioss::Init::Initializer init_db;

    Ioss::DatabaseUsage   db_usage = Ioss::READ_MODEL;
    Ioss::PropertyManager properties;

    properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));

    auto *db_io = new Ioex::DatabaseIO(nullptr, filename, db_usage,
                                       Ioss::ParallelUtils::comm_world(), properties);
    return db_io;
  }

  // BeginDocTest1
  Ioex::DatabaseIO *create_output_db_io(const std::string &filename)
  {
    Ioss::Init::Initializer init_db;
    Ioss::DatabaseUsage     db_usage = Ioss::WRITE_RESULTS;
    Ioss::PropertyManager   properties;

    properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));

    auto *db_io = new Ioex::DatabaseIO(nullptr, filename, db_usage,
                                       Ioss::ParallelUtils::comm_world(), properties);
    return db_io;
  }

  TEST_CASE("Ioex::constructor")
  {
    Ioex::DatabaseIO *db_io = create_input_db_io(input_filename);
    db_io->set_surface_split_type(Ioss::SPLIT_BY_ELEMENT_BLOCK);

    Ioss::Region region(db_io);

    CHECK(db_io->ok());

    const std::vector<Ioss::ElementBlock *> &element_blocks = region.get_element_blocks();
    CHECK(2u == element_blocks.size());

    std::vector<std::string>      gold_strings{"block_2", "block_1"};
    std::vector<std::string>      gold_top_names{Ioss::Hex8::name, Ioss::Shell4::name};
    std::vector<int>              parametric_dim{3, 2};
    std::vector<int>              num_vertices{8, 4};
    std::vector<int>              number_nodes{8, 4};
    std::vector<int>              num_edges{12, 4};
    std::vector<int>              num_faces{6, 2};
    std::vector<int>              num_boundaries{6, 6};
    std::vector<size_t>           gold_conn_size{16, 4};
    std::vector<std::vector<int>> gold_connectivity{
        {1, 2, 3, 4, 5, 6, 7, 8, 5, 6, 7, 8, 9, 10, 11, 12}, {5, 6, 7, 8}};

    std::vector<size_t>               gold_num_elements_per_block{2, 1};
    std::vector<std::vector<int64_t>> gold_ids{{1, 2}, {3}};

    std::vector<bool>   attributes_exist{false, true};
    std::vector<size_t> num_attributes{0, 1};

    // element block testing

    for (size_t i = 0; i < element_blocks.size(); ++i) {
      std::vector<std::string> empty = element_blocks[i]->get_block_adjacencies();
      CHECK(!empty.empty());
      CHECK(gold_strings[i] == empty[0]);

      const Ioss::ElementTopology *topology = element_blocks[i]->topology();
      test_topology(topology, gold_top_names[i], parametric_dim[i], num_vertices[i],
                    number_nodes[i], num_edges[i], num_faces[i], num_boundaries[i]);

      std::vector<int64_t> connectivity;
      element_blocks[i]->get_field_data("connectivity_raw", connectivity);

      CHECK(gold_conn_size[i] == connectivity.size());
      for (size_t j = 0; j < gold_conn_size[i]; ++j) {
        CHECK(gold_connectivity[i][j] == connectivity[j]);
      }

      std::vector<int64_t> ids;
      element_blocks[i]->get_field_data("ids", ids);

      CHECK(gold_num_elements_per_block[i] == ids.size());
      for (size_t j = 0; j < ids.size(); ++j) {
        CHECK(gold_ids[i][j] == ids[j]);
      }

      std::vector<double> attributeValues;
      CHECK(attributes_exist[i] == element_blocks[i]->field_exists("attribute"));
      if (attributes_exist[i]) {
        int num_attr = element_blocks[i]->get_property("attribute_count").get_int();
        CHECK(1 == num_attr);
        element_blocks[i]->get_field_data("attribute", attributeValues);
        CHECK(num_attributes[i] == attributeValues.size());
      }
    }

    // node block testing

    std::vector<double> gold_coordinates{-0.5, -0.5, -0.5, -0.5, 0.5,  0.5, 0.5,  0.5,  1.5,
                                         1.5,  1.5,  1.5,  -0.5, -0.5, 0.5, 0.5,  -0.5, -0.5,
                                         0.5,  0.5,  -0.5, -0.5, 0.5,  0.5, 0.5,  -0.5, -0.5,
                                         0.5,  0.5,  -0.5, -0.5, 0.5,  0.5, -0.5, -0.5, 0.5};

    Ioss::NodeBlock    *nb = region.get_node_blocks()[0];
    std::vector<double> coordinates;
    nb->get_field_data("mesh_model_coordinates", coordinates);
    int64_t num_nodes   = nb->entity_count();
    int64_t spatial_dim = nb->get_property("component_degree").get_int();

    size_t num_coordinates = num_nodes * spatial_dim;

    CHECK(coordinates.size() == num_coordinates);
    for (int i = 0; i < num_nodes; ++i) {
      for (int j = 0; j < spatial_dim; ++j) {
        int index_coordinate      = spatial_dim * i + j;
        int index_gold_coordinate = num_nodes * j + i;
        // gold_coordinates = { all_x, all_y, all_z };
        // coordinates = { x1, y1, z1, x2, y2, z2, ..., };
        CHECK(gold_coordinates[index_gold_coordinate] == coordinates[index_coordinate]);
      }
    }

    // sidesets

    const std::vector<Ioss::SideSet *> &sidesets = region.get_sidesets();
    CHECK(2u == sidesets.size());

    // std::vector<bool> gold_df_exist{true, true};

    // sidesets might be sorted by local element id

    // in exodus file: sideset 1 is elements 1 and 3, sideset 2 is 3 and 2 (however it is reversed
    // in results below)

    std::vector<std::vector<int64_t>> gold_element_ids  = {{1, 3}, {2, 3}};
    std::vector<std::vector<int64_t>> gold_side_ids     = {{6, 2}, {5, 1}};
    std::vector<std::vector<int64_t>> gold_sideset_conn = {{5, 6, 7, 8, 5, 8, 7, 6},
                                                           {5, 8, 7, 6, 5, 6, 7, 8}};

    for (size_t i = 0; i < sidesets.size(); ++i) {
      std::vector<int64_t> element_ids;
      std::vector<int64_t> side_ids;
      std::vector<int64_t> connectivity;

      for (size_t j = 0; j < sidesets[i]->block_count(); ++j) {
        Ioss::SideBlock     *block = sidesets[i]->get_block(j);
        std::vector<int64_t> side_ids_per_block;
        std::vector<int64_t> connectivity_per_block;

        CHECK(block->field_exists("element_side_raw"));
        block->get_field_data("element_side_raw", side_ids_per_block);
        for (size_t k = 0; k < side_ids_per_block.size(); k += 2) {
          element_ids.push_back(side_ids_per_block[k]);
          side_ids.push_back(side_ids_per_block[k + 1]);
        }

        CHECK(block->field_exists("connectivity_raw"));
        block->get_field_data("connectivity_raw", connectivity_per_block);
        CHECK(4u == connectivity_per_block.size());
        connectivity.insert(connectivity.end(), connectivity_per_block.begin(),
                            connectivity_per_block.end());
      }

      CHECK(element_ids.size() == side_ids.size());
      for (size_t j = 0; j < element_ids.size(); ++j) {
        INFO(fmt::format("Sideset Element and Side Ids check for sideset {}, element {}", i + 1,
                         j + 1));
        CHECK(gold_element_ids[i][j] == element_ids[j]);
        CHECK(gold_side_ids[i][j] == side_ids[j]);
      }

      CHECK(gold_sideset_conn[i].size() == connectivity.size());
      for (size_t j = 0; j < connectivity.size(); ++j) {
        INFO(fmt::format("Sideset Connectivity for sideset {}, node {}", i + 1, j + 1));
        CHECK(gold_sideset_conn[i][j] == connectivity[j]);
      }
    }
  }
  // EndDocTest1

  void test_topology(const Ioss::ElementTopology *topology, const std::string &gold_top,
                     const int parameteric_dim, const int num_vertices, const int num_nodes,
                     const int num_edges, const int num_faces, const int num_boundaries)
  {
    const std::string &name = topology->name();
    CHECK(gold_top == name);

    CHECK(topology->is_element());
    CHECK(3 == topology->spatial_dimension());
    CHECK(parameteric_dim == topology->parametric_dimension());
    CHECK(1 == topology->order());

    CHECK(topology->edges_similar());
    CHECK(topology->faces_similar());

    CHECK(num_vertices == topology->number_corner_nodes());
    CHECK(num_nodes == topology->number_nodes());
    CHECK(num_edges == topology->number_edges());
    CHECK(num_faces == topology->number_faces());
    CHECK(num_boundaries == topology->number_boundaries());

    for (int i = 0; i < num_edges; ++i) {
      CHECK(2 == topology->number_nodes_edge(i));
    }

    for (int i = 0; i < num_faces; ++i) {
      CHECK(4 == topology->number_nodes_face(i));
      CHECK(4 == topology->number_edges_face(i));
    }

    std::vector<int> element_connectivity = topology->element_connectivity();
    CHECK((size_t)num_nodes == element_connectivity.size());
  }

  // BeginDocTest2
  TEST_CASE("Ioex::write_file")
  {
    Ioex::DatabaseIO *db_in = create_input_db_io(input_filename);

    Ioss::Region input_region(db_in);

    CHECK(db_in->ok());

    const std::string output_filename = "ADeDA_out.e";

    Ioex::DatabaseIO *db_out = create_output_db_io(output_filename);
    Ioss::Region      output_region(db_out);
    CHECK(db_out->ok());

    output_region.begin_mode(Ioss::STATE_DEFINE_MODEL);

    std::string NodeBlockName = "nodeblock_1";

    auto    nb          = input_region.get_node_blocks()[0];
    int64_t num_nodes   = nb->entity_count();
    int     spatial_dim = nb->get_property("component_degree").get_int();

    auto *output_node_block = new Ioss::NodeBlock(db_out, NodeBlockName, num_nodes, spatial_dim);
    output_region.add(output_node_block);

    const auto &input_element_blocks = input_region.get_element_blocks();
    for (const auto &block : input_element_blocks) {
      const auto *topology = block->topology();
      int         block_id = block->get_property("id").get_int();

      const std::string &name         = block->name();
      std::string        exotype      = topology->name();
      int64_t            num_elements = block->entity_count();

      auto *output_element_block = new Ioss::ElementBlock(db_out, name, exotype, num_elements);
      output_element_block->property_add(Ioss::Property("original_topology_type", exotype));
      output_element_block->property_add(Ioss::Property("id", block_id));
      output_region.add(output_element_block);

      int num_attributes = block->get_property("attribute_count").get_int();
      for (int j = 0; j < num_attributes; ++j) {
        output_element_block->field_add(Ioss::Field("attribute", Ioss::Field::REAL, "scalar",
                                                    Ioss::Field::ATTRIBUTE, num_elements, j + 1));
      }
    }

    const auto &nodesets_input = input_region.get_nodesets();
    for (const auto &nset : nodesets_input) {
      const std::string &nodeset_name            = nset->name();
      int64_t            number_nodes_in_nodeset = nset->entity_count();

      auto *const nodeset = new Ioss::NodeSet(db_out, nodeset_name, number_nodes_in_nodeset);
      output_region.add(nodeset);
      int num_local_owned_nodes = number_nodes_in_nodeset;
      nodeset->property_add(Ioss::Property("locally_owned_count", num_local_owned_nodes));
    }

    const auto &sidesets_input = input_region.get_sidesets();
    for (const auto &sset : sidesets_input) {
      const std::string &sideset_name = sset->name();

      auto *const sideset_output = new Ioss::SideSet(db_out, sideset_name);
      output_region.add(sideset_output);

      const auto &side_blocks = sset->get_side_blocks();
      for (const auto &sblock : side_blocks) {
        const std::string &topo_name        = sblock->topology()->name();
        int64_t            side_count       = sblock->entity_count();
        const std::string &parent_topo_name = sblock->parent_element_topology()->name();
        const std::string &side_block_name  = sblock->name();
        auto              *side_block =
            new Ioss::SideBlock(db_out, side_block_name, topo_name, parent_topo_name, side_count);

        sideset_output->add(side_block);

        int         nodes_per_side = sblock->topology()->number_nodes();
        std::string storage_type   = "Real[";
        storage_type += std::to_string(nodes_per_side);
        storage_type += "]";
        side_block->field_add(Ioss::Field("distribution_factors", Ioss::Field::REAL, storage_type,
                                          Ioss::Field::MESH, side_count));
      }
    }

    output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

    output_region.begin_mode(Ioss::STATE_MODEL);

    {
      auto *input_node_block = input_region.get_node_blocks()[0];

      std::vector<double> coordinates;
      input_node_block->get_field_data("mesh_model_coordinates", coordinates);
      output_node_block->put_field_data("mesh_model_coordinates", coordinates);

      std::vector<int64_t> node_ids;
      input_node_block->get_field_data("ids", node_ids);
      output_node_block->put_field_data("ids", node_ids);
    }

    auto output_element_blocks = output_region.get_element_blocks();
    for (size_t blk = 0; blk < output_element_blocks.size(); blk++) {
      std::vector<int64_t> elem_ids;
      input_element_blocks[blk]->get_field_data("ids", elem_ids);
      output_element_blocks[blk]->put_field_data("ids", elem_ids);

      std::vector<int64_t> connectivity;
      input_element_blocks[blk]->get_field_data("connectivity_raw", connectivity);
      output_element_blocks[blk]->put_field_data("connectivity_raw", connectivity);

      if (input_element_blocks[blk]->field_exists("attribute")) {
        std::vector<double> attributeValues;
        input_element_blocks[blk]->get_field_data("attribute", attributeValues);
        output_element_blocks[blk]->put_field_data("attribute", attributeValues);
      }
    }

    const auto &nodesets_output = output_region.get_nodesets();
    for (size_t i = 0; i < nodesets_output.size(); ++i) {
      std::vector<int64_t> node_ids;
      nodesets_input[i]->get_field_data("ids", node_ids);
      nodesets_output[i]->put_field_data("ids", node_ids);

      std::vector<double> dist_factors;
      nodesets_input[i]->get_field_data("distribution_factors", dist_factors);
      nodesets_output[i]->put_field_data("distribution_factors", dist_factors);
    }

    const auto &sidesets_output = output_region.get_sidesets();
    for (size_t i = 0; i < sidesets_output.size(); ++i) {
      const auto &side_blocks_input  = sidesets_input[i]->get_side_blocks();
      const auto &side_blocks_output = sidesets_output[i]->get_side_blocks();
      for (size_t k = 0; k < side_blocks_output.size(); ++k) {
        std::vector<int64_t> elem_side_ids;
        side_blocks_input[k]->get_field_data("element_side", elem_side_ids);
        side_blocks_output[k]->put_field_data("element_side", elem_side_ids);

        std::vector<double> dist_factors;
        side_blocks_input[k]->get_field_data("distribution_factors", dist_factors);
        side_blocks_output[k]->put_field_data("distribution_factors", dist_factors);
      }
    }

    output_region.end_mode(Ioss::STATE_MODEL);
  }
  // EndDocTest2
} // namespace
int main(IOSS_MAYBE_UNUSED int argc, char **argv)
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(MPI_Finalize);
#endif

  Catch::Session session; // There must be exactly one instance

  // Build a new parser on top of Catch2's
  using namespace Catch::Clara;
  auto cli = session.cli()                     // Get Catch2's command line parser
             | Opt(input_filename, "filename") // bind variable to a new option, with a hint string
                   ["-f"]["--filename"]        // the option names it will respond to
             ("Filename to read mesh from.");  // description string for the help output

  // Now pass the new composite back to Catch2 so it uses that
  session.cli(cli);

  // Let Catch2 (using Clara) parse the command line
  int returnCode = session.applyCommandLine(argc, argv);
  if (returnCode != 0) // Indicates a command line error
    return returnCode;

  fmt::print("'{}'\n", input_filename);
  return session.run();
}
