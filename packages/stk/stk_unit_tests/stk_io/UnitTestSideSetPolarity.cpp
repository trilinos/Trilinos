#include <gtest/gtest.h>
#include <stk_unit_test_utils/ioUtils.hpp>
#include "stk_mesh/base/GetEntities.hpp"
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_io/StkIoUtils.hpp>
#include "stk_unit_test_utils/ReadWriteSidesetTester.hpp"
#include "stk_unit_test_utils/FaceTestingUtils.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"

namespace {

enum SidesetDirection {LEFT = 0, RIGHT, DOUBLE};

stk::mesh::PartVector create_sideset_parts(stk::mesh::MetaData &meta, const std::vector<std::string>&names)
{
  stk::mesh::PartVector parts;

  int id = 1;
  for(const std::string &name : names)
  {
      stk::mesh::Part &part = meta.declare_part(name, meta.side_rank());
      meta.set_part_id(part, id);
      stk::io::put_io_part_attribute(part);
      parts.push_back(&part);
      ++id;
  }

  return parts;
}

void create_AA_mesh(stk::mesh::BulkData &bulk)
{
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5, 6, 7, 8,block_1\n\
                            0,2,HEX_8,5,6,7,8,9,10,11,12,block_1";

    std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
                                        0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                        0,0,2, 1,0,2, 1,1,2, 0,1,2 };

    bulk.initialize_face_adjacent_element_graph();
    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, bulk);
}

void create_AB_mesh(stk::mesh::BulkData &bulk)
{
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5, 6, 7, 8,block_1\n\
                            0,2,HEX_8,5,6,7,8,9,10,11,12,block_2";

    std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
                                        0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                        0,0,2, 1,0,2, 1,1,2, 0,1,2 };

    bulk.initialize_face_adjacent_element_graph();
    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, bulk);
}

stk::mesh::Entity verify_and_get_face(stk::mesh::BulkData &bulk, stk::mesh::Part* ssPart)
{
    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(*ssPart, bulk.buckets(stk::topology::FACE_RANK), faces);
    EXPECT_EQ(1u, faces.size());
    return faces[0];
}

void populate_elem_sides(SidesetDirection direction,
                         stk::mesh::EntityIdVector &elem,
                         std::vector<int> &ordinal)
{
    if(direction == LEFT)
    {
        elem.push_back(1u);
        ordinal.push_back(5);
    }
    else if(direction == RIGHT)
    {
        elem.push_back(2u);
        ordinal.push_back(4);
    }
    else
    {
        elem.push_back(1u);
        ordinal.push_back(5);

        elem.push_back(2u);
        ordinal.push_back(4);
    }
}

void populate_sideset_names(SidesetDirection direction,
                            std::vector<std::string> &names)
{
    if(direction == LEFT)
    {
        names = {"surface_1", "surface_block_1_QUAD4_1"};
    }
    else if(direction == RIGHT)
    {
        names = {"surface_1", "surface_block_2_QUAD4_1"};
    }
    else
    {
        names = {"surface_1", "surface_block_1_QUAD4_1", "surface_block_2_QUAD4_1"};
    }
}

stk::mesh::Part* create_AA_mesh_with_sideset(stk::mesh::BulkData &bulk,
                                             SidesetDirection direction)
{
    stk::mesh::EntityIdVector elem;
    std::vector<int> ordinal;

    stk::mesh::MetaData &meta = bulk.mesh_meta_data();
    stk::mesh::PartVector parts = create_sideset_parts(meta, std::vector<std::string>{"surface_1"});

    create_AA_mesh(bulk);

    populate_elem_sides(direction, elem, ordinal);
    ThrowRequire(elem.size() == ordinal.size());

    stk::mesh::SideSet &sideSet = bulk.create_sideset(*parts[0]);

    for(unsigned i=0; i<elem.size(); ++i) {
        sideSet.add({bulk.get_entity(stk::topology::ELEMENT_RANK, elem[i]), ordinal[i]});
    }

    bulk.create_side_entities(sideSet, parts);

    stk::mesh::Part* block_1 = meta.get_part("block_1");
    EXPECT_TRUE(block_1 != nullptr);
    meta.set_part_id(*block_1, 1);

    std::vector<const stk::mesh::Part*> touchingParts{block_1};

    meta.set_surface_to_block_mapping(parts[0], touchingParts);

    return parts[0];
}

stk::mesh::Part* create_AB_mesh_with_sideset(stk::mesh::BulkData &bulk,
                                             SidesetDirection direction)
{
    stk::mesh::EntityIdVector elem;
    std::vector<int> ordinal;

    stk::mesh::MetaData &meta = bulk.mesh_meta_data();
    std::vector<std::string> sideSetNames;
    populate_sideset_names(direction, sideSetNames);

    stk::mesh::PartVector parts = create_sideset_parts(meta, sideSetNames);

    for(unsigned i = 0; i<parts.size(); ++i) {
        meta.set_part_id(*parts[i], 1);
    }

    for(unsigned i=1; i<parts.size(); ++i) {
        meta.declare_part_subset(*parts[0], *parts[i]);
    }

    create_AB_mesh(bulk);

    populate_elem_sides(direction, elem, ordinal);
    ThrowRequire(elem.size() == ordinal.size());

    stk::mesh::SideSet &sideSet = bulk.create_sideset(*parts[0]);

    for(unsigned i=0; i<elem.size(); ++i) {
        sideSet.add({bulk.get_entity(stk::topology::ELEMENT_RANK, elem[i]), ordinal[i]});
    }

    bulk.create_side_entities(sideSet, parts);

    stk::mesh::Part* block_1 = meta.get_part("block_1");
    EXPECT_TRUE(block_1 != nullptr);
    stk::mesh::Part* block_2 = meta.get_part("block_2");
    EXPECT_TRUE(block_2 != nullptr);

    meta.set_part_id(*block_1, 1);
    meta.set_part_id(*block_2, 2);

    std::vector<const stk::mesh::Part*> touchingParts;

    if(direction == LEFT)
    {
        touchingParts = {block_1};
        meta.set_surface_to_block_mapping(parts[1], touchingParts);
    }
    else if(direction == RIGHT)
    {
        touchingParts = {block_2};
        meta.set_surface_to_block_mapping(parts[1], touchingParts);
    }
    else
    {
        touchingParts = {block_1};
        meta.set_surface_to_block_mapping(parts[1], touchingParts);

        touchingParts = {block_2};
        meta.set_surface_to_block_mapping(parts[2], touchingParts);
    }

    return parts[0];
}

void check_polarity_of_modified_AA(stk::ParallelMachine pm,
                                   SidesetDirection direction,
                                   const std::string& name,
                                   bool isPositivePolarity)
{
    size_t spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData bulk(meta, pm);

    stk::mesh::Part* ssPart = create_AA_mesh_with_sideset(bulk, direction);
    stk::mesh::Entity  face = verify_and_get_face(bulk, ssPart);

    std::pair<bool,bool> sidesetPolarity = stk::io::is_positive_sideset_polarity(bulk, *ssPart, face);
    EXPECT_TRUE(sidesetPolarity.first);
    EXPECT_EQ(isPositivePolarity, sidesetPolarity.second);
}

void check_polarity_of_modified_AB(stk::ParallelMachine pm,
                                   SidesetDirection direction,
                                   const std::string& name,
                                   bool isPositivePolarity)
{
    size_t spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData bulk(meta, pm);

    stk::mesh::Part* ssPart = create_AB_mesh_with_sideset(bulk, direction);
    stk::mesh::Entity  face = verify_and_get_face(bulk, ssPart);

    std::pair<bool,bool> sidesetPolarity = stk::io::is_positive_sideset_polarity(bulk, *ssPart, face);
    EXPECT_TRUE(sidesetPolarity.first);
    EXPECT_EQ(isPositivePolarity, sidesetPolarity.second);

    std::string file = name+".e";
    stk::io::write_mesh(file, bulk);

    sidesetPolarity = stk::io::is_positive_sideset_polarity(bulk, *ssPart, face);
    EXPECT_TRUE(sidesetPolarity.first);
    EXPECT_EQ(isPositivePolarity, sidesetPolarity.second);
    unlink(file.c_str());
}

TEST(StkIo, sideset_polarity_AA)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned p_size = stk::parallel_machine_size(pm);

  if(p_size == 1)
  {
      check_polarity_of_modified_AA(pm, RIGHT, "ARA", false);
      check_polarity_of_modified_AA(pm, LEFT , "ALA", true);
  }
}

TEST(StkIo, sideset_polarity_AB)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned p_size = stk::parallel_machine_size(pm);

  if(p_size == 1)
  {
      check_polarity_of_modified_AB(pm, RIGHT, "ARB", false);
      check_polarity_of_modified_AB(pm, LEFT , "ALB", true);
  }
}

}
