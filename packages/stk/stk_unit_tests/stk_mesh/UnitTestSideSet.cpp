#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/InputFile.hpp>
#include <stk_io/SidesetUpdater.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

namespace {
void move_elems_from_block_to_block(stk::mesh::BulkData& bulk,
                                    const std::vector<stk::mesh::EntityId>& elemIDs,
                                    const std::string& fromBlockName,
                                    const std::string& toBlockName)
{
  stk::mesh::Part& fromBlock = *bulk.mesh_meta_data().get_part(fromBlockName);
  stk::mesh::Part& toBlock = *bulk.mesh_meta_data().get_part(toBlockName);

  stk::mesh::EntityVector elems;
  for(stk::mesh::EntityId elemID : elemIDs) {
    stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, elemID);
    ThrowRequireMsg(bulk.is_valid(elem), "Failed to find element with ID="<<elemID);
    elems.push_back(elem);
  }

  bulk.batch_change_entity_parts(elems, stk::mesh::PartVector{&toBlock}, stk::mesh::PartVector{&fromBlock});
}

void create_sides_between_blocks(stk::mesh::BulkData& bulk,
                                 const std::string& block1Name,
                                 const std::string& block2Name,
                                 const std::string& sidePartName)
{
  stk::mesh::Part& block1 = *bulk.mesh_meta_data().get_part(block1Name);
  stk::mesh::Part& block2 = *bulk.mesh_meta_data().get_part(block2Name);
  stk::mesh::Part& sidePart = *bulk.mesh_meta_data().get_part(sidePartName);

  stk::mesh::Selector blockSelector = block1 | block2;
  stk::mesh::create_interior_block_boundary_sides(bulk, blockSelector, stk::mesh::PartVector{&sidePart});
}

void move_sides_between_blocks_into_sideset_part(stk::mesh::BulkData& bulk,
                                 const std::string& block1Name,
                                 const std::string& block2Name,
                                 const std::string& sidePartName)
{
  stk::mesh::Part& block1 = *bulk.mesh_meta_data().get_part(block1Name);
  stk::mesh::Part& block2 = *bulk.mesh_meta_data().get_part(block2Name);
  stk::mesh::Part& sidePart = *bulk.mesh_meta_data().get_part(sidePartName);

  stk::mesh::Selector betweenBlocks = block1 & block2;
  stk::mesh::EntityVector sidesBetweenBlocks;
  stk::mesh::get_selected_entities(betweenBlocks, bulk.buckets(stk::topology::FACE_RANK), sidesBetweenBlocks);

  bulk.batch_change_entity_parts(sidesBetweenBlocks, stk::mesh::PartVector{&sidePart}, stk::mesh::PartVector{});
}

void create_sideset(stk::mesh::BulkData& bulk,
                    const std::string& surfacePartName,
                    const std::string& blockPartName)
{
  stk::mesh::Part& blockPart = *bulk.mesh_meta_data().get_part(blockPartName);
  stk::mesh::Part& surfacePart = *bulk.mesh_meta_data().get_part(surfacePartName);

  bulk.mesh_meta_data().set_surface_to_block_mapping(&surfacePart, stk::mesh::ConstPartVector{&blockPart});
  bulk.create_sideset(surfacePart);
}

void check_sideset_elems(const stk::mesh::BulkData& bulk,
                         const std::string surfacePartName,
                         const std::vector<stk::mesh::EntityId>& elemIDs)
{
  EXPECT_EQ(1u, bulk.get_number_of_sidesets());

  stk::mesh::Part& surfacePart = *bulk.mesh_meta_data().get_part(surfacePartName);
  const stk::mesh::SideSet& sideset = bulk.get_sideset(surfacePart);
  EXPECT_EQ(4u, sideset.size());

  for(size_t i=0; i<elemIDs.size(); ++i) {
    EXPECT_EQ(elemIDs[i], bulk.identifier(sideset[i].element));
  }
}

void load_mesh_with_specified_splitting(stk::mesh::BulkData &bulk, const std::string &filename, Ioss::SurfaceSplitType splitType)
{
  stk::io::StkMeshIoBroker stkIo;
  stkIo.set_bulk_data(bulk);
  size_t inputIdx = stkIo.add_mesh_database(filename, stk::io::READ_MESH);
  stk::io::InputFile &inputFile = stkIo.get_mesh_database(inputIdx);
  inputFile.set_surface_split_type(splitType);

  stkIo.create_input_mesh();
  stkIo.add_all_mesh_fields_as_input_fields();
  stkIo.populate_bulk_data();
}

void write_mesh_with_specified_splitting(stk::mesh::BulkData &bulk, const std::string &filename, Ioss::SurfaceSplitType splitType)
{
  stk::io::StkMeshIoBroker stkIo;
  stkIo.set_bulk_data(bulk);

  size_t outputIdx = stkIo.create_output_mesh(filename, stk::io::WRITE_RESULTS);
  Ioss::Region* ioRegion = stkIo.get_output_io_region(outputIdx).get();

  Ioss::DatabaseIO *db = ioRegion->get_database();
  assert(db != nullptr);
  db->set_surface_split_type(splitType);

  stkIo.write_output_mesh(outputIdx);
}

void create_two_elem_block_mesh_with_spanning_sidesets(const std::string& filename)
{
  const int dim = 3;
  stk::mesh::MetaData meta(dim);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);


  const std::string meshDesc = "generated:1x1x2|sideset:xX";

  stk::mesh::Part& block_2 = meta.declare_part("block_2", stk::topology::ELEM_RANK);
  stk::mesh::set_topology(block_2, stk::topology::HEX_8);
  meta.set_part_id(block_2, 2);
  stk::io::put_io_part_attribute(block_2);

  load_mesh_with_specified_splitting(bulk, meshDesc, Ioss::SPLIT_BY_TOPOLOGIES);

  stk::mesh::Part* block_1 = meta.get_part("block_1");
  EXPECT_TRUE(block_1 != nullptr);

  stk::mesh::Part* surface_1 = meta.get_part("surface_1");
  stk::mesh::Part* surface_2 = meta.get_part("surface_2");
  stk::mesh::Part* surface_1_quad4 = meta.get_part("surface_1_quad4");
  stk::mesh::Part* surface_2_quad4 = meta.get_part("surface_2_quad4");

  EXPECT_TRUE(surface_1 != nullptr);
  EXPECT_TRUE(surface_2 != nullptr);
  EXPECT_TRUE(surface_1_quad4 != nullptr);
  EXPECT_TRUE(surface_2_quad4 != nullptr);

  std::vector<const stk::mesh::Part*> touchingParts{block_1, &block_2};

  meta.set_surface_to_block_mapping(surface_1, touchingParts);
  meta.set_surface_to_block_mapping(surface_2, touchingParts);
  meta.set_surface_to_block_mapping(surface_1_quad4, touchingParts);
  meta.set_surface_to_block_mapping(surface_2_quad4, touchingParts);

  bulk.modification_begin("");
  stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, 2);
  if(bulk.is_valid(elem2) && bulk.bucket(elem2).owned()) {
    stk::mesh::PartVector addParts    = {&block_2};
    stk::mesh::PartVector removeParts = {block_1};

    bulk.change_entity_parts(elem2, addParts, removeParts);
  }
  bulk.modification_end();

  write_mesh_with_specified_splitting(bulk, filename, Ioss::SPLIT_BY_ELEMENT_BLOCK);
}
}

class TestSideSet : public stk::unit_test_util::MeshFixture
{
protected:
  void setup_2_block_mesh()
  {
    stk::mesh::MetaData& meta = get_meta();
    stk::mesh::Part& block2 = meta.declare_part_with_topology("block_2", stk::topology::HEX_8);
    stk::mesh::Part& surface1 = meta.declare_part_with_topology("surface_1", stk::topology::QUAD_4);
    stk::io::put_io_part_attribute(block2);
    stk::io::put_io_part_attribute(surface1);
    meta.set_part_id(block2, 2);
    meta.set_part_id(surface1, 1);

    allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::io::fill_mesh("generated:2x2x2", get_bulk());

    move_elems_from_block_to_block(get_bulk(), {2, 4, 6, 8}, "block_1", "block_2");
  }
};

TEST_F(TestSideSet, creatingSideOfOneElem_eachProcHasOneSide)
{
  setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::EntityVector localElems;
  stk::mesh::get_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEM_RANK), localElems);
  ASSERT_GE(localElems.size(), 1u);

  stk::mesh::SideSet sideSet(get_bulk());
  sideSet.add(stk::mesh::SideSetEntry(localElems[0], 1));

  EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(get_meta().side_rank())));
  get_bulk().create_side_entities(sideSet, {});
  EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(get_meta().side_rank())));
}

TEST_F(TestSideSet, createSideSetsSpanningMultipleBlocks)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) < 3) {
    const std::string filename = "sideSetMesh.g";
    create_two_elem_block_mesh_with_spanning_sidesets(filename);

    allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);

    load_mesh_with_specified_splitting(*bulkData, filename, Ioss::SPLIT_BY_ELEMENT_BLOCK);

    stk::io::create_bulkdata_sidesets(get_bulk());

    EXPECT_EQ(2u, get_bulk().get_number_of_sidesets());
    stk::unit_test_util::delete_mesh(filename);
  }
}


TEST_F(TestSideSet, createSingleSidedSideSetOnBlock1BetweenBlocks)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    setup_2_block_mesh();

    create_sideset(get_bulk(), "surface_1", "block_1");
    create_sides_between_blocks(get_bulk(), "block_1", "block_2", "surface_1");

    check_sideset_elems(get_bulk(), "surface_1", {1, 3, 5, 7});
  }
}

TEST_F(TestSideSet, createSingleSidedSideSetOnBlock1BetweenBlocks_ChangeSideParts)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    setup_2_block_mesh();

    create_sideset(get_bulk(), "surface_1", "block_1");
    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), stk::mesh::PartVector{}, false);
    move_sides_between_blocks_into_sideset_part(get_bulk(), "block_1", "block_2", "surface_1");

    check_sideset_elems(get_bulk(), "surface_1", {1, 3, 5, 7});
  }
}

TEST_F(TestSideSet, createSingleSidedSideSetOnBlock2BetweenBlocks)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    setup_2_block_mesh();

    create_sideset(get_bulk(), "surface_1", "block_2");
    create_sides_between_blocks(get_bulk(), "block_1", "block_2", "surface_1");

    check_sideset_elems(get_bulk(), "surface_1", {2, 4, 6, 8});
  }
}

TEST_F(TestSideSet, createSingleSidedSideSetOnBlock2BetweenBlocks_ChangeSideParts)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    setup_2_block_mesh();

    create_sideset(get_bulk(), "surface_1", "block_2");
    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), stk::mesh::PartVector{}, false);
    move_sides_between_blocks_into_sideset_part(get_bulk(), "block_1", "block_2", "surface_1");

    check_sideset_elems(get_bulk(), "surface_1", {2, 4, 6, 8});
  }
}

TEST_F(TestSideSet, findSidesetEntry)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    setup_mesh("generated:1x1x1|sideset:x", stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Part* sidePart = get_meta().get_part("surface_1");

    stk::mesh::SideSet& sideset = get_bulk().get_sideset(*sidePart);
    stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
    const stk::mesh::ConnectivityOrdinal expectedSideOrdinal = 3;
    const int numSides = 6;
    for (int i = 0; i < numSides; ++i)
    {
      if (expectedSideOrdinal != stk::mesh::ConnectivityOrdinal(i))
      {
        EXPECT_FALSE(sideset.contains(elem1, stk::mesh::ConnectivityOrdinal(i)));
        EXPECT_FALSE(sideset.contains(stk::mesh::SideSetEntry(elem1, stk::mesh::ConnectivityOrdinal(i))));
      }
      else
      {
        EXPECT_TRUE(sideset.contains(elem1, stk::mesh::ConnectivityOrdinal(i)));
        EXPECT_TRUE(sideset.contains(stk::mesh::SideSetEntry(elem1, stk::mesh::ConnectivityOrdinal(i))));
      }
    }
  }
}
