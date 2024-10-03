#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/SideSetUtil.hpp>
#include <stk_mesh/base/SidesetUpdater.hpp>
#include <stk_mesh/base/DestroyElements.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/InputFile.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/SideSetUtil.hpp>
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp"

namespace {
using stk::unit_test_util::build_mesh;

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
    STK_ThrowRequireMsg(bulk.is_valid(elem), "Failed to find element with ID="<<elemID);
    elems.push_back(elem);
  }

  bulk.batch_change_entity_parts(elems, stk::mesh::PartVector{&toBlock}, stk::mesh::PartVector{&fromBlock});
}

void copy_elems_from_block_to_block(stk::mesh::BulkData& bulk,
                                    const std::vector<stk::mesh::EntityId>& elemIDs,
                                    const std::string& fromBlockName,
                                    const std::string& toBlockName)
{
  stk::mesh::Part& toBlock = *bulk.mesh_meta_data().get_part(toBlockName);

  stk::mesh::EntityVector elems;
  for(stk::mesh::EntityId elemID : elemIDs) {
    stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, elemID);
    STK_ThrowRequireMsg(bulk.is_valid(elem), "Failed to find element with ID="<<elemID);
    elems.push_back(elem);
  }

  bulk.batch_change_entity_parts(elems, stk::mesh::PartVector{&toBlock}, stk::mesh::PartVector{});
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
  bulk.create_sideset(stk::mesh::get_sideset_parent(surfacePart));
}

void check_sideset_elems(const stk::mesh::BulkData& bulk,
                         const std::string surfacePartName,
                         const std::vector<stk::mesh::EntityId>& elemIDs)
{
  EXPECT_EQ(1u, bulk.get_number_of_sidesets());

  stk::mesh::Part& surfacePart = *bulk.mesh_meta_data().get_part(surfacePartName);
  const stk::mesh::SideSet& sideset = bulk.get_sideset(surfacePart);
  ASSERT_EQ(elemIDs.size(), sideset.size());

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
  Ioss::Region* ioRegion = stkIo.get_output_ioss_region(outputIdx).get();

  Ioss::DatabaseIO *db = ioRegion->get_database();
  assert(db != nullptr);
  db->set_surface_split_type(splitType);

  stkIo.write_output_mesh(outputIdx);
}

void create_two_elem_block_mesh_with_spanning_sidesets(const std::string& filename)
{
  const int dim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(dim, MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk = *bulkPtr;


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

TEST(SkinBoundary, check_interior_shared_shell_side)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                      .set_spatial_dimension(3)
                                                      .set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA)
                                                      .create();
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::Part& boundaryPart = meta.declare_part("boundaryPart", meta.side_rank());
  std::string meshDesc = "1,10098,TET_4, 1,2,3,4, block_1\n"
                         "0,320234,SHELL_TRI_3, 1,3,4, block_2|sideset:name=surface_1; data=10098,3";

  std::vector<double> coords = {0,0,0,  0,1,0,  1,0,0, 1,1,1};

  stk::unit_test_util::setup_text_mesh(*bulkPtr, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

  stk::mesh::Entity sharedSide = bulkPtr->get_entity(meta.side_rank(), 100983);
  ASSERT_TRUE(bulkPtr->is_valid(sharedSide));
  ASSERT_TRUE(bulkPtr->bucket(sharedSide).shared());

  int owner = bulkPtr->parallel_owner_rank(sharedSide);
  int otherProc = 1 - owner;
  stk::mesh::EntityProcVec entitiesToChange;
  if (owner == bulkPtr->parallel_rank()) {
    entitiesToChange = {stk::mesh::EntityProc(sharedSide, otherProc)};
  }
  bulkPtr->change_entity_owner(entitiesToChange);

  stk::mesh::create_interior_block_boundary_sides(*bulkPtr, meta.universal_part(), {&boundaryPart});

  sharedSide = bulkPtr->get_entity(meta.side_rank(), 100983);
  EXPECT_TRUE(bulkPtr->bucket(sharedSide).member(boundaryPart));

  EXPECT_TRUE(stk::mesh::check_interior_block_boundary_sides(*bulkPtr, meta.universal_part(), boundaryPart));
}

class TestSideSet : public stk::unit_test_util::MeshFixture
{
protected:
  void setup_2_block_mesh()
  {
    allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::MetaData& meta = get_meta();
    stk::mesh::Part& block2 = meta.declare_part_with_topology("block_2", stk::topology::HEX_8);
    stk::mesh::Part& surface1 = meta.declare_part_with_topology("surface_1", stk::topology::QUAD_4);
    stk::io::put_io_part_attribute(block2);
    stk::io::put_io_part_attribute(surface1);
    meta.set_part_id(block2, 2);
    meta.set_part_id(surface1, 1);

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

    stk::mesh::create_bulkdata_sidesets(get_bulk());

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

TEST_F(TestSideSet, createSingleSidedSideSetOnBlock1BetweenBlocks_WithSubsetSurfaceTouchingBlock1AndSurface1TouchingBothBlocks)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::MetaData& meta = get_meta();
    stk::mesh::Part& block1 = meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
    stk::mesh::Part& block2 = meta.declare_part_with_topology("block_2", stk::topology::HEX_8);
    stk::mesh::Part& surface1 = meta.declare_part_with_topology("surface_1", stk::topology::QUAD_4);
    stk::mesh::Part& surface1TouchingBlock1 = meta.declare_part_with_topology("surface_1_touching_block1", stk::topology::QUAD_4);
    meta.declare_part_subset(surface1, surface1TouchingBlock1);

    stk::io::put_io_part_attribute(block2);
    stk::io::put_io_part_attribute(surface1);
    stk::io::put_io_part_attribute(surface1TouchingBlock1);

    meta.set_surface_to_block_mapping(&surface1, stk::mesh::ConstPartVector{&block1, &block2});

    meta.set_part_id(block2, 2);
    meta.set_part_id(surface1, 1);
    meta.set_part_id(surface1TouchingBlock1, 1);

    stk::io::fill_mesh("generated:2x1x1", get_bulk());

    move_elems_from_block_to_block(get_bulk(), {2}, "block_1", "block_2");

    create_sideset(get_bulk(), surface1TouchingBlock1.name(), "block_1");
    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), stk::mesh::PartVector{}, false);
    move_sides_between_blocks_into_sideset_part(get_bulk(), "block_1", "block_2", surface1TouchingBlock1.name());

    check_sideset_elems(get_bulk(), "surface_1", {1});
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


namespace {
void create_sideset_observer(stk::mesh::BulkData& bulk, stk::mesh::Selector activeSelector = stk::mesh::Selector())
{
  if (!bulk.has_observer_type<stk::mesh::SidesetUpdater>()) {
    if (activeSelector == stk::mesh::Selector()) {
      activeSelector = bulk.mesh_meta_data().universal_part();
    }
    bulk.register_observer(std::make_shared<stk::mesh::IncrementalSidesetUpdater>(bulk, activeSelector));
  }
}
}

class SideSetModification : public stk::unit_test_util::MeshFixture
{
protected:
  stk::mesh::Entity create_hex_solo_side(const stk::mesh::PartVector& parts)
  {
    stk::mesh::PartVector quadParts = parts;
    stk::mesh::Part& quadTopoPart = get_bulk().mesh_meta_data().get_topology_root_part(stk::topology::QUAD_4);
    quadParts.push_back(&quadTopoPart);
    stk::mesh::Entity soloSide = get_bulk().declare_solo_side(quadParts);

    return soloSide;
  }

  stk::mesh::Entity create_hex_solo_side(unsigned nodeOffset, const stk::mesh::PartVector& parts)
  {
    stk::mesh::Part& nodePart = get_bulk().mesh_meta_data().get_topology_root_part(stk::topology::NODE);
    stk::mesh::PartVector nodeParts = {&nodePart};

    stk::mesh::PartVector quadParts = parts;
    stk::mesh::Part& quadTopoPart = get_bulk().mesh_meta_data().get_topology_root_part(stk::topology::QUAD_4);
    quadParts.push_back(&quadTopoPart);
    stk::mesh::Entity soloSide = get_bulk().declare_solo_side(quadParts);

    stk::mesh::Permutation perm = stk::mesh::Permutation::INVALID_PERMUTATION;
    stk::mesh::OrdinalVector scratch1, scratch2, scratch3;

    for(unsigned j = 0; j < 4; ++j)
    {
      unsigned id = nodeOffset + j;
      stk::mesh::Entity node = get_bulk().get_entity(stk::topology::NODE_RANK, id);
      if (!get_bulk().is_valid(node)) {
        node = get_bulk().declare_node(id, nodeParts);
      }
      get_bulk().declare_relation(soloSide, node, j, perm, scratch1, scratch2, scratch3);
    }

    return soloSide;
  }

  void snap_hex_solo_side_to_element(stk::mesh::Entity elem, stk::mesh::ConnectivityOrdinal ordinal, stk::mesh::Entity side)
  {
    stk::mesh::MetaData& meta = get_bulk().mesh_meta_data();
    ASSERT_EQ(stk::topology::ELEM_RANK, get_bulk().entity_rank(elem));
    ASSERT_EQ(meta.side_rank(), get_bulk().entity_rank(side));

    stk::mesh::EntityVector sideNodes;
    stk::mesh::impl::fill_element_side_nodes_from_topology(get_bulk().bucket(elem).topology(), get_bulk().begin_nodes(elem), ordinal, sideNodes);

    stk::mesh::Permutation perm = stk::mesh::Permutation::INVALID_PERMUTATION;
    stk::mesh::OrdinalVector scratch1, scratch2, scratch3;

    for(unsigned j = 0; j < sideNodes.size(); ++j) {
      get_bulk().declare_relation(side, sideNodes[j], j, perm, scratch1, scratch2, scratch3);
    }

    get_bulk().declare_relation(elem, side, ordinal);
  }

  void setup_test(const std::string& meshSpec)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    m_soloPart = &get_meta().declare_part("soloPart", get_meta().side_rank());

    stk::io::fill_mesh(meshSpec, get_bulk());
    create_sideset_observer(get_bulk());
    m_sidePart = get_meta().get_part("surface_1");

    if(m_sidePart == nullptr) {
      m_sidePart = &get_meta().declare_part_with_topology("sidesetPart", stk::topology::QUAD_4);
      get_meta().set_part_id(*m_sidePart, 100u);

      stk::mesh::Part* block_1 = get_meta().get_part("block_1");
      ASSERT_TRUE(block_1 != nullptr);
      get_meta().set_surface_to_block_mapping(m_sidePart, stk::mesh::ConstPartVector{block_1});
    }

    m_sideset = &get_bulk().create_sideset(*m_sidePart);
    m_elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
    ASSERT_TRUE(get_bulk().is_valid(m_elem1));

    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK), faces);

    for(stk::mesh::Entity face : faces) {
      ASSERT_TRUE(get_bulk().is_valid(face));
      ASSERT_TRUE(get_bulk().bucket(face).member(*m_sidePart));
    }
  }

  void create_declared_element_side(stk::mesh::EntityId elemId, stk::mesh::ConnectivityOrdinal ordinal, const stk::mesh::PartVector& addParts)
  {
    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    ASSERT_TRUE(get_bulk().is_valid(elem));

    const stk::mesh::ConnectivityOrdinal* ordinals = get_bulk().begin_ordinals(elem, stk::topology::FACE_RANK);
    for(unsigned i=0; i<get_bulk().num_connectivity(elem, stk::topology::FACE_RANK); i++) {
      ASSERT_TRUE(ordinals[i] != ordinal);
    }

    get_bulk().modification_begin();
    m_declaredSideNotInSideset = get_bulk().declare_element_side(elem, ordinal, addParts);
    get_bulk().modification_end();
  }

  void setup_test_with_declared_element_side(const std::string& meshSpec, stk::mesh::EntityId elemId, stk::mesh::ConnectivityOrdinal ordinal)
  {
    setup_test(meshSpec);

    stk::mesh::PartVector addParts;
    create_declared_element_side(elemId, ordinal, addParts);

    ASSERT_TRUE(get_bulk().is_valid(m_declaredSideNotInSideset));
    ASSERT_FALSE(get_bulk().bucket(m_declaredSideNotInSideset).member(*m_sidePart));
  }

  void setup_test_with_declared_element_side_in_solo_part(const std::string& meshSpec, stk::mesh::EntityId elemId, stk::mesh::ConnectivityOrdinal ordinal)
  {
    setup_test(meshSpec);

    stk::mesh::PartVector addParts;
    addParts.push_back(m_soloPart);
    create_declared_element_side(elemId, ordinal, addParts);

    ASSERT_TRUE(get_bulk().is_valid(m_declaredSideNotInSideset));
    ASSERT_FALSE(get_bulk().bucket(m_declaredSideNotInSideset).member(*m_sidePart));
    ASSERT_TRUE(get_bulk().bucket(m_declaredSideNotInSideset).member(*m_soloPart));
    ASSERT_FALSE(get_bulk().does_sideset_exist(*m_soloPart));
  }

  void create_solo_side_no_mesh_modification(bool insertIntoSoloPart)
  {
    unsigned numNodes = stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::NODE_RANK));

    unsigned nodeOffset = numNodes + 1;

    stk::mesh::PartVector addParts;
    if(insertIntoSoloPart) {
      addParts.push_back(m_soloPart);
    }

    m_soloSide = create_hex_solo_side(nodeOffset, addParts);
  }

  void create_solo_side(bool insertIntoSoloPart)
  {
    get_bulk().modification_begin();
    create_solo_side_no_mesh_modification(insertIntoSoloPart);
    get_bulk().modification_end();
  }

  void verify_test_with_solo_side(bool isSoloSideInSoloPart)
  {
    ASSERT_TRUE(get_bulk().is_valid(m_soloSide));
    ASSERT_FALSE(get_bulk().does_sideset_exist(*m_soloPart));
    ASSERT_FALSE(get_bulk().bucket(m_soloSide).member(*m_sidePart));
    ASSERT_EQ(isSoloSideInSoloPart, get_bulk().bucket(m_soloSide).member(*m_soloPart));
  }

  void setup_test_with_solo_side(const std::string& meshSpec)
  {
    bool isSoloSideInSoloPart = false;
    setup_test(meshSpec);
    create_solo_side(isSoloSideInSoloPart);
    verify_test_with_solo_side(isSoloSideInSoloPart);
  }

  void setup_test_with_solo_side_in_solo_part(const std::string& meshSpec)
  {
    bool isSoloSideInSoloPart = true;
    setup_test(meshSpec);
    create_solo_side(isSoloSideInSoloPart);
    verify_test_with_solo_side(isSoloSideInSoloPart);
  }

  stk::mesh::Part* m_sidePart = nullptr;
  stk::mesh::SideSet* m_sideset = nullptr;
  stk::mesh::Entity m_soloSide;
  stk::mesh::Part* m_soloPart = nullptr;
  stk::mesh::Entity m_declaredSideNotInSideset;
  stk::mesh::Entity m_elem1;
};

TEST_F(SideSetModification, removeSidesetEntry_AfterDestroyElement)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    stk::mesh::ConnectivityOrdinal inputSidesetOrdinal = 3;
    setup_test("generated:1x1x1|sideset:x");
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));

    get_bulk().modification_begin();
    get_bulk().destroy_entity(m_elem1);
    get_bulk().modification_end();

    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(0u, m_sideset->size());
  }
}

TEST_F(SideSetModification, keepSidesetEntry_AfterDestroyElementNotConnectedToSideset)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    stk::mesh::ConnectivityOrdinal inputSidesetOrdinal = 3;
    setup_test("generated:2x1x1|sideset:x");
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));

    stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
    ASSERT_TRUE(get_bulk().is_valid(elem2));
    get_bulk().modification_begin();
    get_bulk().destroy_entity(elem2);
    get_bulk().modification_end();

    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::ELEM_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
  }
}

TEST_F(SideSetModification, removeSidesetEntry_AfterDestroyFace)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    stk::mesh::ConnectivityOrdinal inputSidesetOrdinal = 3;
    setup_test("generated:1x1x1|sideset:x");
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));

    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK), faces);
    get_bulk().modification_begin();
    bool destroyFaceWithUpwardRelationToElement = get_bulk().destroy_entity(faces[0]);
    EXPECT_FALSE(destroyFaceWithUpwardRelationToElement);
    get_bulk().modification_end();

    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
  }
}

TEST_F(SideSetModification, removeSidesetEntry_AfterFaceChangePart)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    stk::mesh::ConnectivityOrdinal inputSidesetOrdinal = 3;
    setup_test("generated:1x1x1|sideset:x");
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));

    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK), faces);
    get_bulk().modification_begin();
    get_bulk().change_entity_parts(faces[0], stk::mesh::PartVector{}, stk::mesh::PartVector{m_sidePart});
    get_bulk().modification_end();

    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(0u, m_sideset->size());
    ASSERT_FALSE(get_bulk().bucket(faces[0]).member(*m_sidePart));
  }
}

TEST_F(SideSetModification, removeSidesetEntry_AfterDestroyRelationFromElementToFace)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    stk::mesh::ConnectivityOrdinal inputSidesetOrdinal = 3;
    setup_test("generated:1x1x1|sideset:x");
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));

    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK), faces);
    get_bulk().modification_begin();
    stk::mesh::ConnectivityOrdinal ordinal = 3;
    EXPECT_TRUE(get_bulk().destroy_relation(m_elem1, faces[0], ordinal));
    EXPECT_TRUE(get_bulk().destroy_entity(faces[0]));
    get_bulk().modification_end();

    EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(0u, m_sideset->size());
  }
}

TEST_F(SideSetModification, removeSidesetEntry_AfterFaceChangePartAndThenDestroyRelationFromElementToFace)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    stk::mesh::ConnectivityOrdinal inputSidesetOrdinal = 3;
    setup_test("generated:1x1x1|sideset:x");
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));

    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK), faces);
    get_bulk().modification_begin();
    get_bulk().change_entity_parts(faces[0], stk::mesh::PartVector{}, stk::mesh::PartVector{m_sidePart});
    stk::mesh::ConnectivityOrdinal ordinal = 3;
    EXPECT_TRUE(get_bulk().destroy_relation(m_elem1, faces[0], ordinal));
    get_bulk().modification_end();

    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(0u, m_sideset->size());
  }
}

TEST_F(SideSetModification, removeSidesetEntry_AfterDestroyRelationFromElementToFaceAndThenFaceChangePart)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    stk::mesh::ConnectivityOrdinal inputSidesetOrdinal = 3;
    setup_test("generated:1x1x1|sideset:x");
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));

    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK), faces);
    get_bulk().modification_begin();
    stk::mesh::ConnectivityOrdinal ordinal = 3;
    EXPECT_TRUE(get_bulk().destroy_relation(m_elem1, faces[0], ordinal));
    get_bulk().change_entity_parts(faces[0], stk::mesh::PartVector{}, stk::mesh::PartVector{m_sidePart});
    get_bulk().modification_end();

    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(0u, m_sideset->size());
  }
}

TEST_F(SideSetModification, keepSidesetEntry_AfterDestroyFaceNotInSideset)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    stk::mesh::ConnectivityOrdinal inputSidesetOrdinal = 3;
    setup_test_with_solo_side("generated:1x1x1|sideset:x");
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));

    get_bulk().modification_begin();
    get_bulk().destroy_entity(m_soloSide);
    get_bulk().modification_end();

    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));
  }
}

TEST_F(SideSetModification, keepSidesetEntry_AfterFaceChangePartNotInSideset)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    stk::mesh::ConnectivityOrdinal inputSidesetOrdinal = 3;
    setup_test_with_solo_side_in_solo_part("generated:1x1x1|sideset:x");
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));

    get_bulk().modification_begin();
    get_bulk().change_entity_parts(m_soloSide, stk::mesh::PartVector{}, stk::mesh::PartVector{m_soloPart});
    get_bulk().modification_end();

    EXPECT_EQ(2u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));
  }
}

TEST_F(SideSetModification, keepSidesetEntry_AfterDestroyRelationFromElementToFaceNotInSideset)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    stk::mesh::ConnectivityOrdinal inputSidesetOrdinal = 3;
    stk::mesh::ConnectivityOrdinal ordinal = 5;
    stk::mesh::EntityId elemId = 1u;
    setup_test_with_declared_element_side("generated:1x1x1|sideset:x", elemId, ordinal);
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));

    get_bulk().modification_begin();
    EXPECT_TRUE(get_bulk().destroy_relation(m_elem1, m_declaredSideNotInSideset, ordinal));
    get_bulk().modification_end();

    EXPECT_EQ(2u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));
  }
}

TEST_F(SideSetModification, keepSidesetEntry_AfterFaceChangePartAndThenDestroyRelationFromElementToFaceNotInSideset)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    stk::mesh::ConnectivityOrdinal inputSidesetOrdinal = 3;
    stk::mesh::ConnectivityOrdinal ordinal = 5;
    stk::mesh::EntityId elemId = 1u;
    setup_test_with_declared_element_side_in_solo_part("generated:1x1x1|sideset:x", elemId, ordinal);
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));

    get_bulk().modification_begin();
    get_bulk().change_entity_parts(m_declaredSideNotInSideset, stk::mesh::PartVector{}, stk::mesh::PartVector{m_soloPart});
    EXPECT_TRUE(get_bulk().destroy_relation(m_elem1, m_declaredSideNotInSideset, ordinal));
    get_bulk().modification_end();

    EXPECT_EQ(2u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));
  }
}

TEST_F(SideSetModification, keepSidesetEntry_AfterDestroyRelationFromElementToFaceAndThenFaceChangePartNotInSideset)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    stk::mesh::ConnectivityOrdinal inputSidesetOrdinal = 3;
    stk::mesh::ConnectivityOrdinal ordinal = 5;
    stk::mesh::EntityId elemId = 1u;
    setup_test_with_declared_element_side_in_solo_part("generated:1x1x1|sideset:x", elemId, ordinal);
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));

    get_bulk().modification_begin();
    EXPECT_TRUE(get_bulk().destroy_relation(m_elem1, m_declaredSideNotInSideset, ordinal));
    get_bulk().change_entity_parts(m_declaredSideNotInSideset, stk::mesh::PartVector{}, stk::mesh::PartVector{m_soloPart});
    get_bulk().modification_end();

    EXPECT_EQ(2u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, inputSidesetOrdinal));
  }
}

TEST_F(SideSetModification, addSidesetEntry_AfterCreateFaceInSideset)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    setup_test("generated:1x1x1");
    EXPECT_EQ(0u, m_sideset->size());

    get_bulk().modification_begin();
    stk::mesh::ConnectivityOrdinal ordinal = 3;
    stk::mesh::Entity sidesetSide = get_bulk().declare_element_side(m_elem1, ordinal, stk::mesh::PartVector{m_sidePart});
    get_bulk().modification_end();

    EXPECT_TRUE(get_bulk().is_valid(sidesetSide));
    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, ordinal));
  }
}

TEST_F(SideSetModification, noAddedSidesetEntry_AfterFaceChangePartWithoutDeclareElementToFaceRelation)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    setup_test_with_solo_side("generated:1x1x1");
    EXPECT_EQ(0u, m_sideset->size());
    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));

    get_bulk().modification_begin();
    get_bulk().change_entity_parts(m_soloSide, stk::mesh::PartVector{m_sidePart}, stk::mesh::PartVector{});
    get_bulk().modification_end();

    EXPECT_EQ(0u, m_sideset->size());
    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
  }
}


TEST_F(SideSetModification, noAddedSidesetEntry_AfterDeclareElementToFaceRelationWithoutFaceChangePart)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    setup_test("generated:1x1x1");
    EXPECT_EQ(0u, m_sideset->size());
    EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));

    get_bulk().modification_begin();
    stk::mesh::ConnectivityOrdinal ordinal = 3;
    stk::mesh::Entity side = get_bulk().declare_element_side(m_elem1, ordinal, stk::mesh::PartVector{});
    get_bulk().modification_end();

    EXPECT_TRUE(get_bulk().is_valid(side));
    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(0u, m_sideset->size());
  }
}

TEST_F(SideSetModification, addSidesetEntry_AfterDeclareElementToFaceRelationAndThenFaceChangePart)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    setup_test("generated:1x1x1");
    EXPECT_EQ(0u, m_sideset->size());
    EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));

    get_bulk().modification_begin();
    stk::mesh::ConnectivityOrdinal ordinal = 3;
    stk::mesh::Entity side = get_bulk().declare_element_side(m_elem1, ordinal, stk::mesh::PartVector{});
    get_bulk().change_entity_parts(side, stk::mesh::PartVector{m_sidePart}, stk::mesh::PartVector{});
    get_bulk().modification_end();

    EXPECT_TRUE(get_bulk().is_valid(side));
    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, ordinal));
  }
}

TEST_F(SideSetModification, addSidesetEntry_AfterFaceChangePartAndThenDeclareElementToFaceRelation)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    setup_test_with_solo_side("generated:1x1x1");
    EXPECT_EQ(0u, m_sideset->size());
    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));

    get_bulk().modification_begin();
    get_bulk().change_entity_parts(m_soloSide, stk::mesh::PartVector{m_sidePart}, stk::mesh::PartVector{});

    stk::mesh::ConnectivityOrdinal ordinal = 3;
    get_bulk().declare_relation(m_elem1, m_soloSide, ordinal);
    get_bulk().modification_end();

    EXPECT_TRUE(get_bulk().is_valid(m_soloSide));
    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, ordinal));
  }
}

TEST_F(SideSetModification, noAddedSidesetEntry_AfterCreateFace_ForFaceNotInSideset)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    setup_test("generated:1x1x1");
    EXPECT_EQ(0u, m_sideset->size());
    EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));

    get_bulk().modification_begin();
    stk::mesh::ConnectivityOrdinal sidesetOrdinal = 3;
    get_bulk().declare_element_side(m_elem1, sidesetOrdinal, stk::mesh::PartVector{m_sidePart});
    bool isInSoloPartWhichTriggersPartChange = false;
    create_solo_side_no_mesh_modification(isInSoloPartWhichTriggersPartChange);
    get_bulk().modification_end();

    EXPECT_EQ(2u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, sidesetOrdinal));
  }
}

TEST_F(SideSetModification, noAddedSidesetEntry_AfterCreateFaceWithFaceChangePart_ForFaceNotInSideset)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    setup_test("generated:1x1x1");
    EXPECT_EQ(0u, m_sideset->size());
    EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));

    get_bulk().modification_begin();
    stk::mesh::ConnectivityOrdinal sidesetOrdinal = 3;
    get_bulk().declare_element_side(m_elem1, sidesetOrdinal, stk::mesh::PartVector{m_sidePart});
    bool isInSoloPartWhichTriggersPartChange = true;
    create_solo_side_no_mesh_modification(isInSoloPartWhichTriggersPartChange);
    get_bulk().modification_end();

    EXPECT_EQ(2u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, sidesetOrdinal));
  }
}

TEST_F(SideSetModification, noAddedSidesetEntry_AfterCreateFaceWithDeclareElementToFaceRelation_ForFaceNotInSideset)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    setup_test("generated:1x1x1");
    EXPECT_EQ(0u, m_sideset->size());
    EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));

    get_bulk().modification_begin();
    stk::mesh::ConnectivityOrdinal sidesetOrdinal = 3;
    get_bulk().declare_element_side(m_elem1, sidesetOrdinal, stk::mesh::PartVector{m_sidePart});
    stk::mesh::Entity soloSide = create_hex_solo_side(stk::mesh::PartVector{});
    stk::mesh::ConnectivityOrdinal soloSideOrdinal = 5;
    snap_hex_solo_side_to_element(m_elem1, soloSideOrdinal, soloSide);
    get_bulk().modification_end();

    EXPECT_EQ(2u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, sidesetOrdinal));
  }
}

TEST_F(SideSetModification, noAddedSidesetEntry_AfterCreateFaceWithFaceChangePartAndThenDeclareElementToFaceRelation_ForFaceNotInSideset)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    setup_test("generated:1x1x1");
    EXPECT_EQ(0u, m_sideset->size());
    EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));

    get_bulk().modification_begin();
    stk::mesh::ConnectivityOrdinal sidesetOrdinal = 3;
    get_bulk().declare_element_side(m_elem1, sidesetOrdinal, stk::mesh::PartVector{m_sidePart});
    stk::mesh::Entity soloSide = create_hex_solo_side(stk::mesh::PartVector{m_soloPart});
    stk::mesh::ConnectivityOrdinal soloSideOrdinal = 5;
    snap_hex_solo_side_to_element(m_elem1, soloSideOrdinal, soloSide);
    get_bulk().modification_end();

    EXPECT_EQ(2u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, sidesetOrdinal));
  }
}

TEST_F(SideSetModification, noAddedSidesetEntry_AfterCreateFaceWithDeclareElementToFaceRelationAndThenFaceChangePart_ForFaceNotInSideset)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    setup_test("generated:1x1x1");
    EXPECT_EQ(0u, m_sideset->size());
    EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));

    get_bulk().modification_begin();
    stk::mesh::ConnectivityOrdinal sidesetOrdinal = 3;
    get_bulk().declare_element_side(m_elem1, sidesetOrdinal, stk::mesh::PartVector{m_sidePart});
    stk::mesh::Entity soloSide = create_hex_solo_side(stk::mesh::PartVector{});
    stk::mesh::ConnectivityOrdinal soloSideOrdinal = 5;
    snap_hex_solo_side_to_element(m_elem1, soloSideOrdinal, soloSide);
    get_bulk().change_entity_parts(soloSide, stk::mesh::PartVector{m_soloPart}, stk::mesh::PartVector{});
    get_bulk().modification_end();

    EXPECT_EQ(2u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(1u, m_sideset->size());
    EXPECT_TRUE(m_sideset->contains(m_elem1, sidesetOrdinal));
  }
}

TEST_F(SideSetModification, emptySideset_AfterRemoveOneConnectedElementFromBlock)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_test("generated:1x1x1");
  EXPECT_EQ(0u, m_sideset->size());
  EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));

  get_bulk().modification_begin();
  stk::mesh::ConnectivityOrdinal sidesetOrdinal = 5;
  stk::mesh::Entity side = get_bulk().declare_element_side(m_elem1, sidesetOrdinal, stk::mesh::PartVector{m_sidePart});
  get_bulk().modification_end();

  EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
  EXPECT_EQ(1u, m_sideset->size());
  EXPECT_TRUE(m_sideset->contains(m_elem1, sidesetOrdinal));
  EXPECT_EQ(1u, get_bulk().num_elements(side));

  stk::mesh::Part* block_1 = get_meta().get_part("block_1");
  ASSERT_TRUE(block_1 != nullptr);
  ASSERT_TRUE(get_bulk().bucket(side).member(*block_1));

  get_bulk().modification_begin();
  get_bulk().change_entity_parts(m_elem1, stk::mesh::PartVector{}, stk::mesh::PartVector{block_1});
  get_bulk().modification_end();

  EXPECT_FALSE(get_bulk().bucket(side).member(*block_1));
  EXPECT_EQ(0u, m_sideset->size());
}

TEST_F(SideSetModification, noChange_AfterRemoveThenAddInSameModCycleOneConnectedElementFromBlock)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_test("generated:1x1x1");
  EXPECT_EQ(0u, m_sideset->size());
  EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));

  get_bulk().modification_begin();
  stk::mesh::ConnectivityOrdinal sidesetOrdinal = 5;
  stk::mesh::Entity side = get_bulk().declare_element_side(m_elem1, sidesetOrdinal, stk::mesh::PartVector{m_sidePart});
  get_bulk().modification_end();

  EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
  EXPECT_EQ(1u, m_sideset->size());
  EXPECT_TRUE(m_sideset->contains(m_elem1, sidesetOrdinal));
  EXPECT_EQ(1u, get_bulk().num_elements(side));

  stk::mesh::Part* block_1 = get_meta().get_part("block_1");
  ASSERT_TRUE(block_1 != nullptr);
  ASSERT_TRUE(get_bulk().bucket(side).member(*block_1));

  get_bulk().modification_begin();
  get_bulk().change_entity_parts(m_elem1, stk::mesh::PartVector{}, stk::mesh::PartVector{block_1});
  get_bulk().change_entity_parts(m_elem1, stk::mesh::PartVector{block_1}, stk::mesh::PartVector{});
  get_bulk().modification_end();

  EXPECT_TRUE(get_bulk().bucket(side).member(*block_1));
  EXPECT_TRUE(get_bulk().bucket(m_elem1).member(*block_1));
  EXPECT_EQ(1u, m_sideset->size());
  EXPECT_TRUE(m_sideset->contains(m_elem1, sidesetOrdinal));
  EXPECT_EQ(1u, get_bulk().num_elements(side));
}

TEST_F(SideSetModification, nonEmptyInternalSideset_AfterRemoveOneConnectedElementFromBlock)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_test("generated:1x1x2");
  EXPECT_EQ(0u, m_sideset->size());
  EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));

  get_bulk().modification_begin();
  stk::mesh::ConnectivityOrdinal sidesetOrdinal = 5;
  stk::mesh::Entity side = get_bulk().declare_element_side(m_elem1, sidesetOrdinal, stk::mesh::PartVector{m_sidePart});
  get_bulk().modification_end();

  EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK)));
  EXPECT_EQ(2u, m_sideset->size());
  EXPECT_TRUE(m_sideset->contains(m_elem1, sidesetOrdinal));
  EXPECT_EQ(2u, get_bulk().num_elements(side));

  stk::mesh::Part* block_1 = get_meta().get_part("block_1");
  ASSERT_TRUE(block_1 != nullptr);
  ASSERT_TRUE(get_bulk().bucket(side).member(*block_1));

  get_bulk().modification_begin();
  get_bulk().change_entity_parts(m_elem1, stk::mesh::PartVector{}, stk::mesh::PartVector{block_1});
  get_bulk().modification_end();

  EXPECT_TRUE(get_bulk().bucket(side).member(*block_1));
  EXPECT_EQ(1u, m_sideset->size());
  EXPECT_EQ(2u, get_bulk().num_elements(side));
}

TEST(CreateAndWrite, DISABLED_textmesh_shell_quad_4_EdgeSides)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
//shell-quad-4 mesh:
//       6
// 3*----*----*9
//  | E2 | E4 |
//  |    |    |
// 2*---5*----*8
//  | E1 | E3 |
//  |    |    |
// 1*----*----*7
//       4
//
  const std::string meshDesc =
       "0,1,SHELL_QUAD_4, 1,4,5,2, block_1\n\
        0,2,SHELL_QUAD_4, 2,5,6,3, block_1\n\
        0,3,SHELL_QUAD_4, 4,7,8,5, block_1\n\
        0,4,SHELL_QUAD_4, 5,8,9,6, block_1|sideset:name=surface_1; data=1,3,3,3, 3,4,4,4, 4,5,2,5, 2,6,1,6";

  std::vector<double> coords = {0,0,0,  0,1,0,  0,2,0,
                                1,0,0,  1,1,0,  1,2,0,
                                2,0,0,  2,1,0,  2,2,0};

  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

  stk::io::write_mesh("shellq4_edge_sides.g", *bulk);
}

TEST(CreateAndWrite, DISABLED_textmesh_shell_quad_4_FullExteriorSkin)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
  const std::string meshDesc =
       "0,1,SHELL_QUAD_4, 1,4,5,2, block_1\n\
        0,2,SHELL_QUAD_4, 2,5,6,3, block_1\n\
        0,3,SHELL_QUAD_4, 4,7,8,5, block_1\n\
        0,4,SHELL_QUAD_4, 5,8,9,6, block_1|sideset:name=surface_1; data=1,1,2,1,3,1,4,1, 1,2,2,2,3,2,4,2, 1,3,3,3, 3,4,4,4, 4,5,2,5, 2,6,1,6";

  std::vector<double> coords = {0,0,0,  0,1,0,  0,2,0,
                                1,0,0,  1,1,0,  1,2,0,
                                2,0,0,  2,1,0,  2,2,0};

  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

  stk::io::write_mesh("shellq4_full_exterior_skin.g", *bulk);
}

TEST(CreateAndWrite, DISABLED_textmesh_shell_quad_8_EdgeSides)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
//shell-quad-8 mesh:
//       6         12
// 3*----*----*9----*----*15
//  |         |          |
//  |         |          |
// 2*         *8         *14
//  |         |          |
//  |         |          |
// 1*----*----*7----*----*13
//       4         10
//
  const std::string meshDesc =
       "0,1,SHELL_QUAD_8, 1,7,9,3,4,8,6,2, block_1\n\
        0,2,SHELL_QUAD_8, 7,13,15,9,10,14,12,8, block_2|sideset:name=surface_1; data=1,1,1,3,2,3,2,6";
//Note: this sideset spec creates a face-side (elem 1, side 1), and 3 edge-sides.
//stk-io can't currently read the resulting exodus mesh file if the face-side is present. (Because
//the surface-df field gets conflicting sides (restrictions)).

  std::vector<double> coords = {
     /*1*/0,0,0,  /*2*/0,1,0,  /*3*/0,2,0,
     /*4*/1,0,0,  /*6*/1,2,0,
     /*7*/2,0,0,  /*8*/2,1,0,  /*9*/2,2,0,
     /*10*/3,0,0,  /*12*/3,2,0,
     /*13*/4,0,0,  /*14*/4,1,0,  /*15*/4,2,0};

  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

  stk::io::write_mesh("shellq8_2blk_face_edge_sides.g", *bulk);
}

TEST(CreateAndWrite, DISABLED_textmesh_shell_tri_3_EdgeSides)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
//shell-tri-3 mesh:
//       4
//  2*---*---*6
//   |E1/|E3/|
//   | / | / |
//   |/E2|/E4|
//  1*---*---*5
//       3
//
  const std::string meshDesc =
       "0,1,SHELL_TRI_3, 1,4,2, block_1\n\
        0,2,SHELL_TRI_3, 1,3,4, block_1\n\
        0,3,SHELL_TRI_3, 3,6,4, block_1\n\
        0,4,SHELL_TRI_3, 3,5,6, block_1|sideset:name=surface_1; data=2,3,4,3, 4,4, 3,4,1,4, 1,5";

  std::vector<double> coords = {0,0,0,  0,1,0,
                                1,0,0,  1,1,0,
                                2,0,0,  2,1,0};

  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

  stk::io::write_mesh("shellt3_edge_sides.g", *bulk);
}

TEST(CreateAndWrite, DISABLED_textmesh_shell_tri_3_FullExteriorSkin)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
  const std::string meshDesc =
       "0,1,SHELL_TRI_3, 1,4,2, block_1\n\
        0,2,SHELL_TRI_3, 1,3,4, block_1\n\
        0,3,SHELL_TRI_3, 3,6,4, block_1\n\
        0,4,SHELL_TRI_3, 3,5,6, block_1|sideset:name=surface_1; data=1,1,2,1,3,1,4,1, 1,2,2,2,3,2,4,2, 2,3,4,3, 4,4, 3,4,1,4, 1,5";

  std::vector<double> coords = {0,0,0,  0,1,0,
                                1,0,0,  1,1,0,
                                2,0,0,  2,1,0};

  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

  stk::io::write_mesh("shellt3_full_exterior_skin.g", *bulk);
}

class InternalSideSet : public stk::unit_test_util::MeshFixture
{
protected:
  void setup_internal_sideset_test()
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    get_bulk().initialize_face_adjacent_element_graph();
    create_sideset_observer(get_bulk());

    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n\
        1,2,HEX_8,5,6,7,8,9,10,11,12,block_1";

        std::vector<double> coordinates = {0,0,0,  1,0,0,   1,1,0,   0,1,0,
        0,0,1,  1,0,1,   1,1,1,   0,1,1,
        0,0,2,  1,0,2,   1,1,2,   0,1,2};

    stk::mesh::Part& sidesetPart = get_meta().declare_part(sidesetName, stk::topology::FACE_RANK);
    get_meta().set_part_id(sidesetPart, 100u);

    stk::unit_test_util::setup_text_mesh(
          get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
    create_sideset(get_bulk(), sidesetName, "block_1");

    std::vector<stk::mesh::SideSet*> sidesets = get_bulk().get_sidesets();
    ASSERT_EQ(1u, sidesets.size());

    sidesets[0]->set_accept_all_internal_non_coincident_entries(false);
  }

  void execute_internal_sideset_test(stk::mesh::EntityId elemId, stk::mesh::ConnectivityOrdinal ordinal)
  {
    std::vector<stk::mesh::SideSet*> sidesets = get_bulk().get_sidesets();
    ASSERT_EQ(1u, sidesets.size());

    stk::mesh::Part* sidesetPart = get_meta().get_part(sidesetName);
    ASSERT_TRUE(sidesetPart != nullptr);

    get_bulk().modification_begin();
    stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK,elemId);

    if(get_bulk().is_valid(element) && get_bulk().bucket(element).owned()) {
      sidesets[0]->add(element, ordinal);
      get_bulk().declare_element_side(element, ordinal, stk::mesh::ConstPartVector{sidesetPart});
    }
    get_bulk().modification_end();
  }

  void test_internal_sideset_test(stk::mesh::EntityId faceId, const std::vector<unsigned>& expectedNumberOfSidesetEntriesPerProc)
  {
    std::vector<stk::mesh::SideSet*> sidesets = get_bulk().get_sidesets();
    ASSERT_EQ(1u, sidesets.size());

    stk::mesh::Entity face = get_bulk().get_entity(stk::topology::FACE_RANK,faceId);
    ASSERT_TRUE(get_bulk().is_valid(face));

    ASSERT_TRUE(expectedNumberOfSidesetEntriesPerProc.size() >= (unsigned)get_bulk().parallel_size());
    EXPECT_EQ(expectedNumberOfSidesetEntriesPerProc[get_bulk().parallel_rank()], sidesets[0]->size());
  }

  void run_internal_sideset_test(stk::mesh::EntityId elemId, stk::mesh::ConnectivityOrdinal ordinal, stk::mesh::EntityId faceId,
                                 const std::vector<unsigned>& expectedNumberOfSidesetEntriesPerProc)
  {
    if(stk::parallel_machine_size(MPI_COMM_WORLD) != 2) {return;}

    setup_internal_sideset_test();
    execute_internal_sideset_test(elemId, ordinal);
    test_internal_sideset_test(faceId, expectedNumberOfSidesetEntriesPerProc);
  }

  std::string sidesetName = "sidesetPart";
};

TEST_F(InternalSideSet, maintainSingleSidedAfterParallelCreationOnProcessorBoundary_onP0)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    stk::mesh::EntityId elemId = 1;
    stk::mesh::ConnectivityOrdinal ordinal = 5;
    stk::mesh::EntityId faceId = 16;
    std::vector<unsigned> expectedNumberOfSidesetEntriesPerProc = {1, 0};
    run_internal_sideset_test(elemId, ordinal, faceId, expectedNumberOfSidesetEntriesPerProc);
  }
}

TEST_F(InternalSideSet, maintainSingleSidedAfterParallelCreationOnProcessorBoundary_onP1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    stk::mesh::EntityId elemId = 2;
    stk::mesh::ConnectivityOrdinal ordinal = 4;
    stk::mesh::EntityId faceId = 16;
    std::vector<unsigned> expectedNumberOfSidesetEntriesPerProc = {0, 1};
    run_internal_sideset_test(elemId, ordinal, faceId, expectedNumberOfSidesetEntriesPerProc);
  }
}

class ParallelCoincidence : public stk::unit_test_util::MeshFixture
{
public:
  struct ParallelCoincidenceEntry
  {
    ParallelCoincidenceEntry(const stk::mesh::EntityId elemId_, const stk::mesh::ConnectivityOrdinal ordinal_, const bool expectedCoincidence_)
      : elemId(elemId_), ordinal(ordinal_), expectedCoincidence(expectedCoincidence_) { }

    stk::mesh::EntityId elemId;
    stk::mesh::ConnectivityOrdinal ordinal;
    bool expectedCoincidence;
  };

protected:
  void initialize_mesh()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    get_bulk().initialize_face_adjacent_element_graph();
    create_sideset_observer(get_bulk());
  }

  void setup_coincident_mesh()
  {
    initialize_mesh();

    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n\
        1,2,HEX_8,1,2,3,4,9,10,11,12,block_1";

        std::vector<double> coordinates = {0,0,0,  1,0,0,   1,1,0,   0,1,0,
        0,0,1,  1,0,1,   1,1,1,   0,1,1,
        0,0,2,  1,0,2,   1,1,2,   0,1,2};

    stk::unit_test_util::setup_text_mesh(
          get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
  }

  void setup_non_coincident_mesh()
  {
    initialize_mesh();

    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n\
        1,2,HEX_8,5,6,7,8,9,10,11,12,block_1";

        std::vector<double> coordinates = {0,0,0,  1,0,0,   1,1,0,   0,1,0,
        0,0,1,  1,0,1,   1,1,1,   0,1,1,
        0,0,2,  1,0,2,   1,1,2,   0,1,2};

    stk::unit_test_util::setup_text_mesh(
          get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
  }

  void test_parallel_coincidence(const std::vector<ParallelCoincidenceEntry>& expectedValues)
  {
    std::ostream* outputStream = nullptr;

#ifndef NDEBUG
    std::ostringstream oss;
    oss << "output." << get_bulk().parallel_rank();
    std::string fileName = oss.str();
    std::ofstream ofs(fileName, std::ofstream::out);
    STK_ThrowRequireMsg(ofs.fail() == false, "Failed to open debug file: " << oss.str());
    outputStream = &ofs;
#endif

    stk::mesh::SideSetHelper helper(get_bulk(), get_bulk().mesh_meta_data().universal_part(), outputStream);

    for(const ParallelCoincidenceEntry& entry : expectedValues) {
      stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, entry.elemId);

      if(get_bulk().is_valid(element) && get_bulk().bucket(element).owned()) {
        bool isParallelCoincident = helper.element_side_has_coincidence_using_elem_elem_graph(stk::mesh::Entity(), element, entry.ordinal);
        EXPECT_EQ(entry.expectedCoincidence, isParallelCoincident);
      }
    }
#ifndef NDEBUG
    unlink(fileName.c_str());
#endif
  }
};

TEST_F(ParallelCoincidence, checkParallelCoincidenceWithElemElemGraph)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    setup_coincident_mesh();

    std::vector<ParallelCoincidenceEntry> expectedValues;
    expectedValues.emplace_back(1u, 4u, true);
    expectedValues.emplace_back(2u, 4u, true);

    test_parallel_coincidence(expectedValues);
  }
}

TEST_F(ParallelCoincidence, checkParallelNonCoincidenceWithElemElemGraph)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    setup_non_coincident_mesh();

    std::vector<ParallelCoincidenceEntry> expectedValues;
    expectedValues.emplace_back(1u, 5u, false);
    expectedValues.emplace_back(2u, 4u, false);

    test_parallel_coincidence(expectedValues);
  }
}

TEST(DeclareElementSide, baseScenario_sidesBetweenTwoTriangles)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                    .set_spatial_dimension(2).create();
  stk::mesh::Part& mySidePart = bulk->mesh_meta_data().declare_part("mySidePart");
  std::string meshDesc =
      "0,1,TRI_3_2D,1,2,3,block_1\n"
      "0,2,TRI_3_2D,2,4,3,block_1\n"
      "|dimension:2";

  std::vector<double> coords = {0,0,  0,1,  1,0,  1,1};

  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  EXPECT_EQ(0u, stk::mesh::count_entities(*bulk, meta.side_rank(), meta.universal_part()));

  stk::mesh::Entity elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity elem2 = bulk->get_entity(stk::topology::ELEM_RANK, 2);

  stk::mesh::PartVector parts = {&mySidePart};
  stk::mesh::Entity side11, side22;
  {
    bulk->modification_begin();
    side11 = bulk->declare_element_side(elem1, 1, parts);
    side22 = bulk->declare_element_side(elem2, 2, parts);
    bulk->modification_end();
  }

  EXPECT_EQ(1u, stk::mesh::count_entities(*bulk, meta.side_rank(), meta.universal_part()));
}

TEST(DeclareElementSide, destroyElemAndReconnectElem_sidesBetweenTwoTriangles)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                    .set_spatial_dimension(2).create();
  stk::mesh::Part& mySidePart = bulk->mesh_meta_data().declare_part("mySidePart");
  std::string meshDesc =
      "0,1,TRI_3_2D,1,2,3,block_1\n"
      "0,2,TRI_3_2D,2,4,3,block_1\n"
      "0,3,TRI_3_2D,1,3,5,block_1\n"
      "0,4,TRI_3_2D,2,6,4,block_1\n"
      "|dimension:2";

  std::vector<double> coords = {0,0,  0,1,  1,0,  1,1,  -0.5,0.5, 1.5,0.5};

  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  EXPECT_EQ(0u, stk::mesh::count_entities(*bulk, meta.side_rank(), meta.universal_part()));

  stk::mesh::Entity elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity elem2 = bulk->get_entity(stk::topology::ELEM_RANK, 2);

  stk::mesh::PartVector parts = {&mySidePart};
  stk::mesh::Entity side11, side22, side42;
  {
    bulk->modification_begin();
    side11 = bulk->declare_element_side(elem1, 1, parts);
    side22 = bulk->declare_element_side(elem2, 2, parts);
    bulk->modification_end();
  }

  EXPECT_EQ(1u, stk::mesh::count_entities(*bulk, meta.side_rank(), meta.universal_part()));

  //now we simulate a "collapse-edge" scenario from the NGS team:
  //element 2 is between elements 1 and 4. We delete element 2 and
  //reconnect element 4 (disconnect node 4 and connect to node 3) so
  //that elements 1 and 4 should now share a "graph edge" in stk-mesh's
  //elem-elem-graph.
  stk::mesh::Entity elem4 = bulk->get_entity(stk::topology::ELEM_RANK, 4);
  {
    bulk->modification_begin();
    EXPECT_TRUE(bulk->destroy_entity(elem2));
    stk::mesh::Entity node4 = bulk->get_entity(stk::topology::NODE_RANK, 4);
    EXPECT_TRUE(bulk->destroy_relation(elem4, node4, 2));
    stk::mesh::Entity node3 = bulk->get_entity(stk::topology::NODE_RANK, 3);
    bulk->declare_relation(elem4, node3, 2);
    bulk->modification_end();
  }

  {
    bulk->modification_begin();
    side11 = bulk->declare_element_side(elem1, 1, parts);
    side42 = bulk->declare_element_side(elem4, 2, parts);
    bulk->modification_end();
  }

  //if the elem-elem-graph correctly knows that elements 1 and 4 share a
  //graph edge, then there should still be just 1 side between them.
  EXPECT_EQ(1u, stk::mesh::count_entities(*bulk, meta.side_rank(), meta.universal_part()));
}

TEST(Skinning, createSidesForBlock)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                            .set_spatial_dimension(3)
                                                            .create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::mesh::Part& block2 = meta.declare_part_with_topology("block_2", stk::topology::HEX_8);
  stk::mesh::Part& surface1 = meta.declare_part_with_topology("surface_1", stk::topology::QUAD_4);
  stk::io::put_io_part_attribute(block2);
  stk::io::put_io_part_attribute(surface1);
  meta.set_part_id(block2, 2);
  meta.set_part_id(surface1, 1);

  stk::io::fill_mesh("generated:2x2x2", bulk);

  copy_elems_from_block_to_block(bulk, {1, 2}, "block_1", "block_2");

  stk::mesh::create_exposed_block_boundary_sides(bulk, block2, stk::mesh::PartVector{&surface1}, (!block2));
  EXPECT_EQ(10u, stk::mesh::count_entities(bulk, stk::topology::FACE_RANK, surface1));
}

TEST(Skinning, createSidesForShellQuad4Block)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  std::unique_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
//shell-quad-4 mesh:
//       6
// 3*----*----*9
//  | E2 | E4 |
//  |    |    |
// 2*---5*----*8
//  | E1 | E3 |
//  |    |    |
// 1*----*----*7
//       4
//
  stk::mesh::Part& skinPart = bulk->mesh_meta_data().declare_part("mySkin");
  const std::string meshDesc =
       "0,1,SHELL_QUAD_4, 1,4,5,2, block_1\n\
        0,2,SHELL_QUAD_4, 2,5,6,3, block_1\n\
        0,3,SHELL_QUAD_4, 4,7,8,5, block_1\n\
        0,4,SHELL_QUAD_4, 5,8,9,6, block_1|sideset:name=surface_1; data=1,3,3,3, 3,4,4,4, 4,5,2,5, 2,6,1,6";

  std::vector<double> coords = {0,0,0,  0,1,0,  0,2,0,
                                1,0,0,  1,1,0,  1,2,0,
                                2,0,0,  2,1,0,  2,2,0};

  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

  stk::mesh::create_exposed_block_boundary_sides(*bulk, bulk->mesh_meta_data().universal_part(), stk::mesh::PartVector{&skinPart});
  EXPECT_EQ(8u, stk::mesh::count_entities(*bulk, stk::topology::FACE_RANK, skinPart));
}

TEST(Skinning, createSidesForShellQuad8Block)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  std::unique_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
//shell-quad-8 mesh:
//       6         12
// 3*----*----*9----*----*15
//  |         |          |
//  |         |          |
// 2*         *8         *14
//  |         |          |
//  |         |          |
// 1*----*----*7----*----*13
//       4         10
//
  stk::mesh::Part& skinPart = bulk->mesh_meta_data().declare_part("mySkin");
  const std::string meshDesc =
       "0,1,SHELL_QUAD_8, 1,7,9,3,4,8,6,2, block_1\n\
        0,2,SHELL_QUAD_8, 7,13,15,9,10,14,12,8, block_2|sideset:name=surface_1; data=1,1,1,3,2,3,2,6";

  std::vector<double> coords = {
     /*1*/0,0,0,  /*2*/0,1,0,  /*3*/0,2,0,
     /*4*/1,0,0,  /*6*/1,2,0,
     /*7*/2,0,0,  /*8*/2,1,0,  /*9*/2,2,0,
     /*10*/3,0,0,  /*12*/3,2,0,
     /*13*/4,0,0,  /*14*/4,1,0,  /*15*/4,2,0};

  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

  stk::mesh::create_exposed_block_boundary_sides(*bulk, bulk->mesh_meta_data().universal_part(), stk::mesh::PartVector{&skinPart});
  EXPECT_EQ(4u, stk::mesh::count_entities(*bulk, stk::topology::FACE_RANK, skinPart));
}

TEST(Skinning, createAllSidesForBlock)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                            .set_spatial_dimension(3)
                                                            .create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::mesh::Part& block2 = meta.declare_part_with_topology("block_2", stk::topology::HEX_8);
  stk::mesh::Part& surface1 = meta.declare_part_with_topology("surface_1", stk::topology::QUAD_4);
  stk::io::put_io_part_attribute(block2);
  stk::io::put_io_part_attribute(surface1);
  meta.set_part_id(block2, 2);
  meta.set_part_id(surface1, 1);

  stk::io::fill_mesh("generated:2x2x2", bulk);

  copy_elems_from_block_to_block(bulk, {1, 2}, "block_1", "block_2");

  stk::mesh::create_all_block_boundary_sides(bulk, meta.universal_part(), stk::mesh::PartVector{&surface1});
  EXPECT_EQ(28u, stk::mesh::count_entities(bulk, stk::topology::FACE_RANK, surface1));
}

TEST(Skinning, createAllSidesForBlock_separatePartForInteriorSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                            .set_spatial_dimension(3)
                                                            .create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::mesh::Part& block2 = meta.declare_part_with_topology("block_2", stk::topology::HEX_8);
  stk::mesh::Part& surface1 = meta.declare_part_with_topology("surface_1", stk::topology::QUAD_4);
  stk::mesh::Part& interiorSkin = meta.declare_part_with_topology("interiorSkin", stk::topology::QUAD_4);
  stk::io::put_io_part_attribute(block2);
  stk::io::put_io_part_attribute(surface1);
  stk::io::put_io_part_attribute(interiorSkin);
  meta.set_part_id(block2, 2);
  meta.set_part_id(surface1, 1);
  meta.set_part_id(interiorSkin, 2);

  stk::io::fill_mesh("generated:2x2x2", bulk);

  copy_elems_from_block_to_block(bulk, {1, 2}, "block_1", "block_2");

  stk::mesh::PartVector interiorSkinParts = {&interiorSkin};
  stk::mesh::create_all_block_boundary_sides(bulk, meta.universal_part(), stk::mesh::PartVector{&surface1}, &interiorSkinParts);
  EXPECT_EQ(24u, stk::mesh::count_entities(bulk, stk::topology::FACE_RANK, surface1));
  EXPECT_EQ(4u, stk::mesh::count_entities(bulk, stk::topology::FACE_RANK, interiorSkin));
}

