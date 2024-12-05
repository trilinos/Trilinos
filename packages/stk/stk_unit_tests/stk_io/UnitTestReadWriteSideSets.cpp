#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include "stk_mesh/base/GetEntities.hpp"
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/SideSetUtil.hpp>
#include "stk_unit_test_utils/ReadWriteSidesetTester.hpp"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include "stk_unit_test_utils/FaceTestingUtils.hpp"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "IOMeshFixture.hpp"
#include <stk_unit_test_utils/BuildMesh.hpp>

namespace stk
{
namespace io
{
namespace unit_test
{

using stk::unit_test_util::build_mesh;
using stk::unit_test_util::build_mesh_no_simple_fields;

class StkIoSubset : public stk::io::unit_test::IOMeshFixture
{
protected:
  void test_write_then_read(const std::vector<size_t>& expectedEntityCounts,
                            stk::mesh::Part* blockToExclude = nullptr)
  {
    const std::string fileName("meshSubset.g");
    stk::mesh::Selector meshSubsetSelector = create_block_subset_selector({blockToExclude});
    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(get_bulk());
    size_t outputFileIndex = stkIo.create_output_mesh(fileName, stk::io::WRITE_RESULTS);
    stkIo.set_output_selector(outputFileIndex, stk::topology::ELEM_RANK, meshSubsetSelector);
    stkIo.write_output_mesh(outputFileIndex);

    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(MPI_COMM_WORLD);
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::io::fill_mesh(fileName, *bulk);

    std::vector<size_t> entityCounts;
    stk::mesh::count_entities(meta.locally_owned_part(), *bulk, entityCounts);

    ASSERT_TRUE(entityCounts.size() <= expectedEntityCounts.size());
    for(size_t i=0; i<entityCounts.size(); ++i) {
       EXPECT_EQ(entityCounts[i], expectedEntityCounts[i]);
    }
    unlink(fileName.c_str());
  }
};

TEST_F(StkIoSubset, outputOneOfTwoBlocks)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { return; }
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  const std::vector<std::string> partNames {"block_1", "block_2"};
  stk::mesh::Part& block1 = create_io_part(partNames[0], 1);
  stk::mesh::Part& block2 = create_io_part(partNames[1], 2);
  stk::mesh::Part& surface1 = create_io_part("surface_1", 1, stk::topology::QUAD_4);
  stk::io::fill_mesh("generated:1x1x2", get_bulk());
  stk::mesh::EntityId elemId = 1;
  move_element(elemId, block1, block2);
  stk::mesh::ConnectivityOrdinal sideOrd = 0;
  create_side(elemId, sideOrd, surface1);
  elemId = 2;
  create_side(elemId, sideOrd, surface1);

  std::vector<size_t> entityCounts(get_meta().entity_rank_count(), 0);
  entityCounts[stk::topology::NODE_RANK] = 8;
  entityCounts[stk::topology::FACE_RANK] = 1;
  entityCounts[stk::topology::ELEM_RANK] = 1;

  test_write_then_read(entityCounts, &block1);
}

class StkIoSideset : public stk::io::unit_test::IOMeshFixture
{
protected:
  using VecField = stk::mesh::Field<double>;

  void set_face_field_data(VecField& ssField, stk::mesh::Entity face)
  {
    const stk::mesh::Entity* faceNodes = get_bulk().begin_nodes(face);
    unsigned numFaceNodes = get_bulk().num_nodes(face);
    EXPECT_EQ(numFaceNodes, stk::mesh::field_scalars_per_entity(ssField, face));
    double* fieldData = stk::mesh::field_data(ssField, face);
    for (unsigned n=0; n<numFaceNodes; ++n) {
      fieldData[n] = static_cast<double>(get_bulk().identifier(faceNodes[n]));
    }
  }

  void verify_face_field_sizes_and_values(const stk::mesh::BulkData& bulk,
                                          const std::string& ssFieldName,
                                          unsigned expectedNumFaces)
  {
    const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    stk::mesh::Part* surface1 = meta.get_part("surface_1");
    ASSERT_TRUE(surface1 != nullptr);

    stk::mesh::FieldBase* ssField = meta.get_field(meta.side_rank(), ssFieldName);
    ASSERT_TRUE(ssField != nullptr);

    stk::mesh::Selector selector(*ssField & meta.locally_owned_part());
    stk::mesh::EntityVector faces;
    stk::mesh::get_entities(bulk, meta.side_rank(), selector, faces);
    EXPECT_EQ(expectedNumFaces, faces.size());

    for(stk::mesh::Entity face : faces) {
      unsigned numFaceNodes = bulk.num_nodes(face);
      unsigned numScalars = stk::mesh::field_scalars_per_entity(*ssField, face);
      ASSERT_TRUE(numFaceNodes <= numScalars);
      const stk::mesh::Entity* faceNodes = bulk.begin_nodes(face);
      const double* fieldData = static_cast<const double*>(stk::mesh::field_data(*ssField, face));
      for(unsigned n=0; n<numFaceNodes; ++n) {
        EXPECT_NEAR(fieldData[n], static_cast<double>(bulk.identifier(faceNodes[n])), 1.e-6);
      }
    }
  }

  void confirm_fmwk_restart_sequence(const std::string& fileName,
                                     const std::string& ssFieldName,
                                     unsigned expectedNumFaces)
  {
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(MPI_COMM_WORLD);
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(*bulk);
    stkIo.add_mesh_database(fileName, stk::io::READ_RESTART);
    stkIo.create_input_mesh();
    stkIo.add_all_mesh_fields_as_input_fields();

    stk::mesh::FieldBase* sidesetField = meta.get_field(meta.side_rank(), ssFieldName);
    ASSERT_TRUE(sidesetField != nullptr);
    stkIo.populate_bulk_data();

    stk::io::MeshField meshField(*sidesetField, sidesetField->name());
    const double time = 1.0;
    meshField.set_read_time(time);
    stkIo.read_input_field(meshField);

    verify_face_field_sizes_and_values(*bulk, ssFieldName, expectedNumFaces);
  }

  void read_and_verify_face_field_sizes_and_values(const std::string& fileName,
                                                   const std::string& ssFieldName,
                                                   unsigned expectedNumFaces,
                                                   bool isRestart = false)
  {
    stk::io::DatabasePurpose dbReadPurpose = isRestart ? stk::io::READ_RESTART : stk::io::READ_MESH;
    {
      std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(get_comm());

      stk::io::fill_mesh_with_fields(fileName, *bulk, dbReadPurpose);

      verify_face_field_sizes_and_values(*bulk, ssFieldName, expectedNumFaces);
    }

    if (isRestart) {
      EXPECT_NO_THROW(confirm_fmwk_restart_sequence(fileName, ssFieldName, expectedNumFaces));
    }
  }

  void set_coords(stk::mesh::BulkData& bulk, VecField& coordField,
                  const stk::mesh::EntityIdVector& nodeIds,
                  const std::vector<double>& coords)
  {
    const unsigned spatialDim = bulk.mesh_meta_data().spatial_dimension();
    ASSERT_EQ(coords.size(), spatialDim*nodeIds.size());
    unsigned offset = 0;
    for(stk::mesh::EntityId id : nodeIds) {
      stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, id);
      ASSERT_TRUE(bulk.is_valid(node));
      double* coordData = stk::mesh::field_data(coordField, node);
      for(unsigned d=0; d<spatialDim; ++d) {
        coordData[d] = coords[offset++];
      }
    }
  }

  void create_parts_and_fields(bool createQuadBlockFirst = true)
  {
    const std::vector<std::string> partNames {"block_1", "block_2"};
    block1 = &create_io_part(partNames[0], 1, stk::topology::HEX_8);
    block2 = &create_io_part(partNames[1], 2, stk::topology::TET_4);
    surface1 = &get_meta().declare_part("surface_1", stk::topology::FACE_RANK);
    stk::io::put_io_part_attribute(*surface1);
    get_meta().set_part_id(*surface1, 1);

    if (createQuadBlockFirst) {
      sideBlock1 = &create_io_part("surface_hex8_quad4_1", 1, stk::topology::QUAD_4);
      get_meta().declare_part_subset(*surface1, *sideBlock1);
      sideBlock2 = &create_io_part("surface_tet4_tri3_1", 1, stk::topology::TRI_3);
      get_meta().declare_part_subset(*surface1, *sideBlock2);
    }
    else {
      sideBlock2 = &create_io_part("surface_tet4_tri3_1", 1, stk::topology::TRI_3);
      get_meta().declare_part_subset(*surface1, *sideBlock2);
      sideBlock1 = &create_io_part("surface_hex8_quad4_1", 1, stk::topology::QUAD_4);
      get_meta().declare_part_subset(*surface1, *sideBlock1);
    }

    coordField = &get_meta().declare_field<double>(stk::topology::NODE_RANK, "coordinates");
    get_meta().set_coordinate_field(coordField);
    ssField = &get_meta().declare_field<double>(stk::topology::FACE_RANK, "ssfield");
    stk::mesh::put_field_on_mesh(*coordField, get_meta().universal_part(), 3, nullptr);
    stk::mesh::put_field_on_mesh(*ssField, *sideBlock1, 4, nullptr);
    stk::mesh::put_field_on_mesh(*ssField, *sideBlock2, 3, nullptr);
  }

  void create_hex_with_face()
  {
    stk::mesh::PartVector elem1Parts{block1};
    stk::mesh::EntityId elemId = 1;
    stk::mesh::EntityIdVector nodeIds{1, 2, 3, 4, 5, 6, 7, 8};
    stk::mesh::Entity elem1 = stk::mesh::declare_element(get_bulk(), elem1Parts, elemId, nodeIds);
    std::vector<double> coords{0,0,0, 1,0,0, 1,1,0, 0,1,0,
                               0,0,1, 1,0,1, 1,1,1, 0,1,1};
    set_coords(get_bulk(), *coordField, nodeIds, coords);

    stk::mesh::PartVector quadFaceParts{surface1, sideBlock1};
    stk::mesh::ConnectivityOrdinal faceOrd = 0;
    stk::mesh::Entity face11 = get_bulk().declare_element_side(elem1, faceOrd, quadFaceParts);
    EXPECT_EQ(stk::topology::QUAD_4, get_bulk().bucket(face11).topology());
    set_face_field_data(*ssField, face11);
  }

  void create_tet_with_face()
  {
    stk::mesh::PartVector elem2Parts{block2};
    stk::mesh::EntityId elemId = 2;
    stk::mesh::EntityIdVector nodeIds = {1, 2, 9, 10};
    std::vector<double> coords = {0,0,0, 1,0,0, 0.5,-0.5,-0.5, 0.5,-0.5,0};
    stk::mesh::Entity elem2 = stk::mesh::declare_element(get_bulk(), elem2Parts,
                                            elemId, nodeIds);
    set_coords(get_bulk(), *coordField, nodeIds, coords);

    stk::mesh::PartVector triFaceParts{surface1, sideBlock2};
    stk::mesh::ConnectivityOrdinal faceOrd = 0;
    stk::mesh::Entity face21 = get_bulk().declare_element_side(elem2, faceOrd, triFaceParts);
    EXPECT_EQ(stk::topology::TRI_3, get_bulk().bucket(face21).topology());
    set_face_field_data(*ssField, face21);
  }

  void setup_mesh_hex_tet_quad_tri()
  {
    get_bulk().modification_begin();

    if (get_bulk().parallel_size() == 1 || get_bulk().parallel_rank() == 0) {
      create_tet_with_face();
    }

    if (get_bulk().parallel_size() == 1 || get_bulk().parallel_rank() == 1) {
      create_hex_with_face();
    }

    get_bulk().modification_end();
  }

  void test_write_then_read(bool isRestart = false)
  {
    std::string fileName("ssfield.exo");
    int step = 1;
    double time = 1.0;
    stk::io::DatabasePurpose dbWritePurpose = isRestart ? stk::io::WRITE_RESTART : stk::io::WRITE_RESULTS;
    stk::io::write_mesh_with_fields(fileName, get_bulk(), step, time, dbWritePurpose);

    unsigned expectedNumFaces = get_bulk().parallel_size()==1 ? 2 : 1;
    read_and_verify_face_field_sizes_and_values(fileName, ssField->name(), expectedNumFaces, isRestart);
    unlink(fileName.c_str());
  }

private:
  stk::mesh::Part* block1;
  stk::mesh::Part* block2;
  stk::mesh::Part* surface1;
  stk::mesh::Part* sideBlock1;
  stk::mesh::Part* sideBlock2;
  VecField* coordField;
  VecField* ssField;
};

TEST_F(StkIoSideset, field_QuadAndTriSides)
{
  if (stk::parallel_machine_size(get_comm()) > 2) { GTEST_SKIP(); }
  unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);

  bool createQuadBlockFirst = true;
  create_parts_and_fields(createQuadBlockFirst);
  setup_mesh_hex_tet_quad_tri();

  test_write_then_read();
}

TEST_F(StkIoSideset, field_TriAndQuadSides)
{
  if (stk::parallel_machine_size(get_comm()) > 2) { GTEST_SKIP(); }
  unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);

  bool createQuadBlockFirst = false;
  create_parts_and_fields(createQuadBlockFirst);
  setup_mesh_hex_tet_quad_tri();

  test_write_then_read();
}

TEST_F(StkIoSideset, field_QuadAndTriSides_restart)
{
  if (stk::parallel_machine_size(get_comm()) > 2) { GTEST_SKIP(); }
  unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);

  bool createQuadBlockFirst = true;
  create_parts_and_fields(createQuadBlockFirst);
  setup_mesh_hex_tet_quad_tri();

  const bool isRestart = true;
  test_write_then_read(isRestart);
}

TEST_F(StkIoSideset, field_TriAndQuadSides_restart)
{
  if (stk::parallel_machine_size(get_comm()) > 2) { GTEST_SKIP(); }
  unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);

  bool createQuadBlockFirst = false;
  create_parts_and_fields(createQuadBlockFirst);
  setup_mesh_hex_tet_quad_tri();

  const bool isRestart = true;
  test_write_then_read(isRestart);
}

TEST(StkIo, read_write_and_compare_exo_files_with_sidesets)
{
    std::vector<std::string> filesToTest = {
                                            "AA.e", "ADeDB.e", "ADeLB.e", "ADReA.e", "AefA.e",  "AL.e",
                                            "ALe.e",    "ALeLB.e", "ALRA.e",  "ARA.e",   "ARReA.e",
                                            "ADA.e",    "ADe.e",   "ADeRA.e", "ADReB.e", "AeL.e",   "ALeDA.e",
                                            "ALeRA.e", "ALReA.e", "ARe.e",   "ARReB.e",
                                            "ADeDA.e",  "ADeLA.e", "ADeRB.e", "AeA.e",   "ALA.e",   "ALeDB.e",
                                            "ALeLA.e",  "ALeRB.e", "ALReB.e", "ARefLA.e"
    };

    stk::ParallelMachine comm = MPI_COMM_WORLD;
    if ( stk::parallel_machine_size(comm) == 2)
    {
        for(std::string& input_file_name : filesToTest)
            EXPECT_NO_THROW(stk::unit_test_util::sideset::test_reading_writing_sideset_from_file(comm, input_file_name, "new.e")) << " for file " << input_file_name;
    }
}

TEST(StkIo, read_write_and_compare_exo_files_with_sidesets_because_PMR_for_coincident_not_implemented_yet)
{
    std::vector<std::string> filesToTest = { "ALefRA.e" };

    stk::ParallelMachine comm = MPI_COMM_WORLD;
    if ( stk::parallel_machine_size(comm) == 2)
    {
        for(std::string& input_file_name : filesToTest)
            EXPECT_NO_THROW(stk::unit_test_util::sideset::test_reading_writing_sideset_from_file(comm, input_file_name, "new.e")) << " for file " << input_file_name;
    }
}

bool does_entity_have_part(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity, stk::mesh::Part& input_part)
{
    for(stk::mesh::Part* part : bulk.bucket(entity).supersets())
    {
        if(part->name()==input_part.name())
            return true;
    }
    return false;
}

TEST(StkIo, test_side_membership)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD)==2)
    {
        const unsigned spatialDim = 3;
        std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD);
        stk::mesh::MetaData& meta = bulk->mesh_meta_data();

        stk::mesh::Part& partProc0 = meta.declare_part("partProc0", meta.side_rank());
        stk::mesh::Part& partProc1 = meta.declare_part("partProc1", meta.side_rank());
        stk::io::fill_mesh("generated:1x1x2", *bulk);

        stk::mesh::EntityId elementId = 1+bulk->parallel_rank();
        stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, elementId);

        std::vector<unsigned> side_ordinal_on_proc = { 5, 4 };
        stk::mesh::PartVector part_on_proc = {&partProc0, &partProc1};

        bulk->modification_begin();
        unsigned side_ord = side_ordinal_on_proc[bulk->parallel_rank()];
        stk::mesh::PartVector parts={part_on_proc[bulk->parallel_rank()]};
        stk::mesh::Entity side = bulk->declare_element_side(element, side_ord, parts);

        if(bulk->parallel_rank()==0)
        {
            EXPECT_TRUE(does_entity_have_part(*bulk, side, partProc0));
        }
        else
        {
            EXPECT_TRUE(does_entity_have_part(*bulk, side, partProc1));
        }

        bulk->modification_end();

        EXPECT_TRUE(bulk->parallel_owner_rank(side)==0);
        EXPECT_TRUE(does_entity_have_part(*bulk, side, partProc0));
        EXPECT_TRUE(does_entity_have_part(*bulk, side, partProc1));
    }
}

void test_output_sideset(stk::unit_test_util::sideset::BulkDataTester &bulk,
                         const std::string &outputFileName,
                         stk::unit_test_util::sideset::ReadMode readMode)
{
  stk::unit_test_util::sideset::write_exo_file(bulk, outputFileName);

  auto meta2 = std::make_shared<stk::mesh::MetaData>();
  stk::unit_test_util::sideset::BulkDataTester bulk2(*meta2, bulk.parallel());
  stk::unit_test_util::sideset::read_exo_file(bulk2, outputFileName, readMode);
  EXPECT_NO_THROW(stk::unit_test_util::sideset::compare_sidesets(outputFileName, bulk, bulk2));
}

void test_mesh_with_single_face_connected_to_two_elements(stk::unit_test_util::sideset::BulkDataTester &bulk)
{
  stk::mesh::EntityVector faces;
  stk::mesh::get_selected_entities(bulk.mesh_meta_data().locally_owned_part(), bulk.buckets(bulk.mesh_meta_data().side_rank()), faces);
  EXPECT_EQ(faces.size(), 1u);
  ASSERT_EQ(2u, bulk.num_elements(faces[0]));
}

void create_new_sideset_and_faces(stk::unit_test_util::sideset::BulkDataTester &bulk,
                                  stk::mesh::PartVector &parts,
                                  const stk::unit_test_util::sideset::ElemIdSideVector &newSideSet)
{
  stk::mesh::SideSet &sideSet = bulk.create_sideset(stk::mesh::get_sideset_parent(*parts[0]));

  for(unsigned i=0; i<newSideSet.size(); ++i)
  {
    stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEMENT_RANK, static_cast<stk::mesh::EntityId>(newSideSet[i].elem_id));
    EXPECT_TRUE(bulk.is_valid(element));
    sideSet.add({element, static_cast<stk::mesh::ConnectivityOrdinal>(newSideSet[i].side_ordinal)});
  }

  bulk.create_side_entities(sideSet, parts);
}

void test_create_and_write_new_sideset(stk::unit_test_util::sideset::BulkDataTester &bulk,
                                       stk::mesh::PartVector &parts,
                                       const stk::unit_test_util::sideset::ElemIdSideVector &newSideSet,
                                       const std::string &outputFileName)
{
  create_new_sideset_and_faces(bulk, parts, newSideSet);
  test_mesh_with_single_face_connected_to_two_elements(bulk);
  test_output_sideset(bulk, outputFileName, stk::unit_test_util::sideset::READ_SERIAL_AND_DECOMPOSE);
}

stk::mesh::PartVector create_parts(stk::mesh::MetaData &meta)
{
  stk::mesh::PartVector parts;
  stk::mesh::Part &part = meta.declare_part("surface_1", meta.side_rank());
  meta.set_part_id(part, 1);
  stk::io::put_io_part_attribute(part);
  parts.push_back(&part);
  return parts;
}

void load_mesh_with_no_sidesets(stk::mesh::BulkData &bulk, const std::string &inputFileName)
{
  stk::unit_test_util::sideset::read_exo_file(bulk, inputFileName, stk::unit_test_util::sideset::READ_SERIAL_AND_DECOMPOSE);
  EXPECT_TRUE(bulk.get_number_of_sidesets() == 0);
}

void test_create_and_write_new_sideset(stk::ParallelMachine pm,
                                       const std::string &inputFileName,
                                       const std::string &outputFileName,
                                       const stk::unit_test_util::sideset::ElemIdSideVector &newSideSet)
{
  auto meta = std::make_shared<stk::mesh::MetaData>(3);
  stk::unit_test_util::sideset::BulkDataTester bulk(*meta, pm);
  stk::mesh::PartVector parts = create_parts(*meta);

  load_mesh_with_no_sidesets(bulk, inputFileName);
  test_create_and_write_new_sideset(bulk, parts, newSideSet, outputFileName);
}

TEST(StkIo, create_and_write_new_sideset)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned p_size = stk::parallel_machine_size(pm);
  if(p_size == 1)
  {
      test_create_and_write_new_sideset(pm, "Ae.e", "new_ADe.e", {{1, 5}, {2, 1}});
      test_create_and_write_new_sideset(pm, "Ae.e", "new_ALe.e", {{1, 5}        });
      test_create_and_write_new_sideset(pm, "Ae.e", "new_ARe.e", {        {2, 1}});

      test_create_and_write_new_sideset(pm, "AA.e", "new_ADA.e", {{1, 5}, {2, 4}});
      test_create_and_write_new_sideset(pm, "AA.e", "new_ALA.e", {{1, 5}        });
      test_create_and_write_new_sideset(pm, "AA.e", "new_ARA.e", {        {2, 4}});
  }
}

void read_and_test_preexisting_sidesets(stk::unit_test_util::sideset::BulkDataTester& bulk,
                                       const std::string& inputFileName,
                                       int inputId,
                                       const stk::unit_test_util::sideset::ElemIdSideVector& expectedInputElemIdSides)
{
    stk::unit_test_util::sideset::read_exo_file(bulk, inputFileName, stk::unit_test_util::sideset::READ_SERIAL_AND_DECOMPOSE);
    stk::unit_test_util::sideset::SideSetIdAndElemIdSidesVector expectedInputSideset;
    expectedInputSideset.push_back( {inputId, expectedInputElemIdSides});
    stk::unit_test_util::sideset::compare_sidesets(inputFileName, bulk, expectedInputSideset);
}

void delete_sides_from_sideset(stk::unit_test_util::sideset::BulkDataTester& bulk,
                               int inputId,
                               const stk::unit_test_util::sideset::ElemIdSideVector& deletedElemIdSides)
{
    stk::mesh::SideSet* deletedSideset = stk::unit_test_util::sideset::get_stk_side_set(bulk, deletedElemIdSides);
    stk::mesh::Part *surface_part = stk::unit_test_util::get_surface_part_with_id(bulk.mesh_meta_data(), inputId);
    STK_ThrowRequire(nullptr != surface_part);
    stk::mesh::SideSet& sideSet = bulk.get_sideset(*surface_part);
    for(const stk::mesh::SideSetEntry &entry : *deletedSideset)
    {
        auto iter = std::find(sideSet.begin(), sideSet.end(), entry);
        EXPECT_TRUE(iter != sideSet.end());
        sideSet.erase(iter);
    }
    delete deletedSideset;
}

void test_read_and_modify_sideset(stk::ParallelMachine pm,
                                  const std::string &inputFileName,
                                  const std::string &outputFileName,
                                  const stk::unit_test_util::sideset::ElemIdSideVector &expectedInputElemIdSides,
                                  const stk::unit_test_util::sideset::ElemIdSideVector &deletedElemIdSides,
                                  const stk::unit_test_util::sideset::ElemIdSideVector &expectedOutputElemIdSides)
{
    auto meta = std::make_shared<stk::mesh::MetaData>();
    stk::unit_test_util::sideset::BulkDataTester bulk(*meta, pm);
    int inputId = 1;

    read_and_test_preexisting_sidesets(bulk, inputFileName, inputId, expectedInputElemIdSides);
    delete_sides_from_sideset(bulk, inputId, deletedElemIdSides);
    EXPECT_NO_THROW(stk::unit_test_util::sideset::compare_sidesets(outputFileName, bulk, {{inputId, expectedOutputElemIdSides}}));

    test_output_sideset(bulk, outputFileName, stk::unit_test_util::sideset::READ_SERIAL_AND_DECOMPOSE);
}

TEST(StkIo, modify_sideset)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned p_size = stk::parallel_machine_size(pm);

  if(p_size == 1)
  {
      test_read_and_modify_sideset(pm, "ADe.e", "new_ARe.e", {{1, 5}, {2, 1}}, {{1, 5}        }, {{2, 1}});
      test_read_and_modify_sideset(pm, "ADe.e", "new_ALe.e", {{1, 5}, {2, 1}}, {        {2, 1}}, {{1, 5}});
      test_read_and_modify_sideset(pm, "ADe.e", "new_Ae.e",  {{1, 5}, {2, 1}}, {{1, 5}, {2, 1}}, {      });
  }
}


TEST(StkIo, skinned_sideset_from_badly_named_element_block)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(pm) > 1) GTEST_SKIP();

  std::string meshDesc = "textmesh: 0,1,HEX_8,1,2,3,4,5,6,7,8, shell_a"
                         "|coordinates: 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1"
                         "|sideset:name=surface_skinned_shell_a; data=1,1";

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh_no_simple_fields(spatialDim, pm, stk::mesh::BulkData::NO_AUTO_AURA);
  EXPECT_NO_THROW(stk::io::fill_mesh(meshDesc, *bulk));

  std::string fileName("missing_skinned_shell.g");
  stk::io::write_mesh(fileName, *bulk);

  std::shared_ptr<stk::mesh::BulkData> bulk2 = build_mesh_no_simple_fields(spatialDim, pm, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::unit_test_util::sideset::read_exo_file(*bulk2, fileName, stk::unit_test_util::sideset::READ_SERIAL_AND_DECOMPOSE);
  stk::unit_test_util::sideset::SideSetIdAndElemIdSidesVector expectedInputSideset{ {1, {{1, 0}}} };
  stk::unit_test_util::sideset::compare_sidesets(fileName, *bulk2, expectedInputSideset);
  
  unlink(fileName.c_str());
}

TEST(StkIo, parallel_transform_AA_to_ADA_to_ARA)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned p_size = stk::parallel_machine_size(pm);

  if(p_size == 2)
  {
      auto meta = std::make_shared<stk::mesh::MetaData>(3);
      stk::unit_test_util::sideset::BulkDataTester bulk(*meta, pm);

      stk::mesh::PartVector parts = create_parts(*meta);
      load_mesh_with_no_sidesets(bulk, "AA.e");

      stk::mesh::SideSet &sideSet = bulk.create_sideset(*parts[0]);

      if(bulk.parallel_rank() == 0)
          sideSet.add({bulk.get_entity(stk::topology::ELEMENT_RANK, 1u), 5});
      else
          sideSet.add({bulk.get_entity(stk::topology::ELEMENT_RANK, 2u), 4});

      bulk.create_side_entities(sideSet, parts);
      test_output_sideset(bulk,"modified_ADA.e", stk::unit_test_util::sideset::READ_ALREADY_DECOMPOSED);

      if(bulk.parallel_rank() == 0)
          sideSet.clear();

      test_output_sideset(bulk,"modified_ARA.e", stk::unit_test_util::sideset::READ_ALREADY_DECOMPOSED);
  }
}

class ShellSidesets : public stk::unit_test_util::MeshFixture
{
protected:

  stk::mesh::Part& create_io_part(const std::string& name,
                                  stk::topology topo,
                                  int id)
  {
    stk::mesh::Part& part = get_meta().declare_part_with_topology(name, topo);
    get_meta().set_part_id(part, id);
    stk::io::put_io_part_attribute(part);
    return part;
  }

  void create_parts_for_shell_and_sidesets()
  {
    get_meta().declare_field<double>(stk::topology::NODE_RANK, "coordinates");
    shellBlock = &create_io_part("block_1", stk::topology::SHELL_TRI_3, 1);

    surface1 = &create_io_part("surface_1", stk::topology::TRI_3, 1);
    ss1blk1 = &create_io_part("surface_trishell3_tri3_1", stk::topology::TRI_3, 1);
    get_meta().declare_part_subset(*surface1, *ss1blk1);

    surface2 = &create_io_part("surface_2", stk::topology::LINE_2, 2);
    ss2blk1 = &create_io_part("surface_trishell3_edge2_2", stk::topology::LINE_2, 2);
    get_meta().declare_part_subset(*surface2, *ss2blk1);

    stk::mesh::ConstPartVector blocks = {shellBlock};
    get_meta().set_surface_to_block_mapping(surface1, blocks);
    get_meta().set_surface_to_block_mapping(ss1blk1, blocks);
    get_meta().set_surface_to_block_mapping(surface2, blocks);
    get_meta().set_surface_to_block_mapping(ss2blk1, blocks);
  }

  void create_shell_and_sides()
  {
    get_bulk().modification_begin();

    stk::mesh::EntityId shellId = 1;
    stk::mesh::EntityIdVector nodeIds = {1, 2, 3};
    stk::mesh::Entity shell1 = stk::mesh::declare_element(get_bulk(), *shellBlock, shellId, nodeIds);

    stk::mesh::ConstPartVector triParts = {ss1blk1};
    get_bulk().declare_element_side(shell1, 0, triParts);

    stk::mesh::ConstPartVector edgeParts = {ss2blk1};
    stk::mesh::EntityId edgeId = 1;
    stk::mesh::Entity edge = get_bulk().declare_edge(edgeId, edgeParts);
    stk::mesh::Entity node1 = get_bulk().get_entity(stk::topology::NODE_RANK, 1);
    stk::mesh::Entity node2 = get_bulk().get_entity(stk::topology::NODE_RANK, 2);
    get_bulk().declare_relation(edge, node1, 0);
    get_bulk().declare_relation(edge, node2, 1);
    get_bulk().declare_relation(shell1, edge, 0);

    get_bulk().modification_end();
  }

  void test_write_then_read()
  {
    std::string fileName("shellss.g");
    stk::io::write_mesh(fileName, get_bulk());

    const unsigned spatialDim = 3;
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh_no_simple_fields(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::io::fill_mesh(fileName, *bulk);

    stk::mesh::Part* test_surface1 = meta.get_part("surface_1");
    ASSERT_TRUE(test_surface1 != nullptr);
    EXPECT_EQ(stk::topology::FACE_RANK, test_surface1->primary_entity_rank());

    stk::mesh::Part* test_surface2 = meta.get_part("surface_2");
    ASSERT_TRUE(test_surface2 != nullptr);
    EXPECT_EQ(stk::topology::FACE_RANK, test_surface2->primary_entity_rank());

    unlink(fileName.c_str());
  }

  stk::mesh::Part* shellBlock;
  stk::mesh::Part* surface1;
  stk::mesh::Part* ss1blk1;
  stk::mesh::Part* surface2;
  stk::mesh::Part* ss2blk1;
};

TEST_F(ShellSidesets, DISABLED_testWriteThenRead)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts_for_shell_and_sidesets();
  create_shell_and_sides();
  test_write_then_read();
}

void check_shell_edge_side(const stk::mesh::BulkData& bulk, stk::mesh::ConnectivityOrdinal sideOrdinal)
{
  stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 1);
  ASSERT_TRUE(bulk.is_valid(elem));
  EXPECT_TRUE(bulk.bucket(elem).topology().is_shell());
  const unsigned numEdges = bulk.num_connectivity(elem, stk::topology::EDGE_RANK);
  ASSERT_EQ(1u, numEdges);
  const stk::mesh::ConnectivityOrdinal* ords = bulk.begin_ordinals(elem, stk::topology::EDGE_RANK);
  const stk::mesh::ConnectivityOrdinal sideOrdinalOffset = 2;
  const stk::mesh::ConnectivityOrdinal expectedZeroBasedOrd = sideOrdinal - sideOrdinalOffset - 1;
  EXPECT_EQ(expectedZeroBasedOrd, ords[0]);
}

TEST(ShellSideSidesets, DISABLED_create3DEdgesFromShellSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "textmesh: 0,1,SHELL_QUAD_4,1,2,3,4 |coordinates: 0,0,0, 1,0,0, 1,1,0, 0,1,0 |sideset:name=surface_1; data=1,3";

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh_no_simple_fields(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  std::cout<<"reading mesh from textmesh"<<std::endl;
  EXPECT_NO_THROW(stk::io::fill_mesh(meshDesc, *bulk));

  std::cout<<"checking mesh created from textmesh:"<<std::endl;
  stk::mesh::ConnectivityOrdinal expectedSideOrdinal = 3;
  check_shell_edge_side(*bulk, expectedSideOrdinal);

  std::string fileName("shell_edge_side.g");
  std::cout<<"writing mesh"<<std::endl;
  stk::io::write_mesh(fileName, *bulk);
 
  std::shared_ptr<stk::mesh::BulkData> bulk2 = build_mesh_no_simple_fields(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  std::cout<<"reading mesh from file"<<std::endl;
  EXPECT_NO_THROW(stk::io::fill_mesh(fileName, *bulk2));

  std::cout<<"checking mesh written-then-read from stk-io:"<<std::endl;
  check_shell_edge_side(*bulk2, expectedSideOrdinal);
}

}  // namespace unit_test
}  // namespace io
}  // namespace stk
