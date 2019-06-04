#include <gtest/gtest.h>
#include <stk_unit_test_utils/ioUtils.hpp>
#include "stk_mesh/base/GetEntities.hpp"
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include "stk_unit_test_utils/ReadWriteSidesetTester.hpp"
#include "stk_unit_test_utils/FaceTestingUtils.hpp"

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
        stk::mesh::MetaData meta(3);
        stk::mesh::Part& partProc0 = meta.declare_part("partProc0", meta.side_rank());
        stk::mesh::Part& partProc1 = meta.declare_part("partProc1", meta.side_rank());
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
        stk::io::fill_mesh("generated:1x1x2", bulk);

        stk::mesh::EntityId elementId = 1+bulk.parallel_rank();
        stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEM_RANK, elementId);

        std::vector<unsigned> side_ordinal_on_proc = { 5, 4 };
        stk::mesh::PartVector part_on_proc = {&partProc0, &partProc1};

        bulk.modification_begin();
        unsigned side_ord = side_ordinal_on_proc[bulk.parallel_rank()];
        stk::mesh::PartVector parts={part_on_proc[bulk.parallel_rank()]};
        stk::mesh::Entity side = bulk.declare_element_side(element, side_ord, parts);

        if(bulk.parallel_rank()==0)
        {
            EXPECT_TRUE(does_entity_have_part(bulk, side, partProc0));
        }
        else
        {
            EXPECT_TRUE(does_entity_have_part(bulk, side, partProc1));
        }

        bulk.modification_end();

        EXPECT_TRUE(bulk.parallel_owner_rank(side)==0);
        EXPECT_TRUE(does_entity_have_part(bulk, side, partProc0));
        EXPECT_TRUE(does_entity_have_part(bulk, side, partProc1));
    }
}

void test_output_sideset(stk::unit_test_util::sideset::BulkDataTester &bulk,
                         const std::string &outputFileName,
                         stk::unit_test_util::sideset::ReadMode readMode)
{
  stk::unit_test_util::sideset::write_exo_file(bulk, outputFileName);

  stk::mesh::MetaData meta2;
  stk::unit_test_util::sideset::BulkDataTester bulk2(meta2, bulk.parallel());
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
  stk::mesh::SideSet &sideSet = bulk.create_sideset(stk::io::get_sideset_parent(*parts[0]));

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

void test_create_and_write_new_sideset(stk::ParallelMachine pm, const std::string &inputFileName, const std::string &outputFileName, const stk::unit_test_util::sideset::ElemIdSideVector &newSideSet)
{
  stk::mesh::MetaData meta(3);
  stk::unit_test_util::sideset::BulkDataTester bulk(meta, pm);
  stk::mesh::PartVector parts = create_parts(meta);

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
    stk::mesh::SideSet deletedSideset = stk::unit_test_util::sideset::get_stk_side_set(bulk, deletedElemIdSides);
    stk::mesh::Part *surface_part = stk::unit_test_util::get_surface_part_with_id(bulk.mesh_meta_data(), inputId);
    ThrowRequire(nullptr != surface_part);
    stk::mesh::SideSet& sideSet = bulk.get_sideset(*surface_part);
    for(const stk::mesh::SideSetEntry &entry : deletedSideset)
    {
        auto iter = std::find(sideSet.begin(), sideSet.end(), entry);
        EXPECT_TRUE(iter != sideSet.end());
        sideSet.erase(iter);
    }
}

void test_read_and_modify_sideset(stk::ParallelMachine pm,
                                  const std::string &inputFileName,
                                  const std::string &outputFileName,
                                  const stk::unit_test_util::sideset::ElemIdSideVector &expectedInputElemIdSides,
                                  const stk::unit_test_util::sideset::ElemIdSideVector &deletedElemIdSides,
                                  const stk::unit_test_util::sideset::ElemIdSideVector &expectedOutputElemIdSides)
{
    stk::mesh::MetaData meta;
    stk::unit_test_util::sideset::BulkDataTester bulk(meta, pm);
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

TEST(StkIo, parallel_transform_AA_to_ADA_to_ARA)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned p_size = stk::parallel_machine_size(pm);

  if(p_size == 2)
  {
      stk::mesh::MetaData meta(3);
      stk::unit_test_util::sideset::BulkDataTester bulk(meta, pm);

      stk::mesh::PartVector parts = create_parts(meta);
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
