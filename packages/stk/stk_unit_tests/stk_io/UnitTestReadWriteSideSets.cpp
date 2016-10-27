#include <gtest/gtest.h>
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_unit_test_utils/ReadWriteSidesetTester.hpp"

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
        for(std::string& input_file_name : filesToTest)
            EXPECT_NO_THROW(stk::unit_test_util::sideset::test_reading_writing_sideset_from_file(comm, input_file_name, "new.e")) << " for file " << input_file_name;
}

TEST(StkIo, DISABLED_read_write_and_compare_exo_files_with_sidesets_because_PMR_for_coincident_not_implemented_yet)
{
    std::vector<std::string> filesToTest = {
                                            "ALefRA.e"
    };

    stk::ParallelMachine comm = MPI_COMM_WORLD;
    if ( stk::parallel_machine_size(comm) == 2)
        for(std::string& input_file_name : filesToTest)
            EXPECT_NO_THROW(stk::unit_test_util::sideset::test_reading_writing_sideset_from_file(comm, input_file_name, "new.e")) << " for file " << input_file_name;
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

        bulk.initialize_face_adjacent_element_graph();

        bulk.modification_begin();
        unsigned side_ord = side_ordinal_on_proc[bulk.parallel_rank()];
        stk::mesh::PartVector parts={part_on_proc[bulk.parallel_rank()]};
        stk::mesh::Entity side = bulk.declare_element_side(element, side_ord, parts);

        if(bulk.parallel_rank()==0)
            EXPECT_TRUE(does_entity_have_part(bulk, side, partProc0));
        else
            EXPECT_TRUE(does_entity_have_part(bulk, side, partProc1));

        bulk.modification_end();

        EXPECT_TRUE(bulk.parallel_owner_rank(side)==0);
        EXPECT_TRUE(does_entity_have_part(bulk, side, partProc0));
        EXPECT_TRUE(does_entity_have_part(bulk, side, partProc1));
    }
}

void test_serial_write_new_sideset(stk::unit_test_util::sideset::BulkDataTester &bulk,
                                   stk::mesh::PartVector &parts,
                                   const std::string &outputFileName)
{
  bulk.initialize_face_adjacent_element_graph();

  stk::mesh::Entity element1 = bulk.get_entity(stk::topology::ELEMENT_RANK, 1u);
  EXPECT_TRUE(bulk.is_valid(element1));

  stk::mesh::Entity element2 = bulk.get_entity(stk::topology::ELEMENT_RANK, 2u);
  EXPECT_TRUE(bulk.is_valid(element2));

  const stk::unit_test_util::sideset::SideSetData newSidesetData = {{1, {{element1, 5}, {element2, 1}}}};
  bulk.create_side_entities(newSidesetData[0].sideSet, parts);

  stk::mesh::EntityVector faces;
  stk::mesh::get_selected_entities(bulk.mesh_meta_data().locally_owned_part(), bulk.buckets(bulk.mesh_meta_data().side_rank()), faces);
  EXPECT_EQ(faces.size(), 1u);
  ASSERT_EQ(2u, bulk.num_elements(faces[0]));

  stk::unit_test_util::sideset::write_exo_file(bulk, outputFileName, newSidesetData);
}

TEST(StkIo, create_and_write_new_sideset)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned p_size = stk::parallel_machine_size(pm);
  if(p_size == 1)
  {
      stk::mesh::MetaData meta;
      stk::unit_test_util::sideset::BulkDataTester bulk(meta, pm);

      stk::unit_test_util::sideset::StkMeshIoBrokerTester stkIo;
      stk::unit_test_util::sideset::setup_io_broker_for_read(stkIo, bulk, "Ae.e", stk::unit_test_util::sideset::READ_SERIAL_AND_DECOMPOSE);

      stk::mesh::PartVector parts;
      stk::mesh::Part &part = meta.declare_part("surface_1", meta.side_rank());
      meta.set_part_id(part, 1);
      stk::io::put_io_part_attribute(part);
      parts.push_back(&part);

      stk::unit_test_util::sideset::SideSetData sidesetData;
      stkIo.populate_bulk_data();
      stkIo.fill_sideset_data(sidesetData);
      EXPECT_TRUE(sidesetData.empty());

      test_serial_write_new_sideset(bulk, parts, "Ae_with_ss.e");
  }
}
