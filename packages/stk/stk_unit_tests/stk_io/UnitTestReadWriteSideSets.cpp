#include <vector>

#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/InputFile.hpp>   // for InputFile for m_input_files
#include <stk_util/parallel/Parallel.hpp> // for stk::parallel_machine_size
#include "Ioss_Region.h"                // for Region, NodeSetContainer, etc
#include "Ioss_SideBlock.h"             // for SideBlock
#include "Ioss_SideSet.h"               // for SideSet, SideBlockContainer
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"

struct IdAndSideSet {
    int id;
    stk::mesh::SideSet sideSet;
};

class StkMeshIoBrokerTester : public stk::io::StkMeshIoBroker
{
public:
    StkMeshIoBrokerTester() {}

    virtual void populate_mesh(bool delay_field_data_allocation = true)
    {
        stk::io::StkMeshIoBroker::populate_mesh(delay_field_data_allocation);
        extract_sideset_data_from_io(sideset_data);
    }

private:
    int convert_to_zero_based_ordinal(int one_based_ordinal) const
    {
        return one_based_ordinal-1;
    }

    stk::mesh::ElemIdSide add_elem_side_pair(int elem_id, int one_based_ordinal)
    {
        int zero_based_side_ordinal = convert_to_zero_based_ordinal(one_based_ordinal);
        return stk::mesh::ElemIdSide{elem_id, zero_based_side_ordinal};
    }

    void convert_elem_sides_pairs_into_sideset(const stk::mesh::BulkData& bulk, const std::vector<int>& elem_side, stk::mesh::SideSet& sideset)
    {
        for(size_t is=0; is<elem_side.size() / 2; ++is)
        {
            stk::mesh::Entity const elem = bulk.get_entity(stk::topology::ELEMENT_RANK, elem_side[is*2]);
            if (bulk.is_valid(elem))
                sideset.push_back(add_elem_side_pair(elem_side[is*2], elem_side[is*2+1]));
        }
    }

    void convert_block_to_sideset(const stk::mesh::BulkData& bulk, const Ioss::SideBlock* block, stk::mesh::SideSet& sideset)
    {
        if (stk::io::include_entity(block))
        {
            std::vector<int> elem_side ;
            block->get_field_data("element_side", elem_side);
            convert_elem_sides_pairs_into_sideset(bulk, elem_side, sideset);
        }
    }

    void convert_ioss_sideset_to_stk_sideset(const stk::mesh::BulkData& bulk, const Ioss::SideSet* sset, stk::mesh::SideSet& sideset)
    {
        if(stk::io::include_entity(sset))
        {
            for (size_t i=0; i < sset->block_count(); i++)
            {
                Ioss::SideBlock *block = sset->get_block(i);
                convert_block_to_sideset(bulk, block, sideset);
            }
        }
    }

    void add_ioss_sideset_to_stk_sideset_with_sideset_id(const stk::mesh::BulkData& bulk, const Ioss::SideSet* sset, std::vector<IdAndSideSet>& sideset_data)
    {
        int sset_id = sset->get_property("id").get_int();
        sideset_data.push_back(IdAndSideSet{sset_id, stk::mesh::SideSet{}});
        stk::mesh::SideSet& currentSideset = sideset_data.back().sideSet;
        convert_ioss_sideset_to_stk_sideset(bulk, sset, currentSideset);
    }

public:
    void fill_sideset_data(std::vector<IdAndSideSet>& sidesets) const
    {
        sidesets = sideset_data;
    }

private:
    void extract_sideset_data_from_io(std::vector<IdAndSideSet> &sidesetData)
    {
        stk::mesh::BulkData &bulk = this->bulk_data();
        Ioss::Region *region = m_input_files[m_active_mesh_index]->get_input_io_region().get();
        for ( const Ioss::SideSet * sset : region->get_sidesets() )
            add_ioss_sideset_to_stk_sideset_with_sideset_id(bulk, sset, sidesetData);
    }

    std::vector<IdAndSideSet> sideset_data;
};

class BulkDataTester : public stk::mesh::BulkData
{
public:
    BulkDataTester(stk::mesh::MetaData & mesh_meta_data
                   , stk::ParallelMachine parallel)
    :   stk::mesh::BulkData(mesh_meta_data, parallel)
    {
    }

    void tester_save_sideset_data(int sideset_id, const stk::mesh::SideSet& data)
    {
        save_sideset_data(sideset_id, data);
    }
};

enum ReadMode { READ_SERIAL_AND_DECOMPOSE, READ_ALREADY_DECOMPOSED };

void read_exo_file( stk::mesh::BulkData &bulkData, std::string filename, std::vector<IdAndSideSet> &sideset_data, ReadMode read_mode)
{
    StkMeshIoBrokerTester stkIo;
    if(read_mode == READ_SERIAL_AND_DECOMPOSE)
        stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
    stkIo.set_bulk_data(bulkData);
    stkIo.add_mesh_database(filename, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.add_all_mesh_fields_as_input_fields();
    stkIo.populate_bulk_data();

    stkIo.fill_sideset_data(sideset_data);
}

void write_exo_file(BulkDataTester &bulkData, const std::string &filename, const std::vector<IdAndSideSet>& sideset_data)
{
    for(const IdAndSideSet& sset : sideset_data)
        bulkData.tester_save_sideset_data(sset.id, sset.sideSet);

    StkMeshIoBrokerTester io_broker;
    io_broker.set_bulk_data(bulkData);
    size_t resultFileIndex = io_broker.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    io_broker.write_output_mesh(resultFileIndex);
}

void fill_sideset_data_from_serial_input_file_and_write_decomposed_file(stk::ParallelMachine comm, const std::string& input_filename, const std::string& output_file_name, std::vector<IdAndSideSet>& sideset_data)
{
    stk::mesh::MetaData meta;
    BulkDataTester bulkData(meta, comm);
    read_exo_file( bulkData, input_filename, sideset_data, READ_SERIAL_AND_DECOMPOSE);
    write_exo_file( bulkData, output_file_name, sideset_data);
}

void fill_sideset_data_from_decomposed_input_file(stk::ParallelMachine comm, const std::string& input_filename, std::vector<IdAndSideSet>& sideset_data)
{
    stk::mesh::MetaData meta;
    BulkDataTester bulkData(meta, comm);
    read_exo_file( bulkData, input_filename, sideset_data, READ_ALREADY_DECOMPOSED);
}

void compare_sideset_read_against_sideset_written(const std::string& input_file_name, const std::vector<IdAndSideSet> &sideset_data1, const std::vector<IdAndSideSet> &sideset_data2)
{
    ASSERT_EQ(sideset_data1.size(), sideset_data2.size()) << "for file: " << input_file_name;
    for(size_t ss=0; ss<sideset_data1.size(); ++ss)
    {
        const stk::mesh::SideSet& sideSet1 = sideset_data1[ss].sideSet;
        const stk::mesh::SideSet& sideSet2 = sideset_data2[ss].sideSet;
        EXPECT_EQ(sideset_data1[ss].id, sideset_data2[ss].id);
        ASSERT_EQ(sideSet1.size(), sideSet2.size()) << "for file: " << input_file_name;

        for(size_t i=0;i<sideSet1.size();++i)
        {
            EXPECT_EQ(sideSet1[i].elem_id, sideSet2[i].elem_id) << "for file: " << input_file_name;
            EXPECT_EQ(sideSet1[i].side_ordinal, sideSet2[i].side_ordinal) << "for file: " << input_file_name;
        }
    }
}

std::vector<IdAndSideSet> get_sideset_data_from_file_read_and_write_new_file(stk::ParallelMachine comm, const std::string& input_file_name, const std::string& output_file_name)
{
    std::vector<IdAndSideSet> sideset_data1;
    fill_sideset_data_from_serial_input_file_and_write_decomposed_file(comm, input_file_name, output_file_name, sideset_data1);
    return sideset_data1;
}

std::vector<IdAndSideSet> get_sideset_data_from_written_file(stk::ParallelMachine comm, const std::string& output_file_name)
{
    std::vector<IdAndSideSet> sideset_data2;
    fill_sideset_data_from_decomposed_input_file(comm, output_file_name, sideset_data2);
    return sideset_data2;
}

void test_reading_writing_sideset_from_file(stk::ParallelMachine comm, const std::string& input_file_name)
{
    std::string output_file_name = "new.e";
    std::vector<IdAndSideSet> sideset_data1 = get_sideset_data_from_file_read_and_write_new_file(comm, input_file_name, output_file_name);
    std::vector<IdAndSideSet> sideset_data2 = get_sideset_data_from_written_file(comm, output_file_name);
    compare_sideset_read_against_sideset_written(input_file_name, sideset_data1, sideset_data2);
    unlink(output_file_name.c_str());
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
        for(std::string& input_file_name : filesToTest)
            EXPECT_NO_THROW(test_reading_writing_sideset_from_file(comm, input_file_name)) << " for file " << input_file_name;
}

TEST(StkIo, DISABLED_read_write_and_compare_exo_files_with_sidesets_because_PMR_for_coincident_not_implemented_yet)
{
    std::vector<std::string> filesToTest = {
                                            "ALefRA.e"
    };

    stk::ParallelMachine comm = MPI_COMM_WORLD;
    if ( stk::parallel_machine_size(comm) == 2)
        for(std::string& input_file_name : filesToTest)
            EXPECT_NO_THROW(test_reading_writing_sideset_from_file(comm, input_file_name)) << " for file " << input_file_name;
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
