#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t, nullptr
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD
#include <stk_io/DatabasePurpose.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkIoUtils.hpp>
#include <Ioss_SideSet.h>
#include <Ioss_SideBlock.h>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/SideSetUtil.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"

using stk::unit_test_util::build_mesh_no_simple_fields;

namespace
{

stk::mesh::EntityVector get_sides(const stk::mesh::BulkData &bulkData, const stk::mesh::Part& sidesetPart)
{
    const stk::mesh::MetaData &meta = bulkData.mesh_meta_data();
    stk::mesh::EntityVector sides;
    stk::mesh::Selector sideSelector = sidesetPart & ( meta.locally_owned_part() | meta.globally_shared_part());
    stk::mesh::get_selected_entities(sideSelector, bulkData.buckets(meta.side_rank()), sides);
    return sides;
}

bool is_sideset_case_supported(const std::string& input_file_name, stk::mesh::BulkData::AutomaticAuraOption auraOption)
{
    bool supported = true;
    if(stk::parallel_machine_size(MPI_COMM_WORLD) < 3)
    {
        stk::ParallelMachine comm = MPI_COMM_WORLD;

        std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh_no_simple_fields(comm, auraOption);
        stk::mesh::MetaData& meta = bulk->mesh_meta_data();
        stk::io::fill_mesh_with_auto_decomp(input_file_name, *bulk);

        const stk::mesh::ElemElemGraph &graph = bulk->get_face_adjacent_element_graph();
        stk::mesh::impl::ParallelPartInfo parallelPartInfo;
        stk::mesh::impl::populate_part_ordinals_for_remote_edges(*bulk, graph, parallelPartInfo);

        const stk::mesh::PartVector& allparts = meta.get_mesh_parts();
        for(auto sidesetPart : allparts)
        {
            if(stk::mesh::is_side_set(*sidesetPart))
                supported = supported & stk::mesh::does_not_contain_internal_sideset(*bulk, get_sides(*bulk, *sidesetPart), parallelPartInfo);
        }
    }
    return supported;
}

void test_supported_sideset_cases_with_aura_option(stk::mesh::BulkData::AutomaticAuraOption auraOption)
{
    std::string exodusFileName = stk::unit_test_util::get_option("-i", "none");

    if(exodusFileName=="none")
    {
        if(stk::parallel_machine_size(MPI_COMM_WORLD) <= 2)
        {
            std::vector<std::string> files ={"ADA.e", "ARA.e", "ALA.e"};
            for(std::string& filename : files)
                EXPECT_FALSE(is_sideset_case_supported(filename, auraOption));
        }
    }
    else
    {
        if(is_sideset_case_supported(exodusFileName, auraOption))
            std::cerr << "Sidesets in file are OK\n";
        else
            std::cerr << "Sidesets in file are not supported\n";
    }
}

TEST(StkIo, unsupported_cases_with_aura)
{
    test_supported_sideset_cases_with_aura_option(stk::mesh::BulkData::AUTO_AURA);
}

TEST(StkIo, unsupported_cases_without_aura)
{
    test_supported_sideset_cases_with_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
}

}
