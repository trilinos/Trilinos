#include <stk_balance/balance.hpp>
#include <gtest/gtest.h>
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_balance/internal/balanceDefaults.hpp>

class ReverseOrderAttributes : public stk::unit_test_util::MeshFixture {};

TEST_F(ReverseOrderAttributes, balance_attributeOrderPreserved)
{
    const std::string inputFile = "reverseOrderAttr.exo";
    const std::string outputDir = "outputDir";
    stk::balance::run_stk_rebalance(outputDir, inputFile, stk::balance::SD_DEFAULTS, MPI_COMM_WORLD);

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::io::StkMeshIoBroker stkIo(get_comm());
    const std::string balancedFile = outputDir + "/" + inputFile;
    stk::io::fill_mesh_preexisting(stkIo, balancedFile, get_bulk());
    stk::mesh::Part *blockPart = get_meta().get_part("block_17");
    stk::mesh::FieldVector balancedAttrFields = stkIo.get_ordered_attribute_fields(blockPart);

    std::vector<std::string> expectedAttrFieldNames = {"j", "i", "h", "g", "f", "e", "d", "c", "b", "a" };
    ASSERT_EQ(expectedAttrFieldNames.size(), balancedAttrFields.size());
    for(size_t i=0; i<expectedAttrFieldNames.size(); i++)
        EXPECT_EQ(expectedAttrFieldNames[i], balancedAttrFields[i]->name());
}
