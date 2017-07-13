#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/balance.hpp>
#include <gtest/gtest.h>
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_unit_test_utils/FaceCreationTestUtils.hpp>

namespace {

class CoincidentElems : public stk::unit_test_util::MeshFixture
{
protected:
    CoincidentElems()
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    }
    void expect_coincidents_on_same_proc(stk::mesh::EntityId elem1Id, stk::mesh::EntityId elem2Id)
    {
        stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, elem1Id);
        if(get_bulk().is_valid(elem1))
        {
            stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, elem2Id);
            EXPECT_TRUE(get_bulk().is_valid(elem2));
        }
    }
    void expect_proper_face_connection(const std::string& filename, size_t numSides, const std::vector<SideTestUtil::Side>& sideSet)
    {
        SideTestUtil::expect_global_num_sides_in_part(get_bulk(), numSides, get_meta().universal_part());
        SideTestUtil::expect_all_sides_exist_for_elem_side(get_bulk(), filename, sideSet);
    }
    void test_decomp_and_balance(const std::string& fileName,
                                 stk::mesh::EntityId elem1Id,
                                 stk::mesh::EntityId elem2Id,
                                 const size_t numGlobalSides,
                                 const std::vector<SideTestUtil::Side>& expectedElemSides)
    {
        EXPECT_NO_THROW(stk::balance::initial_decomp_and_balance(get_bulk(), graphOptions, fileName, "output_" + fileName));
        expect_coincidents_on_same_proc(elem1Id, elem2Id);
        expect_proper_face_connection(fileName, numGlobalSides, expectedElemSides);
    }
    stk::balance::GraphCreationSettings graphOptions;
};

TEST_F(CoincidentElems, balance_coincidentsNotSplit)
{
    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
         1,2,HEX_8,5,6,7,8,9,10,11,12\n\
         1,3,HEX_8,5,6,7,8,9,10,11,12";
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());

    expect_coincidents_on_same_proc(2, 3);
    stk::balance::balanceStkMesh(graphOptions, get_bulk());
    expect_coincidents_on_same_proc(2, 3);
}

TEST_F(CoincidentElems, linearDecompJL_noThrow)
{
    test_decomp_and_balance("JL.e", 1, 2, 1, std::vector<SideTestUtil::Side>{{1,5}, {2,5}});
}

TEST_F(CoincidentElems, linearDecompALJ_noThrow)
{
    test_decomp_and_balance("ALJ.e", 2, 3, 1, std::vector<SideTestUtil::Side>{{1,5}, {2,4}, {3,4}});
}

TEST_F(CoincidentElems, linearDecompARefLA_noThrow)
{
    test_decomp_and_balance("ARefLA.e", 3, 4, 2, std::vector<SideTestUtil::Side>{{1,5}, {2,4}, {3,0}, {3,1}, {4,0}, {4,1}});
}

}
