#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/ElemElemGraph.hpp>
#include <unit_tests/BulkDataTester.hpp>
#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_io/IossBridge.hpp>
#include <unistd.h>
#include "FaceCreatorFixture.hpp"

namespace
{

class FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester : public FaceCreatorFixture
{
protected:

    FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester() : FaceCreatorFixture(2) {}

    virtual void create_faces(stk::mesh::Entity element, stk::mesh::EntityVector& nodes_of_face)
    {
        bulkData->modification_begin();
        create_face_per_proc(element, nodes_of_face);
        test_that_num_sides_is_expected_value(2);
        bulkData->modification_end();
    }

    void setup_mesh(stk::mesh::BulkData::AutomaticAuraOption aura_option)
    {
        bulkData = new stk::mesh::unit_test::BulkDataElemGraphFaceSharingTester(metaData, get_comm(), aura_option);
        unsigned numX = 2, numY = 1;
        convert_quad_fixture_to_my_bulk_data_flavor(numX, numY);
    }

    void convert_quad_fixture_to_my_bulk_data_flavor(unsigned numX, unsigned numY)
    {
        stk::mesh::fixtures::QuadFixture fixture(get_comm(), numX, numY, false);

        stk::mesh::Field<double, stk::mesh::Cartesian2d> &coordField = fixture.m_meta.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian2d>>(stk::topology::NODE_RANK, "model_coordinates");
        stk::mesh::put_field(coordField, fixture.m_meta.universal_part(), fixture.m_meta.spatial_dimension());
        stk::mesh::Part& block_1 = fixture.m_meta.declare_part_with_topology("block_1", stk::topology::QUADRILATERAL_4_2D);
        stk::io::put_io_part_attribute(block_1);

        fixture.m_meta.commit();
        fixture.generate_mesh();

        std::vector<double> x;
        std::vector<double> y;

        for(unsigned j=0;j<=numY;++j)
        {
            for(unsigned i=0;i<=numX;i++)
            {
                x.push_back(i); // 0 1 2, 0 1 2, 0 1 2, ...
                y.push_back(j); // 0 0 0, 1 1 1
            }
        }

        stk::mesh::EntityVector nodes;
        stk::mesh::get_selected_entities(fixture.m_meta.universal_part(), fixture.m_bulk_data.buckets(stk::topology::NODE_RANK), nodes);
        for(stk::mesh::Entity node : nodes )
        {
            double* coords = stk::mesh::field_data(coordField, node);
            unsigned id = fixture.m_bulk_data.identifier(node);
            coords[0] = x[id-1];
            coords[1] = y[id-1];
        }

        fixture.m_bulk_data.modification_begin();
        stk::mesh::EntityVector elems;
        stk::mesh::get_selected_entities(fixture.m_meta.locally_owned_part(), fixture.m_bulk_data.buckets(stk::topology::ELEM_RANK), elems);
        for(stk::mesh::Entity element : elems)
        {
            fixture.m_bulk_data.change_entity_parts(element, {&block_1});
        }
        fixture.m_bulk_data.modification_end();

        std::ostringstream os;
        const std::string file_temp("testadfasdasdfas.exo");
        stk::unit_test_util::write_mesh_using_stk_io(file_temp, fixture.m_bulk_data, get_comm());
        stk::unit_test_util::fill_mesh_using_stk_io(file_temp, *bulkData, get_comm());

        ThrowRequireMsg(fixture.m_bulk_data.parallel_size()<10, "Testing assumption violated.");
        os << file_temp << "." << fixture.m_bulk_data.parallel_size() << "." << fixture.m_bulk_data.parallel_rank();
        unlink(os.str().c_str());
    }

    virtual stk::mesh::EntityVector get_nodes_of_face_for_this_proc()
    {
        std::vector<unsigned> face_node_ids = { 2, 5 };
        return get_nodes_for_proc(face_node_ids);
    }

    virtual unsigned get_permuted_index(unsigned i)
    {
        std::vector<std::vector<unsigned> > index_for_proc = {
                {0, 1},
                {1, 0}
        };
        return index_for_proc[get_bulk().parallel_rank()][i];
    }

    void test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face()
    {
        only_proc_0_makes_a_face();
        test_that_num_sides_is_expected_value(1);
        test_that_each_proc_has_num_sides_with_expected_value(1);
    }

    void only_proc_0_makes_a_face()
    {
        unsigned id = get_bulk().parallel_rank()+1;
        stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, id);
        stk::mesh::EntityVector nodes_of_face = get_nodes_of_face_for_this_proc();
        create_faces(elem, nodes_of_face);
    }

    void test_that_each_proc_has_num_sides_with_expected_value(unsigned expected_num_sides)
    {
        unsigned num_local_sides = stk::mesh::count_selected_entities(get_bulk().mesh_meta_data().globally_shared_part(), get_bulk().buckets(get_bulk().mesh_meta_data().side_rank()));
        EXPECT_EQ(expected_num_sides, num_local_sides);
    }
};


TEST_F(FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateTwoFacesWithAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_mesh(stk::mesh::BulkData::AUTO_AURA);
        test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
    }
}

TEST_F(FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateTwoFacesWithoutAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
    }
}

TEST_F(FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateOneFaceWithAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_mesh(stk::mesh::BulkData::AUTO_AURA);
        test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face();
    }
}

TEST_F(FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateOneFaceWithoutAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face();
    }
}

}
