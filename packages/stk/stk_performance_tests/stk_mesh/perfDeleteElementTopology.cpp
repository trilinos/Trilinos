#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <string>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>

namespace
{

class DestroyElementTopologyPerformanceTest : public stk::unit_test_util::MeshFixture
{
protected:
    DestroyElementTopologyPerformanceTest()
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        create_num_fields_of_rank(10, stk::topology::NODE_RANK);
        create_num_fields_of_rank(10, stk::topology::ELEM_RANK);
        stk::io::fill_mesh("generated:50x50x64", get_bulk());
    }

    void expect_no_entities()
    {
        EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::ELEM_RANK)));
        EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::NODE_RANK)));
    }
private:
    void create_num_fields_of_rank(const unsigned numFields, stk::mesh::EntityRank rank)
    {
        const std::string name = "field_" + std::to_string(rank) + "_";
        for(unsigned i = 0; i < numFields; i++)
            create_vector_field_on_universal_part(rank, name + std::to_string(i));
    }
    void create_vector_field_on_universal_part(stk::mesh::EntityRank rank, const std::string& name)
    {
        auto& field = get_meta().declare_field<double>(rank, name, 3);
        stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
    }
};


class DestroyElementTopologyPerformance : public stk::unit_test_util::PerformanceTester
{
public:
    DestroyElementTopologyPerformance(stk::mesh::BulkData &bulk)
      : stk::unit_test_util::PerformanceTester(bulk.parallel()),
        bulkData(bulk)
    { }
protected:
    virtual void run_algorithm_to_time()
    {
        bulkData.destroy_elements_of_topology(stk::topology::HEX_8);
    }
    virtual size_t get_value_to_output_as_iteration_count() { return 1; }
private:
    stk::mesh::BulkData & bulkData;
};
TEST_F(DestroyElementTopologyPerformanceTest, DestroyElementTopology)
{
    DestroyElementTopologyPerformance(get_bulk()).run_performance_test();
    expect_no_entities();
}


class DestroyAllElementsIndividuallyPerformance : public stk::unit_test_util::PerformanceTester
{
public:
    DestroyAllElementsIndividuallyPerformance(stk::mesh::BulkData &bulk)
      : stk::unit_test_util::PerformanceTester(bulk.parallel()),
        bulkData(bulk)
    { }
protected:
    virtual void run_algorithm_to_time()
    {
        bulkData.modification_begin();
        for(stk::mesh::EntityRank rank : {stk::topology::ELEM_RANK, stk::topology::FACE_RANK ,stk::topology::NODE_RANK})
            destroy_entities_of_rank(rank);
        bulkData.modification_end();
    }
    virtual size_t get_value_to_output_as_iteration_count() { return 1; }
private:
    void destroy_entities_of_rank(stk::mesh::EntityRank rank)
    {
        stk::mesh::EntityVector entities;
        stk::mesh::get_entities(bulkData, rank, entities);
        for(stk::mesh::Entity elem : entities)
            bulkData.destroy_entity(elem);
    }
    stk::mesh::BulkData & bulkData;
};
TEST_F(DestroyElementTopologyPerformanceTest, DestroyElementsIndividually)
{
    DestroyAllElementsIndividuallyPerformance(get_bulk()).run_performance_test();
    expect_no_entities();
}

}
