#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetEntities.hpp>

inline stk::mesh::FieldBase* declareTriStateNodalField(stk::mesh::MetaData &stkMeshMetaData, const std::string &fieldName)
{
    const int numberOfStates = 3;
    stk::mesh::Field<double> &triStateField = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(fieldName, numberOfStates);
    stk::mesh::put_field(triStateField, stk::mesh::Entity::NODE, stkMeshMetaData.universal_part());
    return &triStateField;
}

inline void putDataOnTestField(stk::mesh::BulkData &stkMeshBulkData, const double value, stk::mesh::FieldBase &field)
{
    std::vector<stk::mesh::Entity> nodes;
    stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);
    for(size_t i=0; i<nodes.size(); i++)
    {
        double *fieldDataForNode = reinterpret_cast<double*>(stkMeshBulkData.field_data(field, nodes[i]));
        *fieldDataForNode = value;
    }
}

inline void testDataOnField(stk::mesh::BulkData &stkMeshBulkData, const double goldValue, stk::mesh::FieldBase &field)
{
    std::vector<stk::mesh::Entity> nodes;
    stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);
    for(size_t i=0; i<nodes.size(); i++)
    {
        double *fieldDataForNode = reinterpret_cast<double*>(stkMeshBulkData.field_data(field, nodes[i]));
        EXPECT_DOUBLE_EQ(goldValue, *fieldDataForNode);
    }
}

