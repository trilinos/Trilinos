#ifndef stkMeshTestUtilsHpp
#define stkMeshTestUtilsHpp

#include <gtest/gtest.h>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_topology/topology.hpp>

void testTemperatureFieldSetCorrectly(const stk::mesh::Field<double> &temperatureField, double prescribedTemperatureValue, const std::set<stk::mesh::EntityId> &boundaryNodeIds)
{
    stk::mesh::BulkData &stkMeshBulkData = temperatureField.get_mesh();
    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);
    for(size_t i=0; i<nodes.size(); ++i)
    {
        double *temperature = stk::mesh::field_data(temperatureField, nodes[i]);
        if(boundaryNodeIds.find(stkMeshBulkData.identifier(nodes[i])) != boundaryNodeIds.end())
        {
            EXPECT_EQ(prescribedTemperatureValue, *temperature);
        }
        else
        {
            EXPECT_EQ(0.0, *temperature);
        }
    }
}

#endif
