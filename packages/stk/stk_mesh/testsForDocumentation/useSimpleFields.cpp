#include <gtest/gtest.h>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

namespace {

//BEGIN
TEST(stkMeshHowTo, useSimpleFields)
{
    const unsigned spatialDimension = 3;
    stk::mesh::MetaData metaData(spatialDimension, stk::mesh::entity_rank_names());

    typedef stk::mesh::Field<double> ScalarField;
    typedef stk::mesh::Field<double, stk::mesh::Cartesian3d> VectorField;
    ScalarField& pressureField = metaData.declare_field<ScalarField>(stk::topology::ELEM_RANK, "pressure");
    VectorField& displacementsField = metaData.declare_field<VectorField>(stk::topology::NODE_RANK, "displacements");

    double initialPressureValue = 4.4;
    stk::mesh::put_field_on_entire_mesh_with_initial_value(pressureField, &initialPressureValue);
    stk::mesh::put_field_on_entire_mesh(displacementsField);

    stk::mesh::Part &tetPart = metaData.declare_part_with_topology("tetElementPart", stk::topology::TET_4);

    metaData.commit();
    stk::mesh::BulkData mesh(metaData, MPI_COMM_WORLD);
    mesh.modification_begin();
    stk::mesh::EntityId elem1Id = 1;
    stk::mesh::EntityId elem1Nodes[] = {1, 2, 3, 4};
    stk::mesh::Entity elem1=stk::mesh::declare_element(mesh, tetPart, elem1Id, elem1Nodes);
    stk::mesh::EntityId elem2Id = 2;
    stk::mesh::EntityId elem2Nodes[] = {2, 3, 4, 5};
    stk::mesh::Entity elem2=stk::mesh::declare_element(mesh, tetPart, elem2Id, elem2Nodes);
    mesh.modification_end();

    const stk::mesh::BucketVector& nodeBuckets = mesh.buckets(stk::topology::NODE_RANK);
    EXPECT_TRUE(!nodeBuckets.empty());
    for(size_t bucketIndex=0; bucketIndex<nodeBuckets.size(); bucketIndex++)
    {
        const stk::mesh::Bucket& bucket = *nodeBuckets[bucketIndex];
        double* displacementDataForBucket = stk::mesh::field_data(displacementsField, bucket);
        EXPECT_GT(bucket.size(), 0u);
        for(size_t nodeIndex=0; nodeIndex<bucket.size(); nodeIndex++)
        {
            unsigned numValuesPerNode = stk::mesh::field_scalars_per_entity(displacementsField, bucket);
            const unsigned cartesian_enum_size = stk::mesh::Cartesian3d::Size;
            EXPECT_EQ(cartesian_enum_size, numValuesPerNode);
            for(unsigned i=0; i<numValuesPerNode; i++)
            {
                EXPECT_EQ(0.0, displacementDataForBucket[nodeIndex*numValuesPerNode + i]);
                displacementDataForBucket[nodeIndex*numValuesPerNode + i] = 99.9;
            }
        }
    }

    double* pressureFieldDataForElem1 = stk::mesh::field_data(pressureField, elem1);
    EXPECT_EQ(initialPressureValue, *pressureFieldDataForElem1);

    double* pressureFieldDataForElem2 = stk::mesh::field_data(pressureField, elem2);
    EXPECT_EQ(initialPressureValue, *pressureFieldDataForElem2);
}
//END

}
