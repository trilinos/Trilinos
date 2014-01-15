#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_topology/topology.hpp>

namespace
{
TEST(StkMeshHowTo, betterUnderstandSelectorConstruction)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x1"; // syntax creates a 1x1x1 cube
    stkMeshIoBroker.open_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::mesh::Selector nothingSelector_byDefaultConstruction;
    size_t expectingZeroBuckets = 0;
    EXPECT_EQ(expectingZeroBuckets, stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, nothingSelector_byDefaultConstruction).size());

    stk::mesh::Selector allSelector(!nothingSelector_byDefaultConstruction);
    size_t numberOfAllNodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK).size();
    EXPECT_EQ(numberOfAllNodeBuckets, stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, allSelector).size());
}

TEST(StkMeshHowTo, makeSureYouAreNotIntersectingNothingSelector)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    // syntax creates faces for surface on the positive: 'x-side', 'y-side', and 'z-side'
    // of a 1x1x1 cube, these parts are given the names: 'surface_1', 'surface_2', and 'surface_3'
    // automagically when it is created [create_input_mesh()]
    const std::string generatedMeshSpecification = "generated:1x1x1|sideset:XYZ";
    stkMeshIoBroker.open_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::Part *surface1Part = stkMeshMetaData.get_part("surface_1");
    stk::mesh::Part *surface2Part = stkMeshMetaData.get_part("surface_2");
    stk::mesh::Part *surface3Part = stkMeshMetaData.get_part("surface_3");
    stk::mesh::PartVector allSurfaces;
    allSurfaces.push_back(surface1Part);
    allSurfaces.push_back(surface2Part);
    allSurfaces.push_back(surface3Part);

    stk::mesh::Selector boundaryNodesSelectorWithNothingBeingIntersected;

    std::ostringstream readableSelectorDescription;
    readableSelectorDescription << boundaryNodesSelectorWithNothingBeingIntersected;
    EXPECT_STREQ("NOTHING", readableSelectorDescription.str().c_str());

    for (size_t surfaceIndex = 0; surfaceIndex < allSurfaces.size(); ++surfaceIndex)
    {
        stk::mesh::Part &surfacePart = *(allSurfaces[surfaceIndex]);
        stk::mesh::Selector surfaceSelector(surfacePart);
        boundaryNodesSelectorWithNothingBeingIntersected &= surfacePart;
    }
    size_t expectedNumberOfBucketsWhenIntersectingNothing = 0;
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();
    EXPECT_EQ(expectedNumberOfBucketsWhenIntersectingNothing, stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, boundaryNodesSelectorWithNothingBeingIntersected).size());

    stk::mesh::Selector preferredBoundaryNodesSelector = stk::mesh::selectIntersection(allSurfaces);
    size_t expectedNumberOfNodeBucketsWhenIntersectingAllSurfaces = 1;
    EXPECT_EQ(expectedNumberOfNodeBucketsWhenIntersectingAllSurfaces, stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, preferredBoundaryNodesSelector).size());
}
}
