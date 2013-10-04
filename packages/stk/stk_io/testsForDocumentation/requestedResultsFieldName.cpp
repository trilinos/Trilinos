#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/MeshReadWriteUtils.hpp>
#include <stk_io/IossBridge.hpp>
#include <fieldNameTestUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>

namespace
{
TEST(StkIoTestForDocumentation, OutputtingFieldWithAlternativeName)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::MeshData stkIo(communicator);
    stk::mesh::MetaData &stkMeshMetaData = generateMetaData(stkIo);

    // describe nodal 'displacement' field
    const int numberOfStates = 1;
    const std::string primaryFieldName = "displacement";
    stk::mesh::Field<double> &nodalDisplacement = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(primaryFieldName, numberOfStates);
    stk::mesh::put_field(nodalDisplacement, stk::mesh::Entity::NODE, stkMeshMetaData.universal_part());

    // mark field for output
    std::string alternateFieldName("deformation");
    stkIo.add_results_field(nodalDisplacement, alternateFieldName);

    // field descriptions done, fill field data (done with mesh input)
    stkIo.populate_bulk_data();

    EXPECT_TRUE( fieldWithNameChangedIsOutput(stkIo, communicator, alternateFieldName) );
}
}
