#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <fieldNameTestUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>

namespace
{
TEST(StkMeshIoBrokerHowTo, writeFieldWithAlternativeName)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    size_t index = stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();

    // describe nodal 'displacement' field
    const int numberOfStates = 1;
    const std::string primaryFieldName = "displacement";
    stk::mesh::Field<double> &nodalDisplacement = stkMeshMetaData.
      declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK,
					       primaryFieldName, numberOfStates);
    stk::mesh::put_field(nodalDisplacement, stkMeshMetaData.universal_part());

    const std::string outputFileName = "resultsOutput.exo";
    size_t resultsOutputIndex = stkIo.create_output_mesh(outputFileName,
							 stk::io::WRITE_RESULTS);

    // mark field for output
    std::string alternateFieldName("deformation");
    stkIo.add_field(resultsOutputIndex, nodalDisplacement, alternateFieldName);

    // field descriptions done, fill field data (done with mesh input)
    stkIo.populate_bulk_data();

    EXPECT_TRUE( fieldWithNameChangedIsOutput(stkIo, communicator, resultsOutputIndex,
					      alternateFieldName));
    unlink(outputFileName.c_str());
}
}
