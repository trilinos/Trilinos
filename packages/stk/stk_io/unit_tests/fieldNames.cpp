#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <string>
#include <mpi.h>
#include <stk_io/MeshReadWriteUtils.hpp>
#include <stk_io/IossBridge.hpp>
#include <Ioss_SubSystem.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>

namespace {

STKUNIT_UNIT_TEST(StkIoTest, FieldName)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string exodusFileName = "generated:1x1x1";

    stk::io::MeshData stkIo(communicator);
    stkIo.open_mesh_database(exodusFileName);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
    const int numberOfStates = 1;
    const std::string internalClientFieldName = "Field0";
    stk::mesh::Field<double> &field0 = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(internalClientFieldName, numberOfStates);
    stk::mesh::put_field(field0, stk::mesh::Entity::NODE, stkMeshMetaData.universal_part());
    stk::io::set_field_role(field0, Ioss::Field::TRANSIENT);

    const std::string requestedFieldNameForOutput("jeSSe");
    stk::io::set_results_field_name(field0, requestedFieldNameForOutput);

    stkIo.populate_bulk_data();

    stk::mesh::FieldBase *thisWontBurnUs = stkMeshMetaData.get_field(internalClientFieldName);
    EXPECT_TRUE(NULL != thisWontBurnUs);

    const std::string outputFileName = "ourSillyOwlput.exo";
    stkIo.create_output_mesh(outputFileName);
    stkIo.define_output_fields();

    Ioss::Region *ioRegion = stkIo.output_io_region().get();
    Ioss::NodeBlock *nodeBlockAssociatedWithField0 = ioRegion->get_node_blocks()[0];
    Ioss::NameList fieldNames;
    nodeBlockAssociatedWithField0->field_describe(Ioss::Field::TRANSIENT, &fieldNames);

    ASSERT_EQ(1u, fieldNames.size());
    EXPECT_STREQ(requestedFieldNameForOutput.c_str(), fieldNames[0].c_str());

    unlink(outputFileName.c_str());
}

STKUNIT_UNIT_TEST(StkIoTest, FieldNameWithRestart)
{
    std::string restartFilename = "output.restart";
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string internalClientFieldName = "Field0";
    {
        const std::string exodusFileName = "generated:1x1x1";

        stk::io::MeshData stkIo(communicator);
        stkIo.open_mesh_database(exodusFileName);
        stkIo.create_input_mesh();

        stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
        const int numberOfStates = 1;
        stk::mesh::Field<double> &field0 = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(internalClientFieldName, numberOfStates);
        stk::mesh::put_field(field0, stk::mesh::Entity::NODE, stkMeshMetaData.universal_part());
        stk::io::set_field_role(field0, Ioss::Field::TRANSIENT);

        std::string requestedFieldNameForOutput("notjeSSe");
        stk::io::set_results_field_name(field0, requestedFieldNameForOutput);
        requestedFieldNameForOutput = "jeSSe";
        stk::io::set_results_field_name(field0, requestedFieldNameForOutput);

        stkIo.add_restart_field(field0);

        stkIo.populate_bulk_data();

        stk::mesh::FieldBase *thisWontBurnUs = stkMeshMetaData.get_field(internalClientFieldName);
        EXPECT_TRUE(NULL != thisWontBurnUs);

        stkIo.create_restart_output(restartFilename);
        stkIo.define_restart_fields();

        double time = 0.0;
        stkIo.process_restart_output(time);
    }

    Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", restartFilename, Ioss::READ_RESTART, communicator);
    Ioss::Region iossRegion(iossDb);
    Ioss::NodeBlock *nodeBlockAssociatedWithField0 = iossRegion.get_node_blocks()[0];
    Ioss::NameList fieldNames;
    nodeBlockAssociatedWithField0->field_describe(Ioss::Field::TRANSIENT, &fieldNames);

    ASSERT_EQ(1u, fieldNames.size());
    EXPECT_STRCASEEQ(internalClientFieldName.c_str(), fieldNames[0].c_str()); //lowercased on restart read

    unlink(restartFilename.c_str());
}

// another test, output to disk and read back in
//    const std::string outputFileName = "ourSillyOwlput.exo";
//    stkIo.create_output_mesh(outputFileName);
//    stkIo.define_output_fields();
//    double time = 0.0;
//    stkIo.process_output_request(time);
//
//
//    stk::io::MeshData stkIo(communicator);
//    stkIo.open_mesh_database(outputFileName);
//    stkIo.create_input_mesh();

}
