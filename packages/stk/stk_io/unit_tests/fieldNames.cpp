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

void generateMetaData(stk::io::MeshData &stkIo)
{
    const std::string exodusFileName = "generated:1x1x1";
    stkIo.open_mesh_database(exodusFileName);
    stkIo.create_input_mesh();
}

void createNamedFieldOnMesh(stk::mesh::MetaData &stkMeshMetaData, const std::string &internalClientFieldName)
{
    const int numberOfStates = 1;
    stk::mesh::Field<double> &field0 = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(internalClientFieldName, numberOfStates);
    stk::mesh::put_field(field0, stk::mesh::Entity::NODE, stkMeshMetaData.universal_part());
}

void markFieldForResultsOutputWithNewName(stk::mesh::FieldBase *field0, const std::string &requestedFieldNameForResultsOutput)
{
    stk::io::set_field_role(*field0, Ioss::Field::TRANSIENT);
    stk::io::set_results_field_name(*field0, requestedFieldNameForResultsOutput);
}

void testFieldNamedCorrectly(Ioss::Region *ioRegion, MPI_Comm communicator, const std::string &goldFieldName)
{
    Ioss::NodeBlock *nodeBlockAssociatedWithField0 = ioRegion->get_node_blocks()[0];
    Ioss::NameList fieldNames;
    nodeBlockAssociatedWithField0->field_describe(Ioss::Field::TRANSIENT, &fieldNames);

    ASSERT_EQ(1u, fieldNames.size());
    EXPECT_STREQ(goldFieldName.c_str(), fieldNames[0].c_str());
}

STKUNIT_UNIT_TEST(StkIoTest, FieldName)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string internalClientFieldName = "Field0";
    stk::io::MeshData stkIo(communicator);
    generateMetaData(stkIo);

    stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
    createNamedFieldOnMesh(stkMeshMetaData, internalClientFieldName);

    stk::mesh::FieldBase *field0 = stkMeshMetaData.get_field(internalClientFieldName);
    std::string requestedFieldNameForResultsOutput("jeSSe");
    markFieldForResultsOutputWithNewName(field0, requestedFieldNameForResultsOutput);

    stkIo.populate_bulk_data();

    const std::string outputFileName = "ourSillyOutput.exo";
    stkIo.create_output_mesh(outputFileName);
    stkIo.define_output_fields();

    Ioss::Region *ioRegion = stkIo.output_io_region().get();
    testFieldNamedCorrectly(ioRegion, communicator, requestedFieldNameForResultsOutput);

    unlink(outputFileName.c_str());
}

STKUNIT_UNIT_TEST(StkIoTest, FieldNameRenameTwice)
{
    std::string restartFilename = "output.restart";
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string internalClientFieldName = "Field0";
    stk::io::MeshData stkIo(communicator);
    generateMetaData(stkIo);

    stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
    createNamedFieldOnMesh(stkMeshMetaData, internalClientFieldName);

    stk::mesh::FieldBase *field0 = stkMeshMetaData.get_field(internalClientFieldName);
    std::string requestedFieldNameForResultsOutput("NotjeSSe");
    markFieldForResultsOutputWithNewName(field0, requestedFieldNameForResultsOutput);

    requestedFieldNameForResultsOutput = "jeSSe";
    markFieldForResultsOutputWithNewName(field0, requestedFieldNameForResultsOutput);

    stkIo.populate_bulk_data();

    const std::string outputFileName = "ourSillyOwlput.exo";
    stkIo.create_output_mesh(outputFileName);
    stkIo.define_output_fields();

    Ioss::Region *ioRegion = stkIo.output_io_region().get();
    testFieldNamedCorrectly(ioRegion, communicator, requestedFieldNameForResultsOutput);

    unlink(outputFileName.c_str());
}

STKUNIT_UNIT_TEST(StkIoTest, FieldNameWithRestart)
{
    std::string restartFilename = "output.restart";
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string internalClientFieldName = "Field0";
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);
        
        stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
        createNamedFieldOnMesh(stkMeshMetaData, internalClientFieldName);
        
        stk::mesh::FieldBase *field0 = stkMeshMetaData.get_field(internalClientFieldName);

        stkIo.add_restart_field(*field0);

        stkIo.populate_bulk_data();

        stkIo.create_restart_output(restartFilename);
        stkIo.define_restart_fields();

        double time = 0.0;
        stkIo.process_restart_output(time);
    }

    Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", restartFilename, Ioss::READ_RESTART, communicator);
    Ioss::Region iossRegion(iossDb);
    std::string goldFieldName = internalClientFieldName;
    sierra::make_lower(goldFieldName);
    testFieldNamedCorrectly(&iossRegion, communicator, goldFieldName);

    unlink(restartFilename.c_str());
}

STKUNIT_UNIT_TEST(StkIoTest, FieldNameWithResultsAndRestart)
{
    const std::string restartFilename = "output2.restart";
    const std::string outputFileName = "resultsOutput.exo";
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string internalClientFieldName = "Field0";
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);

        stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
        createNamedFieldOnMesh(stkMeshMetaData, internalClientFieldName);

        stk::mesh::FieldBase *field0 = stkMeshMetaData.get_field(internalClientFieldName);
        std::string requestedFieldNameForResultsOutput("jeSSe");
        markFieldForResultsOutputWithNewName(field0, requestedFieldNameForResultsOutput);

        stkIo.add_restart_field(*field0);

        stkIo.populate_bulk_data();

        stkIo.create_output_mesh(outputFileName);
        stkIo.define_output_fields();

        stkIo.create_restart_output(restartFilename);
        stkIo.define_restart_fields();

        double time = 0.0;
        stkIo.process_output_request(time);
        stkIo.process_restart_output(time);

        Ioss::Region *ioRegion = stkIo.output_io_region().get();
        testFieldNamedCorrectly(ioRegion, communicator, requestedFieldNameForResultsOutput);
    }

    Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", restartFilename, Ioss::READ_RESTART, communicator);
    Ioss::Region iossRegion(iossDb);
    std::string goldFieldName = internalClientFieldName;
    sierra::make_lower(goldFieldName);
    testFieldNamedCorrectly(&iossRegion, communicator, goldFieldName);

    unlink(outputFileName.c_str());
    unlink(restartFilename.c_str());
}

}
