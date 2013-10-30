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

void testFieldNamedCorrectly(Ioss::Region &ioRegion, MPI_Comm communicator, const std::string &goldFieldName)
{
    Ioss::NodeBlock *nodeBlockAssociatedWithField0 = ioRegion.get_node_blocks()[0];
    Ioss::NameList fieldNames;
    nodeBlockAssociatedWithField0->field_describe(Ioss::Field::TRANSIENT, &fieldNames);

    ASSERT_EQ(1u, fieldNames.size());
    EXPECT_STRCASEEQ(goldFieldName.c_str(), fieldNames[0].c_str());
}

STKUNIT_UNIT_TEST(FieldNamesTest, FieldNameRenameTwice)
{
    const std::string outputFilename = "ourSillyOwlput.exo";
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string internalClientFieldName = "Field0";
    std::string requestedFieldNameForResultsOutput("NotjeSSe");
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);

        stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
        createNamedFieldOnMesh(stkMeshMetaData, internalClientFieldName);
        stkIo.populate_bulk_data();

        size_t results_output_index = stkIo.create_output_mesh(outputFilename);

        stk::mesh::FieldBase *field0 = stkMeshMetaData.get_field(internalClientFieldName);
        stkIo.add_results_field(results_output_index, *field0, requestedFieldNameForResultsOutput);

        requestedFieldNameForResultsOutput = "jeSSe";
        stkIo.add_results_field(results_output_index, *field0, requestedFieldNameForResultsOutput);

        //stkIo.define_output_fields();
        stkIo.process_output_request(0.0, results_output_index);
    }

    Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", outputFilename, Ioss::READ_MODEL, communicator);
    Ioss::Region ioRegion(iossDb);
    testFieldNamedCorrectly(ioRegion, communicator, requestedFieldNameForResultsOutput);

    unlink(outputFilename.c_str());
}

STKUNIT_UNIT_TEST(FieldNamesTest, FieldNameWithRestart)
{
    std::string restartFilename = "output.restart";
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string internalClientFieldName = "Field0";
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);
        
        stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
        createNamedFieldOnMesh(stkMeshMetaData, internalClientFieldName);
        stkIo.populate_bulk_data();
        
        stk::mesh::FieldBase *field0 = stkMeshMetaData.get_field(internalClientFieldName);

        stkIo.add_restart_field(*field0);


        stkIo.create_restart_output(restartFilename);
        stkIo.define_restart_fields();

        double time = 0.0;
        stkIo.begin_restart_output_at_time(time);
        stkIo.process_restart_output();
        stkIo.end_current_restart_output();
    }

    Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", restartFilename, Ioss::READ_RESTART, communicator);
    Ioss::Region iossRegion(iossDb);
    std::string goldFieldName = internalClientFieldName;
    sierra::make_lower(goldFieldName);
    testFieldNamedCorrectly(iossRegion, communicator, goldFieldName);

    unlink(restartFilename.c_str());
}

STKUNIT_UNIT_TEST(FieldNamesTest, FieldNameWithResultsAndRestart)
{
    const std::string restartFilename = "FieldNameWithResultsAndRestart.restart";
    const std::string outputFileName = "FieldNameWithResultsAndRestart.exo";
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string internalClientFieldName = "Field0";
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);

        stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
        createNamedFieldOnMesh(stkMeshMetaData, internalClientFieldName);
        stkIo.populate_bulk_data();

        size_t results_output_index = stkIo.create_output_mesh(outputFileName);
        stk::mesh::FieldBase *field0 = stkMeshMetaData.get_field(internalClientFieldName);
        std::string requestedFieldNameForResultsOutput("jeSSe");
        stkIo.add_results_field(results_output_index, *field0, requestedFieldNameForResultsOutput);

        stkIo.add_restart_field(*field0);

        stkIo.create_restart_output(restartFilename);
        stkIo.define_restart_fields();

        double time = 0.0;
        stkIo.begin_results_output_at_time(time, results_output_index);
        stkIo.process_output_request(results_output_index);
        stkIo.end_current_results_output(results_output_index);

        stkIo.begin_restart_output_at_time(time);
        stkIo.process_restart_output();
        stkIo.end_current_restart_output();
    }
    Ioss::DatabaseIO *iossResultDb = Ioss::IOFactory::create("exodus", restartFilename, Ioss::READ_MODEL, communicator);
    Ioss::Region resultRegion(iossResultDb);
    std::string goldFieldName = internalClientFieldName;
    sierra::make_lower(goldFieldName);
    testFieldNamedCorrectly(resultRegion, communicator, goldFieldName);

    Ioss::DatabaseIO *iossRestartDb = Ioss::IOFactory::create("exodus", restartFilename, Ioss::READ_RESTART, communicator);
    Ioss::Region restartRegion(iossRestartDb);
    testFieldNamedCorrectly(restartRegion, communicator, goldFieldName);

    unlink(outputFileName.c_str());
    unlink(restartFilename.c_str());
}

}
