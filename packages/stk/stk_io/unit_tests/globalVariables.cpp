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
    stkIo.populate_bulk_data();
}

void testGlobalVarOnFile(const std::string &outputFileName, const std::vector<std::string> &goldGlobalVarName, const std::vector<double> goldGlobalVarValue, MPI_Comm comm)
{
    Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", outputFileName, Ioss::READ_RESTART, comm);
    Ioss::Region inputRegion(iossDb);
    Ioss::NameList fieldNames;
    //Get all the fields on the ioss region, which are the global variables
    inputRegion.field_describe(Ioss::Field::TRANSIENT, &fieldNames);

    ASSERT_EQ(goldGlobalVarName.size(), fieldNames.size());

    const int stepNumber = 1;
    inputRegion.begin_state(stepNumber);

    const double tolerance = 1e-16;
    for(size_t i=0; i<goldGlobalVarName.size(); i++)
    {
        EXPECT_STRCASEEQ(goldGlobalVarName[i].c_str(), fieldNames[i].c_str());
        double globalValue = -1.0;
        inputRegion.get_field_data(fieldNames[i], &globalValue, sizeof(double));
        EXPECT_NEAR(goldGlobalVarValue[i], globalValue, tolerance);
    }
}

STKUNIT_UNIT_TEST(StkIoTest, OneGlobalDouble)
{
    const std::string outputFileName = "ourSillyOutput.exo";
    const std::string globalVarName = "testGlobal";
    const double globalVarValue = 13.0;
    MPI_Comm communicator = MPI_COMM_WORLD;
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);

        stkIo.create_output_mesh(outputFileName);

        stkIo.add_results_global(globalVarName, Ioss::Field::REAL);

        const double time = 1.0;
        stkIo.begin_results_output_at_time(time);

        stkIo.write_results_global(globalVarName, globalVarValue);

        stkIo.end_current_results_output();
    }

    std::vector<std::string> globalVarNames(1, globalVarName);
    std::vector<double> globalVarValues(1,globalVarValue);
    testGlobalVarOnFile(outputFileName, globalVarNames, globalVarValues, communicator);
    unlink(outputFileName.c_str());
}

STKUNIT_UNIT_TEST(StkIoTest, TwoGlobalDoubles)
{
    const std::string outputFileName = "ourSillyOutput.exo";
    MPI_Comm communicator = MPI_COMM_WORLD;
    std::vector<std::string> globalVarNames;
    globalVarNames.push_back("testGlobal");
    globalVarNames.push_back("testGlobal2");
    std::vector<double> globalVarValues;
    globalVarValues.push_back(13.0);
    globalVarValues.push_back(14.0);
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);

        stkIo.create_output_mesh(outputFileName);

        stkIo.add_results_global(globalVarNames[0], Ioss::Field::REAL);
        stkIo.add_results_global(globalVarNames[1], Ioss::Field::REAL);

        const double time = 1.0;
        stkIo.begin_results_output_at_time(time);

        stkIo.write_results_global(globalVarNames[0], globalVarValues[0]);
        stkIo.write_results_global(globalVarNames[1], globalVarValues[1]);

        stkIo.end_current_results_output();
    }

    testGlobalVarOnFile(outputFileName, globalVarNames, globalVarValues, communicator);
    unlink(outputFileName.c_str());
}

}
