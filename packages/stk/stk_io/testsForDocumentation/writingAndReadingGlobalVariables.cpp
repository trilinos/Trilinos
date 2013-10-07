#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <mpi.h>
#include <stk_io/MeshReadWriteUtils.hpp>
#include <Ioss_SubSystem.h>
#include <fieldNameTestUtils.hpp>

namespace
{
TEST(StkIoTestForDocumentation, WriteAndReadGlobalVariables)
{
    const std::string restartFileName = "OneGlobalDouble.restart";
    const std::string timeStepVarName = "timeStep";
    const double timeStepSize = 1e-6;
    const double currentTime = 1.0;
    MPI_Comm communicator = MPI_COMM_WORLD;

    //Write restart file with time step size as a global variable
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);
        stkIo.populate_bulk_data();

        stkIo.create_restart_output(restartFileName);
        stkIo.add_restart_global(timeStepVarName, Ioss::Field::REAL);
        stkIo.begin_restart_output_at_time(currentTime);
        stkIo.write_restart_global(timeStepVarName, timeStepSize);
        stkIo.end_current_restart_output();
    }

    //Read restart file with time step size as a global variable
    {
        stk::io::MeshData stkIo(communicator);
        stkIo.open_mesh_database(restartFileName);
        stkIo.create_input_mesh();
        stkIo.populate_bulk_data();
        stkIo.process_restart_input(currentTime);
        std::vector<std::string> globalNamesOnFile;
        stkIo.get_global_variable_names(globalNamesOnFile);

        ASSERT_EQ(1u, globalNamesOnFile.size());
        EXPECT_STRCASEEQ(timeStepVarName.c_str(),
                         globalNamesOnFile[0].c_str());
        double timeStepSizeReadFromFile =
                stkIo.get_global(globalNamesOnFile[0]);
        const double tolerance = 1e-16;
        EXPECT_NEAR(timeStepSize, timeStepSizeReadFromFile, tolerance);
    }
    unlink(restartFileName.c_str());
}
}
