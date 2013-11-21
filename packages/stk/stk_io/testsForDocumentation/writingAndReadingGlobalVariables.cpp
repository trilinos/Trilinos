#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Ioss_SubSystem.h>
#include <fieldNameTestUtils.hpp>

namespace
{
TEST(StkMeshIoBrokerHowTo, writeAndReadGlobalVariables)
{
    const std::string restartFileName = "OneGlobalDouble.restart";
    const std::string timeStepVarName = "timeStep";
    const double timeStepSize = 1e-6;
    const double currentTime = 1.0;
    MPI_Comm communicator = MPI_COMM_WORLD;

    //Write restart file with time step size as a global variable
    {
        stk::io::StkMeshIoBroker stkIo(communicator);
        generateMetaData(stkIo);
        stkIo.populate_bulk_data();

        size_t fileIndex =
	  stkIo.create_output_mesh(restartFileName, stk::io::WRITE_RESTART);
        stkIo.add_global(fileIndex, timeStepVarName, Ioss::Field::REAL);
        stkIo.begin_output_step(fileIndex, currentTime);
        stkIo.write_global(fileIndex, timeStepVarName, timeStepSize);
        stkIo.end_output_step(fileIndex);
    }

    //Read restart file with time step size as a global variable
    {
        stk::io::StkMeshIoBroker stkIo(communicator);
        stkIo.open_mesh_database(restartFileName, stk::io::READ_RESTART);
        stkIo.create_input_mesh();
        stkIo.populate_bulk_data();
        stkIo.read_defined_input_fields(currentTime);
        std::vector<std::string> globalNamesOnFile;
        stkIo.get_global_variable_names(globalNamesOnFile);

        ASSERT_EQ(1u, globalNamesOnFile.size());
        EXPECT_STRCASEEQ(timeStepVarName.c_str(),
                         globalNamesOnFile[0].c_str());
        double timeStepSizeReadFromFile = 0.0;
	stkIo.get_global(globalNamesOnFile[0], timeStepSizeReadFromFile);
        const double tolerance = 1e-16;
        EXPECT_NEAR(timeStepSize, timeStepSizeReadFromFile, tolerance);
    }
    unlink(restartFileName.c_str());
}
}
