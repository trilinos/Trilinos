#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>
#include <Ioss_SubSystem.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <fieldNameTestUtils.hpp>
#include <restartTestUtils.hpp>

namespace {

TEST(StkMeshIoBrokerHowTo, restartWithMultistateField)
{
    std::string restartFilename = "output.restart";
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string fieldName = "disp";
    const double stateNp1Value = 1.0;
    const double stateNValue = 2.0;
    const double stateNm1Value = 3.0;
    double time = 0.0;
    {
        stk::io::StkMeshIoBroker stkIo(communicator);
        stk::mesh::MetaData &stkMeshMetaData = generateMetaData(stkIo);
        stk::mesh::FieldBase *triStateField =
                declareTriStateNodalField(stkMeshMetaData, fieldName);
        stkIo.populate_bulk_data();

        putDataOnTriStateField(stkIo.bulk_data(), triStateField,
                stateNp1Value, stateNValue, stateNm1Value);

        size_t fileHandle = stkIo.create_output_mesh(restartFilename);
        stkIo.add_restart_field(fileHandle, *triStateField);

        stkIo.begin_output_step(time, fileHandle);
        stkIo.process_output_request(fileHandle);
        stkIo.end_output_step(fileHandle);
    }

    //code to test that the field wrote correctly
    testMultistateFieldWroteCorrectlyToRestart(restartFilename, time,
            fieldName, stateNp1Value, stateNValue);
    unlink(restartFilename.c_str());
}
}
