#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/MeshReadWriteUtils.hpp>
#include <stk_io/IossBridge.hpp>
#include <Ioss_SubSystem.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <fieldNameTestUtils.hpp>
#include <restartTestUtils.hpp>

namespace {

TEST(StkIoTestForDocumentation, restartingWithMultistateField)
{
    std::string restartFilename = "output.restart";
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string fieldName = "disp";
    const double stateNp1Value = 1.0;
    const double stateNValue = 2.0;
    const double stateNm1Value = 3.0;
    double time = 0.0;
    {
        stk::io::MeshData stkIo(communicator);
        stk::mesh::MetaData &stkMeshMetaData = generateMetaData(stkIo);
        stk::mesh::FieldBase *triStateField =
                declareTriStateNodalField(stkMeshMetaData, fieldName);
        stkIo.populate_bulk_data();

        stk::mesh::FieldBase *statedFieldNp1 =
                triStateField->field_state(stk::mesh::StateNP1);
        putDataOnTestField(stkIo.bulk_data(), stateNp1Value,
                           *statedFieldNp1);
        stk::mesh::FieldBase *statedFieldN =
                triStateField->field_state(stk::mesh::StateN);
        putDataOnTestField(stkIo.bulk_data(), stateNValue,
                           *statedFieldN);
        stk::mesh::FieldBase *statedFieldNm1 =
                triStateField->field_state(stk::mesh::StateNM1);
        putDataOnTestField(stkIo.bulk_data(), stateNm1Value,
                           *statedFieldNm1);

        stkIo.add_restart_field(*triStateField);

        stkIo.create_restart_output(restartFilename);
        stkIo.define_restart_fields();

        stkIo.begin_restart_output_at_time(time);
        stkIo.process_restart_output();
        stkIo.end_current_restart_output();
    }

    {
        stk::io::MeshData stkIo(communicator);
        stkIo.open_mesh_database(restartFilename);
        stkIo.create_input_mesh();

        stk::mesh::MetaData &restartedMetaData = stkIo.meta_data();
        stk::mesh::FieldBase *triStateField =
                declareTriStateNodalField(restartedMetaData, fieldName);

        stkIo.add_restart_field(*triStateField);
        stkIo.populate_bulk_data();
        stkIo.process_restart_input(time);

        stk::mesh::FieldBase *statedFieldNp1 =
                triStateField->field_state(stk::mesh::StateNP1);
        testDataOnField(stkIo.bulk_data(), stateNp1Value, *statedFieldNp1);
        stk::mesh::FieldBase *statedFieldN =
                triStateField->field_state(stk::mesh::StateN);
        testDataOnField(stkIo.bulk_data(), stateNValue, *statedFieldN);
    }
    unlink(restartFilename.c_str());
}
}
