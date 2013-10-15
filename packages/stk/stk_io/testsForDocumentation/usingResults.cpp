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

TEST(StkIoTestForDocumentation, resultsWithMultistateField)
{
    std::string resultsFilename = "output.results";
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string fieldName = "disp";
    const double stateNp1Value = 1.0;
    const double stateNValue = 2.0;
    const double stateNm1Value = 3.0;
    double time = 0.0;

    const std::string np1Name = fieldName+"NP1";
    const std::string nName   = fieldName+"N";
    const std::string nm1Name = fieldName+"Nm1";
    {
        stk::io::MeshData stkIo(communicator);

        stk::mesh::MetaData &stkMeshMetaData = generateMetaData(stkIo);
	// Declare a multi-state field
        stk::mesh::FieldBase *triStateField =
	  declareNodalField(stkMeshMetaData, fieldName, 3);

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

	// Output each state of the multi-state field individually to results file.
	stkIo.add_results_field(*statedFieldNp1, np1Name);
	stkIo.add_results_field(*statedFieldN,   nName);
	stkIo.add_results_field(*statedFieldNm1, nm1Name);
	
        stkIo.create_output_mesh(resultsFilename);

	stkIo.begin_results_output_at_time(time);
	stkIo.process_output_request();
        stkIo.end_current_results_output();
    }

    {
        stk::io::MeshData stkIo(communicator);
        stkIo.open_mesh_database(resultsFilename);
        stkIo.create_input_mesh();

        stk::mesh::MetaData &resultsedMetaData = stkIo.meta_data();
        stk::mesh::FieldBase *FieldNp1 = declareNodalField(resultsedMetaData, np1Name, 1);
        stk::mesh::FieldBase *FieldN   = declareNodalField(resultsedMetaData, nName, 1);
        stk::mesh::FieldBase *FieldNm1 = declareNodalField(resultsedMetaData, nm1Name, 1);

        stkIo.add_restart_field(*FieldNp1, np1Name);
        stkIo.add_restart_field(*FieldN,   nName);
        stkIo.add_restart_field(*FieldNm1, nm1Name);

        stkIo.populate_bulk_data();
        stkIo.process_restart_input(time);

        testDataOnField(stkIo.bulk_data(), stateNp1Value, *FieldNp1);
        testDataOnField(stkIo.bulk_data(), stateNValue,   *FieldN);
        testDataOnField(stkIo.bulk_data(), stateNm1Value, *FieldNm1);
    }
    unlink(resultsFilename.c_str());
}
}
