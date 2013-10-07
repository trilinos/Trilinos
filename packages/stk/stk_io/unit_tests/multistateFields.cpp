#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <string>
#include <mpi.h>
#include <stk_io/MeshReadWriteUtils.hpp>
#include <stk_io/IossBridge.hpp>
#include <Ioss_SubSystem.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <fieldNameTestUtils.hpp>

namespace {

stk::mesh::FieldBase* declareTriStateNodalField(stk::mesh::MetaData &stkMeshMetaData, const std::string &fieldName)
{
    const int numberOfStates = 3;
    stk::mesh::Field<double> &triStateField = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(fieldName, numberOfStates);
    stk::mesh::put_field(triStateField, stk::mesh::Entity::NODE, stkMeshMetaData.universal_part());
    return &triStateField;
}

void putDataOnTestField(stk::mesh::BulkData &stkMeshBulkData, const double value, stk::mesh::FieldBase &field)
{
    std::vector<stk::mesh::Entity> nodes;
    stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);
    for(size_t i=0; i<nodes.size(); i++)
    {
        double *fieldDataForNode = reinterpret_cast<double*>(stkMeshBulkData.field_data(field, nodes[i]));
        *fieldDataForNode = value;
    }
}

void testDataOnField(stk::mesh::BulkData &stkMeshBulkData, const double goldValue, stk::mesh::FieldBase &field)
{
    std::vector<stk::mesh::Entity> nodes;
    stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);
    for(size_t i=0; i<nodes.size(); i++)
    {
        double *fieldDataForNode = reinterpret_cast<double*>(stkMeshBulkData.field_data(field, nodes[i]));
        EXPECT_DOUBLE_EQ(goldValue, *fieldDataForNode);
    }
}

STKUNIT_UNIT_TEST(MultistateFieldTest, SingleNodalMultistatedRestartField)
{
    std::string restartFilename = "output.restart";
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string fieldName = "multistateField";
    double time = 0.0;
    {
        stk::io::MeshData stkIo(communicator);
        stk::mesh::MetaData &stkMeshMetaData = generateMetaData(stkIo);
        stk::mesh::FieldBase *triStateField = declareTriStateNodalField(stkMeshMetaData, fieldName);
        stkIo.populate_bulk_data();

        stk::mesh::FieldBase *statedFieldNp1 = triStateField->field_state(stk::mesh::StateNP1);
        putDataOnTestField(stkIo.bulk_data(), 1.0, *statedFieldNp1);

        stk::mesh::FieldBase *statedFieldN = triStateField->field_state(stk::mesh::StateN);
        putDataOnTestField(stkIo.bulk_data(), 2.0, *statedFieldN);

        stk::mesh::FieldBase *statedFieldNm1 = triStateField->field_state(stk::mesh::StateNM1);
        putDataOnTestField(stkIo.bulk_data(), 3.0, *statedFieldNm1);

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
        stk::mesh::FieldBase *triStateField = declareTriStateNodalField(restartedMetaData, fieldName);
        stkIo.add_restart_field(*triStateField);

        stkIo.populate_bulk_data();

        stkIo.process_restart_input(time);

        stk::mesh::FieldBase *statedFieldN = triStateField->field_state(stk::mesh::StateN);
        testDataOnField(stkIo.bulk_data(), 1.0, *statedFieldN);
        stk::mesh::FieldBase *statedFieldNm1 = triStateField->field_state(stk::mesh::StateNM1);
        testDataOnField(stkIo.bulk_data(), 2.0, *statedFieldNm1);
    }
    unlink(restartFilename.c_str());
}
}
