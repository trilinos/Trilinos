#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetEntities.hpp>

inline void checkFileForNodalVarNames(const std::string &exodusFilename, const std::vector<std::string>& nodalVarNames)
{
   Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", exodusFilename, Ioss::READ_MODEL, MPI_COMM_WORLD);
   Ioss::Region ioRegion(iossDb);
   Ioss::NodeBlock *nodeBlockAssociatedWithField0 = ioRegion.get_node_blocks()[0];
   Ioss::NameList fieldNames;
   nodeBlockAssociatedWithField0->field_describe(Ioss::Field::TRANSIENT, &fieldNames);

   ASSERT_EQ(nodalVarNames.size(), fieldNames.size());
   for (size_t i=0;i<nodalVarNames.size();i++)
   {
       EXPECT_TRUE(nodeBlockAssociatedWithField0->field_exists(nodalVarNames[i])) << nodalVarNames[i];
   }
}

inline void checkFileForGlobal(const std::string &exodusFilename, const std::string &globalVarName, const double value)
{
    Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", exodusFilename, Ioss::READ_MODEL, MPI_COMM_WORLD);
    Ioss::Region ioRegion(iossDb);
    ioRegion.begin_state(1);
    ASSERT_TRUE(ioRegion.field_exists(globalVarName));
    double valueOnFile = 0.0;
    ioRegion.get_field_data(globalVarName, &valueOnFile, sizeof(double));
    EXPECT_EQ(value, valueOnFile);
}

inline stk::mesh::FieldBase* declareNodalField(stk::mesh::MetaData &stkMeshMetaData, const std::string &fieldName,
					       int numberOfStates)
{
    stk::mesh::Field<double> &multiStateField = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(fieldName, numberOfStates);
    stk::mesh::put_field(multiStateField, stk::mesh::Entity::NODE, stkMeshMetaData.universal_part());
    return &multiStateField;
}

inline stk::mesh::FieldBase* declareTriStateNodalField(stk::mesh::MetaData &stkMeshMetaData, const std::string &fieldName)
{
    const int numberOfStates = 3;
    return declareNodalField(stkMeshMetaData, fieldName, numberOfStates);
}

inline void putDataOnTestField(stk::mesh::BulkData &stkMeshBulkData, const double value, stk::mesh::FieldBase &field)
{
    std::vector<stk::mesh::Entity> nodes;
    stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);
    for(size_t i=0; i<nodes.size(); i++)
    {
        double *fieldDataForNode = reinterpret_cast<double*>(stkMeshBulkData.field_data(field, nodes[i]));
        *fieldDataForNode = value;
    }
}

inline void putDataOnTriStateField(stk::mesh::BulkData &bulkData, stk::mesh::FieldBase *triStateField,
        const double stateNp1Value,
        const double stateNValue,
        const double stateNm1Value)
{
    stk::mesh::FieldBase *statedFieldNp1 =
            triStateField->field_state(stk::mesh::StateNP1);
    putDataOnTestField(bulkData, stateNp1Value,
                       *statedFieldNp1);
    stk::mesh::FieldBase *statedFieldN =
            triStateField->field_state(stk::mesh::StateN);
    putDataOnTestField(bulkData, stateNValue,
                       *statedFieldN);
    stk::mesh::FieldBase *statedFieldNm1 =
            triStateField->field_state(stk::mesh::StateNM1);
    putDataOnTestField(bulkData, stateNm1Value,
                       *statedFieldNm1);
}

inline void testDataOnField(stk::mesh::BulkData &stkMeshBulkData, const double goldValue, stk::mesh::FieldBase &field)
{
    std::vector<stk::mesh::Entity> nodes;
    stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);
    for(size_t i=0; i<nodes.size(); i++)
    {
        double *fieldDataForNode = reinterpret_cast<double*>(stkMeshBulkData.field_data(field, nodes[i]));
        EXPECT_DOUBLE_EQ(goldValue, *fieldDataForNode);
    }
}

inline void createExampleInputFile(stk::io::StkMeshIoBroker& stkMeshIoBroker)
{
    const std::string exodusFileName = "generated:1x1x1";
    stkMeshIoBroker.open_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
}

inline void setupMeshAndFieldsForTest(stk::io::StkMeshIoBroker &stkMeshIoBroker,
        const std::string &displacementFieldName,
        const std::string &velocityFieldName)
{
    const double displacementValue = 1.0;
    const double velocityValue = 2.0;

    createExampleInputFile(stkMeshIoBroker);

    const int numberOfStates = 1;
    stk::mesh::FieldBase *displacementField = declareNodalField(stkMeshIoBroker.meta_data(), displacementFieldName, numberOfStates);
    stk::mesh::FieldBase *velocityField = declareNodalField(stkMeshIoBroker.meta_data(), velocityFieldName, numberOfStates);

    stkMeshIoBroker.populate_bulk_data();

    putDataOnTestField(stkMeshIoBroker.bulk_data(), displacementValue, *displacementField);
    putDataOnTestField(stkMeshIoBroker.bulk_data(), velocityValue, *velocityField);
}

inline void testMultistateFieldWroteCorrectlyToRestart(const std::string &restartFilename,
        const double time,
        const std::string &fieldName,
        const double stateNp1Value,
        const double stateNValue)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.open_mesh_database(restartFilename, stk::io::READ_RESTART);
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

inline void testMultistateFieldWroteCorrectly(const std::string &resultsFilename,
        const double time,
        const std::string &np1Name,
        const std::string &nName,
        const std::string &nm1Name,
        const double stateNp1Value,
        const double stateNValue,
        const double stateNm1Value)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.open_mesh_database(resultsFilename, stk::io::READ_RESTART);
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
