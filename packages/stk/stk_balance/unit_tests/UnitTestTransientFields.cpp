#include <stk_unit_test_utils/unittestMeshUtils.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/internal/StkBalanceUtils.hpp>
#include "stk_balance/internal/TransientFieldTransferById.hpp"
#include <stk_balance/internal/balanceDefaults.hpp>

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include "stk_io/FillMesh.hpp"
#include "stk_io/WriteMesh.hpp"
#include <stk_io/StkMeshIoBroker.hpp>

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_io/StkIoUtils.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_mesh/base/FieldParallel.hpp"
#include <stk_tools/mesh_clone/MeshClone.hpp>
#include <iostream>
#include <limits>

namespace
{

void unlink_serial_file(const std::string &baseName)
{
    unlink(baseName.c_str());
}

void unlink_parallel_file(const std::string &baseName)
{
    std::string file = stk::balance::internal::get_parallel_filename(stk::parallel_machine_rank(MPI_COMM_WORLD),
                                                                     stk::parallel_machine_size(MPI_COMM_WORLD),
                                                                     baseName);

    unlink(file.c_str());
}


void verify_expected_time_steps_in_parallel_file(const std::string &parallelMeshName, const std::vector<double>& expectedTimeSteps)
{
    stk::io::StkMeshIoBroker broker;
    stk::mesh::MetaData metaData;
    stk::mesh::BulkData bulkData(metaData, MPI_COMM_WORLD);
    stk::io::fill_mesh_preexisting(broker, parallelMeshName, bulkData);
    std::vector<double> inputTimeSteps = broker.get_time_steps();
    EXPECT_EQ(expectedTimeSteps, inputTimeSteps);
}

void verify_expected_global_variable_in_parallel_file(const std::string &parallelMeshName, const std::vector<double>& expectedTimeSteps,
                                                      const std::string& globalVariableName)
{
    stk::io::StkMeshIoBroker broker;
    stk::mesh::MetaData metaData;
    stk::mesh::BulkData bulkData(metaData, MPI_COMM_WORLD);
    stk::io::fill_mesh_preexisting(broker, parallelMeshName, bulkData);
    std::vector<std::string> names;
    broker.get_global_variable_names(names);
    ASSERT_EQ(3u, names.size());
    EXPECT_EQ(globalVariableName+"_double", names[0]);
    EXPECT_EQ(globalVariableName+"_int", names[1]);
    EXPECT_EQ(globalVariableName+"_real_vec", names[2]);

    std::vector<double> inputTimeSteps = broker.get_time_steps();
    EXPECT_EQ(expectedTimeSteps, inputTimeSteps);
    double testTimeValue = 0;
    const double epsilon = std::numeric_limits<double>::epsilon();
    for(size_t i=0; i<expectedTimeSteps.size(); ++i) {
        int timestep = i+1;
        broker.read_defined_input_fields(timestep);
        bool success = broker.get_global(globalVariableName+"_double", testTimeValue);
        EXPECT_TRUE(success);
        EXPECT_NEAR(expectedTimeSteps[i], testTimeValue, epsilon);

        int goldIntValue = timestep;
        int testIntValue = 0;
        success = broker.get_global(globalVariableName+"_int", testIntValue);
        EXPECT_TRUE(success);
        EXPECT_EQ(goldIntValue, testIntValue);

        std::vector<double> goldVecValue(3, expectedTimeSteps[i]);
        std::vector<double> testVecValue;
        success = broker.get_global(globalVariableName+"_real_vec", testVecValue);
        EXPECT_EQ(goldVecValue, testVecValue);
    }
}

void verify_input_transient_fields(stk::io::StkMeshIoBroker &stkIo, const std::vector<double>& expectedTimeSteps)
{
    EXPECT_EQ(stkIo.get_num_time_steps(), (int)expectedTimeSteps.size());
    EXPECT_EQ(expectedTimeSteps, stkIo.get_time_steps());

    stk::io::FieldNameToPartVector fieldNamePartVector = stkIo.get_nodal_var_names();

    EXPECT_EQ(fieldNamePartVector.size(), 2u);
}

void verify_input_static_fields(stk::io::StkMeshIoBroker &stkIo)
{
    EXPECT_EQ(0, stkIo.get_num_time_steps());

    stk::io::FieldNameToPartVector fieldNamePartVector = stkIo.get_nodal_var_names();

    EXPECT_EQ(fieldNamePartVector.size(), 0u);
}

void verify_transient_field_values(stk::mesh::BulkData& bulk, stk::mesh::FieldBase* field, double timeStep)
{
    const stk::mesh::BucketVector & entityBuckets = bulk.get_buckets(field->entity_rank(),bulk.mesh_meta_data().locally_owned_part());
    for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
        stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
        for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
            stk::mesh::Entity entity = entityBucket[entityIndex];
            double * data = static_cast<double*> (stk::mesh::field_data(*field, entity));
            unsigned numEntriesPerEntity = stk::mesh::field_scalars_per_entity(*field, entity);
            for(unsigned i=0; i<numEntriesPerEntity; i++)
                EXPECT_EQ(i + 100*timeStep + static_cast<double>(bulk.identifier(entity)), data[i]);
        }
    }
}

void verify_transient_data_from_parallel_file(const std::string &inputParallelFile, const std::vector<double>& expectedTimeSteps)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
    stk::io::StkMeshIoBroker broker;
    stk::io::fill_mesh_preexisting(broker, inputParallelFile, bulk);

    stk::mesh::FieldVector transientFields = stk::io::get_transient_fields(meta);

    verify_input_transient_fields(broker, expectedTimeSteps);
    std::vector<double> timeSteps = broker.get_time_steps();
    for(int iStep=0; iStep<broker.get_num_time_steps(); iStep++)
    {
        double readTime = broker.read_defined_input_fields_at_step(iStep+1, nullptr);
        EXPECT_EQ(timeSteps[iStep], readTime);

        for(stk::mesh::FieldBase* field : transientFields)
            verify_transient_field_values(bulk, field, readTime);
    }
}

void verify_static_data_from_parallel_file(const std::string &inputParallelFile)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
    stk::io::StkMeshIoBroker broker;
    stk::io::fill_mesh_preexisting(broker, inputParallelFile, bulk);

    verify_input_static_fields(broker);
}

void verify_static_data_from_serial_file(const std::string &inputSerialFile)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
    stk::io::StkMeshIoBroker broker;
    broker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "LINEAR"));
    stk::io::fill_mesh_preexisting(broker, inputSerialFile, bulk);

    verify_input_static_fields(broker);
}

void verify_transient_field_name(stk::io::StkMeshIoBroker &brokerA, stk::io::StkMeshIoBroker &brokerB, const std::string &fieldBaseName)
{
    EXPECT_EQ(brokerA.meta_data().entity_rank_names(), brokerB.meta_data().entity_rank_names());

    stk::mesh::FieldVector transientFieldsA = stk::io::get_transient_fields(brokerA.meta_data(), stk::topology::NODE_RANK);
    stk::mesh::FieldVector transientFieldsB = stk::io::get_transient_fields(brokerB.meta_data(), stk::topology::NODE_RANK);

    EXPECT_EQ(2u, transientFieldsA.size());
    EXPECT_EQ(2u, transientFieldsB.size());

    std::string scalarFieldName = fieldBaseName + "_scalar";
    std::string vectorFieldName = fieldBaseName + "_vector";

    EXPECT_EQ(0, strcasecmp(scalarFieldName.c_str(), transientFieldsA[0]->name().c_str())) << scalarFieldName << "   " << transientFieldsA[0]->name();
    EXPECT_EQ(0, strcasecmp(vectorFieldName.c_str(), transientFieldsA[1]->name().c_str())) << vectorFieldName << "   " << transientFieldsA[1]->name();

    EXPECT_EQ(0, strcasecmp(scalarFieldName.c_str(), transientFieldsB[0]->name().c_str())) << scalarFieldName << "   " << transientFieldsB[0]->name();
    EXPECT_EQ(0, strcasecmp(vectorFieldName.c_str(), transientFieldsB[1]->name().c_str())) << vectorFieldName << "   " << transientFieldsB[1]->name();
}

void verify_transient_transfer(const std::string& serialOutputMeshName,
                               const std::string& parallelOutputMeshName,
                               const std::string& fieldBaseName,
                               const std::vector<double>& expectedTimeSteps)
{
    stk::mesh::MetaData metaA;
    stk::mesh::BulkData bulkA(metaA, MPI_COMM_WORLD);
    stk::io::StkMeshIoBroker brokerA;
    brokerA.property_add(Ioss::Property("DECOMPOSITION_METHOD", "LINEAR"));
    stk::io::fill_mesh_preexisting(brokerA, serialOutputMeshName, bulkA);

    stk::mesh::MetaData metaB;
    stk::mesh::BulkData bulkB(metaB, MPI_COMM_WORLD);
    stk::io::StkMeshIoBroker brokerB;
    stk::io::fill_mesh_preexisting(brokerB, parallelOutputMeshName, bulkB);

    verify_input_transient_fields(brokerA, expectedTimeSteps);
    verify_transient_field_name(brokerA, brokerB, fieldBaseName);
    verify_transient_data_from_parallel_file(parallelOutputMeshName, expectedTimeSteps);
}

void verify_static_transfer(const std::string& serialOutputMeshName,
                            const std::string& parallelOutputMeshName)
{
    verify_static_data_from_serial_file(serialOutputMeshName);
    verify_static_data_from_parallel_file(parallelOutputMeshName);
}

TEST(TestTransientFieldBalance, verifyNumberOfSteps)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
    {
        std::string outputMeshName = "twenty_hex_transient.e";
        std::string fieldBaseName = "myField";
        std::string globalVariableName = "TestTime";
        const std::vector<double> expectedTimeSteps = {0.0, 1.0, 2.0, 3.0, 4.0};

        stk::unit_test_util::IdAndTimeFieldValueSetter fieldSetter;
        stk::unit_test_util::generated_mesh_with_transient_data_to_file_in_serial("1x1x20", outputMeshName, fieldBaseName, globalVariableName, expectedTimeSteps, fieldSetter);

        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
        stk::io::StkMeshIoBroker broker;
        stk::io::fill_mesh_with_auto_decomp(outputMeshName, bulk, broker);
        verify_input_transient_fields(broker, expectedTimeSteps);

        stk::balance::run_stk_rebalance(".", outputMeshName, stk::balance::SD_DEFAULTS, MPI_COMM_WORLD);

        verify_expected_time_steps_in_parallel_file(outputMeshName, expectedTimeSteps);

        unlink_serial_file(outputMeshName.c_str());
        unlink_parallel_file(outputMeshName.c_str());
    }
}

TEST(TestTransientFieldBalance, verifyTransientDataTransferOnFourProcessors)
{
    std::string serialOutputMeshName = "sixteen_hex_transient.e";
    std::string parallelOutputMeshName = "sixteen_hex_transient_balanced.e";
    std::string fieldBaseName = "myField";
    std::string globalVariableName = "TestTime";
    const std::vector<double> expectedTimeSteps = {0.0, 1.0, 2.0, 3.0, 4.0};

    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 4)
    {
        stk::unit_test_util::IdAndTimeFieldValueSetter fieldSetter;
        stk::unit_test_util::generated_mesh_with_transient_data_to_file_in_serial("1x4x4", serialOutputMeshName, fieldBaseName, globalVariableName, expectedTimeSteps, fieldSetter);

        stk::balance::BasicZoltan2Settings rcbOptions;
        stk::balance::run_stk_balance_with_settings(parallelOutputMeshName, serialOutputMeshName, MPI_COMM_WORLD, rcbOptions);

        verify_transient_transfer(serialOutputMeshName, parallelOutputMeshName, fieldBaseName, expectedTimeSteps);

        unlink_serial_file(serialOutputMeshName);
        unlink_parallel_file(parallelOutputMeshName);
    }
}

TEST(TestTransientFieldBalance, verifyStaticDataTransfer)
{
    std::string serialOutputMeshName = "sixteen_hex_transient.e";
    std::string parallelOutputMeshName = "sixteen_hex_transient_balanced.e";

    if((stk::parallel_machine_size(MPI_COMM_WORLD) ==  2) ||
       (stk::parallel_machine_size(MPI_COMM_WORLD) ==  4) ||
       (stk::parallel_machine_size(MPI_COMM_WORLD) ==  8) ||
       (stk::parallel_machine_size(MPI_COMM_WORLD) == 16))
    {
        stk::unit_test_util::generated_mesh_to_file_in_serial("1x4x4", serialOutputMeshName);

        stk::balance::BasicZoltan2Settings rcbOptions;
        stk::balance::run_stk_balance_with_settings(parallelOutputMeshName, serialOutputMeshName, MPI_COMM_WORLD, rcbOptions);

        verify_static_transfer(serialOutputMeshName, parallelOutputMeshName);

        unlink_serial_file(serialOutputMeshName);
        unlink_parallel_file(parallelOutputMeshName);
    }
}

TEST(TestTransientFieldBalance, verifyGlobalVariable)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
    {
        std::string outputMeshName = "twenty_hex_transient.e";
        std::string fieldBaseName = "myField";
        const std::vector<double> expectedTimeSteps = {0.0, 1.0, 2.0, 3.0, 4.0};
        std::string globalVariableName = "test_time";
        stk::unit_test_util::IdAndTimeFieldValueSetter fieldSetter;
        stk::unit_test_util::generated_mesh_with_transient_data_to_file_in_serial("1x1x20", outputMeshName, fieldBaseName, globalVariableName, expectedTimeSteps, fieldSetter);

        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
        stk::io::StkMeshIoBroker broker;
        stk::io::fill_mesh_with_auto_decomp(outputMeshName, bulk, broker);
        verify_input_transient_fields(broker, expectedTimeSteps);

        stk::balance::run_stk_rebalance(".", outputMeshName, stk::balance::SD_DEFAULTS, MPI_COMM_WORLD);

        verify_expected_global_variable_in_parallel_file(outputMeshName, expectedTimeSteps, globalVariableName);

        unlink_serial_file(outputMeshName.c_str());
        unlink_parallel_file(outputMeshName.c_str());
    }
}

TEST(TestTransientFieldBalance, verifyThrowIfInputFileEqualsOutputFile)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    std::string serialMeshName = "sixteen_hex_transient.e";
    std::string parallelOutputMeshName = "sixteen_hex_transient.e";

    stk::unit_test_util::generated_mesh_to_file_in_serial("1x4x4", serialMeshName);

    stk::balance::BasicZoltan2Settings rcbOptions;
    EXPECT_THROW(stk::balance::run_stk_balance_with_settings(parallelOutputMeshName, serialMeshName, MPI_COMM_WORLD, rcbOptions), std::logic_error);
    unlink_serial_file(serialMeshName);
  }
}

TEST(TestTransientFieldBalance, verifyThrowIfInputFileEqualsDotSlashOutputFile)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    std::string serialMeshName = "sixteen_hex_transient.e";
    std::string parallelOutputMeshName = "./sixteen_hex_transient.e";

    stk::unit_test_util::generated_mesh_to_file_in_serial("1x4x4", serialMeshName);

    stk::balance::BasicZoltan2Settings rcbOptions;
    EXPECT_THROW(stk::balance::run_stk_balance_with_settings(parallelOutputMeshName, serialMeshName, MPI_COMM_WORLD, rcbOptions), std::logic_error);
    unlink_serial_file(serialMeshName);
  }
}

}
