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

namespace {

Ioss::Field::BasicType iossBasicType(double)
{
    return Ioss::Field::REAL;
}

Ioss::Field::BasicType iossBasicType(int)
{
    return Ioss::Field::INTEGER;
}

void generateMetaData(stk::io::MeshData &stkIo)
{
    const std::string exodusFileName = "generated:1x1x1";
    stkIo.open_mesh_database(exodusFileName);
    stkIo.create_input_mesh();
}

template <typename DataType>
void testGlobalVarOnFile(const std::string &outputFileName, const int stepNumber, const std::vector<std::string> &goldGlobalVarName,
                         const std::vector<DataType> goldGlobalVarValue, MPI_Comm comm)
{
    Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", outputFileName, Ioss::READ_RESTART, comm);
    Ioss::Region inputRegion(iossDb);
    Ioss::NameList fieldNames;
    //Get all the fields on the ioss region, which are the global variables
    inputRegion.field_describe(Ioss::Field::TRANSIENT, &fieldNames);

    ASSERT_EQ(goldGlobalVarName.size(), fieldNames.size());

    inputRegion.begin_state(stepNumber);

    const double tolerance = 1e-16;
    for(size_t i=0; i<goldGlobalVarName.size(); i++)
    {
        EXPECT_STRCASEEQ(goldGlobalVarName[i].c_str(), fieldNames[i].c_str());
	std::vector<double> doubleGlobalValue;
        inputRegion.get_field_data(fieldNames[i], doubleGlobalValue);

	Ioss::Field field = inputRegion.get_field(fieldNames[i]);
	size_t num_components = field.transformed_storage()->component_count();
	EXPECT_EQ(num_components, doubleGlobalValue.size());

	for (size_t j=0; j < num_components; j++) {
	  //Workaround exodus only storing doubles
	  DataType globalValue = static_cast<DataType>(doubleGlobalValue[j]);
	  // For a multi-component global variable, the value is 'goldGlobalVarValue' + 0.1 * the component
	  // E.g., a vector_3d with goldGlobalVarValue 42 would have values 42.0 42.1 42.2
	  EXPECT_NEAR(goldGlobalVarValue[i]+0.1*j, globalValue, tolerance);
	}
    }
}

void testNodalFieldOnFile(const std::string &outputFileName, const int stepNumber, const std::string &goldNodalFieldName,
                         const std::vector<double> goldNodalFieldValues, MPI_Comm comm)
{
    Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", outputFileName, Ioss::READ_RESTART, comm);
    Ioss::Region inputRegion(iossDb);
    Ioss::NodeBlock *nodeBlock = inputRegion.get_node_block("nodeblock_1");
    ASSERT_TRUE(nodeBlock->field_exists(goldNodalFieldName));

    inputRegion.begin_state(stepNumber);

    std::vector<double> fieldValues;
    nodeBlock->get_field_data(goldNodalFieldName, fieldValues);
    ASSERT_EQ(goldNodalFieldValues.size(), fieldValues.size());
    const double tolerance = 1e-16;
    for(size_t i=0; i<goldNodalFieldValues.size(); i++)
    {
        EXPECT_NEAR(goldNodalFieldValues[i], fieldValues[i], tolerance);
    }
}

STKUNIT_UNIT_TEST(GlobalVariablesTest, OneGlobalDouble)
{
    const std::string outputFileName = "OneGlobalDouble.exo";
    const std::string globalVarName = "testGlobal";
    const double globalVarValue = 13.0;
    MPI_Comm communicator = MPI_COMM_WORLD;
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);
        stkIo.populate_bulk_data();

        stkIo.create_output_mesh(outputFileName);

        stkIo.add_results_global(globalVarName, Ioss::Field::REAL);

        const double time = 1.0;
        stkIo.begin_results_output_at_time(time);

        stkIo.write_results_global(globalVarName, globalVarValue);

        stkIo.end_current_results_output();
    }

    const int stepNumber = 1;
    std::vector<std::string> globalVarNames(1, globalVarName);
    std::vector<double> globalVarValues(1,globalVarValue);
    testGlobalVarOnFile(outputFileName, stepNumber, globalVarNames, globalVarValues, communicator);
    unlink(outputFileName.c_str());
}

STKUNIT_UNIT_TEST(GlobalVariablesTest, OneGlobalDoubleVector3)
{
    const std::string outputFileName = "OneGlobalDoubleVector3.exo";
    const std::string globalVarName = "testGlobal";
    std::vector<double> globalVarValue;
    globalVarValue.push_back(13.0);
    globalVarValue.push_back(13.1);
    globalVarValue.push_back(13.2);

    MPI_Comm communicator = MPI_COMM_WORLD;
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);
        stkIo.populate_bulk_data();

        stkIo.create_output_mesh(outputFileName);

        stkIo.add_results_global(globalVarName, "vector_3d", Ioss::Field::REAL);

        const double time = 1.0;
        stkIo.begin_results_output_at_time(time);

        stkIo.write_results_global(globalVarName, globalVarValue);

        stkIo.end_current_results_output();
    }

    const int stepNumber = 1;
    std::vector<std::string> globalVarNames(1, globalVarName);
    testGlobalVarOnFile(outputFileName, stepNumber, globalVarNames, globalVarValue, communicator);
    unlink(outputFileName.c_str());
}

STKUNIT_UNIT_TEST(GlobalVariablesTest, OneGlobalDouble10)
{
    const std::string outputFileName = "OneGlobalDouble10.exo";
    const std::string globalVarName = "testGlobal";
    std::vector<double> globalVarValue;
    globalVarValue.push_back(13.0);
    globalVarValue.push_back(13.1);
    globalVarValue.push_back(13.2);
    globalVarValue.push_back(13.3);
    globalVarValue.push_back(13.4);
    globalVarValue.push_back(13.5);
    globalVarValue.push_back(13.6);
    globalVarValue.push_back(13.7);
    globalVarValue.push_back(13.8);
    globalVarValue.push_back(13.9);

    MPI_Comm communicator = MPI_COMM_WORLD;
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);
        stkIo.populate_bulk_data();

        stkIo.create_output_mesh(outputFileName);

        stkIo.add_results_global(globalVarName, globalVarValue.size(), Ioss::Field::REAL);

        const double time = 1.0;
        stkIo.begin_results_output_at_time(time);

        stkIo.write_results_global(globalVarName, globalVarValue);

        stkIo.end_current_results_output();
    }

    const int stepNumber = 1;
    std::vector<std::string> globalVarNames(1, globalVarName);
    testGlobalVarOnFile(outputFileName, stepNumber, globalVarNames, globalVarValue, communicator);
    unlink(outputFileName.c_str());
}

template <typename DataType>
void testTwoGlobals(const std::string &outputFileName, const std::vector<std::string> &globalVarNames)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    std::vector<DataType> globalVarValues;
    globalVarValues.push_back(13);
    globalVarValues.push_back(14);
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);
        stkIo.populate_bulk_data();

        stkIo.create_output_mesh(outputFileName);

        Ioss::Field::BasicType iossDataType = iossBasicType(DataType());
        stkIo.add_results_global(globalVarNames[0], iossDataType);
        stkIo.add_results_global(globalVarNames[1], iossDataType);

        const double time = 1.0;
        stkIo.begin_results_output_at_time(time);

        stkIo.write_results_global(globalVarNames[0], globalVarValues[0]);
        stkIo.write_results_global(globalVarNames[1], globalVarValues[1]);

        stkIo.end_current_results_output();
    }

    const int stepNumber = 1;
    testGlobalVarOnFile(outputFileName, stepNumber, globalVarNames, globalVarValues, communicator);
    unlink(outputFileName.c_str());
}

STKUNIT_UNIT_TEST(GlobalVariablesTest, TwoGlobalIntegers)
{
    std::vector<std::string> globalVarNames;
    globalVarNames.push_back("testGlobal");
    globalVarNames.push_back("testGlobal2");
    testTwoGlobals<int>("TwoGlobalIntegers.exo", globalVarNames);
}

STKUNIT_UNIT_TEST(GlobalVariablesTest, TwoGlobalDoubles)
{
    std::vector<std::string> globalVarNames;
    globalVarNames.push_back("testGlobal");
    globalVarNames.push_back("testGlobal2");
    testTwoGlobals<double>("TwoGlobalDoubles.exo", globalVarNames);
}

STKUNIT_UNIT_TEST(GlobalVariablesTest, TwoGlobalDoublesSameName)
{
    std::vector<std::string> globalVarNames;
    globalVarNames.push_back("testGlobal");
    globalVarNames.push_back("testGlobal");
    EXPECT_THROW(testTwoGlobals<double>("TwoGlobalDoublesSameName.exo", globalVarNames), std::exception);
}

stk::mesh::Field<double> &createNodalTestField(stk::mesh::MetaData &stkMeshMetaData, const std::string &fieldName)
{
    const int numberOfStates = 1;
    stk::mesh::Field<double> &field0 = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(fieldName, numberOfStates);
    stk::mesh::put_field(field0, stk::mesh::Entity::NODE, stkMeshMetaData.universal_part());
    stk::io::set_field_role(field0, Ioss::Field::TRANSIENT);
    return field0;
}
void putDataOnTestField(stk::mesh::BulkData &stkMeshBulkData, stk::mesh::Field<double> &field0, std::vector<double> &nodalFieldValues)
{
    std::vector<stk::mesh::Entity> nodes;
    stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);
    for(size_t i=0; i<nodes.size(); i++)
    {
        double *fieldDataForNode = stkMeshBulkData.field_data(field0, nodes[i]);
        *fieldDataForNode = static_cast<double>(stkMeshBulkData.identifier(nodes[i]));
        nodalFieldValues.push_back(*fieldDataForNode);
    }
}

STKUNIT_UNIT_TEST(GlobalVariablesTest, GlobalDoubleWithFieldMultipleTimeSteps)
{
    const std::string outputFileName = "GlobalDoubleWithFieldMultipleTimeSteps.exo";
    const std::string fieldName = "field0";
    std::vector<double> nodalFieldValues;
    const std::string globalVarName = "testGlobal";
    std::vector<double> globalVarValuesOverTime;
    const int numTimeSteps = 5;
    MPI_Comm communicator = MPI_COMM_WORLD;
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);

        stk::mesh::Field<double> &field0 = createNodalTestField(stkIo.meta_data(), fieldName);

        stkIo.populate_bulk_data();

        putDataOnTestField(stkIo.bulk_data(), field0, nodalFieldValues);

        stkIo.create_output_mesh(outputFileName);

        stkIo.add_results_global(globalVarName, Ioss::Field::REAL);
        stkIo.define_output_fields();

        double time = 1.0;
        const double stepSize = 1.0;
        for(int i=0; i<numTimeSteps; i++)
        {
            stkIo.begin_results_output_at_time(time);

            const double globalVarValue = time;
            stkIo.write_results_global(globalVarName, globalVarValue);
            globalVarValuesOverTime.push_back(globalVarValue);

            stkIo.process_output_request();

            stkIo.end_current_results_output();
            time += stepSize;
        }
    }

    for(int i=0; i<numTimeSteps; i++)
    {
        std::vector<std::string> globalVarNames(1, globalVarName);
        std::vector<double> globalVarValues(1,globalVarValuesOverTime[i]);
        testGlobalVarOnFile(outputFileName, i+1, globalVarNames, globalVarValues, communicator);
        testNodalFieldOnFile(outputFileName, i+1, fieldName, nodalFieldValues, communicator);
    }
    unlink(outputFileName.c_str());
}

STKUNIT_UNIT_TEST(GlobalVariablesTest, OneGlobalDoubleRestart)
{
    const std::string outputFileName = "OneGlobalDouble.restart";
    const std::string globalVarName = "testGlobal";
    const double globalVarValue = 13.0;
    MPI_Comm communicator = MPI_COMM_WORLD;
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);
        stkIo.populate_bulk_data();

        stkIo.create_restart_output(outputFileName);

        stkIo.add_restart_global(globalVarName, Ioss::Field::REAL);

        const double time = 1.0;
        stkIo.begin_restart_output_at_time(time);

        stkIo.write_restart_global(globalVarName, globalVarValue);

        stkIo.end_current_restart_output();
    }

    const int stepNumber = 1;
    std::vector<std::string> globalVarNames(1, globalVarName);
    std::vector<double> globalVarValues(1,globalVarValue);
    testGlobalVarOnFile(outputFileName, stepNumber, globalVarNames, globalVarValues, communicator);
    unlink(outputFileName.c_str());
}

STKUNIT_UNIT_TEST(GlobalVariablesTest, OneGlobalDoubleWithFieldRestart)
{
    const std::string outputFileName = "GlobalDoubleWithFieldMultipleTimeSteps.restart";
    const std::string fieldName = "field0";
    std::vector<double> nodalFieldValues;
    const std::string globalVarName = "testGlobal";
    std::vector<double> globalVarValuesOverTime;
    const int numTimeSteps = 5;
    MPI_Comm communicator = MPI_COMM_WORLD;
    {
        stk::io::MeshData stkIo(communicator);
        generateMetaData(stkIo);

        stk::mesh::Field<double> &field0 = createNodalTestField(stkIo.meta_data(), fieldName);
        stkIo.add_restart_field(field0);

        stkIo.populate_bulk_data();

        putDataOnTestField(stkIo.bulk_data(), field0, nodalFieldValues);

        stkIo.create_restart_output(outputFileName);

        stkIo.add_restart_global(globalVarName, Ioss::Field::REAL);
        stkIo.define_restart_fields();

        double time = 1.0;
        const double stepSize = 1.0;
        for(int i=0; i<numTimeSteps; i++)
        {
            stkIo.begin_restart_output_at_time(time);

            const double globalVarValue = time;
            stkIo.write_restart_global(globalVarName, globalVarValue);
            globalVarValuesOverTime.push_back(globalVarValue);

            stkIo.process_restart_output();

            stkIo.end_current_restart_output();
            time += stepSize;
        }
    }

    for(int i=0; i<numTimeSteps; i++)
    {
        std::vector<std::string> globalVarNames(1, globalVarName);
        std::vector<double> globalVarValues(1,globalVarValuesOverTime[i]);
        testGlobalVarOnFile(outputFileName, i+1, globalVarNames, globalVarValues, communicator);
        testNodalFieldOnFile(outputFileName, i+1, fieldName, nodalFieldValues, communicator);
    }
    unlink(outputFileName.c_str());
}

}
