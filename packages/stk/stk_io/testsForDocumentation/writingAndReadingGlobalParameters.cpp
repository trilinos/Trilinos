#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <mpi.h>
#include <stk_util/util/ParameterList.hpp>
#include <stk_io/MeshReadWriteUtils.hpp>
#include <fieldNameTestUtils.hpp>

namespace
{
TEST(StkIoTestForDocumentation, WriteAndReadGlobalParameters)
{
    const std::string file_name = "GlobalParameters.e";
    MPI_Comm communicator = MPI_COMM_WORLD;

    stk::util::ParameterList parameters;
    stk::util::ParameterList gold_parameters; // To compare values read
    
    // Add some parameters to write and read...
    parameters.set_param("PI", 3.14159);  // Double 
    parameters.set_param("Answer", 42);   // Integer
    gold_parameters.set_param("PI", 3.14159);  // Double 
    gold_parameters.set_param("Answer", 42);   // Integer

    std::vector<double> my_vector;
    my_vector.push_back(2.78);
    my_vector.push_back(5.30);
    my_vector.push_back(6.21);
    parameters.set_param("some_doubles", my_vector);   // Vector of doubles...
    gold_parameters.set_param("some_doubles", my_vector);   // Vector of doubles...
    
    std::vector<int> ages;
    ages.push_back(55);
    ages.push_back(49);
    ages.push_back(21);
    ages.push_back(19);
    
    parameters.set_param("Ages", ages);   // Vector of integers...
    gold_parameters.set_param("Ages", ages);   // Vector of integers...
    
    // Write output file with all parameters in parameters list...
    {
        stk::io::StkMeshIoBroker stkIo(communicator);
        generateMetaData(stkIo);
        stkIo.populate_bulk_data();

        size_t fileIndex = stkIo.create_output_mesh(file_name);

	stk::util::ParameterMapType::const_iterator i = parameters.begin();
	stk::util::ParameterMapType::const_iterator iend = parameters.end();
	for (; i != iend; ++i) {
	  const std::string parameterName = (*i).first;
	  stk::util::Parameter &parameter = parameters.get_param(parameterName);
	  stkIo.add_global(fileIndex, parameterName, parameter);
	}

        stkIo.begin_output_at_time(0.0, fileIndex);

	for (i = parameters.begin(); i != iend; ++i) {
	  const std::string parameterName = (*i).first;
	  stk::util::Parameter &parameter = parameters.get_param(parameterName);
	  stkIo.write_global(fileIndex, parameterName, parameter);
	}

        stkIo.end_current_output(fileIndex);
    }

    // Read parameters from file...
    {
        stk::io::StkMeshIoBroker stkIo(communicator);
        stkIo.open_mesh_database(file_name);
        stkIo.create_input_mesh();
        stkIo.populate_bulk_data();

        stkIo.process_restart_input(0.0);

	size_t param_count = 0;
	stk::util::ParameterMapType::const_iterator i = parameters.begin();
	stk::util::ParameterMapType::const_iterator iend = parameters.end();
	for (; i != iend; ++i) {
	  param_count++;
	  const std::string parameterName = (*i).first;
	  stk::util::Parameter &parameter = parameters.get_param(parameterName);
	  stk::util::Parameter &gold_parameter = gold_parameters.get_param(parameterName);
	  stkIo.get_global(parameterName, parameter);
	  validate_parameters_equal_value(parameter, gold_parameter);
	}

        std::vector<std::string> globalNamesOnFile;
        stkIo.get_global_variable_names(globalNamesOnFile);
        ASSERT_EQ(param_count, globalNamesOnFile.size());

    }
    unlink(file_name.c_str());
}
}
