#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <mpi.h>
#include <stk_util/util/ParameterList.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <fieldNameTestUtils.hpp>

namespace
{
  TEST(StkMeshIoBrokerHowTo, writeAndReadGlobalParametersAuto)
  {
    const std::string file_name = "GlobalParameters.e";
    MPI_Comm communicator = MPI_COMM_WORLD;

    stk::util::ParameterList params;
    
    // Add some parameters to write and read...
    params.set_param("PI", 3.14159);  // Double 
    params.set_param("Answer", 42);   // Integer

    std::vector<double> my_vector;
    my_vector.push_back(2.78);
    my_vector.push_back(5.30);
    my_vector.push_back(6.21);
    params.set_param("some_doubles", my_vector);   // Vector of doubles...
    
    std::vector<int> ages;
    ages.push_back(55);
    ages.push_back(49);
    ages.push_back(21);
    ages.push_back(19);
    
    params.set_param("Ages", ages);   // Vector of integers...
    
    //-BEGIN
    // ... Setup is the same as in the previous example
    // Write output file with all parameters in params list...
    {
      stk::io::StkMeshIoBroker stkIo(communicator);
      const std::string exodusFileName = "generated:1x1x8";
      size_t input_index = stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
      stkIo.set_active_mesh(input_index);
      stkIo.create_input_mesh();
      stkIo.populate_bulk_data();

      size_t idx = stkIo.create_output_mesh(file_name,
					    stk::io::WRITE_RESTART);

      stk::util::ParameterMapType::const_iterator i = params.begin();
      stk::util::ParameterMapType::const_iterator iend = params.end();
      for (; i != iend; ++i) {
	const std::string paramName = (*i).first;
	//+ NOTE: Need a reference to the parameter.
	stk::util::Parameter &param = params.get_param(paramName);
	//+ NOTE: Calling add_global_ref, passing address of value
	stkIo.add_global_ref(idx, paramName, &param.value, param.type);/*@\label{io:global:autoparam}*/
      }

      //+ All writing of the values is handled automatically,
      //+ do not need to call write_global
      stkIo.process_output_request(idx, 0.0);/*@\label{io:global:autowrite}*/
    }
    // ... Reading is the same as in previous example
    //-END

    // Read params from file...
    {
      stk::util::ParameterList gold_params; // To compare values read
      gold_params.set_param("PI", 3.14159);  // Double 
      gold_params.set_param("Answer", 42);   // Integer
      gold_params.set_param("some_doubles", my_vector); 
      gold_params.set_param("Ages", ages);   // Vector of integers...

      stk::io::StkMeshIoBroker stkIo(communicator);
      stkIo.add_mesh_database(file_name, stk::io::READ_MESH);
      stkIo.create_input_mesh();
      stkIo.populate_bulk_data();

      stkIo.read_defined_input_fields(0.0);

      size_t param_count = 0;
      stk::util::ParameterMapType::const_iterator i = params.begin();
      stk::util::ParameterMapType::const_iterator iend = params.end();
      for (; i != iend; ++i) {
	param_count++;
	const std::string paramName = (*i).first;
	stk::util::Parameter &param = params.get_param(paramName);
	stk::util::Parameter &gold_param
	  = gold_params.get_param(paramName);
	stkIo.get_global(paramName, param.value, param.type);
	validate_parameters_equal_value(param, gold_param);
      }

      std::vector<std::string> globalNamesOnFile;
      stkIo.get_global_variable_names(globalNamesOnFile);
      ASSERT_EQ(param_count, globalNamesOnFile.size());

    }
    unlink(file_name.c_str());
  }
}
