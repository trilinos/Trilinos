#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <mpi.h>
#include <stk_util/util/ParameterList.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <fieldNameTestUtils.hpp>

namespace
{
  TEST(StkMeshIoBrokerHowTo, writeAndReadGlobalParameters)
  {
    // ============================================================
    //+ INITIALIZATION
    const std::string file_name = "GlobalParameters.e";
    MPI_Comm communicator = MPI_COMM_WORLD;

    stk::util::ParameterList params;
    stk::util::ParameterList gold_params; // To compare values read
    
    // Add some parameters to write and read...
    params.set_param("PI", 3.14159);  // Double 
    params.set_param("Answer", 42);   // Integer
    gold_params.set_param("PI", 3.14159);  // Double 
    gold_params.set_param("Answer", 42);   // Integer

    std::vector<double> my_vector;
    my_vector.push_back(2.78);
    my_vector.push_back(5.30);
    my_vector.push_back(6.21);
    params.set_param("doubles", my_vector); // Vector of doubles...
    gold_params.set_param("doubles", my_vector);
    
    std::vector<int> ages;
    ages.push_back(55);
    ages.push_back(49);
    ages.push_back(21);
    ages.push_back(19);
    
    params.set_param("Ages", ages);   // Vector of integers...
    gold_params.set_param("Ages", ages);   // Vector of integers...
    

    {
      stk::io::StkMeshIoBroker stkIo(communicator);
      generateMetaData(stkIo);
      stkIo.populate_bulk_data();

      // ============================================================
      //+ EXAMPLE
      //+ Write output file with all parameters in params list...
      size_t idx = stkIo.create_output_mesh(file_name,
					    stk::io::WRITE_RESTART);

      stk::util::ParameterMapType::const_iterator i  = params.begin();
      stk::util::ParameterMapType::const_iterator ie = params.end();
      for (; i != ie; ++i) {
	const std::string parameterName = (*i).first;
	stk::util::Parameter &param = params.get_param(parameterName);
	stkIo.add_global(idx, parameterName, param.value, param.type);
      }

      stkIo.begin_output_step(idx, 0.0);/*@\label{io:global:write_begin}*/

      for (i = params.begin(); i != ie; ++i) {
	const std::string parameterName = (*i).first;
	stk::util::Parameter &param = params.get_param(parameterName);
	stkIo.write_global(idx, parameterName, param.value, param.type);
      }

      stkIo.end_output_step(idx);/*@\label{io:global:write_end}*/
    }

    {
      // ============================================================
      //+ EXAMPLE
      //+ Read parameters from file...
      stk::io::StkMeshIoBroker stkIo(communicator);
      stkIo.open_mesh_database(file_name, stk::io::READ_MESH);
      stkIo.create_input_mesh();
      stkIo.populate_bulk_data();

      stkIo.read_defined_input_fields(0.0);

      stk::util::ParameterMapType::const_iterator i = params.begin();
      stk::util::ParameterMapType::const_iterator ie = params.end();
      for (; i != ie; ++i) {
	const std::string parameterName = (*i).first;
	stk::util::Parameter &param = params.get_param(parameterName);
	stkIo.get_global(parameterName, param.value, param.type);
      }

      // ============================================================
      //+ VALIDATION
      std::vector<std::string> globalNamesOnFile;
      stkIo.get_global_variable_names(globalNamesOnFile);
      ASSERT_EQ(param_count, globalNamesOnFile.size());

      size_t param_count = 0;
      for (i = params.begin(); i != ie; ++i) {
	param_count++;
	const std::string parameterName = (*i).first;
	stk::util::Parameter &param = params.get_param(parameterName);
	stk::util::Parameter &gold_parameter =
	  gold_params.get_param(parameterName);
	validate_parameters_equal_value(param, gold_parameter);
      }

      std::vector<std::string> globalNamesOnFile;
      stkIo.get_global_variable_names(globalNamesOnFile);
      ASSERT_EQ(param_count, globalNamesOnFile.size());

    }
    // ============================================================
    // CLEAN UP
    unlink(file_name.c_str());
  }
}
