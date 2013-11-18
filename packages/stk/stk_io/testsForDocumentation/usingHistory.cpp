#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <mpi.h>
#include <stk_util/util/ParameterList.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <fieldNameTestUtils.hpp>

namespace
{
  TEST(StkMeshIoBrokerHowTo, writeHistory)
  {

    const std::string file_name = "History.e";
    MPI_Comm communicator = MPI_COMM_WORLD;

    stk::util::ParameterList parameters;
    
    {
      // ========================================================================
      // INITIALIZATION...
      // Add some parameters to write and read...
      parameters.set_param("PI", 3.14159);  // Double 
      parameters.set_param("Answer", 42);   // Integer

      std::vector<double> my_vector;
      my_vector.push_back(2.78);
      my_vector.push_back(5.30);
      my_vector.push_back(6.21);
      parameters.set_param("some_doubles", my_vector);   // Vector of doubles...
    
      std::vector<int> ages;
      ages.push_back(55);
      ages.push_back(49);
      ages.push_back(21);
      ages.push_back(19);
    
      parameters.set_param("Ages", ages);   // Vector of integers...
    }

    {
      // ========================================================================
      // EXAMPLE USAGE...
      // Begin use of stk io history file...
      stk::io::StkMeshIoBroker stkIo(communicator);

      // Define the history output.
      size_t history_index = stkIo.add_heartbeat_output(file_name, stk::io::BINARY);

      stk::util::ParameterMapType::const_iterator i = parameters.begin();
      stk::util::ParameterMapType::const_iterator iend = parameters.end();
      for (; i != iend; ++i) {
	const std::string parameterName = (*i).first;
	stk::util::Parameter &parameter = parameters.get_param(parameterName);

	// Tell history database which global variables should be output at each step...
	stkIo.add_heartbeat_global(history_index, parameterName, parameter.value, parameter.type);
      }

      // Now output the global variables...
      int timestep_count = 1;
      double time = 0.0;
      for (int step=1; step <= timestep_count; step++) {
	stkIo.process_heartbeat_output(history_index, step, time);
      }
    }

    {
      // ========================================================================
      // VERIFICATION:
      Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", file_name, Ioss::READ_MODEL, communicator);
      Ioss::Region region(iossDb);

      EXPECT_EQ(region.get_property("state_count").get_int(), 1);
      region.begin_state(1);
      
      Ioss::NameList fields;
      region.field_describe(Ioss::Field::TRANSIENT, &fields);
      EXPECT_EQ(fields.size(), 4u);
      
      std::vector<double> values;
      region.get_field_data("PI", values);
      EXPECT_NEAR(values[0], 3.14159, 1.0e-6);

      std::vector<int> ages;
      region.get_field_data("Ages", ages);
      EXPECT_EQ(ages.size(), 4u);
      EXPECT_EQ(ages[0], 55);
      EXPECT_EQ(ages[1], 49);
      EXPECT_EQ(ages[2], 21);
      EXPECT_EQ(ages[3], 19);
    }

    // ========================================================================
    // CLEANUP:
    unlink(file_name.c_str());
  }
}
