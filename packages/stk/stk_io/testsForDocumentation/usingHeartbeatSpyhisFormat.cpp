#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <mpi.h>
#include <stk_util/util/ParameterList.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <fieldNameTestUtils.hpp>

namespace
{
  TEST(StkMeshIoBrokerHowTo, writeHeartbeatSpyhisFormat)
  {

    const std::string file_name = "Heartbeat-Spyhis.txt";
    MPI_Comm communicator = MPI_COMM_WORLD;

    stk::util::ParameterList parameters;
    
    // Initialization...
    {
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

    // Begin use of stk io heartbeat file...
    stk::io::StkMeshIoBroker stkIo(communicator);

    // ========================================================================
    // Define the heartbeat output.
    Ioss::PropertyManager hb_props;
    hb_props.add(Ioss::Property("FILE_FORMAT", "spyhis"));

    size_t heartbeat_index = stkIo.add_heartbeat_output(file_name, stk::io::TEXT, hb_props);

    stk::util::ParameterMapType::const_iterator i = parameters.begin();
    stk::util::ParameterMapType::const_iterator iend = parameters.end();
    for (; i != iend; ++i) {
      const std::string parameterName = (*i).first;
      stk::util::Parameter &parameter = parameters.get_param(parameterName);

      // Tell heartbeat database which global variables should be output at each step...
      stkIo.add_heartbeat_global(heartbeat_index, parameterName, parameter.value, parameter.type);
    }

    // ========================================================================
    // Now output the global variables...
    int timestep_count = 1;
    double time = 0.0;
    for (int step=1; step <= timestep_count; step++) {
      stkIo.process_heartbeat_output(heartbeat_index, step, time);
    }

    //    unlink(file_name.c_str());
  }
}
