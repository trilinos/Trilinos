#include <gtest/gtest.h>                // for AssertHelper, EXPECT_STREQ, etc
#include <mpi.h>                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <istream>                      // for basic_istream, ifstream
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker, etc
#include <stk_util/util/ParameterList.hpp>  // for ParameterList, etc
#include <stk_util/parallel/Parallel.hpp>
#include <string>                       // for string, getline
#include <utility>                      // for pair
#include <vector>                       // for vector
#include "Ioss_Property.h"              // for Property
#include "Ioss_PropertyManager.h"       // for PropertyManager

namespace
{
  TEST(StkMeshIoBrokerHowTo, writeHeartbeatOverrideSeparator)
  {

    const std::string file_name = "HeartbeatBarSep.txt";
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);

    if (numProcs != 1) {
      return;
    }

    stk::util::ParameterList parameters;
    
    {
      // ========================================================================
      // INITIALIZATION...
      // Add some parameters to write and read...
      parameters.set_param("PI", -3.14159);  // Double 
      parameters.set_param("Answer", 42);   // Integer

      std::vector<double> my_vector;
      my_vector.push_back(2.78);
      my_vector.push_back(5.30);
      my_vector.push_back(6.21);
      parameters.set_param("dbl_vec", my_vector);   // Vector of doubles...
    
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
      // Begin use of stk io heartbeat file...
      stk::io::StkMeshIoBroker stkIo(communicator);

      //-BEGIN
      //+ Use vertical bar as field separator 
      Ioss::PropertyManager hb_props;
      hb_props.add(Ioss::Property("FIELD_SEPARATOR", " | "));
      size_t hb =
	stkIo.add_heartbeat_output(file_name, stk::io::TEXT, hb_props);
      //-END

      stk::util::ParameterMapType::const_iterator i = parameters.begin();
      stk::util::ParameterMapType::const_iterator iend = parameters.end();
      for (; i != iend; ++i) {
	const std::string parameterName = (*i).first;
	stk::util::Parameter &parameter = parameters.get_param(parameterName);

	// Tell heartbeat database which global variables should be output at each step...
	stkIo.add_heartbeat_global(hb, parameterName, &parameter.value, parameter.type);
      }

      // Now output the global variables...
      int timestep_count = 1;
      double time = 0.0;
      for (int step=1; step <= timestep_count; step++) {
	stkIo.process_heartbeat_output(hb, step, time);
      }
    }

    {
      // ========================================================================
      // VERIFICATION:
      // open the heartbeat file...
      std::ifstream heartbeat(file_name.c_str());
      std::string header_line;
      std::string data_line;

      std::string expected_header_line = "        Time |       Ages_1 |       Ages_2 |       Ages_3 |       Ages_4 |       Answer |           PI |    dbl_vec_1 |    dbl_vec_2 |    dbl_vec_3";
      std::string expected_data_line = " 0.00000e+00 |           55 |           49 |           21 |           19 |           42 | -3.14159e+00 |  2.78000e+00 |  5.30000e+00 |  6.21000e+00";

      EXPECT_TRUE(!std::getline(heartbeat, header_line).fail());
      EXPECT_STREQ(header_line.c_str(), expected_header_line.c_str());
      EXPECT_TRUE(!std::getline(heartbeat, data_line).fail());
      EXPECT_STREQ(data_line.c_str(), expected_data_line.c_str());
    }

    // ========================================================================
    // CLEANUP:
    unlink(file_name.c_str());
  }
}
