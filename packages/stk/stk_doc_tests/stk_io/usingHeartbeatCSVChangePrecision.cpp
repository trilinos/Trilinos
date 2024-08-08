// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_STREQ, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <istream>                      // for basic_istream, ifstream
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker, etc
#include <stk_util/util/ParameterList.hpp>  // for ParameterList, etc
#include <string>                       // for string, getline
#include <utility>                      // for pair
#include <vector>                       // for vector
#include "Ioss_Property.h"              // for Property
#include "Ioss_PropertyManager.h"       // for PropertyManager
#include "stk_util/parallel/Parallel.hpp"

namespace
{
TEST(StkMeshIoBrokerHowTo, writeHeartbeatCSVChangePrecision)
{

  const std::string file_name = "Heartbeat-CSV.txt";
  MPI_Comm communicator = MPI_COMM_WORLD;
  int num_procs = stk::parallel_machine_size(communicator);
  if (num_procs != 1) {
    return;
  }
  int my_processor = stk::parallel_machine_rank(communicator);

  stk::util::ParameterList params;

  {
    // ============================================================
    // INITIALIZATION...
    // Add some params to write and read...
    params.set_param("PI", -3.14159);  // Double
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
  }

  {
    // ============================================================
    // EXAMPLE USAGE...
    // Define the heartbeat output...
    stk::io::StkMeshIoBroker stkIo(communicator);

    //-BEGIN
    //+ Output should have 10 digits of precision (1.0123456789e+00)
    //+             default precision is 5 digits (1.012345e+00)
    Ioss::PropertyManager hb_props;
    hb_props.add(Ioss::Property("PRECISION", 10));

    //+ Define the heartbeat output and the format (CSV)
    size_t hb = stkIo.add_heartbeat_output(file_name, stk::io::CSV, hb_props);
    //-END

    stk::util::ParameterMapType::const_iterator i = params.begin();
    stk::util::ParameterMapType::const_iterator iend = params.end();
    for (; i != iend; ++i) {
      const std::string parameterName = (*i).first;
      stk::util::Parameter &parameter = params.get_param(parameterName);

      //+ Tell heartbeat which variables to output at each step...
      stkIo.add_heartbeat_global(hb, parameterName, parameter);
    }

    int timestep_count = 1;
    double time = 0.0;
    for (int step=1; step <= timestep_count; step++) {
      //+ Now output the global variables...
      stkIo.process_heartbeat_output(hb, step, time);
    }
  }

  if (my_processor == 0) { // Heartbeat is only output on processor 0.
    // ============================================================
    // VERIFICATION:
    // Expected output:
    std::string expected_header_line =
        "             Time,            Ages_1,            Ages_2,            Ages_3,            Ages_4,            Answer,                PI,    some_doubles_1,    some_doubles_2,    some_doubles_3";
    std::string expected_data_line =
        " 0.0000000000e+00,                55,                49,                21,                19,                42, -3.1415900000e+00,  2.7800000000e+00,  5.3000000000e+00,  6.2100000000e+00";

    // open the heartbeat file...
    std::ifstream heartbeat(file_name.c_str());
    std::string header_line;
    std::string data_line;

    EXPECT_TRUE(!std::getline(heartbeat, header_line).fail());
    EXPECT_STREQ(header_line.c_str(), expected_header_line.c_str());
    EXPECT_TRUE(!std::getline(heartbeat, data_line).fail());
    EXPECT_STREQ(data_line.c_str(), expected_data_line.c_str());

    // ============================================================
    // CLEANUP:
    unlink(file_name.c_str());
  }
}
}
