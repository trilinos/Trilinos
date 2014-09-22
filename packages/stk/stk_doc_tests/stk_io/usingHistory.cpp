// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <ostream>                      // for basic_ostream::operator<<
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker, etc
#include <stk_util/util/ParameterList.hpp>  // for ParameterList, etc
#include <stk_util/parallel/Parallel.hpp>
#include <string>                       // for string
#include <utility>                      // for pair
#include <vector>                       // for vector
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL
#include "Ioss_ElementTopology.h"       // for NameList
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_Property.h"              // for Property
#include "Ioss_Region.h"                // for Region
namespace Ioss { class DatabaseIO; }

namespace
{
  TEST(StkMeshIoBrokerHowTo, writeHistory)
  {

    const std::string file_name = "History.e";
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if (numProcs != 1) {
      return;
    }

    stk::util::ParameterList params;
    
    {
      // ========================================================================
      // INITIALIZATION...
      // Add some params to write and read...
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
    }

    {
      // ========================================================================
      // EXAMPLE USAGE...
      // Begin use of stk io history file...
      stk::io::StkMeshIoBroker stkIo(communicator);

      //-BEGIN
      //+ Define the heartbeat output and the format (BINARY) 
      size_t hb = stkIo.add_heartbeat_output(file_name, stk::io::BINARY);
      //-END

      stk::util::ParameterMapType::const_iterator i = params.begin();
      stk::util::ParameterMapType::const_iterator iend = params.end();
      for (; i != iend; ++i) {
	const std::string parameterName = (*i).first;
	stk::util::Parameter &parameter = params.get_param(parameterName);

	// Tell history database which global variables should be output at each step...
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
