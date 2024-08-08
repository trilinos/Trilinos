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

#include <gtest/gtest.h>                // for ASSERT_EQ, AssertHelper, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include "parameterTestUtils.hpp"
#include <stk_util/util/ParameterList.hpp>  // for ParameterList, etc
#include <string>                       // for string
#include <utility>                      // for pair
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc

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
      stkIo.add_global_ref(idx, paramName, param);/*@\label{io:global:autoparam}*/
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
      stkIo.get_global(paramName, param);
      validate_parameters_equal_value(param, gold_param);
    }

    std::vector<std::string> globalNamesOnFile;
    stkIo.get_global_variable_names(globalNamesOnFile);
    ASSERT_EQ(param_count, globalNamesOnFile.size());

  }
  unlink(file_name.c_str());
}
}
