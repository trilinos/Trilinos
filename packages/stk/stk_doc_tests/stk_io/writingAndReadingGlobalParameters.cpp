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
//-BEGIN
TEST(StkMeshIoBrokerHowTo, writeAndReadGlobalParameters)
{
  // ============================================================
  //+ INITIALIZATION
  const std::string file_name = "GlobalParameters.e";
  MPI_Comm communicator = MPI_COMM_WORLD;

  // Add some parameters to write and read...
  stk::util::ParameterList params;
  params.set_param("PI", 3.14159);   // Double
  params.set_param("Answer", 42);    // Integer

  std::vector<double> my_vector = { 2.78, 5.30, 6.21 };
  params.set_param("doubles", my_vector); // Vector of doubles...

  std::vector<int> ages = { 55, 49, 21, 19};
  params.set_param("Ages", ages);   // Vector of integers...

  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    size_t index = stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();
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
      stkIo.add_global(idx, parameterName, param);
    }

    stkIo.begin_output_step(idx, 0.0);/*@\label{io:global:write_begin}*/

    for (i = params.begin(); i != ie; ++i) {
      const std::string parameterName = (*i).first;
      stk::util::Parameter &param = params.get_param(parameterName);
      stkIo.write_global(idx, parameterName, param);
    }

    stkIo.end_output_step(idx);/*@\label{io:global:write_end}*/
  }

  {
    // ============================================================
    //+ EXAMPLE
    //+ Read parameters from file...
    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.add_mesh_database(file_name, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    stkIo.read_defined_input_fields(0.0);

    stk::util::ParameterMapType::const_iterator i = params.begin();
    stk::util::ParameterMapType::const_iterator ie = params.end();
    for (; i != ie; ++i) {
      const std::string parameterName = (*i).first;
      stk::util::Parameter &param = params.get_param(parameterName);
      stkIo.get_global(parameterName, param);
    }

    // ============================================================
    //+ VALIDATION
    stk::util::ParameterList gold_params; // To compare values read
    gold_params.set_param("PI", 3.14159);        // Double
    gold_params.set_param("Answer", 42);         // Integer
    gold_params.set_param("doubles", my_vector); // Vector of doubles
    gold_params.set_param("Ages", ages);         // Vector of integers...

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
//-END
}
