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

#include <gtest/gtest.h>                // for AssertHelper, ASSERT_EQ, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <exception>                    // for exception
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <string>                       // for string, basic_string
#include <vector>                       // for vector
#include <limits>                       // for std::numeric_limits
#include "Ioss_Field.h"                 // for Field, etc
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_util/parallel/Parallel.hpp"

namespace
{
//-BEGIN
TEST(StkMeshIoBrokerHowTo, writeAndReadGlobalVariables)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 1) { return; }

  const std::string restartFileName = "OneGlobalDouble.restart";
  const std::string timeStepVarName = "timeStep";
  const double timeStepSize = 1e-6;
  const double currentTime = 1.0;

  //+ Write restart file with time step size as a global variable
  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    size_t fileIndex =
        stkIo.create_output_mesh(restartFileName, stk::io::WRITE_RESTART);
    stkIo.add_global(fileIndex, timeStepVarName, Ioss::Field::REAL);
    stkIo.begin_output_step(fileIndex, currentTime);
    stkIo.write_global(fileIndex, timeStepVarName, timeStepSize);
    stkIo.end_output_step(fileIndex);
  }

  //+ Read restart file with time step size as a global variable
  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.add_mesh_database(restartFileName, stk::io::READ_RESTART);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();
    stkIo.read_defined_input_fields(currentTime);
    std::vector<std::string> globalNamesOnFile;
    stkIo.get_global_variable_names(globalNamesOnFile);

    ASSERT_EQ(1u, globalNamesOnFile.size());
    EXPECT_STRCASEEQ(timeStepVarName.c_str(),
                     globalNamesOnFile[0].c_str());
    double timeStepSizeReadFromFile = 0.0;
    stkIo.get_global(globalNamesOnFile[0], timeStepSizeReadFromFile);

    const double epsilon = std::numeric_limits<double>::epsilon();
    EXPECT_NEAR(timeStepSize, timeStepSizeReadFromFile, epsilon);

    //+ If try to get a global that does not exist, will throw
    //+ an exception by default...
    double value = 0.0;
    EXPECT_ANY_THROW(stkIo.get_global("does_not_exist", value));

    //+ If the application wants to handle the error instead (without a try/catch),
    //+ can pass in an optional boolean:
    bool abort_if_not_found = false;
    bool found = stkIo.get_global("does_not_exist", value, abort_if_not_found);
    ASSERT_FALSE(found);
  }

  unlink(restartFileName.c_str());
}
//-END
}
