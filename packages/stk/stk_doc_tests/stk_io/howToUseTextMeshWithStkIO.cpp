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

#include <gtest/gtest.h>
#include <stk_unit_test_utils/getOption.h>
#include <unistd.h>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>

namespace
{

TEST(StkIoHowTo, useTextMesh)
{
  stk::io::StkMeshIoBroker stkIo(MPI_COMM_WORLD);

  std::string textMeshDesc = "textmesh:0,1,HEX_8,1,2,3,4,5,6,7,8";

  stkIo.add_mesh_database(textMeshDesc, stk::io::READ_MESH);

  Ioss::Region* ioRegion = stkIo.get_input_ioss_region().get();
  EXPECT_TRUE(ioRegion != nullptr);

  Ioss::DatabaseIO* ioDatabase = ioRegion->get_database();
  EXPECT_TRUE(ioDatabase != nullptr);
  EXPECT_TRUE(ioDatabase->ok(true));
  EXPECT_EQ(ioDatabase->get_format(), "TextMesh");
}

TEST(StkIoHowTo, useTextMesh_withAllOptions)
{
  stk::io::StkMeshIoBroker stkIo(MPI_COMM_WORLD);

  std::string textMeshDesc =
      "textmesh:"
      "0,1,QUAD_4_2D,1,2,5,4\n"
      "0,2,QUAD_4_2D,2,3,6,5"
      "|coordinates:0,0, 1,0, 2,0, 0,1, 1,1, 2,1"
      "|dimension:2";

  stkIo.add_mesh_database(textMeshDesc, stk::io::READ_MESH);

  Ioss::Region* ioRegion = stkIo.get_input_ioss_region().get();
  EXPECT_TRUE(ioRegion != nullptr);

  Ioss::DatabaseIO* ioDatabase = ioRegion->get_database();
  EXPECT_TRUE(ioDatabase != nullptr);
  EXPECT_TRUE(ioDatabase->ok(true));
  EXPECT_EQ(ioDatabase->get_format(), "TextMesh");
}

}  // namespace
