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

#include "stk_io/StkIoUtils.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_io/FillMesh.hpp"
#include "stk_io/IossBridge.hpp"
#include <stk_util/parallel/Parallel.hpp>
#include <gtest/gtest.h>
#include <stk_unit_test_utils/BuildMesh.hpp>

namespace {

using stk::unit_test_util::build_mesh;
using stk::unit_test_util::build_mesh_no_simple_fields;

TEST(ParallelFileName, serialNameShouldBeUnchanged)
{
  std::string fileName = "filename.exo";
  EXPECT_EQ(stk::io::construct_filename_for_serial_or_parallel(fileName, 1, 0), "filename.exo");
}

TEST(ParallelFileName, singleDigitNames)
{
  std::string fileName = "filename.exo";
  EXPECT_EQ(stk::io::construct_filename_for_serial_or_parallel(fileName, 5, 0), "filename.exo.5.0");
  EXPECT_EQ(stk::io::construct_filename_for_serial_or_parallel(fileName, 9, 0), "filename.exo.9.0");
  EXPECT_EQ(stk::io::construct_filename_for_serial_or_parallel(fileName, 9, 8), "filename.exo.9.8");
}

TEST(ParallelFileName, multiDigitNames)
{
  std::string fileName = "filename.exo";
  EXPECT_EQ(stk::io::construct_filename_for_serial_or_parallel(fileName, 10, 0),  "filename.exo.10.00");
  EXPECT_EQ(stk::io::construct_filename_for_serial_or_parallel(fileName, 10, 9),  "filename.exo.10.09");
  EXPECT_EQ(stk::io::construct_filename_for_serial_or_parallel(fileName, 99, 0),  "filename.exo.99.00");
  EXPECT_EQ(stk::io::construct_filename_for_serial_or_parallel(fileName, 99, 9),  "filename.exo.99.09");
  EXPECT_EQ(stk::io::construct_filename_for_serial_or_parallel(fileName, 99, 10), "filename.exo.99.10");
  EXPECT_EQ(stk::io::construct_filename_for_serial_or_parallel(fileName, 99, 98), "filename.exo.99.98");
  EXPECT_EQ(stk::io::construct_filename_for_serial_or_parallel(fileName, 100, 98), "filename.exo.100.098");
}

TEST(ParallelFileName, InvalidInputs) {
  std::string fileName = "filename.exo";
  EXPECT_THROW(stk::io::construct_filename_for_serial_or_parallel(fileName,  -1,  0),  std::exception);
  EXPECT_THROW(stk::io::construct_filename_for_serial_or_parallel(fileName, 100, -5),  std::exception);
  EXPECT_THROW(stk::io::construct_filename_for_serial_or_parallel(fileName, 107, 107), std::exception);
  EXPECT_THROW(stk::io::construct_filename_for_serial_or_parallel(fileName, 107, 130), std::exception);

}

TEST(CheckElemBlockTopology, validTopology)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD);
  stk::io::fill_mesh("generated:2x2x2", *bulk);
  EXPECT_NO_THROW(stk::io::throw_if_any_elem_block_has_invalid_topology(bulk->mesh_meta_data(), "test"));
}

TEST(CheckElemBlockTopology, invalidTopology)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::MetaData meta(3);
  stk::mesh::Part& block2 = meta.declare_part("block_2", stk::topology::ELEM_RANK);
  stk::io::put_io_part_attribute(block2);
  EXPECT_THROW(stk::io::throw_if_any_elem_block_has_invalid_topology(meta, "test"), std::runtime_error);
}

}
