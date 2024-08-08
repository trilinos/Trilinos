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

#include <gtest/gtest.h>                // for TEST
#include <sstream>                      // for ostringstream
#include <unistd.h>
#include "mpi.h"                        // for MPI_COMM_WORLD

#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/DumpMeshInfo.hpp"
#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"  // for HexFixture

namespace {

TEST(UnitTestDebugDump, dump_all_meta_info)
{
  const unsigned NX = 3;
  const unsigned NY = 1;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, NX, NY, NZ);
  fixture.m_meta.commit();

  // Doesn't check anything, but at least makes sure it builds and runs
  std::ostringstream oss;
  stk::mesh::impl::dump_all_meta_info(fixture.m_meta, oss);

//  std::cout << oss.str() << std::endl;
}

TEST(UnitTestDebugDump, dump_all_mesh_info)
{
  const unsigned NX = 3;
  const unsigned NY = 1;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD, NX, NY, NZ);
  hf.m_meta.commit();
  hf.generate_mesh();

  // Doesn't check anything, but at least makes sure it builds and runs
  std::ostringstream oss;
  stk::mesh::impl::dump_all_mesh_info(hf.m_bulk_data, oss);

//  std::cout << oss.str() << std::endl;
}

TEST(UnitTestDebugDump, dump_mesh_per_proc)
{
  const unsigned NX = 3;
  const unsigned NY = 1;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD, NX, NY, NZ);
  hf.m_meta.commit();
  hf.generate_mesh();

  // Doesn't check anything, but at least makes sure it builds and runs
  std::ostringstream oss;
  stk::mesh::impl::dump_mesh_per_proc(hf.m_bulk_data, "junk_mesh");

  unlink(("junk_mesh." + std::to_string(stk::parallel_machine_rank(MPI_COMM_WORLD))).c_str());
}

TEST(UnitTestDebugDump, dump_partition_summary)
{
  const unsigned NX = 3;
  const unsigned NY = 1;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD, NX, NY, NZ);
  hf.m_meta.commit();
  hf.generate_mesh();

  // Doesn't check anything, but at least makes sure it builds and runs
  std::ostringstream oss;
  stk::mesh::impl::dump_partition_summary(hf.m_bulk_data, oss);

//  std::cout << oss.str() << std::endl;
}

TEST(UnitTestDebugDump, dump_partition_summary_per_proc)
{
  const unsigned NX = 3;
  const unsigned NY = 1;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD, NX, NY, NZ);
  hf.m_meta.commit();
  hf.generate_mesh();

  // Doesn't check anything, but at least makes sure it builds and runs
  std::ostringstream oss;
  stk::mesh::impl::dump_partition_summary_per_proc(hf.m_bulk_data, "junk_partition_summary");

  unlink(("junk_partition_summary." + std::to_string(stk::parallel_machine_rank(MPI_COMM_WORLD))).c_str());
}

TEST(UnitTestDebugDump, dump_bucket_size_histogram)
{
  const unsigned NX = 3;
  const unsigned NY = 1;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD, NX, NY, NZ);
  hf.m_meta.commit();
  hf.generate_mesh();

  // Doesn't check anything, but at least makes sure it builds and runs
  std::ostringstream oss;
  stk::mesh::impl::dump_bucket_size_histogram(hf.m_bulk_data, oss);

//  std::cout << oss.str() << std::endl;
}

TEST(UnitTestDebugDump, dump_bucket_size_histogram_per_proc)
{
  const unsigned NX = 3;
  const unsigned NY = 1;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD, NX, NY, NZ);
  hf.m_meta.commit();
  hf.generate_mesh();

  // Doesn't check anything, but at least makes sure it builds and runs
  std::ostringstream oss;
  stk::mesh::impl::dump_bucket_size_histogram_per_proc(hf.m_bulk_data, "junk_histogram");

  unlink(("junk_histogram." + std::to_string(stk::parallel_machine_rank(MPI_COMM_WORLD))).c_str());
}

TEST(UnitTestDebugDump, dump_global_bucket_size_histogram)
{
  const unsigned NX = 3;
  const unsigned NY = 1;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD, NX, NY, NZ);
  hf.m_meta.commit();
  hf.generate_mesh();

  // Doesn't check anything, but at least makes sure it builds and runs
  std::ostringstream oss;
  stk::mesh::impl::dump_bucket_size_histogram_per_proc(hf.m_bulk_data, "junk_histogram");

  unlink(("junk_histogram." + std::to_string(stk::parallel_machine_rank(MPI_COMM_WORLD))).c_str());
}

}
