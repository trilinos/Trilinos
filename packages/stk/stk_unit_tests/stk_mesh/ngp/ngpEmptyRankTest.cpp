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
#include <stk_util/stk_config.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <string>
#include <memory>
#include <cstdlib>

namespace {

void write_single_element_mesh(const std::string& fileName)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_rank(comm) == 0) {
    MPI_Comm commSelf = stk::parallel_machine_self();
    stk::mesh::MeshBuilder builder(commSelf);
    std::shared_ptr<stk::mesh::BulkData> bulk =
         builder.set_spatial_dimension(3).create();
    stk::io::fill_mesh("generated:1x1x1",*bulk);
    stk::io::write_mesh(fileName, *bulk);
  }
}

TEST(NgpEmptyRank, get_updated_ngp_mesh)
{
//This test exercises a case that was the subject of a user bug-report
//After reading a single-element mesh on 2 MPI procs, and with
//inconsistent setting (or no setting) of spatial-dimension on the
//mesh-builder, the call to get_updated_ngp_mesh was seg-faulting.
//This issue has now been fixed.
//
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }
  const int myProc = stk::parallel_machine_rank(comm);

  std::string fileName("hex8.exo");
  write_single_element_mesh(fileName);
  stk::parallel_machine_barrier(comm);

  stk::mesh::MeshBuilder builder(comm);
  std::shared_ptr<stk::mesh::BulkData> bulk;
  if (myProc == 0) {
    //Inconsistent setting of spatial-dim (on one proc but not the other)
    //is not really valid. But we don't prevent it, and we can tolerate it.
    //(It was included in the user bug-report.)
    bulk = builder.set_spatial_dimension(3).create();
  }
  else {
    bulk = builder.create();
  }

  stk::io::fill_mesh_with_auto_decomp(fileName,*bulk);

  //get_updated_ngp_mesh should neither throw nor seg-fault.
  //(user bug-report was seg-fault beneath the internal call to
  //DeviceMesh::fill_buckets)
  EXPECT_NO_THROW(stk::mesh::get_updated_ngp_mesh(*bulk));
}

}

