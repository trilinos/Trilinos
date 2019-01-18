/*
//@HEADER
// ************************************************************************
// 
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions Contact  H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#include "mtk_kokkos.h"
#include <stk_ngp/Ngp.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/stk_config.h>
#include <stk_util/parallel/Parallel.hpp>
 
void test_field() {

  unsigned spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  unsigned dimX = 10;
  unsigned dimY = 10;
  unsigned dimZ = 10;
  std::ostringstream mesh_spec;
  mesh_spec << "generated:"<<dimX<<"x"<<dimY<<"x"<<dimZ;
  stk::io::fill_mesh(mesh_spec.str(), bulk);

  ngp::StaticMesh staticMesh(bulk);
  double initialValue = 9.9;
  ngp::StaticField<double> scalarField(stk::topology::NODE_RANK, initialValue, bulk, meta.locally_owned_part());

  stk::mesh::EntityVector nodes;
  stk::mesh::get_selected_entities(meta.locally_owned_part(), bulk.buckets(stk::topology::NODE_RANK), nodes);
  Kokkos::View<stk::mesh::Entity*> device_nodes("device_nodes", nodes.size());
  Kokkos::View<stk::mesh::Entity*,Kokkos::HostSpace> host_nodes("host_nodes", nodes.size());

  for(size_t i=0; i<nodes.size(); ++i) {
      host_nodes(i) = nodes[i];
  }
  Kokkos::deep_copy(device_nodes, host_nodes);

  struct timeval begin,end;

  gettimeofday(&begin,NULL);

  int nrepeat = 500;
  double result = 0;

  for(int n=0; n<nrepeat; ++n) {

      result = 0;
      Kokkos::parallel_reduce(nodes.size(), KOKKOS_LAMBDA(int i, double& update) {
        update += scalarField.get(staticMesh, device_nodes(i), 0);
      }, result);

  }

  gettimeofday(&end,NULL);

  Kokkos::fence();

  double time = 1.0*(end.tv_sec-begin.tv_sec) +
                1.0e-6*(end.tv_usec-begin.tv_usec);

  std::cout << "Time: " << time << std::endl;

  double goldValue = nodes.size()*initialValue;
  EXPECT_NEAR(goldValue, result, 0.001);
}

TEST_F(MTK_Kokkos, field_exp) {
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
      test_field();
  }
}

