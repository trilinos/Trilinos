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
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/NgpFieldBLAS.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_unit_test_utils/timer.hpp>

namespace
{

TEST(StkNgpField, multiStateRotation)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 3000;
  std::string meshSpec = "generated:80x80x80";

  std::cout << "Using mesh-spec: " << meshSpec << std::endl;

  stk::unit_test_util::BatchTimer batchTimer(comm);

  batchTimer.initialize_batch_timer();

  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm)
                                                      .set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA)
                                                      .set_spatial_dimension(3)
                                                      .create();

  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  const int numFieldStates = 3;
  stk::mesh::Field<double>& tensorField1 = meta.declare_field<double>(stk::topology::ELEM_RANK, "tensorField1", numFieldStates);
  stk::mesh::Field<double>& tensorField2 = meta.declare_field<double>(stk::topology::ELEM_RANK, "tensorField2", numFieldStates);
  stk::mesh::Field<double>& vectorField1 = meta.declare_field<double>(stk::topology::ELEM_RANK, "vectorField1", numFieldStates);
  stk::mesh::Field<double>& vectorField2 = meta.declare_field<double>(stk::topology::ELEM_RANK, "vectorField2", numFieldStates);
  stk::mesh::put_field_on_mesh(tensorField1, meta.universal_part(), 9, nullptr);
  stk::mesh::put_field_on_mesh(tensorField2, meta.universal_part(), 9, nullptr);
  stk::mesh::put_field_on_mesh(vectorField1, meta.universal_part(), 3, nullptr);
  stk::mesh::put_field_on_mesh(vectorField2, meta.universal_part(), 3, nullptr);

  stk::io::fill_mesh(meshSpec, *bulkPtr);

  Kokkos::Profiling::pushRegion("get_updated_ngp_mesh");
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);
  EXPECT_FALSE(ngpMesh.need_sync_to_host());
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("initialize fields");
  stk::ngp::ExecSpace execSpace;
  constexpr double initValue1 = 1.14;
  constexpr double initValue2 = 3.14;
  for(int s=0; s<numFieldStates; ++s) {
    const stk::mesh::FieldState state = static_cast<stk::mesh::FieldState>(s);
    stk::mesh::Field<double>& tensorField1_state = tensorField1.field_of_state(state);
    stk::mesh::Field<double>& tensorField2_state = tensorField2.field_of_state(state);
    stk::mesh::Field<double>& vectorField1_state = vectorField1.field_of_state(state);
    stk::mesh::Field<double>& vectorField2_state = vectorField2.field_of_state(state);
    stk::mesh::field_fill(initValue1, tensorField1_state, execSpace);
    stk::mesh::field_fill(initValue2, tensorField2_state, execSpace);
    stk::mesh::field_fill(initValue1, vectorField1_state, execSpace);
    stk::mesh::field_fill(initValue2, vectorField2_state, execSpace);
  }
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("multiStateRotation test");

  for (unsigned j = 0; j < NUM_RUNS; j++) {

    batchTimer.start_batch_timer();

    for(unsigned i=0; i<NUM_ITERS; ++i) {
      Kokkos::Profiling::pushRegion("field_copy");
      if (i%2==0) {
        stk::mesh::field_copy(tensorField1, tensorField2, execSpace);
        stk::mesh::field_copy(vectorField1, vectorField2, execSpace);
      }
      else {
        stk::mesh::field_copy(tensorField2, tensorField1, execSpace);
        stk::mesh::field_copy(vectorField2, vectorField1, execSpace);
      }
      Kokkos::Profiling::popRegion();

      Kokkos::Profiling::pushRegion("update_field_data_states");
      const bool rotateNgpFieldViews = true;
      bulkPtr->update_field_data_states(rotateNgpFieldViews);
      Kokkos::Profiling::popRegion();
    }

    batchTimer.stop_batch_timer();
  }

  Kokkos::Profiling::popRegion();
  batchTimer.print_batch_timing(NUM_ITERS);
}

}
