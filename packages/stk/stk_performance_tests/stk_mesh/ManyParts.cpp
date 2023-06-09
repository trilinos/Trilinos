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

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_unit_test_utils/timer.hpp>

TEST(many_parts, many_parts)
{
  const unsigned NUM_RUNS = 5;
  unsigned num_parts = 8000;
  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);

    stk::mesh::Part& super1 = meta.declare_part("super1");
    stk::mesh::Part& super2 = meta.declare_part("super2");
    stk::mesh::Part& super3 = meta.declare_part("super3");
    stk::mesh::Part& super4 = meta.declare_part("super4");

    stk::mesh::Part& sub1 = meta.declare_part("sub1");
    stk::mesh::Part& sub2 = meta.declare_part("sub2");
    stk::mesh::Part& sub3 = meta.declare_part("sub3");
    stk::mesh::Part& sub4 = meta.declare_part("sub4");

    stk::mesh::PartVector parts;
    for(unsigned i=0; i<num_parts; ++i) {
      std::string partName = "part_" + std::to_string(i);
      stk::mesh::Part& part = meta.declare_part(partName);
      meta.declare_part_subset(super1, part);
      meta.declare_part_subset(super2, part);
      meta.declare_part_subset(super3, part);
      meta.declare_part_subset(part, sub1);
      meta.declare_part_subset(part, sub2);
      meta.declare_part_subset(part, sub3);
      parts.push_back(&part);
    }

    typedef stk::mesh::Field<double> VectorField;
    VectorField& field = meta.declare_field<double>(stk::topology::NODE_RANK, "field");
    for(size_t i=0; i<parts.size(); ++i) {
      const stk::mesh::Part& part = *parts[i];
      stk::mesh::put_field_on_mesh(field, part, 3, nullptr);
    }

    for(size_t i=0; i<parts.size(); ++i) {
      meta.declare_part_subset(super4, *parts[i]);
      meta.declare_part_subset(*parts[i], sub4);
    }

    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(num_parts);
}
