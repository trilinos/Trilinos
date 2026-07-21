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
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_expreval/Eval.hpp>
#include <stk_tools/mesh_tools/MeshEvalFunctor.hpp>
#include <stk_io/FillMesh.hpp>
#include <iostream>
#include <memory>

namespace {

TEST(ExprEval, stkMesh)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  auto meshPtr = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh("generated:2x2x2", *meshPtr);

  stk::mesh::MetaData& meta = meshPtr->mesh_meta_data();
  unsigned numElements = stk::mesh::count_entities(*meshPtr, stk::topology::ELEM_RANK, meta.universal_part());
  EXPECT_EQ(8u, numElements);

  auto eval = std::make_shared<stk::expreval::Eval>("newField = x*2");
  eval->parse();
  EXPECT_TRUE(eval->getSyntaxStatus());
  EXPECT_TRUE(eval->getParseStatus());

  stk::tools::MeshEvalFunctor meshFunc(eval, meta, stk::topology::NODE_RANK);

  stk::mesh::for_each_entity_run(*meshPtr, stk::topology::NODE_RANK, meta.locally_owned_part(), meshFunc);

  const auto* newField = meta.get_field(stk::topology::NODE_RANK, "newField");
  const auto* coordField = meta.coordinate_field();
  EXPECT_TRUE(newField != nullptr);
  EXPECT_TRUE(coordField != nullptr);

  auto newFieldData = newField->data<double,stk::mesh::ReadOnly>();
  auto coordFieldData = coordField->data<double,stk::mesh::ReadOnly>();

  stk::mesh::for_each_entity_run(*meshPtr, stk::topology::NODE_RANK, meta.locally_owned_part(),
    [&](const stk::mesh::BulkData&, const stk::mesh::MeshIndex& mi)
    {
      auto newFieldNodeData = newFieldData.entity_values(mi);
      auto nodeCoord = coordFieldData.entity_values(mi);
      EXPECT_DOUBLE_EQ(2.0*nodeCoord(0_comp), newFieldNodeData(0_comp));
    }
  );
}

} // namespace anonymous

