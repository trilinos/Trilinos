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

#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field, etc
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_bytes_per_entity, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for operator<<, Selector, etc
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc

namespace {

struct FieldTypeChecker {
  template <typename T>
  void operator()(const stk::mesh::FieldBase& fieldBase) {
    EXPECT_EQ(fieldBase.type_is<T>(), true);
  }
};

TEST(UnitTestFieldDatatypeExecute, functor) {
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::mesh::Field<unsigned>& field = meta.declare_field<unsigned>(stk::topology::NODE_RANK, "unsigned_field");

  FieldTypeChecker fieldTypeChecker{};
  stk::mesh::field_datatype_execute(field, fieldTypeChecker);
}

struct FieldTypeCheckerExtraArg {
  template <typename T>
  void operator()(const stk::mesh::FieldBase& fieldBase, int extraArg) {
    EXPECT_EQ(fieldBase.type_is<T>(), true);
    EXPECT_EQ(extraArg, 3);
  }
};

TEST(UnitTestFieldDatatypeExecute, functorExtraArg) {
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::mesh::Field<unsigned>& field = meta.declare_field<unsigned>(stk::topology::NODE_RANK, "unsigned_field");

  FieldTypeCheckerExtraArg fieldTypeCheckerExtraArg{};
  stk::mesh::field_datatype_execute(field, fieldTypeCheckerExtraArg, 3);
}

TEST(UnitTestFieldDatatypeExecute, templatedLambda) {
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::mesh::Field<unsigned>& field = meta.declare_field<unsigned>(stk::topology::NODE_RANK, "unsigned_field");

  auto fieldTypeChecker = [&]<typename T>(const stk::mesh::FieldBase& fieldBase) {
    EXPECT_EQ(fieldBase.type_is<T>(), true);
  };

  stk::mesh::field_datatype_execute(field, fieldTypeChecker);
}

TEST(UnitTestFieldDatatypeExecute, templatedLambdaExtraArgs) {
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::mesh::Field<unsigned>& field = meta.declare_field<unsigned>(stk::topology::NODE_RANK, "unsigned_field");

  auto fieldTypeChecker = [&]<typename T>(const stk::mesh::FieldBase& fieldBase, int extraArg) {
    EXPECT_EQ(fieldBase.type_is<T>(), true);
    EXPECT_EQ(extraArg, 5);
  };

  stk::mesh::field_datatype_execute(field, fieldTypeChecker, 5);
}

}






