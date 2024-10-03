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
#include <stddef.h>
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>
#include <string>
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/ForEachEntity.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_transfer/copy_by_id/TransferCopyById.hpp"
#include "stk_transfer/copy_by_id/TransferCopyByIdStkMeshAdapter.hpp"
#include "stk_transfer/copy_by_id/SearchByIdGeometric.hpp"

namespace
{
using DoubleField = stk::mesh::Field<double>;

void set_field_vals_from_node_ids(const stk::mesh::BulkData& mesh,
                                  const DoubleField& field)
{
  stk::mesh::for_each_entity_run(mesh, stk::topology::NODE_RANK,
                                 mesh.mesh_meta_data().locally_owned_part(),
                                 [&field](const stk::mesh::BulkData& bulkData, const stk::mesh::Entity& node)
  {
    double * scalar = stk::mesh::field_data(field, node);
    *scalar = static_cast<double>(bulkData.identifier(node));
  });
}

void change_mesh_decomposition(stk::mesh::BulkData& mesh)
{
  if (mesh.parallel_size() == 1) { return; }
  stk::mesh::Selector owned = mesh.mesh_meta_data().locally_owned_part();
  stk::mesh::EntityProcVec entityProcPairs;
  int myProc = mesh.parallel_rank();
  int otherProc = 1-myProc;

  stk::mesh::for_each_entity_run_no_threads(mesh, stk::topology::ELEM_RANK,
                                            mesh.mesh_meta_data().locally_owned_part(),
                                            [&entityProcPairs, &otherProc](const stk::mesh::BulkData& bulkData, const stk::mesh::Entity& elem)
  {
    entityProcPairs.emplace_back(elem, otherProc);
  });

  stk::mesh::for_each_entity_run_no_threads(mesh, stk::topology::NODE_RANK,
                                 mesh.mesh_meta_data().locally_owned_part(),
                                 [&entityProcPairs, &otherProc](const stk::mesh::BulkData& bulkData, const stk::mesh::Entity& node)
  {
    entityProcPairs.emplace_back(node, otherProc);
  });

  mesh.change_entity_owner(entityProcPairs);
}

//BEGIN
TEST(StkTransferHowTo, useCopyTransfer)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) > 2) { GTEST_SKIP(); }

  const std::string meshSpec("generated:3x3x4");
  double initVals = std::numeric_limits<double>::max();
  const unsigned spatialDim = 3;

  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> meshA = builder.create();
  stk::mesh::MetaData& metaA = meshA->mesh_meta_data();
  DoubleField & scalarFieldNodeA = metaA.declare_field<double>(stk::topology::NODE_RANK, "Node Scalar Field");
  stk::mesh::put_field_on_mesh(scalarFieldNodeA, metaA.universal_part(), &initVals);
  stk::io::fill_mesh(meshSpec, *meshA);

  std::shared_ptr<stk::mesh::BulkData> meshB = builder.create();
  stk::mesh::MetaData& metaB = meshB->mesh_meta_data();
  DoubleField & scalarFieldNodeB = metaB.declare_field<double>(stk::topology::NODE_RANK, "Node Scalar Field");
  stk::mesh::put_field_on_mesh(scalarFieldNodeB, metaB.universal_part(), &initVals);
  stk::io::fill_mesh(meshSpec, *meshB);

  change_mesh_decomposition(*meshB);

  set_field_vals_from_node_ids(*meshA, scalarFieldNodeA);

  // Set up CopyTransfer
  stk::mesh::EntityVector entitiesA;
  stk::mesh::get_entities(*meshA, stk::topology::NODE_RANK, metaA.locally_owned_part(), entitiesA);
  std::vector<stk::mesh::FieldBase*> fieldsA = {&scalarFieldNodeA};
  stk::transfer::TransferCopyByIdStkMeshAdapter transferMeshA(*meshA,entitiesA,fieldsA);

  stk::mesh::EntityVector entitiesB;
  stk::mesh::get_entities(*meshB, stk::topology::NODE_RANK, metaB.locally_owned_part(), entitiesB);
  std::vector<stk::mesh::FieldBase*> fieldsB = {&scalarFieldNodeB};
  stk::transfer::TransferCopyByIdStkMeshAdapter transferMeshB(*meshB,entitiesB,fieldsB);

  stk::transfer::SearchByIdGeometric copySearch;

  stk::transfer::TransferCopyById copyTransfer(copySearch,transferMeshA,transferMeshB);
  copyTransfer.initialize();
  copyTransfer.apply();

  // Verify nodal fields on meshB are correct
  stk::mesh::Selector owned = metaB.locally_owned_part();
  auto check_nodal_fields = [&scalarFieldNodeB](const stk::mesh::BulkData& mesh, const stk::mesh::Entity& node)
  { 
    const double tolerance = 1.0e-8;
    double * scalar = stk::mesh::field_data(scalarFieldNodeB, node);
    EXPECT_NEAR( static_cast<double>(mesh.identifier(node)), *scalar, tolerance);
  };
  stk::mesh::for_each_entity_run(*meshB, stk::topology::NODE_RANK, owned, check_nodal_fields);
}
//END

} // namespace
