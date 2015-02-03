// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <string>
#include <limits>
#include "gtest/gtest.h"
#include "stk_util/parallel/Parallel.hpp"      // for ParallelMachine, etc
#include "stk_mesh/base/BulkData.hpp"          // for BulkData, etc
#include "stk_mesh/base/MetaData.hpp"          // for MetaData, entity_rank_names, etc
#include "stk_topology/topology.hpp"           // for topology, etc
#include "stk_mesh/base/CoordinateSystems.hpp" // for Cartesian2d, etc.
#include "stk_mesh/base/Part.hpp"              // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"        // for declare_element
#include "stk_mesh/base/GetEntities.hpp"       // for get_selected_entities, etc.

#include "stk_transfer/GeometricTransfer.hpp"
#include "use_cases/STKNode.hpp"
#include "use_cases/LinearInterpolate.hpp"
#include <boost/shared_ptr.hpp>



namespace
{

typedef stk::mesh::Field<double>                         ScalarField;
typedef stk::mesh::Field<double, stk::mesh::Cartesian3d> VectorField;

void build_mesh(stk::mesh::MetaData & meta,
                stk::mesh::BulkData & mesh,
                const size_t num_elements,
                const size_t num_nodes,
                stk::mesh::EntityId element_ids[],
                int element_owner[],
                stk::mesh::EntityId elem_node_ids[][8],
                int node_sharing[],
                double coordinates[][3])
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  double init_vals[] = {std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::max()};

  stk::mesh::Part * elem_part = &meta.declare_part_with_topology("elem_part", stk::topology::HEX_8);
  ScalarField & scalarField = meta.declare_field<ScalarField>(stk::topology::NODE_RANK, "Scalar Field");
  VectorField & vectorField = meta.declare_field<VectorField>(stk::topology::NODE_RANK, "Vector Field");
  VectorField & coordsField = meta.declare_field<VectorField>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field(scalarField, meta.universal_part(), init_vals);
  stk::mesh::put_field(vectorField, meta.universal_part(), init_vals);
  stk::mesh::put_field(coordsField, meta.universal_part(), init_vals);
  meta.commit();

  mesh.modification_begin();
  for (size_t i = 0; i < num_elements; ++i) {
    if (p_rank == element_owner[i]) {
      stk::mesh::declare_element(mesh, *elem_part, element_ids[i], elem_node_ids[i]);
    }
  }
  for (size_t i = 0; i < num_nodes; ++i) {
    if (node_sharing[p_rank*num_nodes+i] != -1) {
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, i+1);
      mesh.add_node_sharing(node, node_sharing[p_rank*num_nodes+i]);
    }
  }
  mesh.modification_end();

  const stk::mesh::BucketVector & nodeBuckets = mesh.get_buckets(stk::topology::NODE_RANK, meta.locally_owned_part());
  for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
      for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
          stk::mesh::Entity node = nodeBucket[nodeIndex];
          double * coordsSource = stk::mesh::field_data(coordsField, node);

          // Assume compact node numbering from 1 so that index into provided coordinates
          // array are a fixed offset from the ID.
          int nodeOffset = mesh.identifier(node) - 1;
          for (size_t i = 0; i < meta.spatial_dimension(); ++i) {
              coordsSource[i] = coordinates[nodeOffset][i];
          }
      }
  }
}

void fill_mesh_values(stk::mesh::BulkData & mesh)
{
    const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    const unsigned spatial_dimension = meta.spatial_dimension();
    ScalarField & scalarField = static_cast<ScalarField&>(*meta.get_field(stk::topology::NODE_RANK, "Scalar Field"));
    VectorField & vectorField = static_cast<VectorField&>(*meta.get_field(stk::topology::NODE_RANK, "Vector Field"));

    const stk::mesh::BucketVector & nodeBuckets = mesh.get_buckets(stk::topology::NODE_RANK, meta.locally_owned_part());
    for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
        stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
        for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
            stk::mesh::Entity node = nodeBucket[nodeIndex];
            double * scalar = stk::mesh::field_data(scalarField, node);
            double * vector = stk::mesh::field_data(vectorField, node);

            *scalar = static_cast<double>(mesh.identifier(node));
            for (unsigned i = 0; i < spatial_dimension; ++i) {
                vector[i] = static_cast<double>((mesh.identifier(node)-1)*spatial_dimension+i);
            }
        }
    }
}


TEST(Transfer, copy0T0)
{
//    ^ Y       ID.owning_proc
//    |
//    |    X         meshA                   meshB
//    .---->    4.0--------3.0          4.0--------3.0
//   /          /|         /|           /|         /|
//  /          / |        / |          / |        / |
// v Z       8.0--------7.0 |        8.0--------7.0 |
//            |  |  1.0  |  |   -->   |  |  1.0  |  |
//            | 1.0------|-2.0        | 1.0------|-2.0
//            | /        | /          | /        | /
//            |/         |/           |/         |/
//           5.0--------6.0          5.0--------6.0
//

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 1) {
    return;
  }

  (void)Intrepid::INTREPID_TOL;  // Get rid of unused variable build error

  const size_t spatial_dimension = 3;
  const size_t num_elements = 1;
  const size_t num_nodes = 8;
  stk::mesh::EntityId element_ids[] = {1};
  int element_owner[] = {0};
  stk::mesh::EntityId elem_node_ids[][8] = { {1, 2, 3, 4, 5, 6, 7, 8} };
  int node_sharing[] = { -1, -1, -1, -1, -1, -1, -1, -1,
                         -1, -1, -1, -1, -1, -1, -1, -1 };
  double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0} };

  // Set up the "source" mesh for the transfer
  //
  stk::mesh::MetaData metaA(spatial_dimension);
  stk::mesh::BulkData meshA(metaA, pm);
  build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_owner, elem_node_ids, node_sharing, coordinates);

  // Set up the "target" mesh for the transfer
  //
  stk::mesh::MetaData metaB(spatial_dimension);
  stk::mesh::BulkData meshB(metaB, pm);
  build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_owner, elem_node_ids, node_sharing, coordinates);

  // Fill "source" fields with valid data
  //
  fill_mesh_values(meshA);

  ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "coordinates"));
  ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "coordinates"));

  // Set up the transfer
  //
  std::vector<stk::mesh::Entity> sourceNodes;
  stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
  std::vector<stk::mesh::Entity> targetNodes;
  stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

  const double radius = 0.25;
  std::vector<stk::mesh::FieldBase*> sourceFields;
  sourceFields.push_back(&scalarSourceField);
  sourceFields.push_back(&vectorSourceField);
  boost::shared_ptr<stk::transfer::STKNode> transferSource(
      new stk::transfer::STKNode(sourceNodes, coordsSourceField, sourceFields, radius));

  std::vector<stk::mesh::FieldBase*> targetFields;
  targetFields.push_back(&scalarTargetField);
  targetFields.push_back(&vectorTargetField);
  boost::shared_ptr<stk::transfer::STKNode> transferTarget(
      new stk::transfer::STKNode(targetNodes, coordsTargetField, targetFields, radius));

  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode,
      class stk::transfer::STKNode
    >
  > transfer(transferSource, transferTarget, "copy0T0 unit test");

  // Do the transfer
  //
  transfer.initialize();
  transfer.apply();

  // Check "target" fields to make sure they hold the expected values
  //
  const double tolerance = 1.e-8;
  const stk::mesh::BucketVector & nodeBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() | metaB.globally_shared_part());
  for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      stk::mesh::Entity node = nodeBucket[nodeIndex];
      double * scalarTarget = stk::mesh::field_data(scalarTargetField, node);
      double * vectorTarget = stk::mesh::field_data(vectorTargetField, node);

      EXPECT_NEAR(static_cast<double>(meshB.identifier(node)), *scalarTarget, tolerance);
      for (size_t i = 0; i < spatial_dimension; ++i) {
        EXPECT_NEAR(static_cast<double>((meshB.identifier(node)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
      }
    }
  }

}

TEST(Transfer, copy0T1)
{
//    ^ Y       ID.owning_proc
//    |
//    |    X         meshA                   meshB
//    .---->    4.0--------3.0          4.1--------3.1
//   /          /|         /|           /|         /|
//  /          / |        / |          / |        / |
// v Z       8.0--------7.0 |        8.1--------7.1 |
//            |  |  1.0  |  |   -->   |  |  1.1  |  |
//            | 1.0------|-2.0        | 1.1------|-2.1
//            | /        | /          | /        | /
//            |/         |/           |/         |/
//           5.0--------6.0          5.1--------6.1
//

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const size_t spatial_dimension = 3;
  const size_t num_elements = 1;
  const size_t num_nodes = 8;
  stk::mesh::EntityId element_ids[] = {1};
  int element_ownerA[] = {0};
  int element_ownerB[] = {1};
  stk::mesh::EntityId elem_node_ids[][8] = { {1, 2, 3, 4, 5, 6, 7, 8} };
  int node_sharing[] = { -1, -1, -1, -1, -1, -1, -1, -1,
                         -1, -1, -1, -1, -1, -1, -1, -1 };
  double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0} };

  // Set up the "source" mesh for the transfer
  //
  stk::mesh::MetaData metaA(spatial_dimension);
  stk::mesh::BulkData meshA(metaA, pm);
  build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_ownerA, elem_node_ids, node_sharing, coordinates);

  // Set up the "target" mesh for the transfer
  //
  stk::mesh::MetaData metaB(spatial_dimension);
  stk::mesh::BulkData meshB(metaB, pm);
  build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_ownerB, elem_node_ids, node_sharing, coordinates);

  // Fill "source" fields with valid data
  //
  fill_mesh_values(meshA);

  ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "coordinates"));
  ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "coordinates"));

  // Set up the transfer
  //
  std::vector<stk::mesh::Entity> sourceNodes;
  stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
  std::vector<stk::mesh::Entity> targetNodes;
  stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

  const double radius = 0.25;
  std::vector<stk::mesh::FieldBase*> sourceFields;
  sourceFields.push_back(&scalarSourceField);
  sourceFields.push_back(&vectorSourceField);
  boost::shared_ptr<stk::transfer::STKNode> transferSource(
      new stk::transfer::STKNode(sourceNodes, coordsSourceField, sourceFields, radius));

  std::vector<stk::mesh::FieldBase*> targetFields;
  targetFields.push_back(&scalarTargetField);
  targetFields.push_back(&vectorTargetField);
  boost::shared_ptr<stk::transfer::STKNode> transferTarget(
      new stk::transfer::STKNode(targetNodes, coordsTargetField, targetFields, radius));

  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode,
      class stk::transfer::STKNode
    >
  > transfer(transferSource, transferTarget, "copy0T1 unit test");

  // Do the transfer
  //
  transfer.initialize();
  transfer.apply();

  // Check "target" fields to make sure they hold the expected values
  //
  const double tolerance = 1.e-8;
  const stk::mesh::BucketVector & nodeBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() | metaB.globally_shared_part());
  for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      stk::mesh::Entity node = nodeBucket[nodeIndex];
      double * scalarTarget = stk::mesh::field_data(scalarTargetField, node);
      double * vectorTarget = stk::mesh::field_data(vectorTargetField, node);

      EXPECT_NEAR(static_cast<double>(meshB.identifier(node)), *scalarTarget, tolerance);
      for (size_t i = 0; i < spatial_dimension; ++i) {
        EXPECT_NEAR(static_cast<double>((meshB.identifier(node)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
      }
    }
  }

}

TEST(Transfer, copy1T0)
{
//    ^ Y       ID.owning_proc
//    |
//    |    X         meshA                   meshB
//    .---->    4.1--------3.1          4.0--------3.0
//   /          /|         /|           /|         /|
//  /          / |        / |          / |        / |
// v Z       8.1--------7.1 |        8.0--------7.0 |
//            |  |  1.1  |  |   -->   |  |  1.0  |  |
//            | 1.1------|-2.1        | 1.0------|-2.0
//            | /        | /          | /        | /
//            |/         |/           |/         |/
//           5.1--------6.1          5.0--------6.0
//

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const size_t spatial_dimension = 3;
  const size_t num_elements = 1;
  const size_t num_nodes = 8;
  stk::mesh::EntityId element_ids[] = {1};
  int element_ownerA[] = {1};
  int element_ownerB[] = {0};
  stk::mesh::EntityId elem_node_ids[][8] = { {1, 2, 3, 4, 5, 6, 7, 8} };
  int node_sharing[] = { -1, -1, -1, -1, -1, -1, -1, -1,
                         -1, -1, -1, -1, -1, -1, -1, -1 };
  double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0} };

  // Set up the "source" mesh for the transfer
  //
  stk::mesh::MetaData metaA(spatial_dimension);
  stk::mesh::BulkData meshA(metaA, pm);
  build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_ownerA, elem_node_ids, node_sharing, coordinates);

  // Set up the "target" mesh for the transfer
  //
  stk::mesh::MetaData metaB(spatial_dimension);
  stk::mesh::BulkData meshB(metaB, pm);
  build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_ownerB, elem_node_ids, node_sharing, coordinates);

  // Fill "source" fields with valid data
  //
  fill_mesh_values(meshA);

  ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "coordinates"));
  ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "coordinates"));

  // Set up the transfer
  //
  std::vector<stk::mesh::Entity> sourceNodes;
  stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
  std::vector<stk::mesh::Entity> targetNodes;
  stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

  const double radius = 0.25;
  std::vector<stk::mesh::FieldBase*> sourceFields;
  sourceFields.push_back(&scalarSourceField);
  sourceFields.push_back(&vectorSourceField);
  boost::shared_ptr<stk::transfer::STKNode> transferSource(
      new stk::transfer::STKNode(sourceNodes, coordsSourceField, sourceFields, radius));

  std::vector<stk::mesh::FieldBase*> targetFields;
  targetFields.push_back(&scalarTargetField);
  targetFields.push_back(&vectorTargetField);
  boost::shared_ptr<stk::transfer::STKNode> transferTarget(
      new stk::transfer::STKNode(targetNodes, coordsTargetField, targetFields, radius));

  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode,
      class stk::transfer::STKNode
    >
  > transfer(transferSource, transferTarget, "copy1T0 unit test");

  // Do the transfer
  //
  transfer.initialize();
  transfer.apply();

  // Check "target" fields to make sure they hold the expected values
  //
  const double tolerance = 1.e-8;
  const stk::mesh::BucketVector & nodeBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() | metaB.globally_shared_part());
  for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      stk::mesh::Entity node = nodeBucket[nodeIndex];
      double * scalarTarget = stk::mesh::field_data(scalarTargetField, node);
      double * vectorTarget = stk::mesh::field_data(vectorTargetField, node);

      EXPECT_NEAR(static_cast<double>(meshB.identifier(node)), *scalarTarget, tolerance);
      for (size_t i = 0; i < spatial_dimension; ++i) {
        EXPECT_NEAR(static_cast<double>((meshB.identifier(node)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
      }
    }
  }
}

TEST(Transfer, copy01T10)
{
//    ^ Y       ID.owning_proc
//    |
//    |    X              meshA                             meshB
//    .---->    4.0--------5.0--------6.1         4.1--------5.0--------6.0
//   /          /|         /|         /|          /|         /|         /|
//  /          / |        / |        / |         / |        / |        / |
// v Z      10.0-------11.0-------12.1 |      10.1-------11.0-------12.0 |
//            |  |  1.0  |  |  2.1  |  |  -->   |  |  1.1  |  |  2.0  |  |
//            | 1.0------|-2.0------|-3.1       | 1.1------|-2.0------|-3.0
//            | /        | /        | /         | /        | /        | /
//            |/         |/         |/          |/         |/         |/
//           7.0--------8.0--------9.1         7.1--------8.0--------9.0
//

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const size_t spatial_dimension = 3;
  const size_t num_elements = 2;
  const size_t num_nodes = 12;
  stk::mesh::EntityId element_ids[] = {1, 2};
  int element_ownerA[] = {0, 1};
  int element_ownerB[] = {1, 0};
  stk::mesh::EntityId elem_node_ids[][8] = { {1, 2, 5, 4, 7, 8, 11, 10},
                                             {2, 3, 6, 5, 8, 9, 12, 11} };
  int node_sharingA[] = { -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1,
                          -1, 0, -1, -1, 0, -1, -1, 0, -1, -1, 0, -1 };
  int node_sharingB[] = { -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1,
                          -1, 0, -1, -1, 0, -1, -1, 0, -1, -1, 0, -1 };
  double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
                              {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0},
                              {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0} };

  // Set up the "source" mesh for the transfer
  //
  stk::mesh::MetaData metaA(spatial_dimension);
  stk::mesh::BulkData meshA(metaA, pm);
  build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_ownerA, elem_node_ids, node_sharingA, coordinates);

  // Set up the "target" mesh for the transfer
  //
  stk::mesh::MetaData metaB(spatial_dimension);
  stk::mesh::BulkData meshB(metaB, pm);
  build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_ownerB, elem_node_ids, node_sharingB, coordinates);

  // Fill "source" fields with valid data
  //
  fill_mesh_values(meshA);

  ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "coordinates"));
  ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "coordinates"));

  // Set up the transfer
  //
  std::vector<stk::mesh::Entity> sourceNodes;
  stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
  std::vector<stk::mesh::Entity> targetNodes;
  stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

  const double radius = 0.25;
  std::vector<stk::mesh::FieldBase*> sourceFields;
  sourceFields.push_back(&scalarSourceField);
  sourceFields.push_back(&vectorSourceField);
  boost::shared_ptr<stk::transfer::STKNode> transferSource(
      new stk::transfer::STKNode(sourceNodes, coordsSourceField, sourceFields, radius));

  std::vector<stk::mesh::FieldBase*> targetFields;
  targetFields.push_back(&scalarTargetField);
  targetFields.push_back(&vectorTargetField);
  boost::shared_ptr<stk::transfer::STKNode> transferTarget(
      new stk::transfer::STKNode(targetNodes, coordsTargetField, targetFields, radius));

  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode,
      class stk::transfer::STKNode
    >
  > transfer(transferSource, transferTarget, "copy01T10 unit test");

  // Do the transfer
  //
  transfer.initialize();
  transfer.apply();

  // Check "target" fields to make sure they hold the expected values.
  //
  const double tolerance = 1.e-8;
  const stk::mesh::BucketVector & nodeBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() | metaB.globally_shared_part());
  for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      stk::mesh::Entity node = nodeBucket[nodeIndex];
      double * scalarTarget = stk::mesh::field_data(scalarTargetField, node);
      double * vectorTarget = stk::mesh::field_data(vectorTargetField, node);

      EXPECT_NEAR(static_cast<double>(meshB.identifier(node)), *scalarTarget, tolerance);
      for (size_t i = 0; i < spatial_dimension; ++i) {
        EXPECT_NEAR(static_cast<double>((meshB.identifier(node)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
      }
    }
  }
}

TEST(Transfer, copy001T011)
{
//    ^ Y       ID.owning_proc
//    |
//    |    X                   meshA                                         meshB
//    .---->    5.0--------6.0--------7.0--------8.1          5.0--------6.0--------7.1--------8.1
//   /          /|         /|         /|         /|           /|         /|         /|         /|
//  /          / |        / |        / |        / |          / |        / |        / |        / |
// v Z      13.0-------14.0-------15.0-------16.1 |       13.0-------14.0-------15.1-------16.1 |
//            |  |  1.0  |  |  2.0  |  |  3.1  |  |   -->   |  |  1.0  |  |  2.1  |  |  3.1  |  |
//            | 1.0------|-2.0------|-3.0------|-4.1        | 1.0------|-2.0------|-3.1------|-4.1
//            | /        | /        | /        | /          | /        | /        | /        | /
//            |/         |/         |/         |/           |/         |/         |/         |/
//           9.0-------10.0-------11.0-------12.1          9.0-------10.0-------11.1-------12.1
//

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const size_t spatial_dimension = 3;
  const size_t num_elements = 3;
  const size_t num_nodes = 16;
  stk::mesh::EntityId element_ids[] = {1, 2, 3};
  int element_ownerA[] = {0, 0, 1};
  int element_ownerB[] = {0, 1, 1};
  stk::mesh::EntityId elem_node_ids[][8] = { {1, 2, 6, 5, 9, 10, 14, 13},
                                             {2, 3, 7, 6, 10, 11, 15, 14},
                                             {3, 4, 8, 7, 11, 12, 16, 15} };
  int node_sharingA[] = { -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1,
                          -1, -1, 0, -1, -1, -1, 0, -1, -1, -1, 0, -1, -1, -1, 0, -1 };
  int node_sharingB[] = { -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1,
                          -1, 0, -1, -1, -1, 0, -1, -1, -1, 0, -1, -1, -1, 0, -1, -1 };
  double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
                              {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0},
                              {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0} };

  // Set up the "source" mesh for the transfer
  //
  stk::mesh::MetaData metaA(spatial_dimension);
  stk::mesh::BulkData meshA(metaA, pm);
  build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_ownerA, elem_node_ids, node_sharingA, coordinates);

  // Set up the "target" mesh for the transfer
  //
  stk::mesh::MetaData metaB(spatial_dimension);
  stk::mesh::BulkData meshB(metaB, pm);
  build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_ownerB, elem_node_ids, node_sharingB, coordinates);

  // Fill "source" fields with valid data
  //
  fill_mesh_values(meshA);

  ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "coordinates"));
  ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "coordinates"));

  // Set up the transfer
  //
  std::vector<stk::mesh::Entity> sourceNodes;
  stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
  std::vector<stk::mesh::Entity> targetNodes;
  stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

  const double radius = 0.25;
  std::vector<stk::mesh::FieldBase*> sourceFields;
  sourceFields.push_back(&scalarSourceField);
  sourceFields.push_back(&vectorSourceField);
  boost::shared_ptr<stk::transfer::STKNode> transferSource(
      new stk::transfer::STKNode(sourceNodes, coordsSourceField, sourceFields, radius));

  std::vector<stk::mesh::FieldBase*> targetFields;
  targetFields.push_back(&scalarTargetField);
  targetFields.push_back(&vectorTargetField);
  boost::shared_ptr<stk::transfer::STKNode> transferTarget(
      new stk::transfer::STKNode(targetNodes, coordsTargetField, targetFields, radius));

  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode,
      class stk::transfer::STKNode
    >
  > transfer(transferSource, transferTarget, "copy001T011 unit test");

  // Do the transfer
  //
  transfer.initialize();
  transfer.apply();

  // Check "target" fields to make sure they hold the expected values.
  //
  const double tolerance = 1.e-8;
  const stk::mesh::BucketVector & nodeBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() | metaB.globally_shared_part());
  for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      stk::mesh::Entity node = nodeBucket[nodeIndex];
      double * scalarTarget = stk::mesh::field_data(scalarTargetField, node);
      double * vectorTarget = stk::mesh::field_data(vectorTargetField, node);

      EXPECT_NEAR(static_cast<double>(meshB.identifier(node)), *scalarTarget, tolerance);
      for (size_t i = 0; i < spatial_dimension; ++i) {
        EXPECT_NEAR(static_cast<double>((meshB.identifier(node)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
      }
    }
  }
}

TEST(Transfer, copy012T000)
{
//    ^ Y       ID.owning_proc
//    |
//    |    X                   meshA                                         meshB
//    .---->    5.0--------6.0--------7.1--------8.2          5.0--------6.0--------7.0--------8.0
//   /          /|         /|         /|         /|           /|         /|         /|         /|
//  /          / |        / |        / |        / |          / |        / |        / |        / |
// v Z      13.0-------14.0-------15.1-------16.2 |       13.0-------14.0-------15.0-------16.0 |
//            |  |  1.0  |  |  2.1  |  |  3.2  |  |   -->   |  |  1.0  |  |  2.0  |  |  3.0  |  |
//            | 1.0------|-2.0------|-3.1------|-4.2        | 1.0------|-2.0------|-3.0------|-4.0
//            | /        | /        | /        | /          | /        | /        | /        | /
//            |/         |/         |/         |/           |/         |/         |/         |/
//           9.0-------10.0-------11.1-------12.2          9.0-------10.0-------11.0-------12.0
//

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 3) {
    return;
  }

  const size_t spatial_dimension = 3;
  const size_t num_elements = 3;
  const size_t num_nodes = 16;
  stk::mesh::EntityId element_ids[] = {1, 2, 3};
  int element_ownerA[] = {0, 1, 2};
  int element_ownerB[] = {0, 0, 0};
  stk::mesh::EntityId elem_node_ids[][8] = { {1, 2, 6, 5, 9, 10, 14, 13},
                                             {2, 3, 7, 6, 10, 11, 15, 14},
                                             {3, 4, 8, 7, 11, 12, 16, 15} };
  int node_sharingA[] = { -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1,
                          -1, 0, 2, -1, -1, 0, 2, -1, -1, 0, 2, -1, -1, 0, 2, -1,
                          -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1 };
  int node_sharingB[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                          -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                          -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
                              {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0},
                              {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0} };

  // Set up the "source" mesh for the transfer
  //
  stk::mesh::MetaData metaA(spatial_dimension);
  stk::mesh::BulkData meshA(metaA, pm);
  build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_ownerA, elem_node_ids, node_sharingA, coordinates);

  // Set up the "target" mesh for the transfer
  //
  stk::mesh::MetaData metaB(spatial_dimension);
  stk::mesh::BulkData meshB(metaB, pm);
  build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_ownerB, elem_node_ids, node_sharingB, coordinates);

  // Fill "source" fields with valid data
  //
  fill_mesh_values(meshA);

  ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "coordinates"));
  ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "coordinates"));

  // Set up the transfer
  //
  std::vector<stk::mesh::Entity> sourceNodes;
  stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
  std::vector<stk::mesh::Entity> targetNodes;
  stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

  const double radius = 0.25;
  std::vector<stk::mesh::FieldBase*> sourceFields;
  sourceFields.push_back(&scalarSourceField);
  sourceFields.push_back(&vectorSourceField);
  boost::shared_ptr<stk::transfer::STKNode> transferSource(
      new stk::transfer::STKNode(sourceNodes, coordsSourceField, sourceFields, radius));

  std::vector<stk::mesh::FieldBase*> targetFields;
  targetFields.push_back(&scalarTargetField);
  targetFields.push_back(&vectorTargetField);
  boost::shared_ptr<stk::transfer::STKNode> transferTarget(
      new stk::transfer::STKNode(targetNodes, coordsTargetField, targetFields, radius));

  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode,
      class stk::transfer::STKNode
    >
  > transfer(transferSource, transferTarget, "copy012T000 unit test");

  // Do the transfer
  //
  transfer.initialize();
  transfer.apply();

  // Check "target" fields to make sure they hold the expected values.
  //
  const double tolerance = 1.e-8;
  const stk::mesh::BucketVector & nodeBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() | metaB.globally_shared_part());
  for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      stk::mesh::Entity node = nodeBucket[nodeIndex];
      double * scalarTarget = stk::mesh::field_data(scalarTargetField, node);
      double * vectorTarget = stk::mesh::field_data(vectorTargetField, node);

      EXPECT_NEAR(static_cast<double>(meshB.identifier(node)), *scalarTarget, tolerance);
      for (size_t i = 0; i < spatial_dimension; ++i) {
        EXPECT_NEAR(static_cast<double>((meshB.identifier(node)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
      }
    }
  }
}

TEST(Transfer, copy000T012)
{
//    ^ Y       ID.owning_proc
//    |
//    |    X                   meshA                                         meshB
//    .---->    5.0--------6.0--------7.0--------8.0         5.0--------6.0--------7.1--------8.2
//   /          /|         /|         /|         /|          /|         /|         /|         /|
//  /          / |        / |        / |        / |         / |        / |        / |        / |
// v Z      13.0-------14.0-------15.0-------16.0 |      13.0-------14.0-------15.1-------16.2 |
//            |  |  1.0  |  |  2.0  |  |  3.0  |  |  -->   |  |  1.0  |  |  2.1  |  |  3.2  |  |
//            | 1.0------|-2.0------|-3.0------|-4.0       | 1.0------|-2.0------|-3.1------|-4.2
//            | /        | /        | /        | /         | /        | /        | /        | /
//            |/         |/         |/         |/          |/         |/         |/         |/
//           9.0-------10.0-------11.0-------12.0         9.0-------10.0-------11.1-------12.2
//

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 3) {
    return;
  }

  const size_t spatial_dimension = 3;
  const size_t num_elements = 3;
  const size_t num_nodes = 16;
  stk::mesh::EntityId element_ids[] = {1, 2, 3};
  int element_ownerA[] = {0, 0, 0};
  int element_ownerB[] = {0, 1, 2};
  stk::mesh::EntityId elem_node_ids[][8] = { {1, 2, 6, 5, 9, 10, 14, 13},
                                             {2, 3, 7, 6, 10, 11, 15, 14},
                                             {3, 4, 8, 7, 11, 12, 16, 15} };
  int node_sharingA[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                          -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                          -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  int node_sharingB[] = { -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1,
                          -1, 0, 2, -1, -1, 0, 2, -1, -1, 0, 2, -1, -1, 0, 2, -1,
                          -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1 };
  double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
                              {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0},
                              {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0} };

  // Set up the "source" mesh for the transfer
  //
  stk::mesh::MetaData metaA(spatial_dimension);
  stk::mesh::BulkData meshA(metaA, pm);
  build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_ownerA, elem_node_ids, node_sharingA, coordinates);

  // Set up the "target" mesh for the transfer
  //
  stk::mesh::MetaData metaB(spatial_dimension);
  stk::mesh::BulkData meshB(metaB, pm);
  build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_ownerB, elem_node_ids, node_sharingB, coordinates);

  // Fill "source" fields with valid data
  //
  fill_mesh_values(meshA);

  ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "coordinates"));
  ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "coordinates"));

  // Set up the transfer
  //
  std::vector<stk::mesh::Entity> sourceNodes;
  stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
  std::vector<stk::mesh::Entity> targetNodes;
  stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

  const double radius = 0.25;
  std::vector<stk::mesh::FieldBase*> sourceFields;
  sourceFields.push_back(&scalarSourceField);
  sourceFields.push_back(&vectorSourceField);
  boost::shared_ptr<stk::transfer::STKNode> transferSource(
      new stk::transfer::STKNode(sourceNodes, coordsSourceField, sourceFields, radius));

  std::vector<stk::mesh::FieldBase*> targetFields;
  targetFields.push_back(&scalarTargetField);
  targetFields.push_back(&vectorTargetField);
  boost::shared_ptr<stk::transfer::STKNode> transferTarget(
      new stk::transfer::STKNode(targetNodes, coordsTargetField, targetFields, radius));

  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode,
      class stk::transfer::STKNode
    >
  > transfer(transferSource, transferTarget, "copy000T012 unit test");

  // Do the transfer
  //
  transfer.initialize();
  transfer.apply();

  // Check "target" fields to make sure they hold the expected values.
  //
  const double tolerance = 1.e-8;
  const stk::mesh::BucketVector & nodeBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() | metaB.globally_shared_part());
  for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      stk::mesh::Entity node = nodeBucket[nodeIndex];
      double * scalarTarget = stk::mesh::field_data(scalarTargetField, node);
      double * vectorTarget = stk::mesh::field_data(vectorTargetField, node);

      EXPECT_NEAR(static_cast<double>(meshB.identifier(node)), *scalarTarget, tolerance);
      for (size_t i = 0; i < spatial_dimension; ++i) {
        EXPECT_NEAR(static_cast<double>((meshB.identifier(node)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
      }
    }
  }
}

TEST(Transfer, copy0011T1010)
{
//    ^ Y       ID.owning_proc
//    |
//    |    X                         meshA                                                     meshB
//    .---->    6.0--------7.0--------8.0--------9.1-------10.1           6.1--------7.0--------8.0--------9.0-------10.0
//   /          /|         /|         /|         /|         /|            /|         /|         /|         /|         /|
//  /          / |        / |        / |        / |        / |           / |        / |        / |        / |        / |
// v Z      16.0-------17.0-------18.0-------19.1-------20.1 |        16.1-------17.0-------18.0-------19.0-------20.0 |
//            |  |  1.0  |  |  2.0  |  |  3.1  |  |  4.1  |  |   -->    |  |  1.1  |  |  2.0  |  |  3.1  |  |  4.0  |  |
//            | 1.0------|-2.0------|-3.0------|-4.1------|-5.1         | 1.1------|-2.0------|-3.0------|-4.0------|-5.0
//            | /        | /        | /        | /        | /           | /        | /        | /        | /        | /
//            |/         |/         |/         |/         |/            |/         |/         |/         |/         |/
//          11.0-------12.0-------13.0-------14.1-------15.1          11.1-------12.0-------13.0-------14.0-------15.0
//

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const size_t spatial_dimension = 3;
  const size_t num_elements = 4;
  const size_t num_nodes = 20;
  stk::mesh::EntityId element_ids[] = {1, 2, 3, 4};
  int element_ownerA[] = {0, 0, 1, 1};
  int element_ownerB[] = {1, 0, 1, 0};
  stk::mesh::EntityId elem_node_ids[][8] = { {1, 2, 7, 6, 11, 12, 17, 16},
                                             {2, 3, 8, 7, 12, 13, 18, 17},
                                             {3, 4, 9, 8, 13, 14, 19, 18},
                                             {4, 5, 10, 9, 14, 15, 20, 19} };
  int node_sharingA[] = { -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1,
                          -1, -1, 0, -1, -1, -1, -1, 0, -1, -1, -1, -1, 0, -1, -1, -1, -1, 0, -1, -1 };
  int node_sharingB[] = { -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1,
                          -1, 0, 0, 0, -1, -1, 0, 0, 0, -1, -1, 0, 0, 0, -1, -1, 0, 0, 0, -1 };
  double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {4.0, 0.0, 0.0},
                              {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0}, {4.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0}, {4.0, 0.0, 1.0},
                              {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0}, {4.0, 1.0, 1.0} };

  // Set up the "source" mesh for the transfer
  //
  stk::mesh::MetaData metaA(spatial_dimension);
  stk::mesh::BulkData meshA(metaA, pm);
  build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_ownerA, elem_node_ids, node_sharingA, coordinates);

  // Set up the "target" mesh for the transfer
  //
  stk::mesh::MetaData metaB(spatial_dimension);
  stk::mesh::BulkData meshB(metaB, pm);
  build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_ownerB, elem_node_ids, node_sharingB, coordinates);

  // Fill "source" fields with valid data
  //
  fill_mesh_values(meshA);

  ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "coordinates"));
  ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "coordinates"));

  // Set up the transfer
  //
  std::vector<stk::mesh::Entity> sourceNodes;
  stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
  std::vector<stk::mesh::Entity> targetNodes;
  stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

  const double radius = 0.25;
  std::vector<stk::mesh::FieldBase*> sourceFields;
  sourceFields.push_back(&scalarSourceField);
  sourceFields.push_back(&vectorSourceField);
  boost::shared_ptr<stk::transfer::STKNode> transferSource(
      new stk::transfer::STKNode(sourceNodes, coordsSourceField, sourceFields, radius));

  std::vector<stk::mesh::FieldBase*> targetFields;
  targetFields.push_back(&scalarTargetField);
  targetFields.push_back(&vectorTargetField);
  boost::shared_ptr<stk::transfer::STKNode> transferTarget(
      new stk::transfer::STKNode(targetNodes, coordsTargetField, targetFields, radius));

  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode,
      class stk::transfer::STKNode
    >
  > transfer(transferSource, transferTarget, "copy0011T1010 unit test");

  // Do the transfer
  //
  transfer.initialize();
  transfer.apply();

  // Check "target" fields to make sure they hold the expected values.
  //
  const double tolerance = 1.e-8;
  const stk::mesh::BucketVector & nodeBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() | metaB.globally_shared_part());
  for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      stk::mesh::Entity node = nodeBucket[nodeIndex];
      double * scalarTarget = stk::mesh::field_data(scalarTargetField, node);
      double * vectorTarget = stk::mesh::field_data(vectorTargetField, node);

      EXPECT_NEAR(static_cast<double>(meshB.identifier(node)), *scalarTarget, tolerance);
      for (size_t i = 0; i < spatial_dimension; ++i) {
        EXPECT_NEAR(static_cast<double>((meshB.identifier(node)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
      }
    }
  }
}

TEST(Transfer, copy0T_)
{
//    ^ Y       ID.owning_proc
//    |
//    |    X         meshA              meshB
//    .---->    4.0--------3.0
//   /          /|         /|
//  /          / |        / |
// v Z       8.0--------7.0 |
//            |  |  1.0  |  |   -->    (empty)
//            | 1.0------|-2.0
//            | /        | /
//            |/         |/
//           5.0--------6.0
//

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 1) {
    return;
  }

  (void)Intrepid::INTREPID_TOL;  // Get rid of unused variable build error

  const size_t spatial_dimension = 3;
  const size_t num_elements = 1;
  const size_t num_nodes = 8;
  stk::mesh::EntityId element_ids[] = {1};
  int element_owner[] = {0};
  stk::mesh::EntityId elem_node_ids[][8] = { {1, 2, 3, 4, 5, 6, 7, 8} };
  int node_sharing[] = { -1, -1, -1, -1, -1, -1, -1, -1,
                         -1, -1, -1, -1, -1, -1, -1, -1 };
  double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0} };

  // Set up the "source" mesh for the transfer
  //
  stk::mesh::MetaData metaA(spatial_dimension);
  stk::mesh::BulkData meshA(metaA, pm);
  build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_owner, elem_node_ids, node_sharing, coordinates);

  // Set up the "target" mesh for the transfer without creating any elements
  //
  stk::mesh::MetaData metaB(spatial_dimension);
  stk::mesh::BulkData meshB(metaB, pm);

  double init_vals[] = {std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::max()};
  ScalarField & scalarTargetField = metaB.declare_field<ScalarField>(stk::topology::NODE_RANK, "Scalar Field");
  VectorField & vectorTargetField = metaB.declare_field<VectorField>(stk::topology::NODE_RANK, "Vector Field");
  VectorField & coordsTargetField = metaB.declare_field<VectorField>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field(scalarTargetField, metaB.universal_part(), init_vals);
  stk::mesh::put_field(vectorTargetField, metaB.universal_part(), init_vals);
  stk::mesh::put_field(coordsTargetField, metaB.universal_part(), init_vals);
  metaB.commit();

  // Fill "source" fields with valid data
  //
  fill_mesh_values(meshA);

  ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "coordinates"));

  // Set up the transfer
  //
  std::vector<stk::mesh::Entity> sourceNodes;
  stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
  std::vector<stk::mesh::Entity> targetNodes;
  stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

  const double radius = 0.25;
  std::vector<stk::mesh::FieldBase*> sourceFields;
  sourceFields.push_back(&scalarSourceField);
  sourceFields.push_back(&vectorSourceField);
  boost::shared_ptr<stk::transfer::STKNode> transferSource(
      new stk::transfer::STKNode(sourceNodes, coordsSourceField, sourceFields, radius));

  std::vector<stk::mesh::FieldBase*> targetFields;
  targetFields.push_back(&scalarTargetField);
  targetFields.push_back(&vectorTargetField);
  boost::shared_ptr<stk::transfer::STKNode> transferTarget(
      new stk::transfer::STKNode(targetNodes, coordsTargetField, targetFields, radius));

  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode,
      class stk::transfer::STKNode
    >
  > transfer(transferSource, transferTarget, "copy0T_ unit test");

  // Do the transfer
  //
  transfer.initialize();
  transfer.apply();

  // Check "target" fields to make sure they hold the expected values
  //
  const double tolerance = 1.e-8;
  const stk::mesh::BucketVector & nodeBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() | metaB.globally_shared_part());
  for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      stk::mesh::Entity node = nodeBucket[nodeIndex];
      double * scalarTarget = stk::mesh::field_data(scalarTargetField, node);
      double * vectorTarget = stk::mesh::field_data(vectorTargetField, node);

      EXPECT_NEAR(static_cast<double>(meshB.identifier(node)), *scalarTarget, tolerance);
      for (size_t i = 0; i < spatial_dimension; ++i) {
        EXPECT_NEAR(static_cast<double>((meshB.identifier(node)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
      }
    }
  }
}

TEST(Transfer, copy_T0)
{
//    ^ Y       ID.owning_proc
//    |
//    |    X      meshA                 meshB
//    .---->                        4.0--------3.0
//   /                              /|         /|
//  /                              / |        / |
// v Z                           8.0--------7.0 |
//               (empty)  -->     |  |  1.0  |  |
//                                | 1.0------|-2.0
//                                | /        | /
//                                |/         |/
//                               5.0--------6.0
//

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 1) {
    return;
  }

  (void)Intrepid::INTREPID_TOL;  // Get rid of unused variable build error

  const size_t spatial_dimension = 3;
  const size_t num_elements = 1;
  const size_t num_nodes = 8;
  stk::mesh::EntityId element_ids[] = {1};
  int element_owner[] = {0};
  stk::mesh::EntityId elem_node_ids[][8] = { {1, 2, 3, 4, 5, 6, 7, 8} };
  int node_sharing[] = { -1, -1, -1, -1, -1, -1, -1, -1,
                         -1, -1, -1, -1, -1, -1, -1, -1 };
  double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0} };

  // Set up the "source" mesh for the transfer
  //
  stk::mesh::MetaData metaA(spatial_dimension);
  stk::mesh::BulkData meshA(metaA, pm);

  double init_vals[] = {std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::max()};
  ScalarField & scalarSourceField = metaA.declare_field<ScalarField>(stk::topology::NODE_RANK, "Scalar Field");
  VectorField & vectorSourceField = metaA.declare_field<VectorField>(stk::topology::NODE_RANK, "Vector Field");
  VectorField & coordsSourceField = metaA.declare_field<VectorField>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field(scalarSourceField, metaA.universal_part(), init_vals);
  stk::mesh::put_field(vectorSourceField, metaA.universal_part(), init_vals);
  stk::mesh::put_field(coordsSourceField, metaA.universal_part(), init_vals);
  metaA.commit();

  // Set up the "target" mesh for the transfer without creating any elements
  //
  stk::mesh::MetaData metaB(spatial_dimension);
  stk::mesh::BulkData meshB(metaB, pm);
  build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_owner, elem_node_ids, node_sharing, coordinates);

  ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "coordinates"));

  // Set up the transfer
  //
  std::vector<stk::mesh::Entity> sourceNodes;
  stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
  std::vector<stk::mesh::Entity> targetNodes;
  stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

  const double radius = 0.25;
  std::vector<stk::mesh::FieldBase*> sourceFields;
  sourceFields.push_back(&scalarSourceField);
  sourceFields.push_back(&vectorSourceField);
  boost::shared_ptr<stk::transfer::STKNode> transferSource(
      new stk::transfer::STKNode(sourceNodes, coordsSourceField, sourceFields, radius));

  std::vector<stk::mesh::FieldBase*> targetFields;
  targetFields.push_back(&scalarTargetField);
  targetFields.push_back(&vectorTargetField);
  boost::shared_ptr<stk::transfer::STKNode> transferTarget(
      new stk::transfer::STKNode(targetNodes, coordsTargetField, targetFields, radius));

  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode,
      class stk::transfer::STKNode
    >
  > transfer(transferSource, transferTarget, "copy_T0 unit test");

  // Do the transfer
  //
  // FIXME: Uncomment when convert to the real copy transfer.  This hangs.
//  transfer.initialize();
//  transfer.apply();
//
//  // Check "target" fields to make sure they hold the expected values
//  //
//  const double tolerance = 1.e-8;
//  const stk::mesh::BucketVector & nodeBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() | metaB.globally_shared_part());
//  for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
//    stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
//    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
//      stk::mesh::Entity node = nodeBucket[nodeIndex];
//      double * scalarTarget = stk::mesh::field_data(scalarTargetField, node);
//      double * vectorTarget = stk::mesh::field_data(vectorTargetField, node);
//
//      std::cout << "scalarTarget[" << meshB.identifier(node) << "]=" << *scalarTarget << std::endl;
//      EXPECT_NEAR(static_cast<double>(meshB.identifier(node)), *scalarTarget, tolerance);
//      for (size_t i = 0; i < spatial_dimension; ++i) {
//        EXPECT_NEAR(static_cast<double>((meshB.identifier(node)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
//      }
//    }
//  }
}

TEST(Transfer, copy00_T_11)
{
//    ^ Y       ID.owning_proc
//    |
//    |    X              meshA                             meshB
//    .---->    5.0--------6.0--------7.0         6.1--------7.1--------8.1
//   /          /|         /|         /|          /|         /|         /|
//  /          / |        / |        / |         / |        / |        / |
// v Z      13.0-------14.0-------15.0 |      14.1-------15.1-------16.1 |
//            |  |  1.0  |  |  2.0  |  |  -->   |  |  2.1  |  |  3.1  |  |
//            | 1.0------|-2.0------|-3.0       | 2.1------|-3.1------|-4.1
//            | /        | /        | /         | /        | /        | /
//            |/         |/         |/          |/         |/         |/
//           9.0-------10.0-------11.0        10.1-------11.1-------12.1
//
//  Test of three-element mesh where the source and target are different (overlapping)
//  subsections of the mesh.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const size_t spatial_dimension = 3;
  const size_t num_elements = 3;
  const size_t num_nodes = 16;
  stk::mesh::EntityId element_ids[] = {1, 2, 3};
  int element_ownerA[] = {0, 0, -1};
  int element_ownerB[] = {-1, 1, 1};
  stk::mesh::EntityId elem_node_ids[][8] = { {1, 2, 6, 5, 9, 10, 14, 13},
                                             {2, 3, 7, 6, 10, 11, 15, 14},
                                             {3, 4, 8, 7, 11, 12, 16, 15} };
  int node_sharingA[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                          -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  int node_sharingB[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                          -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
                              {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0},
                              {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0} };

  // Set up the "source" mesh for the transfer
  //
  stk::mesh::MetaData metaA(spatial_dimension);
  stk::mesh::BulkData meshA(metaA, pm);
  build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_ownerA, elem_node_ids, node_sharingA, coordinates);

  // Set up the "target" mesh for the transfer
  //
  stk::mesh::MetaData metaB(spatial_dimension);
  stk::mesh::BulkData meshB(metaB, pm);
  build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_ownerB, elem_node_ids, node_sharingB, coordinates);

  // Fill "source" fields with valid data
  //
  fill_mesh_values(meshA);

  ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "coordinates"));
  ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "coordinates"));

  // Set up the transfer
  //
  std::vector<stk::mesh::Entity> sourceNodes;
  stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
  std::vector<stk::mesh::Entity> targetNodes;
  stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

  const double radius = 0.25;
  std::vector<stk::mesh::FieldBase*> sourceFields;
  sourceFields.push_back(&scalarSourceField);
  sourceFields.push_back(&vectorSourceField);
  boost::shared_ptr<stk::transfer::STKNode> transferSource(
      new stk::transfer::STKNode(sourceNodes, coordsSourceField, sourceFields, radius));

  std::vector<stk::mesh::FieldBase*> targetFields;
  targetFields.push_back(&scalarTargetField);
  targetFields.push_back(&vectorTargetField);
  boost::shared_ptr<stk::transfer::STKNode> transferTarget(
      new stk::transfer::STKNode(targetNodes, coordsTargetField, targetFields, radius));

  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode,
      class stk::transfer::STKNode
    >
  > transfer(transferSource, transferTarget, "copy01T10 unit test");

  // Do the transfer
  //
  transfer.initialize();
  transfer.apply();

  // Check "target" fields to make sure they hold the expected values.
  //
  const double tolerance = 1.e-8;
  const stk::mesh::BucketVector & nodeBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() | metaB.globally_shared_part());
  for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      stk::mesh::Entity node = nodeBucket[nodeIndex];
      double * scalarTarget = stk::mesh::field_data(scalarTargetField, node);
      double * vectorTarget = stk::mesh::field_data(vectorTargetField, node);

      EXPECT_NEAR(static_cast<double>(meshB.identifier(node)), *scalarTarget, tolerance);
      for (size_t i = 0; i < spatial_dimension; ++i) {
        EXPECT_NEAR(static_cast<double>((meshB.identifier(node)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
      }
    }
  }
}

TEST(Transfer, copy00___T___11)
{
//    ^ Y       ID.owning_proc
//    |
//    |    X              meshA                             meshB
//    .---->    7.0--------8.0--------9.0        10.1-------11.1-------12.1
//   /          /|         /|         /|          /|         /|         /|
//  /          / |        / |        / |         / |        / |        / |
// v Z      19.0-------20.0-------21.0 |      22.1-------23.1-------24.1 |
//            |  |  1.0  |  |  2.0  |  |  -->   |  |  2.1  |  |  3.1  |  |
//            | 1.0------|-2.0------|-3.0       | 4.1------|-5.1------|-6.1
//            | /        | /        | /         | /        | /        | /
//            |/         |/         |/          |/         |/         |/
//          13.0-------14.0-------15.0        16.1-------17.1-------18.1
//
//  Test of three-element mesh where the source and target are different (non-overlapping)
//  subsections of the mesh.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const size_t spatial_dimension = 3;
  const size_t num_elements = 5;
  const size_t num_nodes = 24;
  stk::mesh::EntityId element_ids[] = {1, 2, 3, 4, 5};
  int element_ownerA[] = {0, 0, -1, -1, -1};
  int element_ownerB[] = {-1, -1, -1, 1, 1};
  stk::mesh::EntityId elem_node_ids[][8] = { {1, 2, 8, 7, 13, 14, 20, 19},
                                             {2, 3, 9, 8, 14, 15, 21, 20},
                                             {3, 4, 10, 9, 15, 16, 22, 21},
                                             {4, 5, 11, 10, 16, 17, 23, 22},
                                             {5, 6, 12, 11, 17, 18, 24, 23} };
  int node_sharingA[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                          -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  int node_sharingB[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                          -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {4.0, 0.0, 0.0}, {5.0, 0.0, 0.0},
                              {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0}, {4.0, 1.0, 0.0}, {5.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0}, {4.0, 0.0, 1.0}, {5.0, 0.0, 1.0},
                              {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0}, {4.0, 1.0, 1.0}, {5.0, 1.0, 1.0} };

  // Set up the "source" mesh for the transfer
  //
  stk::mesh::MetaData metaA(spatial_dimension);
  stk::mesh::BulkData meshA(metaA, pm);
  build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_ownerA, elem_node_ids, node_sharingA, coordinates);

  // Set up the "target" mesh for the transfer
  //
  stk::mesh::MetaData metaB(spatial_dimension);
  stk::mesh::BulkData meshB(metaB, pm);
  build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_ownerB, elem_node_ids, node_sharingB, coordinates);

  // Fill "source" fields with valid data
  //
  fill_mesh_values(meshA);

  ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "coordinates"));
  ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Scalar Field"));
  VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Vector Field"));
  VectorField & coordsTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "coordinates"));

  // Set up the transfer
  //
  std::vector<stk::mesh::Entity> sourceNodes;
  stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
  std::vector<stk::mesh::Entity> targetNodes;
  stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

  const double radius = 0.25;
  std::vector<stk::mesh::FieldBase*> sourceFields;
  sourceFields.push_back(&scalarSourceField);
  sourceFields.push_back(&vectorSourceField);
  boost::shared_ptr<stk::transfer::STKNode> transferSource(
      new stk::transfer::STKNode(sourceNodes, coordsSourceField, sourceFields, radius));

  std::vector<stk::mesh::FieldBase*> targetFields;
  targetFields.push_back(&scalarTargetField);
  targetFields.push_back(&vectorTargetField);
  boost::shared_ptr<stk::transfer::STKNode> transferTarget(
      new stk::transfer::STKNode(targetNodes, coordsTargetField, targetFields, radius));

  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode,
      class stk::transfer::STKNode
    >
  > transfer(transferSource, transferTarget, "copy01T10 unit test");

  // Do the transfer
  //
  transfer.initialize();
  transfer.apply();

  // Check "target" fields to make sure they hold the expected values.
  //
  const double tolerance = 1.e-8;
  const stk::mesh::BucketVector & nodeBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() | metaB.globally_shared_part());
  for (size_t bucketIndex = 0; bucketIndex < nodeBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & nodeBucket = * nodeBuckets[bucketIndex];
    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      stk::mesh::Entity node = nodeBucket[nodeIndex];
      double * scalarTarget = stk::mesh::field_data(scalarTargetField, node);
      double * vectorTarget = stk::mesh::field_data(vectorTargetField, node);

      EXPECT_NEAR(static_cast<double>(meshB.identifier(node)), *scalarTarget, tolerance);
      for (size_t i = 0; i < spatial_dimension; ++i) {
        EXPECT_NEAR(static_cast<double>((meshB.identifier(node)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
      }
    }
  }
}

}

