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
#include "stk_mesh/base/FEMHelpers.hpp"        // for declare_element
#include "stk_topology/topology.hpp"           // for topology, etc
#include "stk_mesh/base/CoordinateSystems.hpp" // for Cartesian2d, etc.
#include "stk_mesh/base/Part.hpp"              // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/GetEntities.hpp"       // for get_selected_entities, etc.

#include "CopySearchCommAll.hpp"
#include "CopySearchGeometric.hpp"
#include "CopyTransfer.hpp"
#include "CopyTransferStkMeshAdapter.hpp"



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
                double coordinates[][3],
                stk::mesh::EntityId face_node_ids[][4]=NULL,
                stk::mesh::EntityId elem_face_ids[][6]=NULL)
{
  const int p_rank = mesh.parallel_rank();
  double init_vals[] = {std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::max()};

  stk::mesh::Part * elem_part = &meta.declare_part_with_topology("elem_part", stk::topology::HEX_8);
  stk::mesh::Part * face_part = &meta.declare_part_with_topology("face_part", stk::topology::QUAD_4);
  stk::mesh::Part * shell_part = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  ScalarField & scalarFieldNode = meta.declare_field<ScalarField>(stk::topology::NODE_RANK, "Node Scalar Field");
  VectorField & vectorFieldNode = meta.declare_field<VectorField>(stk::topology::NODE_RANK, "Node Vector Field");
  VectorField & coordsFieldNode = meta.declare_field<VectorField>(stk::topology::NODE_RANK, "coordinates");
  meta.set_coordinate_field(&coordsFieldNode);
  stk::mesh::put_field(scalarFieldNode, meta.universal_part(), init_vals);
  stk::mesh::put_field(vectorFieldNode, meta.universal_part(), init_vals);
  stk::mesh::put_field(coordsFieldNode, meta.universal_part(), init_vals);
  ScalarField & scalarFieldElement = meta.declare_field<ScalarField>(stk::topology::ELEM_RANK, "Element Scalar Field");
  VectorField & vectorFieldElement = meta.declare_field<VectorField>(stk::topology::ELEM_RANK, "Element Vector Field");
  stk::mesh::put_field(scalarFieldElement, meta.universal_part(), init_vals);
  stk::mesh::put_field(vectorFieldElement, meta.universal_part(), init_vals);
  ScalarField & scalarFieldFace = meta.declare_field<ScalarField>(stk::topology::FACE_RANK, "Face Scalar Field");
  VectorField & vectorFieldFace = meta.declare_field<VectorField>(stk::topology::FACE_RANK, "Face Vector Field");
  stk::mesh::put_field(scalarFieldFace, meta.universal_part(), init_vals);
  stk::mesh::put_field(vectorFieldFace, meta.universal_part(), init_vals);
  ScalarField & scalarFieldShell = meta.declare_field<ScalarField>(stk::topology::ELEM_RANK, "Shell Scalar Field");
  VectorField & vectorFieldShell = meta.declare_field<VectorField>(stk::topology::ELEM_RANK, "Shell Vector Field");
  stk::mesh::put_field(scalarFieldShell, *shell_part, init_vals);
  stk::mesh::put_field(vectorFieldShell, *shell_part, init_vals);
  meta.commit();

  mesh.modification_begin();
  for (size_t i = 0; i < num_elements; ++i) {
    if (p_rank == element_owner[i]) {
      stk::mesh::declare_element(mesh, *elem_part, element_ids[i], elem_node_ids[i]);
    }
  }
  if (elem_face_ids != NULL && face_node_ids != NULL) {
    stk::mesh::PartVector add_parts;
    add_parts.push_back(face_part);
    for (size_t i = 0; i < num_elements; ++i) {
      if (p_rank == element_owner[i]) {
        stk::mesh::Entity element = mesh.get_entity(stk::topology::ELEM_RANK,element_ids[i]);
        for (size_t side_id = 0; side_id < 6; ++side_id) {
          stk::mesh::EntityId side_global_id = elem_face_ids[i][side_id];
          stk::mesh::Entity side = mesh.get_entity(stk::topology::FACE_RANK,side_global_id);
          if (!mesh.is_valid(side)) {
            side = mesh.declare_entity(stk::topology::FACE_RANK,side_global_id,add_parts);
            const size_t face_node_id_index = i*6+side_id;
            for (size_t node_index = 0; node_index < 4 ; ++node_index) {
              stk::mesh::EntityId node_global_id = face_node_ids[face_node_id_index][node_index];
              stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK,node_global_id);
              mesh.declare_relation(side,node,node_index);
            }
          }
          mesh.declare_relation(element,side,side_id);
        }
      }
    }
  }
  for (size_t i = 0; i < num_nodes; ++i) {
    if (node_sharing[p_rank*num_nodes+i] != -1) {
      stk::mesh::Entity entity = mesh.get_entity(stk::topology::NODE_RANK, i+1);
      mesh.add_node_sharing(entity, node_sharing[p_rank*num_nodes+i]);
    }
  }
  mesh.modification_end();

  const stk::mesh::BucketVector & entityBuckets = mesh.get_buckets(stk::topology::NODE_RANK, meta.locally_owned_part());
  for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
      for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
          stk::mesh::Entity entity = entityBucket[entityIndex];
          double * coordsSource = stk::mesh::field_data(coordsFieldNode, entity);

          // Assume compact entity numbering from 1 so that index into provided coordinates
          // array are a fixed offset from the ID.
          int nodeOffset = mesh.identifier(entity) - 1;
          for (size_t i = 0; i < meta.spatial_dimension(); ++i) {
              coordsSource[i] = coordinates[nodeOffset][i];
          }
      }
  }
  std::vector<const stk::mesh::FieldBase *> fields_to_communicate;
  fields_to_communicate.push_back(&coordsFieldNode);
  stk::mesh::copy_owned_to_shared(mesh,fields_to_communicate);
}

void add_shells_to_mesh(stk::mesh::MetaData & meta,
                        stk::mesh::BulkData & mesh,
                        const size_t num_elements,
                        stk::mesh::EntityId element_ids[],
                        stk::mesh::EntityId shell_node_ids[][4],
                        stk::mesh::EntityId elem_shell_ids[][6],
                        int shell_owner_by_elem_side[][6],
                        const size_t num_shells)
{
  const int p_rank = mesh.parallel_rank();
  stk::mesh::PartVector add_parts;
  add_parts.push_back(meta.get_part("shell_part"));
  mesh.modification_begin();
  for (size_t i = 0; i < num_elements; ++i) {
    stk::mesh::Entity element = mesh.get_entity(stk::topology::ELEM_RANK,element_ids[i]);
    if (mesh.is_valid(element) && mesh.bucket(element).owned()) {
      for (size_t side_id = 0; side_id < 6; ++side_id) {
        stk::mesh::EntityId shell_global_id = elem_shell_ids[i][side_id];
        stk::mesh::Entity shell = mesh.get_entity(stk::topology::ELEM_RANK,shell_global_id);
        const bool i_own_shell = (p_rank == shell_owner_by_elem_side[i][side_id]);
        if (!mesh.is_valid(shell) && i_own_shell) {
          for (size_t shell_index=0; shell_index<num_shells; ++shell_index) {
            shell = mesh.declare_entity(stk::topology::ELEM_RANK,
                                        shell_global_id+100*shell_index,
                                        add_parts);
            const size_t shell_node_id_index = i*6+side_id;
            for (size_t node_index = 0; node_index < 4 ; ++node_index) {
              stk::mesh::EntityId node_global_id = shell_node_ids[shell_node_id_index][node_index];
              stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK,node_global_id);
              mesh.declare_relation(shell,node,node_index);
            }
          }
        }
      }
    }
  }
  mesh.modification_end();
}

void fill_mesh_values_for_rank(stk::mesh::BulkData & mesh,
                               stk::mesh::EntityRank rank,
                               ScalarField & sF,
                               VectorField & vF,
                               const stk::mesh::Selector & sel)
{
    const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    const unsigned spatial_dimension = meta.spatial_dimension();
    ScalarField & scalarField = sF;
    VectorField & vectorField = vF;

    const stk::mesh::BucketVector & entityBuckets = mesh.get_buckets(rank,sel);
    for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
        stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
        for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
            stk::mesh::Entity entity = entityBucket[entityIndex];
            double * scalar = stk::mesh::field_data(scalarField, entity);
            double * vector = stk::mesh::field_data(vectorField, entity);

            *scalar = static_cast<double>(mesh.identifier(entity));
            for (unsigned i = 0; i < spatial_dimension; ++i) {
                vector[i] = static_cast<double>((mesh.identifier(entity)-1)*spatial_dimension+i);
            }
        }
    }
}

void fill_mesh_values(stk::mesh::BulkData & mesh)
{
    const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    ScalarField & scalarField =
      static_cast<ScalarField&>(*meta.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorField =
      static_cast<VectorField&>(*meta.get_field(stk::topology::NODE_RANK, "Node Vector Field"));

    fill_mesh_values_for_rank(mesh,stk::topology::NODE_RANK,scalarField,vectorField,meta.locally_owned_part());

    ScalarField & scalarFieldElement =
      static_cast<ScalarField&>(*meta.get_field(stk::topology::ELEM_RANK, "Element Scalar Field"));
    VectorField & vectorFieldElement =
      static_cast<VectorField&>(*meta.get_field(stk::topology::ELEM_RANK, "Element Vector Field"));

    fill_mesh_values_for_rank(mesh,stk::topology::ELEM_RANK,scalarFieldElement,vectorFieldElement,meta.locally_owned_part());

    ScalarField & scalarFieldFace =
      static_cast<ScalarField&>(*meta.get_field(stk::topology::FACE_RANK, "Face Scalar Field"));
    VectorField & vectorFieldFace =
      static_cast<VectorField&>(*meta.get_field(stk::topology::FACE_RANK, "Face Vector Field"));

    fill_mesh_values_for_rank(mesh,stk::topology::FACE_RANK,scalarFieldFace,vectorFieldFace,meta.locally_owned_part());

    ScalarField & scalarFieldShell =
      static_cast<ScalarField&>(*meta.get_field(stk::topology::ELEM_RANK, "Shell Scalar Field"));
    VectorField & vectorFieldShell =
      static_cast<VectorField&>(*meta.get_field(stk::topology::ELEM_RANK, "Shell Vector Field"));

    stk::mesh::Part & shell_part = *meta.get_part("shell_part");
    fill_mesh_values_for_rank(mesh,stk::topology::ELEM_RANK,
                              scalarFieldShell,
                              vectorFieldShell,
                              meta.locally_owned_part() & shell_part);
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

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

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Field"));
    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    {
      // Unit test for CopySearchCommAll,
      // also verifies do_search can be called twice
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor::const_iterator map_iter=key_to_target_processor.begin();
      for ( ; map_iter != key_to_target_processor.end()  ; ++map_iter) {
        const int expected_target_processor = 0;
        EXPECT_EQ( expected_target_processor, map_iter->second );
      }

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      const MeshIDSet & remote_keys = copySearch.get_remote_keys();
      EXPECT_TRUE( remote_keys.empty() );
    }

    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    transfer.apply();

    // Check "target" fields to make sure they hold the expected values
    //
    const double tolerance = 1.e-8;
    const stk::mesh::BucketVector & entityBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() );
    for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
      for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
        stk::mesh::Entity entity = entityBucket[entityIndex];
        double * scalarTarget = stk::mesh::field_data(scalarTargetField, entity);
        double * vectorTarget = stk::mesh::field_data(vectorTargetField, entity);

        EXPECT_NEAR(static_cast<double>(meshB.identifier(entity)), *scalarTarget, tolerance);
        for (size_t i = 0; i < spatial_dimension; ++i) {
          EXPECT_NEAR(static_cast<double>((meshB.identifier(entity)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
        }
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;


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

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Field"));
    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    {
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor::const_iterator map_iter=key_to_target_processor.begin();
      for ( ; map_iter != key_to_target_processor.end()  ; ++map_iter) {
        const int expected_target_processor = 1;
        EXPECT_EQ( expected_target_processor, map_iter->second );
      }

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      const MeshIDSet & remote_keys = copySearch.get_remote_keys();
      const int p_rank = stk::parallel_machine_rank( pm );
      if (p_rank == 0) {
        EXPECT_TRUE( remote_keys.empty() );
      } else {
        EXPECT_EQ( 8u, remote_keys.size() );
      }
    }
    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    transfer.apply();

    // Check "target" fields to make sure they hold the expected values
    //
    const double tolerance = 1.e-8;
    const stk::mesh::BucketVector & entityBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() );
    for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
      for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
        stk::mesh::Entity entity = entityBucket[entityIndex];
        double * scalarTarget = stk::mesh::field_data(scalarTargetField, entity);
        double * vectorTarget = stk::mesh::field_data(vectorTargetField, entity);

        EXPECT_NEAR(static_cast<double>(meshB.identifier(entity)), *scalarTarget, tolerance);
        for (size_t i = 0; i < spatial_dimension; ++i) {
          EXPECT_NEAR(static_cast<double>((meshB.identifier(entity)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
        }
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

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

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Field"));
    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    {
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor::const_iterator map_iter=key_to_target_processor.begin();
      for ( ; map_iter != key_to_target_processor.end()  ; ++map_iter) {
        const int expected_target_processor = 0;
        EXPECT_EQ( expected_target_processor, map_iter->second );
      }

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      const MeshIDSet & remote_keys = copySearch.get_remote_keys();
      const int p_rank = stk::parallel_machine_rank( pm );
      if (p_rank == 0) {
        EXPECT_EQ( 8u, remote_keys.size() );
      } else {
        EXPECT_TRUE( remote_keys.empty() );
      }
    }
    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);


    // Do the transfer
    //
    transfer.initialize();
    transfer.apply();

    // Check "target" fields to make sure they hold the expected values
    //
    const double tolerance = 1.e-8;
    const stk::mesh::BucketVector & entityBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() );
    for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
      for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
        stk::mesh::Entity entity = entityBucket[entityIndex];
        double * scalarTarget = stk::mesh::field_data(scalarTargetField, entity);
        double * vectorTarget = stk::mesh::field_data(vectorTargetField, entity);

        EXPECT_NEAR(static_cast<double>(meshB.identifier(entity)), *scalarTarget, tolerance);
        for (size_t i = 0; i < spatial_dimension; ++i) {
          EXPECT_NEAR(static_cast<double>((meshB.identifier(entity)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
        }
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

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

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Field"));
    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    //  GeometricTransfer
    //  stk::transfer::GeometricTransfer<
    //    class stk::transfer::LinearInterpolate<
    //      class stk::transfer::STKNode,
    //      class stk::transfer::STKNode
    //    >
    //  > transfer(transferSource, transferTarget, "copy01T10 unit test");

    {
      const int p_rank = stk::parallel_machine_rank( pm );
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold_map;
      if (0 == p_rank) {
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,1)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,4)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,7)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,10)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,2)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,5)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,8)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,11)] = 0;
      } else {
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,3)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,6)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,9)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,12)] = 0;
      }
      EXPECT_TRUE( gold_map == key_to_target_processor );

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      MeshIDSet gold_remote_keys;
      if (0 == p_rank) {
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,3).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,6).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,9).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,12).m_value);
      } else {
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,1).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,4).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,7).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,10).m_value);
      }
      EXPECT_TRUE( copySearch.get_remote_keys() == gold_remote_keys );
    }
    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    transfer.apply();

    // Check "target" fields to make sure they hold the expected values.
    //
    const double tolerance = 1.e-8;
    const stk::mesh::BucketVector & entityBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part());
    for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
      for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
        stk::mesh::Entity entity = entityBucket[entityIndex];
        double * scalarTarget = stk::mesh::field_data(scalarTargetField, entity);
        double * vectorTarget = stk::mesh::field_data(vectorTargetField, entity);

        EXPECT_NEAR(static_cast<double>(meshB.identifier(entity)), *scalarTarget, tolerance);
        for (size_t i = 0; i < spatial_dimension; ++i) {
          EXPECT_NEAR(static_cast<double>((meshB.identifier(entity)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
        }
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

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

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Field"));
    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    //  GeometricTransfer
    //  stk::transfer::GeometricTransfer<
    //    class stk::transfer::LinearInterpolate<
    //      class stk::transfer::STKNode,
    //      class stk::transfer::STKNode
    //    >
    //  > transfer(transferSource, transferTarget, "copy001T011 unit test");

    {
      const int p_rank = stk::parallel_machine_rank( pm );
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold_map;
      if (0 == p_rank) {
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,1)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,5)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,9)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,13)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,2)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,6)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,10)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,14)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,3)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,7)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,11)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,15)] = 1;
      } else {
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,4)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,8)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,12)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,16)] = 1;
      }
      EXPECT_TRUE( gold_map == key_to_target_processor );

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      MeshIDSet gold_remote_keys;
      if (0 == p_rank) {
      } else {
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,3).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,7).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,11).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,15).m_value);
      }
      EXPECT_TRUE( copySearch.get_remote_keys() == gold_remote_keys );
    }
    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    transfer.apply();

    // Check "target" fields to make sure they hold the expected values.
    //
    const double tolerance = 1.e-8;
    const stk::mesh::BucketVector & entityBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part());
    for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
      for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
        stk::mesh::Entity entity = entityBucket[entityIndex];
        double * scalarTarget = stk::mesh::field_data(scalarTargetField, entity);
        double * vectorTarget = stk::mesh::field_data(vectorTargetField, entity);

        EXPECT_NEAR(static_cast<double>(meshB.identifier(entity)), *scalarTarget, tolerance);
        for (size_t i = 0; i < spatial_dimension; ++i) {
          EXPECT_NEAR(static_cast<double>((meshB.identifier(entity)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
        }
      }
    }
  }
}

TEST(Transfer, copy001T011Element)
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

    const size_t spatial_dimension = 3;
    const size_t num_elements = 3;
    const size_t num_nodes = 16;
    stk::mesh::EntityId element_ids[] = {1, 2, 3};
    int element_ownerA[] = {0, 0, 1};
    int element_ownerB[] = {0, 1, 1};
    stk::mesh::EntityId elem_node_ids[][8] = { {1, 2, 6, 5,  9, 10, 14, 13},
                                               {2, 3, 7, 6, 10, 11, 15, 14},
                                               {3, 4, 8, 7, 11, 12, 16, 15} };
    int node_sharingA[] = { -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1,
                            -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1 };
    int node_sharingB[] = { -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1,
                            -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1 };
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

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::ELEM_RANK, "Element Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::ELEM_RANK, "Element Vector Field"));
    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::ELEM_RANK, "Element Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::ELEM_RANK, "Element Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceElements;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::ELEM_RANK), sourceElements);
    std::vector<stk::mesh::Entity> targetElements;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::ELEM_RANK), targetElements);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceElements, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetElements, targetFields);

    {
      const int p_rank = stk::parallel_machine_rank( pm );
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold_map;
      if (0 == p_rank) {
        gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,1)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,2)] = 1;
      } else {
        gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,3)] = 1;
      }
      EXPECT_TRUE( gold_map == key_to_target_processor );

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      MeshIDSet gold_remote_keys;
      if (1 == p_rank) {
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::ELEM_RANK,2).m_value);
      }
      EXPECT_TRUE( copySearch.get_remote_keys() == gold_remote_keys );
    }
    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    transfer.apply();

    // Check "target" fields to make sure they hold the expected values.
    //
    const double tolerance = 1.e-8;
    const stk::mesh::BucketVector & entityBuckets = meshB.get_buckets(stk::topology::ELEM_RANK, metaB.locally_owned_part());
    for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
      for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
        stk::mesh::Entity entity = entityBucket[entityIndex];
        double * scalarTarget = stk::mesh::field_data(scalarTargetField, entity);
        double * vectorTarget = stk::mesh::field_data(vectorTargetField, entity);

        EXPECT_NEAR(static_cast<double>(meshB.identifier(entity)), *scalarTarget, tolerance);
        for (size_t i = 0; i < spatial_dimension; ++i) {
          EXPECT_NEAR(static_cast<double>((meshB.identifier(entity)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
        }
      }
    }
  }
}

TEST(Transfer, copy001T011Face)
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

    const size_t spatial_dimension = 3;
    const size_t num_elements = 3;
    const size_t num_nodes = 16;
    stk::mesh::EntityId element_ids[] = {1, 2, 3};
    int element_ownerA[] = {0, 0, 1};
    int element_ownerB[] = {0, 1, 1};
    stk::mesh::EntityId elem_node_ids[][8] = { { 9,10,2,1,13,14,6,5},
                                               {10,11,3,2,14,15,7,6},
                                               {11,12,4,3,15,16,8,7} };
    int node_sharingA[] = { -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1,
                            -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1 };
    int node_sharingB[] = { -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1,
                            -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1 };
    double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
                                {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0},
                                {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0},
                                {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0} };
    stk::mesh::EntityId face_node_ids[][4] =
    { { 9,10,14,13}, {10,2,6,14}, {2,1,5,6}, { 9,13,5, 1}, { 9,1,2,10}, {13,14,6,5},
      {10,11,15,14}, {11,3,7,15}, {3,2,6,7}, {10, 2,6,14}, {10,2,3,11}, {14,15,7,6},
      {11,12,16,15}, {12,4,8,16}, {4,3,7,8}, {11, 3,7,15}, {11,3,4,12}, {15,16,8,7} };
    stk::mesh::EntityId elem_face_ids[][6] = { { 1,  2,  3,  4,  5,  6},
                                               { 7,  8,  9,  2, 10, 11},
                                               {12, 13, 14,  8, 15, 16} };


    // Set up the "source" mesh for the transfer
    //
    stk::mesh::MetaData metaA(spatial_dimension);
    stk::mesh::BulkData meshA(metaA, pm);
    build_mesh(metaA,
               meshA,
               num_elements,
               num_nodes,
               element_ids,
               element_ownerA,
               elem_node_ids,
               node_sharingA,
               coordinates,
               face_node_ids,
               elem_face_ids);

    // Set up the "target" mesh for the transfer
    //
    stk::mesh::MetaData metaB(spatial_dimension);
    stk::mesh::BulkData meshB(metaB, pm);
    build_mesh(metaB,
               meshB,
               num_elements,
               num_nodes,
               element_ids,
               element_ownerB,
               elem_node_ids,
               node_sharingB,
               coordinates,
               face_node_ids,
               elem_face_ids);

    // Fill "source" fields with valid data
    //
    fill_mesh_values(meshA);

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::FACE_RANK, "Face Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::FACE_RANK, "Face Vector Field"));
    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::FACE_RANK, "Face Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::FACE_RANK, "Face Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceFaces;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::FACE_RANK), sourceFaces);
    std::vector<stk::mesh::Entity> targetFaces;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::FACE_RANK), targetFaces);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceFaces, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetFaces, targetFields);

    {
      const int p_rank = stk::parallel_machine_rank( pm );
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold_map;
      if (0 == p_rank) {
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,1)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,2)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,3)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,4)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,5)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,6)] = 0;

        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,7)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,8)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,9)] = 1;
        //gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,2)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,10)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,11)] = 1;
      } else {
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,12)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,13)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,14)] = 1;
        //gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,8)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,15)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::FACE_RANK,16)] = 1;
      }
      EXPECT_TRUE( gold_map == key_to_target_processor );

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      MeshIDSet gold_remote_keys;
      if (1 == p_rank) {
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::FACE_RANK,7).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::FACE_RANK,8).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::FACE_RANK,9).m_value);
        //gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::FACE_RANK,2).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::FACE_RANK,10).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::FACE_RANK,11).m_value);
      }
      EXPECT_TRUE( copySearch.get_remote_keys() == gold_remote_keys );
    }
    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    transfer.apply();

    // Check "target" fields to make sure they hold the expected values.
    //
    const double tolerance = 1.e-8;
    const stk::mesh::BucketVector & entityBuckets = meshB.get_buckets(stk::topology::FACE_RANK, metaB.locally_owned_part());
    for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
      for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
        stk::mesh::Entity entity = entityBucket[entityIndex];
        double * scalarTarget = stk::mesh::field_data(scalarTargetField, entity);
        double * vectorTarget = stk::mesh::field_data(vectorTargetField, entity);

        EXPECT_NEAR(static_cast<double>(meshB.identifier(entity)), *scalarTarget, tolerance);
        for (size_t i = 0; i < spatial_dimension; ++i) {
          EXPECT_NEAR(static_cast<double>((meshB.identifier(entity)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
        }
      }
    }
  }
}

TEST(Transfer, copy001T011Shell)
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

    const size_t spatial_dimension = 3;
    const size_t num_elements = 3;
    const size_t num_nodes = 16;
    stk::mesh::EntityId element_ids[] = {1, 2, 3};
    int element_ownerA[] = {0, 0, 1};
    int element_ownerB[] = {0, 1, 1};
    stk::mesh::EntityId elem_node_ids[][8] = { { 9,10,2,1,13,14,6,5},
                                               {10,11,3,2,14,15,7,6},
                                               {11,12,4,3,15,16,8,7} };
    int node_sharingA[] = { -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1,
                            -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1 };
    int node_sharingB[] = { -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1,
                            -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1 };
    double coordinates[][3] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
                                {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0},
                                {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0},
                                {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0} };


    // Set up the "source" mesh for the transfer
    //
    stk::mesh::MetaData metaA(spatial_dimension);
    stk::mesh::BulkData meshA(metaA, pm);
    build_mesh(metaA,
               meshA,
               num_elements,
               num_nodes,
               element_ids,
               element_ownerA,
               elem_node_ids,
               node_sharingA,
               coordinates);

    // Set up the "target" mesh for the transfer
    //
    stk::mesh::MetaData metaB(spatial_dimension);
    stk::mesh::BulkData meshB(metaB, pm);
    build_mesh(metaB,
               meshB,
               num_elements,
               num_nodes,
               element_ids,
               element_ownerB,
               elem_node_ids,
               node_sharingB,
               coordinates);


    stk::mesh::EntityId shell_node_ids[][4] =
    { { 9,10,14,13}, {10,2,6,14}, {2,1,5,6}, { 9,13,5, 1}, { 9,1,2,10}, {13,14,6,5},
      {10,11,15,14}, {11,3,7,15}, {3,2,6,7}, {10, 2,6,14}, {10,2,3,11}, {14,15,7,6},
      {11,12,16,15}, {12,4,8,16}, {4,3,7,8}, {11, 3,7,15}, {11,3,4,12}, {15,16,8,7} };
    stk::mesh::EntityId elem_shell_ids[][6] = { {10, 11, 12, 13, 14, 15},
                                                {16, 17, 18, 11, 19, 20},
                                                {21, 22, 23, 17, 24, 25} };
    int shell_owner_by_elem_sideA[][6] = { {0, 0, 0, 0, 0, 0},
                                           {0, 0, 0, 0, 0, 0},
                                           {1, 1, 1, 0, 1, 1} };
    int shell_owner_by_elem_sideB[][6] = { {0, 0, 0, 0, 0, 0},
                                           {1, 1, 1, 0, 1, 1},
                                           {1, 1, 1, 1, 1, 1} };
    const int num_shells = 10;
    add_shells_to_mesh(metaA,
                       meshA,
                       num_elements,
                       element_ids,
                       shell_node_ids,
                       elem_shell_ids,
                       shell_owner_by_elem_sideA,
                       num_shells);
    add_shells_to_mesh(metaB,
                       meshB,
                       num_elements,
                       element_ids,
                       shell_node_ids,
                       elem_shell_ids,
                       shell_owner_by_elem_sideB,
                       num_shells);

    // Fill "source" fields with valid data
    //
    fill_mesh_values(meshA);

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::ELEM_RANK, "Shell Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::ELEM_RANK, "Shell Vector Field"));
    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::ELEM_RANK, "Shell Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::ELEM_RANK, "Shell Vector Field"));

    // Set up the transfer
    //
    stk::mesh::Part & shell_part = *metaA.get_part("shell_part");

    std::vector<stk::mesh::Entity> sourceShells;
    stk::mesh::get_selected_entities(metaA.locally_owned_part() & shell_part,
                                     meshA.buckets(stk::topology::ELEM_RANK),
                                     sourceShells);

    std::vector<stk::mesh::Entity> targetShells;
    stk::mesh::get_selected_entities(metaB.locally_owned_part() & shell_part,
                                     meshB.buckets(stk::topology::ELEM_RANK),
                                     targetShells);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceShells, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetShells, targetFields);

    {
      const int p_rank = stk::parallel_machine_rank( pm );
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold_map;
      for (int shell_index=0; shell_index<num_shells; ++shell_index) {
        if (0 == p_rank) {
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,10+100*shell_index)] = 0;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,11+100*shell_index)] = 0;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,12+100*shell_index)] = 0;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,13+100*shell_index)] = 0;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,14+100*shell_index)] = 0;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,15+100*shell_index)] = 0;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,16+100*shell_index)] = 1;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,17+100*shell_index)] = 1;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,18+100*shell_index)] = 1;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,19+100*shell_index)] = 1;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,20+100*shell_index)] = 1;
        } else {
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,21+100*shell_index)] = 1;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,22+100*shell_index)] = 1;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,23+100*shell_index)] = 1;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,24+100*shell_index)] = 1;
          gold_map[stk::mesh::EntityKey(stk::topology::ELEM_RANK,25+100*shell_index)] = 1;
        }
      }
      EXPECT_TRUE( gold_map == key_to_target_processor );

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      MeshIDSet gold_remote_keys;
      for (int shell_index=0; shell_index<num_shells; ++shell_index) {
        if (1 == p_rank) {
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::ELEM_RANK,16+100*shell_index).m_value);
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::ELEM_RANK,17+100*shell_index).m_value);
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::ELEM_RANK,18+100*shell_index).m_value);
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::ELEM_RANK,19+100*shell_index).m_value);
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::ELEM_RANK,20+100*shell_index).m_value);
        }
      }
      EXPECT_TRUE( copySearch.get_remote_keys() == gold_remote_keys );
    }
    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    transfer.apply();

    // Check "target" fields to make sure they hold the expected values.
    //
    const double tolerance = 1.e-8;
    const stk::mesh::BucketVector & entityBuckets =
      meshB.get_buckets(stk::topology::ELEM_RANK, metaB.locally_owned_part() & shell_part);
    for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
      for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
        stk::mesh::Entity entity = entityBucket[entityIndex];
        double * scalarTarget = stk::mesh::field_data(scalarTargetField, entity);
        double * vectorTarget = stk::mesh::field_data(vectorTargetField, entity);

        EXPECT_NEAR(static_cast<double>(meshB.identifier(entity)), *scalarTarget, tolerance);
        for (size_t i = 0; i < spatial_dimension; ++i) {
          EXPECT_NEAR(static_cast<double>((meshB.identifier(entity)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
        }
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

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

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Field"));
    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    //  GeometricTransfer
    //  stk::transfer::GeometricTransfer<
    //    class stk::transfer::LinearInterpolate<
    //      class stk::transfer::STKNode,
    //      class stk::transfer::STKNode
    //    >
    //  > transfer(transferSource, transferTarget, "copy012T000 unit test");

    {
      const int p_rank = stk::parallel_machine_rank( pm );
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold_map;
      if (0 == p_rank) {
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,1)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,5)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,9)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,13)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,2)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,6)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,10)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,14)] = 0;
      } else if (1 == p_rank) {
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,3)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,7)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,11)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,15)] = 0;
      } else {
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,4)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,8)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,12)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,16)] = 0;
      }
      EXPECT_TRUE( gold_map == key_to_target_processor );

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      MeshIDSet gold_remote_keys;
      if (0 == p_rank) {
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,3).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,7).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,11).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,15).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,4).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,8).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,12).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,16).m_value);
      }
      EXPECT_TRUE( copySearch.get_remote_keys() == gold_remote_keys );
    }
    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    transfer.apply();

    // Check "target" fields to make sure they hold the expected values.
    //
    const double tolerance = 1.e-8;
    const stk::mesh::BucketVector & entityBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() );
    for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
      for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
        stk::mesh::Entity entity = entityBucket[entityIndex];
        double * scalarTarget = stk::mesh::field_data(scalarTargetField, entity);
        double * vectorTarget = stk::mesh::field_data(vectorTargetField, entity);

        EXPECT_NEAR(static_cast<double>(meshB.identifier(entity)), *scalarTarget, tolerance);
        for (size_t i = 0; i < spatial_dimension; ++i) {
          EXPECT_NEAR(static_cast<double>((meshB.identifier(entity)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
        }
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

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

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Field"));
    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    //  GeometricTransfer
    //  stk::transfer::GeometricTransfer<
    //    class stk::transfer::LinearInterpolate<
    //      class stk::transfer::STKNode,
    //      class stk::transfer::STKNode
    //    >
    //  > transfer(transferSource, transferTarget, "copy000T012 unit test");

    {
      const int p_rank = stk::parallel_machine_rank( pm );
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold_map;
      if (0 == p_rank) {
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,1)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,5)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,9)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,13)] = 0;

        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,2)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,6)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,10)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,14)] = 0;

        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,3)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,7)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,11)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,15)] = 1;

        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,4)] = 2;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,8)] = 2;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,12)] = 2;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,16)] = 2;
      }
      EXPECT_TRUE( gold_map == key_to_target_processor );

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      MeshIDSet gold_remote_keys;
      if (1 == p_rank) {
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,3).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,7).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,11).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,15).m_value);
      } else if (2 == p_rank) {
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,4).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,8).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,12).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,16).m_value);
      }
      EXPECT_TRUE( copySearch.get_remote_keys() == gold_remote_keys );
    }

    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    transfer.apply();

    // Check "target" fields to make sure they hold the expected values.
    //
    const double tolerance = 1.e-8;
    const stk::mesh::BucketVector & entityBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() );
    for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
      for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
        stk::mesh::Entity entity = entityBucket[entityIndex];
        double * scalarTarget = stk::mesh::field_data(scalarTargetField, entity);
        double * vectorTarget = stk::mesh::field_data(vectorTargetField, entity);

        EXPECT_NEAR(static_cast<double>(meshB.identifier(entity)), *scalarTarget, tolerance);
        for (size_t i = 0; i < spatial_dimension; ++i) {
          EXPECT_NEAR(static_cast<double>((meshB.identifier(entity)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
        }
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

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

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Field"));
    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    //  GeometricTransfer
    //  stk::transfer::GeometricTransfer<
    //    class stk::transfer::LinearInterpolate<
    //      class stk::transfer::STKNode,
    //      class stk::transfer::STKNode
    //    >
    //  > transfer(transferSource, transferTarget, "copy0011T1010 unit test");

    {
      const int p_rank = stk::parallel_machine_rank( pm );
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold_map;
      if (0 == p_rank) {
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,1)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,6)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,11)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,16)] = 1;

        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,2)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,7)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,12)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,17)] = 0;

        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,3)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,8)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,13)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,18)] = 0;
      } else if (1 == p_rank) {
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,4)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,9)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,14)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,19)] = 0;

        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,5)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,10)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,15)] = 0;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,20)] = 0;
      }
      EXPECT_TRUE( gold_map == key_to_target_processor );

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      MeshIDSet gold_remote_keys;
      if (0 == p_rank) {
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,4).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,9).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,14).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,19).m_value);

        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,5).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,10).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,15).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,20).m_value);
      } else {
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,1).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,6).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,11).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,16).m_value);
      }
      EXPECT_TRUE( copySearch.get_remote_keys() == gold_remote_keys );
    }
    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    transfer.apply();

    // Check "target" fields to make sure they hold the expected values.
    //
    const double tolerance = 1.e-8;
    const stk::mesh::BucketVector & entityBuckets = meshB.get_buckets(stk::topology::NODE_RANK, metaB.locally_owned_part() );
    for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
      stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
      for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
        stk::mesh::Entity entity = entityBucket[entityIndex];
        double * scalarTarget = stk::mesh::field_data(scalarTargetField, entity);
        double * vectorTarget = stk::mesh::field_data(vectorTargetField, entity);

        EXPECT_NEAR(static_cast<double>(meshB.identifier(entity)), *scalarTarget, tolerance);
        for (size_t i = 0; i < spatial_dimension; ++i) {
          EXPECT_NEAR(static_cast<double>((meshB.identifier(entity)-1)*spatial_dimension+i), vectorTarget[i], tolerance);
        }
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

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
    ScalarField & scalarTargetField = metaB.declare_field<ScalarField>(stk::topology::NODE_RANK, "Node Scalar Field");
    VectorField & vectorTargetField = metaB.declare_field<VectorField>(stk::topology::NODE_RANK, "Node Vector Field");
    VectorField & coordsTargetField = metaB.declare_field<VectorField>(stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field(scalarTargetField, metaB.universal_part(), init_vals);
    stk::mesh::put_field(vectorTargetField, metaB.universal_part(), init_vals);
    stk::mesh::put_field(coordsTargetField, metaB.universal_part(), init_vals);
    metaB.commit();

    // Fill "source" fields with valid data
    //
    fill_mesh_values(meshA);

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    //  GeometricTransfer
    //  stk::transfer::GeometricTransfer<
    //    class stk::transfer::LinearInterpolate<
    //      class stk::transfer::STKNode,
    //      class stk::transfer::STKNode
    //    >
    //  > transfer(transferSource, transferTarget, "copy0T_ unit test");

    {
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold_map;
      EXPECT_TRUE( gold_map == key_to_target_processor );

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      MeshIDSet gold_remote_keys;
      EXPECT_TRUE( copySearch.get_remote_keys() == gold_remote_keys );
    }
    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    EXPECT_NO_THROW(transfer.initialize());
    EXPECT_NO_THROW(transfer.apply());
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

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
    ScalarField & scalarSourceField = metaA.declare_field<ScalarField>(stk::topology::NODE_RANK, "Node Scalar Field");
    VectorField & vectorSourceField = metaA.declare_field<VectorField>(stk::topology::NODE_RANK, "Node Vector Field");
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

    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    {
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold_map;
      EXPECT_TRUE( gold_map == key_to_target_processor );

      EXPECT_EQ( 8u, copySearch.get_remote_keys().size() );
    }
    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    EXPECT_THROW(transfer.apply(), std::runtime_error);
  }
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

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

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Field"));
    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    //  GeometricTransfer
    //  stk::transfer::GeometricTransfer<
    //    class stk::transfer::LinearInterpolate<
    //      class stk::transfer::STKNode,
    //      class stk::transfer::STKNode
    //    >
    //  > transfer(transferSource, transferTarget, "copy01T10 unit test");

    {
      const int p_rank = stk::parallel_machine_rank( pm );
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold_map;
      if (0 == p_rank) {
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,2)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,6)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,10)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,14)] = 1;

        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,3)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,7)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,11)] = 1;
        gold_map[stk::mesh::EntityKey(stk::topology::NODE_RANK,15)] = 1;
      }
      EXPECT_TRUE( gold_map == key_to_target_processor );

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      MeshIDSet gold_remote_keys;
      if (1 == p_rank) {
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,2).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,6).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,10).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,14).m_value);

        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,3).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,7).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,11).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,15).m_value);

        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,4).m_value); // these don't get filled.
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,8).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,12).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,16).m_value);
      }
      EXPECT_TRUE( copySearch.get_remote_keys() == gold_remote_keys );
    }
    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    EXPECT_THROW(transfer.apply(), std::runtime_error);
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
//            |  |  1.0  |  |  2.0  |  |  -->   |  |  4.1  |  |  5.1  |  |
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

  stk::transfer::CopySearchGeometric geometricSearch;
  stk::transfer::CopySearchCommAll commAllSearch;
  stk::transfer::CopySearchBase * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_TRUE( copySearchPtr == &geometricSearch );
    }
    stk::transfer::CopySearchBase & copySearch = *copySearchPtr;

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

    ScalarField & scalarSourceField = static_cast<ScalarField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorSourceField = static_cast<VectorField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Field"));
    ScalarField & scalarTargetField = static_cast<ScalarField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Field"));
    VectorField & vectorTargetField = static_cast<VectorField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(), meshA.buckets(stk::topology::NODE_RANK), sourceNodes);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(), meshB.buckets(stk::topology::NODE_RANK), targetNodes);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceField);
    sourceFields.push_back(&vectorSourceField);
    stk::transfer::CopyTransferStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetField);
    targetFields.push_back(&vectorTargetField);
    stk::transfer::CopyTransferStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    //  GeometricTransfer
    //  stk::transfer::GeometricTransfer<
    //    class stk::transfer::LinearInterpolate<
    //      class stk::transfer::STKNode,
    //      class stk::transfer::STKNode
    //    >
    //  > transfer(transferSource, transferTarget, "copy01T10 unit test");

    {
      const int p_rank = stk::parallel_machine_rank( pm );
      typedef stk::transfer::CopySearchBase::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold_map;
      EXPECT_TRUE( gold_map == key_to_target_processor );

      typedef stk::transfer::CopySearchBase::MeshIDSet MeshIDSet;
      MeshIDSet gold_remote_keys;
      if (1 == p_rank) {
        // none of these get fulfilled
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,4).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,10).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,16).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,22).m_value);

        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,5).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,11).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,17).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,23).m_value);

        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,6).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,12).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,18).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,24).m_value);
      }
      EXPECT_TRUE( copySearch.get_remote_keys() == gold_remote_keys );
    }
    stk::transfer::CopyTransfer transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    EXPECT_THROW(transfer.apply(), std::runtime_error);

  }
}

} // namespace

