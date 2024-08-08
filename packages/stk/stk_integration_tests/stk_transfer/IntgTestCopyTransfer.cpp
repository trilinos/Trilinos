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

#include "stk_mesh/base/BulkData.hpp"          // for BulkData, etc
#include "stk_mesh/base/FEMHelpers.hpp"        // for declare_element
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MetaData.hpp"          // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"              // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_topology/topology.hpp"           // for topology, etc
#include "stk_util/parallel/Parallel.hpp"      // for ParallelMachine, etc
#include "stk_unit_test_utils/TextMesh.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_io/IossBridge.hpp"
#include "gtest/gtest.h"
#include <limits>
#include <memory>
#include <string>

#include "stk_transfer/copy_by_id/SearchByIdCommAll.hpp"
#include "stk_transfer/copy_by_id/SearchByIdGeometric.hpp"
#include "stk_transfer/copy_by_id/TransferCopyById.hpp"
#include "stk_transfer/copy_by_id/TransferCopyByIdMpmdMeshAdapter.hpp"
#include "stk_transfer/copy_by_id/TransferCopyByIdStkMeshAdapter.hpp"
#include "stk_transfer/copy_by_id/TransferCopyTranslator.hpp"

namespace
{
using stk::unit_test_util::build_mesh;

typedef stk::mesh::Field<int>    ScalarIntField;
typedef stk::mesh::Field<double> ScalarDoubleField;
typedef stk::mesh::Field<double> VectorDoubleField;

void build_mesh(stk::mesh::MetaData & meta,
                stk::mesh::BulkData & mesh,
                const size_t num_elements,
                const size_t num_nodes,
                const stk::mesh::EntityIdVector & element_ids,
                std::vector<int> element_owner,
                const stk::mesh::EntityIdVector * elem_node_ids,
                std::vector<int> node_sharing,
                std::vector<std::vector<double>> coordinates,
                bool createFaces = false )
{
  const int p_rank = mesh.parallel_rank();
  int int_init_vals = std::numeric_limits<int>::max();
  double double_init_vals[] = {std::numeric_limits<double>::max(),
                               std::numeric_limits<double>::max(),
                               std::numeric_limits<double>::max()};

  stk::mesh::Part * elem_part = &meta.declare_part_with_topology("elem_part", stk::topology::HEX_8);
  stk::mesh::Part * face_part = &meta.declare_part_with_topology("face_part", stk::topology::QUAD_4);
  stk::mesh::Part * shell_part = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  ScalarIntField    & scalarIntFieldNode    = meta.declare_field<int>(stk::topology::NODE_RANK, "Node Scalar Int Field");
  ScalarDoubleField & scalarDoubleFieldNode = meta.declare_field<double>(stk::topology::NODE_RANK, "Node Scalar Double Field");
  VectorDoubleField & vectorDoubleFieldNode = meta.declare_field<double>(stk::topology::NODE_RANK, "Node Vector Double Field");
  VectorDoubleField & coordsFieldNode       = meta.declare_field<double>(stk::topology::NODE_RANK, "coordinates");
  meta.set_coordinate_field(&coordsFieldNode);
  stk::mesh::put_field_on_mesh(scalarIntFieldNode, meta.universal_part(), &int_init_vals);
  stk::mesh::put_field_on_mesh(scalarDoubleFieldNode, meta.universal_part(), double_init_vals);
  stk::mesh::put_field_on_mesh(vectorDoubleFieldNode, meta.universal_part(), 3, double_init_vals);
  stk::mesh::put_field_on_mesh(coordsFieldNode, meta.universal_part(), 3, double_init_vals);
  ScalarIntField    & scalarIntFieldElement    = meta.declare_field<int>(stk::topology::ELEM_RANK, "Element Scalar Int Field");
  ScalarDoubleField & scalarDoubleFieldElement = meta.declare_field<double>(stk::topology::ELEM_RANK, "Element Scalar Double Field");
  VectorDoubleField & vectorDoubleFieldElement = meta.declare_field<double>(stk::topology::ELEM_RANK, "Element Vector Double Field");
  stk::mesh::put_field_on_mesh(scalarIntFieldElement, meta.universal_part(), &int_init_vals);
  stk::mesh::put_field_on_mesh(scalarDoubleFieldElement, meta.universal_part(), double_init_vals);
  stk::mesh::put_field_on_mesh(vectorDoubleFieldElement, meta.universal_part(), 3, double_init_vals);
  ScalarIntField    & scalarIntFieldFace    = meta.declare_field<int>(stk::topology::FACE_RANK, "Face Scalar Int Field");
  ScalarDoubleField & scalarDoubleFieldFace = meta.declare_field<double>(stk::topology::FACE_RANK, "Face Scalar Double Field");
  VectorDoubleField & vectorDoubleFieldFace = meta.declare_field<double>(stk::topology::FACE_RANK, "Face Vector Double Field");
  stk::mesh::put_field_on_mesh(scalarIntFieldFace, meta.universal_part(), &int_init_vals);
  stk::mesh::put_field_on_mesh(scalarDoubleFieldFace, meta.universal_part(), double_init_vals);
  stk::mesh::put_field_on_mesh(vectorDoubleFieldFace, meta.universal_part(), 3, double_init_vals);
  ScalarIntField    & scalarIntFieldShell    = meta.declare_field<int>(stk::topology::ELEM_RANK, "Shell Scalar Int Field");
  ScalarDoubleField & scalarDoubleFieldShell = meta.declare_field<double>(stk::topology::ELEM_RANK, "Shell Scalar Double Field");
  VectorDoubleField & vectorDoubleFieldShell = meta.declare_field<double>(stk::topology::ELEM_RANK, "Shell Vector Double Field");
  stk::mesh::put_field_on_mesh(scalarIntFieldShell, *shell_part, &int_init_vals);
  stk::mesh::put_field_on_mesh(scalarDoubleFieldShell, *shell_part, double_init_vals);
  stk::mesh::put_field_on_mesh(vectorDoubleFieldShell, *shell_part, 3, double_init_vals);
  meta.commit();

  mesh.initialize_face_adjacent_element_graph();
  mesh.modification_begin();
  for (size_t i = 0; i < num_elements; ++i) {
    if (p_rank == element_owner[i]) {
      stk::mesh::declare_element(mesh, *elem_part, element_ids[i], elem_node_ids[i]);
    }
  }
  for (size_t i = 0; i < num_nodes; ++i) {
    if (node_sharing[p_rank*num_nodes+i] != -1) {
      stk::mesh::Entity entity = mesh.get_entity(stk::topology::NODE_RANK, i+1);
      mesh.add_node_sharing(entity, node_sharing[p_rank*num_nodes+i]);
    }
  }
  mesh.modification_end();

  if (createFaces) {
    mesh.modification_begin();
    stk::mesh::PartVector add_parts;
    add_parts.push_back(face_part);
    for (size_t i = 0; i < num_elements; ++i) {
      if (p_rank == element_owner[i]) {
        stk::mesh::Entity element = mesh.get_entity(stk::topology::ELEM_RANK,element_ids[i]);
        for (size_t side_id = 0; side_id < 6; ++side_id) {
          mesh.declare_element_side(element, side_id, add_parts);
        }
      }
    }
    mesh.modification_end();
  }

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
                        const stk::mesh::EntityIdVector & element_ids,
                        const stk::mesh::EntityIdVector shell_node_ids[],
                        const stk::mesh::EntityIdVector elem_shell_ids[],
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
            shell = mesh.declare_element(shell_global_id+100*shell_index, add_parts);
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
                               ScalarIntField & scalarIntField,
                               ScalarDoubleField & scalarDoubleField,
                               VectorDoubleField & vectorDoubleField,
                               const stk::mesh::Selector & sel)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  const unsigned spatial_dimension = meta.spatial_dimension();

  const stk::mesh::BucketVector & entityBuckets = mesh.get_buckets(rank,sel);
  for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
    for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
      stk::mesh::Entity entity = entityBucket[entityIndex];
      int    * scalarInt    = stk::mesh::field_data(scalarIntField, entity);
      double * scalarDouble = stk::mesh::field_data(scalarDoubleField, entity);
      double * vectorDouble = stk::mesh::field_data(vectorDoubleField, entity);

      *scalarInt = mesh.identifier(entity);
      *scalarDouble = static_cast<double>(mesh.identifier(entity));
      for (unsigned i = 0; i < spatial_dimension; ++i) {
        vectorDouble[i] = static_cast<double>((mesh.identifier(entity)-1)*spatial_dimension+i);
      }
    }
  }
}

void fill_mesh_values(stk::mesh::BulkData & mesh)
{
    const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    ScalarIntField    & scalarIntField    = static_cast<ScalarIntField&>(*meta.get_field(stk::topology::NODE_RANK, "Node Scalar Int Field"));
    ScalarDoubleField & scalarDoubleField = static_cast<ScalarDoubleField&>(*meta.get_field(stk::topology::NODE_RANK, "Node Scalar Double Field"));
    VectorDoubleField & vectorDoubleField = static_cast<VectorDoubleField&>(*meta.get_field(stk::topology::NODE_RANK, "Node Vector Double Field"));

    fill_mesh_values_for_rank(mesh,stk::topology::NODE_RANK,scalarIntField,scalarDoubleField,vectorDoubleField,meta.locally_owned_part());

    ScalarIntField    & scalarIntFieldElement    = static_cast<ScalarIntField&>(*meta.get_field(stk::topology::ELEM_RANK, "Element Scalar Int Field"));
    ScalarDoubleField & scalarDoubleFieldElement = static_cast<ScalarDoubleField&>(*meta.get_field(stk::topology::ELEM_RANK, "Element Scalar Double Field"));
    VectorDoubleField & vectorDoubleFieldElement = static_cast<VectorDoubleField&>(*meta.get_field(stk::topology::ELEM_RANK, "Element Vector Double Field"));

    fill_mesh_values_for_rank(mesh,stk::topology::ELEM_RANK,scalarIntFieldElement,scalarDoubleFieldElement,vectorDoubleFieldElement,meta.locally_owned_part());

    ScalarIntField    & scalarIntFieldFace    = static_cast<ScalarIntField&>(*meta.get_field(stk::topology::FACE_RANK, "Face Scalar Int Field"));
    ScalarDoubleField & scalarDoubleFieldFace = static_cast<ScalarDoubleField&>(*meta.get_field(stk::topology::FACE_RANK, "Face Scalar Double Field"));
    VectorDoubleField & vectorDoubleFieldFace = static_cast<VectorDoubleField&>(*meta.get_field(stk::topology::FACE_RANK, "Face Vector Double Field"));

    fill_mesh_values_for_rank(mesh,stk::topology::FACE_RANK,scalarIntFieldFace,scalarDoubleFieldFace,vectorDoubleFieldFace,meta.locally_owned_part());

    ScalarIntField    & scalarIntFieldShell    = static_cast<ScalarIntField&>(*meta.get_field(stk::topology::ELEM_RANK, "Shell Scalar Int Field"));
    ScalarDoubleField & scalarDoubleFieldShell = static_cast<ScalarDoubleField&>(*meta.get_field(stk::topology::ELEM_RANK, "Shell Scalar Double Field"));
    VectorDoubleField & vectorDoubleFieldShell = static_cast<VectorDoubleField&>(*meta.get_field(stk::topology::ELEM_RANK, "Shell Vector Double Field"));

    stk::mesh::Part & shell_part = *meta.get_part("shell_part");
    fill_mesh_values_for_rank(mesh,stk::topology::ELEM_RANK,scalarIntFieldShell,scalarDoubleFieldShell,vectorDoubleFieldShell,meta.locally_owned_part() & shell_part);
}

struct TwoElemMeshInfo
{
    const size_t num_elements = 1;
    const size_t num_nodes = 8;
    stk::mesh::EntityIdVector element_ids {1};
    std::vector<stk::mesh::EntityIdVector> elem_node_ids { {1, 2, 3, 4, 5, 6, 7, 8} };
    std::vector<int> node_sharingA = { -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1 };
    std::vector<int> node_sharingB = { -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1 };
    std::vector<std::vector<double>> coordinates = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0} };
};

struct FourElemMeshInfo
{
    const size_t spatial_dimension = 3;
    const size_t num_elements = 2;
    const size_t num_nodes = 12;
    stk::mesh::EntityIdVector element_ids {1, 2};
    std::vector<stk::mesh::EntityIdVector> elem_node_ids {
        {1, 2, 5, 4, 7, 8, 11, 10},
        {2, 3, 6, 5, 8, 9, 12, 11}
    };
    std::vector<int> node_sharingA = { -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1,
      -1, 0, -1, -1, 0, -1, -1, 0, -1, -1, 0, -1 };
    std::vector<int> node_sharingB = { -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1,
      -1, 0, -1, -1, 0, -1, -1, 0, -1, -1, 0, -1 };
    std::vector<std::vector<double>> coordinates = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
      {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0},
      {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0} };
};

struct SixElemMeshInfo
{
  const size_t spatial_dimension = 3;
  const size_t num_elements = 3;
  const size_t num_nodes = 16;
  stk::mesh::EntityIdVector element_ids {1, 2, 3};
  std::vector<stk::mesh::EntityIdVector> elem_node_ids {
      {1, 2, 6, 5, 9, 10, 14, 13},
      {2, 3, 7, 6, 10, 11, 15, 14},
      {3, 4, 8, 7, 11, 12, 16, 15} };
  std::vector<int> node_sharingA = { -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1,
      -1, -1, 0, -1, -1, -1, 0, -1, -1, -1, 0, -1, -1, -1, 0, -1 };
  std::vector<int> node_sharingB = { -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1,
      -1, 0, -1, -1, -1, 0, -1, -1, -1, 0, -1, -1, -1, 0, -1, -1 };
  std::vector<std::vector<double>> coordinates = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
      {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0},
      {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0} };
};

struct TwoElemKeyToTargetTest
{
  explicit TwoElemKeyToTargetTest(int expected_target_processor) : m_expected_target_processor(expected_target_processor) {}
  void operator()(const stk::transfer::SearchById::KeyToTargetProcessor & key_to_target) const
  {
    for (auto && elem : key_to_target) {
      EXPECT_EQ( m_expected_target_processor, elem.second );
    }
  }
  int m_expected_target_processor;
};

namespace {
std::string field_prefix(stk::mesh::EntityRank rank)
{
  switch (rank)
  {
  case stk::topology::NODE_RANK:
    return "Node";
  case stk::topology::ELEM_RANK:
    return "Element";
  case stk::topology::FACE_RANK:
    return "Face";
  default:
    throw std::runtime_error("Unknown field prefix");
  }
}
}

using SearchByIdTypes = ::testing::Types<stk::transfer::SearchByIdCommAll, stk::transfer::SearchByIdGeometric>;

#ifdef TYPED_TEST_SUITE
  TYPED_TEST_SUITE(CopyTransferFixture, SearchByIdTypes);
#else
  TYPED_TEST_CASE(CopyTransferFixture, SearchByIdTypes);
#endif

template <class SearchById>
class CopyTransferFixture : public ::testing::Test
{
public:
  typedef stk::transfer::SearchById::KeyToTargetProcessor KeyToTargetProcessor;
  void init(stk::ParallelMachine global_pm, int color, std::vector<int> mesh_color_ownership = {0, 0})
  {
    pm = global_pm;
    MPI_Comm_split(global_pm, color, stk::parallel_machine_rank(global_pm), &pmSub);
    commOwnsMesh.resize(2);
    commOwnsMesh[0] = color == mesh_color_ownership[0];
    commOwnsMesh[1] = color == mesh_color_ownership[1];
  }

  template <class MeshInfo>
  void build_fixture(std::vector<int> element_ownerA, std::vector<int> element_ownerB, MeshInfo info, bool create_faces = false)
  {
    // Set up the "source" mesh for the transfer
    //
    if (commOwnsMesh[0])
    {
      metaA = stk::mesh::MeshBuilder().set_spatial_dimension(spatial_dimension).create_meta_data();
      meshA = stk::mesh::MeshBuilder(pmSub).create(metaA);
      build_mesh(*metaA, *meshA, info.num_elements, info.num_nodes, info.element_ids, element_ownerA, &info.elem_node_ids[0], info.node_sharingA, info.coordinates, create_faces);
    }

    // Set up the "target" mesh for the transfer
    //
    if (commOwnsMesh[1])
    {
      metaB = stk::mesh::MeshBuilder().set_spatial_dimension(spatial_dimension).create_meta_data();
      meshB = stk::mesh::MeshBuilder(pmSub).create(metaB);
      build_mesh(*metaB, *meshB, info.num_elements, info.num_nodes, info.element_ids, element_ownerB, &info.elem_node_ids[0], info.node_sharingB, info.coordinates, create_faces);
    }
  }

  template <class KeyToTargetCheck, class RemoteKeyCheck>
  void run_test(KeyToTargetCheck keyToTargetCheck, RemoteKeyCheck keyCheck, stk::mesh::EntityRank field_rank = stk::topology::NODE_RANK)
  {

    std::unique_ptr<stk::transfer::TransferCopyByIdMeshAdapter> transferSourcePtr;
    std::unique_ptr<stk::transfer::TransferCopyByIdMeshAdapter> transferTargetPtr;
    auto prefix = field_prefix(field_rank);
    if(commOwnsMesh[0])
    {
      fill_mesh_values(*meshA);

      ScalarIntField    & scalarIntSourceField    = static_cast<ScalarIntField&>(*metaA->get_field(field_rank, prefix + " Scalar Int Field"));
      ScalarDoubleField & scalarDoubleSourceField = static_cast<ScalarDoubleField&>(*metaA->get_field(field_rank, prefix + " Scalar Double Field"));
      VectorDoubleField & vectorDoubleSourceField = static_cast<VectorDoubleField&>(*metaA->get_field(field_rank, prefix + " Vector Double Field"));
      std::vector<stk::mesh::Entity> sourceNodes;
      const bool sortById = true;
      stk::mesh::get_entities(*meshA, field_rank, metaA->locally_owned_part(), sourceNodes, sortById);
      std::vector<stk::mesh::FieldBase*> sourceFields;
      sourceFields.push_back(&scalarIntSourceField);
      sourceFields.push_back(&scalarDoubleSourceField);
      sourceFields.push_back(&vectorDoubleSourceField);
      transferSourcePtr.reset(new stk::transfer::TransferCopyByIdStkMeshAdapter(*meshA, sourceNodes, sourceFields, pm));
    }
    else
    {
      transferSourcePtr.reset(new stk::transfer::TransferCopyByIdMpmdMeshAdapter(pm, 3));
    }

    if(commOwnsMesh[1])
    {
      scalarIntTargetField    = static_cast<ScalarIntField*>(metaB->get_field(field_rank, prefix + " Scalar Int Field"));
      scalarDoubleTargetField = static_cast<ScalarDoubleField*>(metaB->get_field(field_rank, prefix + " Scalar Double Field"));
      vectorDoubleTargetField = static_cast<VectorDoubleField*>(metaB->get_field(field_rank, prefix + " Vector Double Field"));

      std::vector<stk::mesh::Entity> targetNodes;
      stk::mesh::Selector receiverSelector = metaB->locally_owned_part();
      if (receiverIncludesSharedNodes) {
        receiverSelector |= metaB->globally_shared_part();
      }
      stk::mesh::get_entities(*meshB, field_rank, receiverSelector, targetNodes);

      std::vector<stk::mesh::FieldBase*> targetFields;
      targetFields.push_back(scalarIntTargetField);
      targetFields.push_back(scalarDoubleTargetField);
      targetFields.push_back(vectorDoubleTargetField);
      transferTargetPtr.reset(new stk::transfer::TransferCopyByIdStkMeshAdapter(*meshB, targetNodes, targetFields, pm));
    }

    else
    {
      transferTargetPtr.reset(new stk::transfer::TransferCopyByIdMpmdMeshAdapter(pm, 3));
    }


    {
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(*transferSourcePtr, *transferTargetPtr,key_to_target_processor);

      stk::util::sort_and_unique(key_to_target_processor);
      keyToTargetCheck(key_to_target_processor);

      typedef stk::transfer::SearchById::MeshIDSet MeshIDSet;
      const MeshIDSet & remote_keys = copySearch.get_remote_keys();
      keyCheck(remote_keys);
    }

    stk::transfer::TransferCopyById transfer(copySearch, *transferSourcePtr, *transferTargetPtr);

    // Do the transfer
    //
    transfer.initialize();
    transfer.apply();
  }

  void check_target_fields(stk::mesh::EntityRank field_rank = stk::topology::NODE_RANK)
  {
    // Check "target" fields to make sure they hold the expected values
    //
    const double tolerance = 1.e-8;
    if (meshB) { //mesh B might not exist in MPMD case
      const stk::mesh::BucketVector & entityBuckets = meshB->get_buckets(field_rank, metaB->locally_owned_part() );
      for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
        stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
        for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
          stk::mesh::Entity entity = entityBucket[entityIndex];
          int * scalarIntTarget       = stk::mesh::field_data(*scalarIntTargetField, entity);
          double * scalarDoubleTarget = stk::mesh::field_data(*scalarDoubleTargetField, entity);
          double * vectorDoubleTarget = stk::mesh::field_data(*vectorDoubleTargetField, entity);

          EXPECT_EQ(static_cast<int>(meshB->identifier(entity)), *scalarIntTarget);
          EXPECT_NEAR(static_cast<double>(meshB->identifier(entity)), *scalarDoubleTarget, tolerance);
          for (size_t i = 0; i < spatial_dimension; ++i) {
            EXPECT_NEAR(static_cast<double>((meshB->identifier(entity)-1)*spatial_dimension+i), vectorDoubleTarget[i], tolerance);
          }
        }
      }
    }
    else
    {
      EXPECT_TRUE (meshA); //if mehsB didn't exist on processor, mehsA should exist
    }
  }

  void add_shared_nodes_to_receiver() { receiverIncludesSharedNodes = true; }

  virtual ~CopyTransferFixture()
  {
    MPI_Comm_free(&pmSub);
  }

protected:
  stk::ParallelMachine pm;
  stk::ParallelMachine pmSub;
//  Does this communicator own the mesh
  std::vector<bool> commOwnsMesh;
  SearchById copySearch;
  bool receiverIncludesSharedNodes = false;
  const size_t spatial_dimension = 3;
  std::shared_ptr<stk::mesh::MetaData> metaA;
  std::unique_ptr<stk::mesh::BulkData> meshA;
  std::shared_ptr<stk::mesh::MetaData> metaB;
  std::unique_ptr<stk::mesh::BulkData> meshB;
  ScalarIntField * scalarIntTargetField = nullptr;
  ScalarDoubleField * scalarDoubleTargetField = nullptr;
  VectorDoubleField * vectorDoubleTargetField = nullptr;
};

TYPED_TEST(CopyTransferFixture, copy0T0)
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

  const int color = 0; //all ranks same communicator
  this->init(MPI_COMM_WORLD, color);

  const int p_size = stk::parallel_machine_size( this->pm );
  if (p_size != 1) {
    return;
  }

  std::vector<int> element_ownerA = {0};
  std::vector<int> element_ownerB = {0};
  const int expected_target_processor = 0;
  this->build_fixture(element_ownerA, element_ownerB, TwoElemMeshInfo());
  this->run_test(TwoElemKeyToTargetTest(expected_target_processor),
      [=](const stk::transfer::SearchById::MeshIDSet & remote_keys){
      EXPECT_TRUE( remote_keys.empty() );
  });

  this->check_target_fields();
}

TYPED_TEST(CopyTransferFixture, copy0T1)
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

  const int color = 0; //all ranks same communicator
  this->init(MPI_COMM_WORLD, color);

  const int p_size = stk::parallel_machine_size( this->pm );
  if (p_size != 2) {
    return;
  }

  std::vector<int> element_ownerA = {0};
  std::vector<int> element_ownerB = {1};
  const int expected_target_processor = 1;
  this->build_fixture(element_ownerA, element_ownerB, TwoElemMeshInfo());
  const int p_rank = stk::parallel_machine_rank( this->pm );
  this->run_test(TwoElemKeyToTargetTest(expected_target_processor),
      [=](const stk::transfer::SearchById::MeshIDSet & remote_keys){
      if (p_rank == 0) {
        EXPECT_TRUE( remote_keys.empty() );
      } else {
        EXPECT_EQ( 8u, remote_keys.size() );
      }
  });

  this->check_target_fields();
}

TYPED_TEST(CopyTransferFixture, copy0T1_MPMD)
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
// color            0                       1
// split comm rank  0                       0

  const int color = stk::parallel_machine_rank(MPI_COMM_WORLD); //subcommunicator based on global rank
  this->init(MPI_COMM_WORLD, color, {0,1});
  const int p_size = stk::parallel_machine_size( this->pm );

  if (p_size != 2) {
    return;
  }

  std::vector<int> element_ownerA_subcomm_rank = {0};
  std::vector<int> element_ownerB_subcomm_rank = {0};
  const int expected_target_processor = 1;
  this->build_fixture(element_ownerA_subcomm_rank, element_ownerB_subcomm_rank, TwoElemMeshInfo());
  const int p_rank = stk::parallel_machine_rank( this->pm );
  this->run_test(TwoElemKeyToTargetTest(expected_target_processor),
      [=](const stk::transfer::SearchById::MeshIDSet & remote_keys){
      if (p_rank == 0) {
        EXPECT_TRUE( remote_keys.empty() );
      } else {
        EXPECT_EQ( 8u, remote_keys.size() );
      }
  });

  this->check_target_fields();
}

TYPED_TEST(CopyTransferFixture, copy1T0)
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

  const int color = 0; //all ranks same communicator
  this->init(MPI_COMM_WORLD, color);
  const int p_size = stk::parallel_machine_size( this->pm );

  if (p_size != 2) {
    return;
  }

  std::vector<int> element_ownerA = {1};
  std::vector<int> element_ownerB = {0};
  const int expected_target_processor = 0;
  this->build_fixture(element_ownerA, element_ownerB, TwoElemMeshInfo());
  stk::ParallelMachine comm = this->pm;
  this->run_test(TwoElemKeyToTargetTest(expected_target_processor), [=](const stk::transfer::SearchById::MeshIDSet & remote_keys){
      const int p_rank = stk::parallel_machine_rank( comm );
      if (p_rank == 0) {
        EXPECT_EQ( 8u, remote_keys.size() );
      } else {
        EXPECT_TRUE( remote_keys.empty() );
      }
  });

  this->check_target_fields();
}

TYPED_TEST(CopyTransferFixture, copy1T0_MPMD)
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
// color            1                       0
// split comm rank  0                       0

  int color = stk::parallel_machine_rank(MPI_COMM_WORLD); //subcommunicator based on global rank
  this->init(MPI_COMM_WORLD, color, {1,0});
  const int p_size = stk::parallel_machine_size( this->pm );

  if (p_size != 2) {
    return;
  }

  std::vector<int> element_ownerA_subcomm_rank = {0};
  std::vector<int> element_ownerB_subcomm_rank = {0};
  const int expected_target_processor = 0;
  this->build_fixture(element_ownerA_subcomm_rank, element_ownerB_subcomm_rank, TwoElemMeshInfo());
  stk::ParallelMachine comm = this->pm;
  this->run_test(TwoElemKeyToTargetTest(expected_target_processor),
      [=](const stk::transfer::SearchById::MeshIDSet & remote_keys){
      const int p_rank = stk::parallel_machine_rank(comm);
      if (p_rank == 0) {
        EXPECT_EQ( 8u, remote_keys.size() );
      } else {
        EXPECT_TRUE( remote_keys.empty() );
      }
  });

  this->check_target_fields();
}

namespace {
stk::transfer::SearchById::KeyToTargetProcessor get_01T10_key_to_target_processor_gold(stk::ParallelMachine pm)
{
  const int p_rank = stk::parallel_machine_rank( pm );
  stk::transfer::SearchById::KeyToTargetProcessor gold;
  if (0 == p_rank) {
    gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,1), 1);
    gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,4), 1);
    gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,7), 1);
    gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,10), 1);
    gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,2), 0);
    gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,5), 0);
    gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,8), 0);
    gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,11), 0);
  } else {
    gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,3), 0);
    gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,6), 0);
    gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,9), 0);
    gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,12), 0);
  }
  stk::util::sort_and_unique(gold);
  return gold;
}


stk::transfer::SearchById::MeshIDSet get_01T10_remote_key_gold(stk::ParallelMachine pm)
{
  const int p_rank = stk::parallel_machine_rank( pm );
  stk::transfer::SearchById::MeshIDSet gold_remote_keys;
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
  return gold_remote_keys;
}

}

TYPED_TEST(CopyTransferFixture, copy01T10)
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

  const int color = 0; //all ranks same communicator
  this->init(MPI_COMM_WORLD, color);

  const int p_size = stk::parallel_machine_size( this->pm );
  if (p_size != 2) {
    return;
  }

  std::vector<int> element_ownerA = {0, 1};
  std::vector<int> element_ownerB = {1, 0};
  this->build_fixture(element_ownerA, element_ownerB, FourElemMeshInfo());
  stk::ParallelMachine comm = this->pm;
  this->run_test([=](const stk::transfer::SearchById::KeyToTargetProcessor & key_to_target_processor)
      {
    auto gold = get_01T10_key_to_target_processor_gold(comm);
    EXPECT_EQ(gold, key_to_target_processor);
      },
      [=](const stk::transfer::SearchById::MeshIDSet & remote_keys){
        auto gold_remote_keys = get_01T10_remote_key_gold(comm);
        EXPECT_EQ(remote_keys, gold_remote_keys);
      });

  this->check_target_fields();
}

TYPED_TEST(CopyTransferFixture, copy01T32_MPMD)
{
//    ^ Y       ID.owning_proc
//    |
//    |    X              meshA                             meshB
//    .---->    4.0--------5.0--------6.1         4.3--------5.2--------6.2
//   /          /|         /|         /|          /|         /|         /|
//  /          / |        / |        / |         / |        / |        / |
// v Z      10.0-------11.0-------12.1 |      10.3-------11.2-------12.2 |
//            |  |  1.0  |  |  2.1  |  |  -->   |  |  1.3  |  |  2.2  |  |
//            | 1.0------|-2.0------|-3.1       | 1.3------|-2.2------|-3.2
//            | /        | /        | /         | /        | /        | /
//            |/         |/         |/          |/         |/         |/
//           7.0--------8.0--------9.1         7.3--------8.2--------9.2
//
  //   color       0         0                        1         1

  //global ranks 0, 1 are color 0, global ranks 2, 3 are color 1
  //subcommunicator ranks are the same as in the non mpmd 01T10 case above

  int color = stk::parallel_machine_rank(MPI_COMM_WORLD) / 2;
  this->init(MPI_COMM_WORLD, color, {0,1});
  const int p_size = stk::parallel_machine_size( this->pm );

  if (p_size != 4) {
    return;
  }

  std::vector<int> element_ownerA_subcomm_rank = {0, 1};
  std::vector<int> element_ownerB_subcomm_rank = {1, 0};

  this->build_fixture(element_ownerA_subcomm_rank , element_ownerB_subcomm_rank, FourElemMeshInfo());
  stk::ParallelMachine comm = this->pm;
  stk::ParallelMachine commSub = this->pmSub;
  this->run_test([=](const stk::transfer::SearchById::KeyToTargetProcessor & key_to_target_processor)
      {
        if (color == 1) return;

        //Map wrt subcommunicators is same as non mpmd 01T10 case above
        auto gold = get_01T10_key_to_target_processor_gold(commSub);

        //adjust from subcommunicator to global
        for (auto && elem : gold)
        {
          elem.second+=2;
        }
        EXPECT_EQ(gold, key_to_target_processor);
      },
      [=](const stk::transfer::SearchById::MeshIDSet & remote_keys){
        if (color == 0) return;
        auto gold_remote_keys = get_01T10_remote_key_gold(commSub);
        //Add nodes that were on the same proc in the non mpmd case
        if (stk::parallel_machine_rank(comm) == 2)
        {
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,2).m_value);
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,5).m_value);
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,8).m_value);
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,11).m_value);
        }

        EXPECT_EQ(remote_keys, gold_remote_keys);
      });

  this->check_target_fields();
}

TYPED_TEST(CopyTransferFixture, copy001T011)
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

  const int color = 0; //all ranks same communicator
  this->init(MPI_COMM_WORLD, color);

  const int p_size = stk::parallel_machine_size( this->pm );
  if (p_size != 2) {
    return;
  }

  std::vector<int> element_ownerA = {0, 0, 1};
  std::vector<int> element_ownerB = {0, 1, 1};
  this->build_fixture(element_ownerA, element_ownerB, SixElemMeshInfo());
  const int p_rank = stk::parallel_machine_rank( this->pm );
  this->run_test([=](const stk::transfer::SearchById::KeyToTargetProcessor & key_to_target_processor)
      {
      stk::transfer::SearchById::KeyToTargetProcessor gold;
      if (0 == p_rank) {
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,1), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,5), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,9), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,13), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,2), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,6), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,10), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,14), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,3), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,7), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,11), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,15), 1);
      } else {
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,4), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,8), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,12), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,16), 1);
      }
      stk::util::sort_and_unique(gold);
      EXPECT_EQ(gold, key_to_target_processor);
      },

      [=](const stk::transfer::SearchById::MeshIDSet & remote_keys){
        typedef stk::transfer::SearchById::MeshIDSet MeshIDSet;
        MeshIDSet gold_remote_keys;
        if (0 == p_rank) {
        } else {
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,3).m_value);
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,7).m_value);
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,11).m_value);
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,15).m_value);
        }
        EXPECT_EQ(remote_keys, gold_remote_keys);
      });

  this->check_target_fields();
}

TYPED_TEST(CopyTransferFixture, copy001T011Element)
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
  const int color = 0; //all ranks same communicator
  this->init(MPI_COMM_WORLD, color);

  const int p_size = stk::parallel_machine_size( this->pm );
  if (p_size != 2) {
    return;
  }

  std::vector<int> element_ownerA = {0, 0, 1};
  std::vector<int> element_ownerB = {0, 1, 1};
  this->build_fixture(element_ownerA, element_ownerB, SixElemMeshInfo());
  const int p_rank = stk::parallel_machine_rank( this->pm );
  this->run_test([=](const stk::transfer::SearchById::KeyToTargetProcessor & key_to_target_processor)
      {
        stk::transfer::SearchById::KeyToTargetProcessor gold;
        if (0 == p_rank) {
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,1), 0);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,2), 1);
        } else {
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,3), 1);
        }
        stk::util::sort_and_unique(gold);
        EXPECT_EQ(gold, key_to_target_processor);
      },

      [=](const stk::transfer::SearchById::MeshIDSet & remote_keys){
        typedef stk::transfer::SearchById::MeshIDSet MeshIDSet;
        MeshIDSet gold_remote_keys;
        if (1 == p_rank) {
          gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::ELEM_RANK,2).m_value);
        }
        EXPECT_EQ(remote_keys, gold_remote_keys);
      }, stk::topology::ELEM_RANK);

  this->check_target_fields(stk::topology::ELEMENT_RANK);
}

TYPED_TEST(CopyTransferFixture, copy001T011Face)
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

  const int color = 0; //all ranks same communicator
  this->init(MPI_COMM_WORLD, color);

  const int p_size = stk::parallel_machine_size( this->pm );
  if (p_size != 2) {
    return;
  }

  std::vector<int> element_ownerA = {0, 0, 1};
  std::vector<int> element_ownerB = {0, 1, 1};

  auto meshInfo = SixElemMeshInfo();
  meshInfo.elem_node_ids = {
        { 9,10,2,1,13,14,6,5},
        {10,11,3,2,14,15,7,6},
        {11,12,4,3,15,16,8,7} };
  meshInfo.node_sharingA = { -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1,
                            -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1 };
  meshInfo.node_sharingB = { -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1,
                            -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1 };
  meshInfo.coordinates = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
                                {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0},
                                {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0},
                                {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0} };

  const bool create_faces = true;
  this->build_fixture(element_ownerA, element_ownerB, meshInfo, create_faces);
  const int p_rank = stk::parallel_machine_rank( this->pm );
  auto & mesh_a = *this->meshA;
  auto & mesh_b = *this->meshB;
  this->run_test([&, p_rank](const stk::transfer::SearchById::KeyToTargetProcessor & key_to_target_processor)
      {
      typedef stk::transfer::SearchById::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor gold;
      if (0 == p_rank) {
        stk::mesh::Entity elem1 = mesh_a.get_entity(stk::topology::ELEM_RANK, 1);
        stk::mesh::Entity elem2 = mesh_a.get_entity(stk::topology::ELEM_RANK, 2);
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem1, 0)), 0);
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem1, 1)), 0);
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem1, 2)), 0);
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem1, 3)), 0);
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem1, 4)), 0);
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem1, 5)), 0);

        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem2, 0)), 1);
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem2, 1)), 1);
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem2, 2)), 1);
        //gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(meshA, elem2, 3)), 0);  // Already in map from elem1
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem2, 4)), 1);
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem2, 5)), 1);
      } else {
        stk::mesh::Entity elem3 = mesh_a.get_entity(stk::topology::ELEM_RANK, 3);
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem3, 0)), 1);
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem3, 1)), 1);
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem3, 2)), 1);
        //gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(meshA, elem3, 3)), 0);  // Not owned by this proc
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem3, 4)), 1);
        gold.emplace_back(mesh_a.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_a, elem3, 5)), 1);
      }
      stk::util::sort_and_unique(gold);
      EXPECT_EQ(gold, key_to_target_processor);
      },

      [&, p_rank](const stk::transfer::SearchById::MeshIDSet & remote_keys){
      typedef stk::transfer::SearchById::MeshIDSet MeshIDSet;
      MeshIDSet gold_remote_keys;
      if (1 == p_rank) {
        stk::mesh::Entity elem2 = mesh_b.get_entity(stk::topology::ELEM_RANK, 2);
        gold_remote_keys.insert(mesh_b.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_b, elem2, 0)).m_value);
        gold_remote_keys.insert(mesh_b.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_b, elem2, 1)).m_value);
        gold_remote_keys.insert(mesh_b.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_b, elem2, 2)).m_value);
        //gold_remote_keys.insert(meshB.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(meshB, elem2, 3)).m_value);  // Not received because not owned
        gold_remote_keys.insert(mesh_b.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_b, elem2, 4)).m_value);
        gold_remote_keys.insert(mesh_b.entity_key(stk::mesh::get_side_entity_for_elem_side_pair(mesh_b, elem2, 5)).m_value);
      }
      EXPECT_EQ(remote_keys, gold_remote_keys);
      }, stk::topology::FACE_RANK);

  this->check_target_fields(stk::topology::FACE_RANK);
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

  stk::transfer::SearchByIdGeometric geometricSearch;
  stk::transfer::SearchByIdCommAll commAllSearch;
  stk::transfer::SearchById * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_EQ(copySearchPtr, &geometricSearch);
    }
    stk::transfer::SearchById & copySearch = *copySearchPtr;

    const size_t spatial_dimension = 3;
    const size_t num_elements = 3;
    const size_t num_nodes = 16;
    stk::mesh::EntityIdVector element_ids {1, 2, 3};
    std::vector<int> element_ownerA = {0, 0, 1};
    std::vector<int> element_ownerB = {0, 1, 1};
    std::vector<stk::mesh::EntityIdVector> elem_node_ids {
        { 9,10,2,1,13,14,6,5},
        {10,11,3,2,14,15,7,6},
        {11,12,4,3,15,16,8,7} };
    std::vector<int> node_sharingA = { -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1,
                                       -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1 };
    std::vector<int> node_sharingB = { -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1, -1, -1,
                                       -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1 };
    std::vector<std::vector<double>> coordinates = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
                                                     {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0},
                                                     {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0},
                                                     {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0} };


    // Set up the "source" mesh for the transfer
    //
    std::shared_ptr<stk::mesh::BulkData> meshA = build_mesh(spatial_dimension, pm);
    stk::mesh::MetaData& metaA = meshA->mesh_meta_data();
    build_mesh(metaA,
               *meshA,
               num_elements,
               num_nodes,
               element_ids,
               element_ownerA,
               &elem_node_ids[0],
               node_sharingA,
               coordinates);

    // Set up the "target" mesh for the transfer
    //
    std::shared_ptr<stk::mesh::BulkData> meshBPtr = build_mesh(spatial_dimension, pm);
    stk::mesh::MetaData& metaB = meshBPtr->mesh_meta_data();
    stk::mesh::BulkData& meshB = *meshBPtr;
    build_mesh(metaB,
               meshB,
               num_elements,
               num_nodes,
               element_ids,
               element_ownerB,
               &elem_node_ids[0],
               node_sharingB,
               coordinates);


    stk::mesh::EntityIdVector shell_node_ids[] {
        { 9,10,14,13}, {10,2,6,14}, {2,1,5,6}, { 9,13,5, 1}, { 9,1,2,10}, {13,14,6,5},
        {10,11,15,14}, {11,3,7,15}, {3,2,6,7}, {10, 2,6,14}, {10,2,3,11}, {14,15,7,6},
        {11,12,16,15}, {12,4,8,16}, {4,3,7,8}, {11, 3,7,15}, {11,3,4,12}, {15,16,8,7} };
    stk::mesh::EntityIdVector elem_shell_ids[] {
        {10, 11, 12, 13, 14, 15},
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
                       *meshA,
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
    fill_mesh_values(*meshA);

    ScalarIntField    & scalarSourceIntField    = static_cast<ScalarIntField&>(*metaA.get_field(stk::topology::ELEM_RANK, "Shell Scalar Int Field"));
    ScalarDoubleField & scalarSourceDoubleField = static_cast<ScalarDoubleField&>(*metaA.get_field(stk::topology::ELEM_RANK, "Shell Scalar Double Field"));
    VectorDoubleField & vectorSourceDoubleField = static_cast<VectorDoubleField&>(*metaA.get_field(stk::topology::ELEM_RANK, "Shell Vector Double Field"));
    ScalarIntField    & scalarTargetIntField    = static_cast<ScalarIntField&>(*metaB.get_field(stk::topology::ELEM_RANK, "Shell Scalar Int Field"));
    ScalarDoubleField & scalarTargetDoubleField = static_cast<ScalarDoubleField&>(*metaB.get_field(stk::topology::ELEM_RANK, "Shell Scalar Double Field"));
    VectorDoubleField & vectorTargetDoubleField = static_cast<VectorDoubleField&>(*metaB.get_field(stk::topology::ELEM_RANK, "Shell Vector Double Field"));

    // Set up the transfer
    //
    stk::mesh::Part & shell_part = *metaA.get_part("shell_part");

    std::vector<stk::mesh::Entity> sourceShells;
    const bool sortById = true;
    stk::mesh::get_entities(*meshA, stk::topology::ELEM_RANK,
                            metaA.locally_owned_part() & shell_part,
                            sourceShells, sortById);

    std::vector<stk::mesh::Entity> targetShells;
    stk::mesh::get_entities(meshB, stk::topology::ELEM_RANK,
                            metaB.locally_owned_part() & shell_part,
                            targetShells, sortById);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceIntField);
    sourceFields.push_back(&scalarSourceDoubleField);
    sourceFields.push_back(&vectorSourceDoubleField);
    stk::transfer::TransferCopyByIdStkMeshAdapter transferSource(*meshA, sourceShells, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetIntField);
    targetFields.push_back(&scalarTargetDoubleField);
    targetFields.push_back(&vectorTargetDoubleField);
    stk::transfer::TransferCopyByIdStkMeshAdapter transferTarget(meshB, targetShells, targetFields);

    {
      const int p_rank = stk::parallel_machine_rank( pm );
      typedef stk::transfer::SearchById::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);
      stk::util::sort_and_unique(key_to_target_processor);

      KeyToTargetProcessor gold;
      for (int shell_index=0; shell_index<num_shells; ++shell_index) {
        if (0 == p_rank) {
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,10+100*shell_index), 0);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,11+100*shell_index), 0);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,12+100*shell_index), 0);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,13+100*shell_index), 0);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,14+100*shell_index), 0);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,15+100*shell_index), 0);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,16+100*shell_index), 1);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,17+100*shell_index), 1);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,18+100*shell_index), 1);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,19+100*shell_index), 1);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,20+100*shell_index), 1);
        } else {
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,21+100*shell_index), 1);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,22+100*shell_index), 1);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,23+100*shell_index), 1);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,24+100*shell_index), 1);
          gold.emplace_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK,25+100*shell_index), 1);
        }
      }
      stk::util::sort_and_unique(gold);
      EXPECT_EQ(gold, key_to_target_processor);

      typedef stk::transfer::SearchById::MeshIDSet MeshIDSet;
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
      EXPECT_EQ(copySearch.get_remote_keys(), gold_remote_keys);
    }
    stk::transfer::TransferCopyById transfer(copySearch, transferSource, transferTarget);

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
        int    * scalarIntTarget    = stk::mesh::field_data(scalarTargetIntField, entity);
        double * scalarDoubleTarget = stk::mesh::field_data(scalarTargetDoubleField, entity);
        double * vectorDoubleTarget = stk::mesh::field_data(vectorTargetDoubleField, entity);

        EXPECT_EQ(static_cast<int>(meshB.identifier(entity)), *scalarIntTarget);
        EXPECT_NEAR(static_cast<double>(meshB.identifier(entity)), *scalarDoubleTarget, tolerance);
        for (size_t i = 0; i < spatial_dimension; ++i) {
          EXPECT_NEAR(static_cast<double>((meshB.identifier(entity)-1)*spatial_dimension+i), vectorDoubleTarget[i], tolerance);
        }
      }
    }
  }
}

namespace {

std::vector<int> get_node_sharing_000()
{
    return { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
}

std::vector<int> get_node_sharing_012()
{
  return { -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1,
      -1, 0, 2, -1, -1, 0, 2, -1, -1, 0, 2, -1, -1, 0, 2, -1,
      -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1 };
}

SixElemMeshInfo create_SixElemMeshInfo_012T000()
{
  auto meshInfo = SixElemMeshInfo();
  meshInfo.node_sharingA = get_node_sharing_012();
  meshInfo.node_sharingB = get_node_sharing_000();
  return meshInfo;
}

SixElemMeshInfo create_SixElemMeshInfo_000T012()
{
  auto meshInfo = SixElemMeshInfo();
  meshInfo.node_sharingA = get_node_sharing_000();
  meshInfo.node_sharingB = get_node_sharing_012();
  return meshInfo;
}
}

TYPED_TEST(CopyTransferFixture, copy012T000)
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

  const int color = 0; //all ranks same communicator
  this->init(MPI_COMM_WORLD, color);

  const int p_size = stk::parallel_machine_size( this->pm );
  if (p_size != 3) {
    return;
  }

  std::vector<int> element_ownerA = {0, 1, 2};
  std::vector<int> element_ownerB = {0, 0, 0};

  auto meshInfo = create_SixElemMeshInfo_012T000();

  this->build_fixture(element_ownerA, element_ownerB, meshInfo);
  const int p_rank = stk::parallel_machine_rank( this->pm );
  this->run_test([=](const stk::transfer::SearchById::KeyToTargetProcessor & key_to_target_processor)
    {
      stk::transfer::SearchById::KeyToTargetProcessor gold;
      if (0 == p_rank) {
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,1), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,5), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,9), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,13), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,2), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,6), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,10), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,14), 0);
      } else if (1 == p_rank) {
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,3), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,7), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,11), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,15), 0);
      } else {
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,4), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,8), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,12), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,16), 0);
      }
      stk::util::sort_and_unique(gold);
      EXPECT_EQ(gold, key_to_target_processor);
      },

      [=](const stk::transfer::SearchById::MeshIDSet & remote_keys){
      typedef stk::transfer::SearchById::MeshIDSet MeshIDSet;
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
      EXPECT_EQ(remote_keys, gold_remote_keys);
    });

  this->check_target_fields();

}

TYPED_TEST(CopyTransferFixture, copy000T012)
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

  const int color = 0; //all ranks same communicator
  this->init(MPI_COMM_WORLD, color);

  const int p_size = stk::parallel_machine_size( this->pm );
  if (p_size != 3) {
    return;
  }


  auto meshInfo = create_SixElemMeshInfo_000T012();
  std::vector<int> element_ownerA = {0, 0, 0};
  std::vector<int> element_ownerB = {0, 1, 2};
  this->build_fixture(element_ownerA, element_ownerB, meshInfo);
  const int p_rank = stk::parallel_machine_rank( this->pm );
  this->run_test([=](const stk::transfer::SearchById::KeyToTargetProcessor & key_to_target_processor)
    {
      typedef stk::transfer::SearchById::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor gold;
      if (0 == p_rank) {
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,1), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,5), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,9), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,13), 0);

        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,2), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,6), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,10), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,14), 0);

        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,3), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,7), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,11), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,15), 1);

        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,4), 2);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,8), 2);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,12), 2);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,16), 2);
      }
      stk::util::sort_and_unique(gold);
      EXPECT_EQ(gold, key_to_target_processor);
      },

      [=](const stk::transfer::SearchById::MeshIDSet & remote_keys){
      typedef stk::transfer::SearchById::MeshIDSet MeshIDSet;
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
      EXPECT_EQ(remote_keys, gold_remote_keys);
    });

  this->check_target_fields();
}

TYPED_TEST(CopyTransferFixture, copy000T012_multipleTargets)
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
//  This is the same as the previous test, but with the destination receiving values
//  on both locally-owned and shared instead of just locally-owned.  So, every node
//  that isn't ghosted will receive some data.

  const int color = 0; //all ranks same communicator
  this->init(MPI_COMM_WORLD, color);

  const int p_size = stk::parallel_machine_size( this->pm );
  if (p_size != 3) {
    return;
  }


  auto meshInfo = create_SixElemMeshInfo_000T012();
  std::vector<int> element_ownerA = {0, 0, 0};
  std::vector<int> element_ownerB = {0, 1, 2};
  this->build_fixture(element_ownerA, element_ownerB, meshInfo);
  const int p_rank = stk::parallel_machine_rank( this->pm );
  this->add_shared_nodes_to_receiver();

  this->run_test([=](const stk::transfer::SearchById::KeyToTargetProcessor & key_to_target_processor)
    {
      typedef stk::transfer::SearchById::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor gold;
      if (0 == p_rank) {
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,1), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,5), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,9), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,13), 0);

        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,2), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,2), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,6), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,6), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,10), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,10), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,14), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,14), 1);

        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,3), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,3), 2);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,7), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,7), 2);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,11), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,11), 2);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,15), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,15), 2);

        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,4), 2);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,8), 2);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,12), 2);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,16), 2);
      }
      stk::util::sort_and_unique(gold);
      EXPECT_EQ(gold, key_to_target_processor);
      },

      [=](const stk::transfer::SearchById::MeshIDSet & remote_keys){
      typedef stk::transfer::SearchById::MeshIDSet MeshIDSet;
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
      } else if (2 == p_rank) {
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,3).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,7).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,11).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,15).m_value);

        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,4).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,8).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,12).m_value);
        gold_remote_keys.insert(stk::mesh::EntityKey(stk::topology::NODE_RANK,16).m_value);
      }
      EXPECT_EQ(remote_keys, gold_remote_keys);
    });

  this->check_target_fields();
}

TYPED_TEST(CopyTransferFixture, copy0011T1010)
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
  const int color = 0; //all ranks same communicator
  this->init(MPI_COMM_WORLD, color);

  const int p_size = stk::parallel_machine_size( this->pm );
  if (p_size != 2) {
    return;
  }

  struct EightElemMeshInfo
  {
    const size_t spatial_dimension = 3;
    const size_t num_elements = 4;
    const size_t num_nodes = 20;
    stk::mesh::EntityIdVector element_ids = {1, 2, 3, 4};
    std::vector<stk::mesh::EntityIdVector> elem_node_ids = {
        {1, 2, 7, 6, 11, 12, 17, 16},
        {2, 3, 8, 7, 12, 13, 18, 17},
        {3, 4, 9, 8, 13, 14, 19, 18},
        {4, 5, 10, 9, 14, 15, 20, 19} };
    std::vector<int> node_sharingA = { -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1,
      -1, -1, 0, -1, -1, -1, -1, 0, -1, -1, -1, -1, 0, -1, -1, -1, -1, 0, -1, -1 };
    std::vector<int> node_sharingB = { -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1,
      -1, 0, 0, 0, -1, -1, 0, 0, 0, -1, -1, 0, 0, 0, -1, -1, 0, 0, 0, -1 };
    std::vector<std::vector<double>> coordinates = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {4.0, 0.0, 0.0},
      {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0}, {4.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0}, {4.0, 0.0, 1.0},
      {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0}, {4.0, 1.0, 1.0} };
  };

  std::vector<int> element_ownerA = {0, 0, 1, 1};
  std::vector<int> element_ownerB = {1, 0, 1, 0};
  this->build_fixture(element_ownerA, element_ownerB, EightElemMeshInfo());
  const int p_rank = stk::parallel_machine_rank( this->pm );
  this->run_test([=](const stk::transfer::SearchById::KeyToTargetProcessor & key_to_target_processor)
    {
      typedef stk::transfer::SearchById::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor gold;
      if (0 == p_rank) {
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,1), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,6), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,11), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,16), 1);

        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,2), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,7), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,12), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,17), 0);

        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,3), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,8), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,13), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,18), 0);
      } else if (1 == p_rank) {
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,4), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,9), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,14), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,19), 0);

        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,5), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,10), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,15), 0);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,20), 0);
      }
      stk::util::sort_and_unique(gold);
      EXPECT_EQ(gold, key_to_target_processor);
      },

      [=](const stk::transfer::SearchById::MeshIDSet & remote_keys){
      typedef stk::transfer::SearchById::MeshIDSet MeshIDSet;
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
      EXPECT_EQ(remote_keys, gold_remote_keys);
    });

  this->check_target_fields();
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

  stk::transfer::SearchByIdGeometric geometricSearch;
  stk::transfer::SearchByIdCommAll commAllSearch;
  stk::transfer::SearchById * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_EQ(copySearchPtr, &geometricSearch);
    }
    stk::transfer::SearchById & copySearch = *copySearchPtr;

    const size_t spatial_dimension = 3;
    const size_t num_elements = 1;
    const size_t num_nodes = 8;
    stk::mesh::EntityIdVector element_ids {1};
    std::vector<int> element_owner = {0};
    stk::mesh::EntityIdVector elem_node_ids[] { {1, 2, 3, 4, 5, 6, 7, 8} };
    std::vector<int> node_sharing = { -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1 };
    std::vector<std::vector<double>> coordinates = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0} };

    // Set up the "source" mesh for the transfer
    //
    std::shared_ptr<stk::mesh::BulkData> meshAPtr = build_mesh(spatial_dimension, pm);
    stk::mesh::MetaData& metaA = meshAPtr->mesh_meta_data();
    stk::mesh::BulkData& meshA = *meshAPtr;
    build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_owner, elem_node_ids, node_sharing, coordinates);

    // Set up the "target" mesh for the transfer without creating any elements
    //
    std::shared_ptr<stk::mesh::BulkData> meshBPtr = build_mesh(spatial_dimension, pm);
    stk::mesh::MetaData& metaB = meshBPtr->mesh_meta_data();
    stk::mesh::BulkData& meshB = *meshBPtr;

    int int_init_vals = std::numeric_limits<int>::max();
    double double_init_vals[] = {std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max()};
    ScalarIntField    & scalarTargetIntField    = metaB.declare_field<int>(stk::topology::NODE_RANK, "Node Scalar Int Field");
    ScalarDoubleField & scalarTargetDoubleField = metaB.declare_field<double>(stk::topology::NODE_RANK, "Node Scalar Double Field");
    VectorDoubleField & vectorTargetDoubleField = metaB.declare_field<double>(stk::topology::NODE_RANK, "Node Vector Double Field");
    VectorDoubleField & coordsTargetField = metaB.declare_field<double>(stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field_on_mesh(scalarTargetIntField, metaB.universal_part(), &int_init_vals);
    stk::mesh::put_field_on_mesh(scalarTargetDoubleField, metaB.universal_part(), double_init_vals);
    stk::mesh::put_field_on_mesh(vectorTargetDoubleField, metaB.universal_part(), 3, double_init_vals);
    stk::mesh::put_field_on_mesh(coordsTargetField, metaB.universal_part(), 3, double_init_vals);
    metaB.commit();

    // Fill "source" fields with valid data
    //
    fill_mesh_values(meshA);

    ScalarIntField    & scalarSourceIntField    = static_cast<ScalarIntField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Int Field"));
    ScalarDoubleField & scalarSourceDoubleField = static_cast<ScalarDoubleField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Double Field"));
    VectorDoubleField & vectorSourceDoubleField = static_cast<VectorDoubleField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Double Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    const bool sortById = true;
    stk::mesh::get_entities(meshA, stk::topology::NODE_RANK, metaA.locally_owned_part(), sourceNodes, sortById);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_entities(meshB, stk::topology::NODE_RANK, metaB.locally_owned_part(), targetNodes, sortById);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceIntField);
    sourceFields.push_back(&scalarSourceDoubleField);
    sourceFields.push_back(&vectorSourceDoubleField);
    stk::transfer::TransferCopyByIdStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetIntField);
    targetFields.push_back(&scalarTargetDoubleField);
    targetFields.push_back(&vectorTargetDoubleField);
    stk::transfer::TransferCopyByIdStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    //  GeometricTransfer
    //  stk::transfer::GeometricTransfer<
    //    class stk::transfer::LinearInterpolate<
    //      class stk::transfer::STKNode,
    //      class stk::transfer::STKNode
    //    >
    //  > transfer(transferSource, transferTarget, "copy0T_ unit test");

    {
      typedef stk::transfer::SearchById::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold;
      EXPECT_EQ(gold, key_to_target_processor);

      typedef stk::transfer::SearchById::MeshIDSet MeshIDSet;
      MeshIDSet gold_remote_keys;
      EXPECT_EQ(copySearch.get_remote_keys(), gold_remote_keys);
    }
    stk::transfer::TransferCopyById transfer(copySearch, transferSource, transferTarget);

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

  stk::transfer::SearchByIdGeometric geometricSearch;
  stk::transfer::SearchByIdCommAll commAllSearch;
  stk::transfer::SearchById * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_EQ(copySearchPtr, &geometricSearch);
    }
    stk::transfer::SearchById & copySearch = *copySearchPtr;

    const size_t spatial_dimension = 3;
    const size_t num_elements = 1;
    const size_t num_nodes = 8;
    stk::mesh::EntityIdVector element_ids {1};
    std::vector<int> element_owner = {0};
    stk::mesh::EntityIdVector elem_node_ids[] { {1, 2, 3, 4, 5, 6, 7, 8} };
    std::vector<int> node_sharing = { -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1 };
    std::vector<std::vector<double>> coordinates = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0} };

    // Set up the "source" mesh for the transfer
    //
    std::shared_ptr<stk::mesh::BulkData> meshAPtr = build_mesh(spatial_dimension, pm);
    stk::mesh::MetaData& metaA = meshAPtr->mesh_meta_data();
    stk::mesh::BulkData& meshA = *meshAPtr;

    int int_init_vals = std::numeric_limits<int>::max();
    double double_init_vals[] = {std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max()};
    ScalarIntField    & scalarSourceIntField    = metaA.declare_field<int>(stk::topology::NODE_RANK, "Node Scalar Int Field");
    ScalarDoubleField & scalarSourceDoubleField = metaA.declare_field<double>(stk::topology::NODE_RANK, "Node Scalar Double Field");
    VectorDoubleField & vectorSourceDoubleField = metaA.declare_field<double>(stk::topology::NODE_RANK, "Node Vector Double Field");
    VectorDoubleField & coordsSourceField = metaA.declare_field<double>(stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field_on_mesh(scalarSourceIntField, metaA.universal_part(), &int_init_vals);
    stk::mesh::put_field_on_mesh(scalarSourceDoubleField, metaA.universal_part(), double_init_vals);
    stk::mesh::put_field_on_mesh(vectorSourceDoubleField, metaA.universal_part(), 3, double_init_vals);
    stk::mesh::put_field_on_mesh(coordsSourceField, metaA.universal_part(), 3, double_init_vals);
    metaA.commit();

    // Set up the "target" mesh for the transfer without creating any elements
    //
    std::shared_ptr<stk::mesh::BulkData> meshBPtr = build_mesh(spatial_dimension, pm);
    stk::mesh::MetaData& metaB = meshBPtr->mesh_meta_data();
    stk::mesh::BulkData& meshB = *meshBPtr;
    build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_owner, elem_node_ids, node_sharing, coordinates);

    ScalarIntField & scalarTargetIntField       = static_cast<ScalarIntField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Int Field"));
    ScalarDoubleField & scalarTargetDoubleField = static_cast<ScalarDoubleField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Double Field"));
    VectorDoubleField & vectorTargetDoubleField = static_cast<VectorDoubleField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Double Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    const bool sortById = true;
    stk::mesh::get_entities(meshA, stk::topology::NODE_RANK, metaA.locally_owned_part(), sourceNodes, sortById);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_entities(meshB, stk::topology::NODE_RANK, metaB.locally_owned_part(), targetNodes, sortById);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceIntField);
    sourceFields.push_back(&scalarSourceDoubleField);
    sourceFields.push_back(&vectorSourceDoubleField);
    stk::transfer::TransferCopyByIdStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetIntField);
    targetFields.push_back(&scalarTargetDoubleField);
    targetFields.push_back(&vectorTargetDoubleField);
    stk::transfer::TransferCopyByIdStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    {
      typedef stk::transfer::SearchById::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold;
      EXPECT_EQ(gold, key_to_target_processor);

      EXPECT_EQ( 8u, copySearch.get_remote_keys().size() );
    }
    stk::transfer::TransferCopyById transfer(copySearch, transferSource, transferTarget);

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

  stk::transfer::SearchByIdGeometric geometricSearch;
  stk::transfer::SearchByIdCommAll commAllSearch;
  stk::transfer::SearchById * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_EQ(copySearchPtr, &geometricSearch);
    }
    stk::transfer::SearchById & copySearch = *copySearchPtr;

    const size_t spatial_dimension = 3;
    const size_t num_elements = 3;
    const size_t num_nodes = 16;
    stk::mesh::EntityIdVector element_ids {1, 2, 3};
    std::vector<int> element_ownerA = {0, 0, -1};
    std::vector<int> element_ownerB = {-1, 1, 1};
    stk::mesh::EntityIdVector elem_node_ids[] {
        {1, 2, 6, 5, 9, 10, 14, 13},
        {2, 3, 7, 6, 10, 11, 15, 14},
        {3, 4, 8, 7, 11, 12, 16, 15} };
    std::vector<int> node_sharingA = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
    std::vector<int> node_sharingB = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
    std::vector<std::vector<double>> coordinates = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
                                                     {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0},
                                                     {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0},
                                                     {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0} };

    // Set up the "source" mesh for the transfer
    //
    std::shared_ptr<stk::mesh::BulkData> meshAPtr = build_mesh(spatial_dimension, pm);
    stk::mesh::MetaData& metaA = meshAPtr->mesh_meta_data();
    stk::mesh::BulkData& meshA = *meshAPtr;
    build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_ownerA, elem_node_ids, node_sharingA, coordinates);

    // Set up the "target" mesh for the transfer
    //
    std::shared_ptr<stk::mesh::BulkData> meshBPtr = build_mesh(spatial_dimension, pm);
    stk::mesh::MetaData& metaB = meshBPtr->mesh_meta_data();
    stk::mesh::BulkData& meshB = *meshBPtr;
    build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_ownerB, elem_node_ids, node_sharingB, coordinates);

    // Fill "source" fields with valid data
    //
    fill_mesh_values(meshA);

    ScalarIntField    & scalarSourceIntField    = static_cast<ScalarIntField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Int Field"));
    ScalarDoubleField & scalarSourceDoubleField = static_cast<ScalarDoubleField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Double Field"));
    VectorDoubleField & vectorSourceDoubleField = static_cast<VectorDoubleField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Double Field"));
    ScalarIntField    & scalarTargetIntField    = static_cast<ScalarIntField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Int Field"));
    ScalarDoubleField & scalarTargetDoubleField = static_cast<ScalarDoubleField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Double Field"));
    VectorDoubleField & vectorTargetDoubleField = static_cast<VectorDoubleField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Double Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    const bool sortById = true;
    stk::mesh::get_entities(meshA, stk::topology::NODE_RANK, metaA.locally_owned_part(), sourceNodes, sortById);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_entities(meshB, stk::topology::NODE_RANK, metaB.locally_owned_part(), targetNodes, sortById);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceIntField);
    sourceFields.push_back(&scalarSourceDoubleField);
    sourceFields.push_back(&vectorSourceDoubleField);
    stk::transfer::TransferCopyByIdStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetIntField);
    targetFields.push_back(&scalarTargetDoubleField);
    targetFields.push_back(&vectorTargetDoubleField);
    stk::transfer::TransferCopyByIdStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    //  GeometricTransfer
    //  stk::transfer::GeometricTransfer<
    //    class stk::transfer::LinearInterpolate<
    //      class stk::transfer::STKNode,
    //      class stk::transfer::STKNode
    //    >
    //  > transfer(transferSource, transferTarget, "copy01T10 unit test");

    {
      const int p_rank = stk::parallel_machine_rank( pm );
      typedef stk::transfer::SearchById::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold;
      if (0 == p_rank) {
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,2), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,6), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,10), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,14), 1);

        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,3), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,7), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,11), 1);
        gold.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,15), 1);
      }
      stk::util::sort_and_unique(gold);
      EXPECT_EQ(gold, key_to_target_processor);

      typedef stk::transfer::SearchById::MeshIDSet MeshIDSet;
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
      EXPECT_EQ(copySearch.get_remote_keys(), gold_remote_keys);
    }
    stk::transfer::TransferCopyById transfer(copySearch, transferSource, transferTarget);

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

  stk::transfer::SearchByIdGeometric geometricSearch;
  stk::transfer::SearchByIdCommAll commAllSearch;
  stk::transfer::SearchById * copySearchPtr = &commAllSearch;
  for (int search_index=0 ; search_index<2 ; ++search_index)
  {
    if (1 == search_index) {
      copySearchPtr = &geometricSearch;
      EXPECT_EQ(copySearchPtr, &geometricSearch);
    }
    stk::transfer::SearchById & copySearch = *copySearchPtr;

    const size_t spatial_dimension = 3;
    const size_t num_elements = 5;
    const size_t num_nodes = 24;
    stk::mesh::EntityIdVector element_ids {1, 2, 3, 4, 5};
    std::vector<int> element_ownerA = {0, 0, -1, -1, -1};
    std::vector<int> element_ownerB = {-1, -1, -1, 1, 1};
    stk::mesh::EntityIdVector elem_node_ids[] {
        {1, 2, 8, 7, 13, 14, 20, 19},
        {2, 3, 9, 8, 14, 15, 21, 20},
        {3, 4, 10, 9, 15, 16, 22, 21},
        {4, 5, 11, 10, 16, 17, 23, 22},
        {5, 6, 12, 11, 17, 18, 24, 23} };
    std::vector<int> node_sharingA = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
    std::vector<int> node_sharingB = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
    std::vector<std::vector<double>> coordinates = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {4.0, 0.0, 0.0}, {5.0, 0.0, 0.0},
      {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0}, {4.0, 1.0, 0.0}, {5.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {3.0, 0.0, 1.0}, {4.0, 0.0, 1.0}, {5.0, 0.0, 1.0},
      {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}, {3.0, 1.0, 1.0}, {4.0, 1.0, 1.0}, {5.0, 1.0, 1.0} };

    // Set up the "source" mesh for the transfer
    //
    std::shared_ptr<stk::mesh::BulkData> meshAPtr = build_mesh(spatial_dimension, pm);
    stk::mesh::MetaData& metaA = meshAPtr->mesh_meta_data();
    stk::mesh::BulkData& meshA = *meshAPtr;
    build_mesh(metaA, meshA, num_elements, num_nodes, element_ids, element_ownerA, elem_node_ids, node_sharingA, coordinates);

    // Set up the "target" mesh for the transfer
    //
    std::shared_ptr<stk::mesh::BulkData> meshBPtr = build_mesh(spatial_dimension, pm);
    stk::mesh::MetaData& metaB = meshBPtr->mesh_meta_data();
    stk::mesh::BulkData& meshB = *meshBPtr;
    build_mesh(metaB, meshB, num_elements, num_nodes, element_ids, element_ownerB, elem_node_ids, node_sharingB, coordinates);

    // Fill "source" fields with valid data
    //
    fill_mesh_values(meshA);

    ScalarIntField    & scalarSourceIntField    = static_cast<ScalarIntField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Int Field"));
    ScalarDoubleField & scalarSourceDoubleField = static_cast<ScalarDoubleField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Scalar Double Field"));
    VectorDoubleField & vectorSourceDoubleField = static_cast<VectorDoubleField&>(*metaA.get_field(stk::topology::NODE_RANK, "Node Vector Double Field"));
    ScalarIntField    & scalarTargetIntField    = static_cast<ScalarIntField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Int Field"));
    ScalarDoubleField & scalarTargetDoubleField = static_cast<ScalarDoubleField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Scalar Double Field"));
    VectorDoubleField & vectorTargetDoubleField = static_cast<VectorDoubleField&>(*metaB.get_field(stk::topology::NODE_RANK, "Node Vector Double Field"));

    // Set up the transfer
    //
    std::vector<stk::mesh::Entity> sourceNodes;
    const bool sortById = true;
    stk::mesh::get_entities(meshA, stk::topology::NODE_RANK, metaA.locally_owned_part(), sourceNodes, sortById);
    std::vector<stk::mesh::Entity> targetNodes;
    stk::mesh::get_entities(meshB, stk::topology::NODE_RANK, metaB.locally_owned_part(), targetNodes, sortById);

    std::vector<stk::mesh::FieldBase*> sourceFields;
    sourceFields.push_back(&scalarSourceIntField);
    sourceFields.push_back(&scalarSourceDoubleField);
    sourceFields.push_back(&vectorSourceDoubleField);
    stk::transfer::TransferCopyByIdStkMeshAdapter transferSource(meshA, sourceNodes, sourceFields);

    std::vector<stk::mesh::FieldBase*> targetFields;
    targetFields.push_back(&scalarTargetIntField);
    targetFields.push_back(&scalarTargetDoubleField);
    targetFields.push_back(&vectorTargetDoubleField);
    stk::transfer::TransferCopyByIdStkMeshAdapter transferTarget(meshB, targetNodes, targetFields);

    //  GeometricTransfer
    //  stk::transfer::GeometricTransfer<
    //    class stk::transfer::LinearInterpolate<
    //      class stk::transfer::STKNode,
    //      class stk::transfer::STKNode
    //    >
    //  > transfer(transferSource, transferTarget, "copy01T10 unit test");

    {
      const int p_rank = stk::parallel_machine_rank( pm );
      typedef stk::transfer::SearchById::KeyToTargetProcessor KeyToTargetProcessor;
      KeyToTargetProcessor key_to_target_processor;
      copySearch.do_search(transferSource,transferTarget,key_to_target_processor);

      KeyToTargetProcessor gold;
      EXPECT_EQ(gold, key_to_target_processor);

      typedef stk::transfer::SearchById::MeshIDSet MeshIDSet;
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
      EXPECT_EQ(copySearch.get_remote_keys(), gold_remote_keys);
    }
    stk::transfer::TransferCopyById transfer(copySearch, transferSource, transferTarget);

    // Do the transfer
    //
    transfer.initialize();
    EXPECT_THROW(transfer.apply(), std::runtime_error);

  }
}

TEST(Transfer, dataTypeTranslation)
{
  unsigned expectedLength = 10;

  for(stk::transfer::DataTypeKey::data_t i = stk::transfer::DataTypeKey::data_t::UNSIGNED_INTEGER; i <= stk::transfer::DataTypeKey::data_t::DOUBLE; 
      i = static_cast<stk::transfer::DataTypeKey::data_t>(i+1)) {
    stk::transfer::DataTypeKey dataTranslator(i, expectedLength);
    EXPECT_EQ(expectedLength, dataTranslator.get_data_length());
    EXPECT_EQ(i, dataTranslator.get_data_type());
  }
}

TEST(Transfer, dataLengthIsTooLong)
{
  unsigned expectedLength = stk::transfer::DataTypeKey::key_t::MAX_LENGTH + 1;

  for(stk::transfer::DataTypeKey::data_t i = stk::transfer::DataTypeKey::data_t::UNSIGNED_INTEGER; i <= stk::transfer::DataTypeKey::data_t::DOUBLE; 
      i = static_cast<stk::transfer::DataTypeKey::data_t>(i+1)) {
    EXPECT_THROW(stk::transfer::DataTypeKey dataTranslator(i, expectedLength), std::logic_error);
  }
}

template<typename SENDTYPE, typename RECVTYPE>
void test_mismatched_data_type_copy_transfer(unsigned sendBufferSize, unsigned recvBufferSize,
                                             stk::transfer::DataTypeKey::data_t sendType, stk::transfer::DataTypeKey::data_t recvType)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { return; }
  stk::CommSparse commSparse(MPI_COMM_WORLD);

  std::vector<SENDTYPE> sendBuffer;
  std::vector<RECVTYPE> recvBuffer;

  for(unsigned i = 0; i < sendBufferSize; i++) {
    sendBuffer.push_back(i);
  }
  for(unsigned i = 0; i < recvBufferSize; i++) {
    recvBuffer.push_back(std::numeric_limits<RECVTYPE>::max());
  }

  unsigned sendBufferDataSize = sendBufferSize * sizeof(SENDTYPE);
  unsigned recvBufferDataSize = recvBufferSize * sizeof(RECVTYPE);
  int communicatingProc = 1 - stk::parallel_machine_rank(MPI_COMM_WORLD);
  const uint8_t* sendData = reinterpret_cast<const uint8_t*>(sendBuffer.data());

  for(int phase = 0; phase < 2; phase++) {
    stk::transfer::DataTypeKey dataTranslator(sendType, sendBufferDataSize);
    commSparse.send_buffer(communicatingProc).pack<unsigned>(dataTranslator.m_value);
    for(unsigned index = 0; index < sendBufferDataSize; index++) {
      commSparse.send_buffer(communicatingProc).pack<uint8_t>(sendData[index]);
    }
    if(phase == 0) {
      commSparse.allocate_buffers();
    } else {
      commSparse.communicate();
    }
  }

  std::vector<uint8_t> tmpBuffer(sendBufferDataSize);
  unsigned newBufferSize = 0;
  unsigned dataKey = 0;
  commSparse.recv_buffer(communicatingProc).unpack<unsigned>(dataKey);

  stk::transfer::DataTypeKey unpackedData(dataKey);
  stk::transfer::DataTypeKey::data_t sentDataType = unpackedData.get_data_type();
  newBufferSize = unpackedData.get_data_length();

  EXPECT_EQ(sendType, sentDataType);
  EXPECT_EQ(sendBufferDataSize, newBufferSize);

  for(unsigned index = 0; index < sendBufferDataSize; index++) {
    commSparse.recv_buffer(communicatingProc).unpack<uint8_t>(tmpBuffer[index]);
  }

  stk::transfer::DataTypeTranslator<unsigned> translateUnsigned;
  stk::transfer::DataTypeTranslator<int64_t> translateInt64;
  stk::transfer::DataTypeTranslator<uint64_t> translateUInt64;
  stk::transfer::DataTypeTranslator<int> translateInt;
  stk::transfer::DataTypeTranslator<long double> translateLongDouble;
  stk::transfer::DataTypeTranslator<double> translateDouble;
  std::vector<stk::transfer::TranslatorBase*> dataTranslators = { &translateUnsigned, 
                                                                  &translateInt64,
                                                                  &translateUInt64,
                                                                  &translateInt,
                                                                  &translateLongDouble,
                                                                  &translateDouble };

  dataTranslators[sentDataType]->translate(tmpBuffer.data(), newBufferSize, recvType, recvBuffer.data(), recvBufferDataSize);

  unsigned numEntries = std::min(sendBufferSize, recvBufferSize);

  if(stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    for(unsigned i = 0; i < numEntries; i++) {
      EXPECT_EQ(recvBuffer[i], sendBuffer[i]);
    }
    for(unsigned i = numEntries; i < recvBufferSize; i++) {
      EXPECT_EQ(std::numeric_limits<RECVTYPE>::max(), recvBuffer[i]);
    }
  }
}

TEST(Transfer, mismatchedDataTypeCopyTransfer)
{
  using DataType = stk::transfer::DataTypeKey::data_t;
  test_mismatched_data_type_copy_transfer<unsigned, unsigned long>(2, 2, DataType::UNSIGNED_INTEGER, DataType::UNSIGNED_INTEGER_64);
  test_mismatched_data_type_copy_transfer<unsigned, unsigned long>(4, 2, DataType::UNSIGNED_INTEGER, DataType::UNSIGNED_INTEGER_64);
  test_mismatched_data_type_copy_transfer<unsigned, unsigned long>(2, 4, DataType::UNSIGNED_INTEGER, DataType::UNSIGNED_INTEGER_64);

  test_mismatched_data_type_copy_transfer<int, double>(2, 2, DataType::INTEGER, DataType::DOUBLE);
  test_mismatched_data_type_copy_transfer<int, double>(4, 2, DataType::INTEGER, DataType::DOUBLE);
  test_mismatched_data_type_copy_transfer<int, double>(2, 4, DataType::INTEGER, DataType::DOUBLE);

  test_mismatched_data_type_copy_transfer<double, unsigned>(2, 2, DataType::DOUBLE, DataType::UNSIGNED_INTEGER);
  test_mismatched_data_type_copy_transfer<double, unsigned>(4, 2, DataType::DOUBLE, DataType::UNSIGNED_INTEGER);
  test_mismatched_data_type_copy_transfer<double, unsigned>(2, 4, DataType::DOUBLE, DataType::UNSIGNED_INTEGER);

  test_mismatched_data_type_copy_transfer<uint64_t, double>(2, 2, DataType::UNSIGNED_INTEGER_64, DataType::DOUBLE);
  test_mismatched_data_type_copy_transfer<uint64_t, double>(4, 2, DataType::UNSIGNED_INTEGER_64, DataType::DOUBLE);
  test_mismatched_data_type_copy_transfer<uint64_t, double>(2, 4, DataType::UNSIGNED_INTEGER_64, DataType::DOUBLE);
}

namespace {
  stk::mesh::Part&
  create_block_part(stk::mesh::MetaData& meta, const stk::topology topology, const std::string& blockName, int64_t blockId)
  {
    stk::mesh::Part& part = meta.declare_part_with_topology(blockName, topology);
    stk::io::put_io_part_attribute(part);
    meta.set_part_id(part, blockId);
    return part;
  }
}

TEST(Transfer, mismatchedFieldDataTypeCopyTransfer)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) { return; }

  std::string fieldName = "node_field";
  int intInitVals = std::numeric_limits<int>::max();
  double doubleInitVals = std::numeric_limits<double>::max();
  std::vector<double> coords = {0,0, 1,0, 0,1, 1,1};

  std::shared_ptr<stk::mesh::BulkData> bulkAPtr = build_mesh(2, MPI_COMM_WORLD);
  stk::mesh::MetaData& metaA = bulkAPtr->mesh_meta_data();
  stk::mesh::BulkData& bulkA = *bulkAPtr;
  create_block_part(metaA, stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::FieldBase* fieldBaseA = &metaA.declare_field<int>(stk::topology::NODE_RANK, fieldName);
  stk::mesh::put_field_on_mesh(*fieldBaseA, metaA.universal_part(), &intInitVals);

  std::string meshDescA = "0,1,QUAD_4_2D,1,2,4,3,block_1";
  stk::unit_test_util::setup_text_mesh(bulkA, stk::unit_test_util::get_full_text_mesh_desc(meshDescA, coords));

  std::shared_ptr<stk::mesh::BulkData> bulkBPtr = build_mesh(2, MPI_COMM_WORLD);
  stk::mesh::MetaData& metaB = bulkBPtr->mesh_meta_data();
  stk::mesh::BulkData& bulkB = *bulkBPtr;
  create_block_part(metaB, stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::FieldBase* fieldBaseB = &metaB.declare_field<double>(stk::topology::NODE_RANK, fieldName);
  stk::mesh::put_field_on_mesh(*fieldBaseB, metaB.universal_part(), &doubleInitVals);

  std::string meshDescB;
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    meshDescB = "0,1,QUAD_4_2D,1,2,4,3,block_1";
  } else {
    meshDescB = "1,1,QUAD_4_2D,1,2,4,3,block_1";
  }
  stk::unit_test_util::setup_text_mesh(bulkB, stk::unit_test_util::get_full_text_mesh_desc(meshDescB, coords));

  // Set up CopyTransfer
  stk::mesh::EntityVector entitiesA;
  const bool sortById = true;
  stk::mesh::get_entities(bulkA, stk::topology::NODE_RANK, metaA.locally_owned_part(),
                          entitiesA, sortById);
  std::vector<stk::mesh::FieldBase*> fieldsA;
  fieldsA.push_back(fieldBaseA);
  stk::transfer::TransferCopyByIdStkMeshAdapter transferMeshA(bulkA,entitiesA,fieldsA);

  for(stk::mesh::Entity node : entitiesA) {
    int* scalar = reinterpret_cast<int*>(stk::mesh::field_data(*fieldBaseA, node));
    *scalar = static_cast<int>(0);
  }

  stk::mesh::EntityVector entitiesB;
  stk::mesh::get_entities(bulkB, stk::topology::NODE_RANK, metaB.locally_owned_part(),
                          entitiesB, sortById);
  std::vector<stk::mesh::FieldBase*> fieldsB;
  fieldsB.push_back(fieldBaseB);
  stk::transfer::TransferCopyByIdStkMeshAdapter transferMeshB(bulkB,entitiesB,fieldsB);

  for(stk::mesh::Entity node : entitiesB) {
    double* scalar = reinterpret_cast<double*>(stk::mesh::field_data(*fieldBaseB, node));
    *scalar = static_cast<double>(std::numeric_limits<double>::max());
  }

  stk::transfer::SearchByIdCommAll copySearch;

  stk::transfer::TransferCopyById copyTransfer(copySearch,transferMeshA,transferMeshB);
  copyTransfer.initialize();

  // Apply CopyTransfer
  copyTransfer.apply();

  if(bulkA.parallel_rank() == 1) {
    for(stk::mesh::Entity node : entitiesB) {
      double* scalar = reinterpret_cast<double*>(stk::mesh::field_data(*fieldBaseB, node));
      EXPECT_DOUBLE_EQ(0.0, *scalar);
    }
  }
}

} // namespace

