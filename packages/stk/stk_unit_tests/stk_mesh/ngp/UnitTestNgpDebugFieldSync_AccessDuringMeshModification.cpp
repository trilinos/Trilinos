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
#include <stk_mesh/base/Ngp.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/NgpFieldSyncDebugger.hpp>
#include <stk_util/stk_config.h>
#include "NgpDebugFieldSync_Fixtures.hpp"
#include <string>
#include <sstream>
#include <vector>
#include <numeric>

namespace {

class NgpDebugFieldSync_AccessDuringMeshModification : public NgpDebugFieldSyncFixture
{
public:
  template <typename T>
  void write_scalar_field_on_host_using_entity_vector(stk::mesh::Field<T> & stkField, stk::mesh::EntityVector entities, T value)
  {
    for (const stk::mesh::Entity & entity : entities) {
      if (get_bulk().is_valid(entity)) {
        const stk::mesh::MeshIndex & meshIndex = get_bulk().mesh_index(entity);
        const unsigned bucketId = meshIndex.bucket->bucket_id();
        const stk::mesh::Bucket::size_type bucketOrd = meshIndex.bucket_ordinal;
        T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, bucketId, bucketOrd);
        fieldData[0] = value;
      }
    }
  }

  template <typename T>
  void write_vector_field_on_host_using_bucket(stk::mesh::Field<T> & stkField, T value)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, *bucket);
      for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
        const size_t yComponent = iEntity*3 + 1;
        fieldData[yComponent] = value;
      }
    }
  }

  template <typename T>
  void read_vector_field_on_host_using_bucket(stk::mesh::Field<T> & stkField)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      const T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, *bucket);
      access_for_memory_checking_tool(fieldData, bucket->size()*3);
    }
  }

  template <typename T>
  void read_vector_field_on_device(stk::mesh::Field<T> & stkField)
  {
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::mesh::NgpField<T, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<T, NgpDebugger>(stkField);
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), meta.locally_owned_part());
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned i = 0; i < bucketIds.size(); ++i) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(i));
                             for (unsigned j = 0; j < bucket.size(); ++j) {
                               stk::mesh::FastMeshIndex index = ngpMesh.fast_mesh_index(bucket[j]);
                               unsigned numComponents = ngpField.get_num_components_per_entity(index);
                               for (unsigned component = 0; component < numComponents; ++component) {
#if defined(DEVICE_USE_LOCATION_BUILTINS)
                                 access_for_memory_checking_tool(&ngpField(index, component));
#else
                                 access_for_memory_checking_tool(&ngpField(index, component, __FILE__, __LINE__));
#endif
                               }
                             }
                           }
                         });
  }


  template <typename T>
  void create_element_with_scalar_field_write_using_entity(
           const std::vector<std::pair<stk::mesh::EntityId, std::string>> & elemParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(get_bulk(), counts);

    get_bulk().modification_begin();

    write_scalar_field_on_host_using_entity(stkField, value);

    for (const auto & elemPart : elemParts) {
      stk::mesh::PartVector parts {get_meta().get_part(elemPart.second), get_meta().get_part("block_1")};
      stk::mesh::EntityIdVector nodeIds;
      stk::mesh::EntityId nodeId = counts[stk::topology::ELEM_RANK] * 4 + 1;
      for (unsigned i = 0; i < 8; ++i) {
        nodeIds.push_back(nodeId++);
      }

      stk::mesh::Entity newElem = stk::mesh::declare_element(get_bulk(), parts, elemPart.first, nodeIds);
      T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, newElem);
      fieldData[0] = value;
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void create_element_with_scalar_field_write_using_bucket(
           const std::vector<std::pair<stk::mesh::EntityId, std::string>> & elemParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(get_bulk(), counts);

    get_bulk().modification_begin();

    write_scalar_field_on_host_using_bucket(stkField, value);

    for (const auto & elemPart : elemParts) {
      stk::mesh::PartVector parts {get_meta().get_part(elemPart.second), get_meta().get_part("block_1")};
      stk::mesh::EntityIdVector nodeIds;
      stk::mesh::EntityId nodeId = counts[stk::topology::ELEM_RANK] * 4 + 1;
      for (unsigned i = 0; i < 8; ++i) {
        nodeIds.push_back(nodeId++);
      }

      stk::mesh::Entity newElem = stk::mesh::declare_element(get_bulk(), parts, elemPart.first, nodeIds);

      const stk::mesh::Bucket & bucket = get_bulk().bucket(newElem);
      T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, bucket);
      for(size_t iEntity = 0; iEntity < bucket.size(); ++iEntity) {
        fieldData[iEntity] = value;
      }
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void create_element_with_scalar_field_write_using_bucket_id(
           const std::vector<std::pair<stk::mesh::EntityId, std::string>> & elemParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(get_bulk(), counts);

    get_bulk().modification_begin();

    write_scalar_field_on_host_using_bucket_id(stkField, value);

    for (const auto & elemPart : elemParts) {
      stk::mesh::PartVector parts {get_meta().get_part(elemPart.second), get_meta().get_part("block_1")};
      stk::mesh::EntityIdVector nodeIds;
      stk::mesh::EntityId nodeId = counts[stk::topology::ELEM_RANK] * 4 + 1;
      for (unsigned i = 0; i < 8; ++i) {
        nodeIds.push_back(nodeId++);
      }

      stk::mesh::Entity newElem = stk::mesh::declare_element(get_bulk(), parts, elemPart.first, nodeIds);

      const stk::mesh::Bucket & bucket = get_bulk().bucket(newElem);
      const unsigned bucketId = bucket.bucket_id();
      T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, bucketId);
      for(size_t iEntity = 0; iEntity < bucket.size(); ++iEntity) {
        fieldData[iEntity] = value;
      }
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void create_element_with_scalar_field_write_using_bucket_id_and_ordinal(
           const std::vector<std::pair<stk::mesh::EntityId, std::string>> & elemParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(get_bulk(), counts);

    get_bulk().modification_begin();

    write_scalar_field_on_host_using_bucket_id_and_ordinal(stkField, value);

    for (const auto & elemPart : elemParts) {
      stk::mesh::PartVector parts {get_meta().get_part(elemPart.second), get_meta().get_part("block_1")};
      stk::mesh::EntityIdVector nodeIds;
      stk::mesh::EntityId nodeId = counts[stk::topology::ELEM_RANK] * 4 + 1;
      for (unsigned i = 0; i < 8; ++i) {
        nodeIds.push_back(nodeId++);
      }

      stk::mesh::Entity newElem = stk::mesh::declare_element(get_bulk(), parts, elemPart.first, nodeIds);

      const stk::mesh::MeshIndex & meshIndex = get_bulk().mesh_index(newElem);
      const unsigned bucketId = meshIndex.bucket->bucket_id();
      const stk::mesh::Bucket::size_type bucketOrd = meshIndex.bucket_ordinal;
      T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, bucketId, bucketOrd);
      fieldData[0] = value;
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void create_element_with_vector_field_write_using_entity(
           const std::vector<std::pair<stk::mesh::EntityId, std::string>> & elemParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(get_bulk(), counts);

    get_bulk().modification_begin();

    write_vector_field_on_host_using_entity(stkField, value);

    for (const auto & elemPart : elemParts) {
      stk::mesh::PartVector parts {get_meta().get_part(elemPart.second), get_meta().get_part("block_1")};
      stk::mesh::EntityIdVector nodeIds;
      stk::mesh::EntityId nodeId = counts[stk::topology::ELEM_RANK] * 4 + 1;
      for (unsigned i = 0; i < 8; ++i) {
        nodeIds.push_back(nodeId++);
      }

      stk::mesh::Entity newElem = stk::mesh::declare_element(get_bulk(), parts, elemPart.first, nodeIds);
      T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, newElem);
      fieldData[1] = value;
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void create_element_with_vector_field_write_using_bucket(
           const std::vector<std::pair<stk::mesh::EntityId, std::string>> & elemParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(get_bulk(), counts);

    get_bulk().modification_begin();

    write_vector_field_on_host_using_bucket(stkField, value);

    for (const auto & elemPart : elemParts) {
      stk::mesh::PartVector parts {get_meta().get_part(elemPart.second), get_meta().get_part("block_1")};
      stk::mesh::EntityIdVector nodeIds;
      stk::mesh::EntityId nodeId = counts[stk::topology::ELEM_RANK] * 4 + 1;
      for (unsigned i = 0; i < 8; ++i) {
        nodeIds.push_back(nodeId++);
      }

      stk::mesh::Entity newElem = stk::mesh::declare_element(get_bulk(), parts, elemPart.first, nodeIds);

      const stk::mesh::Bucket & bucket = get_bulk().bucket(newElem);
      T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, bucket);
      for(size_t iEntity = 0; iEntity < bucket.size(); ++iEntity) {
        const size_t yComponent = iEntity*3 + 1;
        fieldData[yComponent] = value;
      }
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void create_element_with_scalar_field_read_using_entity(
           const std::vector<std::pair<stk::mesh::EntityId, std::string>> & elemParts,
           stk::mesh::Field<T> & stkField)
  {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(get_bulk(), counts);

    get_bulk().modification_begin();

    for (const auto & elemPart : elemParts) {
      stk::mesh::PartVector parts {get_meta().get_part(elemPart.second), get_meta().get_part("block_1")};
      stk::mesh::EntityIdVector nodeIds;
      stk::mesh::EntityId nodeId = counts[stk::topology::ELEM_RANK] * 4 + 1;
      for (unsigned i = 0; i < 8; ++i) {
        nodeIds.push_back(nodeId++);
      }

      stk::mesh::declare_element(get_bulk(), parts, elemPart.first, nodeIds);
    }

    read_scalar_field_on_host_using_entity(stkField);

    get_bulk().modification_end();
  }

  template <typename T>
  void create_element_with_scalar_field_read_using_bucket(
           const std::vector<std::pair<stk::mesh::EntityId, std::string>> & elemParts,
           stk::mesh::Field<T> & stkField)
  {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(get_bulk(), counts);

    get_bulk().modification_begin();

    for (const auto & elemPart : elemParts) {
      stk::mesh::PartVector parts {get_meta().get_part(elemPart.second), get_meta().get_part("block_1")};
      stk::mesh::EntityIdVector nodeIds;
      stk::mesh::EntityId nodeId = counts[stk::topology::ELEM_RANK] * 4 + 1;
      for (unsigned i = 0; i < 8; ++i) {
        nodeIds.push_back(nodeId++);
      }

      stk::mesh::declare_element(get_bulk(), parts, elemPart.first, nodeIds);
    }

    read_scalar_field_on_host_using_bucket(stkField);

    get_bulk().modification_end();
  }

  template <typename T>
  void create_element_with_vector_field_read_using_entity(
           const std::vector<std::pair<stk::mesh::EntityId, std::string>> & elemParts,
           stk::mesh::Field<T> & stkField)
  {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(get_bulk(), counts);

    get_bulk().modification_begin();

    for (const auto & elemPart : elemParts) {
      stk::mesh::PartVector parts {get_meta().get_part(elemPart.second), get_meta().get_part("block_1")};
      stk::mesh::EntityIdVector nodeIds;
      stk::mesh::EntityId nodeId = counts[stk::topology::ELEM_RANK] * 4 + 1;
      for (unsigned i = 0; i < 8; ++i) {
        nodeIds.push_back(nodeId++);
      }

      stk::mesh::declare_element(get_bulk(), parts, elemPart.first, nodeIds);
    }

    read_vector_field_on_host_using_entity<T>(stkField);

    get_bulk().modification_end();
  }

  template <typename T>
  void create_element_with_vector_field_read_using_bucket(
           const std::vector<std::pair<stk::mesh::EntityId, std::string>> & elemParts,
           stk::mesh::Field<T> & stkField)
  {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(get_bulk(), counts);

    get_bulk().modification_begin();

    for (const auto & elemPart : elemParts) {
      stk::mesh::PartVector parts {get_meta().get_part(elemPart.second), get_meta().get_part("block_1")};
      stk::mesh::EntityIdVector nodeIds;
      stk::mesh::EntityId nodeId = counts[stk::topology::ELEM_RANK] * 4 + 1;
      for (unsigned i = 0; i < 8; ++i) {
        nodeIds.push_back(nodeId++);
      }

      stk::mesh::declare_element(get_bulk(), parts, elemPart.first, nodeIds);
    }

    read_vector_field_on_host_using_bucket(stkField);

    get_bulk().modification_end();
  }


  template <typename T>
  void delete_element_with_scalar_field_write_using_entity(const std::vector<stk::mesh::EntityId> & elemIds,
                                                           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_scalar_field_on_host_using_entity(stkField, value);

    for (const stk::mesh::EntityId & elemId : elemIds) {
      get_bulk().destroy_entity(get_bulk().get_entity(stk::topology::ELEM_RANK, elemId));
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void delete_element_with_scalar_field_write_using_bucket(const std::vector<stk::mesh::EntityId> & elemIds,
                                                           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_scalar_field_on_host_using_bucket(stkField, value);

    for (const stk::mesh::EntityId & elemId : elemIds) {
      get_bulk().destroy_entity(get_bulk().get_entity(stk::topology::ELEM_RANK, elemId));
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void delete_element_with_scalar_field_write_using_bucket_id(const std::vector<stk::mesh::EntityId> & elemIds,
                                                              stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_scalar_field_on_host_using_bucket_id(stkField, value);

    for (const stk::mesh::EntityId & elemId : elemIds) {
      get_bulk().destroy_entity(get_bulk().get_entity(stk::topology::ELEM_RANK, elemId));
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void delete_element_with_scalar_field_write_using_bucket_id_and_ordinal(const std::vector<stk::mesh::EntityId> & elemIds,
                                                                          stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    stk::mesh::EntityVector processedEntities = write_scalar_field_on_host_using_bucket_id_and_ordinal(stkField, value);

    for (const stk::mesh::EntityId & elemId : elemIds) {
      get_bulk().destroy_entity(get_bulk().get_entity(stk::topology::ELEM_RANK, elemId));
    }

    write_scalar_field_on_host_using_entity_vector(stkField, processedEntities, value);


    get_bulk().modification_end();
  }

  template <typename T>
  void delete_element_with_vector_field_write_using_entity(const std::vector<stk::mesh::EntityId> & elemIds,
                                                           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_vector_field_on_host_using_entity(stkField, value);

    for (const stk::mesh::EntityId & elemId : elemIds) {
      get_bulk().destroy_entity(get_bulk().get_entity(stk::topology::ELEM_RANK, elemId));
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void delete_element_with_vector_field_write_using_bucket(const std::vector<stk::mesh::EntityId> & elemIds,
                                                           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_vector_field_on_host_using_bucket(stkField, value);

    for (const stk::mesh::EntityId & elemId : elemIds) {
      get_bulk().destroy_entity(get_bulk().get_entity(stk::topology::ELEM_RANK, elemId));
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void delete_element_with_scalar_field_read_using_entity(const std::vector<stk::mesh::EntityId> & elemIds,
                                                          stk::mesh::Field<T> & stkField)
  {
    get_bulk().modification_begin();

    for (const stk::mesh::EntityId & elemId : elemIds) {
      get_bulk().destroy_entity(get_bulk().get_entity(stk::topology::ELEM_RANK, elemId));
    }

    read_scalar_field_on_host_using_entity(stkField);

    get_bulk().modification_end();
  }

  template <typename T>
  void delete_element_with_scalar_field_read_using_bucket(const std::vector<stk::mesh::EntityId> & elemIds,
                                                          stk::mesh::Field<T> & stkField)
  {
    get_bulk().modification_begin();

    for (const stk::mesh::EntityId & elemId : elemIds) {
      get_bulk().destroy_entity(get_bulk().get_entity(stk::topology::ELEM_RANK, elemId));
    }

    read_scalar_field_on_host_using_bucket(stkField);

    get_bulk().modification_end();
  }

  template <typename T>
  void delete_element_with_vector_field_read_using_entity(const std::vector<stk::mesh::EntityId> & elemIds,
                                                          stk::mesh::Field<T> & stkField)
  {
    get_bulk().modification_begin();

    for (const stk::mesh::EntityId & elemId : elemIds) {
      get_bulk().destroy_entity(get_bulk().get_entity(stk::topology::ELEM_RANK, elemId));
    }

    read_vector_field_on_host_using_entity<T>(stkField);

    get_bulk().modification_end();
  }

  template <typename T>
  void delete_element_with_vector_field_read_using_bucket(const std::vector<stk::mesh::EntityId> & elemIds,
                                                          stk::mesh::Field<T> & stkField)
  {
    get_bulk().modification_begin();

    for (const stk::mesh::EntityId & elemId : elemIds) {
      get_bulk().destroy_entity(get_bulk().get_entity(stk::topology::ELEM_RANK, elemId));
    }

    read_vector_field_on_host_using_bucket(stkField);

    get_bulk().modification_end();
  }


  template <typename T>
  void modify_element_part_membership_with_scalar_field_write_using_bucket_id_and_ordinal(
           const std::vector<EntityIdAddRemovePart> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    stk::mesh::EntityVector processedEntities = write_scalar_field_on_host_using_bucket_id_and_ordinal(stkField, value);

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, elemAddRemovePart.id)};
      stk::mesh::PartVector addParts {get_meta().get_part(elemAddRemovePart.addPart)};
      stk::mesh::PartVector removeParts {get_meta().get_part(elemAddRemovePart.removePart)};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    write_scalar_field_on_host_using_entity_vector(stkField, processedEntities, value);

    get_bulk().modification_end();
  }

  template <typename T>
  void modify_element_part_membership_with_scalar_field_write_using_bucket(
           const std::vector<EntityIdAddRemovePart> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_scalar_field_on_host_using_bucket(stkField, value);

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, elemAddRemovePart.id)};
      stk::mesh::PartVector addParts {get_meta().get_part(elemAddRemovePart.addPart)};
      stk::mesh::PartVector removeParts {get_meta().get_part(elemAddRemovePart.removePart)};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void modify_element_part_membership_with_scalar_field_write_using_bucket_id(
           const std::vector<EntityIdAddRemovePart> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_scalar_field_on_host_using_bucket_id(stkField, value);

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, elemAddRemovePart.id)};
      stk::mesh::PartVector addParts {get_meta().get_part(elemAddRemovePart.addPart)};
      stk::mesh::PartVector removeParts {get_meta().get_part(elemAddRemovePart.removePart)};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void modify_element_part_membership_with_vector_field_write_using_entity(
           const std::vector<EntityIdAddRemovePart> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_vector_field_on_host_using_entity(stkField, value);

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, elemAddRemovePart.id)};
      stk::mesh::PartVector addParts {get_meta().get_part(elemAddRemovePart.addPart)};
      stk::mesh::PartVector removeParts {get_meta().get_part(elemAddRemovePart.removePart)};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void modify_element_part_membership_with_vector_field_write_using_bucket(
           const std::vector<EntityIdAddRemovePart> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_vector_field_on_host_using_bucket(stkField, value);

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, elemAddRemovePart.id)};
      stk::mesh::PartVector addParts {get_meta().get_part(elemAddRemovePart.addPart)};
      stk::mesh::PartVector removeParts {get_meta().get_part(elemAddRemovePart.removePart)};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    get_bulk().modification_end();
  }

  template <typename T>
  void modify_element_part_membership_with_scalar_field_read_using_entity(
           const std::vector<EntityIdAddRemovePart> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField)
  {
    get_bulk().modification_begin();

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, elemAddRemovePart.id)};
      stk::mesh::PartVector addParts {get_meta().get_part(elemAddRemovePart.addPart)};
      stk::mesh::PartVector removeParts {get_meta().get_part(elemAddRemovePart.removePart)};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    read_scalar_field_on_host_using_entity(stkField);

    get_bulk().modification_end();
  }

  template <typename T>
  void modify_element_part_membership_with_scalar_field_read_using_bucket(
           const std::vector<EntityIdAddRemovePart> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField)
  {
    get_bulk().modification_begin();

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, elemAddRemovePart.id)};
      stk::mesh::PartVector addParts {get_meta().get_part(elemAddRemovePart.addPart)};
      stk::mesh::PartVector removeParts {get_meta().get_part(elemAddRemovePart.removePart)};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    read_scalar_field_on_host_using_bucket(stkField);

    get_bulk().modification_end();
  }

  template <typename T>
  void modify_element_part_membership_with_vector_field_read_using_entity(
           const std::vector<EntityIdAddRemovePart> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField)
  {
    get_bulk().modification_begin();

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, elemAddRemovePart.id)};
      stk::mesh::PartVector addParts {get_meta().get_part(elemAddRemovePart.addPart)};
      stk::mesh::PartVector removeParts {get_meta().get_part(elemAddRemovePart.removePart)};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    read_vector_field_on_host_using_entity<T>(stkField);

    get_bulk().modification_end();
  }

  template <typename T>
  void modify_element_part_membership_with_vector_field_read_using_bucket(
           const std::vector<EntityIdAddRemovePart> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField)
  {
    get_bulk().modification_begin();

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, elemAddRemovePart.id)};
      stk::mesh::PartVector addParts {get_meta().get_part(elemAddRemovePart.addPart)};
      stk::mesh::PartVector removeParts {get_meta().get_part(elemAddRemovePart.removePart)};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    read_vector_field_on_host_using_bucket(stkField);

    get_bulk().modification_end();
  }
};

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_entity({{3, "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucketIdAndOrdinal_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_scalar_field_write_using_bucket_id_and_ordinal({{2, "Part2", "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucketIdAndOrdinal_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_bucket_id_and_ordinal({{3, "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucketIdAndOrdinal_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_bucket_id_and_ordinal({2}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingEntity_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_vector_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingEntity_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  create_element_with_vector_field_write_using_entity({{3, "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingEntity_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  delete_element_with_vector_field_write_using_entity({2}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucket_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_scalar_field_write_using_bucket({{2, "Part2", "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucket_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_bucket({{3, "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucket_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_bucket({2}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucketId_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_scalar_field_write_using_bucket_id({{2, "Part2", "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucketId_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_bucket_id({{3, "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucketId_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_bucket_id({2}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingBucket_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_vector_field_write_using_bucket({{2, "Part2", "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingBucket_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  create_element_with_vector_field_write_using_bucket({{3, "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingBucket_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  delete_element_with_vector_field_write_using_bucket({2}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_ChangeBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_CreateBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_DeleteBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingEntity_ChangeBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_vector_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 3.14);
  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingEntity_CreateBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  create_element_with_vector_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingEntity_DeleteBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  delete_element_with_vector_field_write_using_entity({2}, stkField, 3.14);
  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_bucket({{2, "Part2", "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucket_CreateBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  create_element_with_scalar_field_write_using_bucket({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucket_DeleteBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  delete_element_with_scalar_field_write_using_bucket({2}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucketId_ChangeBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_bucket_id({{2, "Part2", "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucketId_CreateBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  create_element_with_scalar_field_write_using_bucket_id({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucketId_DeleteBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  delete_element_with_scalar_field_write_using_bucket_id({2}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_vector_field_write_using_bucket({{2, "Part2", "Part1"}}, stkField, 3.14);
  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingBucket_CreateBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  create_element_with_vector_field_write_using_bucket({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingBucket_DeleteBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  delete_element_with_vector_field_write_using_bucket({2}, stkField, 3.14);
  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_ChangeBucket_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 3.14);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_CreateBucket_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);

  const stk::mesh::EntityId maxEntityIdInOldField = 1;
  read_old_scalar_field_on_device(stkField, ngpField, maxEntityIdInOldField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_DeleteBucket_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_ModifyBucket_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 3.14);

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_old_scalar_field_on_device(stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_CreateBucket_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  stkField.modify_on_host();
  stkField.sync_to_device();

  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation
  read_old_scalar_field_on_device(stkField, ngpFieldCopy, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_DeleteBucket_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_old_scalar_field_on_device(stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership_with_scalar_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element_with_scalar_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element_with_scalar_field_read_using_entity({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingEntity_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership_with_vector_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingEntity_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element_with_vector_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingEntity_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element_with_vector_field_read_using_entity({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucket_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership_with_scalar_field_read_using_bucket({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucket_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element_with_scalar_field_read_using_bucket({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucket_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element_with_scalar_field_read_using_bucket({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingBucket_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership_with_vector_field_read_using_bucket({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingBucket_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element_with_vector_field_read_using_bucket({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingBucket_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element_with_vector_field_read_using_bucket({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership_with_scalar_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  create_element_with_scalar_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  delete_element_with_scalar_field_read_using_entity({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingEntity_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  modify_element_part_membership_with_vector_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Lost Device values for Field doubleVectorField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingEntity_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  create_element_with_vector_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleVectorField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingEntity_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  delete_element_with_vector_field_read_using_entity({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleVectorField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucket_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership_with_scalar_field_read_using_bucket({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucket_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  create_element_with_scalar_field_read_using_bucket({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingBucket_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  delete_element_with_scalar_field_read_using_bucket({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingBucket_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  modify_element_part_membership_with_vector_field_read_using_bucket({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Lost Device values for Field doubleVectorField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingBucket_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  create_element_with_vector_field_read_using_bucket({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleVectorField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, VectorAccessUsingBucket_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_vector_field<double>("doubleVectorField", 3, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  delete_element_with_vector_field_read_using_bucket({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleVectorField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{3, "Part2", "Part1"}}, stkField, 3.14);
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 2.18);

  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  create_element_with_scalar_field_write_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField, 2.18);

  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  delete_element_with_scalar_field_write_using_entity({3}, stkField, 2.18);

  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{3, "Part2", "Part1"}}, stkField, 3.14);
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  create_element_with_scalar_field_write_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_MissingAllModifySyncCalls_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  delete_element_with_scalar_field_write_using_entity({3}, stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{3, "Part2", "Part1"}}, stkField, 3.14);
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 2.18);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  create_element_with_scalar_field_write_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField, 2.18);

  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation
  read_old_scalar_field_on_device(stkField, ngpField, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  delete_element_with_scalar_field_write_using_entity({3}, stkField, 2.18);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership_with_scalar_field_read_using_entity({{3, "Part2", "Part1"}}, stkField);
  modify_element_part_membership_with_scalar_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element_with_scalar_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);
  create_element_with_scalar_field_read_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element_with_scalar_field_read_using_entity({2}, stkField);
  delete_element_with_scalar_field_read_using_entity({3}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  modify_element_part_membership_with_scalar_field_read_using_entity({{3, "Part2", "Part1"}}, stkField);
  modify_element_part_membership_with_scalar_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 8, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  create_element_with_scalar_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);
  create_element_with_scalar_field_read_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 10, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  delete_element_with_scalar_field_read_using_entity({2}, stkField);
  delete_element_with_scalar_field_read_using_entity({3}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 5, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}



TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoMods_ChangeBucket_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{3, "Part2", "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoMods_CreateBucket_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  create_element_with_scalar_field_write_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, TwoMods_DeleteBucket_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  delete_element_with_scalar_field_write_using_entity({3}, stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}



TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_TwoMods_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{3, "Part2", "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 2.18);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_TwoMods_CreateBucket_CreateBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  create_element_with_scalar_field_write_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField, 2.18);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_TwoMods_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  delete_element_with_scalar_field_write_using_entity({3}, stkField, 2.18);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}



TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_TwoMods_ChangeBucket_ChangeBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{3, "Part2", "Part1"}}, stkField, 3.14);
  read_old_scalar_field_on_device(stkField, ngpField);

  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 2.18);
  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 8, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_TwoMods_CreateBucket_CreateBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation

  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  read_old_scalar_field_on_device(stkField, ngpField, maxIdToRead);

  create_element_with_scalar_field_write_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField, 2.18);
  read_old_scalar_field_on_device(stkField, ngpField, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_TwoMods_DeleteBucket_DeleteBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();

  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  read_old_scalar_field_on_device(stkField, ngpField);

  delete_element_with_scalar_field_write_using_entity({3}, stkField, 2.18);
  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 5, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}



TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_TwoMods_ChangeBucket_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();
  modify_element_part_membership_with_scalar_field_read_using_entity({{3, "Part2", "Part1"}}, stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();
  modify_element_part_membership_with_scalar_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_TwoMods_CreateBucket_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();
  create_element_with_scalar_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();
  create_element_with_scalar_field_read_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_TwoMods_DeleteBucket_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();
  delete_element_with_scalar_field_read_using_entity({2}, stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();
  delete_element_with_scalar_field_read_using_entity({3}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_TwoMods_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership_with_scalar_field_read_using_entity({{3, "Part2", "Part1"}}, stkField);

  write_scalar_field_on_device(stkField, 2.18);
  modify_element_part_membership_with_scalar_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 8, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_TwoMods_CreateBucket_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  create_element_with_scalar_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);

  write_scalar_field_on_device(stkField, 2.18);
  create_element_with_scalar_field_read_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 10, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_AccessDuringMeshModification, ScalarAccessUsingEntity_TwoMods_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  delete_element_with_scalar_field_read_using_entity({2}, stkField);

  write_scalar_field_on_device(stkField, 2.18);
  delete_element_with_scalar_field_read_using_entity({3}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 5, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

}
