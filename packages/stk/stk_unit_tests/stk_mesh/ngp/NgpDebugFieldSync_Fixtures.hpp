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

#ifndef NGPDEBUGFIELDSYNCFIXTURES_HPP
#define NGPDEBUGFIELDSYNCFIXTURES_HPP

#include <gtest/gtest.h>
#include <stk_mesh/base/NgpFieldSyncDebugger.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <string>

template <typename T> using NgpDebugger = stk::mesh::NgpFieldSyncDebugger<T>;
template <typename T> using StkDebugger = typename NgpDebugger<T>::StkFieldSyncDebuggerType;

void extract_warning(std::string & stdoutString, int numExpectedOccurrences, const std::string & warningString);

void check_no_warnings(const std::string & stdoutString);

void check_contains_file_name(const std::string & stdoutString, const std::string & fileName);

void check_contains_a_line_number(const std::string & stdoutString);

template <typename T>
STK_INLINE_FUNCTION
void access_for_memory_checking_tool(T valuePtr, unsigned numValues = 1)
{
  bool temp = true;
  for (unsigned i = 0; i < numValues; ++i) {
    temp &= (valuePtr[i] == 0);
  }
}

struct EntityIdAddRemovePart {
  stk::mesh::EntityId id;
  std::string addPart;
  std::string removePart;
};

class NgpDebugFieldSyncFixture : public stk::unit_test_util::MeshFixture
{
public:
  template <typename T>
  void declare_scalar_field(const std::string & fieldName,
                            const std::vector<std::string> & partsForField)
  {
    stk::mesh::Selector fieldParts = create_parts(partsForField);
    create_scalar_field<T>(fieldName, stk::topology::ELEM_RANK, fieldParts);
  }

  template <typename T>
  void declare_vector_field(const std::string & fieldName,
                            unsigned numComponents,
                            const std::vector<std::string> & partsForField)
  {
    stk::mesh::Selector fieldParts = create_parts(partsForField);
    create_vector_field<T>(fieldName, stk::topology::ELEM_RANK, numComponents, fieldParts);
  }

  stk::mesh::Selector
  create_parts(const std::vector<std::string> & partNames)
  {
    stk::mesh::Selector allParts;
    for (const std::string & partName : partNames) {
      allParts |= get_meta().declare_part_with_topology(partName, stk::topology::HEX_8);
    }
    return allParts;
  }

  template <typename T>
  stk::mesh::Field<T> & create_scalar_field(const std::string & name,
                                            stk::topology::rank_t rank,
                                            stk::mesh::Selector & fieldParts)
  {
    unsigned numStates = 1;
    const T init = 1;
    stk::mesh::Field<T> & field = get_meta().declare_field<stk::mesh::Field<T>>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, fieldParts, &init);
    return field;
  }

  template <typename T>
  stk::mesh::Field<T,stk::mesh::Cartesian> &
                        create_vector_field(const std::string & name,
                                            stk::topology::rank_t rank,
                                            unsigned numComponents,
                                            stk::mesh::Selector & fieldParts)
  {
    unsigned numStates = 1;
    const std::vector<T> init(numComponents, 1);
    stk::mesh::Field<T,stk::mesh::Cartesian> & field = get_meta().declare_field<stk::mesh::Field<T,stk::mesh::Cartesian>>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, fieldParts, numComponents, init.data());
    return field;
  }

  struct PartConfiguration {
    std::string partName;
    unsigned numElements;
  };

  void build_mesh(const std::vector<PartConfiguration> & partConfiguration,
                  unsigned bucketCapacity = stk::mesh::impl::BucketRepository::default_bucket_capacity)
  {
    std::vector<std::string> allMeshPartNames;
    unsigned numElems = 0;
    for (const auto & config : partConfiguration) {
      allMeshPartNames.push_back(config.partName);
      numElems += config.numElements;
    }

    create_parts(allMeshPartNames);
    setup_mesh("generated:1x1x" + std::to_string(numElems), stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);

    set_initial_part_membership(partConfiguration);
  }

  void set_initial_part_membership(const std::vector<PartConfiguration> & partConfiguration)
  {
    get_bulk().modification_begin();
    stk::mesh::EntityId elemId = 0;
    for (const auto & config : partConfiguration) {
      stk::mesh::PartVector addParts {get_meta().get_part(config.partName)};
      stk::mesh::EntityVector elemsToChange;
      for (unsigned elemNum = 0; elemNum < config.numElements; ++elemNum) {
        stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, ++elemId);
        ThrowRequireMsg(get_bulk().is_valid(element), "Invalid element in fixture!");
        elemsToChange.push_back(element);
      }
      get_bulk().change_entity_parts(elemsToChange, addParts);
    }
    get_bulk().modification_end();
  }

  template <typename T>
  stk::mesh::Field<T> & initialized_field(const std::string & fieldName)
  {
    stk::mesh::Field<T> & stkField = *static_cast<stk::mesh::Field<T>*>(get_meta().get_field(stk::topology::ELEM_RANK, fieldName));
    fill_initial_field<T>(stkField);
    initialize_ngp_field<T>(stkField);
    return stkField;
  }

  template <typename T>
  stk::mesh::Field<T,stk::mesh::Cartesian> & initialized_vector_field(const std::string & fieldName)
  {
    stk::mesh::Field<T,stk::mesh::Cartesian> & stkField = *static_cast<stk::mesh::Field<T,stk::mesh::Cartesian>*>(get_meta().get_field(stk::topology::ELEM_RANK, fieldName));
    fill_initial_field<T>(stkField);
    initialize_ngp_field<T>(stkField);
    return stkField;
  }

  template <typename T>
  void fill_initial_field(stk::mesh::FieldBase & stkField)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        const stk::mesh::EntityId id = get_bulk().identifier(entity);
        T * fieldData = static_cast<T*>(stk::mesh::field_data<stk::mesh::FieldBase, stk::mesh::EmptyStkFieldSyncDebugger>(stkField, entity));
        const unsigned numComponents = stk::mesh::field_scalars_per_entity(stkField, *bucket);
        for (unsigned component = 0; component < numComponents; ++component) {
          fieldData[component] = 10*id + component;
        }
      }
    }
  }

  template <typename T>
  void initialize_ngp_field(stk::mesh::Field<T> & stkField)
  {
    stk::mesh::get_updated_ngp_field<T, NgpDebugger>(stkField);
  }

  template <typename T>
  void initialize_ngp_field(stk::mesh::FieldBase & stkField)
  {
    stk::mesh::get_updated_ngp_field<T, NgpDebugger>(stkField);
  }

  void change_element_parts(const std::vector<EntityIdAddRemovePart>& elemAddRemoveParts)
  {
    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, elemAddRemovePart.id)};
      stk::mesh::PartVector addParts {get_meta().get_part(elemAddRemovePart.addPart)};
      stk::mesh::PartVector removeParts {get_meta().get_part(elemAddRemovePart.removePart)};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }
  }

  void modify_element_part_membership(const std::vector<EntityIdAddRemovePart>& elemAddRemoveParts)
  {
    get_bulk().modification_begin();
    change_element_parts(elemAddRemoveParts);
    get_bulk().modification_end();
  }

  template <typename T>
  void modify_element_part_membership_with_scalar_field_write_using_entity(
           const std::vector<EntityIdAddRemovePart> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_scalar_field_on_host_using_entity(stkField, value);

    change_element_parts(elemAddRemoveParts);

    get_bulk().modification_end();
  }

  template <typename T>
  void write_scalar_field_on_host_using_entity(stk::mesh::Field<T> & stkField, T value)
  {
    write_scalar_field_on_host_using_entity(stkField, stkField, value);
  }

  template <typename T>
  void write_vector_field_on_host_using_entity(stk::mesh::FieldBase & stkField, const stk::mesh::Selector& selector, T value)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().get_buckets(stkField.entity_rank(), selector);
    for (stk::mesh::Bucket * bucket : buckets) {
      unsigned numComponents = stk::mesh::field_scalars_per_entity(stkField, *bucket);
      for (const stk::mesh::Entity & entity : *bucket) {
        T * fieldData = static_cast<T*>(stk::mesh::field_data<stk::mesh::FieldBase, StkDebugger<T>>(stkField, entity));
        for (unsigned d=0; d<numComponents; ++d) {
          fieldData[d] = value;
        }
      }
    }
  }

  template <typename T>
  void write_vector_field_on_host_using_entity(stk::mesh::FieldBase & stkField, T value)
  {
    write_vector_field_on_host_using_entity(stkField, stkField, value);
  }

  template <typename T>
  void write_scalar_field_on_host_using_entity(stk::mesh::Field<T> & stkField, const stk::mesh::Selector& selector, T value)
  {
    write_vector_field_on_host_using_entity(stkField, selector, value);
  }

  template <typename T>
  void write_scalar_field_on_host_using_bucket(stk::mesh::Field<T> & stkField, T value)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().get_buckets(stkField.entity_rank(), stkField);
    for (stk::mesh::Bucket * bucket : buckets) {
      T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, *bucket);
      for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
        fieldData[iEntity] = value;
      }
    }
  }

  template <typename T>
  void write_scalar_field_on_host_using_bucket_id(stk::mesh::Field<T> & stkField, T value)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().get_buckets(stkField.entity_rank(), stkField);
    for (stk::mesh::Bucket * bucket : buckets) {
      const unsigned bucketId = bucket->bucket_id();
      T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, bucketId);
      for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
        fieldData[iEntity] = value;
      }
    }
  }

  template <typename T>
  stk::mesh::EntityVector write_scalar_field_on_host_using_bucket_id_and_ordinal(stk::mesh::Field<T> & stkField, T value)
  {
    stk::mesh::EntityVector entities;
    const stk::mesh::BucketVector& buckets = get_bulk().get_buckets(stkField.entity_rank(), stkField);
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        entities.push_back(entity);
        const stk::mesh::MeshIndex & meshIndex = get_bulk().mesh_index(entity);
        const unsigned bucketId = meshIndex.bucket->bucket_id();
        const stk::mesh::Bucket::size_type bucketOrd = meshIndex.bucket_ordinal;
        T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, bucketId, bucketOrd);
        fieldData[0] = value;
      }
    }
    return entities;
  }

  template <typename T>
  void write_scalar_field_on_host_using_bucket_id_and_ordinal_and_size(stk::mesh::Field<T> & stkField, T value)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().get_buckets(stkField.entity_rank(), stkField);
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        const stk::mesh::MeshIndex & meshIndex = get_bulk().mesh_index(entity);
        const unsigned bucketId = meshIndex.bucket->bucket_id();
        const stk::mesh::Bucket::size_type bucketOrd = meshIndex.bucket_ordinal;
        const unsigned numBytesPerEntity = stk::mesh::field_bytes_per_entity(stkField, *bucket);
        T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, bucketId, bucketOrd, numBytesPerEntity);
        fieldData[0] = value;
      }
    }
  }

  template <typename T>
  void write_scalar_field_on_device(stk::mesh::Field<T> & stkField, T value)
  {
    const int component = 0;
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    stk::mesh::NgpField<T, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<T, NgpDebugger>(stkField);
    stk::mesh::Selector fieldSelector(stkField);

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, fieldSelector,
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
                                     ngpField(entity, component) = value;
                                   });
  }

  template <typename T>
  void device_field_set_all(stk::mesh::Field<T> & stkField, T value)
  {
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    stk::mesh::NgpField<T, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<T, NgpDebugger>(stkField);

    ngpField.set_all(ngpMesh, value);
  }

  template <typename T>
  void read_scalar_field_on_host_using_entity(stk::mesh::Field<T> & stkField)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().get_buckets(stkField.entity_rank(), stkField);
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        const T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, entity);
        access_for_memory_checking_tool(fieldData);
      }
    }
  }

  template <typename T>
  void read_scalar_field_on_host_using_bucket(stk::mesh::Field<T> & stkField)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().get_buckets(stkField.entity_rank(), stkField);
    for (stk::mesh::Bucket * bucket : buckets) {
      const T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, *bucket);
      access_for_memory_checking_tool(fieldData, bucket->size());
    }
  }

  template <typename T>
  void read_field_on_device(stk::mesh::FieldBase & stkField, const stk::mesh::Selector& selector)
  {
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    stk::mesh::NgpField<T, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<T, NgpDebugger>(stkField);
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), selector);
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned i = 0; i < bucketIds.size(); ++i) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(i));
                             for (unsigned j = 0; j < bucket.size(); ++j) {
                               stk::mesh::FastMeshIndex index = ngpMesh.fast_mesh_index(bucket[j]);
                               unsigned numComponents = ngpField.get_num_components_per_entity(index);
                               for (unsigned component=0; component<numComponents; ++component) {
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
  void read_scalar_field_on_device(stk::mesh::Field<T> & stkField, const stk::mesh::Selector& selector)
  {
    read_field_on_device<T>(stkField, selector);
  }

  template <typename T>
  void read_scalar_field_on_device(stk::mesh::Field<T> & stkField)
  {
    read_scalar_field_on_device(stkField, stkField);
  }

  template <typename T>
  void read_vector_field_on_device(stk::mesh::FieldBase & stkField, const stk::mesh::Selector& selector)
  {
    read_field_on_device<T>(stkField, selector);
  }

  template <typename T>
  void read_vector_field_on_device(stk::mesh::FieldBase & stkField)
  {
    read_vector_field_on_device<T>(stkField, stkField);
  }

  template <typename T, typename NGPFIELD>
  void read_old_scalar_field_on_device(stk::mesh::Field<T> & stkField, NGPFIELD & ngpField,
                                       stk::mesh::EntityId maxIdToRead = std::numeric_limits<stk::mesh::EntityId>::max())
  {
    const int component = 0;
    const stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), stkField);
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned i = 0; i < bucketIds.size(); ++i) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(i));
                             for (unsigned j = 0; j < bucket.size(); ++j) {
                               stk::mesh::FastMeshIndex index = ngpMesh.fast_mesh_index(bucket[j]);
                               const stk::mesh::Entity elem = ngpMesh.get_entity(stk::topology::ELEM_RANK, index);
                               const stk::mesh::EntityId elemId = ngpMesh.identifier(elem);
                               if (elemId <= maxIdToRead) {
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
  void read_field_on_device_using_entity_field_data(stk::mesh::Field<T> & stkField)
  {
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    stk::mesh::NgpField<T, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<T, NgpDebugger>(stkField);
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), stkField);
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned bucketId = 0; bucketId < bucketIds.size(); ++bucketId) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(bucketId));
                             for (unsigned offset = 0; offset < bucket.size(); ++offset) {
                               stk::mesh::FastMeshIndex index = ngpMesh.fast_mesh_index(bucket[offset]);
#if defined(DEVICE_USE_LOCATION_BUILTINS)
                               stk::mesh::EntityFieldData<double> vals = ngpField(index);
#else
                               stk::mesh::EntityFieldData<double> vals = ngpField(index, __FILE__, __LINE__);
#endif
                               access_for_memory_checking_tool(vals, vals.size());
                             }
                           }
                         });
  }

  template <typename T>
  void read_field_on_device_using_mesh_index(stk::mesh::Field<T> & stkField)
  {
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    stk::mesh::NgpField<T, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<T, NgpDebugger>(stkField);
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), stkField);
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned i = 0; i < bucketIds.size(); ++i) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(i));
                             for (unsigned j = 0; j < bucket.size(); ++j) {
                               stk::mesh::NgpMesh::MeshIndex meshIndex{&bucket, static_cast<unsigned>(j)};
                               stk::mesh::FastMeshIndex fastMeshIndex{bucket.bucket_id(), static_cast<unsigned>(j)};
                               const unsigned numComponents = ngpField.get_num_components_per_entity(fastMeshIndex);
                               for (unsigned component = 0; component < numComponents; ++component) {
#if defined(DEVICE_USE_LOCATION_BUILTINS)
                                 access_for_memory_checking_tool(&ngpField(meshIndex, component));
#else
                                 access_for_memory_checking_tool(&ngpField(meshIndex, component, __FILE__, __LINE__));
#endif
                               }
                             }
                           }
                         });
  }

};

#endif // NGPDEBUGFIELDSYNCFIXTURES_HPP
