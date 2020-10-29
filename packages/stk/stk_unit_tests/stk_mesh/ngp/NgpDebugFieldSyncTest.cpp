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
#include <stk_util/stk_config.h>
#include <string>
#include <sstream>
#include <vector>
#include <numeric>
#include <tuple>

namespace {

template <typename T>
STK_INLINE_FUNCTION
bool dummy_check(T value)
{
  return (value == 0);
}

class NgpDebugFieldSync : public stk::unit_test_util::MeshFixture
{
public:
  template <typename T>
  void initialize_ngp_field(stk::mesh::Field<T> & stkField)
  {
    stk::mesh::get_updated_ngp_field<T>(stkField);
  }

  void set_initial_part_membership(const std::vector<std::pair<unsigned, std::string>> & numElemsInEachPart)
  {
    get_bulk().modification_begin();
    stk::mesh::EntityId elemId = 0;
    for (const auto & elemCountAndPart : numElemsInEachPart) {
      stk::mesh::PartVector addParts {get_meta().get_part(elemCountAndPart.second)};
      stk::mesh::EntityVector elemsToChange;
      for (unsigned elemNum = 0; elemNum < elemCountAndPart.first; ++elemNum) {
        stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, ++elemId);
        ThrowRequireMsg(get_bulk().is_valid(element), "Invalid element in fixture!");
        elemsToChange.push_back(element);
      }
      get_bulk().change_entity_parts(elemsToChange, addParts);
    }
    get_bulk().modification_end();
  }

  template <typename T>
  void device_field_set_all(stk::mesh::Field<T> & stkField, T value)
  {
    stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
    stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);

    ngpField.set_all(ngpMesh, value);
  }

  template <typename T>
  void write_scalar_field_on_host_using_entity(stk::mesh::Field<T> & stkField, T value)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        T * fieldData = stk::mesh::field_data(stkField, entity);
        fieldData[0] = value;
      }
    }
  }

  template <typename T>
  void write_vector_field_on_host_using_entity(stk::mesh::Field<T> & stkField, T value)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        T * fieldData = stk::mesh::field_data(stkField, entity);
        fieldData[1] = value;  // Write to second component only
      }
    }
  }

  template <typename T>
  void write_scalar_field_on_host_using_bucket(stk::mesh::Field<T> & stkField, T value)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      T * fieldData = stk::mesh::field_data(stkField, *bucket);
      for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
        fieldData[iEntity] = value;
      }
    }
  }

  template <typename T>
  void write_vector_field_on_host_using_bucket(stk::mesh::Field<T> & stkField, T value)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      T * fieldData = stk::mesh::field_data(stkField, *bucket);
      for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
        const size_t yComponent = iEntity*3 + 1;
        fieldData[yComponent] = value;
      }
    }
  }

  template <typename T>
  void read_scalar_field_on_host_using_entity(stk::mesh::Field<T> & stkField)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        const T * fieldData = stk::mesh::field_data(stkField, entity);
        dummy_check(fieldData[0]);
      }
    }
  }

  template <typename T>
  void read_scalar_field_on_host_using_bucket(stk::mesh::Field<T> & stkField)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      const T * fieldData = stk::mesh::field_data(stkField, *bucket);
      for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
        dummy_check(fieldData[iEntity]);
      }
    }
  }

  template <typename T>
  void read_vector_field_on_host_using_entity(stk::mesh::Field<T> & stkField)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        const T * fieldData = stk::mesh::field_data(stkField, entity);
        for (unsigned component = 0; component < 3; ++component) {
          dummy_check(fieldData[component]);
        }
      }
    }
  }

  template <typename T>
  void read_vector_field_on_host_using_bucket(stk::mesh::Field<T> & stkField)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      const T * fieldData = stk::mesh::field_data(stkField, *bucket);
      for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
        for (unsigned component = 0; component < 3; ++component) {
          dummy_check(fieldData[iEntity*3 + component]);
        }
      }
    }
  }

  template <typename T>
  void write_scalar_field_on_device(stk::mesh::Field<T> & stkField, T value)
  {
    const int component = 0;
    stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, meta.locally_owned_part(),
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
                                     ngpField(entity, component) = value;
                                   });
  }

  template <typename T>
  void write_vector_field_on_device(stk::mesh::Field<T> & stkField, T value)
  {
    const int component = 1;  // Just write to the second component
    stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, meta.locally_owned_part(),
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
                                     ngpField(entity, component) = value;
                                   });
  }

  template <typename T>
  void write_vector_field_on_device_using_mesh_index(stk::mesh::Field<T> & stkField, T value)
  {
    const int component = 1;  // Just write to the second component
    stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), meta.locally_owned_part());
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned i = 0; i < bucketIds.size(); ++i) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(i));
                             for (unsigned j = 0; j < bucket.size(); ++j) {
                               stk::mesh::NgpMesh::MeshIndex index{&bucket, static_cast<unsigned>(j)};
                               ngpField(index, component) = value;
                             }
                           }
                         });
  }

  template <typename T>
  void write_vector_field_on_device_using_entity_field_data(stk::mesh::Field<T> & stkField, T value)
  {
    const int component = 1;  // Just write to the second component
    stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, meta.locally_owned_part(),
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
#if defined(STK_DEBUG_FIELD_SYNC) && !defined(DEVICE_USE_LOCATION_BUILTINS)
                                     stk::mesh::EntityFieldData<double> vals = ngpField(entity, __FILE__, __LINE__);
#else
                                     stk::mesh::EntityFieldData<double> vals = ngpField(entity);
#endif
                                     vals[component] = value;
                                   });
  }

  template <typename T>
  void read_scalar_field_on_device(stk::mesh::Field<T> & stkField)
  {
    const int component = 0;
    stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), meta.locally_owned_part());
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned i = 0; i < bucketIds.size(); ++i) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(i));
                             for (unsigned j = 0; j < bucket.size(); ++j) {
                               stk::mesh::FastMeshIndex index = ngpMesh.fast_mesh_index(bucket[j]);
#if defined(STK_DEBUG_FIELD_SYNC) && !defined(DEVICE_USE_LOCATION_BUILTINS)
                               dummy_check(ngpField(index, component, __FILE__, __LINE__));
#else
                               dummy_check(ngpField(index, component));
#endif
                             }
                           }
                         });
  }

  template <typename T, typename NGPFIELD>
  void read_old_scalar_field_on_device(stk::mesh::Field<T> & stkField, NGPFIELD & ngpField,
                                       stk::mesh::EntityId maxIdToRead = std::numeric_limits<stk::mesh::EntityId>::max())
  {
    const int component = 0;
    const stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), meta.locally_owned_part());
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned i = 0; i < bucketIds.size(); ++i) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(i));
                             for (unsigned j = 0; j < bucket.size(); ++j) {
                               stk::mesh::FastMeshIndex index = ngpMesh.fast_mesh_index(bucket[j]);
                               const stk::mesh::Entity elem = ngpMesh.get_entity(stk::topology::ELEM_RANK, index);
                               const stk::mesh::EntityId elemId = ngpMesh.identifier(elem);
                               if (elemId <= maxIdToRead) {
#if defined(STK_DEBUG_FIELD_SYNC) && !defined(DEVICE_USE_LOCATION_BUILTINS)
                                 dummy_check(ngpField(index, component, __FILE__, __LINE__));
#else
                                 dummy_check(ngpField(index, component));
#endif
                               }
                             }
                           }
                         });
  }

  template <typename T>
  void read_vector_field_on_device(stk::mesh::Field<T> & stkField)
  {
    stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), meta.locally_owned_part());
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned i = 0; i < bucketIds.size(); ++i) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(i));
                             for (unsigned j = 0; j < bucket.size(); ++j) {
                               stk::mesh::FastMeshIndex index = ngpMesh.fast_mesh_index(bucket[j]);
                               for (int component = 0; component < 3; ++component) {
#if defined(STK_DEBUG_FIELD_SYNC) && !defined(DEVICE_USE_LOCATION_BUILTINS)
                                 dummy_check(ngpField(index, component, __FILE__, __LINE__));
#else
                                 dummy_check(ngpField(index, component));
#endif
                               }
                             }
                           }
                         });
  }

  template <typename T>
  void read_vector_field_on_device_using_mesh_index(stk::mesh::Field<T> & stkField)
  {
    stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), meta.locally_owned_part());
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned i = 0; i < bucketIds.size(); ++i) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(i));
                             for (unsigned j = 0; j < bucket.size(); ++j) {
                               stk::mesh::NgpMesh::MeshIndex index{&bucket, static_cast<unsigned>(j)};
                               for (int component = 0; component < 3; ++component) {
#if defined(STK_DEBUG_FIELD_SYNC) && !defined(DEVICE_USE_LOCATION_BUILTINS)
                                 dummy_check(ngpField(index, component, __FILE__, __LINE__));
#else
                                 dummy_check(ngpField(index, component));
#endif
                               }
                             }
                           }
                         });
  }

  template <typename T>
  void read_vector_field_on_device_using_entity_field_data(stk::mesh::Field<T> & stkField)
  {
    stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), meta.locally_owned_part());
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned bucketId = 0; bucketId < bucketIds.size(); ++bucketId) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(bucketId));
                             for (unsigned offset = 0; offset < bucket.size(); ++offset) {
                               stk::mesh::FastMeshIndex index = ngpMesh.fast_mesh_index(bucket[offset]);
#if defined(STK_DEBUG_FIELD_SYNC) && !defined(DEVICE_USE_LOCATION_BUILTINS)
                               stk::mesh::EntityFieldData<double> vals = ngpField(index, __FILE__, __LINE__);
#else
                               stk::mesh::EntityFieldData<double> vals = ngpField(index);
#endif
                               for (unsigned i = 0; i < vals.size(); ++i) {
                                 dummy_check(vals[i]);
                               }
                             }
                           }
                         });
  }

  void modify_element_part_membership(const std::vector<std::tuple<stk::mesh::EntityId, std::string, std::string>> & elemAddRemoveParts)
  {
    get_bulk().modification_begin();
    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, std::get<0>(elemAddRemovePart))};
      stk::mesh::PartVector addParts {get_meta().get_part(std::get<1>(elemAddRemovePart))};
      stk::mesh::PartVector removeParts {get_meta().get_part(std::get<2>(elemAddRemovePart))};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }
    get_bulk().modification_end();
  }

  template <typename T>
  void create_element(const std::vector<std::pair<stk::mesh::EntityId, std::string>> & elemParts,
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
    get_bulk().modification_end();

    fill_initial_field<T>(stkField);
  }

  void delete_element(const std::vector<stk::mesh::EntityId> & elemIds)
  {
    get_bulk().modification_begin();
    for (const stk::mesh::EntityId & elemId : elemIds) {
      get_bulk().destroy_entity(get_bulk().get_entity(stk::topology::ELEM_RANK, elemId));
    }
    get_bulk().modification_end();
  }

  template <typename T>
  void modify_element_part_membership_with_scalar_field_write_using_entity(
           const std::vector<std::tuple<stk::mesh::EntityId, std::string, std::string>> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_scalar_field_on_host_using_entity(stkField, value);

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, std::get<0>(elemAddRemovePart))};
      stk::mesh::PartVector addParts {get_meta().get_part(std::get<1>(elemAddRemovePart))};
      stk::mesh::PartVector removeParts {get_meta().get_part(std::get<2>(elemAddRemovePart))};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    get_bulk().modification_end();
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
      T * fieldData = stk::mesh::field_data(stkField, newElem);
      fieldData[0] = value;
    }

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
  void modify_element_part_membership_with_vector_field_write_using_entity(
           const std::vector<std::tuple<stk::mesh::EntityId, std::string, std::string>> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_vector_field_on_host_using_entity(stkField, value);

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, std::get<0>(elemAddRemovePart))};
      stk::mesh::PartVector addParts {get_meta().get_part(std::get<1>(elemAddRemovePart))};
      stk::mesh::PartVector removeParts {get_meta().get_part(std::get<2>(elemAddRemovePart))};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
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
      T * fieldData = stk::mesh::field_data(stkField, newElem);
      fieldData[1] = value;
    }

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
  void modify_element_part_membership_with_scalar_field_write_using_bucket(
           const std::vector<std::tuple<stk::mesh::EntityId, std::string, std::string>> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_scalar_field_on_host_using_bucket(stkField, value);

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, std::get<0>(elemAddRemovePart))};
      stk::mesh::PartVector addParts {get_meta().get_part(std::get<1>(elemAddRemovePart))};
      stk::mesh::PartVector removeParts {get_meta().get_part(std::get<2>(elemAddRemovePart))};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
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
      T * fieldData = stk::mesh::field_data(stkField, bucket);
      for(size_t iEntity = 0; iEntity < bucket.size(); ++iEntity) {
        fieldData[iEntity] = value;
      }
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
  void modify_element_part_membership_with_vector_field_write_using_bucket(
           const std::vector<std::tuple<stk::mesh::EntityId, std::string, std::string>> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField, T value)
  {
    get_bulk().modification_begin();

    write_vector_field_on_host_using_bucket(stkField, value);

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, std::get<0>(elemAddRemovePart))};
      stk::mesh::PartVector addParts {get_meta().get_part(std::get<1>(elemAddRemovePart))};
      stk::mesh::PartVector removeParts {get_meta().get_part(std::get<2>(elemAddRemovePart))};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
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
      T * fieldData = stk::mesh::field_data(stkField, bucket);
      for(size_t iEntity = 0; iEntity < bucket.size(); ++iEntity) {
        const size_t yComponent = iEntity*3 + 1;
        fieldData[yComponent] = value;
      }
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
  void modify_element_part_membership_with_scalar_field_read_using_entity(
           const std::vector<std::tuple<stk::mesh::EntityId, std::string, std::string>> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField)
  {
    get_bulk().modification_begin();

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, std::get<0>(elemAddRemovePart))};
      stk::mesh::PartVector addParts {get_meta().get_part(std::get<1>(elemAddRemovePart))};
      stk::mesh::PartVector removeParts {get_meta().get_part(std::get<2>(elemAddRemovePart))};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    read_scalar_field_on_host_using_entity(stkField);

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
  void modify_element_part_membership_with_vector_field_read_using_entity(
           const std::vector<std::tuple<stk::mesh::EntityId, std::string, std::string>> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField)
  {
    get_bulk().modification_begin();

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, std::get<0>(elemAddRemovePart))};
      stk::mesh::PartVector addParts {get_meta().get_part(std::get<1>(elemAddRemovePart))};
      stk::mesh::PartVector removeParts {get_meta().get_part(std::get<2>(elemAddRemovePart))};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    read_vector_field_on_host_using_entity(stkField);

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

    read_vector_field_on_host_using_entity(stkField);

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

    read_vector_field_on_host_using_entity(stkField);

    get_bulk().modification_end();
  }

  template <typename T>
  void modify_element_part_membership_with_scalar_field_read_using_bucket(
           const std::vector<std::tuple<stk::mesh::EntityId, std::string, std::string>> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField)
  {
    get_bulk().modification_begin();

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, std::get<0>(elemAddRemovePart))};
      stk::mesh::PartVector addParts {get_meta().get_part(std::get<1>(elemAddRemovePart))};
      stk::mesh::PartVector removeParts {get_meta().get_part(std::get<2>(elemAddRemovePart))};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    read_scalar_field_on_host_using_bucket(stkField);

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
  void modify_element_part_membership_with_vector_field_read_using_bucket(
           const std::vector<std::tuple<stk::mesh::EntityId, std::string, std::string>> & elemAddRemoveParts,
           stk::mesh::Field<T> & stkField)
  {
    get_bulk().modification_begin();

    for (const auto & elemAddRemovePart : elemAddRemoveParts) {
      stk::mesh::EntityVector elemsToChange {get_bulk().get_entity(stk::topology::ELEM_RANK, std::get<0>(elemAddRemovePart))};
      stk::mesh::PartVector addParts {get_meta().get_part(std::get<1>(elemAddRemovePart))};
      stk::mesh::PartVector removeParts {get_meta().get_part(std::get<2>(elemAddRemovePart))};
      get_bulk().change_entity_parts(elemsToChange, addParts, removeParts);
    }

    read_vector_field_on_host_using_bucket(stkField);

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
  stk::mesh::Field<T> & create_scalar_field(stk::topology::rank_t rank, const std::string & name)
  {
    unsigned numStates = 1;
    const T init = 1;
    stk::mesh::Field<T> & field = get_meta().declare_field<stk::mesh::Field<T>>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
    return field;
  }

  template <typename T>
  stk::mesh::Field<T> & create_scalar_multistate_field(stk::topology::rank_t rank, const std::string & name)
  {
    unsigned numStates = 2;
    const T init = 1;
    stk::mesh::Field<T> & field = get_meta().declare_field<stk::mesh::Field<T>>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
    return field;
  }

  template <typename T>
  stk::mesh::Field<T> & create_vector_field(stk::topology::rank_t rank, const std::string & name)
  {
    unsigned numStates = 1;
    unsigned numScalarsPerEntity = 3;
    const T init[] = {1, 2, 3};
    stk::mesh::Field<T> & field = get_meta().declare_field<stk::mesh::Field<T>>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), numScalarsPerEntity, init);
    return field;
  }

  template <typename T>
  void fill_initial_field(stk::mesh::Field<T> & stkField)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        const stk::mesh::EntityId id = get_bulk().identifier(entity);
#ifdef STK_DEBUG_FIELD_SYNC
        T * fieldData = stk::mesh::ngp_debug_field_data(stkField, entity);
#else
        T * fieldData = stk::mesh::field_data(stkField, entity);
#endif
        const unsigned numComponents = stk::mesh::field_scalars_per_entity(stkField, *bucket);
        for (unsigned component = 0; component < numComponents; ++component) {
          fieldData[component] = 10*id + component;
        }
      }
    }
  }

  void
  create_parts(const std::vector<std::pair<unsigned, std::string>> & numElemsInEachPart)
  {
    for (const auto & elemCountAndPart : numElemsInEachPart) {
      get_meta().declare_part_with_topology(elemCountAndPart.second, stk::topology::HEX_8);
    }
  }

  template <typename T>
  stk::mesh::Field<T> & build_mesh_with_scalar_field(const std::string & fieldName,
                                                     const std::vector<std::pair<unsigned, std::string>> & numElemsInEachPart,
                                                     unsigned bucketCapacity = stk::mesh::impl::BucketRepository::default_bucket_capacity)
  {
    stk::mesh::Field<T> & stkField = create_scalar_field<T>(stk::topology::ELEM_RANK, fieldName);
    create_parts(numElemsInEachPart);

    unsigned numElems = 0;
    for (const auto & elemCountAndPart : numElemsInEachPart) {
      numElems += elemCountAndPart.first;
    }
    setup_mesh("generated:1x1x" + std::to_string(numElems), stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);

    set_initial_part_membership(numElemsInEachPart);
    fill_initial_field<T>(stkField);
    initialize_ngp_field(stkField);
    return stkField;
  }

  template <typename T>
  stk::mesh::Field<T> & build_mesh_with_scalar_multistate_field(const std::string & fieldName,
                                                                const std::vector<std::pair<unsigned, std::string>> & numElemsInEachPart,
                                                                unsigned bucketCapacity = stk::mesh::impl::BucketRepository::default_bucket_capacity)
  {
    stk::mesh::Field<T> & stkField = create_scalar_multistate_field<T>(stk::topology::ELEM_RANK, fieldName);
    create_parts(numElemsInEachPart);

    unsigned numElems = 0;
    for (const auto & elemCountAndPart : numElemsInEachPart) {
      numElems += elemCountAndPart.first;
    }
    setup_mesh("generated:1x1x" + std::to_string(numElems), stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);

    set_initial_part_membership(numElemsInEachPart);
    fill_initial_field<T>(stkField);
    initialize_ngp_field(stkField);
    return stkField;
  }

  template <typename T>
  stk::mesh::Field<T> & build_mesh_with_vector_field(const std::string & fieldName,
                                                     const std::vector<std::pair<unsigned, std::string>> & numElemsInEachPart,
                                                     unsigned bucketCapacity = stk::mesh::impl::BucketRepository::default_bucket_capacity)
  {
    stk::mesh::Field<T> & stkField = create_vector_field<T>(stk::topology::ELEM_RANK, fieldName);
    create_parts(numElemsInEachPart);

    unsigned numElems = 0;
    for (const auto & elemCountAndPart : numElemsInEachPart) {
      numElems += elemCountAndPart.first;
    }
    setup_mesh("generated:1x1x" + std::to_string(numElems), stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);

    set_initial_part_membership(numElemsInEachPart);
    fill_initial_field<T>(stkField);
    initialize_ngp_field(stkField);
    return stkField;
  }

};


template <typename Out>
void split_lines(const std::string &s, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item)) {
        *result++ = item;
    }
}

std::vector<std::string> split_lines(const std::string &s) {
    std::vector<std::string> elems;
    split_lines(s, std::back_inserter(elems));
    return elems;
}

void extract_warning(std::string & stdoutString, int numExpectedOccurrences, const std::string & warningString)
{
  std::vector<std::string> warningLines = split_lines(stdoutString);
  std::string newStdoutString;
  int numFound = 0;

  for (const std::string & line : warningLines) {
    const size_t loc = line.find(warningString);
    if (loc != std::string::npos) {
      ++numFound;
    }
    else {
      newStdoutString += line + "\n";
    }
  }

#if defined(STK_DEBUG_FIELD_SYNC) && defined(STK_USE_DEVICE_MESH)
  if (numFound != numExpectedOccurrences) {
  std::cout << "Warning string found " << numFound << " times when expecting " << numExpectedOccurrences << " occurrences: \""
            << warningString << "\"" << std::endl;
    ADD_FAILURE();
  }
#endif

  stdoutString = newStdoutString;
}

void check_no_warnings(const std::string & stdoutString)
{
  std::vector<std::string> warningLines = split_lines(stdoutString);

  for (const std::string & line : warningLines) {
    if (!line.empty()) {
      std::cout << "Found unexpected warning: \"" << line << "\"" << std::endl;
      ADD_FAILURE();
    }
  }
}

void check_contains_file_name(const std::string & stdoutString, const std::string & fileName)
{
#if defined(STK_DEBUG_FIELD_SYNC) && defined(STK_USE_DEVICE_MESH) && defined(HOST_USE_LOCATION_BUILTINS)
  const size_t fileNameLoc = stdoutString.find(fileName);
  EXPECT_NE(fileNameLoc, std::string::npos);
#endif
}

void check_contains_a_line_number(const std::string & stdoutString)
{
#if defined(STK_DEBUG_FIELD_SYNC) && defined(STK_USE_DEVICE_MESH) && defined(HOST_USE_LOCATION_BUILTINS)
  const size_t colonLoc = stdoutString.find(":");
  ASSERT_NE(colonLoc, std::string::npos);
  const size_t lineNumberLoc = colonLoc + 1;
  int lineNumber = std::stoi(stdoutString.data()+lineNumberLoc);
  EXPECT_GT(lineNumber, 0);
#endif
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_contains_file_name(stdoutString, "NgpDebugFieldSyncTest.cpp");
  check_contains_a_line_number(stdoutString);
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MissingSyncToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MissingModifyOnHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3.14);

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceEntityFieldDataAccess_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device_using_entity_field_data(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceEntityFieldDataAccess_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3.14);

  read_vector_field_on_device_using_entity_field_data(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceMeshIndexAccess_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device_using_mesh_index(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceMeshIndexAccess_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3.14);

  read_vector_field_on_device_using_mesh_index(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarIntAccessUsingEntity_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_scalar_field<int>("intScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarIntAccessUsingEntity_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_scalar_field<int>("intScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field intScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingEntity_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field<int>("intVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingEntity_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field<int>("intVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3);

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field intVectorField[1]=11.000000");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket(stkField, 3.14);

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field<int>("intScalarField", {{1, "Part1"},
                                                                                          {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket(stkField, 3);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field<int>("intVectorField", {{1, "Part1"},
                                                                                          {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket(stkField, 3);

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field intVectorField[1]=11.000000");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_contains_file_name(stdoutString, "NgpDebugFieldSyncTest.cpp");
  check_contains_a_line_number(stdoutString);
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MissingSyncToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MissingModifyOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_vector_field_on_host_using_entity(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device(stkField, 3.14);

  read_vector_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=21");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceEntityFieldDataAccess_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device_using_entity_field_data(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_vector_field_on_host_using_entity(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceEntityFieldDataAccess_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device_using_entity_field_data(stkField, 3.14);

  read_vector_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=21");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceMeshIndexAccess_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device_using_mesh_index(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_vector_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceMeshIndexAccess_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device_using_mesh_index(stkField, 3.14);

  read_vector_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=21");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingEntity_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field<int>("intVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device(stkField, 3);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_vector_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingEntity_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field<int>("intVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device(stkField, 3);

  read_vector_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[1]=21");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);

  read_scalar_field_on_host_using_bucket(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_vector_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device(stkField, 3.14);

  read_vector_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=21");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field<int>("intVectorField", {{1, "Part1"},
                                                                                          {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device(stkField, 3);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_vector_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field<int>("intVectorField", {{1, "Part1"},
                                                                                          {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device(stkField, 3);

  read_vector_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[1]=21");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarWriteOnHost_ProperlyMarkAsModified_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarWriteOnHost_MissingMarkAsModified_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarWriteOnHost_ProperlyMarkAsModified_ClearHostSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_host_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarWriteOnHost_MissingMarkAsModified_ClearHostSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_host_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarWriteOnHost_ProperlySyncToDevice_ClearDeviceSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();
  stkField.clear_device_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarWriteOnHost_ProperlyMarkAsModified_ClearDeviceSyncState_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_device_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarWriteOnHost_MissingMarkAsModified_ClearDeviceSyncState_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_device_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarWriteOnDevice_ProperlyMarkAsModified_ClearSyncState_AccessOnHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarWriteOnDevice_MissingMarkAsModified_ClearSyncState_AccessOnHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_sync_state();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarWriteOnDevice_ProperlySyncToHost_ClearHostSyncState_AccessOnHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();
  stkField.clear_host_sync_state();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarWriteOnDevice_ProperlyMarkAsModified_ClearHostSyncState_AccessOnHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_host_sync_state();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarWriteOnDevice_MissingMarkAsModified_ClearHostSyncState_AccessOnHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_host_sync_state();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarWriteOnDevice_ProperlyMarkAsModified_ClearDeviceSyncState_AccessOnHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_device_sync_state();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarWriteOnDevice_MissingMarkAsModified_ClearDeviceSyncState_AccessOnHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_device_sync_state();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarDeviceSetAll_AccessOnHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  device_field_set_all(stkField, 2.18);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarDeviceSetAll_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  device_field_set_all(stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarDeviceSetAll_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);

  device_field_set_all(stkField, 2.18);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarDeviceSetAll_MissingAllModifySyncCallsToHost_AccessOnHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);

  device_field_set_all(stkField, 2.18);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarDeviceSetAll_MissingAllModifySyncCallsToHost_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);

  device_field_set_all(stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleTimestep_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_entity(stkField, 3.14+timeStep);
    stkField.modify_on_host();
    stkField.sync_to_device();
    read_scalar_field_on_device(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleTimestep_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_entity(stkField, 3.14+timeStep);
    read_scalar_field_on_device(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleTimestep_ProperlyMarkAsModified_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_entity(stkField, 3.14+timeStep);
    stkField.modify_on_host();
    stkField.clear_sync_state();
    read_scalar_field_on_device(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleTimestep_MissingMarkAsModified_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_entity(stkField, 3.14+timeStep);
    stkField.clear_sync_state();
    read_scalar_field_on_device(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleAccesses_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleStaleAccesses_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_scalar_field_on_device(stkField);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleWrites_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleWrites_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  write_scalar_field_on_host_using_entity(stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_MultipleTimestep_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_vector_field_on_host_using_entity(stkField, 3.14+timeStep);
    stkField.modify_on_host();
    stkField.sync_to_device();
    read_vector_field_on_device(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_MultipleTimestep_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_vector_field_on_host_using_entity(stkField, 3.14+timeStep);
    read_vector_field_on_device(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=21.000000");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleTimestep_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_bucket(stkField, 3.14+timeStep);
    stkField.modify_on_host();
    stkField.sync_to_device();
    read_scalar_field_on_device(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleTimestep_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_bucket(stkField, 3.14+timeStep);
    read_scalar_field_on_device(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleAccesses_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_bucket(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleStaleAccesses_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_bucket(stkField, 3.14);

  read_scalar_field_on_device(stkField);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleWrites_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_bucket(stkField, 3.14);
  write_scalar_field_on_host_using_bucket(stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleWrites_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_bucket(stkField, 3.14);
  write_scalar_field_on_host_using_bucket(stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20.000000");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleTimestep_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(stkField, 3.14+timeStep);
    stkField.modify_on_device();
    stkField.sync_to_host();
    read_scalar_field_on_host_using_entity(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleTimestep_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(stkField, 3.14+timeStep);
    read_scalar_field_on_host_using_entity(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleTimestep_ProperlyMarkAsModified_ClearSyncState_AccessOnHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(stkField, 3.14+timeStep);
    stkField.modify_on_device();
    stkField.clear_sync_state();
    read_scalar_field_on_host_using_entity(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleTimestep_MissingMarkAsModified_ClearSyncState_AccessOnHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(stkField, 3.14+timeStep);
    stkField.clear_sync_state();
    read_scalar_field_on_host_using_entity(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleAccesses_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleStaleAccesses_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  read_scalar_field_on_host_using_entity(stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleWrites_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleWrites_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  write_scalar_field_on_device(stkField, 2.18);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleTimestep_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(stkField, 3.14+timeStep);
    stkField.modify_on_device();
    stkField.sync_to_host();
    read_scalar_field_on_host_using_bucket(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleTimestep_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(stkField, 3.14+timeStep);
    read_scalar_field_on_host_using_bucket(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleAccesses_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_bucket(stkField);
  read_scalar_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleStaleAccesses_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  read_scalar_field_on_host_using_bucket(stkField);
  read_scalar_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleWrites_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleWrites_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  write_scalar_field_on_device(stkField, 2.18);

  read_scalar_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  delete_element({2});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ChangeBucket_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_CreateBucket_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_DeleteBucket_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  delete_element({2});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_entity({{3, "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_DuringMeshModification_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_vector_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_DuringMeshModification_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  create_element_with_vector_field_write_using_entity({{3, "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_DuringMeshModification_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  delete_element_with_vector_field_write_using_entity({2}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_DuringMeshModification_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_scalar_field_write_using_bucket({{2, "Part2", "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_DuringMeshModification_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_bucket({{3, "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_DuringMeshModification_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_bucket({2}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_DuringMeshModification_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_vector_field_write_using_bucket({{2, "Part2", "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_DuringMeshModification_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  create_element_with_vector_field_write_using_bucket({{3, "Part1"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_DuringMeshModification_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  delete_element_with_vector_field_write_using_bucket({2}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ChangeBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.14");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_CreateBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.14");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_DeleteBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ChangeBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_CreateBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_DeleteBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_ChangeBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.14");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_CreateBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.14");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_DeleteBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_DuringMeshModification_ChangeBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_vector_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 3.14);
  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=3.14");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_DuringMeshModification_CreateBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  create_element_with_vector_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=21");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=3.14");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_DuringMeshModification_DeleteBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  delete_element_with_vector_field_write_using_entity({2}, stkField, 3.14);
  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_DuringMeshModification_ChangeBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_bucket({{2, "Part2", "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.14");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_DuringMeshModification_CreateBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  create_element_with_scalar_field_write_using_bucket({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.14");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_DuringMeshModification_DeleteBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  delete_element_with_scalar_field_write_using_bucket({2}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_DuringMeshModification_ChangeBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_vector_field_write_using_bucket({{2, "Part2", "Part1"}}, stkField, 3.14);
  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=3.14");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_DuringMeshModification_CreateBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  create_element_with_vector_field_write_using_bucket({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=21");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=3.14");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_DuringMeshModification_DeleteBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  delete_element_with_vector_field_write_using_bucket({2}, stkField, 3.14);
  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ChangeBucket_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part2", "Part1"}});

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_CreateBucket_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);

  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation
  read_old_scalar_field_on_device(stkField, ngpField, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_DeleteBucket_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  delete_element({2});

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_ChangeBucket_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 3.14);

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_CreateBucket_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.

  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation
  read_old_scalar_field_on_device(stkField, ngpField, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_DeleteBucket_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ModifyBucket_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});
  stk::mesh::NgpField<double> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part2", "Part1"}});

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_old_scalar_field_on_device(stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_CreateBucket_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation
  read_old_scalar_field_on_device(stkField, ngpFieldCopy, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_DeleteBucket_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  delete_element({2});

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_old_scalar_field_on_device(stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ModifyBucket_StaleDeviceFieldCopy_ClearSyncState_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});
  stk::mesh::NgpField<double> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part2", "Part1"}});

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_old_scalar_field_on_device(stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_CreateBucket_StaleDeviceFieldCopy_ClearSyncState_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation
  read_old_scalar_field_on_device(stkField, ngpFieldCopy, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_DeleteBucket_StaleDeviceFieldCopy_ClearSyncState_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  delete_element({2});

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_old_scalar_field_on_device(stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_ModifyBucket_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});
  stk::mesh::NgpField<double> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double>(stkField);

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

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_CreateBucket_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double>(stkField);

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

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_DeleteBucket_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double>(stkField);

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


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership({{2, "Part2", "Part1"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element({{3, "Part1"}}, stkField);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element({2});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ChangeBucket_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  modify_element_part_membership({{2, "Part2", "Part1"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_CreateBucket_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  create_element({{3, "Part1"}}, stkField);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_DeleteBucket_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  delete_element({2});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership_with_scalar_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element_with_scalar_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element_with_scalar_field_read_using_entity({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_DuringMeshModification_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership_with_vector_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_DuringMeshModification_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element_with_vector_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_DuringMeshModification_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element_with_vector_field_read_using_entity({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_DuringMeshModification_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership_with_scalar_field_read_using_bucket({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_DuringMeshModification_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element_with_scalar_field_read_using_bucket({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_DuringMeshModification_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element_with_scalar_field_read_using_bucket({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_DuringMeshModification_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership_with_vector_field_read_using_bucket({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_DuringMeshModification_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element_with_vector_field_read_using_bucket({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_DuringMeshModification_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element_with_vector_field_read_using_bucket({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership({{2, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  create_element({{3, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  delete_element({2});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ChangeBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership({{2, "Part2", "Part1"}});
  stkField.clear_sync_state();
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_CreateBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  create_element({{3, "Part1"}}, stkField);
  stkField.clear_sync_state();
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_DeleteBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  delete_element({2});
  stkField.clear_sync_state();
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership_with_scalar_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  create_element_with_scalar_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringMeshModification_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  delete_element_with_scalar_field_read_using_entity({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_DuringMeshModification_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  modify_element_part_membership_with_vector_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Lost Device values for Field doubleVectorField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_DuringMeshModification_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  create_element_with_vector_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleVectorField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_DuringMeshModification_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  delete_element_with_vector_field_read_using_entity({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleVectorField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_DuringMeshModification_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership_with_scalar_field_read_using_bucket({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_DuringMeshModification_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  create_element_with_scalar_field_read_using_bucket({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_DuringMeshModification_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  delete_element_with_scalar_field_read_using_bucket({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_DuringMeshModification_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  modify_element_part_membership_with_vector_field_read_using_bucket({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Lost Device values for Field doubleVectorField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_DuringMeshModification_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  create_element_with_vector_field_read_using_bucket({{3, "Part1"}, {4, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleVectorField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_DuringMeshModification_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_vector_field_on_device(stkField, 3.14);
  delete_element_with_vector_field_read_using_bucket({2}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleVectorField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_CreateBucket_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  delete_element({2});
  delete_element({3});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_ModifyOnHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_CreateBucket_CreateBucket_ModifyOnHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_ModifyOnHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  delete_element({2});
  delete_element({3});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{3, "Part2", "Part1"}}, stkField, 3.14);
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 2.18);

  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_CreateBucket_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  create_element_with_scalar_field_write_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField, 2.18);

  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  delete_element_with_scalar_field_write_using_entity({3}, stkField, 2.18);

  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.140000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_CreateBucket_CreateBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.140000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  delete_element({2});
  delete_element({3});

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_CreateBucket_CreateBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  delete_element({2});
  delete_element({3});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{3, "Part2", "Part1"}}, stkField, 3.14);
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=2.18");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_CreateBucket_CreateBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  create_element_with_scalar_field_write_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 3, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=2.18");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_MissingAllModifySyncCalls_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  delete_element_with_scalar_field_write_using_entity({3}, stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_CreateBucket_CreateBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation
  read_old_scalar_field_on_device(stkField, ngpField, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  delete_element({2});
  delete_element({3});

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{3, "Part2", "Part1"}}, stkField, 3.14);
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 2.18);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_CreateBucket_CreateBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  create_element_with_scalar_field_write_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField, 2.18);

  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation
  read_old_scalar_field_on_device(stkField, ngpField, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  delete_element_with_scalar_field_write_using_entity({3}, stkField, 2.18);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_CreateBucket_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element({2});
  delete_element({3});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_ModifyOnDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_CreateBucket_CreateBucket_ModifyOnDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_ModifyOnDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  delete_element({2});
  delete_element({3});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership_with_scalar_field_read_using_entity({{3, "Part2", "Part1"}}, stkField);
  modify_element_part_membership_with_scalar_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_CreateBucket_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element_with_scalar_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);
  create_element_with_scalar_field_read_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element_with_scalar_field_read_using_entity({2}, stkField);
  delete_element_with_scalar_field_read_using_entity({3}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_CreateBucket_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  delete_element({2});
  delete_element({3});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_sync_state();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_CreateBucket_CreateBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_sync_state();

  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_sync_state();

  delete_element({2});
  delete_element({3});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  modify_element_part_membership_with_scalar_field_read_using_entity({{3, "Part2", "Part1"}}, stkField);
  modify_element_part_membership_with_scalar_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 8, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_CreateBucket_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  create_element_with_scalar_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);
  create_element_with_scalar_field_read_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 10, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoConsecutiveMeshModifications_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  delete_element_with_scalar_field_read_using_entity({2}, stkField);
  delete_element_with_scalar_field_read_using_entity({3}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 5, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}



TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_ChangeBucket_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_CreateBucket_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  create_element({{4, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_DeleteBucket_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  delete_element({3});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_ChangeBucket_ChangeBucket_ProperlySyncToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_CreateBucket_CreateBucket_ProperlySyncToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  create_element({{4, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_DeleteBucket_DeleteBucket_ProperlySyncToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  delete_element({3});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_ChangeBucket_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

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

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_CreateBucket_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

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

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_DeleteBucket_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

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



TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.14");
  extract_warning(stdoutString, 3, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=2.18");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_CreateBucket_CreateBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  create_element({{4, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.140000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=2.180000");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  delete_element({3});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.140000");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_CreateBucket_CreateBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  create_element({{4, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  delete_element({3});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{3, "Part2", "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 2.18);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.14");
  extract_warning(stdoutString, 3, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=2.18");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_CreateBucket_CreateBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  create_element_with_scalar_field_write_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  create_element_with_scalar_field_write_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField, 2.18);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20.000000");
  extract_warning(stdoutString, 3, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.140000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=2.180000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  delete_element_with_scalar_field_write_using_entity({3}, stkField, 2.18);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=3.140000");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}



TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_ChangeBucket_ChangeBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_old_scalar_field_on_device(stkField, ngpField);

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 8, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_CreateBucket_CreateBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_old_scalar_field_on_device(stkField, ngpField, maxIdToRead);

  create_element({{4, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  read_old_scalar_field_on_device(stkField, ngpField, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_DeleteBucket_DeleteBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_old_scalar_field_on_device(stkField, ngpField);

  delete_element({3});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 5, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_ChangeBucket_ChangeBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{3, "Part2", "Part1"}}, stkField, 3.14);
  read_old_scalar_field_on_device(stkField, ngpField);

  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part2", "Part1"}}, stkField, 2.18);
  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 8, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_CreateBucket_CreateBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

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

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_DeleteBucket_DeleteBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();

  delete_element_with_scalar_field_write_using_entity({2}, stkField, 3.14);
  read_old_scalar_field_on_device(stkField, ngpField);

  delete_element_with_scalar_field_write_using_entity({3}, stkField, 2.18);
  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 5, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}



TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_ChangeBucket_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();
  modify_element_part_membership({{3, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();
  modify_element_part_membership({{2, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_CreateBucket_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();
  create_element({{3, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();
  create_element({{4, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_DeleteBucket_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();
  delete_element({2});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();
  delete_element({3});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_ChangeBucket_ChangeBucket_ModifyOnDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();
  modify_element_part_membership({{3, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.clear_sync_state();
  modify_element_part_membership({{2, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_CreateBucket_CreateBucket_ModifyOnDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();
  create_element({{3, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.clear_sync_state();
  create_element({{4, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_DeleteBucket_DeleteBucket_ModifyOnDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();
  delete_element({2});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.clear_sync_state();
  delete_element({3});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_ChangeBucket_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

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

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_CreateBucket_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

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

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_DeleteBucket_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

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


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership({{3, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  modify_element_part_membership({{2, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 8, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_CreateBucket_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  create_element({{3, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  create_element({{4, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 7, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  delete_element({2});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  delete_element({3});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 5, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_sync_state();
  modify_element_part_membership({{3, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.clear_sync_state();
  modify_element_part_membership({{2, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_CreateBucket_CreateBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_sync_state();
  create_element({{3, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.clear_sync_state();
  create_element({{4, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_sync_state();
  delete_element({2});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.clear_sync_state();
  delete_element({3});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership_with_scalar_field_read_using_entity({{3, "Part2", "Part1"}}, stkField);

  write_scalar_field_on_device(stkField, 2.18);
  modify_element_part_membership_with_scalar_field_read_using_entity({{2, "Part2", "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 8, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_CreateBucket_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  create_element_with_scalar_field_read_using_entity({{3, "Part1"}, {4, "Part1"}}, stkField);

  write_scalar_field_on_device(stkField, 2.18);
  create_element_with_scalar_field_read_using_entity({{5, "Part1"}, {6, "Part1"}}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 10, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_DuringTwoMeshModifications_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{3, "Part1"},
                                                                                                   {1, "Part2"}}, bucketCapacity);

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  delete_element_with_scalar_field_read_using_entity({2}, stkField);

  write_scalar_field_on_device(stkField, 2.18);
  delete_element_with_scalar_field_read_using_entity({3}, stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 5, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_FieldStateRotation_Throw)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_multistate_field<double>("doubleScalarField", {{2, "Part1"}});

#if defined(STK_DEBUG_FIELD_SYNC) && defined(STK_USE_DEVICE_MESH)
  EXPECT_THROW(get_bulk().update_field_data_states(), std::logic_error);
#endif
}

}
