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
#include <stk_util/stk_config.h>
#include <string>
#include <sstream>
#include <vector>

namespace {

class NgpDebugFieldSync : public stk::unit_test_util::MeshFixture
{
public:
  template <typename T>
  void initialize_ngp_field(stk::mesh::Field<T> & stkField)
  {
    stk::mesh::get_updated_ngp_field<T>(stkField);
  }

  void set_initial_part_membership(stk::mesh::BulkData & bulk, stk::mesh::Part & part1, stk::mesh::Part & part2)
  {
    stk::mesh::PartVector addParts1(1, &part1);
    stk::mesh::PartVector addParts2(1, &part2);
    stk::mesh::EntityVector allElements;

    bulk.get_entities(stk::topology::ELEM_RANK, bulk.mesh_meta_data().universal_part(), allElements);

    bulk.modification_begin();
    bulk.change_entity_parts(allElements[0], addParts1);
    bulk.change_entity_parts(allElements[1], addParts2);
    bulk.modification_end();
  }

  void modify_mesh_part_membership(stk::mesh::BulkData & bulk, const std::string & partName)
  {
    stk::mesh::Part & meshModPart = *bulk.mesh_meta_data().get_part(partName);
    stk::mesh::PartVector addParts(1, &meshModPart);
    stk::mesh::EntityVector allElements;

    bulk.get_entities(stk::topology::ELEM_RANK, bulk.mesh_meta_data().universal_part(), allElements);

    bulk.modification_begin();
    bulk.change_entity_parts(allElements[1], addParts);
    bulk.modification_end();
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
        T * fieldData = stk::mesh::field_data(stkField, entity);
        const unsigned numComponents = stk::mesh::field_scalars_per_entity(stkField, *bucket);
        for (unsigned component = 0; component < numComponents; ++component) {
          fieldData[component] = 10*id + component;
        }
      }
    }
  }

  template <typename T>
  stk::mesh::Field<T> & build_mesh_with_scalar_field(const std::string & fieldName)
  {
    stk::mesh::Field<T> & stkField = create_scalar_field<T>(stk::topology::ELEM_RANK, fieldName);
    get_meta().declare_part("meshModPart1", stk::topology::ELEM_RANK);
    get_meta().declare_part("meshModPart2", stk::topology::ELEM_RANK);
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);
    fill_initial_field<T>(stkField);
    initialize_ngp_field(stkField);
    return stkField;
  }

  template <typename T>
  stk::mesh::Field<T> & build_mesh_with_scalar_multistate_field(const std::string & fieldName)
  {
    stk::mesh::Field<T> & stkField = create_scalar_multistate_field<T>(stk::topology::ELEM_RANK, fieldName);
    get_meta().declare_part("meshModPart1", stk::topology::ELEM_RANK);
    get_meta().declare_part("meshModPart2", stk::topology::ELEM_RANK);
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);
    fill_initial_field<T>(stkField);
    initialize_ngp_field(stkField);
    return stkField;
  }

  template <typename T>
  stk::mesh::Field<T> & build_mesh_with_vector_field(const std::string & fieldName)
  {
    stk::mesh::Field<T> & stkField = create_vector_field<T>(stk::topology::ELEM_RANK, fieldName);
    get_meta().declare_part("meshModPart1", stk::topology::ELEM_RANK);
    get_meta().declare_part("meshModPart2", stk::topology::ELEM_RANK);
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);
    fill_initial_field<T>(stkField);
    initialize_ngp_field(stkField);
    return stkField;
  }

  template <typename T>
  stk::mesh::Field<T> & build_mesh_with_scalar_field_two_parts(const std::string & fieldName)
  {
    stk::mesh::Field<T> & stkField = create_scalar_field<T>(stk::topology::ELEM_RANK, fieldName);
    stk::mesh::Part & dummyPart1 = get_meta().declare_part("DummyPart1", stk::topology::ELEM_RANK);
    stk::mesh::Part & dummyPart2 = get_meta().declare_part("DummyPart2", stk::topology::ELEM_RANK);
    get_meta().declare_part("meshModPart1", stk::topology::ELEM_RANK);
    get_meta().declare_part("meshModPart2", stk::topology::ELEM_RANK);
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);
    set_initial_part_membership(get_bulk(), dummyPart1, dummyPart2);
    fill_initial_field<T>(stkField);
    initialize_ngp_field(stkField);
    return stkField;
  }

  template <typename T>
  stk::mesh::Field<T> & build_mesh_with_vector_field_two_parts(const std::string & fieldName)
  {
    stk::mesh::Field<T> & stkField = create_vector_field<T>(stk::topology::ELEM_RANK, fieldName);
    stk::mesh::Part & dummyPart1 = get_meta().declare_part("DummyPart1", stk::topology::ELEM_RANK);
    stk::mesh::Part & dummyPart2 = get_meta().declare_part("DummyPart2", stk::topology::ELEM_RANK);
    get_meta().declare_part("meshModPart1", stk::topology::ELEM_RANK);
    get_meta().declare_part("meshModPart2", stk::topology::ELEM_RANK);
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);
    set_initial_part_membership(get_bulk(), dummyPart1, dummyPart2);
    fill_initial_field<T>(stkField);
    initialize_ngp_field(stkField);
    return stkField;
  }
};


template <typename T>
void device_field_set_all(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, T value)
{
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);

  ngpField.set_all(ngpMesh, value);
}

template <typename T>
void write_scalar_field_on_host_using_entity(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, T value)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stkField.entity_rank());
  for (stk::mesh::Bucket * bucket : buckets) {
    for (const stk::mesh::Entity & entity : *bucket) {
      T * fieldData = stk::mesh::field_data(stkField, entity);
      fieldData[0] = value;
    }
  }
}

template <typename T>
void write_vector_field_on_host_using_entity(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, T value)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stkField.entity_rank());
  for (stk::mesh::Bucket * bucket : buckets) {
    for (const stk::mesh::Entity & entity : *bucket) {
      T * fieldData = stk::mesh::field_data(stkField, entity);
      fieldData[1] = value;  // Write to second component only
    }
  }
}

template <typename T>
void write_scalar_field_on_host_using_bucket(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, T value)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stkField.entity_rank());
  for (stk::mesh::Bucket * bucket : buckets) {
    T * fieldData = stk::mesh::field_data(stkField, *bucket);
    for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
      fieldData[iEntity] = value;
    }
  }
}

template <typename T>
void write_vector_field_on_host_using_bucket(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, T value)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stkField.entity_rank());
  for (stk::mesh::Bucket * bucket : buckets) {
    T * fieldData = stk::mesh::field_data(stkField, *bucket);
    for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
      const size_t yComponent = iEntity*3 + 1;
      fieldData[yComponent] = value;
    }
  }
}

template <typename T>
STK_INLINE_FUNCTION
bool dummy_check(T value)
{
  return (value == 0);
}

template <typename T>
void read_scalar_field_on_host_using_entity(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stkField.entity_rank());
  for (stk::mesh::Bucket * bucket : buckets) {
    for (const stk::mesh::Entity & entity : *bucket) {
    const T * fieldData = stk::mesh::field_data(stkField, entity);
      dummy_check(fieldData[0]);
    }
  }
}

template <typename T>
void read_scalar_field_on_host_using_bucket(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stkField.entity_rank());
  for (stk::mesh::Bucket * bucket : buckets) {
    const T * fieldData = stk::mesh::field_data(stkField, *bucket);
    for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
      dummy_check(fieldData[iEntity]);
    }
  }
}

template <typename T>
void read_vector_field_on_host_using_entity(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stkField.entity_rank());
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
void read_vector_field_on_host_using_bucket(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stkField.entity_rank());
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
void write_scalar_field_on_device(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, T value)
{
  const int component = 0;
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, meta.locally_owned_part(),
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
                                   ngpField(entity, component) = value;
                                 });
}

template <typename T>
void write_vector_field_on_device(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, T value)
{
  const int component = 1;  // Just write to the second component
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, meta.locally_owned_part(),
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
                                   ngpField(entity, component) = value;
                                 });
}

template <typename T>
void write_vector_field_on_device_using_mesh_index(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, T value)
{
  const int component = 1;  // Just write to the second component
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
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
void write_vector_field_on_device_using_entity_field_data(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, T value)
{
  const int component = 1;  // Just write to the second component
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
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
void read_scalar_field_on_device(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField)
{
  const int component = 0;
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
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
void read_scalar_field_on_device(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, NGPFIELD & ngpField)
{
  const int component = 0;
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
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

template <typename T>
void read_vector_field_on_device(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField)
{
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
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
void read_vector_field_on_device_using_mesh_index(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField)
{
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
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
void read_vector_field_on_device_using_entity_field_data(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField)
{
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
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
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_contains_file_name(stdoutString, "NgpDebugFieldSyncTest.cpp");
  check_contains_a_line_number(stdoutString);
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MissingSyncToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MissingModifyOnHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  stkField.sync_to_device();

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(get_bulk(), stkField, 3.14);

  read_vector_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceEntityFieldDataAccess_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device_using_entity_field_data(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceEntityFieldDataAccess_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(get_bulk(), stkField, 3.14);

  read_vector_field_on_device_using_entity_field_data(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceMeshIndexAccess_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device_using_mesh_index(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceMeshIndexAccess_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(get_bulk(), stkField, 3.14);

  read_vector_field_on_device_using_mesh_index(get_bulk(), stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarIntAccessUsingEntity_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_scalar_field<int>("intScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarIntAccessUsingEntity_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_scalar_field<int>("intScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3);

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field intScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingEntity_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field<int>("intVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(get_bulk(), stkField, 3);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingEntity_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field<int>("intVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(get_bulk(), stkField, 3);

  read_vector_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field intVectorField[1]=11.000000");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket(get_bulk(), stkField, 3.14);

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(get_bulk(), stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field_two_parts<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket(get_bulk(), stkField, 3.14);

  read_vector_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field_two_parts<int>("intScalarField");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket(get_bulk(), stkField, 3);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(get_bulk(), stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field_two_parts<int>("intVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket(get_bulk(), stkField, 3);

  read_vector_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field intVectorField[1]=11.000000");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(get_bulk(), stkField, 3.14);

  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

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
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  stkField.modify_on_device();

  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MissingModifyOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_device(get_bulk(), stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_vector_field_on_host_using_entity(get_bulk(), stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_device(get_bulk(), stkField, 3.14);

  read_vector_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=21");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceEntityFieldDataAccess_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_device_using_entity_field_data(get_bulk(), stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_vector_field_on_host_using_entity(get_bulk(), stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceEntityFieldDataAccess_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_device_using_entity_field_data(get_bulk(), stkField, 3.14);

  read_vector_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=21");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceMeshIndexAccess_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_device_using_mesh_index(get_bulk(), stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_vector_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceMeshIndexAccess_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_device_using_mesh_index(get_bulk(), stkField, 3.14);

  read_vector_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=21");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingEntity_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field<int>("intVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_device(get_bulk(), stkField, 3);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_vector_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingEntity_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field<int>("intVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_device(get_bulk(), stkField, 3);

  read_vector_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[1]=21");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_bucket(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(get_bulk(), stkField, 3.14);

  read_scalar_field_on_host_using_bucket(get_bulk(), stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field_two_parts<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_device(get_bulk(), stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_vector_field_on_host_using_bucket(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field_two_parts<double>("doubleVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_device(get_bulk(), stkField, 3.14);

  read_vector_field_on_host_using_bucket(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=21");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field_two_parts<int>("intVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_device(get_bulk(), stkField, 3);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_vector_field_on_host_using_bucket(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field_two_parts<int>("intVectorField");

  testing::internal::CaptureStdout();
  write_vector_field_on_device(get_bulk(), stkField, 3);

  read_vector_field_on_host_using_bucket(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[1]=21");
  check_no_warnings(stdoutString);
}



TEST_F(NgpDebugFieldSync, ScalarDeviceSetAll_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  device_field_set_all(get_bulk(), stkField, 2.18);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarDeviceSetAll_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);

  device_field_set_all(get_bulk(), stkField, 2.18);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarDeviceSetAll_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  device_field_set_all(get_bulk(), stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarDeviceSetAll_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  device_field_set_all(get_bulk(), stkField, 3.14);

  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleTimestep_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14+timeStep);
    stkField.modify_on_host();
    stkField.sync_to_device();
    read_scalar_field_on_device(get_bulk(), stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleTimestep_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14+timeStep);
    read_scalar_field_on_device(get_bulk(), stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleAccesses_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(get_bulk(), stkField);
  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleStaleAccesses_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);

  read_scalar_field_on_device(get_bulk(), stkField);
  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleWrites_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleWrites_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 2.18);

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_MultipleTimestep_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_vector_field_on_host_using_entity(get_bulk(), stkField, 3.14+timeStep);
    stkField.modify_on_host();
    stkField.sync_to_device();
    read_vector_field_on_device(get_bulk(), stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_MultipleTimestep_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_vector_field_on_host_using_entity(get_bulk(), stkField, 3.14+timeStep);
    read_vector_field_on_device(get_bulk(), stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=21.000000");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleTimestep_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_bucket(get_bulk(), stkField, 3.14+timeStep);
    stkField.modify_on_host();
    stkField.sync_to_device();
    read_scalar_field_on_device(get_bulk(), stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleTimestep_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_bucket(get_bulk(), stkField, 3.14+timeStep);
    read_scalar_field_on_device(get_bulk(), stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleAccesses_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_bucket(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(get_bulk(), stkField);
  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleStaleAccesses_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_bucket(get_bulk(), stkField, 3.14);

  read_scalar_field_on_device(get_bulk(), stkField);
  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleWrites_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_bucket(get_bulk(), stkField, 3.14);
  write_scalar_field_on_host_using_bucket(get_bulk(), stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleWrites_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_bucket(get_bulk(), stkField, 3.14);
  write_scalar_field_on_host_using_bucket(get_bulk(), stkField, 2.18);

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20.000000");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleTimestep_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(get_bulk(), stkField, 3.14+timeStep);
    stkField.modify_on_device();
    stkField.sync_to_host();
    read_scalar_field_on_host_using_entity(get_bulk(), stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleTimestep_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(get_bulk(), stkField, 3.14+timeStep);
    read_scalar_field_on_host_using_entity(get_bulk(), stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleAccesses_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(get_bulk(), stkField);
  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleStaleAccesses_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);

  read_scalar_field_on_host_using_entity(get_bulk(), stkField);
  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleWrites_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  write_scalar_field_on_device(get_bulk(), stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MultipleWrites_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  write_scalar_field_on_device(get_bulk(), stkField, 2.18);

  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleTimestep_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(get_bulk(), stkField, 3.14+timeStep);
    stkField.modify_on_device();
    stkField.sync_to_host();
    read_scalar_field_on_host_using_bucket(get_bulk(), stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleTimestep_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(get_bulk(), stkField, 3.14+timeStep);
    read_scalar_field_on_host_using_bucket(get_bulk(), stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleAccesses_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_bucket(get_bulk(), stkField);
  read_scalar_field_on_host_using_bucket(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleStaleAccesses_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);

  read_scalar_field_on_host_using_bucket(get_bulk(), stkField);
  read_scalar_field_on_host_using_bucket(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleWrites_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  write_scalar_field_on_device(get_bulk(), stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_bucket(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucket_MultipleWrites_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_two_parts<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  write_scalar_field_on_device(get_bulk(), stkField, 2.18);

  read_scalar_field_on_host_using_bucket(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  modify_mesh_part_membership(get_bulk(), "meshModPart1");

  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ImplicitSyncToDeviceWithUpdate_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  modify_mesh_part_membership(get_bulk(), "meshModPart1");

  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  modify_mesh_part_membership(get_bulk(), "meshModPart1");

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);

  read_scalar_field_on_device(get_bulk(), stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing unsynchronized Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");
  stk::mesh::NgpField<double> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  modify_mesh_part_membership(get_bulk(), "meshModPart1");

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(get_bulk(), stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing unsynchronized Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_mesh_part_membership(get_bulk(), "meshModPart1");

  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_MeshModification_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  modify_mesh_part_membership(get_bulk(), "meshModPart1");
  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoAdjacentMeshModifications_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_mesh_part_membership(get_bulk(), "meshModPart1");
  modify_mesh_part_membership(get_bulk(), "meshModPart2");

  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoAdjacentMeshModifications_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();
  modify_mesh_part_membership(get_bulk(), "meshModPart1");
  modify_mesh_part_membership(get_bulk(), "meshModPart2");

  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);

  read_scalar_field_on_device(get_bulk(), stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing unsynchronized Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoAdjacentMeshModifications_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_mesh_part_membership(get_bulk(), "meshModPart1");
  modify_mesh_part_membership(get_bulk(), "meshModPart2");

  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoAdjacentMeshModifications_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);

  modify_mesh_part_membership(get_bulk(), "meshModPart1");
  modify_mesh_part_membership(get_bulk(), "meshModPart2");

  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_mesh_part_membership(get_bulk(), "meshModPart1");
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(get_bulk(), stkField);

  modify_mesh_part_membership(get_bulk(), "meshModPart2");
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");
  stk::mesh::NgpField<double> & ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);

  testing::internal::CaptureStdout();

  modify_mesh_part_membership(get_bulk(), "meshModPart1");
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 3.14);
  read_scalar_field_on_device(get_bulk(), stkField, ngpField);

  modify_mesh_part_membership(get_bulk(), "meshModPart2");
  write_scalar_field_on_host_using_entity(get_bulk(), stkField, 2.18);
  read_scalar_field_on_device(get_bulk(), stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Accessing unsynchronized Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();
  modify_mesh_part_membership(get_bulk(), "meshModPart1");
  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  write_scalar_field_on_device(get_bulk(), stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();
  modify_mesh_part_membership(get_bulk(), "meshModPart2");
  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_TwoMeshModifications_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(get_bulk(), stkField, 3.14);
  modify_mesh_part_membership(get_bulk(), "meshModPart1");
  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  write_scalar_field_on_device(get_bulk(), stkField, 2.18);
  modify_mesh_part_membership(get_bulk(), "meshModPart2");
  read_scalar_field_on_host_using_entity(get_bulk(), stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync, ScalarAccessUsingEntity_FieldStateRotation_Throw)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_multistate_field<double>("doubleScalarField");

#if defined(STK_DEBUG_FIELD_SYNC) && defined(STK_USE_DEVICE_MESH)
  EXPECT_THROW(get_bulk().update_field_data_states(), std::logic_error);
#endif
}

}
