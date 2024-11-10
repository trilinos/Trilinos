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

class NgpDebugFieldSync : public NgpDebugFieldSyncFixture
{
public:
  template <typename T>
  void initialize_ngp_field_default_debugger(stk::mesh::Field<T> & stkField)
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
        STK_ThrowRequireMsg(get_bulk().is_valid(element), "Invalid element in fixture!");
        elemsToChange.push_back(element);
      }
      get_bulk().change_entity_parts(elemsToChange, addParts);
    }
    get_bulk().modification_end();
  }

  template <typename T>
  void write_scalar_field_on_host_using_entity_default_debugger(stk::mesh::Field<T> & stkField, T value)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>>(stkField, entity);
        fieldData[0] = value;
      }
    }
  }

  template <typename T>
  void write_vector_field_on_host_using_bucket_id_and_ordinal(stk::mesh::Field<T> & stkField, T value)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, entity);
        fieldData[1] = value;  // Write to second component only
      }
    }
  }

  template <typename T>
  void write_vector_field_on_host_using_bucket_id_and_ordinal_and_size(stk::mesh::Field<T> & stkField, T value)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        const stk::mesh::MeshIndex & meshIndex = get_bulk().mesh_index(entity);
        const unsigned bucketId = meshIndex.bucket->bucket_id();
        const stk::mesh::Bucket::size_type bucketOrd = meshIndex.bucket_ordinal;
        const unsigned numBytesPerEntity = stk::mesh::field_bytes_per_entity(stkField, *bucket);
        T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, bucketId, bucketOrd, numBytesPerEntity);
        fieldData[1] = value;  // Write to second component only
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
  void write_vector_field_on_host_using_bucket_id(stk::mesh::Field<T> & stkField, T value)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      const unsigned bucketId = bucket->bucket_id();
      T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>, StkDebugger<T>>(stkField, bucketId);
      for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
        const size_t yComponent = iEntity*3 + 1;
        fieldData[yComponent] = value;
      }
    }
  }

  template <typename T>
  void read_scalar_field_on_host_using_entity_default_debugger(stk::mesh::Field<T> & stkField)
  {
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stkField.entity_rank());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        const T * fieldData = stk::mesh::field_data<stk::mesh::Field<T>>(stkField, entity);
        access_for_memory_checking_tool(fieldData);
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

  template <typename T, template <typename> class NgpDebugger = stk::mesh::DefaultNgpFieldSyncDebugger>
  void write_scalar_host_field_on_device(stk::mesh::HostField<T, NgpDebugger> & hostField, T value)
  {
    const int component = 0;
    stk::mesh::HostMesh hostMesh(get_bulk());
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();

    stk::mesh::for_each_entity_run(hostMesh, stk::topology::ELEM_RANK, meta.locally_owned_part(),
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
                                     hostField(entity, component) = value;
                                   });
  }

  template <typename T>
  void write_scalar_field_on_device_default_debugger(stk::mesh::Field<T> & stkField, T value)
  {
    const int component = 0;
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
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
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::mesh::NgpField<T, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<T, NgpDebugger>(stkField);
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), meta.locally_owned_part());
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned i = 0; i < bucketIds.size(); ++i) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(i));
                             for (unsigned j = 0; j < bucket.size(); ++j) {
                               stk::mesh::NgpMesh::MeshIndex index{bucket.bucket_id(), static_cast<unsigned>(j)};
                               ngpField(index, component) = value;
                             }
                           }
                         });
  }

  template <typename T>
  void write_vector_field_on_device_using_entity_field_data(stk::mesh::Field<T> & stkField, T value)
  {
    const int component = 1;  // Just write to the second component
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::mesh::NgpField<T, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<T, NgpDebugger>(stkField);

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, meta.locally_owned_part(),
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
#if defined(DEVICE_USE_LOCATION_BUILTINS)
                                     stk::mesh::EntityFieldData<double> vals = ngpField(entity);
#else
                                     stk::mesh::EntityFieldData<double> vals = ngpField(entity, __FILE__, __LINE__);
#endif
                                     vals[component] = value;
                                   });
  }

  template <typename T>
  void read_scalar_field_on_device_default_debugger(stk::mesh::Field<T> & stkField)
  {
    const int component = 0;
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), meta.locally_owned_part());
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned i = 0; i < bucketIds.size(); ++i) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(i));
                             for (unsigned j = 0; j < bucket.size(); ++j) {
                               stk::mesh::FastMeshIndex index = ngpMesh.fast_mesh_index(bucket[j]);
#if defined(DEVICE_USE_LOCATION_BUILTINS)
                               access_for_memory_checking_tool(&ngpField(index, component));
#else
                               access_for_memory_checking_tool(&ngpField(index, component, __FILE__, __LINE__));
#endif
                             }
                           }
                         });
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

  template <typename T, typename NGPFIELD>
  void read_old_scalar_field_on_device(stk::mesh::Field<T> & stkField, NGPFIELD & ngpField,
                                       stk::mesh::EntityId maxIdToRead = std::numeric_limits<stk::mesh::EntityId>::max())
  {
    const int component = 0;
    const stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    const stk::mesh::MetaData & meta = get_bulk().mesh_meta_data();
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), meta.locally_owned_part());
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(unsigned ) {
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

  void delete_element(const std::vector<stk::mesh::EntityId> & elemIds)
  {
    get_bulk().modification_begin();
    for (const stk::mesh::EntityId & elemId : elemIds) {
      get_bulk().destroy_entity(get_bulk().get_entity(stk::topology::ELEM_RANK, elemId));
    }
    get_bulk().modification_end();
  }



  template <typename T>
  stk::mesh::Field<T> & create_scalar_field(stk::topology::rank_t rank, const std::string & name)
  {
    unsigned numStates = 1;
    const T init = 1;
    stk::mesh::Field<T> & field = get_meta().declare_field<T>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
    return field;
  }

  template <typename T>
  stk::mesh::Field<T> & create_scalar_multistate_field(stk::topology::rank_t rank, const std::string & name)
  {
    unsigned numStates = 2;
    const T init = 1;
    stk::mesh::Field<T> & field = get_meta().declare_field<T>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
    return field;
  }

  template <typename T>
  stk::mesh::Field<T> & create_vector_field(stk::topology::rank_t rank, const std::string & name)
  {
    unsigned numStates = 1;
    unsigned numScalarsPerEntity = 3;
    const T init[] = {1, 2, 3};
    stk::mesh::Field<T> & field = get_meta().declare_field<T>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), numScalarsPerEntity, init);
    return field;
  }

  void
  create_parts(const std::vector<std::pair<unsigned, std::string>> & numElemsInEachPart)
  {
    for (const auto & elemCountAndPart : numElemsInEachPart) {
      get_meta().declare_part_with_topology(elemCountAndPart.second, stk::topology::HEX_8);
    }
  }

  template <typename T>
  stk::mesh::Field<T> & build_mesh_with_scalar_field_default_debugger(const std::string & fieldName,
                                                                      const std::vector<std::pair<unsigned, std::string>> & numElemsInEachPart)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Field<T> & stkField = create_scalar_field<T>(stk::topology::ELEM_RANK, fieldName);
    create_parts(numElemsInEachPart);

    unsigned numElems = 0;
    for (const auto & elemCountAndPart : numElemsInEachPart) {
      numElems += elemCountAndPart.first;
    }
    stk::io::fill_mesh("generated:1x1x" + std::to_string(numElems), get_bulk());

    set_initial_part_membership(numElemsInEachPart);
    fill_initial_field<T>(stkField);
    initialize_ngp_field_default_debugger(stkField);
    return stkField;
  }

  template <typename T>
  stk::mesh::Field<T> & build_mesh_with_scalar_field(const std::string & fieldName,
                                                     const std::vector<std::pair<unsigned, std::string>> & numElemsInEachPart)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Field<T> & stkField = create_scalar_field<T>(stk::topology::ELEM_RANK, fieldName);
    create_parts(numElemsInEachPart);

    unsigned numElems = 0;
    for (const auto & elemCountAndPart : numElemsInEachPart) {
      numElems += elemCountAndPart.first;
    }

    stk::io::fill_mesh("generated:1x1x" + std::to_string(numElems), get_bulk());

    set_initial_part_membership(numElemsInEachPart);
    fill_initial_field<T>(stkField);
    initialize_ngp_field(stkField);
    return stkField;
  }

  template <typename T>
  stk::mesh::Field<T> & build_mesh_with_vector_field(const std::string & fieldName,
                                                     const std::vector<std::pair<unsigned, std::string>> & numElemsInEachPart)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Field<T> & stkField = create_vector_field<T>(stk::topology::ELEM_RANK, fieldName);
    create_parts(numElemsInEachPart);

    unsigned numElems = 0;
    for (const auto & elemCountAndPart : numElemsInEachPart) {
      numElems += elemCountAndPart.first;
    }
    stk::io::fill_mesh("generated:1x1x" + std::to_string(numElems), get_bulk());

    set_initial_part_membership(numElemsInEachPart);
    fill_initial_field<T>(stkField);
    initialize_ngp_field(stkField);
    return stkField;
  }

};



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
  check_contains_file_name(stdoutString, "NgpDebugFieldSync_Fixtures.hpp");
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
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[2]=12.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntityFieldData_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_field_on_device_using_entity_field_data(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntityFieldData_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3.14);

  read_field_on_device_using_entity_field_data(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[2]=12.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingMeshIndex_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_field_on_device_using_mesh_index(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingMeshIndex_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3.14);

  read_field_on_device_using_mesh_index(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[2]=12.000000");
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
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field intVectorField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field intVectorField[1]=11.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field intVectorField[2]=12.000000");
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

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucketId_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket_id(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucketId_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket_id(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucketId_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket_id(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucketId_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket_id(stkField, 3.14);

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

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucketIdAndOrdinal_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket_id_and_ordinal(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucketIdAndOrdinal_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket_id_and_ordinal(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucketIdAndOrdinal_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket_id_and_ordinal(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucketIdAndOrdinal_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket_id_and_ordinal(stkField, 3.14);

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucketIdAndOrdinalAndSize_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket_id_and_ordinal_and_size(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ScalarAccessUsingBucketIdAndOrdinalAndSize_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket_id_and_ordinal_and_size(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucketIdAndOrdinalAndSize_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleScalarField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket_id_and_ordinal_and_size(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_vector_field_on_device(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingBucketIdAndOrdinalAndSize_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{1, "Part1"},
                                                                                                   {1, "Part2"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_bucket_id_and_ordinal_and_size(stkField, 3.14);

  read_vector_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
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
  check_contains_file_name(stdoutString, "NgpDebugFieldSync_Fixtures.hpp");
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

  read_vector_field_on_host_using_entity<double>(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorAccessUsingEntity_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device(stkField, 3.14);

  read_vector_field_on_host_using_entity<double>(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[2]=12");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=21");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[2]=22");
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

  read_vector_field_on_host_using_entity<double>(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceEntityFieldDataAccess_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device_using_entity_field_data(stkField, 3.14);

  read_vector_field_on_host_using_entity<double>(stkField);

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

  read_vector_field_on_host_using_entity<double>(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DeviceMeshIndexAccess_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_vector_field<double>("doubleVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device_using_mesh_index(stkField, 3.14);

  read_vector_field_on_host_using_entity<double>(stkField);

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

  read_vector_field_on_host_using_entity<double>(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, VectorIntAccessUsingEntity_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<int> & stkField = build_mesh_with_vector_field<int>("intVectorField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_vector_field_on_device(stkField, 3);

  read_vector_field_on_host_using_entity<int>(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[2]=12");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[1]=21");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[2]=22");
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
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[2]=12");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[1]=21");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleVectorField[2]=22");
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
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[1]=11");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[2]=12");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[1]=21");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field intVectorField[2]=22");
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
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleVectorField[0]=10.000000");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleVectorField[2]=12.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[0]=20.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=21.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[2]=22.000000");
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


#ifndef STK_DEBUG_FIELD_SYNC
TEST_F(NgpDebugFieldSync, DefaultDebugger_ScalarAccessUsingEntity_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_default_debugger<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity_default_debugger(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device_default_debugger(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DefaultDebugger_ScalarAccessUsingEntity_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_default_debugger<double>("doubleScalarField", {{2, "Part1"}});

  testing::internal::CaptureStdout();
  write_scalar_field_on_device_default_debugger(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity_default_debugger(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DefaultDebugger_ScalarAccessUsingEntity_MeshModification_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_default_debugger<double>("doubleScalarField", {{2, "Part1"},
                                                                                                                    {1, "Part2"}});

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity_default_debugger(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device_default_debugger(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DefaultDebugger_ScalarAccessUsingEntity_MeshModification_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_default_debugger<double>("doubleScalarField", {{2, "Part1"},
                                                                                                                    {1, "Part2"}});

  testing::internal::CaptureStdout();

  write_scalar_field_on_device_default_debugger(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership({{2, "Part2", "Part1"}});

  read_scalar_field_on_host_using_entity_default_debugger(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, DefaultDebugger_HostField_UsageNotProblematic)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field_default_debugger<double>("doubleScalarField", {{2, "Part1"}});
  stk::mesh::HostField<double> hostField(get_bulk(), stkField);

  testing::internal::CaptureStdout();
  write_scalar_host_field_on_device(hostField, 3.14);
  read_scalar_field_on_host_using_entity_default_debugger(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}
#endif

TEST_F(NgpDebugFieldSync, ForcedDebugger_HostField_UsageNotProblematic_UsingEntity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});
  stk::mesh::HostField<double, NgpDebugger> hostField(get_bulk(), stkField);

  testing::internal::CaptureStdout();
  write_scalar_host_field_on_device<double, NgpDebugger>(hostField, 3.14);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync, ForcedDebugger_HostField_UsageNotProblematic_UsingBucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  stk::mesh::Field<double> & stkField = build_mesh_with_scalar_field<double>("doubleScalarField", {{2, "Part1"}});
  stk::mesh::HostField<double, NgpDebugger> hostField(get_bulk(), stkField);

  testing::internal::CaptureStdout();
  write_scalar_host_field_on_device<double, NgpDebugger>(hostField, 3.14);
  read_scalar_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

class NgpDebugFieldSync_SeparateFieldRestrictions : public NgpDebugFieldSyncFixture
{
public:
  void setup_mesh_and_field_with_multiple_restrictions(const std::string& fieldName)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Part& part1 = get_meta().declare_part_with_topology("Part1", stk::topology::HEX_8);
    stk::mesh::Part& part2 = get_meta().declare_part_with_topology("Part2", stk::topology::HEX_8);

    const unsigned numStates = 1;
    stk::mesh::Field<double> & field = get_meta().declare_field<double>(stk::topology::ELEM_RANK, fieldName, numStates);

    stk::mesh::put_field_on_mesh(field, part1, nullptr);
    stk::mesh::put_field_on_mesh(field, part2, nullptr);

    const std::vector<PartConfiguration> part1FullPart2Empty = {{"Part1", 2}};
    build_mesh(part1FullPart2Empty);
  }

  void force_creation_of_last_mod_location_field()
  {
    get_bulk().modification_begin();
    get_bulk().modification_end();
  }
};

TEST_F(NgpDebugFieldSync_SeparateFieldRestrictions, MeshMod_EmptyPart_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  
  std::string fieldName("doubleScalarField");
  setup_mesh_and_field_with_multiple_restrictions(fieldName);
  stk::mesh::Field<double> & stkField = initialized_field<double>(fieldName);

  testing::internal::CaptureStdout();
  force_creation_of_last_mod_location_field();

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

}
