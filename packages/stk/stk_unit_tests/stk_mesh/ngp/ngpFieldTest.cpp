#include <gtest/gtest.h>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
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
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include <stk_util/stk_config.h>
#include "NgpUnitTestUtils.hpp"
#include <Kokkos_Core.hpp>
#include <string>

namespace ngp_field_test {

template<typename T>
class NgpFieldTester : public stk::mesh::NgpField<T>
{
public:
  bool test_need_sync_to_host() const { return this->need_sync_to_host(); }
  bool test_need_sync_to_device() const { return this->need_sync_to_device(); }
  unsigned test_get_contiguous_bucket_offset_end(const stk::mesh::BucketVector& buckets, unsigned i)
  {
    return this->get_contiguous_bucket_offset_end(buckets, i);
  }
};

class NgpFieldFixture : public stk::unit_test_util::MeshFixture
{
public:
  template <typename T>
  stk::mesh::Field<T> & create_field(stk::topology::rank_t rank, const std::string & name, unsigned numComponent = 1)
  {
    unsigned numStates = 1;
    const std::vector<T> init(numComponent, 1);
    stk::mesh::Field<T> & field = get_meta().declare_field<stk::mesh::Field<T>>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), numComponent, init.data());
    return field;
  }

  void setup_two_field_two_element_mesh()
  {
    const unsigned bucketCapacity = 1;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);

    stk::mesh::Field<int>& stkField1 = get_meta().declare_field<stk::mesh::Field<int>>(
                                        stk::topology::ELEM_RANK, "variableLengthField");
    stk::mesh::Field<int>& stkField2 = get_meta().declare_field<stk::mesh::Field<int>>(
                                        stk::topology::ELEM_RANK, "potentiallyOverwrittenField");
    stk::mesh::Part& block1 = get_meta().declare_part_with_topology("block_1", stk::topology::SHELL_QUAD_4);
    stk::mesh::Part& block2 = get_meta().declare_part_with_topology("block_2", stk::topology::SHELL_QUAD_4);

    const int init1[2] = {1, 1};
    stk::mesh::put_field_on_mesh(stkField1, block1, 1, init1);
    stk::mesh::put_field_on_mesh(stkField1, block2, 2, init1);

    const int init2 = 1;
    stk::mesh::put_field_on_mesh(stkField2, block1, 1, &init2);
    stk::mesh::put_field_on_mesh(stkField2, block2, 1, &init2);

    const std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,5,6,block_1\n"
                                 "0,2,SHELL_QUAD_4,2,3,4,5,block_2\n";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void setup_two_element_mesh_field_on_each_element(stk::mesh::Field<int>& intField, stk::mesh::Field<double>& doubleField)
  {
    stk::mesh::Part& block1 = get_meta().declare_part_with_topology("block_1", stk::topology::SHELL_QUAD_4);
    stk::mesh::Part& block2 = get_meta().declare_part_with_topology("block_2", stk::topology::SHELL_QUAD_4);

    int initIntVal = 1;
    double initDoubleVal = 2.0;
    stk::mesh::put_field_on_mesh(intField, block1, 1, &initIntVal);
    stk::mesh::put_field_on_mesh(doubleField, block2, 1, &initDoubleVal);

    const std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,4,3,block_1\n"
                                 "0,2,SHELL_QUAD_4,2,5,6,4,block_2\n";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void setup_mesh_4hex_4block(unsigned bucketCapacity)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
    ngp_unit_test_utils::setup_mesh_4hex_4block(get_bulk(), bucketCapacity);
  }

  void setup_mesh_3hex_3block(unsigned bucketCapacity)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
    ngp_unit_test_utils::setup_mesh_3hex_3block(get_bulk(), bucketCapacity);
  }

  void setup_mesh_3hex_2block(unsigned bucketCapacity)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
    ngp_unit_test_utils::setup_mesh_3hex_2block(get_bulk(), bucketCapacity);
  }

  void setup_mesh_2hex_3block(unsigned bucketCapacity)
  {
    get_meta().declare_part("block_3", stk::topology::ELEMENT_RANK);
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
    ngp_unit_test_utils::setup_mesh_2hex_2block(get_bulk(), bucketCapacity);
  }

  void setup_mesh_2hex_2block(unsigned bucketCapacity)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
    ngp_unit_test_utils::setup_mesh_2hex_2block(get_bulk(), bucketCapacity);
  }

  template<typename T>
  void check_field_data_on_host(stk::mesh::Field<T>& stkField, unsigned multiplier)
  {
    stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;
    stk::mesh::Selector fieldSelector(stkField);
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), rank, elements);

    for(stk::mesh::Entity element : elements) {
      T* data = stk::mesh::field_data<stk::mesh::Field<T>>(stkField, element);
      unsigned numComponents = stk::mesh::field_scalars_per_entity(stkField, element);
      for(unsigned j = 0; j < numComponents; j++) {
        int expectedVal = get_bulk().identifier(element) * multiplier + j;
        EXPECT_EQ(data[j], expectedVal);
      }
    }
  }

  template<typename T, typename FieldData, typename FieldDataMirror>
  void get_field_data_from_device(stk::mesh::EntityVector& elements, stk::mesh::NgpField<T>& ngpField,
                                  FieldData& deviceData, FieldDataMirror& hostData)
  {
    unsigned numElems = elements.size();
    stk::NgpVector<stk::mesh::Entity> ngpElements(numElems);

    for(unsigned i = 0; i < numElems; i++) {
      ngpElements[i] = elements[i];
    }
    ngpElements.copy_host_to_device();

    stk::mesh::NgpMesh& ngpMesh = get_bulk().get_updated_ngp_mesh();

    Kokkos::parallel_for(numElems,
      KOKKOS_LAMBDA(const int& elemIdx) {
        stk::mesh::Entity elem = ngpElements.device_get(elemIdx);
        auto meshIndex = ngpMesh.fast_mesh_index(elem);
        unsigned numScalarsPerEntity = ngpField.get_num_components_per_entity(meshIndex);
  
        for(unsigned i = 0; i < numScalarsPerEntity; i++) { 
          deviceData(elemIdx,i) = ngpField(meshIndex,i);
        }
      }
    );
    Kokkos::deep_copy(hostData, deviceData);
  }

  template<typename T, typename FieldDataMirror, typename Func>
  void verify_field_data_on_device(const stk::mesh::EntityVector& elements, const stk::mesh::Field<T>& stkField,
                                   const FieldDataMirror& hostData, Func&& checkFunc)
  {
    for(unsigned i = 0; i < elements.size(); i++) {
      T* data = stk::mesh::field_data<stk::mesh::Field<T>>(stkField, elements[i]);
      unsigned numComponents = stk::mesh::field_scalars_per_entity(stkField, elements[i]);
      for(unsigned j = 0; j < numComponents; j++) {
        checkFunc(hostData(i,j), data[j]);
      }
    }
  }

  template<typename T, typename Func>
  void check_field_data_equality_on_device(stk::mesh::EntityVector& elements, stk::mesh::EntityRank rank,
                                           stk::mesh::NgpField<T>& ngpField, stk::mesh::Field<T>& stkField,
                                           Func&& checkFunc)
  {
    using FieldData = Kokkos::View<T**, Kokkos::LayoutRight, stk::mesh::MemSpace>;

    unsigned numElems = elements.size();
    unsigned numPerEntity = stkField.max_size(rank);
    FieldData deviceData = FieldData("deviceData", numElems, numPerEntity);
    typename FieldData::HostMirror hostData = Kokkos::create_mirror_view(deviceData);

    get_field_data_from_device<T, FieldData, typename FieldData::HostMirror>(elements, ngpField, deviceData, hostData);
    verify_field_data_on_device<T, typename FieldData::HostMirror>(elements, stkField, hostData, checkFunc);
  }

  template<typename T>
  void check_mismatched_field_data_on_device(stk::mesh::NgpField<T>& ngpField, stk::mesh::Field<T>& stkField, stk::mesh::Selector& syncSelector)
  {
    using FieldData = Kokkos::View<T**, Kokkos::LayoutRight, stk::mesh::MemSpace>;

    stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;
    stk::mesh::Selector fieldSelector(stkField);
    fieldSelector &= syncSelector;
    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(fieldSelector, get_bulk().buckets(rank), elements);
    auto checkFunc = [](T hostData, T stkData)
                     {
                       EXPECT_NE(hostData, stkData);
                     };
    check_field_data_equality_on_device<T>(elements, rank, ngpField, stkField, checkFunc);
  }

  template<typename T>
  void check_field_data_on_device(stk::mesh::NgpField<T>& ngpField, stk::mesh::Field<T>& stkField, stk::mesh::Selector& syncSelector)
  {
    using FieldData = Kokkos::View<T**, Kokkos::LayoutRight, stk::mesh::MemSpace>;

    stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;
    stk::mesh::Selector fieldSelector(stkField);
    fieldSelector &= syncSelector;
    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(fieldSelector, get_bulk().buckets(rank), elements);
    auto checkFunc = [](T hostData, T stkData)
                     {
                       EXPECT_EQ(hostData, stkData);
                     };
    check_field_data_equality_on_device<T>(elements, rank, ngpField, stkField, checkFunc);
  }

  template<typename T>
  void check_field_data_on_device(stk::mesh::NgpField<T>& ngpField, stk::mesh::Field<T>& stkField)
  {
    using FieldData = Kokkos::View<T**, Kokkos::LayoutRight, stk::mesh::MemSpace>;

    stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), rank, elements);
    auto checkFunc = [](T hostData, T stkData)
                     {
                       EXPECT_EQ(hostData, stkData);
                     };
    check_field_data_equality_on_device<T>(elements, rank, ngpField, stkField, checkFunc);
  }

  void set_element_field_data(stk::mesh::FieldBase* field)
  {
    stk::mesh::EntityVector elements;
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::ELEM_RANK);
    stk::mesh::get_selected_entities(get_meta().universal_part(), buckets, elements);

    for(stk::mesh::Entity elem : elements) {
      int* data = static_cast<int*>(stk::mesh::field_data(*field, elem));
      unsigned numComponents = stk::mesh::field_scalars_per_entity(*field, get_bulk().bucket(elem));
      for(unsigned i = 0; i < numComponents; i++) {
        data[i] = get_bulk().identifier(elem) * 10 + i;
      }
    }
  }
};

class OptimizedNgpFieldFixture : public NgpFieldFixture
{

public:
  void set_element_field_data(stk::mesh::Field<int>& stkIntField, stk::mesh::Selector selector, unsigned multiplier)
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(selector, get_bulk().buckets(stk::topology::ELEM_RANK), elements);

    for(stk::mesh::Entity elem : elements) {
      int* data = reinterpret_cast<int*>(stk::mesh::field_data(stkIntField, elem));
      unsigned numComponents = stk::mesh::field_scalars_per_entity(stkIntField, elem);
      for(unsigned j = 0; j < numComponents; j++) {
        data[j] = get_bulk().identifier(elem) * multiplier + j;
      }
    }
  }

  void set_element_field_data_on_device(stk::mesh::NgpMesh& ngpMesh, stk::mesh::Field<int>& stkIntField,
                                        const stk::mesh::Selector& selector, unsigned multiplier)
  {
    stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selector,
                                  KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entityIndex) {
                                    const int numScalarsPerEntity = ngpField.get_num_components_per_entity(entityIndex);
                                    for (int component=0; component<numScalarsPerEntity; component++) {
                                      stk::mesh::Entity entity = ngpMesh.get_entity(stk::topology::ELEM_RANK, entityIndex);
                                      ngpField(entityIndex, component) = ngpMesh.identifier(entity) * multiplier + component;
                                    }
                                  });
  }

  void setup_4hex_4block_mesh_with_field(unsigned bucketCapacity, stk::mesh::Field<int>& stkIntField)
  {
    setup_mesh_4hex_4block(bucketCapacity);
    set_element_field_data(stkIntField, get_meta().universal_part(), 10u);
  }

  void setup_3hex_3block_mesh_with_field(unsigned bucketCapacity, stk::mesh::Field<int>& stkIntField)
  {
    setup_mesh_3hex_3block(bucketCapacity);
    set_element_field_data(stkIntField, get_meta().universal_part(), 10u);
  }

  void setup_3hex_2block_mesh_with_field(unsigned bucketCapacity, stk::mesh::Field<int>& stkIntField)
  {
    setup_mesh_3hex_2block(bucketCapacity);
    set_element_field_data(stkIntField, get_meta().universal_part(), 10u);
  }

  void setup_2hex_3block_mesh_with_field(unsigned bucketCapacity, stk::mesh::Field<int>& stkIntField)
  {
    setup_mesh_2hex_3block(bucketCapacity);
    set_element_field_data(stkIntField, get_meta().universal_part(), 10u);
  }

  void setup_2hex_2block_mesh_with_field(unsigned bucketCapacity, stk::mesh::Field<int>& stkIntField)
  {
    setup_mesh_2hex_2block(bucketCapacity);
    set_element_field_data(stkIntField, get_meta().universal_part(), 10u);
  }

  void run_modify_on_host_and_device_using_part_to_select()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    ThrowRequire(part1 != nullptr && part2 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);

    stk::mesh::get_updated_ngp_field<int>(stkIntField);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 1u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 100u);

    ngpIntField.modify_on_device(*part2);
    ngpIntField.sync_to_host();

    check_field_data_on_device<int>(ngpIntField, stkIntField, selector2);
#ifdef STK_USE_DEVICE_MESH
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector1);
#endif
  }

  void run_check_contiguous_bucket_offset_copy()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_4hex_4block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    NgpFieldTester<int>& ngpIntField = static_cast<NgpFieldTester<int>&>(ngpField);

    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    stk::mesh::Part* part4 = get_meta().get_part("block_4");
    ThrowRequire(part1 != nullptr && part2 != nullptr && part4 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);
    stk::mesh::Selector selector4(*part4);

    stk::mesh::Selector blockSelector = selector1 | selector2 | selector4;
    stk::mesh::Selector selector = stk::mesh::selectField(stkIntField) & blockSelector;
    
    const stk::mesh::BucketVector& buckets = get_bulk().get_buckets(stk::topology::ELEM_RANK, selector);

    unsigned startIndex = 0;
    unsigned endIndex = ngpIntField.test_get_contiguous_bucket_offset_end(buckets, startIndex);
    unsigned numBucketsCopied = endIndex - startIndex + 1;
#ifdef STK_USE_DEVICE_MESH
    EXPECT_EQ(numBucketsCopied, 2u);
#else
    EXPECT_EQ(numBucketsCopied, 1u);
#endif
    startIndex = 1;
    endIndex = ngpIntField.test_get_contiguous_bucket_offset_end(buckets, startIndex);
    numBucketsCopied = endIndex - startIndex + 1;
    EXPECT_EQ(numBucketsCopied, 1u);

    startIndex = 2;
    endIndex = ngpIntField.test_get_contiguous_bucket_offset_end(buckets, startIndex);
    numBucketsCopied = endIndex - startIndex + 1;
    EXPECT_EQ(numBucketsCopied, 1u);
  }

  void run_sync_to_host_after_mesh_mod_and_modify_with_selector()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    ThrowRequire(part1 != nullptr && part2 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);

    stk::mesh::Entity entity = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 1u);

    get_bulk().modification_begin();
    get_bulk().change_entity_parts(entity, stk::mesh::PartVector{part2}, stk::mesh::PartVector{part1});
    get_bulk().modification_end();
    stk::mesh::get_updated_ngp_field<int>(stkIntField);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 1u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 100u);

    ngpIntField.modify_on_device(selector2);
    ngpIntField.sync_to_host();

    check_field_data_on_device<int>(ngpIntField, stkIntField, selector2);
#ifdef STK_USE_DEVICE_MESH
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector1);
#endif
  }

  void run_sync_to_host_with_selector_only_selecting_middle_block()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    stk::mesh::Part* part3 = get_meta().get_part("block_3");
    ThrowRequire(part1 != nullptr && part2 != nullptr && part3 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);
    stk::mesh::Selector selector3(*part3);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 12u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 12u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector3, 12u);

    ngpIntField.modify_on_device(selector2);
    ngpIntField.sync_to_host();

    check_field_data_on_device<int>(ngpIntField, stkIntField, selector2);
#ifdef STK_USE_DEVICE_MESH
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector1);
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector3);
#endif
  }

  void run_sync_to_host_with_selector_not_selecting_middle_block()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    stk::mesh::Part* part3 = get_meta().get_part("block_3");
    ThrowRequire(part1 != nullptr && part2 != nullptr && part3 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);
    stk::mesh::Selector selector3(*part3);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 1u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 100u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector3, 1000u);

    ngpIntField.modify_on_device(selector1);
    ngpIntField.modify_on_device(selector3);
    ngpIntField.sync_to_host();

    check_field_data_on_device<int>(ngpIntField, stkIntField, selector1);
    check_field_data_on_device<int>(ngpIntField, stkIntField, selector3);
#ifdef STK_USE_DEVICE_MESH
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector2);
#endif
  }

  void run_sync_to_host_with_selector_with_multiple_modified_on_device_sync_all()
  {
    unsigned numComponents = 2;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    ThrowRequire(part1 != nullptr && part2 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 1u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 100u);

    ngpIntField.modify_on_device(selector2);
    ngpIntField.modify_on_device();
    ngpIntField.sync_to_host();

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_sync_to_host_with_selector_with_multiple_modified_on_device4()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    ThrowRequire(part1 != nullptr && part2 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 1u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 100u);

    ngpIntField.modify_on_device();
    ngpIntField.clear_device_sync_state();
    ngpIntField.modify_on_device(selector2);
    ngpIntField.sync_to_host();

    check_field_data_on_device<int>(ngpIntField, stkIntField, selector2);
#ifdef STK_USE_DEVICE_MESH
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector1);
#endif
  }

  void run_sync_to_host_with_selector_with_multiple_modified_on_device3()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    ThrowRequire(part1 != nullptr && part2 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 1u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 100u);

    ngpIntField.modify_on_device();
    ngpIntField.clear_sync_state();
    ngpIntField.modify_on_device(selector2);
    ngpIntField.sync_to_host();

    check_field_data_on_device<int>(ngpIntField, stkIntField, selector2);
#ifdef STK_USE_DEVICE_MESH
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector1);
#endif
  }

  void run_sync_to_host_with_selector_with_multiple_modified_on_device2()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    ThrowRequire(part1 != nullptr && part2 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 1u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 100u);

    ngpIntField.modify_on_device();
    ngpIntField.modify_on_device(selector2);
    ngpIntField.sync_to_host();

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_sync_to_host_with_selector_with_multiple_modified_on_device()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    stk::mesh::Part* part3 = get_meta().get_part("block_3");
    ThrowRequire(part1 != nullptr && part2 != nullptr && part3 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);
    stk::mesh::Selector selector3(*part3);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 1u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 100u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector3, 1000u);

    ngpIntField.modify_on_device(selector1);
    ngpIntField.modify_on_device(selector2);
    ngpIntField.sync_to_host();

    stk::mesh::Selector unionSelector = selector1 | selector2;
    check_field_data_on_device<int>(ngpIntField, stkIntField, unionSelector);
#ifdef STK_USE_DEVICE_MESH
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector3);
#endif
  }

  void run_sync_to_device_after_mesh_mod_and_modify_with_selector()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    ThrowRequire(part1 != nullptr && part2 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);

    get_bulk().modification_begin();
    get_bulk().modification_end();
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 1u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 100u);

    ngpIntField.modify_on_host(selector2);
    ngpIntField.sync_to_device();

    check_field_data_on_device<int>(ngpIntField, stkIntField, selector2);
#ifdef STK_USE_DEVICE_MESH
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector1);
#endif
  }

  void run_sync_to_device_with_selector_only_selecting_middle_block()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    stk::mesh::Part* part3 = get_meta().get_part("block_3");
    ThrowRequire(part1 != nullptr && part2 != nullptr && part3 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);
    stk::mesh::Selector selector3(*part3);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 12u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 12u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector3, 12u);

    ngpIntField.modify_on_host(selector2);
    ngpIntField.sync_to_device();

    check_field_data_on_device<int>(ngpIntField, stkIntField, selector2);
#ifdef STK_USE_DEVICE_MESH
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector1);
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector3);
#endif
  }

  void run_sync_to_device_with_selector_not_selecting_middle_block()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    stk::mesh::Part* part3 = get_meta().get_part("block_3");
    ThrowRequire(part1 != nullptr && part2 != nullptr && part3 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);
    stk::mesh::Selector selector3(*part3);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 12u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 12u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector3, 12u);

    ngpIntField.modify_on_host(selector1);
    ngpIntField.modify_on_host(selector3);
    ngpIntField.sync_to_device();

    check_field_data_on_device<int>(ngpIntField, stkIntField, selector1);
    check_field_data_on_device<int>(ngpIntField, stkIntField, selector3);
#ifdef STK_USE_DEVICE_MESH
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector2);
#endif
  }

  void run_sync_to_device_with_selector_with_multiple_modified_on_host_sync_all()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    ThrowRequire(part1 != nullptr && part2 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 12u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 12u);

    ngpIntField.modify_on_host(selector2);
    ngpIntField.modify_on_host();
    ngpIntField.sync_to_device();

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_sync_to_device_with_selector_with_multiple_modified_on_host4()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    ThrowRequire(part1 != nullptr && part2 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 12u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 12u);

    ngpIntField.modify_on_host();
    ngpIntField.clear_host_sync_state();
    ngpIntField.modify_on_host(selector2);
    ngpIntField.sync_to_device();

    check_field_data_on_device<int>(ngpIntField, stkIntField, selector2);
#ifdef STK_USE_DEVICE_MESH
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector1);
#endif
  }

  void run_sync_to_device_with_selector_with_multiple_modified_on_host3()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    ThrowRequire(part1 != nullptr && part2 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 12u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 12u);

    ngpIntField.modify_on_host();
    ngpIntField.clear_sync_state();
    ngpIntField.modify_on_host(selector2);
    ngpIntField.sync_to_device();

    check_field_data_on_device<int>(ngpIntField, stkIntField, selector2);
#ifdef STK_USE_DEVICE_MESH
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector1);
#endif
  }

  void run_sync_to_device_with_selector_with_multiple_modified_on_host2()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    ThrowRequire(part1 != nullptr && part2 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 12u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 12u);

    ngpIntField.modify_on_host();
    ngpIntField.modify_on_host(selector2);
    ngpIntField.sync_to_device();

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_sync_to_device_with_selector_with_multiple_modified_on_host()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    ThrowRequire(part1 != nullptr && part2 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 12u);
    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 12u);

    ngpIntField.modify_on_host(selector1);
    ngpIntField.modify_on_host(selector2);
    ngpIntField.sync_to_device();

    stk::mesh::Selector unionSelector = selector1 | selector2;
    check_field_data_on_device<int>(ngpIntField, stkIntField, unionSelector);
    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_sync_to_device_with_selector()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part = get_meta().get_part("block_2");
    ThrowRequire(part != nullptr);
    stk::mesh::Selector selector(*part);

    modify_host_bucket_value_with_selector(stkIntField, ngpIntField, selector, 12u);

    stkIntField.modify_on_host(selector);
    stkIntField.sync_to_device();

    check_field_data_on_device<int>(ngpIntField, stkIntField, selector);
  }

  void run_sync_to_host_with_selector()
  {
    unsigned numComponents = 1;
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 3;

    setup_2hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    stk::mesh::Part* part1 = get_meta().get_part("block_1");
    stk::mesh::Part* part2 = get_meta().get_part("block_2");
    ThrowRequire(part1 != nullptr && part2 != nullptr);
    stk::mesh::Selector selector1(*part1);
    stk::mesh::Selector selector2(*part2);

    modify_device_bucket_value_with_selector(stkIntField, ngpIntField, selector1, 1u);
    modify_device_bucket_value_with_selector(stkIntField, ngpIntField, selector2, 100u);
    stkIntField.modify_on_device(selector2);
    stkIntField.sync_to_host();

    check_field_data_on_device<int>(ngpIntField, stkIntField, selector2);
#ifdef STK_USE_DEVICE_MESH
    check_mismatched_field_data_on_device<int>(ngpIntField, stkIntField, selector1);
#endif
  }

  void run_add_and_delete_bucket3(unsigned numComponents)
  { 
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 2;

    setup_2hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    modify_mesh_add_and_delete_bucket3(stkIntField, ngpIntField);

    ngpIntField.modify_on_host();
    ngpIntField.update_field();

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_add_and_delete_bucket2(unsigned numComponents)
  { 
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 2;

    setup_2hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    modify_mesh_add_and_delete_bucket2(stkIntField, ngpIntField);

    ngpIntField.modify_on_host();
    ngpIntField.update_field();

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_add_and_delete_bucket(unsigned numComponents)
  { 
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 2;

    setup_2hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    modify_mesh_add_and_delete_bucket(stkIntField, ngpIntField);

    ngpIntField.modify_on_host();
    ngpIntField.update_field();

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }


  void run_delete_bucket_in_middle(unsigned numComponents)
  {
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 2;

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    modify_mesh_delete_bucket_in_middle(stkIntField, ngpIntField);

    ngpIntField.modify_on_host();
    ngpIntField.update_field();

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_add_new_element(unsigned numComponents, unsigned bucketCapacity)
  {
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

    modify_and_test_add_element(stkIntField, ngpIntField, bucketCapacity);
  }

  void run_add_new_element_and_modify_from_device(unsigned numComponents)
  {
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 2;

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    {
      get_bulk().get_updated_ngp_mesh();
      stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

      modify_and_test_add_element(stkIntField, ngpIntField, bucketCapacity);
    }

    {
      stk::mesh::NgpMesh& ngpMesh = get_bulk().get_updated_ngp_mesh();
      stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

      check_field_data_on_host<int>(stkIntField, 10u);

      ngpIntField.sync_to_device();    
      ngpIntField.clear_sync_state();
      ngpIntField.modify_on_device();
      
      stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEMENT_RANK, get_meta().universal_part(),
        KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
          ngpIntField(elem, 0) *= 2;
        }
      );

      ngpIntField.sync_to_host();

      check_field_data_on_host<int>(stkIntField, 20u);
    }
  }

  void run_add_bucket_in_middle_copy(unsigned numComponents)
  {
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 1;

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int> ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

    modify_and_test_add_bucket_in_middle(stkIntField, ngpIntField);
  }

  void run_add_bucket_in_middle_external(unsigned numComponents)
  {
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 1;

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int> ngpIntField(get_bulk(), stkIntField);

    modify_and_test_add_bucket_in_middle(stkIntField, ngpIntField);
  }

  void run_add_bucket_in_middle_internal(unsigned numComponents)
  {
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 1;

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

    modify_and_test_add_bucket_in_middle(stkIntField, ngpIntField);
  }

  void modify_and_test_add_bucket_in_middle(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    modify_mesh_add_bucket_in_middle(stkIntField, ngpIntField);

    ngpIntField.modify_on_host();
    ngpIntField.update_field();

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void modify_and_test_add_element(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField, unsigned bucketCapacity)
  {
    modify_mesh_add_element(stkIntField, ngpIntField, bucketCapacity);

    ngpIntField.modify_on_host();
    ngpIntField.update_field();

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_change_bucket_content_by_mesh_modification(unsigned numComponents)
  {
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 2;

    setup_3hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    check_field_data_on_device<int>(ngpIntField, stkIntField);

    modify_mesh_change_bucket_content(stkIntField, ngpIntField);

    ngpIntField.clear_sync_state();
    ngpIntField.modify_on_host();
    ngpIntField.update_field();

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_change_bucket_content_by_user(unsigned numComponents)
  {
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
    const unsigned bucketCapacity = 2;

    setup_3hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    check_field_data_on_device<int>(ngpIntField, stkIntField);

    set_element_field_data(stkIntField, get_meta().universal_part(), 20u);

    ngpIntField.modify_on_host();
    ngpIntField.sync_to_device();

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void modify_mesh_add_and_delete_bucket(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    stk::mesh::NgpMesh& ngpMesh = get_bulk().get_updated_ngp_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);
    get_bulk().modification_begin();
    stk::mesh::PartVector addParts{get_meta().get_part("block_3")};
    stk::mesh::PartVector removeParts{get_meta().get_part("block_1")};
    get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 1), addParts, removeParts);
    get_bulk().modification_end();
    ngpMesh.update_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_3", {1}}, {"block_2", {2}}});
  }

  void fill_nodes(const stk::mesh::Entity element, unsigned numNodes, stk::mesh::EntityVector& nodes)
  {
    const stk::mesh::Entity* elementNodes = get_bulk().begin_nodes(element);

    for(unsigned i = 0; i < numNodes; i++) {
      nodes.push_back(elementNodes[i]);
    }
  }

  stk::mesh::Entity add_element_and_place_in_block(const std::string& newBlockName)
  {
    stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 1);
    stk::mesh::EntityVector nodes;

    unsigned numNodes = get_bulk().num_nodes(element);

    fill_nodes(element, numNodes, nodes);

    stk::mesh::Part* part = get_meta().get_part(newBlockName);
    stk::mesh::PartVector parts = {part};
    parts.push_back(&get_meta().get_topology_root_part(stk::topology::HEX_8));
    stk::mesh::Entity newElement = get_bulk().declare_element(4, parts);

    for(unsigned i = 0; i < numNodes; i++) {
      get_bulk().declare_relation(newElement, nodes[i], i);
    }

    return newElement;
  }

  void replace_element_and_place_in_block(const std::string& newBlockName)
  {
    add_element_and_place_in_block(newBlockName);

    stk::mesh::EntityVector nodes;
    stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 1);
    unsigned numNodes = get_bulk().num_nodes(element);

    fill_nodes(element, numNodes, nodes);

    for(unsigned i = 0; i < numNodes; i++) {
      get_bulk().destroy_relation(element, nodes[i], i);
    }
    get_bulk().destroy_entity(element);
  }

  void modify_host_bucket_value_with_selector(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField,
                                              stk::mesh::Selector selector, int newMultiplier)
  {
    stk::mesh::NgpMesh& ngpMesh = get_bulk().get_updated_ngp_mesh();
    check_field_data_on_device<int>(ngpIntField, stkIntField, selector);
    set_element_field_data(stkIntField, selector, newMultiplier);
    ngpMesh.update_mesh();
  }

  void modify_device_bucket_value_with_selector(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField,
                                                stk::mesh::Selector selector, int newMultiplier)
  {
    stk::mesh::NgpMesh& ngpMesh = get_bulk().get_updated_ngp_mesh();
    check_field_data_on_device<int>(ngpIntField, stkIntField, selector);
    set_element_field_data_on_device(ngpMesh, stkIntField, selector, newMultiplier);
    // ngpMesh.update_mesh();
  }

  void modify_mesh_add_and_delete_bucket3(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    stk::mesh::NgpMesh& ngpMesh = get_bulk().get_updated_ngp_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);
    get_bulk().modification_begin();
    replace_element_and_place_in_block("block_3");
    get_bulk().modification_end();
    ngpMesh.update_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_3", {4}}, {"block_2", {2}}});
  }

  void modify_mesh_add_and_delete_bucket2(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    stk::mesh::NgpMesh& ngpMesh = get_bulk().get_updated_ngp_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);
    get_bulk().modification_begin();
    stk::mesh::PartVector addParts{get_meta().get_part("block_2")};
    stk::mesh::PartVector removeParts{get_meta().get_part("block_1")};
    get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 1), addParts, removeParts);
    addParts[0] = get_meta().get_part("block_3");
    removeParts[0] = get_meta().get_part("block_2");
    get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
    get_bulk().modification_end();
    ngpMesh.update_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_3", {2}}, {"block_2", {1}}});
  }

  void modify_mesh_delete_bucket_in_middle(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    stk::mesh::NgpMesh& ngpMesh = get_bulk().get_updated_ngp_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);
    get_bulk().modification_begin();
    stk::mesh::PartVector addParts{get_meta().get_part("block_1")};
    stk::mesh::PartVector removeParts{get_meta().get_part("block_2")};
    get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
    get_bulk().modification_end();
    ngpMesh.update_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1,2}}, {"block_3", {3}}});
  }

  void modify_mesh_add_bucket_in_middle(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);
    get_bulk().modification_begin();
    stk::mesh::PartVector addParts{get_meta().get_part("block_1")};
    stk::mesh::PartVector removeParts{get_meta().get_part("block_3")};
    get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 3), addParts, removeParts);
    get_bulk().modification_end();
    ngpMesh.update_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_1", {3}}, {"block_2", {2}}});
  }

  void modify_mesh_add_element(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField, unsigned bucketCapacity)
  {
    stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);

    get_bulk().modification_begin();
    stk::mesh::Entity newElement = add_element_and_place_in_block("block_3");
    get_bulk().modification_end();

    int* data = stk::mesh::field_data(stkIntField, newElement);
    *data = get_bulk().identifier(newElement) * 10u;
    ngpMesh.update_mesh();

    if(bucketCapacity == 1) {
      ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}, {"block_3", {4}}});
    }
    else if(bucketCapacity == 2) {
      ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3,4}}});
    }
  }

  void modify_mesh_change_bucket_content(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    stk::mesh::NgpMesh& ngpMesh = get_bulk().get_updated_ngp_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), { {"block_1", {1, 2}}, {"block_3", {3}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);
    get_bulk().modification_begin();
    stk::mesh::PartVector addParts {get_meta().get_part("block_3")};
    stk::mesh::PartVector removeParts {get_meta().get_part("block_1")};
    get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
    get_bulk().modification_end();
    ngpMesh.update_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), { {"block_1", {1}}, {"block_3", {2, 3}}});
  }
};

void move_data_between_fields_on_host(const stk::mesh::BulkData & bulk,
                                      const stk::mesh::Field<int>& source,
                                      stk::mesh::Field<int>& dest)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::ELEMENT_RANK);
  stk::mesh::NgpField<int>& ngpSource = stk::mesh::get_updated_ngp_field<int>(source);
  ngpSource.sync_to_host();

  for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
  {
    const stk::mesh::Bucket &bucket = *buckets[iBucket];

    int* sourceData = static_cast<int*>(stk::mesh::field_data(source, bucket));
    int* destData   = static_cast<int*>(stk::mesh::field_data(dest, bucket));
    for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
    {
      *destData = *sourceData;
    }
  }

  stk::mesh::NgpField<int>& ngpDest = stk::mesh::get_updated_ngp_field<int>(dest);
  ngpDest.modify_on_host();
}

void test_field_values_on_device(stk::mesh::BulkData& bulk,
                                 const stk::mesh::Field<int>& stkField,
                                 const stk::mesh::Part& part,
                                 const int expectedFieldValue)
{
  stk::mesh::NgpField<int> & ngpField = stk::mesh::get_updated_ngp_field<int>(stkField);
  ngpField.sync_to_device();

  stk::mesh::Selector selection = bulk.mesh_meta_data().locally_owned_part() & part;
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stkField.entity_rank(), selection);
  const unsigned numScalarsPerEntity = stk::mesh::field_scalars_per_entity(stkField, *buckets[0]);
  stk::mesh::NgpMesh& ngpMesh = bulk.get_updated_ngp_mesh();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selection,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
                                   for (unsigned component=0; component<numScalarsPerEntity; component++) {
                                     NGP_ThrowRequire(ngpField(entity, component) == expectedFieldValue);
                                   }
                                 });
}

void test_field_values_on_device(stk::mesh::BulkData &bulk,
                                 const stk::mesh::Field<int> & stkField,
                                 const int expectedFieldValue)
{
  test_field_values_on_device(bulk, stkField, bulk.mesh_meta_data().locally_owned_part(), expectedFieldValue);
}

void test_field_values_on_host(const stk::mesh::BulkData& bulk,
                               const stk::mesh::Field<int>& stkField,
                               const stk::mesh::Part& part,
                               const int expectedFieldValue)
{
  stkField.sync_to_host();

  stk::mesh::Selector selection = bulk.mesh_meta_data().locally_owned_part() & part;
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stkField.entity_rank(), selection);
  for (size_t iBucket=0; iBucket<buckets.size(); iBucket++) {
    const stk::mesh::Bucket &bucket = *buckets[iBucket];
    const unsigned numScalarsPerEntity = stk::mesh::field_scalars_per_entity(stkField, bucket);

    int* fieldData = reinterpret_cast<int*>(stk::mesh::field_data(stkField, bucket));
    for (size_t iEntity=0; iEntity<bucket.size(); iEntity++) {
      for (unsigned component=0; component<numScalarsPerEntity; component++) {
        EXPECT_EQ(expectedFieldValue, fieldData[component]);
      }
    }
  }
}

void test_field_values_on_host(const stk::mesh::BulkData& bulk,
                               const stk::mesh::Field<int>& stkField,
                               const int expectedFieldValue)
{
  test_field_values_on_host(bulk, stkField, bulk.mesh_meta_data().locally_owned_part(), expectedFieldValue);
}

template <typename T>
void initialize_ngp_field(stk::mesh::Field<T> & stkField)
{
  stk::mesh::get_updated_ngp_field<T>(stkField);
}

template <typename T>
void modify_field_on_host(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & field, int multiplier)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(field.entity_rank());
  for (stk::mesh::Bucket * bucket : buckets) {
    T * fieldData = stk::mesh::field_data(field, *bucket);
    for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
      fieldData[iEntity] *= multiplier;
    }
  }
  field.modify_on_host();
}

void modify_mesh_part_membership(stk::mesh::BulkData & bulk, stk::mesh::Part & part)
{
  stk::mesh::PartVector addParts(1, &part);
  stk::mesh::EntityVector allElements;

  stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, allElements);

  bulk.modification_begin();
  bulk.change_entity_parts(allElements[0], addParts);
  bulk.modification_end();
}

template <typename T>
void sync_field_to_device(stk::mesh::Field<T> & field)
{
  field.sync_to_device();
}

template <typename T>
void modify_field_on_device(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, int multiplier)
{
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);
  ngpField.sync_to_device();

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, meta.locally_owned_part(),
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
                                   const int numScalarsPerEntity = ngpField.get_num_components_per_entity(entity);
                                   for (int component=0; component<numScalarsPerEntity; component++) {
                                     ngpField(entity, component) *= multiplier;
                                   }
                                 });
  stkField.modify_on_device();
}

template <typename T>
void modify_field_on_device(stk::mesh::BulkData& bulk, stk::mesh::Field<T>& stkField, stk::mesh::Part& part, int value)
{
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);
  ngpField.sync_to_device();

  stk::mesh::Selector selection = meta.locally_owned_part() & part;
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stkField.entity_rank(), selection);
  const unsigned numScalarsPerEntity = stk::mesh::field_scalars_per_entity(stkField, *buckets[0]);
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selection,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
                                   for (unsigned component=0; component<numScalarsPerEntity; component++) {
                                     ngpField(entity, component) = value;
                                   }
                                 });
  stkField.modify_on_device();
}

template <typename T>
void sync_field_to_host(stk::mesh::Field<T> & stkField)
{
  stkField.sync_to_host();
}

template <typename T>
void check_field_on_host(const stk::mesh::BulkData & bulk,
                         const stk::mesh::Field<T>& stkField,
                         int expectedValue)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stkField.entity_rank());
  for (stk::mesh::Bucket * bucket : buckets) {
    T * fieldData = stk::mesh::field_data(stkField, *bucket);
    for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
      EXPECT_EQ(fieldData[iEntity], expectedValue);
    }
  }
}

TEST_F(NgpFieldFixture, noOverwriteInVariableLengthFields)
{
  if (get_parallel_size() != 1) return;

  setup_two_field_two_element_mesh();

  stk::mesh::Field<int>& variableLengthField = dynamic_cast<stk::mesh::Field<int>&>(
                          *get_meta().get_field(stk::topology::ELEM_RANK, "variableLengthField"));
  stk::mesh::Field<int>& potentiallyOverwrittenField = dynamic_cast<stk::mesh::Field<int>&>(
                          *get_meta().get_field(stk::topology::ELEM_RANK, "potentiallyOverwrittenField"));

  modify_field_on_device(get_bulk(), variableLengthField, 2);
  modify_field_on_host(get_bulk(), potentiallyOverwrittenField, 3);
  variableLengthField.sync_to_host();

  test_field_values_on_host(get_bulk(), potentiallyOverwrittenField, 3);
}

void check_expected_num_elements(const stk::mesh::BulkData & bulk, unsigned numElements)
{
  std::vector<size_t> counts;
  const stk::mesh::Part & locallyOwnedPart = bulk.mesh_meta_data().locally_owned_part();
  stk::mesh::count_entities(locallyOwnedPart, bulk, counts);
  ASSERT_EQ(counts[stk::topology::ELEM_RANK], numElements);
}

void fill_gold_num_scalars_per_entity(const stk::mesh::BulkData & bulk, const stk::mesh::FieldBase & variableLengthField,
                                      unsigned numElements, unsigned goldNumScalarsPerEntity[])
{
  stk::mesh::EntityVector elements;
  stk::mesh::get_selected_entities(bulk.mesh_meta_data().locally_owned_part(), bulk.buckets(stk::topology::ELEM_RANK), elements);
  for (const stk::mesh::Entity & element : elements) {
    const unsigned index = bulk.identifier(element) - 1;
    ASSERT_LT(index, numElements);
    goldNumScalarsPerEntity[index] = stk::mesh::field_scalars_per_entity(variableLengthField, element);
  }
}

void test_num_scalars_per_entity(stk::mesh::BulkData & bulk, const stk::mesh::FieldBase & variableLengthField)
{
  const unsigned numElements = 2;
  check_expected_num_elements(bulk, numElements);

  unsigned goldNumScalarsPerEntity[numElements];
  fill_gold_num_scalars_per_entity(bulk, variableLengthField, numElements, goldNumScalarsPerEntity);

  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  stk::mesh::NgpField<int> ngpVariableLengthField = stk::mesh::get_updated_ngp_field<int>(variableLengthField);

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, bulk.mesh_meta_data().locally_owned_part(),
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex & index) {
                                   const stk::mesh::Entity element = ngpMesh.get_entity(stk::topology::ELEM_RANK, index);
                                   const unsigned goldIndex = ngpMesh.identifier(element) - 1;
                                   NGP_EXPECT_EQ(ngpVariableLengthField.get_num_components_per_entity(index),
                                                 goldNumScalarsPerEntity[goldIndex]);
                                 });
}

TEST_F(NgpFieldFixture, NumScalarsPerEntityOnDevice)
{
  if (get_parallel_size() != 1) return;

  setup_two_field_two_element_mesh();

  stk::mesh::Field<int>& variableLengthField = dynamic_cast<stk::mesh::Field<int>&>(
                                               *get_meta().get_field(stk::topology::ELEM_RANK, "variableLengthField"));

  test_num_scalars_per_entity(get_bulk(), variableLengthField);
}

TEST_F(NgpFieldFixture, TestAriaAlgorithm)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & f1 = create_field<int>(stk::topology::ELEM_RANK, "f1");
  stk::mesh::Field<int> & f2 = create_field<int>(stk::topology::ELEM_RANK, "f2");

  setup_mesh("generated:1x1x1", stk::mesh::BulkData::AUTO_AURA);

  stk::mesh::NgpMesh ngpMesh = get_bulk().get_updated_ngp_mesh();

  int multiplier = 2;
  modify_field_on_device(get_bulk(), f1, multiplier);

  move_data_between_fields_on_host(get_bulk(), f1, f2);

  test_field_values_on_device(get_bulk(), f1, multiplier);
  test_field_values_on_device(get_bulk(), f2, multiplier);

  test_field_values_on_host(get_bulk(), f1, multiplier);
  test_field_values_on_host(get_bulk(), f2, multiplier);
}

TEST_F(NgpFieldFixture, GetNgpField)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::NgpField<int> & ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  EXPECT_EQ(ngpIntField.get_ordinal(), stkIntField.mesh_meta_data_ordinal());

#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE((std::is_same<decltype(ngpIntField), stk::mesh::DeviceField<int>&>::value));
#else
  EXPECT_TRUE((std::is_same<decltype(ngpIntField), stk::mesh::HostField<int>&>::value));
#endif
}

TEST_F(NgpFieldFixture, GetMultipleNgpFields)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & stkIntField1 = create_field<int>(stk::topology::ELEM_RANK, "intField1");
  stk::mesh::Field<int> & stkIntField2 = create_field<int>(stk::topology::ELEM_RANK, "intField2");
  stk::mesh::Field<double> & stkDoubleField1 = create_field<double>(stk::topology::ELEM_RANK, "doubleField1");
  stk::mesh::Field<double> & stkDoubleField2 = create_field<double>(stk::topology::ELEM_RANK, "doubleField2");

  setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::NgpField<int> & ngpIntField1 = stk::mesh::get_updated_ngp_field<int>(stkIntField1);
  stk::mesh::NgpField<int> & ngpIntField2 = stk::mesh::get_updated_ngp_field<int>(stkIntField2);
  stk::mesh::NgpField<double> & ngpDoubleField1 = stk::mesh::get_updated_ngp_field<double>(stkDoubleField1);
  stk::mesh::NgpField<double> & ngpDoubleField2 = stk::mesh::get_updated_ngp_field<double>(stkDoubleField2);

  EXPECT_EQ(ngpIntField1.get_ordinal(), stkIntField1.mesh_meta_data_ordinal());
  EXPECT_EQ(ngpIntField2.get_ordinal(), stkIntField2.mesh_meta_data_ordinal());
  EXPECT_EQ(ngpDoubleField1.get_ordinal(), stkDoubleField1.mesh_meta_data_ordinal());
  EXPECT_EQ(ngpDoubleField2.get_ordinal(), stkDoubleField2.mesh_meta_data_ordinal());
}

TEST_F(NgpFieldFixture, ModifyAndSync)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_ngp_field(stkIntField);

  int multiplier = 2;
  modify_field_on_host(get_bulk(), stkIntField, multiplier);
  check_field_on_host(get_bulk(), stkIntField, multiplier);

  sync_field_to_device(stkIntField);
  modify_field_on_device(get_bulk(), stkIntField, multiplier);

  sync_field_to_host(stkIntField);
  check_field_on_host(get_bulk(), stkIntField, multiplier*multiplier);

  size_t expectedSyncsToDevice = 1;
  size_t expectedSyncsToHost = 1;

  EXPECT_EQ(expectedSyncsToDevice, stkIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, stkIntField.num_syncs_to_host());

  stk::mesh::NgpField<int>& deviceNgpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  stk::mesh::HostField<int> hostNgpIntField(get_bulk(), stkIntField);

  EXPECT_EQ(expectedSyncsToDevice, deviceNgpIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToDevice, hostNgpIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, deviceNgpIntField.num_syncs_to_host());
  EXPECT_EQ(expectedSyncsToHost, hostNgpIntField.num_syncs_to_host());
}

TEST_F(NgpFieldFixture, UpdateNgpFieldAfterMeshMod_WithMostCurrentDataOnHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::mesh::Part & dummyPart = get_meta().declare_part("DummyPart", stk::topology::ELEM_RANK);
  setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_ngp_field(stkIntField);

  int multiplier = 2;
  modify_field_on_host(get_bulk(), stkIntField, multiplier);
  modify_mesh_part_membership(get_bulk(), dummyPart);
  check_field_on_host(get_bulk(), stkIntField, multiplier);

  sync_field_to_device(stkIntField);
  modify_field_on_device(get_bulk(), stkIntField, multiplier);

  sync_field_to_host(stkIntField);
  check_field_on_host(get_bulk(), stkIntField, multiplier*multiplier);

  size_t expectedSyncsToDevice = 1;
  size_t expectedSyncsToHost = 1;

  EXPECT_EQ(expectedSyncsToDevice, stkIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, stkIntField.num_syncs_to_host());

  stk::mesh::NgpField<int>& deviceNgpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  stk::mesh::HostField<int> hostNgpIntField(get_bulk(), stkIntField);

  EXPECT_EQ(expectedSyncsToDevice, deviceNgpIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToDevice, hostNgpIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, deviceNgpIntField.num_syncs_to_host());
  EXPECT_EQ(expectedSyncsToHost, hostNgpIntField.num_syncs_to_host());
}

TEST_F(NgpFieldFixture, UpdateNgpFieldAfterMeshMod_WithMostCurrentDataOnDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::mesh::Part & dummyPart = get_meta().declare_part("DummyPart", stk::topology::ELEM_RANK);
  setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_ngp_field(stkIntField);

  int multiplier = 2;
  modify_field_on_device(get_bulk(), stkIntField, multiplier);
  modify_mesh_part_membership(get_bulk(), dummyPart);

  sync_field_to_host(stkIntField);
}

TEST_F(NgpFieldFixture, ConsistentNeedToSyncForAllCopies)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  get_meta().declare_part("DummyPart", stk::topology::ELEM_RANK);
  setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_ngp_field(stkIntField);
  stk::mesh::NgpField<int> copyNgpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

  int multiplier = 2;
  modify_field_on_host(get_bulk(), stkIntField, multiplier);
  check_field_on_host(get_bulk(), stkIntField, multiplier);

  copyNgpField.sync_to_device();
  modify_field_on_device(get_bulk(), stkIntField, multiplier);

  copyNgpField.sync_to_host();
  check_field_on_host(get_bulk(), stkIntField, multiplier*multiplier);

  size_t expectedSyncsToDevice = 1;
  size_t expectedSyncsToHost = 1;

  EXPECT_EQ(expectedSyncsToDevice, stkIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, stkIntField.num_syncs_to_host());

  EXPECT_EQ(stkIntField.num_syncs_to_device(), copyNgpField.num_syncs_to_device());
  EXPECT_EQ(stkIntField.num_syncs_to_host(), copyNgpField.num_syncs_to_host());
}

TEST_F(NgpFieldFixture, ConsistentNeedToSyncForAllCopiesNoModifyCall)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  get_meta().declare_part("DummyPart", stk::topology::ELEM_RANK);
  setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_ngp_field(stkIntField);
  stk::mesh::NgpField<int> copyNgpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

  int multiplier = 1;
  check_field_on_host(get_bulk(), stkIntField, multiplier);
  copyNgpField.sync_to_device();
  copyNgpField.sync_to_host();
  check_field_on_host(get_bulk(), stkIntField, multiplier);

  size_t expectedSyncsToDevice = 0;
  size_t expectedSyncsToHost = 0;

  EXPECT_EQ(expectedSyncsToDevice, stkIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, stkIntField.num_syncs_to_host());

  EXPECT_EQ(stkIntField.num_syncs_to_device(), copyNgpField.num_syncs_to_device());
  EXPECT_EQ(stkIntField.num_syncs_to_host(), copyNgpField.num_syncs_to_host());
}

TEST_F(NgpFieldFixture, RequireSyncBetweenModifyOnHostAndDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  get_meta().declare_part("DummyPart", stk::topology::ELEM_RANK);
  setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_ngp_field(stkIntField);
  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

  int multiplier = 2;
  modify_field_on_host(get_bulk(), stkIntField, multiplier);
  check_field_on_host(get_bulk(), stkIntField, multiplier);
#ifdef KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK
  EXPECT_THROW(ngpField.modify_on_device(), std::logic_error);
#endif

  stkIntField.clear_sync_state();
  ngpField.modify_on_device();
#ifdef KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK
  EXPECT_THROW(ngpField.modify_on_host(), std::logic_error);
#endif
}

TEST_F(NgpFieldFixture, ClearSyncStateAfterModifyOnDevice)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  NgpFieldTester<int>& testNgpField = static_cast<NgpFieldTester<int>&>(ngpField);

  testNgpField.modify_on_device();

  EXPECT_TRUE(testNgpField.test_need_sync_to_host());

  stkIntField.clear_sync_state();
  EXPECT_FALSE(testNgpField.test_need_sync_to_host());
}

TEST_F(NgpFieldFixture, ClearHostSyncState)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  NgpFieldTester<int>& testNgpField = static_cast<NgpFieldTester<int>&>(ngpField);

  testNgpField.modify_on_host();

  EXPECT_TRUE(testNgpField.test_need_sync_to_device());

  testNgpField.clear_host_sync_state();

  EXPECT_FALSE(testNgpField.test_need_sync_to_device());
}

TEST_F(NgpFieldFixture, ClearDeviceSyncState)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  NgpFieldTester<int>& testNgpField = static_cast<NgpFieldTester<int>&>(ngpField);

  testNgpField.modify_on_device();

  EXPECT_TRUE(testNgpField.test_need_sync_to_host());

  testNgpField.clear_device_sync_state();

  EXPECT_FALSE(testNgpField.test_need_sync_to_host());
}

TEST_F(NgpFieldFixture, ClearHostSyncState_doesntClearDeviceMod)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  NgpFieldTester<int>& testNgpField = static_cast<NgpFieldTester<int>&>(ngpField);

  testNgpField.modify_on_device();

  EXPECT_TRUE(testNgpField.test_need_sync_to_host());

  testNgpField.clear_host_sync_state();

  EXPECT_TRUE(testNgpField.test_need_sync_to_host());
}

TEST_F(NgpFieldFixture, ClearDeviceSyncState_doesntClearHostMod)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  NgpFieldTester<int>& testNgpField = static_cast<NgpFieldTester<int>&>(ngpField);

  testNgpField.modify_on_host();

  EXPECT_TRUE(testNgpField.test_need_sync_to_device());

  testNgpField.clear_device_sync_state();

  EXPECT_TRUE(testNgpField.test_need_sync_to_device());
}

class NgpFieldSwapFixture : public NgpFieldFixture {
public:
  void setup_fields_for_swap() {
    stk::mesh::Field<int>& stkIntField1 = create_field<int>(stk::topology::ELEM_RANK, "intField1");
    stk::mesh::Field<int>& stkIntField2 = create_field<int>(stk::topology::ELEM_RANK, "intField2");
    setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);

    stk::mesh::NgpField<int>& ngpField1 = stk::mesh::get_updated_ngp_field<int>(stkIntField1);
    stk::mesh::NgpField<int>& ngpField2 = stk::mesh::get_updated_ngp_field<int>(stkIntField2);
    testNgpField1 = static_cast<NgpFieldTester<int>&>(ngpField1);
    testNgpField2 = static_cast<NgpFieldTester<int>&>(ngpField2);
  }

protected:
  NgpFieldTester<int> testNgpField1;
  NgpFieldTester<int> testNgpField2;

};

TEST_F(NgpFieldSwapFixture, SwapSyncState_ModFlagsUnset)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  setup_fields_for_swap();

  EXPECT_FALSE(testNgpField1.test_need_sync_to_host());
  EXPECT_FALSE(testNgpField2.test_need_sync_to_host());
  EXPECT_FALSE(testNgpField1.test_need_sync_to_device());
  EXPECT_FALSE(testNgpField2.test_need_sync_to_device());

  testNgpField1.swap(testNgpField2);

  EXPECT_FALSE(testNgpField1.test_need_sync_to_host());
  EXPECT_FALSE(testNgpField2.test_need_sync_to_host());
  EXPECT_FALSE(testNgpField1.test_need_sync_to_device());
  EXPECT_FALSE(testNgpField2.test_need_sync_to_device());
}

TEST_F(NgpFieldSwapFixture, SwapSyncState_ModFlagsSetModDevice)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  setup_fields_for_swap();

  testNgpField1.modify_on_device();

  EXPECT_TRUE(testNgpField1.test_need_sync_to_host());
  EXPECT_FALSE(testNgpField2.test_need_sync_to_host());
  EXPECT_FALSE(testNgpField1.test_need_sync_to_device());
  EXPECT_FALSE(testNgpField2.test_need_sync_to_device());

  testNgpField1.swap(testNgpField2);

  EXPECT_FALSE(testNgpField1.test_need_sync_to_host());
  EXPECT_TRUE(testNgpField2.test_need_sync_to_host());
  EXPECT_FALSE(testNgpField1.test_need_sync_to_device());
  EXPECT_FALSE(testNgpField2.test_need_sync_to_device());
}

TEST_F(NgpFieldSwapFixture, SwapSyncState_ModFlagsSetModHost)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  setup_fields_for_swap();

  testNgpField2.modify_on_host();

  EXPECT_FALSE(testNgpField1.test_need_sync_to_host());
  EXPECT_FALSE(testNgpField2.test_need_sync_to_host());
  EXPECT_FALSE(testNgpField1.test_need_sync_to_device());
  EXPECT_TRUE(testNgpField2.test_need_sync_to_device());

  testNgpField1.swap(testNgpField2);

  EXPECT_FALSE(testNgpField1.test_need_sync_to_host());
  EXPECT_FALSE(testNgpField2.test_need_sync_to_host());
  EXPECT_TRUE(testNgpField1.test_need_sync_to_device());
  EXPECT_FALSE(testNgpField2.test_need_sync_to_device());
}

TEST_F(NgpFieldSwapFixture, SwapSyncState_ModFlagsSetModHostDevice)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  setup_fields_for_swap();

  testNgpField1.modify_on_host();
  testNgpField2.modify_on_device();

  EXPECT_FALSE(testNgpField1.test_need_sync_to_host());
  EXPECT_TRUE(testNgpField2.test_need_sync_to_host());
  EXPECT_TRUE(testNgpField1.test_need_sync_to_device());
  EXPECT_FALSE(testNgpField2.test_need_sync_to_device());

  testNgpField1.swap(testNgpField2);

  EXPECT_TRUE(testNgpField1.test_need_sync_to_host());
  EXPECT_FALSE(testNgpField2.test_need_sync_to_host());
  EXPECT_FALSE(testNgpField1.test_need_sync_to_device());
  EXPECT_TRUE(testNgpField2.test_need_sync_to_device());
}

TEST_F(OptimizedNgpFieldFixture, ChangeBucketContentsByUserWithSingleComponent)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_change_bucket_content_by_user(numComponents);
}

//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |   3   |
//   |block_1|block_1|block_3|        |block_1|block_3|block_3|
//   -------------------------        -------------------------
//
TEST_F(OptimizedNgpFieldFixture, ChangeBucketContentsUsingLayoutModificationWithSingleComponent)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_change_bucket_content_by_mesh_modification(numComponents);
}

TEST_F(OptimizedNgpFieldFixture, ChangeBucketContentsUsingLayoutModificationWithThreeComponents)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 3;

  run_change_bucket_content_by_mesh_modification(numComponents);
}

void construct_ngp_field(stk::mesh::Field<int>& stkIntField)
{
  EXPECT_NO_THROW(stk::mesh::NgpField<int> ngpField(stkIntField.get_mesh(), stkIntField));
}

TEST_F(OptimizedNgpFieldFixture, DISABLED_CreateConsecutiveNgpFields)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);
  const unsigned bucketCapacity = 2;

  setup_3hex_2block_mesh_with_field(bucketCapacity, stkIntField);

  construct_ngp_field(stkIntField);
  construct_ngp_field(stkIntField);
}

//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |   3   |
//   |block_1|block_2|block_3|        |block_1|block_2|block_1|
//   -------------------------        -------------------------
//

TEST_F(OptimizedNgpFieldFixture, AddBucketInMiddleWithSingleComponentInternal)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_add_bucket_in_middle_internal(numComponents);
}

TEST_F(OptimizedNgpFieldFixture, AddBucketInMiddleWithSingleComponentCopy)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_add_bucket_in_middle_copy(numComponents);
}

//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |  3,4  |
//   |block_1|block_2|block_3|        |block_1|block_1|block_3|
//   -------------------------        -------------------------
//

TEST_F(OptimizedNgpFieldFixture, AddNewElementWithBucketCapacity2)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;
  unsigned bucketCapacity = 2;

  run_add_new_element(numComponents, bucketCapacity);
}

TEST_F(OptimizedNgpFieldFixture, AddNewElementAndModifyFromDevice)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_add_new_element_and_modify_from_device(numComponents);
}

TEST_F(OptimizedNgpFieldFixture, AddNewElementWithBucketCapacity1)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;
  unsigned bucketCapacity = 1;

  run_add_new_element(numComponents, bucketCapacity);
}

//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |   3   |
//   |block_1|block_2|block_3|        |block_1|block_1|block_3|
//   -------------------------        -------------------------
//
TEST_F(OptimizedNgpFieldFixture, DeleteBucketInMiddleWithSingleComponent)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_delete_bucket_in_middle(numComponents);
}

//   -----------------        -----------------
//   |       |       |        |       |       |
//   |   1   |   2   |  ===>  |   1   |   2   |
//   |block_1|block_2|        |block_3|block_1|
//   -----------------        -----------------
//

TEST_F(OptimizedNgpFieldFixture, AddAndDeleteBucketWithSingleComponent)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_add_and_delete_bucket(numComponents);
}

//   -----------------        -----------------
//   |       |       |        |       |       |
//   |   1   |   2   |  ===>  |   1   |   2   |
//   |block_1|block_2|        |block_2|block_3|
//   -----------------        -----------------
//

TEST_F(OptimizedNgpFieldFixture, AddAndDeleteBucketWithSingleComponent2)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_add_and_delete_bucket2(numComponents);
}

//   -----------------        -----------------
//   |       |       |        |       |       |
//   |   1   |   2   |  ===>  |   2   |   3   |
//   |block_1|block_2|        |block_2|block_3|
//   -----------------        -----------------
//

TEST_F(OptimizedNgpFieldFixture, AddAndDeleteBucketWithSingleComponent3)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_add_and_delete_bucket3(numComponents);
}

TEST_F(OptimizedNgpFieldFixture, SyncToDeviceWithSelector)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_device_with_selector();
}

TEST_F(OptimizedNgpFieldFixture, SyncToDeviceWithSelectorMultipleModifiedOnHost)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_device_with_selector_with_multiple_modified_on_host();
}

TEST_F(OptimizedNgpFieldFixture, SyncToDeviceWithSelectorMultipleModifiedOnHost2)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_device_with_selector_with_multiple_modified_on_host2();
}

TEST_F(OptimizedNgpFieldFixture, SyncToDeviceWithSelectorMultipleModifiedOnHost3)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_device_with_selector_with_multiple_modified_on_host3();
}

TEST_F(OptimizedNgpFieldFixture, SyncToDeviceWithSelectorMultipleModifiedOnHost4)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_device_with_selector_with_multiple_modified_on_host4();
}

TEST_F(OptimizedNgpFieldFixture, SyncToDeviceWithSelectorMultipleModifiedOnHostSyncAll)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_device_with_selector_with_multiple_modified_on_host_sync_all();
}

TEST_F(OptimizedNgpFieldFixture, SyncToDeviceWithSelectorNotSelectingMiddleBlock)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_device_with_selector_not_selecting_middle_block();
}

TEST_F(OptimizedNgpFieldFixture, SyncToDeviceWithSelectorOnlySelectingMiddleBlock)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_device_with_selector_only_selecting_middle_block();
}

TEST_F(OptimizedNgpFieldFixture, SyncToDeviceAfterMeshModAndModifyWithSelector)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_device_after_mesh_mod_and_modify_with_selector();
}

TEST_F(OptimizedNgpFieldFixture, SyncToHostWithSelector)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_host_with_selector();
}

TEST_F(OptimizedNgpFieldFixture, SyncToHostWithSelectorMultipleModifiedOnDevice)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_host_with_selector_with_multiple_modified_on_device();
}

TEST_F(OptimizedNgpFieldFixture, SyncToHostWithSelectorMultipleModifiedOnDevice2)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_host_with_selector_with_multiple_modified_on_device2();
}

TEST_F(OptimizedNgpFieldFixture, SyncToHostWithSelectorMultipleModifiedOnDevice3)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_host_with_selector_with_multiple_modified_on_device3();
}

TEST_F(OptimizedNgpFieldFixture, SyncToHostWithSelectorMultipleModifiedOnDevice4)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_host_with_selector_with_multiple_modified_on_device4();
}

TEST_F(OptimizedNgpFieldFixture, SyncToHostWithSelectorMultipleModifiedOnDeviceSyncAll)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_host_with_selector_with_multiple_modified_on_device_sync_all();
}

TEST_F(OptimizedNgpFieldFixture, SyncToHostWithSelectorNotSelectingMiddleBlock)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_host_with_selector_not_selecting_middle_block();
}

TEST_F(OptimizedNgpFieldFixture, SyncToHostWithSelectorOnlySelectingMiddleBlock)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_host_with_selector_only_selecting_middle_block();
}

TEST_F(OptimizedNgpFieldFixture, SyncToHostAfterMeshModAndModifyWithSelector)
{
  if (get_parallel_size() != 1) return;

  run_sync_to_host_after_mesh_mod_and_modify_with_selector();
}

TEST_F(OptimizedNgpFieldFixture, CheckContiguousBucketOffsetCopy)
{
  if (get_parallel_size() != 1) return;

  run_check_contiguous_bucket_offset_copy();
}

TEST_F(OptimizedNgpFieldFixture, ModifyOnHostAndDeviceUsingPartToSelect)
{
  if (get_parallel_size() != 1) return;

  run_modify_on_host_and_device_using_part_to_select();
}

}
