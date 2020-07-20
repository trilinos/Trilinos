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
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include <stk_util/stk_config.h>
#include "NgpUnitTestUtils.hpp"
#include <Kokkos_Core.hpp>
#include <string>

namespace ngp_field_test {

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
    stk::mesh::Part& block1 = get_meta().declare_part("block_1", stk::topology::ELEM_RANK);
    stk::mesh::Part& block2 = get_meta().declare_part("block_2", stk::topology::ELEM_RANK);

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
    ngp_unit_test_utils::setup_mesh_2hex_3block(get_bulk(), bucketCapacity);
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

  template<typename T>
  void check_field_data_on_device(stk::mesh::NgpField<T>& ngpField, stk::mesh::Field<T>& stkField)
  {
    using FieldData = Kokkos::View<T**, Kokkos::LayoutRight, stk::mesh::MemSpace>;

    stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;
    stk::mesh::Selector fieldSelector(stkField);
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), rank, elements);
    unsigned numElems = elements.size();
    stk::NgpVector<stk::mesh::Entity> ngpElements(elements.size());

    for(unsigned i = 0; i < elements.size(); i++) {
      ngpElements[i] = elements[i];
    }
    ngpElements.copy_host_to_device();

    unsigned numPerEntity = stkField.max_size(rank);

    FieldData deviceData = FieldData("deviceData", numElems, numPerEntity);
    typename FieldData::HostMirror hostData = Kokkos::create_mirror_view(deviceData);

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

    for(unsigned i = 0; i < numElems; i++) {
      T* data = stk::mesh::field_data<stk::mesh::Field<T>>(stkField, elements[i]);
      unsigned numComponents = stk::mesh::field_scalars_per_entity(stkField, elements[i]);
      for(unsigned j = 0; j < numComponents; j++) {
        EXPECT_EQ(hostData(i,j), data[j]) << "i: " << i << " at (" << get_bulk().identifier(elements[i]) << "," << j << ") hostdata: " << hostData(i,j) << " dev data: " << data[j] << std::endl;
      }
    }
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
    stk::mesh::NgpField<int> ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

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
    int* data = stk::mesh::field_data(stkIntField, newElement);
    *data = get_bulk().identifier(newElement) * 10u;
    get_bulk().modification_end();
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

  bulk.get_entities(stk::topology::ELEM_RANK, bulk.mesh_meta_data().universal_part(), allElements);

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

#ifdef STK_USE_DEVICE_MESH
  size_t expectedSyncsToDevice = 1;
  size_t expectedSyncsToHost = 1;
#else
  size_t expectedSyncsToDevice = 0;
  size_t expectedSyncsToHost = 0;
#endif

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

#ifdef STK_USE_DEVICE_MESH
  size_t expectedSyncsToDevice = 1;
  size_t expectedSyncsToHost = 1;
#else
  size_t expectedSyncsToDevice = 0;
  size_t expectedSyncsToHost = 0;
#endif

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
  check_field_on_host(get_bulk(), stkIntField, multiplier);
}

template<typename T>
class NgpFieldTester : public stk::mesh::NgpField<T>
{
public:
  bool test_need_sync_to_host() const { return this->need_sync_to_host(); }
  bool test_need_sync_to_device() const { return this->need_sync_to_device(); }
};

TEST_F(NgpFieldFixture, ClearSyncStateAfterModifyOnDevice)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  NgpFieldTester<int>& testNgpField = static_cast<NgpFieldTester<int>&>(ngpField);

  testNgpField.modify_on_device();
#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(testNgpField.test_need_sync_to_host());
#else
  EXPECT_FALSE(testNgpField.test_need_sync_to_host());
#endif

  stkIntField.clear_sync_state();
  EXPECT_FALSE(testNgpField.test_need_sync_to_host());
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
}
