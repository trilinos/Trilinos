#ifndef NgpFieldTestUtils_hpp
#define NgpFieldTestUtils_hpp

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace ngp_field_test_utils {

template<typename T>
void check_field_data_on_device(stk::mesh::NgpMesh& ngpMesh,
                                stk::mesh::NgpField<T>& ngpField,
                                const stk::mesh::Selector& selector,
                                T expectedValue,
                                int component = -1,
                                T componentValue = 0)
{
  stk::mesh::for_each_entity_run(ngpMesh, ngpField.get_rank(), selector,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      const int numComponents = ngpField.get_num_components_per_entity(entity);
      for(int i=0; i<numComponents; ++i) {
        if (i == component) {
          STK_NGP_ThrowRequire(ngpField(entity, i) == componentValue);
        }
        else {
          STK_NGP_ThrowRequire(ngpField(entity, i) == expectedValue);
        }
      }
    }
  );
}

template<typename T>
void check_field_data_on_host(const stk::mesh::HostMesh& stkMesh,
                              const stk::mesh::FieldBase& stkField,
                              const stk::mesh::Selector& selector,
                              T expectedValue,
                              int component = -1,
                              T componentValue = 0)
{
  stk::mesh::for_each_entity_run(stkMesh, stkField.entity_rank(), selector,
    [&](const stk::mesh::FastMeshIndex& fastMeshIndex) {
      stk::mesh::Entity entity = stkMesh.get_entity(stkField.entity_rank(), fastMeshIndex);
      const int numComponents = stk::mesh::field_scalars_per_entity(stkField, entity);
      const T* fieldData = reinterpret_cast<const T*>(stk::mesh::field_data(stkField, entity));
      for(int i=0; i<numComponents; ++i) {
        if (i == component) {
          STK_ThrowRequireMsg(fieldData[i] == componentValue, "expected "<<componentValue<<" but actual value is "<<fieldData[i]<<"; i=="<<i<<", component=="<<component);
        }
        else {
          STK_ThrowRequireMsg(fieldData[i] == expectedValue, "expected "<<expectedValue<<" but actual value is "<<fieldData[i]<<"; i=="<<i<<", entity="<<stkField.entity_rank()<<","<<stkMesh.identifier(entity));
        }
      }
    }
  );
}

template<typename T>
void check_field_data_on_host(const stk::mesh::BulkData& stkMesh,
                              const stk::mesh::FieldBase& stkField,
                              const stk::mesh::Selector& selector,
                              T expectedValue,
                              int component = -1,
                              T componentValue = 0)
{
  stk::mesh::for_each_entity_run(stkMesh, stkField.entity_rank(), selector,
    [&](const stk::mesh::BulkData& bulk, const stk::mesh::Entity entity) {
      const int numComponents = stk::mesh::field_scalars_per_entity(stkField, entity);
      const T* fieldData = reinterpret_cast<const T*>(stk::mesh::field_data(stkField, entity));
      for(int i=0; i<numComponents; ++i) {
        if (i == component) {
          STK_ThrowRequireMsg(fieldData[i] == componentValue, "expected "<<componentValue<<" but actual value is "<<fieldData[i]<<"; i=="<<i<<", component=="<<component);
        }
        else {
          STK_ThrowRequireMsg(fieldData[i] == expectedValue, "expected "<<expectedValue<<" but actual value is "<<fieldData[i]<<"; i=="<<i<<", entity="<<bulk.entity_key(entity));
        }
      }
    }
  );
}

inline void set_field_data_on_host(const stk::mesh::BulkData& stkMesh,
                            const stk::mesh::FieldBase& stkField,
                            const stk::mesh::Selector& selector,
                            std::function<std::vector<double>(const double*)> func)
{
  const stk::mesh::FieldBase* coordField = stkMesh.mesh_meta_data().coordinate_field();
  stk::mesh::for_each_entity_run(stkMesh, stkField.entity_rank(), selector,
    [&](const stk::mesh::BulkData& bulk, const stk::mesh::Entity entity) {
      double* entityCoords = static_cast<double*>(stk::mesh::field_data(*coordField, entity));
      auto expectedValues = func(entityCoords);
      const int numComponents = stk::mesh::field_scalars_per_entity(stkField, entity);
      double* fieldData = static_cast<double*>(stk::mesh::field_data(stkField, entity));
      for(int i=0; i<numComponents; ++i) {
          fieldData[i] = expectedValues[i];
        }
      }
  );
}

inline void check_field_data_on_host_func(const stk::mesh::BulkData& stkMesh,
                            const stk::mesh::FieldBase& stkField,
                            const stk::mesh::Selector& selector,
                            const std::vector<const stk::mesh::FieldBase*>& otherFields,
                            std::function<std::vector<double>(const double*)> func)
{
  const stk::mesh::FieldBase* coordField = stkMesh.mesh_meta_data().coordinate_field();
  stk::mesh::for_each_entity_run(stkMesh, stkField.entity_rank(), selector,
    [&](const stk::mesh::BulkData& bulk, const stk::mesh::Entity entity) {
      double* entityCoords = static_cast<double*>(stk::mesh::field_data(*coordField, entity));
      auto expectedValues = func(entityCoords);

      unsigned int numComponents = stk::mesh::field_scalars_per_entity(stkField, entity);
      for (const stk::mesh::FieldBase* otherField : otherFields)
      {
        numComponents = std::min(numComponents, stk::mesh::field_scalars_per_entity(*otherField, entity));
      }
      const double* fieldData = reinterpret_cast<const double*>(stk::mesh::field_data(stkField, entity));
      for(unsigned int i=0; i<numComponents; ++i) {
          EXPECT_FLOAT_EQ(fieldData[i], expectedValues[i]);
        }
      }
  );
}

} // ngp_field_test_utils

#endif

