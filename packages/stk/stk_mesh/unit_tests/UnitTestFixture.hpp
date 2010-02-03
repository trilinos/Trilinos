#ifndef stk_mesh_unit_tests_UnitTestFixture_hpp
#define stk_mesh_unit_tests_UnitTestFixture_hpp

#include <Teuchos_RCP.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>

namespace sunit {

  class Stk_Mesh_Fixture
  {
    public:
      Stk_Mesh_Fixture();
      ~Stk_Mesh_Fixture();

      void createMetaData(const std::vector<std::string> & entity_type_names );

      Teuchos::RCP<const stk::mesh::MetaData> get_MetaData() const;
      Teuchos::RCP<stk::mesh::MetaData> get_NonconstMetaData();

      void createBulkData();

      Teuchos::RCP<const stk::mesh::BulkData> get_BulkData() const;
      Teuchos::RCP<stk::mesh::BulkData> get_NonconstBulkData();

    private:
      Teuchos::RCP<stk::mesh::MetaData> m_MetaData;
      Teuchos::RCP<stk::mesh::BulkData> m_BulkData;
  };
  
  // Nonmember constructor
  Teuchos::RCP<Stk_Mesh_Fixture> stk_Mesh_Fixture();
  

  stk::mesh::PartVector getPartVector(
      const Stk_Mesh_Fixture & fix,
      const std::vector<std::string> & names
      );

  stk::mesh::Part* getPart(
      const Stk_Mesh_Fixture & fix,
      const std::string name
      );

  const stk::mesh::Bucket & getBucketContainingEntity(
      const Stk_Mesh_Fixture & fix,
      stk::mesh::EntityType ent_type,
      stk::mesh::EntityId ent_id
      );


  // Example Mesh primarily used for testing Selectors:
  //
  // PartA, PartB, PartC, PartD
  // Entity1, Entity2, Entity3, Entity4
  // All entities are rank 0
  //
  // PartA contains Entity1, Entity2
  // PartB contains Entity2, Entity3
  // PartC contains Entity3, Entity4
  // PartD contains no entities
  // Entity5 is not contained in any user-defined Part
  // 
  // |----------|--|-------|--|----------|    |-------------|
  // |<--PartA---->|       |<--PartC---->|    |   PartD     |
  // |          |<---PartB--->|          |    |             |
  // |  1       |2 |       |3 |       4  | 5  |             |
  // |          |  |       |  |          |    |             |
  // |          |  |       |  |          |    |             | 
  // |----------|--|-------|--|----------|    |-------------|   
  //
  Teuchos::RCP<Stk_Mesh_Fixture> getExampleFixture();

  stk::mesh::Part* getExamplePart(
      const sunit::Stk_Mesh_Fixture & fix,
      const std::string name
      );

  const stk::mesh::Bucket & getExampleBucket(
      const Stk_Mesh_Fixture & fix,
      stk::mesh::EntityId ent_id
      );

} // namespace sunit


#endif // stk_mesh_unit_tests_UnitTestFixture_hpp
