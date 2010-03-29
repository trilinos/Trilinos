/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_unit_tests_UnitTestFixture_hpp
#define stk_mesh_unit_tests_UnitTestFixture_hpp

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>

namespace sunit {

  class Stk_Mesh_Fixture
  {
    public:
      Stk_Mesh_Fixture( const std::vector<std::string> & entity_type_names );
      virtual ~Stk_Mesh_Fixture();

      void createMetaData(const std::vector<std::string> & entity_type_names );

      const stk::mesh::MetaData & get_MetaData() const { return m_MetaData ; }
      stk::mesh::MetaData       & get_NonconstMetaData() { return m_MetaData ; }

      void createBulkData();

      const stk::mesh::BulkData & get_BulkData() const { return m_BulkData ; }
      stk::mesh::BulkData       & get_NonconstBulkData() { return m_BulkData ; }

    private:
      stk::mesh::MetaData m_MetaData;
      stk::mesh::BulkData m_BulkData;

      Stk_Mesh_Fixture();
      Stk_Mesh_Fixture( const Stk_Mesh_Fixture & );
      Stk_Mesh_Fixture & operator = ( const Stk_Mesh_Fixture & );
  };
  
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

  class ExampleFixture : public Stk_Mesh_Fixture {
  public:
    ExampleFixture();
    ~ExampleFixture();
  };

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
