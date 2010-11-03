/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_fixture_SelectorFixture_hpp
#define stk_mesh_fixture_SelectorFixture_hpp

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk {
namespace mesh {
namespace fixtures {

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

class SelectorFixture {
 public:
  SelectorFixture();
  ~SelectorFixture();

  const stk::mesh::MetaData & get_MetaData() const { return m_meta_data ; }
  stk::mesh::MetaData       & get_NonconstMetaData() { return m_meta_data ; }

  const stk::mesh::BulkData & get_BulkData() const { return m_bulk_data ; }
  stk::mesh::BulkData       & get_NonconstBulkData() { return m_bulk_data ; }

  stk::mesh::MetaData m_meta_data ;
  stk::mesh::BulkData m_bulk_data ;

  stk::mesh::Part & m_partA ;
  stk::mesh::Part & m_partB ;
  stk::mesh::Part & m_partC ;
  stk::mesh::Part & m_partD ;

  stk::mesh::Entity * m_entity1 ;
  stk::mesh::Entity * m_entity2 ;
  stk::mesh::Entity * m_entity3 ;
  stk::mesh::Entity * m_entity4 ;
  stk::mesh::Entity * m_entity5 ;

  void generate_mesh();

 private:
  SelectorFixture( const SelectorFixture & );
  SelectorFixture & operator = ( const SelectorFixture & );
};

class VariableSelectorFixture {
 public:
  VariableSelectorFixture(int NumParts);
  ~VariableSelectorFixture();

  stk::mesh::MetaData m_MetaData ;
  stk::mesh::BulkData m_BulkData ;

  stk::mesh::PartVector m_declared_part_vector;
  std::vector<stk::mesh::Entity*> m_entities ;

 private:
  VariableSelectorFixture( const VariableSelectorFixture & );
  VariableSelectorFixture & operator = ( const VariableSelectorFixture & );
};

} // fixtures
} // mesh
} // stk

#endif // stk_mesh_fixture_UnitTestFixture_hpp
