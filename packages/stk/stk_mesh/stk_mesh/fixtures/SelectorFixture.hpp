// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef stk_mesh_fixture_SelectorFixture_hpp
#define stk_mesh_fixture_SelectorFixture_hpp

#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Types.hpp>      // for PartVector
#include <vector>                       // for vector
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
namespace stk { namespace mesh { class Part; } }

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
// FieldA is defined on PartA
// FieldABC is defined on Parts A, B, C
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

  stk::mesh::Entity m_entity1 ;
  stk::mesh::Entity m_entity2 ;
  stk::mesh::Entity m_entity3 ;
  stk::mesh::Entity m_entity4 ;
  stk::mesh::Entity m_entity5 ;

  stk::mesh::Field<double>& m_fieldA;
  stk::mesh::Field<double>& m_fieldABC;

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
  std::vector<stk::mesh::Entity> m_entities ;

 private:
  VariableSelectorFixture( const VariableSelectorFixture & );
  VariableSelectorFixture & operator = ( const VariableSelectorFixture & );
};

} // fixtures
} // mesh
} // stk

#endif // stk_mesh_fixture_SelectorFixture_hpp
