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

#include <stk_unit_test_utils/stk_mesh_fixtures/SelectorFixture.hpp>
#include <sstream>                      // for ostringstream, etc
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, put_field
#include "stk_mesh/base/Types.hpp"      // for EntityRank, EntityId, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/BuildMesh.hpp"
namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {
namespace fixtures {
using stk::unit_test_util::build_mesh;
using stk::mesh::MeshBuilder;

SelectorFixture::~SelectorFixture() {}

SelectorFixture::SelectorFixture()
  : m_meta_data_ptr(MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create_meta_data()),
    m_meta_data( *m_meta_data_ptr ),
    m_bulk_data( m_meta_data, MPI_COMM_WORLD ),
    m_partA( m_meta_data.declare_part( "PartA" , stk::topology::NODE_RANK ) ),
    m_partB( m_meta_data.declare_part( "PartB" , stk::topology::NODE_RANK ) ),
    m_partC( m_meta_data.declare_part( "PartC" , stk::topology::NODE_RANK ) ),
    m_partD( m_meta_data.declare_part( "PartD" , stk::topology::NODE_RANK ) ),
    m_entity1( ),
    m_entity2( ),
    m_entity3( ),
    m_entity4( ),
    m_entity5( )
{
  m_fieldA = &m_meta_data.declare_field<double>(stk::topology::NODE_RANK, "FieldA");
  m_fieldABC = &m_meta_data.declare_field<double>(stk::topology::NODE_RANK, "FieldABC");

  stk::mesh::put_field_on_mesh(*m_fieldA, m_partA, nullptr);

  stk::mesh::put_field_on_mesh(*m_fieldABC, m_partA, nullptr);
  stk::mesh::put_field_on_mesh(*m_fieldABC, m_partB, nullptr);
  stk::mesh::put_field_on_mesh(*m_fieldABC, m_partC, nullptr);
}

void SelectorFixture::generate_mesh()
{
  const unsigned entity_count = 5 ;

  // Create Entities and assign to parts:
  stk::mesh::EntityId ent_id =
    1 + entity_count * m_bulk_data.parallel_rank(); // Unique ID

  std::vector<stk::mesh::Part*> partMembership;

  // Entity1 is contained in PartA
  partMembership.clear();
  partMembership.push_back( & m_partA );
  m_entity1 = m_bulk_data.declare_node(ent_id, partMembership);
  ++ent_id;

  // Entity2 is contained in PartA and PartB
  partMembership.clear();
  partMembership.push_back( & m_partA );
  partMembership.push_back( & m_partB );
  m_entity2 = m_bulk_data.declare_node(ent_id, partMembership);
  ++ent_id;

  // Entity3 is contained in PartB and PartC
  partMembership.clear();
  partMembership.push_back( & m_partB );
  partMembership.push_back( & m_partC );
  m_entity3 = m_bulk_data.declare_node(ent_id, partMembership);
  ++ent_id;

  // Entity4 is contained in PartC
  partMembership.clear();
  partMembership.push_back( & m_partC );
  m_entity4 = m_bulk_data.declare_node(ent_id, partMembership);
  ++ent_id;

  // Entity5 is not contained in any Part
  partMembership.clear();
  m_entity5 = m_bulk_data.declare_node(ent_id, partMembership);
}

//--------------------------------------------------------------------------

VariableSelectorFixture::~VariableSelectorFixture() {}

VariableSelectorFixture::VariableSelectorFixture(int NumParts)
: m_MetaDataPtr(MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).set_entity_rank_names(std::vector<std::string>(4, std::string("MyEntityRank"))).create_meta_data())
, m_MetaData( *m_MetaDataPtr )
, m_BulkDataPtr( MeshBuilder(MPI_COMM_WORLD).create(m_MetaDataPtr) )
, m_BulkData( *m_BulkDataPtr )
, m_declared_part_vector()
{

  // Create Parts and commit:
  std::string myPartName;
  stk::mesh::EntityRank myRank = stk::topology::NODE_RANK;

  std::string partName = "Part_";
  for (int part_i=0 ; part_i<NumParts; ++part_i) {
    std::ostringstream localPartName(partName);
    localPartName << part_i;
    stk::mesh::Part * part =
      & m_MetaData.declare_part(localPartName.str(),myRank);
    m_declared_part_vector.push_back( part );
  }

  m_MetaData.commit();

  // Create Entities and assign to parts:

  m_BulkData.modification_begin();

  stk::mesh::EntityId ent_id =
    1 + NumParts * m_BulkData.parallel_rank(); // Unique ID

  for (int part_i = 0 ; part_i < NumParts ; ++part_i) {
    std::vector<stk::mesh::Part*> partMembership;
    partMembership.push_back(m_declared_part_vector[part_i]);
    stk::mesh::Entity e = m_BulkData.declare_node(ent_id, partMembership);
    m_entities.push_back( e );
    ++ent_id;
  }

  m_BulkData.modification_end();
}

namespace simple_fields {

SelectorFixture::~SelectorFixture() {}

SelectorFixture::SelectorFixture()
  : m_meta_data_ptr(MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create_meta_data()),
    m_meta_data( *m_meta_data_ptr ),
    m_bulk_data( m_meta_data, MPI_COMM_WORLD ),
    m_partA( m_meta_data.declare_part( "PartA" , stk::topology::NODE_RANK ) ),
    m_partB( m_meta_data.declare_part( "PartB" , stk::topology::NODE_RANK ) ),
    m_partC( m_meta_data.declare_part( "PartC" , stk::topology::NODE_RANK ) ),
    m_partD( m_meta_data.declare_part( "PartD" , stk::topology::NODE_RANK ) ),
    m_entity1( ),
    m_entity2( ),
    m_entity3( ),
    m_entity4( ),
    m_entity5( )
{
  m_fieldA = &m_meta_data.declare_field<double>(stk::topology::NODE_RANK, "FieldA");
  m_fieldABC = &m_meta_data.declare_field<double>(stk::topology::NODE_RANK, "FieldABC");

  stk::mesh::put_field_on_mesh(*m_fieldA, m_partA, nullptr);

  stk::mesh::put_field_on_mesh(*m_fieldABC, m_partA, nullptr);
  stk::mesh::put_field_on_mesh(*m_fieldABC, m_partB, nullptr);
  stk::mesh::put_field_on_mesh(*m_fieldABC, m_partC, nullptr);
}

void SelectorFixture::generate_mesh()
{
  const unsigned entity_count = 5 ;

  // Create Entities and assign to parts:
  stk::mesh::EntityId ent_id =
    1 + entity_count * m_bulk_data.parallel_rank(); // Unique ID

  std::vector<stk::mesh::Part*> partMembership;

  // Entity1 is contained in PartA
  partMembership.clear();
  partMembership.push_back( & m_partA );
  m_entity1 = m_bulk_data.declare_node(ent_id, partMembership);
  ++ent_id;

  // Entity2 is contained in PartA and PartB
  partMembership.clear();
  partMembership.push_back( & m_partA );
  partMembership.push_back( & m_partB );
  m_entity2 = m_bulk_data.declare_node(ent_id, partMembership);
  ++ent_id;

  // Entity3 is contained in PartB and PartC
  partMembership.clear();
  partMembership.push_back( & m_partB );
  partMembership.push_back( & m_partC );
  m_entity3 = m_bulk_data.declare_node(ent_id, partMembership);
  ++ent_id;

  // Entity4 is contained in PartC
  partMembership.clear();
  partMembership.push_back( & m_partC );
  m_entity4 = m_bulk_data.declare_node(ent_id, partMembership);
  ++ent_id;

  // Entity5 is not contained in any Part
  partMembership.clear();
  m_entity5 = m_bulk_data.declare_node(ent_id, partMembership);
}

//--------------------------------------------------------------------------

VariableSelectorFixture::~VariableSelectorFixture() {}

VariableSelectorFixture::VariableSelectorFixture(int NumParts)
: m_MetaDataPtr(MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).set_entity_rank_names(std::vector<std::string>(4, std::string("MyEntityRank"))).create_meta_data())
, m_MetaData( *m_MetaDataPtr )
, m_BulkDataPtr( MeshBuilder(MPI_COMM_WORLD).create(m_MetaDataPtr) )
, m_BulkData( *m_BulkDataPtr )
, m_declared_part_vector()
{

  // Create Parts and commit:
  std::string myPartName;
  stk::mesh::EntityRank myRank = stk::topology::NODE_RANK;

  std::string partName = "Part_";
  for (int part_i=0 ; part_i<NumParts; ++part_i) {
    std::ostringstream localPartName(partName);
    localPartName << part_i;
    stk::mesh::Part * part =
      & m_MetaData.declare_part(localPartName.str(),myRank);
    m_declared_part_vector.push_back( part );
  }

  m_MetaData.commit();

  // Create Entities and assign to parts:

  m_BulkData.modification_begin();

  stk::mesh::EntityId ent_id =
    1 + NumParts * m_BulkData.parallel_rank(); // Unique ID

  for (int part_i = 0 ; part_i < NumParts ; ++part_i) {
    std::vector<stk::mesh::Part*> partMembership;
    partMembership.push_back(m_declared_part_vector[part_i]);
    stk::mesh::Entity e = m_BulkData.declare_node(ent_id, partMembership);
    m_entities.push_back( e );
    ++ent_id;
  }

  m_BulkData.modification_end();
}

} // namespace simple_fields

} // fixtures
} // mesh
} // stk
